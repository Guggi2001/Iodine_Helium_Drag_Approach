"""Emit the Tier-2 I+He_n size-distribution scoreboard for finished campaign runs.

This is Slice F3: the deferred Phase-E assembler. It is a **pure scorer** -- it
composes only accepted A-E modules into **one scoreboard row per finished
campaign run dir** and runs nothing (Slice F2's ``gen_tier2_runs.py`` already
wrote each run's ``cfg.json`` / ``ion.npz`` / ``relaxation.npz``). Per run it
reports:

* **W1_matched** -- E1 terminal size distribution from ``relaxation.npz`` (the
  E2 matched-time result, ``source_tag="relaxed"``) scored by E4 against the
  experimental abundance reference;
* **W1_simend_upper** -- E1 from ``ion.npz`` (the R5 sim-end *upper bound*,
  ``source_tag="sim_end"``) scored by E4. Both are reported; never over-read the
  absolute terminal n (MASS doc Sec.R5);
* **E5 diagnostics** -- ``reconstruct_diagnostics(ion, cfg)`` -> the ensemble
  crossing time t*, the shell-retaining <-> total-strip regime label, and
  total-strip reachability. These are ion-stage quantities (t*/Pi need the ion
  trajectory), so the regime label reads the *ion-end* (upper-bound) n, distinct
  from the relaxed ``n_terminal_*`` columns below;
* **ledger_max_resid_eV** -- the 5-term ``ion_ledger_closure`` residual on
  ``ion.npz`` (a wiring diagnostic, not a physics verdict);
* **n_terminal_mean / n_terminal_spread** -- moments of the *matched* (relaxed)
  distribution.

**Reported, not auto-adjudicated** (the Tier-1a reporting-gate stance): no code
asserts a fidelity verdict. Adjudication is F4's / the user's job.

Run enumeration is **glob-discovery**: every complete ``*_tier2_*`` run dir under
``data/runs`` is scored, with the row's knob columns read from each authoritative
``cfg.json`` (never parsed from the tag). This keeps the scoreboard in sync with
the staged campaign (Stage-1 grid, then the Stage-2 tau sweep) without a settings
block to hand-resync between stages. ``budget_eV`` filters by the run's own
scenario stamp.

Usage::

    python scripts/post_processing/tier2_size_distribution_table.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Any

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =============================================================================
# USER SETTINGS
# =============================================================================

RUNS_ROOT = PROJECT_ROOT / "data" / "runs"
REFERENCE_PATH = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"

# Score only runs at this scenario budget [eV] (None -> every complete tier2 run).
BUDGET_EV: float | None = None

SAVE_CSV_PATH = None  # e.g. PROJECT_ROOT / "data" / "reference" / "drag" / "9A" / "tier2_scoreboard.csv"


from i2_helium_md.physics.constants import N_STAR  # noqa: E402
from i2_helium_md.postprocess import (  # noqa: E402
    HeAbundanceReference,
    compare_size_distributions,
    compute_terminal_shell_distribution,
    ion_ledger_closure,
    load_he_abundance_reference,
    reconstruct_diagnostics,
)
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.checkpoint import load_ion_checkpoint  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402

from scripts.tier2_common import TIER2_BRIDGE_TAG, VALIDATION_BUDGET_EV  # noqa: E402


# The 16 pinned plan-F3 columns, then the three confirmed advisory columns.
TIER2_TABLE_COLUMNS = [
    "case",
    "budget_eV",
    "picture",
    "kappa",
    "f_int",
    "f_ret",
    "tau_ps",
    "N",
    "W1_matched",
    "W1_simend_upper",
    "t_cross_ps",
    "regime_label",
    "total_strip_reachable",
    "ledger_max_resid_eV",
    "n_terminal_mean",
    "n_terminal_spread",
    # advisory (additive; strictly beyond the plan's pinned set):
    "t_cross_spread_ps",
    "t_cross_all_agree",
    "sanity_flags",
]

# Glob marker for a campaign run dir: the Slice-F1 run tag starts with ``tier2_``
# (``<case>_drag_<variant>_N<n>_tier2_...``). The ``total_strip`` variant appends
# a suffix but still carries the marker.
_TIER2_DIR_GLOB = "*_tier2_*"

# A campaign run dir is scorable only with these three artifacts. ``neutral.npz``
# is written by F2 but not read here, so it is not required for scoring.
_REQUIRED_ARTIFACTS: tuple[str, ...] = ("cfg.json", "ion.npz", "relaxation.npz")

_RELAXATION_FILENAME = "relaxation.npz"


@dataclass(frozen=True)
class Tier2RunRecord:
    """A scored campaign run: its dir, loaded cfg, and one scoreboard row."""

    label: str
    run_dir: Path
    cfg: SimConfig
    row: dict[str, Any]


def _freeze_side_expected(cfg: SimConfig, *, total_strip: bool = False) -> bool:
    """Whether the run sits at a freeze-side condition (E5 ``freeze_side_expected``).

    Decision (2026-07-05): key E5's ``freeze_side_expected`` off the run's own
    scenario stamp -- ``True`` at the 0.80 eV validation budget (where ``Pi > 1``
    is a units/wiring smell, MASS doc Sec.6.11) and ``False`` at the 2.70 eV
    production budget (where ``Pi > 1`` is legitimate shedding physics and the
    freeze-side reading must not be carried). Compared with a tolerance against
    the sanctioned validation constant, not a raw float literal.

    A ``total_strip`` variant is **always** shed-side: it is engineered (f_int
    high / tau long) so the gate never self-binds in-window, so ``Pi > 1`` is
    legitimate there regardless of budget -- flagging it would be spurious.
    """
    if total_strip:
        return False
    return bool(np.isclose(cfg.coulomb_available_eV, VALIDATION_BUDGET_EV))


def _is_total_strip_run(run_dir: Path) -> bool:
    """True for a ``total_strip`` regime-axis variant (the F1 ``_totalstrip`` tag).

    The variant is recorded only in the run tag (``tier2_run_tag(total_strip=
    True)`` appends ``_totalstrip``); there is no ``SimConfig`` field for it, so
    the marker is read from the dir-name suffix (symmetric to the ``case``
    prefix parse -- a structural marker, not a numeric knob).
    """
    return run_dir.name.endswith("_totalstrip")


def _case_from_run_dir(run_dir: Path) -> str:
    """Recover the droplet-geometry case from the run-dir basename.

    The Slice-F1 convention is ``<case>_drag_<variant>_N<n>_<run_tag>`` (see
    ``tier0_common.run_dir_name``); the case token is everything before the fixed
    ``_drag_`` separator. Only the case is parsed from the name -- every other row
    column is read from the authoritative ``cfg.json``.
    """
    name = run_dir.name
    sep = "_drag_"
    if sep not in name:
        raise ValueError(
            f"run dir {name!r} does not follow the '<case>_drag_...' convention; "
            "cannot recover the case token."
        )
    return name.split(sep, 1)[0]


def _distribution_moments(dist) -> tuple[float, float]:
    """Return ``(mean, spread)`` of a :class:`ShellDistribution` over integer n.

    Mean ``= sum n*fraction``; spread ``= sqrt(sum fraction*(n-mean)^2)`` (the
    population standard deviation of the size distribution). Consistent with the
    scored object rather than re-reading the raw checkpoint column.
    """
    n = dist.n_values.astype(float)
    p = dist.fraction.astype(float)
    mean = float(np.sum(n * p))
    spread = float(np.sqrt(np.sum(p * (n - mean) ** 2)))
    return mean, spread


def _run_is_complete(run_dir: Path) -> bool:
    """True when the run dir holds every artifact F3 reads."""
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def discover_run_dirs(runs_root: str | Path = RUNS_ROOT) -> list[Path]:
    """List every complete ``*_tier2_*`` campaign run dir under ``runs_root``.

    Discovery is by dir marker + artifact completeness only (no budget filter --
    that needs the loaded cfg and lives in :func:`collect_tier2_records`). A dir
    missing any of ``cfg.json`` / ``ion.npz`` / ``relaxation.npz`` (a crashed or
    in-flight run) is skipped, not scored.
    """
    root = Path(runs_root)
    if not root.is_dir():
        return []
    return sorted(
        p
        for p in root.glob(_TIER2_DIR_GLOB)
        if p.is_dir()
        and TIER2_BRIDGE_TAG not in p.name  # the reserved Phase-D bridge run
        and _run_is_complete(p)
    )


def _score_run(
    run_dir: Path,
    ref: HeAbundanceReference,
    *,
    n_max: int = N_STAR,
) -> dict[str, Any]:
    """Score one finished campaign run dir and return one scoreboard row."""
    run = RunDirectory(run_dir)
    cfg = run.load_cfg()
    ion = run.load_ion()
    relax = load_ion_checkpoint(run_dir / _RELAXATION_FILENAME)

    # E1 terminal size distributions: matched-time (relaxed) + sim-end upper
    # bound. A reloaded relaxation.npz is a v7 IonCheckpoint, so the matched
    # distribution MUST carry the explicit source_tag override (E1's documented
    # reload pitfall) or F3 would pair the two scores wrongly.
    sim_dist = compute_terminal_shell_distribution(
        ion, n_max=n_max, source_tag="sim_end"
    )
    relaxed_dist = compute_terminal_shell_distribution(
        relax, n_max=n_max, source_tag="relaxed"
    )

    # E4 integer-support Wasserstein via the config metric dispatch.
    metric = cfg.validation_histogram_metric
    w1_matched = compare_size_distributions(relaxed_dist, ref, metric=metric)
    w1_simend = compare_size_distributions(sim_dist, ref, metric=metric)

    # E5 regime diagnostics on the ion stage (t*/Pi need the ion trajectory);
    # freeze-side keyed off the run's own budget stamp (a total_strip variant is
    # shed-side regardless of budget).
    freeze_side = _freeze_side_expected(cfg, total_strip=_is_total_strip_run(run_dir))
    diag = reconstruct_diagnostics(ion, cfg, freeze_side_expected=freeze_side)

    closure = ion_ledger_closure(ion)
    n_mean, n_spread = _distribution_moments(relaxed_dist)

    return {
        "case": _case_from_run_dir(run_dir),
        "budget_eV": float(cfg.coulomb_available_eV),
        "picture": cfg.ladder_electronic_picture,
        "kappa": float(cfg.ladder_steepness),
        "f_int": _optional_float(cfg.internal_energy_partition_fraction),
        "f_ret": _optional_float(cfg.internal_energy_retained_fraction),
        "tau_ps": _optional_float(cfg.internal_energy_cooling_tau_ps),
        "N": int(cfg.num_molecules),
        "W1_matched": float(w1_matched),
        "W1_simend_upper": float(w1_simend),
        "t_cross_ps": diag.t_cross_summary.median_ps,
        "regime_label": diag.regime_label,
        "total_strip_reachable": bool(diag.total_strip_reachable),
        "ledger_max_resid_eV": float(closure.max_abs_residual_eV),
        "n_terminal_mean": n_mean,
        "n_terminal_spread": n_spread,
        "t_cross_spread_ps": diag.t_cross_summary.spread_ps,
        "t_cross_all_agree": bool(diag.t_cross_summary.all_agree),
        "sanity_flags": "; ".join(diag.sanity_flags),
    }


def _optional_float(value: float | None) -> float | None:
    """Coerce a set knob to ``float``; pass ``None`` (unset) through unchanged."""
    return None if value is None else float(value)


def score_tier2_run(
    run: RunDirectory | str | Path,
    ref: HeAbundanceReference,
    *,
    n_max: int = N_STAR,
) -> dict[str, Any]:
    """Score one finished campaign run and return one scoreboard row.

    ``W1_matched`` is the matched-time (relaxed) Wasserstein; ``W1_simend_upper``
    is the R5 sim-end upper bound. The ledger residual is a wiring diagnostic,
    not a physics verdict. Composes only accepted E1/E4/E5 modules.
    """
    run_dir = run.root if isinstance(run, RunDirectory) else Path(run)
    return _score_run(run_dir, ref, n_max=n_max)


def collect_tier2_records(
    *,
    runs_root: str | Path = RUNS_ROOT,
    reference_path: str | Path = REFERENCE_PATH,
    budget_eV: float | None = BUDGET_EV,
    n_max: int = N_STAR,
) -> tuple[list[Tier2RunRecord], HeAbundanceReference]:
    """Discover, load, and score every campaign run under ``runs_root``.

    The abundance reference is loaded once. ``budget_eV`` (when set) keeps only
    runs whose ``cfg.coulomb_available_eV`` matches (the run's own scenario
    stamp). Returns the scored records and the loaded reference.
    """
    ref = load_he_abundance_reference(reference_path)
    records: list[Tier2RunRecord] = []
    for run_dir in discover_run_dirs(runs_root):
        run = RunDirectory(run_dir)
        cfg = run.load_cfg()
        if budget_eV is not None and not np.isclose(
            cfg.coulomb_available_eV, budget_eV
        ):
            continue
        row = _score_run(run_dir, ref, n_max=n_max)
        records.append(
            Tier2RunRecord(
                label=run_dir.name,
                run_dir=run_dir,
                cfg=cfg,
                row=row,
            )
        )
    return records, ref


def collect_tier2_rows(
    *,
    runs_root: str | Path = RUNS_ROOT,
    reference_path: str | Path = REFERENCE_PATH,
    budget_eV: float | None = BUDGET_EV,
    n_max: int = N_STAR,
) -> list[dict[str, Any]]:
    """Discover, score, and return one scoreboard row per campaign run."""
    records, _ = collect_tier2_records(
        runs_root=runs_root,
        reference_path=reference_path,
        budget_eV=budget_eV,
        n_max=n_max,
    )
    return [record.row for record in records]


def _format_value(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def format_table(rows: list[dict[str, Any]]) -> str:
    """Return a plain-text table with the Tier-2 scoreboard columns."""
    if not rows:
        return "(no Tier-2 rows)"
    widths = {
        col: max(len(col), *(len(_format_value(row[col])) for row in rows))
        for col in TIER2_TABLE_COLUMNS
    }
    header = "  ".join(col.ljust(widths[col]) for col in TIER2_TABLE_COLUMNS)
    sep = "  ".join("-" * widths[col] for col in TIER2_TABLE_COLUMNS)
    body = [
        "  ".join(
            _format_value(row[col]).ljust(widths[col])
            for col in TIER2_TABLE_COLUMNS
        )
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def write_rows_csv(path: str | Path, rows: list[dict[str, Any]]) -> Path:
    """Write scoreboard rows to CSV using the stable Tier-2 column order."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=TIER2_TABLE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return p


def main() -> int:
    """Print the Tier-2 size-distribution scoreboard and optionally save CSV."""
    records, _ = collect_tier2_records(
        runs_root=RUNS_ROOT,
        reference_path=REFERENCE_PATH,
        budget_eV=BUDGET_EV,
    )
    rows = [record.row for record in records]
    print(
        "Tier-2 size-distribution scoreboard: W1_matched is the matched-time "
        "(relaxed) Wasserstein, W1_simend_upper the R5 sim-end upper bound; "
        "the ledger residual is a wiring diagnostic, not a physics verdict. "
        "Reported, not auto-adjudicated."
    )
    print(f"Scored {len(rows)} campaign run(s) under {RUNS_ROOT}.")
    print(format_table(rows))
    if SAVE_CSV_PATH is not None:
        path = write_rows_csv(SAVE_CSV_PATH, rows)
        print(f"Wrote Tier-2 scoreboard -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
