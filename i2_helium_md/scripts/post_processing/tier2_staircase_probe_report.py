"""Score the Tier-2 staircase capability probe against the anchored 21->19->14.

Pure scorer for the pre-F5 probe (``TIER2_STAIRCASE_PROBE_PLAN.md`` §4): it
runs nothing (``gen_tier2_staircase_probe.py`` already wrote each run's
``cfg.json`` / ``ion.npz`` / ``relaxation.npz``) and emits one row per finished
probe run dir:

* **staircase metrics** vs the anchored comparator ``build_shell_schedule``
  (zero artifact dependence -- the Phase-D bridge route): mean in-window sheds
  ``dn_mean_window`` (target ~ 7 = NUM_SHED_EVENTS), ion-end mean n
  ``n_ion_end_mean`` (target 14; the R5 sim-end read), mean first-shed time
  ``t_first_shed_ps`` vs the anchored first event, the shed fraction, and the
  time-averaged trajectory deviation ``n_traj_mad = <|mean_n(t) - n_anchor(t)|>``;
* **flexibility metrics** from the E2 relaxed checkpoint (E1 matched-time,
  ``source_tag="relaxed"`` -- the documented reload pitfall):
  ``n_relaxed_mean`` / ``n_relaxed_spread``. **Optional since Slice DS:**
  a skip-path run (relaxation stage disabled, detection seeding from
  ``ion.npz``) has no ``relaxation.npz`` and scores ``-`` here;
* **detected metrics** (Slice DS / Wave 7) from ``detection.npz`` **when
  present** (the artifact is optional -- probe dirs predating the detection
  stage score with ``-`` in these columns): ``n_detect_mean`` /
  ``n_detect_spread`` / ``n_detect_min`` plus the per-reason state fractions
  (``frozen`` / ``suppressed`` / ``time_exhausted``). The detected read is
  the Tier-2 arbitration observable; ``n_relaxed_*`` and ``frac_frozen``
  stay as convergence diagnostics beside it
  (``TIER2_DETECTION_STAGE_DESIGN.md`` §3.5);
* the 5-term ``ion_ledger_closure`` residual (wiring diagnostic, not physics).

**Reported, not auto-adjudicated**: the headline states factual minima (the
lowest-MAD grid point) and per-picture reachable relaxed terminal-n ranges;
no code asserts a capability verdict. TDDFT is not ground truth -- a grid-wide
miss is the RRK-dof mechanism-level OQ to surface (bridge findings §2), and a
landing region points the N=500 campaign; neither outcome discharges the F5
gate by itself.

Figures (mean-n(t) small multiples + relaxed terminal-n heat map) are gated
behind ``SAVE_FIGURES``/``SHOW_FIGURES`` and lazy-import matplotlib -- never
touched by pytest.

Usage::

    python scripts/post_processing/tier2_staircase_probe_report.py
"""

from __future__ import annotations

import csv
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

# Anchored comparator: schedule onset t* [ps] (the Tier-1a primary is 5.0).
T_STAR_PS = 5.0

SAVE_CSV_PATH = None  # e.g. PROJECT_ROOT / "data" / "reference" / "drag" / "9A" / "tier2_staircase_probe.csv"
SAVE_FIGURES = False  # figure output dir: RUNS_ROOT.parent / "tier2_probe_figures"
SHOW_FIGURES = True


from i2_helium_md.physics.constants import N_STAR  # noqa: E402
from i2_helium_md.physics.shell_schedule import (  # noqa: E402
    ANCHOR_N_END,
    NUM_SHED_EVENTS,
    ShellSchedule,
    build_shell_schedule,
)
from i2_helium_md.postprocess import (  # noqa: E402
    compute_terminal_shell_distribution,
    ion_ledger_closure,
)
from i2_helium_md.postprocess.derived_diagnostics import mean_shell_count  # noqa: E402
from i2_helium_md.simulation.checkpoint import load_ion_checkpoint  # noqa: E402
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    load_detection_result,
)
from i2_helium_md.simulation.ion_propagation_step import (  # noqa: E402
    ion_state_from_checkpoint_column,
)
from i2_helium_md.simulation.relaxation_stage import _freeze_mask  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


PROBE_TABLE_COLUMNS = [
    "case",
    "budget_eV",
    "picture",
    "kappa",
    "tau_ps",
    "s_eff",
    "cooling_gate",
    "f_int",
    "f_ret",
    "N",
    # staircase metrics (ion window, vs the anchored schedule)
    "dn_mean_window",
    "n_ion_end_mean",
    "t_first_shed_ps",
    "frac_ions_shed",
    "n_traj_mad",
    # flexibility / total-strip metrics (E2 relaxed, matched-time) -- a
    # convergence diagnostic once a detected read exists (design §3.5)
    "n_relaxed_mean",
    "n_relaxed_spread",
    "n_relaxed_min",     # deepest reachable shell (total-strip probe headline)
    "frac_frozen",       # E2 freeze completeness (n==0 or E_int<D_0(n))
    # detected metrics (Slice DS; from optional detection.npz -- "-" when absent)
    "n_detect_mean",
    "n_detect_spread",
    "n_detect_min",
    "frac_det_frozen",
    "frac_det_suppressed",
    "frac_det_time_exhausted",
    # wiring diagnostic
    "ledger_max_resid_eV",
]

# Glob marker for a probe run dir (disjoint from the F3 campaign ``*_tier2_*``
# glob in both directions -- the namespace lock).
_PROBE_DIR_GLOB = "*_tier2probe_*"

# A probe run dir is scorable with these two artifacts (``neutral.npz`` is
# written but not read here). ``relaxation.npz`` became OPTIONAL at Slice DS
# (review fix 2026-07-07): a skip-path run (relaxation_stage_enabled=False,
# detection stage seeding from ion.npz) carries no relaxation artifact and
# must still be discovered and scored -- its n_relaxed_*/frac_frozen columns
# read "-".
_REQUIRED_ARTIFACTS: tuple[str, ...] = ("cfg.json", "ion.npz")

_RELAXATION_FILENAME = "relaxation.npz"

# Optional detected-read artifact (Slice DS): scored when present, "-" when not
# -- deliberately NOT in _REQUIRED_ARTIFACTS so every pre-DS probe dir stays
# scorable unchanged (the Wave-7 back-compat criterion).
_DETECTION_FILENAME = "detection.npz"


def staircase_metrics(
    time_ps: np.ndarray,
    n_shell: np.ndarray,
    schedule: ShellSchedule,
) -> dict[str, float]:
    """Staircase-fidelity metrics of one run's ``n_shell`` vs the anchored schedule.

    Parameters
    ----------
    time_ps : np.ndarray, shape (Nt,)
        Stored ion-stage times [ps].
    n_shell : np.ndarray, shape (num_ions, Nt)
        Per-ion integer-valued shell counts (the v7 column).
    schedule : ShellSchedule
        The anchored comparator (``build_shell_schedule(t*)``).

    Returns
    -------
    dict
        ``dn_mean_window`` -- mean sheds in-window (target ``NUM_SHED_EVENTS``);De
        ``n_ion_end_mean`` -- ensemble-mean terminal (ion-end) n (target 14);
        ``t_first_shed_ps`` -- mean first-shed time over the ions that shed
        (NaN when no ion sheds -- never a crash);
        ``frac_ions_shed`` -- fraction of ions with at least one shed;
        ``n_traj_mad`` -- time-averaged ``|mean_n(t) - n_anchor(t)|`` [He].
    """
    time_ps = np.asarray(time_ps, dtype=float)
    n = np.asarray(n_shell, dtype=float)
    if n.ndim != 2 or n.shape[1] != time_ps.size:
        raise ValueError(
            f"n_shell must be (num_ions, Nt={time_ps.size}); got {n.shape}."
        )

    mean_n = mean_shell_count(n)
    n_anchor = np.asarray(schedule.n_of_t(time_ps), dtype=float)

    shed_mask = n < n[:, :1]  # below the ion's own initial count
    any_shed = shed_mask.any(axis=1)
    first_idx = np.argmax(shed_mask, axis=1)
    t_first = np.where(any_shed, time_ps[first_idx], np.nan)

    return {
        "dn_mean_window": float(mean_n[0] - mean_n[-1]),
        "n_ion_end_mean": float(mean_n[-1]),
        "t_first_shed_ps": (
            float(np.nanmean(t_first)) if any_shed.any() else float("nan")
        ),
        "frac_ions_shed": float(any_shed.mean()),
        "n_traj_mad": float(np.mean(np.abs(mean_n - n_anchor))),
    }


def _case_from_run_dir(run_dir: Path) -> str:
    """Recover the droplet-geometry case token (``<case>_drag_...`` convention)."""
    name = run_dir.name
    sep = "_drag_"
    if sep not in name:
        raise ValueError(
            f"run dir {name!r} does not follow the '<case>_drag_...' convention; "
            "cannot recover the case token."
        )
    return name.split(sep, 1)[0]


def _run_is_complete(run_dir: Path) -> bool:
    """True when the run dir holds every artifact the scorer reads."""
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def _s_eff_label(s_eff: float | None) -> str:
    """Render the RRK-dof knob for display: ``None`` -> ``per-n``, else compact."""
    return "per-n" if s_eff is None else f"{float(s_eff):g}"


def _cfg_matches_row(cfg, row: dict[str, Any]) -> bool:
    """True when a loaded cfg is the grid point a scored row describes.

    Used to re-find the best-case run's ion checkpoint for the overlay figure
    (rows carry only knobs, not the run dir). Matches the tag-encoded knobs:
    picture (exact), kappa/tau/budget (isclose), s_eff (both per-n, or isclose),
    and cooling_gate (exact). This mirrors the run-dir uniqueness the probe tag
    guarantees -- cooling_gate is a distinct tag dimension (``_cgds``), so without
    it an A/B pair at a fixed knob point would alias to the wrong arm's run dir.
    """
    row_s = row.get("s_eff")
    cfg_s = cfg.evap_rrk_dof
    if (row_s is None) != (cfg_s is None):
        return False
    if row_s is not None and not np.isclose(float(cfg_s), float(row_s)):
        return False
    if cfg.cooling_spatial_gate != row.get("cooling_gate", "none"):
        return False
    return (
        cfg.ladder_electronic_picture == row["picture"]
        and np.isclose(float(cfg.ladder_steepness), float(row["kappa"]))
        and np.isclose(float(cfg.internal_energy_cooling_tau_ps), float(row["tau_ps"]))
        and np.isclose(float(cfg.coulomb_available_eV), float(row["budget_eV"]))
    )


def discover_probe_run_dirs(runs_root: str | Path = RUNS_ROOT) -> list[Path]:
    """List every complete ``*_tier2probe_*`` run dir under ``runs_root``.

    The probe namespace is disjoint from the campaign glob by construction
    (``tier2_probe_run_tag``), so no bridge/campaign exclusion is needed here;
    an incomplete (crashed or in-flight) dir is skipped, not scored.
    """
    root = Path(runs_root)
    if not root.is_dir():
        return []
    return sorted(
        p for p in root.glob(_PROBE_DIR_GLOB) if p.is_dir() and _run_is_complete(p)
    )


def score_probe_run(
    run: RunDirectory | str | Path,
    schedule: ShellSchedule,
    *,
    n_max: int = N_STAR,
) -> dict[str, Any]:
    """Score one finished probe run dir into one report row.

    The staircase metrics read the **ion** stage (the anchored staircase is an
    in-window quantity; R5 sim-end read); the flexibility moments read the E2
    **relaxed** checkpoint (matched-time, ``source_tag="relaxed"``). Knob
    columns come from the authoritative ``cfg.json``, never the dir tag.
    """
    run_dir = run.root if isinstance(run, RunDirectory) else Path(run)
    run = RunDirectory(run_dir)
    cfg = run.load_cfg()
    ion = run.load_ion()

    metrics = staircase_metrics(ion.time_ps, ion.n_shell, schedule)

    # Relaxed read (E2) -- optional since Slice DS: a skip-path run carries no
    # relaxation.npz; its relaxed columns read "-" (review fix 2026-07-07).
    n_relaxed_mean = n_relaxed_spread = n_relaxed_min = frac_frozen = None
    relaxation_path = run_dir / _RELAXATION_FILENAME
    if relaxation_path.exists():
        relax = load_ion_checkpoint(relaxation_path)
        relaxed_dist = compute_terminal_shell_distribution(
            relax, n_max=n_max, source_tag="relaxed"
        )
        n_relaxed_mean, n_relaxed_spread = relaxed_dist.moments()
        # Deepest reachable shell = smallest occupied n in the relaxed
        # distribution (the total-strip A/B headline; read from the SAME
        # distribution as the moments).
        occupied = relaxed_dist.n_values[relaxed_dist.counts > 0]
        n_relaxed_min = int(occupied.min()) if occupied.size else int(n_max)

        # E2 freeze completeness: the fraction of ions at the relaxed cascade
        # floor (n==0 or E_int<D_0(n)); reuse relaxation_stage._freeze_mask on
        # the final relaxed column so a non-terminating gated cascade (never
        # all-frozen within the relaxation cap) is surfaced rather than
        # silently read as a converged terminal n.
        relaxed_state = ion_state_from_checkpoint_column(relax, -1)
        frozen = _freeze_mask(
            relaxed_state,
            picture=cfg.ladder_electronic_picture,
            kappa=float(cfg.ladder_steepness),
        )
        frac_frozen = float(np.mean(frozen))

    # Detected read (Slice DS / Wave 7): optional -- probe dirs predating the
    # detection stage carry no detection.npz and score "-" in these columns.
    detect_cols: dict[str, Any] = {
        "n_detect_mean": None, "n_detect_spread": None, "n_detect_min": None,
        "frac_det_frozen": None, "frac_det_suppressed": None,
        "frac_det_time_exhausted": None,
    }
    detection_path = run_dir / _DETECTION_FILENAME
    if detection_path.exists():
        detected = load_detection_result(detection_path)
        # Stale-artifact coherence guard (review fix 2026-07-07): a probe dir
        # regenerated in place with a leftover detection.npz would otherwise
        # score the OLD ensemble's detected read next to the NEW run's other
        # columns -- fail loud on the ensemble-size mismatch instead.
        if detected.num_molecules != int(cfg.num_molecules):
            raise ValueError(
                f"stale detection.npz in {run_dir}: it holds "
                f"num_molecules={detected.num_molecules} but cfg.json says "
                f"{cfg.num_molecules}. Re-run the detection stage for this "
                "dir (stale-artifact policy)."
            )
        detected_dist = compute_terminal_shell_distribution(detected, n_max=n_max)
        n_detect_mean, n_detect_spread = detected_dist.moments()
        det_occupied = detected_dist.n_values[detected_dist.counts > 0]
        fractions = detected.reason_fractions()
        detect_cols = {
            "n_detect_mean": n_detect_mean,
            "n_detect_spread": n_detect_spread,
            "n_detect_min": (
                int(det_occupied.min()) if det_occupied.size else int(n_max)
            ),
            "frac_det_frozen": fractions["frozen"],
            "frac_det_suppressed": fractions["suppressed"],
            "frac_det_time_exhausted": fractions["time_exhausted"],
        }

    closure = ion_ledger_closure(ion)

    return {
        "case": _case_from_run_dir(run_dir),
        "budget_eV": float(cfg.coulomb_available_eV),
        "picture": cfg.ladder_electronic_picture,
        "kappa": float(cfg.ladder_steepness),
        "tau_ps": float(cfg.internal_energy_cooling_tau_ps),
        "s_eff": None if cfg.evap_rrk_dof is None else float(cfg.evap_rrk_dof),
        "cooling_gate": cfg.cooling_spatial_gate,
        "f_int": float(cfg.internal_energy_partition_fraction),
        "f_ret": float(cfg.internal_energy_retained_fraction),
        "N": int(cfg.num_molecules),
        **metrics,
        "n_relaxed_mean": n_relaxed_mean,
        "n_relaxed_spread": n_relaxed_spread,
        "n_relaxed_min": n_relaxed_min,
        "frac_frozen": frac_frozen,
        **detect_cols,
        "ledger_max_resid_eV": float(closure.max_abs_residual_eV),
    }


def collect_probe_rows(
    *,
    runs_root: str | Path = RUNS_ROOT,
    t_star_ps: float = T_STAR_PS,
) -> list[dict[str, Any]]:
    """Discover, score, and return one row per complete probe run."""
    schedule = build_shell_schedule(t_star_ps)
    return [
        score_probe_run(run_dir, schedule)
        for run_dir in discover_probe_run_dirs(runs_root)
    ]


def _format_value(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def format_table(rows: list[dict[str, Any]]) -> str:
    """Return a plain-text table with the probe report columns."""
    if not rows:
        return "(no staircase-probe rows)"
    widths = {
        col: max(len(col), *(len(_format_value(row[col])) for row in rows))
        for col in PROBE_TABLE_COLUMNS
    }
    header = "  ".join(col.ljust(widths[col]) for col in PROBE_TABLE_COLUMNS)
    sep = "  ".join("-" * widths[col] for col in PROBE_TABLE_COLUMNS)
    body = [
        "  ".join(
            _format_value(row[col]).ljust(widths[col]) for col in PROBE_TABLE_COLUMNS
        )
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def select_best_case_row(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Return the probe row with the lowest trajectory MAD (the 'best case').

    This is the single definition of "best" the report uses -- both the text
    headline's factual argmin line and the best-case overlay figure draw the
    row it returns, so they agree by construction. "Best" = closest match to
    the anchored 21->19->14 staircase (lowest ``n_traj_mad``); at the mini-probe
    landing that row also has ``n_ion_end_mean`` ~ 14. Rows with a non-finite
    MAD (no ion trajectory scored) are ignored; ``None`` when nothing is finite.
    """
    finite = [row for row in rows if np.isfinite(row["n_traj_mad"])]
    if not finite:
        return None
    return min(finite, key=lambda row: float(row["n_traj_mad"]))


def format_headline(rows: list[dict[str, Any]]) -> str:
    """Factual probe headline: per-picture reach + the lowest-MAD grid point.

    Factual minima only (the F4 stance): the argmin is *named*, never declared
    a calibration success -- adjudication is the user's job.
    """
    lines = [
        "Staircase capability probe -- reported, not auto-adjudicated "
        "(existence probe; TDDFT is not ground truth).",
        f"targets: {NUM_SHED_EVENTS} sheds to n = {ANCHOR_N_END} in-window.",
    ]
    if not rows:
        lines.append("(no probe rows)")
        return "\n".join(lines)

    lines.append("reachable relaxed terminal-n range per picture:")
    for picture in sorted({row["picture"] for row in rows}):
        values = [
            float(row["n_relaxed_mean"])
            for row in rows
            if row["picture"] == picture and np.isfinite(row["n_relaxed_mean"])
        ]
        if values:
            lines.append(
                f"  {picture}: n_relaxed_mean in [{min(values):.1f}, {max(values):.1f}]"
            )
    best = select_best_case_row(rows)
    if best is not None:
        lines.append(
            "lowest trajectory MAD (factual argmin, not a verdict): "
            f"picture={best['picture']} k={float(best['kappa']):.2f} "
            f"tau={float(best['tau_ps']):.2f} s_eff={_s_eff_label(best.get('s_eff'))} -> "
            f"MAD={float(best['n_traj_mad']):.2f} He"
        )
    return "\n".join(lines)


def write_rows_csv(path: str | Path, rows: list[dict[str, Any]]) -> Path:
    """Write probe rows to CSV using the stable column order."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=PROBE_TABLE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    return p


def emit_figures(rows: list[dict[str, Any]], out_dir: Path) -> list[Path]:
    """Write the probe overlay figures (non-test path; lazy matplotlib import).

    (a) mean-n(t) small multiples: one panel per (picture, tau) cell, one curve
    per (kappa, s_eff), the anchored staircase underlaid. Requires re-loading
    each run's ion checkpoint, so this walks the run dirs again (figures are an
    operator artifact, not part of the scored rows).
    (b) relaxed terminal-n heat map (one panel per picture). The column axis is
    the **swept dof knob**: ``kappa`` when it varies (the delivered 45-run
    probe), else ``s_eff`` (the Addendum-A mini-probe pins ``kappa`` and sweeps
    ``s_eff``). Keying on the actually-swept knob stops an ``s_eff`` sweep at a
    fixed ``(kappa, tau)`` collapsing into one cell. If both ``kappa`` and
    ``s_eff`` vary (no probe does), ``kappa`` wins the axis and the scored
    **table / CSV** is the authoritative read-out.
    (c) best-case mean-n(t) overlay: the single lowest-MAD grid point
    (:func:`select_best_case_row` -- the same 'best' the headline names; at the
    mini-probe landing this is the n_end ~ 14 case), drawn like the Phase-D
    bridge ``bridge_mean_n_overlay`` -- the generative mean n(t) against the
    anchored 21->19->14 step schedule.
    """
    import matplotlib

    if not SHOW_FIGURES:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    schedule = build_shell_schedule(T_STAR_PS)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    run_dirs = discover_probe_run_dirs(RUNS_ROOT)
    pictures = sorted({row["picture"] for row in rows})
    taus = sorted({float(row["tau_ps"]) for row in rows})
    kappas = sorted({float(row["kappa"]) for row in rows})
    # Cooling-gate A/B axis: figures collapse to the pre-arm single-gate layout when
    # only one gate is present (so delivered non-A/B figures stay byte-identical).
    gates = sorted({row.get("cooling_gate", "none") for row in rows})

    # (a) mean-n(t) small multiples.
    fig, axes = plt.subplots(
        len(pictures),
        len(taus),
        figsize=(4.0 * len(taus), 3.0 * len(pictures)),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    for run_dir in run_dirs:
        run = RunDirectory(run_dir)
        cfg = run.load_cfg()
        ion = run.load_ion()
        try:
            i = pictures.index(cfg.ladder_electronic_picture)
            j = taus.index(float(cfg.internal_energy_cooling_tau_ps))
        except ValueError:
            continue
        s_label = (
            "s=per-n" if cfg.evap_rrk_dof is None else f"s={cfg.evap_rrk_dof:.0f}"
        )
        cg = cfg.cooling_spatial_gate
        cg_label = "" if len(gates) <= 1 else f" cg={cg}"
        cg_ls = "--" if (len(gates) > 1 and cg == "density_scaled") else "-"
        axes[i][j].plot(
            ion.time_ps,
            mean_shell_count(ion.n_shell),
            lw=1.2,
            ls=cg_ls,
            label=f"k={cfg.ladder_steepness:.2f} {s_label}{cg_label}",
        )
    for i, picture in enumerate(pictures):
        for j, tau_ps in enumerate(taus):
            ax = axes[i][j]
            t_grid = np.linspace(0.0, 30.0, 601)
            ax.step(
                t_grid,
                schedule.n_of_t(t_grid),
                where="post",
                color="k",
                ls="--",
                lw=1.5,
                label=f"anchored t*={T_STAR_PS:g}",
            )
            ax.set_title(f"{picture}, tau={tau_ps:g} ps", fontsize=9)
            if i == len(pictures) - 1:
                ax.set_xlabel("time [ps]")
            if j == 0:
                ax.set_ylabel("mean n")
            ax.legend(fontsize=6)
    fig.suptitle("Staircase probe: mean n(t) per grid cell vs anchored 21->19->14")
    fig.tight_layout()
    path = out_dir / "probe_mean_n_small_multiples.png"
    fig.savefig(path, dpi=150)
    saved.append(path)

    # (b) relaxed terminal-n heat map: tau x (swept dof knob) per picture. The
    # column knob is kappa when it varies, else s_eff -- so the mini-probe's
    # s_eff sweep is not collapsed into one (kappa, tau) cell.
    if len(kappas) > 1:
        col_values: list = kappas
        col_of = lambda row: float(row["kappa"])
        col_ticklabels = [f"{k:g}" for k in kappas]
        col_label = "kappa"
    else:
        col_values = sorted(
            {row["s_eff"] for row in rows},
            key=lambda s: float("inf") if s is None else float(s),  # per-n last
        )
        col_of = lambda row: row["s_eff"]
        col_ticklabels = [_s_eff_label(s) for s in col_values]
        col_label = "s_eff"

    # Rows facet the cooling-gate A/B (single row when only one gate is present, so
    # the two arms no longer collide into one (tau, col) cell and overwrite).
    fig, axes = plt.subplots(
        len(gates),
        len(pictures),
        figsize=(4.5 * len(pictures), 3.6 * len(gates)),
        squeeze=False,
    )
    for g_idx, gate in enumerate(gates):
        for p_idx, picture in enumerate(pictures):
            grid = np.full((len(taus), len(col_values)), np.nan)
            for row in rows:
                if row["picture"] != picture:
                    continue
                if row.get("cooling_gate", "none") != gate:
                    continue
                grid[
                    taus.index(float(row["tau_ps"])), col_values.index(col_of(row))
                ] = row["n_relaxed_mean"]
            ax = axes[g_idx][p_idx]
            im = ax.imshow(grid, origin="lower", aspect="auto", cmap="viridis")
            ax.set_xticks(range(len(col_values)), col_ticklabels)
            ax.set_yticks(range(len(taus)), [f"{t:g}" for t in taus])
            if g_idx == len(gates) - 1:
                ax.set_xlabel(col_label)
            if p_idx == 0:
                ax.set_ylabel("tau [ps]")
            ax.set_title(
                picture if len(gates) == 1 else f"{picture} / {gate}", fontsize=9
            )
            fig.colorbar(im, ax=ax, label="relaxed mean n")
    fig.suptitle("Staircase probe: relaxed terminal n across the grid")
    fig.tight_layout()
    path = out_dir / "probe_relaxed_n_heatmap.png"
    fig.savefig(path, dpi=150)
    saved.append(path)

    # (c) best-case single-panel overlay (the n_end ~ 14 landing), mirroring the
    # Phase-D bridge bridge_mean_n_overlay: generative mean n(t) vs the anchored
    # 21->19->14 step schedule. Re-load just the best run's ion checkpoint.
    best = select_best_case_row(rows)
    best_ion = None
    if best is not None:
        for run_dir in run_dirs:
            run = RunDirectory(run_dir)
            if _cfg_matches_row(run.load_cfg(), best):
                best_ion = run.load_ion()
                break
    if best_ion is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        t_grid = np.linspace(0.0, 30.0, 601)
        ax.step(
            t_grid, schedule.n_of_t(t_grid), where="post",
            color="k", ls="--", lw=1.5,
            label=f"anchored 21→19→14 (t*={T_STAR_PS:g} ps)",
        )
        ax.plot(
            best_ion.time_ps, mean_shell_count(best_ion.n_shell),
            "b-", lw=2.0, label="generative mean n(t)",
        )
        ax.set_xlabel("time [ps]")
        ax.set_ylabel("He-shell count n")
        ax.set_title(
            f"Best case: {best['picture']}, k={float(best['kappa']):.2f}, "
            f"tau={float(best['tau_ps']):.2f} ps, s_eff={_s_eff_label(best.get('s_eff'))}\n"
            f"n_end={float(best['n_ion_end_mean']):.1f} (target {ANCHOR_N_END}), "
            f"trajectory MAD={float(best['n_traj_mad']):.2f} He",
            fontsize=10,
        )
        ax.legend()
        fig.tight_layout()
        path = out_dir / "probe_best_case_mean_n_overlay.png"
        fig.savefig(path, dpi=150)
        saved.append(path)

    if SHOW_FIGURES:
        plt.show()
    plt.close("all")
    return saved


def main() -> int:
    """Print the probe table + headline; optionally save CSV and figures."""
    rows = collect_probe_rows(runs_root=RUNS_ROOT, t_star_ps=T_STAR_PS)
    print(
        f"Tier-2 staircase capability probe: scored {len(rows)} run(s) under "
        f"{RUNS_ROOT} (anchored comparator t*={T_STAR_PS:g} ps)."
    )
    print(format_table(rows))
    print()
    print(format_headline(rows))
    if SAVE_CSV_PATH is not None:
        path = write_rows_csv(SAVE_CSV_PATH, rows)
        print(f"Wrote probe report -> {path}")
    if (SAVE_FIGURES or SHOW_FIGURES) and rows:
        saved = emit_figures(rows, RUNS_ROOT.parent / "tier2_probe_figures")
        for p in saved:
            print(f"Wrote figure -> {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
