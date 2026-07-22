"""Tier-2 MD-confirmation scorer report (§I.11.4 Stage 0 — T4-absorbed).

Scores finished confirmation-namespace run dirs (``*_tier2probe_conf*``,
the T3/T9/re-pilot matrix) on the two §I.11.4.1 surfaces:

1. **Leg / twin-parity table** (the §4m–§4r columns): trapped and
   suppressed fractions, ``n̄_det``, ``n₁``, ``W₁(MD, twin)`` against the
   committed pre-registration CSVs, and the bin-1 mean detected KE for MD
   and twin.
2. **Experimental table** (the T4 scoring surface): the RQ8 solvated
   (n >= 1 renormalized) histogram vs the committed abundance reference
   (``W₁``, ``n₁``, ``n₂``, ``n₁/n₂``), and the per-n mean-KE curve vs the
   committed ihe_ked reference under the **committed error model**
   (per-point sqrt(statErr² + sysErr²) + the calib/condition correlated
   bands profiled as unit-normal nuisances). The n = 1 bin scores under
   the ``N1_ANCHOR`` convention (I77; default ``"median"``), with the
   mean-legacy chi² printed alongside on every row.

Pure scorer — runs nothing, mutates nothing. Every physics convention
lives in ``i2_helium_md.postprocess.tier2_confirmation``; knob columns are
read from the authoritative ``cfg.json`` (the F3 convention), never parsed
from the tag.

Invocation: edit USER SETTINGS, then

    python scripts/post_processing/tier2_confirmation_score.py

Output: two text tables on stdout; CSVs when the ``SAVE_*_CSV_PATH``
settings are set.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    load_confirmation_run,
    load_twin_ke_curve,
    load_twin_prediction,
    score_histogram_vs_reference,
    score_ke_curve_vs_reference,
    wasserstein_between,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

RUNS_ROOT = PROJECT_ROOT / "data" / "runs"

# label -> (run dir name under RUNS_ROOT, twin config label). Default: the
# T9 leg-D matrix (the §4r certified read).
RUN_SELECTION: dict[str, tuple[str, str]] = {
    f"dc{i}": (
        f"9A_drag_shared_pure_cubic_N50_tier2probe_conf270_dc{i}",
        f"c{i}",
    )
    for i in (1, 2, 3, 4)
}

# Committed twin pre-registration CSVs + the leg key inside them.
TWIN_LEG = "d"
TWIN_PREDICTIONS_CSV = (
    RUNS_ROOT / "h2b_forward_model" / "h2b_leg_d_predictions.csv"
)
TWIN_KE_CSV = RUNS_ROOT / "h2b_forward_model" / "h2b_leg_d_ke.csv"

# Committed experimental references.
ABUNDANCE_REFERENCE_CSV = (
    PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
)
IHE_KED_REFERENCE_CSV = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked" / "IHe_KED_reference.csv"
)

# KE-curve scoring knobs (§I.11.4.1): solvated branch only, thin bins need
# >= MIN_KE_BIN_COUNT ions, sim-side SE widens the per-point sigma.
KE_N_MIN = 1
MIN_KE_BIN_COUNT = 2
INCLUDE_SIM_SE = True

# The n = 1 KE comparison anchor (I77 convention; design frozen 2026-07-21):
# "median" (operative default) scores MD's solvated core against the
# reference median_KE_eV at n = 1; "mean" is the pre-I77 legacy; "exclude"
# drops the n = 1 bin. The report always prints BOTH chi^2 columns —
# ke_chi2_prof (this anchor) and ke_chi2_mean_legacy — so the §4s/§4t
# recorded numbers stay reproducible on every row.
N1_ANCHOR = "median"

SAVE_LEG_CSV_PATH = None  # e.g. PROJECT_ROOT / "data" / "runs" / "conf_leg_table.csv"
SAVE_EXPERIMENT_CSV_PATH = None

# ---------------------------------------------------------------------------


def _nan_if_none(value: float | None) -> float:
    """None -> NaN missing-value stand-in; 0.0 is a value, not missing."""
    return float("nan") if value is None else value


def _fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def format_table(rows: list[dict[str, Any]]) -> str:
    """Fixed-width text table (the F3 idiom)."""
    if not rows:
        return "(no rows)"
    cols = list(rows[0].keys())
    cells = [[_fmt(r[c]) for c in cols] for r in rows]
    widths = [
        max(len(c), *(len(row[i]) for row in cells))
        for i, c in enumerate(cols)
    ]
    lines = [
        "  ".join(c.ljust(w) for c, w in zip(cols, widths)),
        "  ".join("-" * w for w in widths),
    ]
    lines += ["  ".join(v.ljust(w) for v, w in zip(row, widths)) for row in cells]
    return "\n".join(lines)


def write_rows_csv(path: str | Path, rows: list[dict[str, Any]]) -> Path:
    """Write the rows as CSV (the F3 idiom); returns the resolved path."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return p.resolve()


def _knob_columns(run_dir: Path) -> dict[str, Any]:
    """Knob columns from the authoritative cfg.json (never from the tag)."""
    cfg = json.loads((run_dir / "cfg.json").read_text(encoding="utf-8"))
    drag = cfg.get("drag_coefficients", {})
    coeff = drag.get("coefficients", {}) if isinstance(drag, dict) else {}
    return {
        "drag_form": cfg.get("drag_form", "-"),
        "v_c": coeff.get("v_c", "-"),
        "p_tail": coeff.get("p_tail", "-"),
        "tau_ps": cfg.get("internal_energy_cooling_tau_ps", "-"),
        "ladder": cfg.get("dissociation_ladder", "-"),
        "f_int": cfg.get("internal_energy_partition_fraction", "-"),
        "prior": cfg.get("droplet_size_prior", "-"),
    }


def experimental_row(
    label: str,
    twin_config: str,
    read,
    abundance_ref,
    ked_ref,
    *,
    n1_anchor: str = N1_ANCHOR,
    ke_n_min: int = KE_N_MIN,
    min_ke_bin_count: int = MIN_KE_BIN_COUNT,
    include_sim_se: bool = INCLUDE_SIM_SE,
) -> dict[str, Any]:
    """One experimental-table row: histogram score + the KE-curve score
    under the operative ``n1_anchor``, plus the mean-legacy chi^2 column
    (byte-identical to the pre-I77 scorer) for §4s/§4t reproducibility."""
    hist = score_histogram_vs_reference(read, abundance_ref)
    ke_score = score_ke_curve_vs_reference(
        read,
        ked_ref,
        n_min=ke_n_min,
        min_count=min_ke_bin_count,
        include_sim_se=include_sim_se,
        n1_anchor=n1_anchor,
    )
    if n1_anchor == "mean":
        ke_legacy_chi2 = ke_score.chi2_profiled
    else:
        ke_legacy_chi2 = score_ke_curve_vs_reference(
            read,
            ked_ref,
            n_min=ke_n_min,
            min_count=min_ke_bin_count,
            include_sim_se=include_sim_se,
            n1_anchor="mean",
        ).chi2_profiled
    return {
        "label": label,
        "config": twin_config,
        "w1_solv": hist.w1_solvated,
        "n1_solv": hist.sim_n1,
        "n2_solv": hist.sim_n2,
        "ratio_n1_n2": hist.ratio_n1_over_n2,
        "ref_n1": hist.ref_n1,
        "ref_n2": hist.ref_n2,
        "n1_anchor": ke_score.n1_anchor,
        "ke_chi2_prof": ke_score.chi2_profiled,
        "ke_chi2_mean_legacy": ke_legacy_chi2,
        "ke_npts": ke_score.n_points,
        "ke_a_calib": ke_score.a_calib,
        "ke_b_cond": ke_score.b_condition,
        "ke_within_1p25": ke_score.n_within_1p25,
        "ke_chi2_raw": ke_score.chi2_unprofiled,
    }


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    leg_rows: list[dict[str, Any]] = []
    exp_rows: list[dict[str, Any]] = []

    for label, (dir_name, twin_config) in RUN_SELECTION.items():
        run_dir = RUNS_ROOT / dir_name
        read = load_confirmation_run(run_dir, label=label)
        knobs = _knob_columns(run_dir)

        twin = load_twin_prediction(
            TWIN_PREDICTIONS_CSV, leg=TWIN_LEG, config=twin_config
        )
        twin_ke = load_twin_ke_curve(
            TWIN_KE_CSV, leg=TWIN_LEG, config=twin_config
        )
        ke = read.ke_by_n()
        n1_hits = ke.n == 1
        n1_ke_md = float(ke.mean_eV[n1_hits][0]) if n1_hits.any() else float("nan")

        leg_rows.append(
            {
                "label": label,
                "config": twin_config,
                **knobs,
                "num_scored": read.num_scored,
                "trapped": read.trapped_frac,
                "supp": read.suppressed_frac,
                "twin_supp": twin.suppressed_frac,
                "nbar_det": read.n_mean,
                "twin_nbar": twin.nbar_det,
                "n1": read.n1_frac,
                "twin_n1": float(twin.h[1]),
                "w1_md_twin": wasserstein_between(
                    read.n_values, read.fraction, read.n_values, twin.h
                ),
                "n1_ke_md_eV": n1_ke_md,
                "n1_ke_twin_eV": _nan_if_none(twin_ke.mean_ke_for(1)),
                "twin_trapped": twin.trapped_frac,
            }
        )

        exp_rows.append(
            experimental_row(label, twin_config, read, abundance_ref, ked_ref)
        )

    print("=== Leg / twin-parity table (conventions: §4r) ===")
    print(format_table(leg_rows))
    print()
    print(
        "=== Experimental table (RQ8 solvated branch; committed error "
        "model) ==="
    )
    print(format_table(exp_rows))

    if SAVE_LEG_CSV_PATH is not None:
        print(f"\nleg CSV -> {write_rows_csv(SAVE_LEG_CSV_PATH, leg_rows)}")
    if SAVE_EXPERIMENT_CSV_PATH is not None:
        print(
            f"experiment CSV -> "
            f"{write_rows_csv(SAVE_EXPERIMENT_CSV_PATH, exp_rows)}"
        )


if __name__ == "__main__":
    main()
