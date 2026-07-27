"""Atlas Axis A (G1) geometry-grid scorer report.

Scores the finished `gen_tier2atlas_geometry.py` cells (11 fixed-geometry
N = 500 runs) on the committed observable vector, plus the two geometry columns
the axis is about (realized droplet radius and **measured** mean birth depth).

Pure scorer — runs nothing, mutates nothing. Every physics convention lives in
`i2_helium_md.postprocess.tier2_confirmation`; the table helpers are imported
from `tier2_confirmation_score.py` rather than re-implemented (rule 1). Knob
columns come from the authoritative `cfg.json`, never from the run tag.

**Oracle first (§1.4).** The pooled N = 5000 battery row is scored before the
grid and printed against its recorded §4cc / Axis-A-pre-read values; a
mismatch means scorer drift and invalidates the session.

**Reading rule (plan §3.1, mandatory).** The landing is a *mixture* property:
pooled W₁ 0.571 beats every single geometry bin. These fixed-geometry cells are
therefore compared **to each other**, never directly against the committed
acceptance. Per-cell W₁ additionally carries ≈ 0.1 of single-seed scatter at
N = 500 (§6.7 item-2 finding 4), so it is reported but is **not** a
discriminator below |Δ| ≈ 0.15; the discriminating columns are `supp`, `trap`,
`nbar_det`, `n1_solv`, `midHot` and `deepKE`.

Invocation: edit USER SETTINGS, then

    python scripts/post_processing/tier2atlas_geometry_table.py

Output: the oracle row + the geometry table on stdout; CSV when
``SAVE_CSV_PATH`` is set.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    DEEP_KE_BAND,
    MIDHOT_BAND,
    deep_bin_ke_ratio,
    load_confirmation_run,
    midhot_ratio,
    score_histogram_vs_reference,
    score_ke_curve_vs_reference,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

RUNS_ROOT = PROJECT_ROOT / "data" / "runs"

# The grid cells, in grid order (label -> run dir basename).
CELL_LABELS: tuple[str, ...] = (
    "r1l1", "r1l2", "r1l3",
    "r2l1", "r2l2", "r2l3",
    "r3l1", "r3l2", "r3l3",
    "r4l2", "r4l3",
)
RUN_SELECTION: dict[str, str] = {
    label: f"9A_drag_shared_pure_cubic_N500_tier2atlas_conf270_geo{label}"
    for label in CELL_LABELS
}

# Scorer-drift oracle: the pooled standing battery + its recorded values.
ORACLE_RUN = "9A_drag_shared_pure_cubic_N5000_tier2probe_conf270_bigc1v725pooled"
ORACLE_RECORDED: dict[str, float] = {
    "trap": 0.067, "supp": 0.187, "nbar_det": 4.068,
    "n1_solv": 0.243, "w1_solv": 0.571, "midHot": 1.0110, "deepKE": 0.631,
}
ORACLE_TOL = 0.002

ABUNDANCE_REFERENCE_CSV = (
    PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
)
IHE_KED_REFERENCE_CSV = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked" / "IHe_KED_reference.csv"
)

# KE-curve scoring knobs — the committed §I.11.4.1 / I77 conventions.
KE_N_MIN = 1
MIN_KE_BIN_COUNT = 2
INCLUDE_SIM_SE = True
N1_ANCHOR = "median"

SAVE_CSV_PATH = None  # e.g. RUNS_ROOT / "atlas_geometry_grid.csv"

# ---------------------------------------------------------------------------


def _cfg(run_dir: Path) -> dict[str, Any]:
    return json.loads((run_dir / "cfg.json").read_text(encoding="utf-8"))


def geometry_columns(run_dir: Path) -> dict[str, Any]:
    """Cell geometry: the configured law + the **measured** birth statistics.

    ``R`` is exact (fixed size per cell), while ``depth_mean`` is read from the
    neutral checkpoint's stored ``r0`` — the quantity the pre-read identified as
    the physics knob, so it is measured rather than assumed. Under a center-pin
    cell ``r0`` is sampled and then zeroed by the driver, so the reported depth
    is R, matching what the propagation actually saw.
    """
    cfg = _cfg(run_dir)
    from i2_helium_md.simulation.checkpoint import load_neutral_checkpoint

    ckpt = load_neutral_checkpoint(run_dir / "neutral.npz")
    radii = np.asarray(ckpt.droplet_radii, dtype=float)
    R = float(radii[0])
    if not np.allclose(radii, R):
        raise ValueError(
            f"{run_dir.name}: droplet radii are not constant "
            f"({radii.min():.3f}–{radii.max():.3f} Å) — this report assumes the "
            f"fixed-size grid cells."
        )
    r0 = np.zeros(int(ckpt.num_molecules)) if cfg["single_initial_position"] \
        else np.asarray(ckpt.r0, dtype=float)
    law = cfg["birth_position_law"]
    if cfg["single_initial_position"]:
        law = "center"
    return {
        "R_A": R,
        "law": law,
        "margin_A": cfg["initial_position_margin_angstrom"],
        "well_K": cfg["binding_energy_molecule_K"],
        "N_he": cfg["single_droplet_size"],
        "r0_mean_A": float(r0.mean()),
        "depth_mean_A": float((R - r0).mean()),
    }


def observable_row(label: str, run_dir: Path, abundance_ref, ked_ref,
                   *, with_geometry: bool = True) -> dict[str, Any]:
    """One row of the committed observable vector for a finished run."""
    read = load_confirmation_run(run_dir, label=label)
    hist = score_histogram_vs_reference(read, abundance_ref)
    ke = score_ke_curve_vs_reference(
        read, ked_ref, n_min=KE_N_MIN, min_count=MIN_KE_BIN_COUNT,
        include_sim_se=INCLUDE_SIM_SE, n1_anchor=N1_ANCHOR,
    )
    mid = midhot_ratio(read, ked_ref)
    deep = deep_bin_ke_ratio(read, ked_ref)
    row: dict[str, Any] = {"label": label}
    if with_geometry:
        row.update(geometry_columns(run_dir))
    row.update({
        "num_scored": read.num_scored,
        "trap": read.trapped_frac,
        # Decomposed on purpose (plan §3.5b item 2): trap_bound is the exact
        # conservative verdict (physics — the ion cannot escape at any
        # relaxation length), trap_marginal is a MODELLING exclusion under
        # `exclude_all_coupled`, conditional on the drag law and the window.
        # Never merge them in a table; `trap` is their sum.
        "trap_bound": read.trap_bound_frac,
        "trap_marg": read.trap_marginal_frac,
        "supp": read.suppressed_frac,
        "nbar_det": read.n_mean,
        "n1_solv": hist.sim_n1,
        "w1_solv": hist.w1_solvated,
        "midHot": mid.value,
        "midHot_bins": mid.n_bins,
        "deepKE": deep.value,
        "deepKE_bins": deep.n_bins,
        "chi2_med": ke.chi2_profiled,
        # Reported next to chi2 on purpose: the deep-birth cells vacate the
        # low-n bins, so their chi2 is a sum over far fewer terms and a smaller
        # value there is NOT a better KE landing.
        "ke_npts": ke.n_points,
    })
    return row


def check_oracle(row: dict[str, Any]) -> list[str]:
    """Return the list of oracle columns that drifted beyond ORACLE_TOL."""
    return [
        f"{key}: got {row[key]:.4f} vs recorded {expected:.4f}"
        for key, expected in ORACLE_RECORDED.items()
        if not abs(float(row[key]) - expected) <= ORACLE_TOL
    ]


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    print(f"midHot = geometric mean of sim/ref-mean over n = "
          f"{MIDHOT_BAND[0]}–{MIDHOT_BAND[1]}; deepKE = arithmetic mean over "
          f"n = {DEEP_KE_BAND[0]}–{DEEP_KE_BAND[1]} (committed conventions).")
    print()

    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    print("=== Scorer-drift oracle: the pooled standing battery ===")
    print(format_table([oracle_row]))
    drift = check_oracle(oracle_row)
    if drift:
        print("\n*** ORACLE DRIFT — do not read the grid below ***")
        for line in drift:
            print(f"    {line}")
    else:
        print(f"oracle OK: all {len(ORACLE_RECORDED)} recorded columns "
              f"reproduce within {ORACLE_TOL}.")
    print()

    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for label, dir_name in RUN_SELECTION.items():
        run_dir = RUNS_ROOT / dir_name
        if not (run_dir / "detection.npz").exists():
            missing.append(label)
            continue
        rows.append(observable_row(label, run_dir, abundance_ref, ked_ref))

    print("=== Axis A geometry grid (fixed R × birth law, N = 500, one seed) ===")
    print(format_table(rows) if rows else "(no finished cells)")
    if missing:
        print(f"\nnot yet run: {missing}")
    print(
        "\nReading rule: compare cells to EACH OTHER (and to a re-weighted "
        "mixture reconstruction), never directly to the committed acceptance — "
        "the landing is a mixture property (plan §3.1). Per-cell w1_solv is "
        "not a discriminator below |Δ| ~ 0.15 at this N and single seed; "
        "supp / trap / nbar_det / n1_solv / midHot / deepKE are."
    )
    if rows:
        thin = [r["label"] for r in rows if r["num_scored"] < 0.6 * 2 * 500]
        if thin:
            print(
                f"Resolution re-check (plan §3.6): cells with a heavily "
                f"depleted scored ensemble — {thin} — lost >40% of ions to the "
                f"trapped channel; their per-cell scatter exceeds the "
                f"pre-registered yardstick."
            )
        marginal = [r["label"] for r in rows if r["trap_marg"] > 0.0]
        if marginal:
            print(
                f"Marginal-class caveat (plan §3.5b item 6): cells {marginal} "
                f"exclude a non-empty `trap_marg` class — ions energetically "
                f"able to escape but still helium-coupled at handover. That "
                f"fraction is a MODELLING exclusion conditional on a cubic law "
                f"extrapolated ~3× beyond its 9/18 Å calibration band, not a "
                f"prediction of nature. Read it with the bracket diagnostic "
                f"(tier2atlas_retained_bracket.py)."
            )
    if SAVE_CSV_PATH is not None and rows:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
