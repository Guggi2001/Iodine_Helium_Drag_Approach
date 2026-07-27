"""Atlas G3 Step 3 — MD confirmation-ring scorer report (plan §3.5d).

Scores the finished `gen_tier2atlas_g3ring.py` cells (14 corrected-geometry
N = 500 runs) on the committed observable vector, prints each next to its
frozen twin row (Δ = MD − twin), and evaluates the pre-registered
predictions GR-P2..GR-P6.

Pure scorer — runs nothing, mutates nothing. Physics conventions live in
`i2_helium_md.postprocess.tier2_confirmation`; the per-run observable row is
imported from `tier2atlas_geometry_table.py` (rule 1), which also carries the
pooled-battery scorer-drift oracle this report re-runs first (§1.4).

MD-level acceptance (§3.5d, bias-free): n1_solv in [0.19, 0.30] AND
nbar_det in [3.77, 4.37]. W1 reported, not gating; KE direction-only;
supp soft.

Invocation:

    python scripts/post_processing/tier2atlas_g3ring_table.py
"""

from __future__ import annotations

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
from i2_helium_md.simulation.checkpoint import load_neutral_checkpoint  # noqa: E402
from scripts.gen_tier2atlas_g3ring import (  # noqa: E402
    RING_MATRIX,
    TWIN_ROWS,
    ring_run_dir_name,
    verify_twin_preregistration,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUNS_ROOT,
    check_oracle,
    observable_row,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

# MD-level acceptance (plan §3.5d, pre-registered).
GATE_N1 = (0.19, 0.30)
GATE_NBAR = (3.77, 4.37)

# GR-P3/P4 bands (plan §3.5d).
NBAR_BIAS_BRACKET = (0.2, 3.2)
N1_TRANSFER_TOL = 0.05
N1_MIN_SOLVATED = 100

SAVE_CSV_PATH = None  # e.g. RUNS_ROOT / "atlas_g3ring.csv"

# ---------------------------------------------------------------------------


def sampled_geometry_columns(run_dir: Path) -> dict[str, Any]:
    """Measured sampled-geometry columns (the ring has per-molecule sizes)."""
    ckpt = load_neutral_checkpoint(run_dir / "neutral.npz")
    radii = np.asarray(ckpt.droplet_radii, dtype=float)
    r0 = np.asarray(ckpt.r0, dtype=float)
    return {
        "R_q50": float(np.quantile(radii, 0.50)),
        "depth_mean_A": float((radii - r0).mean()),
    }


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # §1.4 oracles: the frozen twin rows, then the scorer-drift oracle.
    verify_twin_preregistration()
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** SCORER-DRIFT ORACLE FAILED — do not read the ring ***")
        for line in drift:
            print(f"    {line}")
        return
    print("scorer-drift oracle OK (pooled battery reproduces).")
    print()

    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for cell in RING_MATRIX:
        run_dir = RUNS_ROOT / ring_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            missing.append(cell.label)
            continue
        row = observable_row(cell.label, run_dir, abundance_ref, ked_ref,
                             with_geometry=False)
        row.update(sampled_geometry_columns(run_dir))
        twin_trap, twin_supp, twin_nbar, twin_n1, twin_w1, twin_gate = (
            TWIN_ROWS[cell.label]
        )
        row.update({
            "v_c": cell.v_c, "well": cell.eb_tag, "tau": cell.tau_ps,
            "E0": cell.e0_eV,
            "twin_nbar": float(twin_nbar), "twin_n1": float(twin_n1),
            "twin_trap": float(twin_trap), "twin_gate": int(twin_gate),
            "d_nbar": float(row["nbar_det"]) - float(twin_nbar),
            "d_n1": float(row["n1_solv"]) - float(twin_n1),
            "d_trap": float(row["trap"]) - float(twin_trap),
            "md_gate": int(
                GATE_N1[0] <= float(row["n1_solv"]) <= GATE_N1[1]
                and GATE_NBAR[0] <= float(row["nbar_det"]) <= GATE_NBAR[1]
            ),
        })
        rows.append(row)

    print("=== G3 Step 3 MD ring (corrected geometry, N = 500, seed 20260728) ===")
    print(format_table(rows) if rows else "(no finished cells)")
    if missing:
        print(f"\nnot yet run: {missing}")
    if not rows:
        return
    by = {r["label"]: r for r in rows}

    # ---- Pre-registered predictions (GR-P2..P6) --------------------------
    print("\n=== Pre-registered predictions (plan §3.5d) ===")
    ab = [r for r in rows if r["label"].startswith(("a0", "b0"))]
    landed = [r["label"] for r in ab if r["md_gate"]]
    print(f"GR-P2 (basin transfer): arms A/B cells landing the MD acceptance:"
          f" {landed or 'NONE'} -> "
          f"{'CONFIRMED' if landed else 'REFUTED (basin fails at MD level)'}")

    paired = [r for r in rows]
    hot = [r for r in paired if r["d_nbar"] < 0.0]
    in_bracket = [r for r in paired
                  if NBAR_BIAS_BRACKET[0] <= -r["d_nbar"] <= NBAR_BIAS_BRACKET[1]]
    print(f"GR-P3 (twin-hot n̄ bias): MD < twin at {len(hot)}/{len(paired)} "
          f"cells; |Δ| inside [0.2, 3.2] at {len(in_bracket)}/{len(paired)} "
          f"(need >= 70%): "
          f"{'CONFIRMED' if len(hot) == len(paired) and 10 * len(in_bracket) >= 7 * len(paired) else 'CHECK'}")

    p4_cells = [r for r in rows if not r["label"].startswith(("c", "f"))]
    p4_bad = [r["label"] for r in p4_cells
              if abs(r["d_n1"]) > N1_TRANSFER_TOL]
    print(f"GR-P4 (n1 transfer <= {N1_TRANSFER_TOL}): violations "
          f"{p4_bad or 'none'} (cells under {N1_MIN_SOLVATED} solvated ions "
          f"read with judgment)")

    p5_msgs = []
    if "c50" in by:
        p5_msgs.append(
            f"c50 nbar={by['c50']['nbar_det']:.2f} "
            f"(predicted < {GATE_NBAR[0]}) gate={by['c50']['md_gate']}")
    if "c65" in by:
        p5_msgs.append(
            f"c65 n1={by['c65']['n1_solv']:.3f} (predicted < 0.19) "
            f"gate={by['c65']['md_gate']}")
    if "f725" in by:
        f = by["f725"]
        ok = (f["n1_solv"] < 0.05 and f["nbar_det"] > 10 and f["trap"] >= 0.17)
        p5_msgs.append(
            f"f725 n1={f['n1_solv']:.3f} nbar={f['nbar_det']:.2f} "
            f"trap={f['trap']:.3f} broken-landing reproduction: "
            f"{'CONFIRMED' if ok else 'CHECK'}")
    print("GR-P5 (controls): " + "; ".join(p5_msgs))

    if {"e0482", "e154"} <= set(by):
        # trap is chord-level (E0-independent), so any arm-A cell serves as
        # the eb1168 comparator on the (v_c 5.5, tau 4.8) chord.
        mid = by.get("a035") or by.get("a037")
        t0482, t154 = by["e0482"]["trap"], by["e154"]["trap"]
        t1168 = mid["trap"] if mid else float("nan")
        order_ok = mid is not None and t0482 <= t1168 <= t154
        print(f"GR-P6 (well trap ordering): eb0482 {t0482:.3f} <= eb1168 "
              f"{t1168:.3f} <= eb154 {t154:.3f}: "
              f"{'CONFIRMED' if order_ok else 'REFUTED/incomplete'}; "
              f"bound/marginal decomposition in the table above "
              f"(G4 retained-policy evidence).")

    if SAVE_CSV_PATH is not None and rows:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
