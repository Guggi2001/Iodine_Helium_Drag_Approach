"""Atlas G4 Step 2 Block V — the h405 pooled-battery scorer (plan §3.5f).

Scores the five `gen_tier2atlas_g4step2_battery.py` members (h405 × fresh
seeds 20260730–34, N = 1000 each) plus their in-memory §4cc pool, and
evaluates the pre-registered GV verdicts.

Pure scorer — runs nothing, mutates nothing. Physics conventions live in
`i2_helium_md.postprocess.tier2_confirmation`; the observable vector is
`observable_columns` from `tier2atlas_geometry_table.py` (rule 1), the pool
is `pool_confirmation_reads` (no pooled container dir), and the joint score
is the finals `S` convention (licensed axes W₁ + midHot, provisional norms).

Order of operations (§1.4 — oracles first, hard-fail before any new number):

1. **Scorer-drift oracle** — the standing pooled battery reproduces its 7
   recorded columns within 0.002.
2. **Committed-row oracle** — the committed Block-3 `g4fh405` run rescored
   through this script's own code path equals the committed
   `atlas_g4finals_table.csv` h405 row to 4 decimals on
   nbar_det / n1_solv / w1_solv / midHot / deepKE.

Pre-registered (plan §3.5f):

* **GV-P1 (gate):** pooled n1_solv ∈ [0.19, 0.30] ∧ pooled nbar_det ∈
  [3.77, 4.37].
* **GV-P2 (floor consistency):** pooled w1_solv ∈ [0.64, 0.78]; per-seed
  SD printed next to the standing 0.095 (informational, not a hard clause).
* **GV-P3 (score):** pooled S < 1.683 (a037's re-measured N = 1000 S);
  distance to the Block-3 single-seed 1.559 printed.

GV-P1/P2 failure ⇒ h405 not seed-robust ⇒ no adjudication fires. All three
pass ⇒ the G4 adjudications become FIREABLE (user calls) — nothing adopts
here.

Invocation:

    python scripts/post_processing/tier2atlas_g4step2_battery_table.py
"""

from __future__ import annotations

import csv
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
    load_confirmation_run,
    pool_confirmation_reads,
)
from scripts.gen_tier2atlas_g4finals import finals_run_dir_name  # noqa: E402
from scripts.gen_tier2atlas_g4step2_battery import (  # noqa: E402
    MEMBER_SEEDS,
    member_run_dir_name,
    verify_h405_twin_row,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUNS_ROOT,
    check_oracle,
    observable_columns,
    observable_row,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)
from scripts.post_processing.tier2atlas_g4_transfer import (  # noqa: E402
    PROVISIONAL_NORM,
    joint_score,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_g4step2_battery.csv"

FINALS_TABLE_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_g4finals_table.csv"

# The finals S convention (Block-0-licensed axes, provisional norms).
LICENSED_AXES = ("w1", "midhot")
INCUMBENT = {"w1_solv": 0.5791, "midHot": 1.0108, "deepKE": 0.6327}
W1_SEED_SD_STANDING = 0.0954    # per-seed SD at N = 1000 (Block 0)
A037_N1000_S = 1.683            # GV-P3 bar (finals table, a037 at N = 1000)
H405_BLOCK3_S = 1.559           # proximity reference (finals table)

# Pre-registered GV bands (plan §3.5f).
GATE_N1 = (0.19, 0.30)
GATE_NBAR = (3.77, 4.37)
GV_P2_W1 = (0.64, 0.78)

# Committed-row oracle columns (4-decimal comparison).
ORACLE_ROW_COLS = ("nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE")

# The pre-registered low-n KE axis (plan §3.5f, user adjudication
# 2026-07-28): KE1 vs the reference n = 1 MEDIAN (I77 anchor), KE2 vs the
# reference n = 2 MEAN. Reported here per member + pooled; the hard-gate
# band is frozen from this battery's measured per-seed SD.
LOWKE_REF = {1: ("median", 1.128), 2: ("mean", 0.706)}

# ---------------------------------------------------------------------------


def _s(row: dict[str, Any]) -> float:
    """Joint score on the licensed axes, provisional norms (plan §3.5e/§3.5f)."""
    return joint_score(float(row["w1_solv"]), float(row["midHot"]),
                       float(row["deepKE"]), PROVISIONAL_NORM,
                       axes=LICENSED_AXES)


def lowke_columns(read) -> dict[str, Any]:
    """The pre-registered low-n KE observables (mean/median/mode per bin).

    Mode = midpoint of the tallest 40-bin histogram cell on [0, 3.5] eV
    (display diagnostic; the scored quantities are KE1_mean / KE2_mean).
    """
    out: dict[str, Any] = {}
    for n in (1, 2):
        sel = read.ke_scored_eV[read.n_scored == n]
        if sel.size < 5:
            out[f"KE{n}_mean"] = float("nan")
            out[f"KE{n}_med"] = float("nan")
            out[f"KE{n}_mode"] = float("nan")
            continue
        hist, edges = np.histogram(sel, bins=40, range=(0.0, 3.5))
        i = int(np.argmax(hist))
        out[f"KE{n}_mean"] = float(sel.mean())
        out[f"KE{n}_med"] = float(np.median(sel))
        out[f"KE{n}_mode"] = float(0.5 * (edges[i] + edges[i + 1]))
    return out


def committed_h405_row() -> dict[str, str]:
    """The committed finals-table h405 row (string fields)."""
    with open(FINALS_TABLE_CSV, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["label"] == "h405":
                return row
    raise AssertionError(f"no h405 row in {FINALS_TABLE_CSV}")


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # §1.4 oracles — all three before any new number is read.
    verify_h405_twin_row()
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** SCORER-DRIFT ORACLE FAILED — do not read the battery ***")
        for line in drift:
            print(f"    {line}")
        return
    print("scorer-drift oracle OK (pooled standing battery reproduces).")

    committed = committed_h405_row()
    rescored = observable_columns(
        load_confirmation_run(RUNS_ROOT / finals_run_dir_name("h405"),
                              label="g4fh405"),
        abundance_ref, ked_ref,
    )
    bad = [
        f"{c}: rescored {float(rescored[c]):.4f} vs committed "
        f"{float(committed[c]):.4f}"
        for c in ORACLE_ROW_COLS
        if f"{float(rescored[c]):.4f}" != f"{float(committed[c]):.4f}"
    ]
    if bad:
        print("*** COMMITTED-ROW ORACLE FAILED — do not read the battery ***")
        for line in bad:
            print(f"    {line}")
        return
    print("committed-row oracle OK (Block-3 h405 rescored to 4 decimals).\n")

    # Members + pool.
    rows: list[dict[str, Any]] = []
    reads = []
    missing: list[str] = []
    for member in MEMBER_SEEDS:
        run_dir = RUNS_ROOT / member_run_dir_name(member)
        if not (run_dir / "detection.npz").exists():
            missing.append(member)
            continue
        read = load_confirmation_run(run_dir, label=member)
        reads.append(read)
        row: dict[str, Any] = {"label": member, "seed": MEMBER_SEEDS[member]}
        row.update(observable_columns(read, abundance_ref, ked_ref))
        row.update(lowke_columns(read))
        row["md_gate"] = int(
            GATE_N1[0] <= float(row["n1_solv"]) <= GATE_N1[1]
            and GATE_NBAR[0] <= float(row["nbar_det"]) <= GATE_NBAR[1]
        )
        row["S"] = _s(row)
        rows.append(row)
    if missing:
        print(f"not yet run: {missing} — battery incomplete, stopping before "
              "any pooled number.")
        return

    pooled_read = pool_confirmation_reads(reads, label="h405pooled")
    pooled: dict[str, Any] = {"label": "pooled", "seed": "-"}
    pooled.update(observable_columns(pooled_read, abundance_ref, ked_ref))
    pooled.update(lowke_columns(pooled_read))
    pooled["md_gate"] = int(
        GATE_N1[0] <= float(pooled["n1_solv"]) <= GATE_N1[1]
        and GATE_NBAR[0] <= float(pooled["nbar_det"]) <= GATE_NBAR[1]
    )
    pooled["S"] = _s(pooled)
    rows.append(pooled)

    print("=== G4 Step 2 Block V — h405 battery (5 × N = 1000, fresh seeds) ===")
    print(format_table(rows))

    w1s = np.asarray([float(r["w1_solv"]) for r in rows if r["label"] != "pooled"])
    print(f"\nper-seed W₁: mean {w1s.mean():.4f} ± SD {w1s.std(ddof=1):.4f} "
          f"(standing battery per-seed SD {W1_SEED_SD_STANDING:.4f})")

    print("\n=== Pre-registered GV verdicts (plan §3.5f) ===")
    p_n1 = GATE_N1[0] <= float(pooled["n1_solv"]) <= GATE_N1[1]
    p_nbar = GATE_NBAR[0] <= float(pooled["nbar_det"]) <= GATE_NBAR[1]
    gv1 = p_n1 and p_nbar
    print(f"GV-P1 (pooled gate): n1 {float(pooled['n1_solv']):.4f} in "
          f"{GATE_N1} -> {p_n1}; nbar {float(pooled['nbar_det']):.3f} in "
          f"{GATE_NBAR} -> {p_nbar} => "
          f"{'CONFIRMED' if gv1 else 'REFUTED'}")
    gv2 = GV_P2_W1[0] <= float(pooled["w1_solv"]) <= GV_P2_W1[1]
    print(f"GV-P2 (floor consistency): pooled W₁ "
          f"{float(pooled['w1_solv']):.4f} in {GV_P2_W1} => "
          f"{'CONFIRMED' if gv2 else 'REFUTED'} "
          f"(per-seed SD {w1s.std(ddof=1):.4f} vs standing "
          f"{W1_SEED_SD_STANDING:.4f}, informational)")
    gv3 = float(pooled["S"]) < A037_N1000_S
    print(f"GV-P3 (score): pooled S {float(pooled['S']):.3f} < a037's "
          f"{A037_N1000_S} => {'CONFIRMED' if gv3 else 'REFUTED'} "
          f"(Block-3 single-seed h405 S {H405_BLOCK3_S}; "
          f"incumbent finc1v725 S "
          f"{joint_score(INCUMBENT['w1_solv'], INCUMBENT['midHot'], INCUMBENT['deepKE'], PROVISIONAL_NORM, axes=LICENSED_AXES):.3f})")

    print("\n=== Pre-registered low-n KE axis (plan §3.5f; reported, "
          "norms frozen from this scatter) ===")
    for n in (1, 2):
        kind, ref_val = LOWKE_REF[n]
        vals = np.asarray([float(r[f"KE{n}_mean"]) for r in rows
                           if r["label"] != "pooled"])
        pv = float(pooled[f"KE{n}_mean"])
        print(f"KE{n}: pooled mean {pv:.3f} eV vs ref {kind} {ref_val:.3f} "
              f"(ratio {pv / ref_val:.3f}); per-seed mean "
              f"{vals.mean():.3f} ± SD {vals.std(ddof=1):.3f} "
              f"-> frozen gate band ratio 1 ± {2 * vals.std(ddof=1) / ref_val:.3f} "
              f"(2×SD/ref)")
        print(f"     pooled median {float(pooled[f'KE{n}_med']):.3f} / mode "
              f"{float(pooled[f'KE{n}_mode']):.3f} eV "
              f"(incumbent pooled n={n} mean: {1.034 if n == 1 else 0.754})")
    s_lowke = (abs(np.log(float(pooled["KE1_mean"]) / LOWKE_REF[1][1]))
               + abs(np.log(float(pooled["KE2_mean"]) / LOWKE_REF[2][1]))) \
        / np.log(1.15)
    print(f"low-n KE score terms (|ln ratio|/ln 1.15, additive to S in "
          f"future rankings): {s_lowke:.3f}")

    if gv1 and gv2 and gv3:
        print("\nAll GV verdicts pass: the G4 adjudications are FIREABLE "
              "(user calls) — nothing adopted here.")
    else:
        print("\nGV failure: h405 is not verified at battery level — no "
              "adjudication fires; back to the ridge with this evidence.")

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
