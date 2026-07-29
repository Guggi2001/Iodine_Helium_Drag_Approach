"""Atlas §3.5k — the budget-slope probe scorer (pre-registered, plan §3.5k).

Scores the finished `gen_tier2atlas_budget_probe.py` cells
(bud226k/bud411k/bud411, N = 1000, seed 20260729) against the committed
h405 finals row (the CRN-paired scale-1 baseline) and evaluates the
registered predictions BP-P1..BP-P4 and the BP-KILL slope criterion.

Pure scorer — runs nothing, mutates nothing. Physics conventions imported
from the committed modules (rule 1); the per-bin KE SD columns reuse the
§3.5i scorer's `sn_ke_columns` (the lowke selection convention).

Order of operations (§1.4 — oracles first, hard-fail before any new number):

1. **O1 scorer-drift** — the standing pooled battery reproduces its recorded
   observable row (`check_oracle`).
2. **O2 committed-row** — the committed `g4fh405` run rescored through this
   script's code path equals the committed finals h405 row to 4 decimals.

The kinematic slope S_k is a least-squares fit of KE1 against the
*kinematic* per-fragment budget (2.26 / 2.700616 / 4.11 — the baseline sits
at scale 1, i.e. 14.39964548/2.666/2) over the E0-pinned cells plus the
baseline, using only cells whose n = 1 bin meets the lowke >= 5-member
floor. bud411 (the wiring cell) is excluded from the fit by registration —
its E0 moves too.

Invocation:

    python scripts/post_processing/tier2atlas_budget_table.py
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
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    load_confirmation_run,
)
from scripts.gen_tier2atlas_g4finals import finals_run_dir_name  # noqa: E402
from scripts.gen_tier2atlas_budget_probe import (  # noqa: E402
    E_FRAG_REF_EV,
    RING,
    budget_run_dir_name,
)
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
    ORACLE_ROW_COLS,
    committed_h405_row,
    lowke_columns,
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
from scripts.post_processing.tier2atlas_ke_lown_scan import (  # noqa: E402
    GATE_N1,
    GATE_NBAR,
    KE1_ANCHOR_EV,
    KE2_REF_EV,
    lowke_counts,
)
from scripts.post_processing.tier2atlas_sn_table import sn_ke_columns  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS — the registered BP magnitudes (plan §3.5k; do not tune)
# ---------------------------------------------------------------------------

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_budget_probe.csv"

BP1_BANDS = {"bud226k": (0.29, 0.48), "bud411k": (1.13, 1.70)}
BP2_SLOPE_BAND = (0.35, 0.75)     # registered S_k band [eV per eV budget]
BP4_NEEDLE_SD_CEIL = 0.08         # k-cells: the needle must NOT break
BP_KILL_SLOPE = 0.2               # S_k below this => honest-residual branch

MIN_KE1_BIN = 5                   # the lowke n=1 bin floor for slope points


def scored_row(label: str, budget_eV: float, run_dir: Path,
               abundance_ref, ked_ref) -> dict[str, Any]:
    """One probe row on the §3.5g/h column conventions + the needle columns."""
    read = load_confirmation_run(run_dir, label=label)
    obs = observable_columns(read, abundance_ref, ked_ref)
    ke = lowke_columns(read)
    row: dict[str, Any] = {"label": label, "budget_eV": budget_eV}
    row.update({
        "num_scored": obs["num_scored"],
        "trap": obs["trap"], "supp": obs["supp"],
        "nbar_det": obs["nbar_det"], "n1_solv": obs["n1_solv"],
        "w1_solv": obs["w1_solv"], "midHot": obs["midHot"],
        "deepKE": obs["deepKE"], "chi2_med": obs["chi2_med"],
    })
    row.update(ke)
    row.update(lowke_counts(read))
    row.update(sn_ke_columns(read))
    row["KE1_r100"] = float(ke["KE1_mean"]) / KE1_ANCHOR_EV
    row["KE2_r"] = float(ke["KE2_mean"]) / KE2_REF_EV
    row["md_gate"] = int(
        GATE_N1[0] <= float(obs["n1_solv"]) <= GATE_N1[1]
        and GATE_NBAR[0] <= float(obs["nbar_det"]) <= GATE_NBAR[1]
    )
    return row


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # ---- O1: scorer drift.
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** O1 SCORER-DRIFT ORACLE FAILED — do not read the probe ***")
        for line in drift:
            print(f"    {line}")
        return
    print("O1 scorer-drift oracle OK (pooled standing battery reproduces).")

    # ---- O2: committed finals h405 row through this code path.
    committed = committed_h405_row()
    baseline = scored_row("h405", E_FRAG_REF_EV,
                          RUNS_ROOT / finals_run_dir_name("h405"),
                          abundance_ref, ked_ref)
    bad = [
        f"{c}: rescored {float(baseline[c]):.4f} vs committed "
        f"{float(committed[c]):.4f}"
        for c in ORACLE_ROW_COLS
        if f"{float(baseline[c]):.4f}" != f"{float(committed[c]):.4f}"
    ]
    if bad:
        print("*** O2 COMMITTED-ROW ORACLE FAILED — do not read the probe ***")
        for line in bad:
            print(f"    {line}")
        return
    print("O2 committed-row oracle OK (finals h405 rescored to 4 decimals).\n")

    # ---- The probe (new numbers — read only after O1/O2).
    rows: list[dict[str, Any]] = [baseline]
    missing: list[str] = []
    for cell in RING:
        run_dir = RUNS_ROOT / budget_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            missing.append(cell.label)
            continue
        rows.append(scored_row(cell.label, cell.budget_eV, run_dir,
                               abundance_ref, ked_ref))
    if missing:
        print(f"not yet run: {missing} — probe incomplete, stopping before "
              "any verdict.")
        return

    print("=== §3.5k budget probe (h405 pins, seed 20260729, CRN-paired "
          "baseline first) ===")
    print(format_table(rows))

    by_label = {r["label"]: r for r in rows}
    base_ke1 = float(baseline["KE1_mean"])

    print("\n=== Registered BP verdicts (plan §3.5k) ===")
    # BP-P1: placement bands on the E0-pinned cells.
    p1_parts = []
    for lb, band in BP1_BANDS.items():
        ke1 = float(by_label[lb]["KE1_mean"])
        ok = band[0] <= ke1 <= band[1]
        p1_parts.append(f"{lb} {ke1:.3f} in {band} -> {ok}")
    rising = (float(by_label["bud226k"]["KE1_mean"]) < base_ke1
              < float(by_label["bud411k"]["KE1_mean"]))
    print(f"BP-P1: {'; '.join(p1_parts)}; "
          f"rising along bud226k < h405 < bud411k -> {rising}")

    # BP-P2 + BP-KILL: the kinematic slope over baseline + E0-pinned cells.
    pts = []
    for lb in ("bud226k", "h405", "bud411k"):
        r = by_label[lb]
        if int(r["KE1_n"]) >= MIN_KE1_BIN:
            pts.append((float(r["budget_eV"]), float(r["KE1_mean"])))
        else:
            print(f"BP-P2 note: {lb} n=1 bin has {int(r['KE1_n'])} < "
                  f"{MIN_KE1_BIN} members — excluded from the slope fit.")
    if len(pts) < 2:
        print("BP-P2: fewer than 2 valid points — the slope is not "
              "measurable on this probe (record and stop; user adjudication).")
        s_k = float("nan")
    else:
        xs, ys = zip(*pts)
        s_k = float(np.polyfit(xs, ys, 1)[0])
        in_band = BP2_SLOPE_BAND[0] <= s_k <= BP2_SLOPE_BAND[1]
        print(f"BP-P2: kinematic slope S_k = {s_k:.3f} eV/eV over "
              f"{[f'{x:.3f}' for x in xs]} -> KE1 {[f'{y:.3f}' for y in ys]}; "
              f"registered band {BP2_SLOPE_BAND} -> "
              f"{'CONFIRMED' if in_band else 'OUTSIDE BAND'}")
    if not np.isnan(s_k):
        if s_k < BP_KILL_SLOPE:
            print(f"BP-KILL: S_k = {s_k:.3f} < {BP_KILL_SLOPE} => the "
                  "source-side lever is DEAD in-model; per the registration "
                  "the HONEST-RESIDUAL branch is taken (no (C) design).")
        else:
            print(f"BP-KILL: S_k = {s_k:.3f} >= {BP_KILL_SLOPE} => not "
                  "killed; the (C) design discussion proceeds with S_k as "
                  "its measured authority.")

    # BP-P3: wiring contrast (direction only).
    supp_w = float(by_label["bud411"]["supp"])
    supp_k = float(by_label["bud411k"]["supp"])
    print(f"BP-P3: supp(bud411) {supp_w:.3f} < supp(bud411k) {supp_k:.3f} "
          f"-> {supp_w < supp_k} (E0 0.6165 vs 0.405 at the same kinematics; "
          "magnitudes exploratory)")

    # BP-P4: needle persistence at the k-cells.
    p4_parts = []
    for lb in ("bud226k", "bud411k"):
        sd = float(by_label[lb]["KE1_sd"])
        p4_parts.append(f"{lb} KE1_sd {sd:.3f} < {BP4_NEEDLE_SD_CEIL} -> "
                        f"{sd < BP4_NEEDLE_SD_CEIL}")
    print("BP-P4 (needle persists at single-valued budget): "
          + "; ".join(p4_parts))

    # Exploratory reads — measured, not predicted (uniform-budget caveat:
    # midHot/deepKE here are NOT the channel-weighted (C) forecast).
    print("Exploratory (uniform-budget, NOT the (C) forecast): "
          + " | ".join(
              f"{lb}: nbar {float(by_label[lb]['nbar_det']):.2f}, "
              f"n1 {float(by_label[lb]['n1_solv']):.3f}, "
              f"supp {float(by_label[lb]['supp']):.3f}, "
              f"trap {float(by_label[lb]['trap']):.3f}, "
              f"W1 {float(by_label[lb]['w1_solv']):.3f}, "
              f"midHot {float(by_label[lb]['midHot']):.3f}, "
              f"deepKE {float(by_label[lb]['deepKE']):.3f}"
              for lb in ("bud226k", "bud411k", "bud411")))

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
