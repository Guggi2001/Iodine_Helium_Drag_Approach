"""Atlas §3.5i — the s(n) drag-state-coupling probe scorer (design §8, pre-registered).

Scores the finished `gen_tier2atlas_sn_probe.py` cells (sa22/sa30/sa44,
N = 1000, seed 20260729) against the committed h405 finals row (the
CRN-paired s ≡ 1 baseline) and evaluates the registered predictions
SC-P1..SC-P5 (`TIER2_DRAG_STATE_COUPLING_DESIGN.md` §8).

Pure scorer — runs nothing, mutates nothing. Physics conventions imported
from the committed modules (rule 1); the per-bin KE SD and the n = 3 bin
(needed by SC-P2/SC-P3, not part of the committed column set) are computed
locally on exactly the `lowke_columns` selection convention.

Order of operations (§1.4 — oracles first, hard-fail before any new number):

1. **O1 scorer-drift** — the standing pooled battery reproduces its recorded
   observable row (`check_oracle`).
2. **O2 committed-row** — the committed `g4fh405` run rescored through this
   script's code path equals the committed `atlas_g4finals_table.csv` h405
   row to 4 decimals.

SC-P1's "3× the seed scatter" floor is realised as 3× the CRN baseline's
per-bin standard error (KE1_sd/√KE1_n — the only scatter scale measurable
from the paired rows themselves; the h405 point has no committed multi-seed
battery). The registered magnitudes (KE₁ bands, the 0.08 eV needle-break
floor, ±0.10 deepKE) are hard-coded from design §8.

Invocation:

    python scripts/post_processing/tier2atlas_sn_table.py
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
from scripts.gen_tier2atlas_sn_probe import (  # noqa: E402
    RING,
    sn_run_dir_name,
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
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS — the registered SC magnitudes (design §8; do not tune)
# ---------------------------------------------------------------------------

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_sn_probe.csv"

SC1_SA30_BAND = (0.80, 1.10)   # SC-P1: the prior cell's KE1 landing band
SC2_NEEDLE_SD_FLOOR = 0.08     # SC-P2: per-bin n=1 KE SD floor [eV] (baseline 0.033)
SC3_DEEPKE_TOL = 0.10          # SC-P3: |deepKE - baseline| band (registered baseline 0.603)
SC4_KE1_REACH = 0.95           # SC-P4: kill evaluates on cells reaching this KE1
SC4_MIDHOT_CEILING = 1.15      # SC-P4: the PT-P3-form midHot ceiling
SUCCESS_KE1 = (0.90, 1.15)     # success shape (design §8)

# The CRN-paired baseline row (committed finals h405, s == 1).
H405_BASELINE = {"label": "h405", "rho_shell": float("nan")}

# ---------------------------------------------------------------------------


def sn_ke_columns(read) -> dict[str, Any]:
    """Per-bin KE SD (n = 1..3) + the n = 3 mean/count, `lowke_columns` convention.

    ``KE{n}_sd`` is the ddof=1 sample SD of the scored fragment KE in the
    exact-n bin (the SC-P2 needle observable); bins under 5 members report
    NaN (the `lowke_columns` floor convention).
    """
    out: dict[str, Any] = {}
    n_arr = np.asarray(read.n_scored)
    ke_arr = np.asarray(read.ke_scored_eV)
    for n in (1, 2, 3):
        sel = ke_arr[n_arr == n]
        if n == 3:
            out["KE3_n"] = int(sel.size)
            out["KE3_mean"] = float(sel.mean()) if sel.size >= 5 else float("nan")
        out[f"KE{n}_sd"] = (
            float(sel.std(ddof=1)) if sel.size >= 5 else float("nan")
        )
    return out


def scored_row(label: str, rho_shell: float, run_dir: Path,
               abundance_ref, ked_ref) -> dict[str, Any]:
    """One probe row on the §3.5g/h column conventions + the SC columns."""
    read = load_confirmation_run(run_dir, label=label)
    obs = observable_columns(read, abundance_ref, ked_ref)
    ke = lowke_columns(read)
    row: dict[str, Any] = {"label": label, "rho_shell": rho_shell}
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
    baseline = scored_row("h405", float("nan"),
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
        run_dir = RUNS_ROOT / sn_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            missing.append(cell.label)
            continue
        rows.append(scored_row(cell.label, cell.rho_shell_per_A3, run_dir,
                               abundance_ref, ked_ref))
    if missing:
        print(f"not yet run: {missing} — probe incomplete, stopping before "
              "any verdict.")
        return

    print("=== §3.5i s(n) probe (h405 pins, seed 20260729, CRN-paired "
          "baseline first) ===")
    print(format_table(rows))

    probe_rows = rows[1:]                       # sa22, sa30, sa44 (rho order)
    by_label = {r["label"]: r for r in probe_rows}
    base_ke1 = float(baseline["KE1_mean"])
    base_sd = float(baseline["KE1_sd"])
    base_n1ct = int(baseline["KE1_n"])
    scatter = base_sd / np.sqrt(base_n1ct) if base_n1ct >= 5 else float("nan")

    print("\n=== Registered SC verdicts (design §8) ===")
    # SC-P1: strict coupling ordering sa44 < sa30 < sa22; every cell above
    # baseline by >= 3x the seed scatter; sa30 in the registered band.
    ke1 = {lb: float(by_label[lb]["KE1_mean"]) for lb in by_label}
    ordered = ke1["sa44"] < ke1["sa30"] < ke1["sa22"]
    above = {lb: (v - base_ke1) >= 3.0 * scatter for lb, v in ke1.items()}
    sa30_in = SC1_SA30_BAND[0] <= ke1["sa30"] <= SC1_SA30_BAND[1]
    p1 = ordered and all(above.values()) and sa30_in
    print(f"SC-P1: KE1 baseline {base_ke1:.3f} -> "
          + ", ".join(f"{lb} {ke1[lb]:.3f}" for lb in ("sa22", "sa30", "sa44"))
          + f"; ordered sa44<sa30<sa22 -> {ordered}; all above baseline by "
          f">= 3x scatter ({3.0 * scatter:.3f}) -> {all(above.values())}; "
          f"sa30 in {SC1_SA30_BAND} -> {sa30_in} "
          f"=> {'CONFIRMED' if p1 else 'REFUTED'}")

    # SC-P2: the needle breaks — per-bin n=1 KE SD >= 0.08 at every cell.
    sds = {lb: float(by_label[lb]["KE1_sd"]) for lb in by_label}
    p2 = all(v >= SC2_NEEDLE_SD_FLOOR for v in sds.values())
    print(f"SC-P2 (SIGNATURE): n=1 KE SD baseline {base_sd:.3f} -> "
          + ", ".join(f"{lb} {sds[lb]:.3f}" for lb in ("sa22", "sa30", "sa44"))
          + f"; floor {SC2_NEEDLE_SD_FLOOR} "
          f"=> {'CONFIRMED (the needle breaks)' if p2 else 'REFUTED (kills the central claim regardless of means)'}")

    # SC-P3: grading dKE1 > dKE2 > dKE3 per cell; deepKE within +-0.10 of
    # the CRN baseline.
    base_deep = float(baseline["deepKE"])
    grading = {}
    for lb, r in by_label.items():
        d1 = float(r["KE1_mean"]) - base_ke1
        d2 = float(r["KE2_mean"]) - float(baseline["KE2_mean"])
        d3 = float(r["KE3_mean"]) - float(baseline["KE3_mean"])
        grading[lb] = (d1, d2, d3, d1 > d2 > d3)
    deep_ok = {lb: abs(float(r["deepKE"]) - base_deep) <= SC3_DEEPKE_TOL
               for lb, r in by_label.items()}
    p3 = all(g[3] for g in grading.values()) and all(deep_ok.values())
    for lb in ("sa22", "sa30", "sa44"):
        d1, d2, d3, ok = grading[lb]
        print(f"SC-P3 [{lb}]: dKE1 {d1:+.3f} > dKE2 {d2:+.3f} > dKE3 "
              f"{d3:+.3f} -> {ok}; deepKE {float(by_label[lb]['deepKE']):.3f} "
              f"(baseline {base_deep:.3f} +- {SC3_DEEPKE_TOL}) -> {deep_ok[lb]}")
    print(f"SC-P3 => {'CONFIRMED' if p3 else 'REFUTED'}")

    # SC-P4: kill criterion (PT-P3 form).
    reaching = [r for r in probe_rows if float(r["KE1_mean"]) >= SC4_KE1_REACH]
    if not reaching:
        print(f"SC-P4: no cell reaches KE1 >= {SC4_KE1_REACH} — kill "
              "criterion not evaluable on this probe.")
    else:
        killed = all(float(r["midHot"]) > SC4_MIDHOT_CEILING for r in reaching)
        print(f"SC-P4 (kill): cells with KE1 >= {SC4_KE1_REACH}: "
              + ", ".join(f"{r['label']} midHot {float(r['midHot']):.3f}"
                          for r in reaching)
              + (" => KILLED (every reaching cell exceeds midHot "
                 f"{SC4_MIDHOT_CEILING}; the axis stops — user adjudication)"
                 if killed else " => NOT killed"))

    # SC-P5: exploratory reads — measured, not predicted (the PT-P4 lesson).
    print("SC-P5 (exploratory, no registered signs): "
          + " | ".join(
              f"{lb}: nbar {float(by_label[lb]['nbar_det']):.2f}, "
              f"n1 {float(by_label[lb]['n1_solv']):.3f}, "
              f"supp {float(by_label[lb]['supp']):.3f}, "
              f"trap {float(by_label[lb]['trap']):.3f}, "
              f"W1 {float(by_label[lb]['w1_solv']):.3f}, "
              f"chi2 {float(by_label[lb]['chi2_med']):.0f}"
              for lb in ("sa22", "sa30", "sa44"))
          + f" (baseline: nbar {float(baseline['nbar_det']):.2f}, "
          f"n1 {float(baseline['n1_solv']):.3f}, "
          f"supp {float(baseline['supp']):.3f}, "
          f"trap {float(baseline['trap']):.3f}, "
          f"W1 {float(baseline['w1_solv']):.3f})")

    # Success shape (design §8): KE1 in band, midHot <= ceiling, gate intact.
    winners = [r["label"] for r in probe_rows
               if SUCCESS_KE1[0] <= float(r["KE1_mean"]) <= SUCCESS_KE1[1]
               and float(r["midHot"]) <= SC4_MIDHOT_CEILING
               and int(r["md_gate"]) == 1]
    near = [r["label"] for r in probe_rows
            if SUCCESS_KE1[0] <= float(r["KE1_mean"]) <= SUCCESS_KE1[1]
            and float(r["midHot"]) <= SC4_MIDHOT_CEILING
            and int(r["md_gate"]) == 0]
    print(f"\nSuccess-shape cells (KE1 in {SUCCESS_KE1}, midHot <= "
          f"{SC4_MIDHOT_CEILING}, gated): {winners if winners else 'NONE'}"
          + (f"; gate-broken but recenter-candidates: {near}" if near else "")
          + (" -> pooled 5-seed battery candidate (GV precedent); nothing "
             "adopts at probe level." if winners or near else ""))

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
