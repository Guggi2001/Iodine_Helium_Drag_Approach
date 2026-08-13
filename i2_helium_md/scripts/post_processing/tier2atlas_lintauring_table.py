"""Short-τ decisive pair scorer (plan §6.8 frozen reads TR-P1..TR-P3).

Scores the two `gen_tier2atlas_lintauring.py` cells against the committed
`h405p` bar, with the §6.6 context rows joined for continuity. Every cell
here shares seed 20260731 with the whole 14-cell linear family, so the
comparisons are CRN-paired.

Frozen reads (plan §6.8, pre-registered before launch):

* **TR-P1** — `t1` beats h405 on **both** W₁ (< 0.7675) and KE₁
  (> 0.6371). The user's actual question, at the best available cell.
* **TR-P2** — `t2` lands W₁ < 0.7675 while holding midHot ≤ 1.15. A
  histogram-only claim; its KE₁ is *expected* below h405 and TR-P2 is
  not failed by that.
* **TR-P3** — predictor adjudication, lands either way: |pred − MD| on
  W₁ against the LOO RMSE 0.081. Both cells missing by ≫ RMSE refutes
  the (n₁, tail) predictor and the τ extrapolation with it.
* **KILL (user pre-committed):** both cells at W₁ ≥ 0.7675 ⇒ the
  free-form linear family is closed for real.
* **midHot: REPORTED, never vetoing** (user, 2026-08-13).

Pure scorer — runs nothing, mutates nothing. Full vector to
`atlas_lintauring_table.csv`.

Invocation::

    python scripts/post_processing/tier2atlas_lintauring_table.py
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
from scripts.gen_tier2atlas_lintauring import (  # noqa: E402
    PREDICTOR_LOO_RMSE,
    TAU_MATRIX,
    TAU_PREDICTED,
    tau_run_dir_name,
    verify_tau_preregistration,
)
from scripts.gen_tier2atlas_linring import ring_run_dir_name  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUNS_ROOT,
    check_oracle,
    observable_row,
)
from scripts.post_processing.tier2atlas_linclone_table import (  # noqa: E402
    band_verdict,
)
from scripts.post_processing.tier2atlas_linjoint_table import (  # noqa: E402
    SEED_SD,
    anchor_oracle,
    score_run,
)
from scripts.post_processing.tier2atlas_linring_table import (  # noqa: E402
    GATE_N1,
    GATE_NBAR,
    ke_path_oracle,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

BAR_W1 = 0.7675          # committed h405p
BAR_KE1 = 0.6371
MIDHOT_BAND_TOP = 1.15
SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_lintauring_table.csv"


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # ---- oracles before any number ---------------------------------------
    verify_tau_preregistration()                     # TR-P0 + LR-P1
    drift = check_oracle(observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False))
    if drift:
        print("*** SCORER-DRIFT ORACLE FAILED — do not read ***")
        for line in drift:
            print(f"    {line}")
        return
    ke_drift = ke_path_oracle(abundance_ref, ked_ref)
    if ke_drift:
        print("*** KE-PATH ORACLE FAILED — do not read ***")
        for line in ke_drift:
            print(f"    {line}")
        return
    partner = score_run("h405p", RUNS_ROOT / ring_run_dir_name("h405p"),
                        abundance_ref, ked_ref)
    problems = anchor_oracle(partner, "h405p")
    if problems:
        print("*** PARTNER-ANCHOR ORACLE FAILED — do not read ***")
        for line in problems:
            print(f"    {line}")
        return
    print("oracles OK (TR-P0 + LR-P1 + drift + KE-path + committed h405p).")
    print()

    rows: list[dict[str, Any]] = []
    for cell in TAU_MATRIX:
        run_dir = RUNS_ROOT / tau_run_dir_name(cell.label)
        if not (run_dir / "detection.npz").exists():
            print(f"not yet run: {cell.label}")
            continue
        row = score_run(cell.label, run_dir, abundance_ref, ked_ref)
        row.update({"a": cell.a, "tau": cell.tau_ps, "E0": cell.e0_eV,
                    "role": cell.role})
        p = TAU_PREDICTED[cell.label]
        row.update({f"pred_{k}": v for k, v in p.items()})
        row["w1_err"] = float(row["w1_solv"]) - p["w1"]
        row["d_w1_h405"] = float(row["w1_solv"]) - BAR_W1
        row["d_KE1_h405"] = float(row["KE1_mean"]) - BAR_KE1
        rows.append(row)
    if not rows:
        return

    partner.update({"a": "-", "tau": 4.4, "E0": 0.405, "role": "bar"})
    table = rows + [partner]
    keys = set().union(*(r.keys() for r in table))
    for r in table:
        for k in keys:
            r.setdefault(k, "-")
    print("=== §6.8 short-τ decisive pair (N = 500, seed 20260731, CRN) ===")
    print(format_table(table))
    print()

    by = {r["label"]: r for r in rows}
    print("=== Verdicts (plan §6.8, pre-registered before launch) ===")

    t1 = by.get("t1")
    if t1 is not None:
        w, k = float(t1["w1_solv"]), float(t1["KE1_mean"])
        ok = w < BAR_W1 and k > BAR_KE1
        print(f"TR-P1 (t1 beats h405 on BOTH): W1 {w:.4f} < {BAR_W1} -> "
              f"{w < BAR_W1}; KE1 {k:.4f} > {BAR_KE1} -> {k > BAR_KE1}"
              f"  => {'*** PASS ***' if ok else 'FAIL'}")
        print(f"   midHot {float(t1['midHot']):.3f} (REPORTED, never vetoing), "
              f"tail {float(t1['tail_ge10']):.4f}, n1 {float(t1['n1_solv']):.4f}, "
              f"nbar {float(t1['nbar_det']):.3f}")

    t2 = by.get("t2")
    if t2 is not None:
        w, mh = float(t2["w1_solv"]), float(t2["midHot"])
        ok = w < BAR_W1 and mh <= MIDHOT_BAND_TOP
        print(f"TR-P2 (t2 histogram-only, in band): W1 {w:.4f} < {BAR_W1} -> "
              f"{w < BAR_W1}; midHot {mh:.3f} <= {MIDHOT_BAND_TOP} -> "
              f"{mh <= MIDHOT_BAND_TOP}  => {'*** PASS ***' if ok else 'FAIL'}")
        print(f"   KE1 {float(t2['KE1_mean']):.4f} (expected below h405; does "
              f"NOT fail TR-P2), tail {float(t2['tail_ge10']):.4f}")

    print(f"TR-P3 (predictor adjudication, LOO RMSE {PREDICTOR_LOO_RMSE}):")
    errs = []
    for r in rows:
        e = float(r["w1_err"])
        errs.append(abs(e))
        print(f"   {r['label']}: predicted W1 {r['pred_w1']:.3f}, measured "
              f"{float(r['w1_solv']):.4f}, err {e:+.4f} "
              f"({abs(e) / PREDICTOR_LOO_RMSE:.1f}x RMSE)")
    if errs and min(errs) > 2 * PREDICTOR_LOO_RMSE:
        print("   => the (n1, tail) predictor is REFUTED at both cells; the "
              "tau extrapolation fails and the §6.7 candidate list retires "
              "wholesale.")
    elif errs and max(errs) <= 2 * PREDICTOR_LOO_RMSE:
        print("   => the predictor holds within 2x its LOO RMSE at every "
              "cell — the licensed-leg route is validated OUT of its "
              "training range in tau.")
    else:
        print("   => mixed: one cell within 2x RMSE, one outside — reported, "
              "no wholesale claim either way.")

    if len(rows) == 2 and all(float(r["w1_solv"]) >= BAR_W1 for r in rows):
        print("KILL FIRES (user pre-committed): both cells at W1 >= "
              f"{BAR_W1} — the free-form linear family is CLOSED for real "
              "and no further scan is warranted.")
    else:
        print("KILL does NOT fire — at least one cell lands below the h405 "
              "W1 bar.")

    print("Gate (REPORTED, never selecting; three-way at the measured SDs):")
    for r in rows:
        print(f"   {r['label']}: n1 {float(r['n1_solv']):.4f} -> "
              f"{band_verdict(float(r['n1_solv']), GATE_N1, SEED_SD['n1_solv'])}"
              f"; nbar {float(r['nbar_det']):.3f} -> "
              f"{band_verdict(float(r['nbar_det']), GATE_NBAR, SEED_SD['nbar_det'])}"
              + ("  [policy-blocked]" if r.get("policy_blocked") == 1 else ""))

    print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, table)}")


if __name__ == "__main__":
    main()
