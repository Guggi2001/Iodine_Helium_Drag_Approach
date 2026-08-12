"""Atlas §6.5 Step 3 scorer — the E_bind MD confirmation pair (capped arm).

Scores the CRN pair

* ``linrh405p``  — the committed free-form-ring h405 cell at the bundle
  well (E_bind = 0.11675778 eV); **not re-run**, read from disk,
* ``ebmdh405s``  — the same cell with the well moved to 0.0482 eV and
  nothing else (``gen_tier2atlas_ebindmd.py``; cfg diff is exactly two
  keys),

and reports the measured MD transfer function ``dKE1/dE_bind`` against
the committed §6.5 Step-2 twin forecast, plus the pre-registered MD-P1..P6
verdicts.

Oracle (gate-on-committed-artifacts, runs before any new number): the
partner run re-scored here must reproduce its committed
``atlas_linring_table.csv`` row on every scored observable.

Band conventions are the committed ones (``midHot`` = geometric mean of
sim/ref per-bin means over n = 2–8, ``deepKE`` = arithmetic mean over
n = 10–17); this scorer additionally reports both bands as **absolute
mean KE in eV**, since the slope question is an energy question and the
ratio columns hide it behind a cold reference denominator.

Output: ``data/runs/h2b_forward_model/atlas_ebind_md.csv``.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import (  # noqa: E402
    load_ihe_ked_reference,
)
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    load_confirmation_run,
)
from scripts.gen_tier2atlas_ebindmd import (  # noqa: E402
    EBIND_EV,
    LABEL,
    MDP1_SLOPE_BAND,
    MDP2_SLOPE_BAND,
    MDP3_KE2_MAX_GAP,
    MDP4_TRAP_MAX,
    MDP5_GRADING_BAND,
    MDP6_KE1_CEILING,
    PARTNER_RUN,
    TWIN_KE1_EB0482_CAPPED,
    TWIN_SLOPE_CAPPED,
    WELL_TAG,
    md_run_dir_name,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    RUNS_ROOT,
    observable_row,
)
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
    lowke_columns,
)
from scripts.post_processing.tier2atlas_linring_table import (  # noqa: E402
    ke_width_columns,
)

OUT_CSV = PROJECT_ROOT / "data" / "runs" / "h2b_forward_model" / \
    "atlas_ebind_md.csv"
RING_TABLE_CSV = PROJECT_ROOT / "data" / "runs" / "h2b_forward_model" / \
    "atlas_linring_table.csv"
TWIN_CSV = PROJECT_ROOT / "data" / "runs" / "h2b_forward_model" / \
    "atlas_ebind_twin.csv"

# The committed ring columns the partner must reproduce exactly.
ORACLE_COLS = ("trap", "supp", "nbar_det", "n1_solv", "w1_solv", "midHot",
               "deepKE", "KE1_mean", "KE2_mean")
ORACLE_TOL = 1e-9
MIN_BIN_COUNT = 20          # the committed g3_score / scorer bin rule


def band_mean_ke_eV(read, lo: int, hi: int, aggregate: str
                    ) -> tuple[float, int]:
    """Absolute band mean KE [eV] over n in ``[lo, hi]``.

    Same present-bin rule as the committed scorers (a bin counts only with
    >= ``MIN_BIN_COUNT`` scored fragments), so the bin count is directly
    comparable to ``midHot_bins`` / ``deepKE_bins``.
    """
    vals = []
    for n in range(lo, hi + 1):
        sel = read.ke_scored_eV[read.n_scored == n]
        if sel.size >= MIN_BIN_COUNT:
            vals.append(float(sel.mean()))
    if not vals:
        return float("nan"), 0
    if aggregate == "geometric":
        return float(np.exp(np.mean(np.log(vals)))), len(vals)
    return float(np.mean(vals)), len(vals)


def committed_ring_row(label: str) -> dict[str, str]:
    with open(RING_TABLE_CSV, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["label"] == label:
                return row
    raise AssertionError(f"no {label!r} row in {RING_TABLE_CSV.name}")


def committed_twin_row(arm: str, tag: str) -> dict[str, str]:
    with open(TWIN_CSV, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["arm"] == arm and row["E_bind_tag"] == tag:
                return row
    raise AssertionError(f"no {arm}/{tag} row in {TWIN_CSV.name}")


def score_cell(label: str, run_dir: Path, abundance_ref, ked_ref
               ) -> dict[str, Any]:
    read = load_confirmation_run(run_dir, label=label)
    row = observable_row(label, run_dir, abundance_ref, ked_ref,
                         with_geometry=False)
    row.update(lowke_columns(read))
    row.update(ke_width_columns(read))
    mid, mid_bins = band_mean_ke_eV(read, 2, 8, "geometric")
    deep, deep_bins = band_mean_ke_eV(read, 10, 17, "arithmetic")
    row.update({
        "KE1_n": int((read.n_scored == 1).sum()),
        "KE2_n": int((read.n_scored == 2).sum()),
        "ke_mid_geo_eV": mid, "ke_mid_bins": mid_bins,
        "ke_deep_eV": deep, "ke_deep_bins": deep_bins,
    })
    return row


def verify_partner_oracle(row: dict[str, Any]) -> None:
    """The partner's re-scored row vs its committed ring row."""
    committed = committed_ring_row("h405p")
    drift = [
        f"{col}: got {float(row[col]):.6f} vs committed {float(committed[col]):.6f}"
        for col in ORACLE_COLS
        if not abs(float(row[col]) - float(committed[col])) <= ORACLE_TOL
    ]
    if drift:
        raise AssertionError(
            "partner oracle FAILED (the committed ring row did not "
            "reproduce): " + "; ".join(drift))
    print(f"[ebindmd] partner oracle PASSED: {len(ORACLE_COLS)} committed "
          "atlas_linring_table.csv observables reproduced", flush=True)


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    deep_well = EBIND_EV["eb1168"]
    shallow = EBIND_EV[WELL_TAG]
    d_well = shallow - deep_well          # negative: shallower

    partner = score_cell("h405p", RUNS_ROOT / PARTNER_RUN,
                         abundance_ref, ked_ref)
    verify_partner_oracle(partner)
    cell = score_cell(LABEL, RUNS_ROOT / md_run_dir_name(),
                      abundance_ref, ked_ref)

    partner.update({"E_bind_tag": "eb1168", "E_bind_eV": deep_well})
    cell.update({"E_bind_tag": WELL_TAG, "E_bind_eV": shallow})

    def slope(col: str) -> float:
        return (float(cell[col]) - float(partner[col])) / d_well

    s_ke1, s_ke2 = slope("KE1_mean"), slope("KE2_mean")
    s_mid, s_deep = slope("ke_mid_geo_eV"), slope("ke_deep_eV")
    s_trap = slope("trap")
    grading = abs(s_mid) - abs(s_ke1)

    twin_deep = committed_twin_row("H", "eb1168")
    twin_shal = committed_twin_row("H", WELL_TAG)
    twin_ke1_deep = float(twin_deep["n1_ke_eV"])
    twin_ke1_shal = float(twin_shal["n1_ke_eV"])
    twin_grading = (abs((float(twin_shal["ke_mid_geo_eV"])
                         - float(twin_deep["ke_mid_geo_eV"])) / d_well)
                    - abs((twin_ke1_shal - twin_ke1_deep) / d_well))

    verdicts = {
        "MD_P1": "PASS" if MDP1_SLOPE_BAND[0] <= s_ke1 <= MDP1_SLOPE_BAND[1]
                 else "FAIL",
        "MD_P2": "PASS" if MDP2_SLOPE_BAND[0] <= s_ke1 <= MDP2_SLOPE_BAND[1]
                 else "FAIL",
        "MD_P3": "PASS" if abs(s_ke2 - s_ke1) <= MDP3_KE2_MAX_GAP else "FAIL",
        "MD_P4": "PASS" if float(cell["trap"]) <= MDP4_TRAP_MAX else "FAIL",
        "MD_P5": "PASS" if MDP5_GRADING_BAND[0] <= grading
                 <= MDP5_GRADING_BAND[1] else "FAIL",
        "MD_P6": "PASS" if float(cell["KE1_mean"]) < MDP6_KE1_CEILING
                 else "FAIL",
    }

    rows = []
    for r in (partner, cell):
        out = dict(r)
        out.update({
            "d_well_eV": round(d_well, 6),
            "slope_KE1": round(s_ke1, 4), "slope_KE2": round(s_ke2, 4),
            "slope_mid_geo": round(s_mid, 4), "slope_deep": round(s_deep, 4),
            "slope_trap_per_eV": round(s_trap, 4),
            "grading_mid_minus_n1": round(grading, 4),
            "twin_slope_KE1": TWIN_SLOPE_CAPPED,
            "twin_KE1_at_well": round(
                twin_ke1_shal if out["E_bind_tag"] == WELL_TAG
                else twin_ke1_deep, 4),
            "twin_grading": round(twin_grading, 4),
            "md_minus_twin_KE1": round(
                float(out["KE1_mean"])
                - (twin_ke1_shal if out["E_bind_tag"] == WELL_TAG
                   else twin_ke1_deep), 4),
        })
        out.update(verdicts)
        rows.append(out)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    print(f"\n=== §6.5 Step 3: the MD transfer function (capped arm) ===")
    print(f"  well {deep_well:.5f} -> {shallow:.4f} eV "
          f"(d_well {d_well:+.5f}), CRN seed shared, N = 500")
    print(f"  KE1  {float(partner['KE1_mean']):.4f} -> "
          f"{float(cell['KE1_mean']):.4f}   slope {s_ke1:+.4f} eV/eV")
    print(f"  KE2  {float(partner['KE2_mean']):.4f} -> "
          f"{float(cell['KE2_mean']):.4f}   slope {s_ke2:+.4f}")
    print(f"  mid  {float(partner['ke_mid_geo_eV']):.4f} -> "
          f"{float(cell['ke_mid_geo_eV']):.4f}   slope {s_mid:+.4f} "
          f"(grading {grading:+.4f}, twin {twin_grading:+.4f})")
    print(f"  deep {float(partner['ke_deep_eV']):.4f} -> "
          f"{float(cell['ke_deep_eV']):.4f}   slope {s_deep:+.4f}")
    print(f"  trap {float(partner['trap']):.4f} -> "
          f"{float(cell['trap']):.4f}   lever {s_trap:+.4f}/eV")
    print(f"  nbar {float(partner['nbar_det']):.3f} -> "
          f"{float(cell['nbar_det']):.3f} ; n1 "
          f"{float(partner['n1_solv']):.4f} -> {float(cell['n1_solv']):.4f}")
    print(f"  twin forecast: slope {TWIN_SLOPE_CAPPED:+.4f}, "
          f"KE1({WELL_TAG}) {TWIN_KE1_EB0482_CAPPED:.4f} vs MD "
          f"{float(cell['KE1_mean']):.4f}")
    print("  pre-registered verdicts: "
          + ", ".join(f"{k}={v}" for k, v in verdicts.items()))
    print(f"written -> {OUT_CSV}")


if __name__ == "__main__":
    main()
