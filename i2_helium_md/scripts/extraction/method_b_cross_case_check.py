"""Method-B cross-case held-out check (METHOD_B §4 axis 2, plan B5).

The load-bearing held-out validation that the full-window calibration leaves
in scope before the VMI tier: a Method-B fit must look like the *same
transport physics* across the two cases, not like per-case curve bending.
Three checks, each with a NAMED, PROVISIONAL threshold (per the §6.10 rule,
numbers are set at first-runs and recorded here with their rationale):

(i)   **18 A A<->B coefficient agreement.** The 18 A case is clean-radial, so
      its Method-B fit should land in the same coefficient region as the
      Method-A force-balance fit; a large drift is the overfit signature.
      The 9 A A<->B check is DROPPED: the 9 A Method-A artifact is
      provenance-broken (hand-edited; recorded in drag_migration_log.md) per
      the user decision of 2026-06-11.
(ii)  **Residual-quality asymmetry.** The 9 A in-window RMSE is expected to
      exceed clean-radial 18 A (the model-dimensionality cost), but an
      *extreme* ratio means the transverse contamination dominated the fit.
(iii) **Cross-case E_bind consistency.** The effective binding stands in for
      the same I+/He dynamical-ejection physics in both cases; a large
      disagreement is a red flag independent of the trajectory fit.

Writes ``held_out_validation.json`` next to each case's Method-B bundle
(keeping ``fit_parameters.json`` immutable post-fit), including the
``"vmi": "pending"`` stamp -- the VMI tier remains the ultimate arbiter and
Tier 1 stays gated until it runs (plan B7).

Usage: run after method_b_extraction.py has produced both cases::

    python scripts/extraction/method_b_cross_case_check.py
"""

from __future__ import annotations

import datetime
import json
from pathlib import Path
import sys

# =============================================================================
# USER SETTINGS
# =============================================================================
# Provisional named thresholds (first-runs rule, §6.10): recorded with the
# verdict so a later re-derivation can re-judge the same numbers.
COEFF_RATIO_BAND = (0.5, 2.0)   # (i): a_B/a_A and b_B/b_A inside this band
RESIDUAL_RATIO_MAX = 3.0        # (ii): rmse_9A / rmse_18A at most this
E_BIND_MAX_DIFF_EV = 0.05       # (iii): |E_9A - E_18A| at most this [eV]

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import _git_branch  # noqa: E402
from i2_helium_md.presets import REFERENCE_DRAG_ROOT  # noqa: E402


def _load(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def main() -> int:
    b_18 = REFERENCE_DRAG_ROOT / "18A" / "trajectory_matching" / "fit_parameters.json"
    b_9 = REFERENCE_DRAG_ROOT / "9A" / "trajectory_matching" / "fit_parameters.json"
    a_18 = REFERENCE_DRAG_ROOT / "18A" / "linear_and_cubic" / "fit_parameters.json"

    if not b_18.is_file():
        print(f"missing 18A Method-B bundle: {b_18}; run the extraction first")
        return 1
    fit18 = _load(b_18)
    methodA18 = _load(a_18)
    fit9 = _load(b_9) if b_9.is_file() else None
    if fit9 is None:
        print(f"NOTE: 9A Method-B bundle absent ({b_9}); checks (ii)/(iii) "
              "and the 9A verdict are deferred.")

    print("===== Method-B cross-case held-out check =====")

    # --- (i) 18A A<->B agreement ---
    ratio_a = fit18["a"] / methodA18["a"]
    ratio_b = fit18["b"] / methodA18["b"]
    lo, hi = COEFF_RATIO_BAND
    check_i = (lo <= ratio_a <= hi) and (lo <= ratio_b <= hi)
    print(f"  (i)  18A A<->B: a_B/a_A = {ratio_a:.3f}, b_B/b_A = {ratio_b:.3f} "
          f"(band [{lo}, {hi}]) -> {'PASS' if check_i else 'FAIL'}")
    print("       (9A A<->B dropped: Method-A artifact provenance-broken, "
          "2026-06-11 decision)")

    # --- (ii) residual-quality asymmetry ---
    check_ii = None
    residual_ratio = None
    if fit9 is not None:
        residual_ratio = fit9["objective_rmse_Aps"] / fit18["objective_rmse_Aps"]
        check_ii = residual_ratio <= RESIDUAL_RATIO_MAX
        print(f"  (ii) residual asymmetry: rmse 9A {fit9['objective_rmse_Aps']:.4f}"
              f" / 18A {fit18['objective_rmse_Aps']:.4f} = {residual_ratio:.2f}"
              f" (max {RESIDUAL_RATIO_MAX}) -> {'PASS' if check_ii else 'FAIL'}")

    # --- (iii) cross-case E_bind consistency ---
    check_iii = None
    e_diff = None
    if fit9 is not None:
        e_diff = abs(fit9["effective_binding_energy_I_ion_eV"]
                     - fit18["effective_binding_energy_I_ion_eV"])
        check_iii = e_diff <= E_BIND_MAX_DIFF_EV
        print(f"  (iii) E_bind: 9A {fit9['effective_binding_energy_I_ion_eV']:.4f}"
              f" vs 18A {fit18['effective_binding_energy_I_ion_eV']:.4f} eV, "
              f"|diff| = {e_diff:.4f} (max {E_BIND_MAX_DIFF_EV}) -> "
              f"{'PASS' if check_iii else 'FAIL'}")

    # --- verdicts (one line per case, METHOD_B §8) ---
    verdict_18 = (
        "PASS: 18A B-fit agrees with its Method-A fit (clean-radial case); "
        "VMI held-out still pending."
        if check_i else
        "FAIL: 18A B-fit drifted out of the Method-A coefficient region -- "
        "overfit signature on the clean case; do not wire into presets."
    )
    if fit9 is None:
        verdict_9 = "DEFERRED: no 9A Method-B bundle yet."
    elif check_ii and check_iii:
        verdict_9 = (
            "PASS (cross-case axes): residual asymmetry and E_bind "
            "consistency inside provisional bands; the standing "
            "transverse-contamination flag remains, VMI held-out pending."
        )
    else:
        verdict_9 = (
            "FAIL: 9A fit violates the cross-case bands -- consistent with "
            "the coefficients absorbing the transverse discrepancy; record "
            "as not-yet-usable for the radial MD (METHOD_B §8), do not wire "
            "into presets."
        )
    print(f"  18A verdict: {verdict_18}")
    print(f"  9A  verdict: {verdict_9}")

    # --- write held_out_validation.json per case ---
    stamp = {
        "date": datetime.date.today().isoformat(),
        "branch": _git_branch(PROJECT_ROOT),
        "thresholds": {
            "coeff_ratio_band": list(COEFF_RATIO_BAND),
            "residual_ratio_max": RESIDUAL_RATIO_MAX,
            "e_bind_max_diff_eV": E_BIND_MAX_DIFF_EV,
            "provenance": "provisional first-runs thresholds (§6.10)",
        },
        "vmi": "pending",  # ultimate arbiter; Tier 1 stays gated until run
    }
    out18 = {
        **stamp,
        "checks": {
            "18A_A_vs_B_coefficient_agreement": {
                "ratio_a": ratio_a, "ratio_b": ratio_b, "pass": check_i,
            },
        },
        "verdict": verdict_18,
    }
    path18 = b_18.parent / "held_out_validation.json"
    with open(path18, "w", encoding="utf-8") as fh:
        json.dump(out18, fh, indent=2)
    print(f"  wrote {path18}")

    if fit9 is not None:
        out9 = {
            **stamp,
            "checks": {
                "9A_A_vs_B_coefficient_agreement": (
                    "dropped: Method-A artifact provenance-broken "
                    "(hand-edited); see drag_migration_log.md 2026-06-11"
                ),
                "residual_quality_asymmetry": {
                    "rmse_ratio_9A_over_18A": residual_ratio, "pass": check_ii,
                },
                "cross_case_e_bind_consistency": {
                    "abs_diff_eV": e_diff, "pass": check_iii,
                },
            },
            "verdict": verdict_9,
        }
        path9 = b_9.parent / "held_out_validation.json"
        with open(path9, "w", encoding="utf-8") as fh:
            json.dump(out9, fh, indent=2)
        print(f"  wrote {path9}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
