"""Atlas D4 Step-1 oracle: verify the Method-B shared-fit machinery.

Sensitivity-atlas stage 2a (`TIER2_SENSITIVITY_ATLAS_PLAN.md` section 6.2):
before any artifact-derived table is trusted, the machinery must reproduce
the committed shared pure-cubic result. Two parts, both **read-only** with
respect to ``data/reference/`` (nothing is written anywhere):

* **Part A -- point evaluation.** Evaluate the joint objective at the stored
  ``shared_pure_cubic`` point (a = 0, stored b, stored shared E_bind) and
  compare objective + per-case RMSE against the committed
  ``verdict.json`` stage-2 record. Expected: bit-for-bit (the drag ion path
  is RNG-free and the neutral stage is seeded); any drift is reported in
  full ``repr`` precision.
* **Part B -- refit oracle.** Re-run the section-9 joint Nelder-Mead fit
  (``shared_pure_cubic`` variant, driver settings from
  ``method_b_shared_refit.py``) and require the best-fit b to reproduce the
  locked b = 2.5154 (PASS band |db| <= 5e-4, the rounding band of the locked
  four-decimal quote; the exact drift vs the stored artifact is printed).

Usage: edit USER SETTINGS, then::

    python scripts/extraction/atlas_d4_step1_oracle.py

Exit code 0 = all enabled parts PASS, 1 = any FAIL.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

# =============================================================================
# USER SETTINGS
# =============================================================================
N = 50                  # ensemble size (committed Tier-0 evidence point)
SEED = 20260604         # neutral-stage seed (ion drag path is RNG-free)
ION_TIME_PS = 20.0
DT_ION_PS = 0.01
MAXFEV_PER_START = 400
RUN_REFIT = True        # Part B (Part A alone is a quick machinery check)
B_LOCKED = 2.5154       # the locked headline value (four-decimal quote)
B_PASS_BAND = 5.0e-4    # |b_refit - B_LOCKED| PASS band (rounding half-width)

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    SHARED_PURE_CUBIC,
    build_case_setup,
    evaluate_joint_objective,
    fit_shared_trajectory_matching,
)
from i2_helium_md.presets import REFERENCE_DRAG_ROOT  # noqa: E402

CASES = ("18A", "9A")
VERDICT_PATH = (
    REFERENCE_DRAG_ROOT / "shared" / "trajectory_matching" / "verdict.json"
)


def _report(name: str, value: float, stored: float) -> bool:
    """Print one comparison line; bit-for-bit match is reported as such."""
    exact = value == stored
    drift = value - stored
    tag = "BIT-EXACT" if exact else f"drift {drift:+.3e}"
    print(f"    {name}: got {value!r} vs stored {stored!r} -> {tag}")
    return exact


def main() -> int:
    stored = json.loads(VERDICT_PATH.read_text(encoding="utf-8"))
    rec = stored["stage2"]["shared_pure_cubic"]
    b_stored = float(rec["b"])
    e_stored = {c: float(rec["e_bind_eV"][c]) for c in CASES}
    obj_stored = float(rec["objective_mean_Aps"])

    print(f"oracle: building case setups (N={N}, seed={SEED}) ...")
    t0 = time.perf_counter()
    setups = {
        c: build_case_setup(
            c, num_molecules=N, seed=SEED,
            ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        )
        for c in CASES
    }
    print(f"  setups built in {time.perf_counter() - t0:.0f} s")

    all_pass = True

    # --- Part A: point evaluation at the stored minimum -------------------
    print("Part A: joint objective at the stored shared_pure_cubic point")
    t0 = time.perf_counter()
    joint = evaluate_joint_objective(0.0, b_stored, e_stored, setups)
    dt_eval = time.perf_counter() - t0
    print(f"  one joint evaluation: {dt_eval:.1f} s")
    ok = _report("objective_mean_Aps", joint.objective, obj_stored)
    for c in CASES:
        ok &= _report(
            f"rmse_{c}_Aps",
            joint.per_case[c].rmse_Aps,
            float(rec["per_case"][c]["rmse_Aps"]),
        )
    print(f"  Part A: {'PASS (bit-exact)' if ok else 'FAIL (drift above)'}")
    all_pass &= ok

    # --- Part B: full refit oracle ----------------------------------------
    if RUN_REFIT:
        print("Part B: shared_pure_cubic joint refit "
              f"(maxfev/start={MAXFEV_PER_START}) ...")
        # Section-9.3 anchors come from the 18A **Method-A** bundle. The
        # original driver read them off the then-wired preset; the presets
        # are now re-wired to the shared bundle (a = 0), and the verdict
        # provenance forbids letting that condition a fit -- so the oracle
        # reads the Method-A artifact directly.
        ma = json.loads(
            (REFERENCE_DRAG_ROOT / "18A" / "linear_and_cubic"
             / "fit_parameters.json").read_text(encoding="utf-8")
        )
        anchors = (float(ma["a"]), float(ma["b"]))
        print(f"  anchors (18A Method-A only): a0={anchors[0]:.4f}, "
              f"b0={anchors[1]:.4f}")
        t0 = time.perf_counter()
        fit = fit_shared_trajectory_matching(
            setups, variant=SHARED_PURE_CUBIC, anchors=anchors,
            maxfev_per_start=MAXFEV_PER_START,
        )
        print(f"  refit done in {time.perf_counter() - t0:.0f} s, "
              f"{fit.n_evaluations} joint evals, converged={fit.converged}")
        best = fit.best
        _report("b_refit_vs_stored", best.b, b_stored)
        _report("e_bind_refit_vs_stored_eV",
                best.e_bind_eV["18A"], e_stored["18A"])
        _report("objective_refit_vs_stored", best.objective, obj_stored)
        db = abs(best.b - B_LOCKED)
        ok = db <= B_PASS_BAND
        print(f"  b_refit = {best.b!r} vs locked {B_LOCKED} "
              f"(|db| = {db:.2e}, band {B_PASS_BAND:.0e}) -> "
              f"{'PASS' if ok else 'FAIL'}")
        all_pass &= ok

    print(f"ORACLE VERDICT: {'PASS' if all_pass else 'FAIL'}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
