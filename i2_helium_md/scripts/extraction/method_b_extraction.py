"""Method-B trajectory-matching extraction driver (METHOD_B doc, plan B4).

Runs the joint ``{a, b, E_bind}`` fit for one case against the same-smoothed
long reference (``cleaned_data_long.csv``) over the full post-dynamic-start
window, with the escape penalty, multi-start degeneracy probing, a
neutral-seed sweep, and RMSE-sensitivity half-widths; writes the provenance-
stamped bundle to ``data/reference/drag/<case>/trajectory_matching/``.

The §6.5.1 pairing is exact inside the loop (each trial stamps its candidate
binding into the bundle), so no guard escape hatch is in play. The held-out
validation (cross-case + VMI) is SEPARATE -- see
``scripts/extraction/method_b_cross_case_check.py``; a fit produced here is
*calibrated, not validated*.

DRY_RUN mode evaluates the objective twice at the Method-A-era start point,
asserts the two results are bitwise identical (the ion drag path must be
RNG-free; any stochasticity leaking in fails loudly here), prints both, and
writes nothing.

Usage: edit USER SETTINGS, then::

    python scripts/extraction/method_b_extraction.py
"""

from __future__ import annotations

from pathlib import Path
import sys
import time

# =============================================================================
# USER SETTINGS
# =============================================================================
CASE = "9A"             # "9A" or "18A"; fit 18A first (the trustworthy case)
N = 50                  # ensemble size per forward integration (Tier-0 evidence)
SEED = 20260604         # neutral-stage seed (ion drag path is RNG-free)
SEED_SWEEP_SEEDS = [20260605, 20260606, 20260607, 20260608]  # K-1 extra seeds
ION_TIME_PS = 20.0      # covers the [t*, ~14 ps] window with escape margin
DT_ION_PS = 0.01
DRY_RUN = False         # True: 1 eval x2 + determinism assert, write nothing
MAXFEV_PER_START = 400
SEED_SWEEP_MAXFEV = 120  # warm-started refits need fewer evaluations

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402
from scipy.optimize import Bounds, minimize  # noqa: E402

from i2_helium_md.extraction.trajectory_matching import (  # noqa: E402
    E_BIND_STATIC_EV,
    build_case_setup,
    evaluate_objective,
    fit_trajectory_matching,
    sensitivity_halfwidths,
    write_fit_parameters,
)
from i2_helium_md.presets import REFERENCE_DRAG_ROOT  # noqa: E402


def dry_run(setup) -> None:
    theta0 = (setup.a0, setup.b0, 0.2)
    t0 = time.perf_counter()
    r1 = evaluate_objective(theta0, setup)
    t1 = time.perf_counter()
    r2 = evaluate_objective(theta0, setup)
    print(f"  eval at (a0={setup.a0:.4f}, b0={setup.b0:.4f}, E=0.2 eV), "
          f"{t1 - t0:.2f} s/eval:")
    print(f"    objective = {r1.objective:.6f} A/ps  "
          f"(rmse {r1.rmse_Aps:.6f}, escape {r1.escape_fraction:.3f}, "
          f"penalty {r1.penalty_Aps:.6f})")
    if r1.objective != r2.objective or r1.rmse_Aps != r2.rmse_Aps:
        raise AssertionError(
            "objective is NOT deterministic: repeated evaluation differs "
            f"({r1.objective!r} vs {r2.objective!r}) -- stochasticity has "
            "leaked into the ion drag path"
        )
    print("    repeated evaluation bitwise identical -> deterministic OK")
    print("  DRY_RUN: nothing written.")


def warm_start_refit(case: str, seed: int, x_best_normalized, *, scale) -> dict:
    """One warm-started NM refit on a fresh neutral seed (seed-sweep point)."""
    setup_k = build_case_setup(
        case, num_molecules=N, seed=seed,
        ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
    )
    lower = np.array([0.1, 1.0e-3, 0.01 / E_BIND_STATIC_EV])
    upper = np.array([10.0, 10.0, 1.0])

    def objective_x(x):
        return evaluate_objective(np.asarray(x) * scale, setup_k).objective

    res = minimize(
        objective_x,
        np.clip(np.asarray(x_best_normalized, dtype=float), lower, upper),
        method="Nelder-Mead",
        bounds=Bounds(lower, upper),
        options={"fatol": 1e-3, "xatol": 1e-3, "maxfev": SEED_SWEEP_MAXFEV},
    )
    theta = np.asarray(res.x, dtype=float) * scale
    final = evaluate_objective(theta, setup_k)
    return {
        "seed": seed,
        "a": final.a,
        "b": final.b,
        "e_bind_eV": final.e_bind_eV,
        "objective_Aps": final.objective,
        "rmse_Aps": final.rmse_Aps,
        "escape_fraction": final.escape_fraction,
        "converged": bool(res.success),
    }


def main() -> int:
    print(f"===== Method-B trajectory-matching extraction: {CASE} =====")
    print(f"  N={N}, seed={SEED}, ion {ION_TIME_PS} ps @ dt={DT_ION_PS} ps")
    print("  building case setup (one neutral run) ...")
    setup = build_case_setup(
        CASE, num_molecules=N, seed=SEED,
        ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
    )
    print(f"  reference: {setup.reference_path.name}, window "
          f"[{setup.window[0]:.2f}, {setup.window[1]:.2f}] ps; "
          f"anchors a0={setup.a0:.4f}, b0={setup.b0:.4f}")

    if DRY_RUN:
        dry_run(setup)
        return 0

    # --- main fit (multi-start NM with E_bind pre-scan) ---
    print("  fitting (3 starts, E_bind pre-scan) ...")
    t0 = time.perf_counter()
    fit = fit_trajectory_matching(setup, maxfev_per_start=MAXFEV_PER_START)
    t1 = time.perf_counter()
    best = fit.best
    print(f"  done in {t1 - t0:.0f} s, {fit.n_evaluations} evaluations, "
          f"converged={fit.converged}")
    for i, r in enumerate(fit.starts):
        print(f"    start {i}: a={r.a:.4f}, b={r.b:.4f}, "
              f"E={r.e_bind_eV:.4f} eV, objective={r.objective:.4f}")
    print(f"  BEST: a={best.a:.4f} amu/ps, b={best.b:.4f} amu*ps/A^2, "
          f"E_bind={best.e_bind_eV:.4f} eV")
    print(f"        rmse={best.rmse_Aps:.4f} A/ps, "
          f"escape_fraction={best.escape_fraction:.3f}, "
          f"penalty={best.penalty_Aps:.4f}")

    if best.penalty_Aps != 0.0:
        print("  ABORT: best fit still traps ions (penalty != 0); not "
              "writing a bundle. Inspect the binding bound / drag scale.")
        return 1

    # --- seed sweep (warm-started refits on fresh neutral seeds) ---
    scale = np.array([setup.a0, setup.b0, E_BIND_STATIC_EV])
    x_best = np.array([best.a / setup.a0, best.b / setup.b0,
                       best.e_bind_eV / E_BIND_STATIC_EV])
    sweep = []
    for seed in SEED_SWEEP_SEEDS:
        print(f"  seed-sweep refit @ seed={seed} ...")
        sweep.append(warm_start_refit(CASE, seed, x_best, scale=scale))
    sweep_all = sweep + [{
        "seed": SEED, "a": best.a, "b": best.b, "e_bind_eV": best.e_bind_eV,
        "objective_Aps": best.objective, "rmse_Aps": best.rmse_Aps,
        "escape_fraction": best.escape_fraction, "converged": fit.converged,
    }]
    a_std = float(np.std([r["a"] for r in sweep_all], ddof=1))
    b_std = float(np.std([r["b"] for r in sweep_all], ddof=1))
    e_std = float(np.std([r["e_bind_eV"] for r in sweep_all], ddof=1))
    print(f"  seed-sweep std: a {a_std:.4f}, b {b_std:.4f}, E {e_std:.4f}")

    # --- RMSE-sensitivity half-widths around the optimum ---
    print("  sensitivity scan (10% objective rise) ...")
    hw_a, hw_b, hw_e = sensitivity_halfwidths(setup, best)
    print(f"  sensitivity halfwidths: a {hw_a:.4f}, b {hw_b:.4f}, E {hw_e:.4f}")

    # Band = max(seed spread, sensitivity) per parameter -- a deliberately
    # conservative SENSITIVITY band, not a statistical CI (so labeled).
    fit.a_err = max(a_std, hw_a)
    fit.b_err = max(b_std, hw_b)
    fit.e_bind_err_eV = max(e_std, hw_e)
    fit.uncertainty_model = "seed_sweep_std_plus_rmse_sensitivity"
    fit.seed_sweep_seeds = [SEED] + list(SEED_SWEEP_SEEDS)
    fit.seed_sweep_results = sweep_all

    # Multi-start agreement check against the seed band (degeneracy ridge).
    spread_a = max(r.a for r in fit.starts) - min(r.a for r in fit.starts)
    spread_b = max(r.b for r in fit.starts) - min(r.b for r in fit.starts)
    spread_e = (max(r.e_bind_eV for r in fit.starts)
                - min(r.e_bind_eV for r in fit.starts))
    ridge = (spread_a > 3 * fit.a_err or spread_b > 3 * fit.b_err
             or spread_e > 3 * fit.e_bind_err_eV)
    if ridge:
        print("  WARNING: multi-start minima disagree beyond 3x the "
              "uncertainty band -- drag<->binding degeneracy ridge is live. "
              "Bundle is still written (spread is recorded in optimizer."
              "start_minima) but flag this for review before any preset "
              "re-wiring.")

    out_dir = REFERENCE_DRAG_ROOT / CASE / "trajectory_matching"
    out = write_fit_parameters(
        fit, out_dir, transverse_contaminated=(CASE == "9A"),
    )
    print(f"  wrote {out}")
    print("  NOTE: this fit is calibrated, NOT validated. Run "
          "method_b_cross_case_check.py (and eventually the VMI tier) "
          "before wiring it into presets.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
