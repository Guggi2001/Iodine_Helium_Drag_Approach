"""Statistical comparison of Python vs MATLAB hard-sphere collision results.

Run after run_python_collisions.py and run_matlab_collisions.m.

We can't compare individual atoms or events because the RNG streams differ.
Instead we compare distributional / statistical properties:

  1. Total collision count -- should agree to Poisson noise (~ sqrt(N) ~ 70)
  2. Mean E_kin(t), mean E_dissip(t) -- ensemble averages over 1000 atoms
  3. ΔE/E0 per-collision distribution -- should match in shape; we compare
     mean, std, and bin counts
  4. Direct comparison to the theoretical prediction
        <ΔE/E0> = 2ρ / (1+ρ)²    (uniform cosΘ)
        max ΔE/E0 = 1 - ((1-ρ)/(1+ρ))²
     for ρ = m_iodine/m_He = 127/4 = 31.75

Run:
    python compare_collisions.py
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


SCRIPT_DIR = Path(__file__).parent
PY_SUMMARY = SCRIPT_DIR / "python_summary.csv"
ML_SUMMARY = SCRIPT_DIR / "matlab_summary.csv"
PY_EVENTS = SCRIPT_DIR / "python_collision_events.csv"
ML_EVENTS = SCRIPT_DIR / "matlab_collision_events.csv"

# Theoretical predictions for ρ = 127/4 = 31.75
RHO = 127.0 / 4.0
PRED_MAX_FRAC_LOSS = 1.0 - ((1.0 - RHO) / (1.0 + RHO)) ** 2
PRED_MEAN_FRAC_LOSS = 2.0 * RHO / (1.0 + RHO) ** 2


def load_summary(path):
    return np.loadtxt(path, delimiter=",", skiprows=1)


def load_events(path):
    if not path.exists() or path.stat().st_size <= 30:  # only header
        return np.empty((0, 4))
    return np.loadtxt(path, delimiter=",", skiprows=1)


def main():
    for p in [PY_SUMMARY, ML_SUMMARY, PY_EVENTS, ML_EVENTS]:
        if not p.exists():
            print(f"ERROR: {p} not found. Run both run_python_collisions.py "
                  "and run_matlab_collisions.m first.")
            return

    py_sum = load_summary(PY_SUMMARY)
    ml_sum = load_summary(ML_SUMMARY)
    py_ev = load_events(PY_EVENTS)
    ml_ev = load_events(ML_EVENTS)

    print("=" * 78)
    print("Hard-sphere collision: Python vs MATLAB statistical comparison")
    print("=" * 78)

    # ---------------------------------------------------------------
    # 1. Collision counts
    # ---------------------------------------------------------------
    n_py = py_ev.shape[0]
    n_ml = ml_ev.shape[0]
    n_atoms = 1000  # from generate_init_state.py
    poisson_err = np.sqrt(0.5 * (n_py + n_ml))   # Poisson std for the average
    diff = abs(n_py - n_ml)
    sigma_dev = diff / max(poisson_err, 1)
    print(f"\n1. Total collision count")
    print(f"   Python: {n_py}")
    print(f"   MATLAB: {n_ml}")
    print(f"   Diff:   {diff}  (Poisson noise ~ {poisson_err:.0f}, "
          f"so {sigma_dev:.2f}σ)")
    pass1 = sigma_dev < 3
    print(f"   {'PASS' if pass1 else 'FAIL'}: agreement within 3σ")

    # ---------------------------------------------------------------
    # 2. Mean E_kin(t) and E_dissip(t)
    # ---------------------------------------------------------------
    print(f"\n2. Time-resolved ensemble means (1000 atoms each)")
    sample_t = [0.5, 1.0, 1.5, 2.0]   # ps
    print(f"   {'t [ps]':>8s}  {'Python <E_kin>':>16s}  "
          f"{'MATLAB <E_kin>':>16s}  {'diff [σ]':>10s}")
    pass2_kin = True
    for t_target in sample_t:
        i_py = int(round(t_target / 0.01))
        if i_py >= len(py_sum) or i_py >= len(ml_sum):
            continue
        E_py = py_sum[i_py, 1]
        E_ml = ml_sum[i_py, 1]
        var_py = py_sum[i_py, 2]
        var_ml = ml_sum[i_py, 2]
        # Standard error of the mean for ensemble of 1000:
        sem = np.sqrt((var_py + var_ml) / 2 / n_atoms)
        sigma_dev = abs(E_py - E_ml) / max(sem, 1e-30)
        print(f"   {t_target:>8.2f}  {E_py:>16.6f}  {E_ml:>16.6f}  "
              f"{sigma_dev:>10.2f}")
        if sigma_dev > 3:
            pass2_kin = False

    print(f"\n   {'t [ps]':>8s}  {'Python <Edis>':>16s}  "
          f"{'MATLAB <Edis>':>16s}  {'diff [σ]':>10s}")
    pass2_dis = True
    for t_target in sample_t:
        i_py = int(round(t_target / 0.01))
        if i_py >= len(py_sum) or i_py >= len(ml_sum):
            continue
        D_py = py_sum[i_py, 3]
        D_ml = ml_sum[i_py, 3]
        var_py = py_sum[i_py, 4]
        var_ml = ml_sum[i_py, 4]
        sem = np.sqrt((var_py + var_ml) / 2 / n_atoms)
        sigma_dev = abs(D_py - D_ml) / max(sem, 1e-30)
        print(f"   {t_target:>8.2f}  {D_py:>16.6f}  {D_ml:>16.6f}  "
              f"{sigma_dev:>10.2f}")
        if sigma_dev > 3:
            pass2_dis = False
    pass2 = pass2_kin and pass2_dis
    print(f"   {'PASS' if pass2 else 'FAIL'}: trajectory means agree within 3σ")

    # ---------------------------------------------------------------
    # 3. Per-collision ΔE/E0 distribution
    # ---------------------------------------------------------------
    if py_ev.shape[0] == 0 or ml_ev.shape[0] == 0:
        print("\n3. SKIPPED -- no collision events recorded")
        pass3 = False
    else:
        py_frac = py_ev[:, 2] / py_ev[:, 3]
        ml_frac = ml_ev[:, 2] / ml_ev[:, 3]
        print(f"\n3. Per-collision ΔE/E0 distribution")
        print(f"   {'metric':>12s}  {'Python':>14s}  {'MATLAB':>14s}  {'theory':>14s}")
        print(f"   {'mean':>12s}  {py_frac.mean():>14.6f}  "
              f"{ml_frac.mean():>14.6f}  {PRED_MEAN_FRAC_LOSS:>14.6f}")
        print(f"   {'min':>12s}  {py_frac.min():>14.6f}  "
              f"{ml_frac.min():>14.6f}  {0:>14.6f}")
        print(f"   {'max':>12s}  {py_frac.max():>14.6f}  "
              f"{ml_frac.max():>14.6f}  {PRED_MAX_FRAC_LOSS:>14.6f}")
        # SEM for the mean: σ(distribution) / sqrt(N)
        sem_py = py_frac.std(ddof=1) / np.sqrt(len(py_frac))
        sem_ml = ml_frac.std(ddof=1) / np.sqrt(len(ml_frac))
        sem = np.sqrt(sem_py ** 2 + sem_ml ** 2)
        sigma_dev = abs(py_frac.mean() - ml_frac.mean()) / max(sem, 1e-30)
        print(f"\n   Mean diff = {abs(py_frac.mean()-ml_frac.mean()):.6f}  "
              f"(SEM ~ {sem:.6f}, so {sigma_dev:.2f}σ)")
        pass3a = sigma_dev < 3

        # KS-style histogram test
        bins = np.linspace(0, PRED_MAX_FRAC_LOSS * 1.05, 25)
        h_py, _ = np.histogram(py_frac, bins=bins)
        h_ml, _ = np.histogram(ml_frac, bins=bins)
        # Normalize to compare shape
        n_py_total = h_py.sum()
        n_ml_total = h_ml.sum()
        # Chi-square style metric: per-bin expected variance ~ count
        chi2 = 0.0
        for i in range(len(bins) - 1):
            f_py = h_py[i] / n_py_total
            f_ml = h_ml[i] / n_ml_total
            # Variance in the count fractions
            var_combined = (f_py * (1 - f_py) / n_py_total
                            + f_ml * (1 - f_ml) / n_ml_total)
            if var_combined > 0:
                chi2 += (f_py - f_ml) ** 2 / var_combined
        chi2_per_bin = chi2 / (len(bins) - 1)
        # If consistent, chi2 / dof should be ~ 1
        print(f"   Chi-square per bin (normalized): {chi2_per_bin:.2f}  "
              f"(<5 is good agreement)")
        pass3b = chi2_per_bin < 5
        pass3 = pass3a and pass3b
        print(f"   {'PASS' if pass3 else 'FAIL'}: ΔE distribution agrees")

    # ---------------------------------------------------------------
    # 4. Comparison to theoretical predictions
    # ---------------------------------------------------------------
    if py_ev.shape[0] > 0:
        print(f"\n4. Theoretical predictions for ρ = m_I/m_He = 127/4 = {RHO}")
        print(f"   Predicted <ΔE/E0> = 2ρ/(1+ρ)² = {PRED_MEAN_FRAC_LOSS*100:.4f}%")
        print(f"   Predicted max ΔE/E0 = 1 - ((1-ρ)/(1+ρ))² = "
              f"{PRED_MAX_FRAC_LOSS*100:.4f}%")
        print(f"   Python mean: {py_frac.mean()*100:.4f}%, "
              f"max: {py_frac.max()*100:.4f}%")
        print(f"   MATLAB mean: {ml_frac.mean()*100:.4f}%, "
              f"max: {ml_frac.max()*100:.4f}%")
        # Sanity: are they BOTH close to theory?
        py_mean_off = abs(py_frac.mean() - PRED_MEAN_FRAC_LOSS) / PRED_MEAN_FRAC_LOSS
        ml_mean_off = abs(ml_frac.mean() - PRED_MEAN_FRAC_LOSS) / PRED_MEAN_FRAC_LOSS
        print(f"   Python deviation from theory: {py_mean_off*100:.3f}%")
        print(f"   MATLAB deviation from theory: {ml_mean_off*100:.3f}%")

    # ---------------------------------------------------------------
    # Verdict
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    overall_pass = pass1 and pass2 and pass3
    if overall_pass:
        print("OVERALL VERDICT: PASS")
        print("Python and MATLAB collision implementations agree statistically.")
    else:
        print("OVERALL VERDICT: FAIL")
        print("One or more tests rejected at 3σ. Check the disagreement above.")


if __name__ == "__main__":
    main()
