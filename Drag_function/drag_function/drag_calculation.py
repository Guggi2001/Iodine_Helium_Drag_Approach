from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.signal import savgol_filter
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import curve_fit
import pandas as pd
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Load Config and Data
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
from config_utils_local import config as C
from drag_function import io

def reconstruct_R_from_v(v_w, t_w, R0):
    """
    Reconstruct R from smoothed velocity data using Savitzky-Golay filter.
    """
    return (R0 + cumulative_trapezoid(v_w, t_w, initial=0))
t_exit_9 = 8.5

test = True
if test:
    dict9 = io.load_data(C.PATH9A)
    dict18 = io.load_data(C.PATH18A)

    t = dict18["t"]
    v = dict18["v2"]
    v1 = dict18["v1"]
    mask = (t >= 4.543) & (t <= 8)
    t_w = t[mask]
    v_w = v[mask]
    v1_w = v1[mask]
    wl = 3401
    polyorder = 1
    v_s_18 = savgol_filter(v_w, window_length=wl, polyorder=polyorder, deriv=0, mode="interp")
    v1_s_18 = savgol_filter(v1_w, window_length=wl, polyorder=polyorder, deriv=0, mode="interp")

    plt.figure(figsize=(10, 4))
    plt.plot(t_w, v_w, color='gray', alpha=0.4, label='Original')
    plt.plot(t_w, v_s_18, color='#0072BD', lw=1.5, label='Smoothed')
    plt.xlabel('t / ps')
    plt.ylabel('v / Å/ps')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

    data = pd.read_csv('cleaned_data.csv')
    t_9 = data['time'].values
    v_9_SG = data['cleaned_SG'].values
    v_9_IMF = data['IMF_cleaned'].values
    t_9_full = dict9["t"]
    v_9 = dict9["v2"]
    v_9_1 = dict9["v1"]
    R_9 = dict9["R"]
    mask_9 = (t_9_full >= 2.67) & (t_9_full <= t_exit_9)
    t_9_w = t_9_full[mask_9]
    v_9_w = v_9[mask_9]
    R_9_w = R_9[mask_9]
    v_9_1_w = v_9_1[mask_9]

    wl = 4901
    polyorder = 1
    v9_test = savgol_filter(v_9_IMF, window_length=wl, polyorder=polyorder, deriv=0, mode="interp")



    plt.figure(figsize=(7, 4))
    plt.plot(t_9_w, v_9_1_w,  alpha=0.2, label='Given v Iodine 1')
    plt.plot(t_9_w, v_9_w, alpha=0.2, label='Given v Iodine 2')
    plt.plot(t_9, v_9_SG, lw=1.5, label='Smoothed (IMF + SG wl: 3901)')
    plt.plot(t_9, v_9_IMF,  lw=1.5, label='Smoothed (IMF)')
    plt.plot(t_9, v9_test, lw=1.5, label=f'Smoothed (IMF + SG wl: {wl})')
    plt.xlabel('t / ps')
    plt.ylabel('v / Å/ps')
    plt.title('Velocity curves of both atoms with smoothed versions for the 9Å case')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()


    # plt.figure(figsize=(7, 4))
    # plt.plot(t_9_w, R_9_w, label='Given R')
    # plt.plot(t_9_w, reconstruct_R_from_v(2 * v_9_IMF, t_9_w, R0=R_9_w[0]), label='Reconstructed R')
    # plt.title('Given R and Reconstructed R in relevant time window for 9Å case')
    # plt.xlabel('Time (ps)')
    # plt.ylabel('Distance Å')
    # plt.grid(True, alpha=0.3)
    # plt.legend()
    # plt.show()

    t18 = dict18["t"]
    R18 = dict18["R"]
    v18 = dict18["v2"]
    v18_1 = dict18['v1']
    mask = (t18 >= 4.543) & (t18 <= 8)
    t_18_w = t18[mask]
    R_18_w = R18[mask]
    v_18_w = v18[mask]
    v_18_1_w = v18_1[mask]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

    # 9Å case
    ax1.plot(t_9_w, v_9_1_w, alpha=0.2, label='Given v Iodine 1')
    ax1.plot(t_9_w, v_9_w, alpha=0.2, label='Given v Iodine 2')
    ax1.plot(t_9, v_9_SG, lw=1.5, label='Smoothed velocity')
    ax1.set_xlabel('t / ps')
    ax1.set_ylabel('v / Å/ps')
    ax1.set_title('9Å case: Velocity curves with smoothed version')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 18Å case
    ax2.plot(t_18_w, v_18_1_w, alpha=0.4, label='Given v Iodine 1')
    ax2.plot(t_18_w, v_18_w, alpha=0.4, label='Given v Iodine 2')
    ax2.plot(t_18_w, v_s_18, label='Smoothed velocity Iodine 2')
    ax2.plot(t_18_w, v1_s_18, ls = '--', label='Smoothed velocity Iodine 1')
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Velocity (Å/ps)')
    ax2.set_title('18Å case: Velocity curves with smoothed version')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    pass
#/////////////////////////////////////////////////////////////////////////////////////////////
#                                     Actual Run Tests
#////////////////////////////////////////////////////////////////////////////////////////////

"""
Instantaneous drag extraction from TD-DFT trajectories using spline differentiation.

Pipeline (trusted interior only)
--------------------------------
1) Inputs:
   - t(t): time [ps]
   - R(t): separation [Å]
   - v(t): speed (or longitudinal component) [Å/ps]  (already smoothed in your chosen window)

2) Compute Coulomb driving force:
   F_C(R) = k_e * q1*q2 / R^2   in [amu*Å/ps^2]

3) Compute acceleration via cubic spline:
   a(t) = dv/dt from spline derivative  in [Å/ps^2]

4) Force balance:
   F_drag(t) = F_C(t) - m_eff * a(t)    in [amu*Å/ps^2]

5) Characterize drag law:
   Plot F_drag vs v and fit |F_drag| ≈ gamma * v^n

Notes
-----
- Spline derivatives are unreliable near edges => we compute on a trusted interior,
  truncating first/last N points.
- Units are explicit on every variable, per your request.
"""




# =============================================================================
# Constants & unit conversions
# =============================================================================

# Masses [amu]
M_IODINE_AMU = 126.90447        # iodine atom mass [amu]
M_HELIUM_AMU = 4.002602         # helium atom mass [amu]

# Effective mass (your definition) [amu]
# meff = m_I + 21*m_He  (no backflow term here; add later if needed)
N_HE_SHELL = 21
M_EFF_AMU = M_IODINE_AMU + N_HE_SHELL * M_HELIUM_AMU  # [amu]

# Coulomb constant in eV·Å for elementary charges (+1e, +1e)
K_E_EV_A = 14.3996454784255     # [eV*Å]

# Given conversion:
# 1 eV/Å ≈ 9648.5 amu*Å/ps^2  (user-provided; use as authoritative here)
EV_PER_A_TO_AMU_A_PER_PS2 = 9648.5  # [(amu*Å/ps^2) / (eV/Å)]


# =============================================================================
# Settings
# =============================================================================

@dataclass(frozen=True)
class DragExtractionSettings:
    q1: int = 1                              # charge number for ion 1 [+e]
    q2: int = 1                              # charge number for ion 2 [+e]
    meff_amu: float = M_EFF_AMU              # effective mass [amu]

    # Spline fit settings
    spline_k: int = 3                        # cubic spline
    spline_s: Optional[float] = None         # smoothing parameter s (None => auto heuristic)

    # Trusted interior handling
    truncate_points: int = 500               # drop first/last N points (dt=0.001 ps => 0.5 ps)
    # If you prefer time-based truncation, set truncate_points=None and use truncate_time_ps
    truncate_time_ps: Optional[float] = None # [ps]

    # Diagnostics
    show_plots: bool = True
    fit_power_law: bool = True               # fit |F_drag| = gamma * v^n
    fit_variant: int = 1                    # 1: fit |F_drag| = gamma * v^n, 2: fit F_drag = a v + b v^3
    case: int = None,                        # Optional for title in plots (e.g. 9 or 18 Å case)
    plot_only_drag_with_fit: bool = False,




# =============================================================================
# Core helpers
# =============================================================================

def check_uniform_dt(t_ps: np.ndarray, rtol: float = 1e-3) -> float:
    """Return dt [ps] and validate uniform sampling."""
    t_ps = np.asarray(t_ps, float)
    dt = np.diff(t_ps)
    dt_mean = float(np.mean(dt))
    if not np.allclose(dt, dt_mean, rtol=rtol, atol=1e-12):
        raise ValueError("Non-uniform dt detected. Resample before spline differentiation.")
    return dt_mean


def coulomb_force_amuAps2(R_A: np.ndarray, q1: int = 1, q2: int = 1) -> np.ndarray:
    """
    Coulomb force magnitude for two +1e charges separated by R.

    Inputs
    ------
    R_A : separation [Å]

    Returns
    -------
    F_C : Coulomb force magnitude [amu*Å/ps^2]
    """
    R_A = np.asarray(R_A, float)
    # Force in [eV/Å]
    F_eV_per_A = (K_E_EV_A * q1 * q2) / (R_A ** 2)  # [eV/Å]
    # Convert to [amu*Å/ps^2]
    F_amuAps2 = F_eV_per_A * EV_PER_A_TO_AMU_A_PER_PS2  # [amu*Å/ps^2]
    return F_amuAps2


def choose_spline_s(v_Aps: np.ndarray) -> float:
    """
    Heuristic for UnivariateSpline smoothing parameter s.

    UnivariateSpline minimizes sum(w_i*(y_i - s(x_i))^2) <= s.
    Larger s => smoother spline.

    We estimate noise scale from a robust MAD of first differences.
    This keeps spline close to data while damping small-scale chatter.
    """
    v = np.asarray(v_Aps, float)
    dv = np.diff(v)
    mad = np.median(np.abs(dv - np.median(dv)))
    sigma = 1.4826 * mad  # robust std estimate in [Å/ps]
    # scale s with N * sigma^2  (dimension: (Å/ps)^2)
    N = len(v)
    s = float(N * (sigma ** 2))
    # Guard: if sigma ~0 (already super smooth), allow very small s
    return max(s, 1e-12)


def trusted_interior_mask(t_ps: np.ndarray, settings: DragExtractionSettings) -> np.ndarray:
    """
    Build a boolean mask selecting the trusted interior region.
    This avoids spline edge artifacts in derivatives.
    """
    t_ps = np.asarray(t_ps, float)
    n = len(t_ps)

    if settings.truncate_time_ps is not None:
        t0 = t_ps[0] + settings.truncate_time_ps
        t1 = t_ps[-1] - settings.truncate_time_ps
        return (t_ps >= t0) & (t_ps <= t1)

    # point-based truncation (default)
    N = int(settings.truncate_points)
    if 2 * N >= n:
        raise ValueError("truncate_points too large for the data length.")
    mask = np.zeros(n, dtype=bool)
    mask[N: n - N] = True
    return mask


# =============================================================================
# Main calculation
# =============================================================================

def main_calculation(
    t_ps: np.ndarray,         # [ps]
    R_A: np.ndarray,          # [Å]
    v_A_per_ps: np.ndarray,        # [Å/ps]   (already smoothed velocity in chosen window)
    settings: DragExtractionSettings = DragExtractionSettings(),
) -> Dict[str, np.ndarray]:
    """
    Compute F_C(t), acceleration a(t) from spline, and F_drag(t).

    Returns a dict with arrays (all same length as input t), but note that only
    the trusted interior is meaningful for derivatives and drag law fitting.
    """
    t_ps = np.asarray(t_ps, float)
    R_A = np.asarray(R_A, float)
    v_A_per_ps = np.asarray(v_A_per_ps, float)

    if not (len(t_ps) == len(R_A) == len(v_A_per_ps)):
        raise ValueError("t, R, v must have the same length.")

    dt_ps = check_uniform_dt(t_ps)  # [ps]

    # --- spline smoothing parameter ---
    spline_s = settings.spline_s
    if spline_s is None:
        spline_s = choose_spline_s(v_A_per_ps)

    # --- spline fit: v(t) ---
    # Units: v_A_per_ps [Å/ps], t_ps [ps] => derivative dv/dt [Å/ps^2]
    spline = UnivariateSpline(t_ps, v_A_per_ps, k=settings.spline_k, s=spline_s)

    v_spline_Aps = spline(t_ps)                           # [Å/ps]
    a_spline_Aps2 = spline.derivative(1)(t_ps)            # [Å/ps^2]

    # --- Coulomb force ---
    F_C_amuAps2 = coulomb_force_amuAps2(R_A, settings.q1, settings.q2)  # [amu*Å/ps^2]

    # --- drag force from force balance ---
    # meff [amu] * a [Å/ps^2] => [amu*Å/ps^2]
    F_inert_amuAps2 = settings.meff_amu * a_spline_Aps2                # [amu*Å/ps^2]
    F_drag_amuAps2 = F_C_amuAps2 - F_inert_amuAps2                     # [amu*Å/ps^2]

    # --- trusted interior mask ---
    mask = trusted_interior_mask(t_ps, settings)

    # --- diagnostics / sanity printouts ---
    # Compute means on trusted interior only
    mean_FC = float(np.mean(F_C_amuAps2[mask]))
    mean_inert = float(np.mean(F_inert_amuAps2[mask]))
    mean_drag = float(np.mean(F_drag_amuAps2[mask]))

    print("\n=== Force magnitude sanity check (trusted interior) ===")
    print(f"dt                       = {dt_ps:.6f} ps")
    print(f"meff                     = {settings.meff_amu:.3f} amu")
    print(f"mean(F_C)                = {mean_FC:.3e} amu*Å/ps^2")
    print(f"mean(meff * vdot)        = {mean_inert:.3e} amu*Å/ps^2")
    print(f"mean(F_drag = F_C - ...) = {mean_drag:.3e} amu*Å/ps^2")
    print(f"spline smoothing s       = {spline_s:.3e} (Å/ps)^2\n")

    # --- diagnostic plots ---
    if settings.show_plots:
        # Residual plot: v_data - v_spline
        resid = v_A_per_ps - v_spline_Aps  # [Å/ps]

        if not settings.plot_only_drag_with_fit:
            fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
            if settings.case is not None:
                fig.suptitle(f"Velocity and residual for {settings.case} Å case", fontweight='bold')
            ax[0].plot(t_ps, v_A_per_ps, "k", lw=0.8, label="v data [Å/ps]")
            ax[0].plot(t_ps, v_spline_Aps, "r--", lw=1.5, label="v spline [Å/ps]")
            ax[0].fill_between(t_ps, np.min(v_A_per_ps), np.max(v_A_per_ps), where=mask, alpha=0.08, label="trusted interior")
            ax[0].set_ylabel("v [Å/ps]")
            ax[0].grid(True)
            ax[0].legend(fontsize="small")

            ax[1].plot(t_ps, resid, lw=1.0, label="residual v_data - v_spline [Å/ps]")
            ax[1].axhline(0.0, color="k", lw=0.8, alpha=0.7)
            ax[1].fill_between(t_ps, np.min(resid), np.max(resid), where=mask, alpha=0.08)
            ax[1].set_xlabel("t [ps]")
            ax[1].set_ylabel("residual [Å/ps]")
            ax[1].grid(True)
            ax[1].legend(fontsize="small")

            plt.tight_layout()
            plt.show()

            # Acceleration and forces
            fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
            if settings.case is not None:
                fig.suptitle(f"Acceleration and forces for {settings.case} Å case", fontweight='bold')
            ax[0].plot(t_ps, a_spline_Aps2, lw=1.2, label="vdot from spline [Å/ps²]")
            ax[0].fill_between(t_ps, np.min(a_spline_Aps2), np.max(a_spline_Aps2), where=mask, alpha=0.08)
            ax[0].set_ylabel("vdot [Å/ps²]")
            ax[0].grid(True)
            ax[0].legend(fontsize="small")

            ax[1].plot(t_ps, F_C_amuAps2, lw=1.2, label="F_C [amu*Å/ps²]")
            ax[1].plot(t_ps, F_inert_amuAps2, lw=1.2, label="meff*vdot [amu*Å/ps²]")
            ax[1].plot(t_ps, F_drag_amuAps2, lw=1.2, label="F_drag [amu*Å/ps²]")
            ax[1].fill_between(t_ps, np.min(F_drag_amuAps2), np.max(F_drag_amuAps2), where=mask, alpha=0.08)
            ax[1].set_xlabel("t [ps]")
            ax[1].set_ylabel("Force [amu*Å/ps²]")
            ax[1].grid(True)
            ax[1].legend(fontsize="small")
            plt.tight_layout()
            plt.show()

            if settings.case is not None:
                # Acceleration forces and velocity
                fig, ax = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
                fig.suptitle(f"Velocity, Acceleration and forces for {settings.case} Å case", fontweight='bold')

                if settings.case == 9:
                    dict = io.load_data(C.PATH9A)
                elif settings.case == 18:
                    dict = io.load_data(C.PATH18A)
                else:
                    raise ValueError(f"Unknown case {settings.case}")
                v_orig = dict["v2"]
                t_full = dict["t"]
                mask_v = (t_full >= t_ps[0]) & (t_full <= t_ps[-1])
                v_orig_w = v_orig[mask_v]

                # data = pd.read_csv('cleaned_data.csv')
                # v_9_IMF = data['IMF_cleaned'].values
                # a_9_SG = savgol_filter(v_9_w, window_length=1201, polyorder=3, deriv=1, mode="interp")


                ax[0].plot(t_ps, v_A_per_ps, lw=1.2, color = 'red' , label="v from IMF+SG smoothing [Å/ps]")
                ax[0].plot(t_ps, v_orig_w, lw=1.2, alpha =0.3 , label="Original v [Å/ps]")
                ax[0].fill_between(t_ps, np.min(v_A_per_ps), np.max(v_A_per_ps), where=mask, alpha=0.08)
                ax[0].set_ylabel("v [Å/ps]")
                ax[0].grid(True)
                ax[0].legend(fontsize="small")

                ax[1].plot(t_ps, a_spline_Aps2, lw=1.2, label="vdot from spline [Å/ps²]")
                #ax[1].plot(t_ps, a_9_SG, lw=1.2, label="vdot from IMF+SG smoothing [Å/ps²]")
                ax[1].fill_between(t_ps, np.min(a_spline_Aps2), np.max(a_spline_Aps2), where=mask, alpha=0.08)
                ax[1].set_ylabel("vdot [Å/ps²]")
                ax[1].grid(True)
                ax[1].legend(fontsize="small")

                ax[2].plot(t_ps, F_C_amuAps2, lw=1.2, label="F_C [amu*Å/ps²]")
                ax[2].plot(t_ps, F_inert_amuAps2, lw=1.2, label="meff*vdot [amu*Å/ps²]")
                ax[2].plot(t_ps, F_drag_amuAps2, lw=1.2, label="F_drag [amu*Å/ps²]")
                ax[2].fill_between(t_ps, np.min(F_drag_amuAps2), np.max(F_drag_amuAps2), where=mask, alpha=0.08)
                ax[2].set_xlabel("t [ps]")
                ax[2].set_ylabel("Force [amu*Å/ps²]")
                ax[2].grid(True)
                ax[2].legend(fontsize="small")
                plt.tight_layout()
                plt.show()

        # --- Drag law plot and fit on trusted interior ---
        v_f = v_spline_Aps[mask]  # [Å/ps]
        F_f = F_drag_amuAps2[mask]  # [amu*Å/ps^2]
        # Keep only finite values (and optionally v>0 if v is a speed magnitude)
        mask_2 = np.isfinite(v_f) & np.isfinite(F_f) & (v_f > 0.0)
        v_fit = v_f[mask_2]
        F_fit = F_f[mask_2]

        # Many models assume drag opposes motion; if v is speed magnitude, F_drag should be positive magnitude.
        # If occasional negatives appear, they indicate either remaining oscillations, wrong meff, or non-Markovian effects.
        # For a power-law fit we fit |F_drag| = gamma * v^n on positive v.
    if settings.show_plots and not settings.plot_only_drag_with_fit:
        plt.figure(figsize=(8, 5))
        plt.scatter(v_fit, F_fit, s=10, alpha=0.6, label="F_drag(t) samples")
        plt.xlabel("v [Å/ps]")
        plt.ylabel("F_drag [amu*Å/ps²]")
        if settings.case is not None:
            plt.title(f"Instantaneous drag law: F_drag(v) (trusted interior) for {settings.case} Å case", fontweight='bold')
        else:
            plt.title("Instantaneous drag law: F_drag(v) (trusted interior)", fontweight='bold')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    fit_result = {}
    if settings.fit_power_law:
        # ---------------------------------------------------------------------------------------
        #---------------------------------------------------------------------------------------
        #           gamma*v^n fit (fit_variant=1) is more robust and interpretable for drag
        #           laws, but if you want to try a linear+v^3 fit, set fit_variant=2
        #---------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------
        if settings.fit_variant == 1:
            # Fit on positive finite values
            m = np.isfinite(v_fit) & np.isfinite(F_fit) & (v_fit > 0)
            v_pos = v_fit[m]
            F_pos = np.abs(F_fit[m])

            # Avoid log(0)
            m2 = F_pos > 0
            v_pos = v_pos[m2]
            F_pos = F_pos[m2]

            if len(v_pos) >= 20:
                # log-linear fit: log F = log gamma + n log v
                X = np.log(v_pos)
                Y = np.log(F_pos)
                A = np.vstack([np.ones_like(X), X]).T
                c0, n = np.linalg.lstsq(A, Y, rcond=None)[0]
                gamma = float(np.exp(c0))
                n = float(n)

                fit_result = {"gamma": gamma, "n": n}

                if settings.show_plots:
                    v_grid = np.linspace(np.min(v_pos), np.max(v_pos), 200)
                    F_grid = gamma * (v_grid ** n)
                    plt.figure(figsize=(8, 5))
                    plt.scatter(v_pos, F_pos, s=10, alpha=0.5, label="|F_drag| data")
                    plt.plot(v_grid, F_grid, lw=2.0, label=f"fit: |F| = γ v^n, γ={gamma:.3e}, n={n:.2f}")
                    plt.xlabel("v [Å/ps]")
                    plt.ylabel("|F_drag| [amu*Å/ps²]")
                    if settings.case is not None:
                        plt.title(f"Power-law fit of drag magnitude for {settings.case} Å case", fontweight='bold')
                    else:
                        plt.title("Power-law fit of drag magnitude", fontweight='bold')
                    plt.grid(True)
                    plt.legend()
                    plt.tight_layout()
                    plt.show()

                print("=== Drag law fit (trusted interior) ===")
                print(f"|F_drag| ≈ gamma * v^n")
                print(f"gamma = {gamma:.3e} [amu*Å/ps² / (Å/ps)^n]")
                print(f"n     = {n:.3f}\n")
            else:
                print("Not enough valid samples for power-law fit.")

        # ---------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------
        #                   A linear+v^3 fit, set fit_variant=2
        # ---------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------
        elif settings.fit_variant == 2:

            def drag_model(v, a, b):
                """F_drag(v) = a*v + b*v^3"""
                return a * v + b * (v ** 3)

            # Sensible initial guesses:
            # - a from small-v slope (robust: fit first ~10% of v-range)
            # - b initially 0 or small
            idx = np.argsort(v_fit)
            v_sorted = v_fit[idx]
            F_sorted = F_fit[idx]
            n_small = max(10, int(0.1 * len(v_sorted)))
            v_small = v_sorted[:n_small]
            F_small = F_sorted[:n_small]

            # Linear slope estimate for a (units: [amu*Å/ps^2] / [Å/ps] = [amu/ps])
            a0 = float(np.linalg.lstsq(v_small.reshape(-1, 1), F_small, rcond=None)[0][0]) if np.any(v_small) else 0.0

            # Crude b guess from large-v region if available (units: [amu*Å/ps^2] / [Å/ps]^3 = [amu*ps/Å^2])
            # Use median of (F - a0*v)/v^3 over top 10% velocities; fallback to 0.0
            n_large = max(10, int(0.1 * len(v_sorted)))
            v_large = v_sorted[-n_large:]
            F_large = F_sorted[-n_large:]
            b_candidates = (F_large - a0 * v_large) / np.maximum(v_large ** 3, 1e-30)
            b0 = float(np.median(b_candidates[np.isfinite(b_candidates)])) if np.any(np.isfinite(b_candidates)) else 0.0


            p0 = (a0, b0)

            # Fit with curve_fit (nonlinear least squares)
            popt, pcov = curve_fit(drag_model, v_fit, F_fit, p0=p0, maxfev=20000)
            a_hat, b_hat = popt
            perr = np.sqrt(np.diag(pcov))  # 1σ parameter uncertainties
            a_err, b_err = perr

            # R^2 calculation
            F_pred = drag_model(v_fit, a_hat, b_hat)
            ss_res = float(np.sum((F_fit - F_pred) ** 2))
            ss_tot = float(np.sum((F_fit - np.mean(F_fit)) ** 2))
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
            fit_result = {"a": a_hat, "b": b_hat, "a_err": a_err, "b_err": b_err}
            print("=== av + b v^3 drag fit ===")
            print(f"a = {a_hat:.6e} ± {a_err:.6e}  [amu/ps]  (since F/[v] = (amu*Å/ps^2)/(Å/ps))")
            print(f"b = {b_hat:.6e} ± {b_err:.6e}  [amu*ps/Å^2]  (since F/[v^3] = (amu*Å/ps^2)/(Å^3/ps^3))")
            print(f"R^2 = {r2:.6f}")

            # Plot overlay: raw samples + fitted curve
            v_grid = np.linspace(np.min(v_fit), np.max(v_fit), 400)
            F_grid = drag_model(v_grid, a_hat, b_hat)

            if settings.show_plots:
                plt.figure(figsize=(8, 5))
                plt.scatter(v_fit, F_fit, s=10, alpha=0.6, label="F_drag samples")
                plt.plot(v_grid, F_grid, lw=2.0, label= r"fit: $F_{\rm drag}(v)" + f"=a v + b v^3$, a={a_hat:.2f}, b={b_hat:.2f}")
                plt.xlabel("v [Å/ps]")
                plt.ylabel("F_drag [amu·Å/ps²]")
                plt.title("Drag law fit: av + b v³")
                plt.grid(True)
                plt.legend()
                plt.tight_layout()
                plt.show()


    return {
        "t_ps": t_ps,                        # [ps]
        "R_A": R_A,                          # [Å]
        "v_data_Aps": v_A_per_ps,                 # [Å/ps]
        "v_spline_Aps": v_spline_Aps,        # [Å/ps]
        "a_spline_Aps2": a_spline_Aps2,      # [Å/ps^2]
        "F_C_amuAps2": F_C_amuAps2,          # [amu*Å/ps^2]
        "F_inert_amuAps2": F_inert_amuAps2,  # [amu*Å/ps^2]
        "F_drag_amuAps2": F_drag_amuAps2,    # [amu*Å/ps^2]
        "trusted_mask": mask,                # [bool]
        "fit_result": fit_result,            # dict
        "spline_s": np.array([spline_s]),    # [(Å/ps)^2]
    }


# =============================================================================
# Example usage (adapt to your dict9/dict18 arrays)
# =============================================================================
def drag_model(v,gamma,n):
    """F_drag(v) = gamma * v^n"""
    return gamma * (v ** n)

def drag_model_tert(v, a, b):
    """F_drag(v) = a*v + b*v^3"""
    return a * v + b * (v ** 3)


dict9 = io.load_data(C.PATH9A)
dict18 = io.load_data(C.PATH18A)

data = pd.read_csv('cleaned_data.csv')
t_9 = data['time'].values
v_9_SG = data['cleaned_SG'].values
v_9_IMF = data['IMF_cleaned'].values
t_9_full = dict9["t"]
v_9 = dict9["v2"]
mask_9 = (t_9_full >= 2.67) & (t_9_full <= t_exit_9)
t_9_w = t_9_full[mask_9]
v_9_w = v_9[mask_9]
R_9 = dict9["R"]
R_9_w = R_9[mask_9]
v_9_SG_orig = savgol_filter(v_9_w, window_length=3901, polyorder=1, deriv=0, mode="interp")
v_9_SG_varying_wl = savgol_filter(v_9_IMF, window_length=2501, polyorder=1, deriv=0, mode="interp")
R_recon = reconstruct_R_from_v(2 * v_9_IMF, t_9_w, R0=R_9_w[0])

out = main_calculation(t_9_w, R_9_w, v_9_SG, DragExtractionSettings(case = 9,
                              truncate_points=500, plot_only_drag_with_fit = True))
out_fit_2 = main_calculation(t_9_w, R_9_w, v_9_SG, DragExtractionSettings(case = 9,
                              truncate_points=500, fit_variant = 2, plot_only_drag_with_fit = True))


t18 = dict18["t"]
v18 = dict18["v2"]
R18 = dict18["R"]
mask = (t18 >= 4.54) & (t18 <= 8)
t_18_w = t18[mask]
v_18_w = v18[mask]
R_18_w = R18[mask]
wl = 3401
polyorder = 1
v_18_SG = savgol_filter(v_18_w, window_length=wl, polyorder=polyorder, deriv=0, mode="interp")

out_18 = main_calculation(t_18_w, R_18_w, v_18_SG, DragExtractionSettings(case = 18, truncate_points=500,
                                                  plot_only_drag_with_fit = True))
out_18_fit_2 = main_calculation(t_18_w, R_18_w, v_18_SG, DragExtractionSettings(case = 18, truncate_points=500,
                                 fit_variant = 2, plot_only_drag_with_fit = True))

a = 3

def compare_drag_two_cases(v_9, F_drag_9, gamma_9, n_9, v_18, F_drag_18, gamma_18, n_18, v_fit_lims = None):
    plt.figure(figsize=(8, 5))
    plt.plot(v_9, F_drag_9, label="F_drag 9Å case")
    plt.plot(v_18, F_drag_18, label="F_drag 18Å case")
    if v_fit_lims is not None:
        v_fit = np.linspace(v_fit_lims[0], v_fit_lims[1], 10000)
    else:
        v_fit = np.linspace(np.min(v_18), np.max(v_9), 10000)
    plt.plot(v_fit, drag_model(v_fit, gamma_9, n_9), ls = "--", alpha = 0.4, label=f"Fit_9, gamma = {gamma_9:.2f}, n = {n_9:.2f}")
    plt.plot(v_fit, drag_model(v_fit, gamma_18, n_18), ls = "--", alpha = 0.4, label=f"Fit_18, gamma = {gamma_18:.2f}, n = {n_18:.2f}")
    plt.xlabel("v [Å/ps]")
    plt.ylabel("F_drag [amu*Å/ps²]")
    plt.title("Comparison of F_drag between 9Å and 18Å cases with $\\gamma$$v^n$ fit")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def compare_drag_two_cases_tert(v_9, F_drag_9, a_9, b_9, v_18, F_drag_18, a_18, b_18, v_fit_lims = None):
    plt.figure(figsize=(8, 5))
    plt.plot(v_9, F_drag_9, label="F_drag 9Å case")
    plt.plot(v_18, F_drag_18, label="F_drag 18Å case")
    if v_fit_lims is not None:
        v_fit = np.linspace(v_fit_lims[0], v_fit_lims[1], 10000)
    else:
        v_fit = np.linspace(np.min(v_18), np.max(v_9), 10000)
    plt.plot(v_fit, drag_model_tert(v_fit, a_9, b_9), ls = "--", alpha = 0.4, label=f"Fit_9, a = {a_9:.2f}, b = {b_9:.2f}")
    plt.plot(v_fit, drag_model_tert(v_fit, a_18, b_18), ls = "--", alpha = 0.4, label=f"Fit_18, a = {a_18:.2f}, b = {b_18:.2f}")
    plt.xlabel("v [Å/ps]")
    plt.ylabel("F_drag [amu*Å/ps²]")
    plt.title("Comparison of F_drag between 9Å and 18Å cases with av+bv^3")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()





compare_drag_two_cases_tert(out_fit_2["v_spline_Aps"], out_fit_2["F_drag_amuAps2"], out_fit_2["fit_result"]["a"], out_fit_2["fit_result"]["b"],
                       out_18_fit_2["v_spline_Aps"], out_18_fit_2["F_drag_amuAps2"], out_18_fit_2["fit_result"]["a"], out_18_fit_2["fit_result"]["b"], v_fit_lims = (0.5, 15))

compare_drag_two_cases(out["v_spline_Aps"], out["F_drag_amuAps2"], out["fit_result"]["gamma"], out["fit_result"]["n"],
                       out_18["v_spline_Aps"], out_18["F_drag_amuAps2"], out_18["fit_result"]["gamma"], out_18["fit_result"]["n"], v_fit_lims = (0.5, 15))