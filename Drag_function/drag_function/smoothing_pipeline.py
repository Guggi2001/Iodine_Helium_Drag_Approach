#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Load Config and Data
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
from __future__ import annotations
from config_utils_local import config as C
from drag_function import io


dict9 = io.load_data(C.PATH9A)
dict18 = io.load_data(C.PATH18A)


"""
Robustness runner for a CEEMDAN + SG hybrid denoising pipeline.

Goals
-----
1) Quantify sensitivity to CEEMDAN hyperparameters (noise_width, trials) and SG window.
2) Quantify sensitivity to temporal windowing (t_min/t_max variations).
3) Estimate uncertainty bands via repeated seeds (stochastic decomposition).

Assumptions
-----------
- Time array t is uniformly sampled (critical for period estimation and stable CEEMDAN behavior).
- Signal y is a scalar time-series (e.g., v1 magnitude in Å/ps).
- You reconstruct a "trend" by dropping IMF(s) near a target period (e.g., 1.2 ps) plus optional rules.

Dependencies
------------
pip install EMD-signal scipy numpy matplotlib pandas
"""



from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from PyEMD import CEEMDAN
from scipy.signal import savgol_filter


# -----------------------------
# Utilities: sampling & stats
# -----------------------------
def check_uniform_dt(t_ps: np.ndarray, rtol: float = 1e-3) -> float:
    """
    CEEMDAN and period estimators assume a constant sampling interval dt.
    Non-uniform sampling can cause:
      - distorted period estimates
      - unstable sifting/IMF extraction
    """
    t_ps = np.asarray(t_ps, dtype=float)
    dt = np.diff(t_ps)
    dt_mean = float(np.mean(dt))
    if not np.allclose(dt, dt_mean, rtol=rtol, atol=1e-12):
        raise ValueError("Time axis is not uniformly sampled enough for this pipeline.")
    return dt_mean


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return float(np.sqrt(np.mean((a - b) ** 2)))


def max_abs_dev(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return float(np.max(np.abs(a - b)))


# -----------------------------
# Period estimation (robust)
# -----------------------------
def estimate_period_zero_crossings(t_ps: np.ndarray, x: np.ndarray) -> Optional[float]:
    """
    Estimate average period from zero crossings.
    For an oscillatory mode, consecutive zero crossings are ~T/2 apart.

    Returns
    -------
    period_ps : float or None
        None if insufficient zero crossings.
    """
    t_ps = np.asarray(t_ps, float)
    x = np.asarray(x, float)

    # Sign changes mark zero crossings; interpolate crossing time for accuracy.
    s = np.sign(x)
    s[s == 0] = 1.0
    zc_idx = np.where(np.diff(s) != 0)[0]
    if zc_idx.size < 4:
        return None

    t_crossings: List[float] = []
    for i in zc_idx:
        x0, x1 = x[i], x[i + 1]
        if x1 == x0:
            continue
        frac = -x0 / (x1 - x0)
        t_crossings.append(float(t_ps[i] + frac * (t_ps[i + 1] - t_ps[i])))

    if len(t_crossings) < 4:
        return None

    t_crossings = np.asarray(t_crossings, float)
    half_periods = np.diff(t_crossings)
    half_period = float(np.median(half_periods))
    return 2.0 * half_period


# -----------------------------
# CEEMDAN + SG pipeline
# -----------------------------
@dataclass(frozen=True)
class CeemdanParams:
    trials: int = 200
    noise_width: float = 0.2
    seed: int = 0
    max_imfs: Optional[int] = None  # None => let CEEMDAN decide


@dataclass(frozen=True)
class SgParams:
    window_length: int = 1401  # samples; must be odd
    polyorder: int = 1


@dataclass(frozen=True)
class ModeSelectionParams:
    target_period_ps: float = 1.2
    rel_tol: float = 0.25
    # If mode is split across two IMFs (mode mixing), allow dropping a second IMF
    # if its amplitude is non-negligible relative to the strongest target IMF.
    allow_second_imf: bool = True
    second_imf_amp_ratio: float = 0.35


def run_ceemdan(
    t_ps: np.ndarray,
    y: np.ndarray,
    params: CeemdanParams,
) -> np.ndarray:
    """
    Run CEEMDAN decomposition and return IMFs: shape (n_imfs, n_samples).

    Reproducibility: set both global numpy seed and CEEMDAN internal seed.
    """
    _ = check_uniform_dt(t_ps)
    t_ps = np.asarray(t_ps, float)
    y = np.asarray(y, float)

    # Ensure reproducible ensemble noise
    np.random.seed(params.seed)

    ce = CEEMDAN(trials=params.trials, noise_width=params.noise_width)
    ce.random_seed = params.seed
    if params.max_imfs is not None:
        ce.MAX_IMF = int(params.max_imfs)

    imfs = ce.ceemdan(y, t_ps)
    return imfs


def identify_target_imfs(
    t_ps: np.ndarray,
    imfs: np.ndarray,
    sel: ModeSelectionParams,
) -> Tuple[List[int], List[Optional[float]], List[float]]:
    """
    Identify which IMF(s) correspond to the target mode (~1.2 ps).

    Uses:
      - period estimate via zero crossings
      - amplitude (RMS) ranking to avoid picking tiny IMFs
      - optional second IMF if mode mixing splits the target mode
    """
    periods: List[Optional[float]] = []
    amps: List[float] = []

    for k in range(imfs.shape[0]):
        T = estimate_period_zero_crossings(t_ps, imfs[k])
        periods.append(T)
        amps.append(float(np.sqrt(np.mean(imfs[k] ** 2))))

    lo = (1.0 - sel.rel_tol) * sel.target_period_ps
    hi = (1.0 + sel.rel_tol) * sel.target_period_ps
    candidates = [k for k, T in enumerate(periods) if (T is not None and lo <= T <= hi)]

    if not candidates:
        return [], periods, amps

    candidates_sorted = sorted(candidates, key=lambda k: amps[k], reverse=True)
    drop = [candidates_sorted[0]]

    if sel.allow_second_imf and len(candidates_sorted) > 1:
        k0 = candidates_sorted[0]
        k1 = candidates_sorted[1]
        if amps[k1] >= sel.second_imf_amp_ratio * amps[k0]:
            drop.append(k1)

    return drop, periods, amps


def reconstruct_trend(
    imfs: np.ndarray,
    drop_indices: Sequence[int],
) -> np.ndarray:
    """
    Reconstruct signal by summing all IMFs except those in drop_indices.
    """
    drop_set = set(drop_indices)
    keep = [k for k in range(imfs.shape[0]) if k not in drop_set]
    if not keep:
        raise ValueError("All IMFs dropped; cannot reconstruct.")
    return np.sum(imfs[keep, :], axis=0)


def apply_sg_smoother(y: np.ndarray, sg: SgParams) -> np.ndarray:
    """
    SG is used as a mild post-processor to remove residual 'ensemble chatter'
    without imposing a sharp frequency-domain cutoff (which can ring).
    """
    wl = int(sg.window_length)
    if wl % 2 == 0:
        wl += 1
    if wl >= len(y):
        raise ValueError("SG window too long for the available segment.")
    if sg.polyorder >= wl:
        raise ValueError("SG polyorder must be < window_length.")
    return savgol_filter(y, window_length=wl, polyorder=sg.polyorder, mode="interp")


def ceemdan_sg_trend(
    t_ps: np.ndarray,
    y: np.ndarray,
    ce_params: CeemdanParams,
    sg_params: SgParams,
    sel_params: ModeSelectionParams,
) -> Dict[str, object]:
    """
    Full pipeline on a given (t,y) segment:
      1) CEEMDAN -> IMFs
      2) drop IMF(s) near target period
      3) reconstruct trend (sum remaining IMFs)
      4) SG smooth the reconstructed trend

    Returns dict with trend, drop indices, IMF periods/amps, and config.
    """
    imfs = run_ceemdan(t_ps, y, ce_params)
    drop_idx, periods, amps = identify_target_imfs(t_ps, imfs, sel_params)
    trend_raw = reconstruct_trend(imfs, drop_idx)
    trend_sg = apply_sg_smoother(trend_raw, sg_params)

    return {
        "trend": trend_sg,
        "trend_raw": trend_raw,
        "drop_idx": drop_idx,
        "imf_periods": periods,
        "imf_amps": amps,
        "n_imfs": imfs.shape[0],
        "ce_params": ce_params,
        "sg_params": sg_params,
        "sel_params": sel_params,
    }


# -----------------------------
# Windowing helpers
# -----------------------------
def select_window(
    t_ps: np.ndarray,
    y: np.ndarray,
    t_min: float,
    t_max: float,
) -> Tuple[np.ndarray, np.ndarray]:
    t_ps = np.asarray(t_ps, float)
    y = np.asarray(y, float)
    mask = (t_ps >= t_min) & (t_ps <= t_max)
    return t_ps[mask], y[mask]


def overlap_mask(
    t_ps: np.ndarray,
    t_min: float,
    t_max: float,
) -> np.ndarray:
    t_ps = np.asarray(t_ps, float)
    return (t_ps >= t_min) & (t_ps <= t_max)


# -----------------------------
# Robustness runner
# -----------------------------
@dataclass(frozen=True)
class RobustnessSweep:
    # noise_widths: Tuple[float, ...] = (0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)
    # sg_windows: Tuple[int, ...] = (201, 401, 601, 801)
    noise_widths: Tuple[float, ...] = (0.15,0.25,)
    sg_windows: Tuple[int, ...] = (1401,)
    # Optional: trials sweep for convergence check
    trials_list: Tuple[int, ...] = (100, 200, 400)
    seeds: Tuple[int, ...] = tuple(range(10))  # for uncertainty / CI
    # Windowing sweep
    t_min_list: Tuple[float, ...] = (2.5, 2.7, 3.0, 3.5, 4.0)
    t_max_list: Tuple[float, ...] = (8.0, 8.5, 9.0)


def run_hyperparameter_sweep(
    t_ps: np.ndarray,
    y: np.ndarray,
    base_t_min: float,
    base_t_max: float,
    base_ce: CeemdanParams,
    base_sg: SgParams,
    sel: ModeSelectionParams,
    sweep: RobustnessSweep,
) -> pd.DataFrame:
    """
    Sweep noise_width and SG window length, measure deviation from a baseline trend.

    Baseline is defined at:
      noise_width = base_ce.noise_width
      trials      = base_ce.trials
      seed        = base_ce.seed
      SG window   = base_sg.window_length
    """
    t0, y0 = select_window(t_ps, y, base_t_min, base_t_max)

    baseline = ceemdan_sg_trend(t0, y0, base_ce, base_sg, sel)["trend"]

    rows = []
    for nw in sweep.noise_widths:
        for wl in sweep.sg_windows:
            ce_params = CeemdanParams(
                trials=base_ce.trials,
                noise_width=nw,
                seed=base_ce.seed,
                max_imfs=base_ce.max_imfs,
            )
            sg_params = SgParams(window_length=wl, polyorder=base_sg.polyorder)

            out = ceemdan_sg_trend(t0, y0, ce_params, sg_params, sel)
            trend = out["trend"]

            rows.append({
                "noise_width": nw,
                "trials": base_ce.trials,
                "seed": base_ce.seed,
                "sg_window": wl,
                "sg_polyorder": base_sg.polyorder,
                "t_min": base_t_min,
                "t_max": base_t_max,
                "n_imfs": out["n_imfs"],
                "drop_idx": str(out["drop_idx"]),
                "rmsd": rmsd(trend, baseline),
                "max_abs_dev": max_abs_dev(trend, baseline),
            })

    return pd.DataFrame(rows)


def run_windowing_sweep(
    t_ps: np.ndarray,
    y: np.ndarray,
    base_ce: CeemdanParams,
    base_sg: SgParams,
    sel: ModeSelectionParams,
    sweep: RobustnessSweep,
    compare_on_overlap: Tuple[float, float] = (4.0, 8.0),
) -> pd.DataFrame:
    """
    Sweep (t_min, t_max) and compare trends on a fixed overlap interval.
    This isolates end-effect sensitivity.

    compare_on_overlap defines where you compute RMSD/Δmax for all runs,
    ensuring the metric is evaluated on the same time region.
    """
    t_overlap_min, t_overlap_max = compare_on_overlap

    # Baseline window = median-ish choice (first entries by default)
    base_t_min = 2.67
    base_t_max = 8.5
    tb, yb = select_window(t_ps, y, base_t_min, base_t_max)
    baseline = ceemdan_sg_trend(tb, yb, base_ce, base_sg, sel)["trend"]

    # Map baseline to full t for overlap comparison
    # (since tb is a subwindow, we compare only on intersection with overlap)
    rows = []
    for t_min in sweep.t_min_list:
        for t_max in sweep.t_max_list:
            if t_max <= t_min:
                continue
            tt, yy = select_window(t_ps, y, t_min, t_max)
            if len(tt) < 10:
                continue

            out = ceemdan_sg_trend(tt, yy, base_ce, base_sg, sel)
            trend = out["trend"]

            # Evaluate on overlap in absolute time
            m_base = overlap_mask(tb, t_overlap_min, t_overlap_max)
            m_run = overlap_mask(tt, t_overlap_min, t_overlap_max)

            # Need same number of points to compare => interpolate run trend onto baseline overlap times
            t_base_ov = tb[m_base]
            if t_base_ov.size < 10:
                raise ValueError("Overlap interval too small on baseline window.")

            trend_base_ov = baseline[m_base]
            trend_run_ov = np.interp(t_base_ov, tt[m_run], trend[m_run])

            rows.append({
                "t_min": t_min,
                "t_max": t_max,
                "seed": base_ce.seed,
                "noise_width": base_ce.noise_width,
                "trials": base_ce.trials,
                "sg_window": base_sg.window_length,
                "rmsd_overlap": rmsd(trend_run_ov, trend_base_ov),
                "max_abs_dev_overlap": max_abs_dev(trend_run_ov, trend_base_ov),
                "n_imfs": out["n_imfs"],
                "drop_idx": str(out["drop_idx"]),
            })

    return pd.DataFrame(rows)


def compute_confidence_band(
    t_ps: np.ndarray,
    y: np.ndarray,
    t_min: float,
    t_max: float,
    base_ce: CeemdanParams,
    base_sg: SgParams,
    sel: ModeSelectionParams,
    seeds: Sequence[int],
) -> Dict[str, np.ndarray]:
    """
    Run multiple independent CEEMDAN realizations (different seeds) and compute:
      mean trend μ(t) and pointwise std σ(t).

    Useful for shaded confidence intervals: μ(t) ± 2σ(t).
    """
    tt, yy = select_window(t_ps, y, t_min, t_max)
    dt = check_uniform_dt(tt)

    trends = []
    for s in seeds:
        ce_params = CeemdanParams(
            trials=base_ce.trials,
            noise_width=base_ce.noise_width,
            seed=int(s),
            max_imfs=base_ce.max_imfs,
        )
        out = ceemdan_sg_trend(tt, yy, ce_params, base_sg, sel)
        trends.append(out["trend"])

    trends = np.asarray(trends, float)  # shape (n_seeds, n_time)
    mu = np.mean(trends, axis=0)
    sigma = np.std(trends, axis=0, ddof=1) if trends.shape[0] > 1 else np.zeros_like(mu)

    return {
        "t": tt,
        "mu": mu,
        "sigma": sigma,
        "dt": np.array([dt]),
        "trends": trends,
    }


# -----------------------------
# Plotting helpers
# -----------------------------
def plot_confidence_band(t: np.ndarray, y: np.ndarray, mu: np.ndarray, sigma: np.ndarray, title: str) -> None:
    plt.figure(figsize=(10, 4))
    plt.plot(t, y, color="k", lw=0.8, alpha=0.5, label="original")
    plt.plot(t, mu, lw=2.0, label="trend mean (CEEMDAN+SG)")
    plt.fill_between(t, mu - 2 * sigma, mu + 2 * sigma, alpha=0.25, label=r"$\mu \pm 2\sigma$")
    plt.xlabel("t (ps)")
    plt.ylabel("signal")
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_windowing_sweep_heatmap(df_win: pd.DataFrame) -> None:
    """
    Visualize windowing sweep results as heatmaps showing RMSD and max_abs_dev
    for different combinations of t_min and t_max.

    Parameters
    ----------
    df_win : pd.DataFrame
        Output from run_windowing_sweep() with columns:
        t_min, t_max, rmsd_overlap, max_abs_dev_overlap
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Pivot data for heatmap visualization
    pivot_rmsd = df_win.pivot(index='t_min', columns='t_max', values='rmsd_overlap')
    pivot_max = df_win.pivot(index='t_min', columns='t_max', values='max_abs_dev_overlap')

    # RMSD heatmap
    im1 = axes[0].imshow(pivot_rmsd.values, aspect='auto', cmap='viridis', origin='lower')
    axes[0].set_xticks(range(len(pivot_rmsd.columns)))
    axes[0].set_yticks(range(len(pivot_rmsd.index)))
    axes[0].set_xticklabels([f'{x:.1f}' for x in pivot_rmsd.columns])
    axes[0].set_yticklabels([f'{y:.1f}' for y in pivot_rmsd.index])
    axes[0].set_xlabel('t_max (ps)')
    axes[0].set_ylabel('t_min (ps)')
    axes[0].set_title('RMSD on Overlap Region')
    cbar1 = plt.colorbar(im1, ax=axes[0])
    cbar1.set_label('RMSD')

    # Add text annotations
    for i in range(len(pivot_rmsd.index)):
        for j in range(len(pivot_rmsd.columns)):
            val = pivot_rmsd.values[i, j]
            if not np.isnan(val):
                axes[0].text(j, i, f'{val:.3f}', ha='center', va='center',
                           color='white' if val > pivot_rmsd.values[~np.isnan(pivot_rmsd.values)].mean() else 'black',
                           fontsize=8)

    # Max absolute deviation heatmap
    im2 = axes[1].imshow(pivot_max.values, aspect='auto', cmap='plasma', origin='lower')
    axes[1].set_xticks(range(len(pivot_max.columns)))
    axes[1].set_yticks(range(len(pivot_max.index)))
    axes[1].set_xticklabels([f'{x:.1f}' for x in pivot_max.columns])
    axes[1].set_yticklabels([f'{y:.1f}' for y in pivot_max.index])
    axes[1].set_xlabel('t_max (ps)')
    axes[1].set_ylabel('t_min (ps)')
    axes[1].set_title('Max Absolute Deviation on Overlap Region')
    cbar2 = plt.colorbar(im2, ax=axes[1])
    cbar2.set_label('Max |Δ|')

    # Add text annotations
    for i in range(len(pivot_max.index)):
        for j in range(len(pivot_max.columns)):
            val = pivot_max.values[i, j]
            if not np.isnan(val):
                axes[1].text(j, i, f'{val:.3f}', ha='center', va='center',
                           color='white' if val > pivot_max.values[~np.isnan(pivot_max.values)].mean() else 'black',
                           fontsize=8)

    plt.tight_layout()
    plt.show()


def plot_windowing_sweep_summary(df_win: pd.DataFrame) -> None:
    """
    Create summary plots showing how windowing choices affect robustness metrics.

    Parameters
    ----------
    df_win : pd.DataFrame
        Output from run_windowing_sweep()
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. RMSD vs t_min for each t_max
    for t_max in sorted(df_win['t_max'].unique()):
        subset = df_win[df_win['t_max'] == t_max]
        axes[0, 0].plot(subset['t_min'], subset['rmsd_overlap'],
                       marker='o', label=f't_max={t_max:.1f}')
    axes[0, 0].set_xlabel('t_min (ps)')
    axes[0, 0].set_ylabel('RMSD (overlap)')
    axes[0, 0].set_title('RMSD Sensitivity to t_min')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # 2. RMSD vs t_max for each t_min
    for t_min in sorted(df_win['t_min'].unique()):
        subset = df_win[df_win['t_min'] == t_min]
        axes[0, 1].plot(subset['t_max'], subset['rmsd_overlap'],
                       marker='o', label=f't_min={t_min:.1f}')
    axes[0, 1].set_xlabel('t_max (ps)')
    axes[0, 1].set_ylabel('RMSD (overlap)')
    axes[0, 1].set_title('RMSD Sensitivity to t_max')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # 3. Max abs dev vs t_min for each t_max
    for t_max in sorted(df_win['t_max'].unique()):
        subset = df_win[df_win['t_max'] == t_max]
        axes[1, 0].plot(subset['t_min'], subset['max_abs_dev_overlap'],
                       marker='s', label=f't_max={t_max:.1f}')
    axes[1, 0].set_xlabel('t_min (ps)')
    axes[1, 0].set_ylabel('Max |Δ| (overlap)')
    axes[1, 0].set_title('Max Deviation Sensitivity to t_min')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)

    # 4. Max abs dev vs t_max for each t_min
    for t_min in sorted(df_win['t_min'].unique()):
        subset = df_win[df_win['t_min'] == t_min]
        axes[1, 1].plot(subset['t_max'], subset['max_abs_dev_overlap'],
                       marker='s', label=f't_min={t_min:.1f}')
    axes[1, 1].set_xlabel('t_max (ps)')
    axes[1, 1].set_ylabel('Max |Δ| (overlap)')
    axes[1, 1].set_title('Max Deviation Sensitivity to t_max')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


# -----------------------------
# Example usage (adapt to your dict)
# -----------------------------

if __name__ == "__main__":
    t_ps = dict9["t"]
    y = dict9["v2"]

    # --- baseline configuration ---
    base_ce = CeemdanParams(trials=200, noise_width=0.2, seed=0, max_imfs=None)
    base_sg = SgParams(window_length=401, polyorder=2)
    sel = ModeSelectionParams(target_period_ps=1.2, rel_tol=0.25)

    sweep = RobustnessSweep()

    # --- define your working window ---
    t_min = 2.67
    t_max = 8.5

    # 1) Hyperparameter sweep (noise_width x SG window)
    df_hp = run_hyperparameter_sweep(t_ps, y, t_min, t_max, base_ce, base_sg, sel, sweep)
    print(df_hp.sort_values("rmsd").head(10))

#%%
    # 2) Windowing sweep (t_min/t_max sensitivity, evaluated on overlap region)
    df_win = run_windowing_sweep(t_ps, y, base_ce, base_sg, sel, sweep, compare_on_overlap=(4.0, 8.0))
    print(df_win.sort_values("rmsd_overlap").head(10))

#%%
    plot_windowing_sweep_heatmap(df_win)
    plot_windowing_sweep_summary(df_win)





#%%
# # 3) Confidence interval band from 10 independent seeds
# band = compute_confidence_band(t_ps, y, t_min, t_max, base_ce, base_sg, sel, seeds=tuple(range(10)))
# plot_confidence_band(band["t"], y[(t_ps >= t_min) & (t_ps <= t_max)], band["mu"], band["sigma"],
#                      title="CEEMDAN+SG trend uncertainty (10 seeds)")