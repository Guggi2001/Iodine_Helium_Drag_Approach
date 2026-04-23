#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Load Config and Data
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#%%
from config_utils_local import config as C
from drag_function import io


dict9 = io.load_data(C.PATH9A)
dict18 = io.load_data(C.PATH18A)


import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Optional, Tuple, List

from PyEMD import CEEMDAN
from scipy.signal import find_peaks


# -----------------------------
# Period estimation utilities
# -----------------------------
def estimate_period_zero_crossings(t_ps: np.ndarray, x: np.ndarray) -> Optional[float]:
    """
    Estimate average period from zero-crossings.
    For an oscillatory IMF, consecutive zero-crossings are ~T/2 apart.

    Returns:
      period in ps (float) or None if insufficient crossings.
    """
    x = np.asarray(x, float)
    t_ps = np.asarray(t_ps, float)

    # Find sign changes
    s = np.sign(x)
    s[s == 0] = 1  # treat exact zeros as positive
    zc_idx = np.where(np.diff(s) != 0)[0]
    if zc_idx.size < 4:
        return None

    # Linear interpolation for better crossing times
    t_zc = []
    for i in zc_idx:
        x0, x1 = x[i], x[i + 1]
        t0, t1 = t_ps[i], t_ps[i + 1]
        # x(t) assumed linear between samples
        if (x1 - x0) == 0:
            continue
        frac = -x0 / (x1 - x0)
        t_cross = t0 + frac * (t1 - t0)
        t_zc.append(t_cross)

    t_zc = np.array(t_zc, float)
    if t_zc.size < 4:
        return None

    # Half-period estimates from consecutive crossings
    half_periods = np.diff(t_zc)
    # Robust average: median is less sensitive to outliers
    half_period = np.median(half_periods)
    return 2.0 * half_period


def estimate_period_peaks(t_ps: np.ndarray, x: np.ndarray) -> Optional[float]:
    """
    Estimate average period from peak-to-peak distance.
    Works well if the IMF has clear maxima.
    """
    x = np.asarray(x, float)
    t_ps = np.asarray(t_ps, float)

    peaks, _ = find_peaks(x)
    if peaks.size < 3:
        return None

    dt_peaks = np.diff(t_ps[peaks])
    return float(np.median(dt_peaks))


# -----------------------------
# CEEMDAN workflow
# -----------------------------
@dataclass
class CeemdanConfig:
    trials: int = 50
    noise_width: float = 0.15
    random_seed: int = 0
    max_imfs: Optional[int] = None  # None = let CEEMDAN decide


def ceemdan_decompose(t_ps: np.ndarray, y: np.ndarray, cfg: CeemdanConfig) -> np.ndarray:
    """
    Compute CEEMDAN IMFs.

    Note: CEEMDAN expects uniform sampling for meaningful period estimates.
    """
    t_ps = np.asarray(t_ps, float)
    y = np.asarray(y, float)

    dt = np.mean(np.diff(t_ps))
    if not np.allclose(np.diff(t_ps), dt, rtol=1e-3, atol=1e-12):
        raise ValueError("Time array is not uniformly sampled; resample before CEEMDAN.")

    ce = CEEMDAN(trials=cfg.trials, noise_width=cfg.noise_width)
    ce.random_seed = cfg.random_seed

    # Optional: limit number of IMFs
    if cfg.max_imfs is not None:
        ce.MAX_IMF = int(cfg.max_imfs)

    imfs = ce.ceemdan(y, t_ps)
    return imfs


def identify_imfs_by_period(
    t_ps: np.ndarray,
    imfs: np.ndarray,
    target_period_ps: float = 1.2,
    rel_tol: float = 0.20,
    method: str = "zero_crossings",
) -> Tuple[List[int], List[Optional[float]]]:
    """
    Identify which IMFs match a target period within relative tolerance.

    rel_tol=0.20 means accept periods in [0.8*target, 1.2*target].

    Returns:
      indices: list of IMF indices to remove
      periods: list of estimated periods for all IMFs
    """
    periods = []
    for k in range(imfs.shape[0]):
        x = imfs[k]
        if method == "zero_crossings":
            T = estimate_period_zero_crossings(t_ps, x)
        elif method == "peaks":
            T = estimate_period_peaks(t_ps, x)
        else:
            raise ValueError("method must be 'zero_crossings' or 'peaks'.")
        periods.append(T)

    lo = (1.0 - rel_tol) * target_period_ps
    hi = (1.0 + rel_tol) * target_period_ps

    idx = []
    for k, T in enumerate(periods):
        if T is None:
            continue
        if lo <= T <= hi:
            idx.append(k)

    return idx, periods


def reconstruct_without_imfs(imfs: np.ndarray, drop_indices: List[int]) -> np.ndarray:
    """
    Reconstruct signal by summing all IMFs except those in drop_indices.
    """
    keep = [k for k in range(imfs.shape[0]) if k not in set(drop_indices)]
    if not keep:
        raise ValueError("All IMFs were dropped; cannot reconstruct.")
    return np.sum(imfs[keep, :], axis=0)

from scipy.signal import savgol_filter
def sg_smooth_v(t, v, window_length, polyorder=3):
    t = np.asarray(t, float)
    v = np.asarray(v, float)

    # Ensure window_length is a plain int for savgol_filter and for parity checks
    window_length = int(window_length)

    dt = np.mean(np.diff(t))
    if not np.allclose(np.diff(t), dt, rtol=1e-3, atol=1e-12):
        raise ValueError("t is not uniformly spaced enough for SG derivative.")

    if window_length % 2 == 0:
        window_length += 1
    if window_length >= len(v):
        raise ValueError("window_length must be smaller than data length.")
    if polyorder >= window_length:
        raise ValueError("polyorder must be < window_length.")

    v_s = savgol_filter(v, window_length=window_length, polyorder=polyorder, deriv=0, delta=dt, mode="interp")
    resid = v - v_s
    return v_s, resid, dt

# -----------------------------
# High-quality visualization
# -----------------------------
def plot_ceemdan_result(
    t_ps: np.ndarray,
    y: np.ndarray,
    imfs: np.ndarray,
    drop_indices: List[int],
    periods: List[Optional[float]],
    target_period_ps: float = 1.2,
):
    cleaned = reconstruct_without_imfs(imfs, drop_indices)
    extracted_mode = np.sum(imfs[drop_indices, :], axis=0) if drop_indices else np.zeros_like(y)

    plt.figure(figsize=(10, 5))
    plt.plot(t_ps, y, lw=1.0, label="original signal")
    plt.plot(t_ps, cleaned, lw=2.0, label="reconstructed baseline (IMFs excluding target)")
    if drop_indices:
        plt.plot(t_ps, extracted_mode, lw=1.0, label=f"isolated ~{target_period_ps:.2f} ps mode (sum of selected IMF(s))")
    plt.xlabel("t (ps)")
    plt.ylabel("signal (arb.)")
    plt.title("CEEMDAN filtering: original vs isolated mode vs reconstructed baseline")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Optional: show IMF periods
    print("IMF period estimates (ps):")
    for k, T in enumerate(periods):
        tag = " <== dropped" if k in drop_indices else ""
        if T is None:
            print(f"  IMF {k:02d}: period ~ None{tag}")
        else:
            print(f"  IMF {k:02d}: period ~ {T:.3f}{tag}")





# -----------------------------
# Example run (dummy or real)
# -----------------------------
if __name__ == "__main__":
    t = dict9["t"]
    v = dict9["v2"]
    mask = (t >= 2.67) & (t <= 8.5)
    t_ps = t[mask]
    y = v[mask]

    cfg = CeemdanConfig(trials=80, noise_width=0.2, random_seed=0, max_imfs=None)
    imfs = ceemdan_decompose(t_ps, y, cfg)

    drop_idx, periods = identify_imfs_by_period(
        t_ps, imfs, target_period_ps=1.2, rel_tol=0.20, method="zero_crossings"
    )

    plot_ceemdan_result(t_ps, y, imfs, drop_idx, periods, target_period_ps=1.2)

    cleaned = reconstruct_without_imfs(imfs, drop_idx)
    #%%
    wls = [1001, 2501,  3901]


    plt.figure(figsize=(10, 5))
    plt.plot(t_ps, y, lw=1.0, label="original signal")
    plt.plot(t_ps, cleaned, lw=2.0, label="IMF filtering")
    for wl in wls:
        cleaned_SG, _, _ = sg_smooth_v(t_ps, cleaned, window_length=wl, polyorder=1)
        plt.plot(t_ps, cleaned_SG, lw=2.0, ls = '--', label="IMF + SG filtering (wl={})".format(wl))
    plt.xlabel("t (ps)")
    plt.ylabel("signal (arb.)")
    plt.title(f"Originial vs IMF filtering vs IMF + SG filtering with various window lengths")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    #%%
    # Export data to CSV
    import pandas as pd

    wl = 3901
    cleaned_SG, _, _ = sg_smooth_v(t_ps, cleaned, window_length=wl, polyorder=1)

    # Create DataFrame with time, cleaned_SG, and IMF cleaned data
    export_data = pd.DataFrame({
        'time': t_ps,
        'cleaned_SG': cleaned_SG,
        'IMF_cleaned': cleaned
    })

    # Export to CSV
    export_data.to_csv('cleaned_data.csv', index=False)
    print(f"Data exported to cleaned_data.csv")

