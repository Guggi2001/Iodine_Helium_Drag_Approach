from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import numpy as np

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Load Config and Data
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
from config_utils_local import config as C
from drag_function import io


dict9 = io.load_data(C.PATH9A)
dict18 = io.load_data(C.PATH18A)
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#                                                         Actual Run Tests
#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def rolling_rms(x, win):
    """Centered rolling RMS using convolution; returns array same length with NaN near edges."""
    x = np.asarray(x, float)
    win = int(win)
    if win < 3:
        raise ValueError("win must be >= 3")
    if win % 2 == 0:
        win += 1
    w = np.ones(win) / win
    x2 = x**2
    rms = np.sqrt(np.convolve(x2, w, mode="same"))
    # mark edges as NaN where window incomplete (optional but honest)
    half = win // 2
    rms[:half] = np.nan
    rms[-half:] = np.nan
    return rms

def find_t_star_stationary_residual(
    t, v, distance_case,
    dt=None,
    t_out=9.0,
    # mild trend removal:
    trend_win_ps=0.8,
    trend_poly=2,
    # rolling RMS window:
    rms_win_ps=0.5,
    # late-time band definition:
    late_window=(7.0, 8.8),
    k=2.5,                 # band width in MAD units
    sustain_ps=0.8,        # must remain in band for this long
    show=True
):
    """
    Determine transient end t* by stationarity of rolling RMS of residual.

    Returns:
      t_star (float or None), diagnostics dict
    """
    t = np.asarray(t, float)
    v = np.asarray(v, float)

    if dt is None:
        dt = float(np.mean(np.diff(t)))

    # --- 1) mild slow trend removal (SG) ---
    trend_win = int(round(trend_win_ps / dt))
    if trend_win % 2 == 0:
        trend_win += 1
    trend_win = max(trend_win, trend_poly + 3)  # ensure valid

    v_trend = savgol_filter(v, window_length=trend_win, polyorder=trend_poly, mode="interp")
    resid = v - v_trend

    # --- 2) rolling RMS of residual ---
    rms_win = int(round(rms_win_ps / dt))
    if rms_win % 2 == 0:
        rms_win += 1
    rms_win = max(rms_win, 5)
    rrms = rolling_rms(resid, rms_win)

    # --- 3) define late-time baseline band via robust stats ---
    lw0, lw1 = late_window
    lw1 = min(lw1, t_out)
    m_late = (t >= lw0) & (t <= lw1) & np.isfinite(rrms)
    if np.sum(m_late) < 20:
        raise ValueError("Late window too small; adjust late_window or rms/trend windows.")

    late_med = float(np.nanmedian(rrms[m_late]))
    late_mad = float(1.4826 * np.nanmedian(np.abs(rrms[m_late] - late_med)))
    band_lo = late_med - k * late_mad
    band_hi = late_med + k * late_mad

    # ensure nonnegative lower bound
    band_lo = max(0.0, band_lo)

    # --- 4) find earliest time where rrms stays within band for sustain_ps ---
    sustain_n = int(round(sustain_ps / dt))
    sustain_n = max(sustain_n, 5)

    # search only up to t_out
    m_search = (t <= t_out) & np.isfinite(rrms)
    idx = np.where(m_search)[0]

    t_star = None
    for i in idx:
        j = i + sustain_n
        if j >= len(t):
            break
        # require we stay within the band and within t_out
        if t[j-1] > t_out:
            break
        window_ok = np.all((rrms[i:j] >= band_lo) & (rrms[i:j] <= band_hi))
        if window_ok:
            t_star = float(t[i])
            break

    diagnostics = {
        "dt": dt,
        "trend_win": trend_win,
        "rms_win": rms_win,
        "late_med": late_med,
        "late_mad": late_mad,
        "band_lo": band_lo,
        "band_hi": band_hi,
        "t_star": t_star,
    }

    if show:
        fig, ax = plt.subplots(2, 1, figsize=(7, 4), sharex=True)
        fig.suptitle(f'{distance_case} Å case for top iodine atom', fontsize=16, fontweight='bold')
        ax[0].plot(t, v, 'k', lw=0.7, label="v_top raw")
        ax[0].plot(t, v_trend, 'r--', lw=1.2, label=f"trend (SG {trend_win_ps} ps)")
        ax[0].axvline(t_out, color='k', ls=':', label="t_out")
        if t_star is not None:
            ax[0].axvline(t_star, color='g', ls='--', label=f"t*={t_star:.2f} ps")
        ax[0].set_ylabel("v (Å/ps)")
        ax[0].grid(True)
        ax[0].legend(fontsize="small")

        ax[1].plot(t, rrms, lw=1.0, label=f"rolling RMS(resid) ({rms_win_ps} ps)")
        ax[1].axhline(band_lo, color='r', ls='--', lw=1.0, label="stationary band")
        ax[1].axhline(band_hi, color='r', ls='--', lw=1.0)
        ax[1].axvspan(lw0, lw1, color='gray', alpha=0.15, label="late window")
        ax[1].axvline(t_out, color='k', ls=':')
        if t_star is not None:
            ax[1].axvline(t_star, color='g', ls='--')
        ax[1].set_xlabel("t (ps)")
        ax[1].set_ylabel("RMS of residual")
        ax[1].grid(True)
        ax[1].legend(fontsize="small")

        plt.tight_layout()
        plt.show()

    return t_star, diagnostics

t18 = dict18["t"]
v18 = dict18["v2"]

t_star18, diag18 = find_t_star_stationary_residual(
    t18, v18, 18,
    t_out=8.5,
    trend_win_ps=0.8,
    rms_win_ps=0.5,
    late_window=(5, 8),
    k=2.5,
    sustain_ps=0.8,
    show=True
)

print("t*18 =", t_star18)
print(diag18)

t9 = dict9["t"]
v9 = dict9["v2"]
t_star9, diag9 = find_t_star_stationary_residual(
    t9, v9, 9,
    t_out=9.0,
    trend_win_ps=0.8,
    rms_win_ps=0.5,
    late_window=(5, 8.5),
    k=2.5,
    sustain_ps=0.8,
    show=True
)

print("t*9 =", t_star9)
print(diag9)