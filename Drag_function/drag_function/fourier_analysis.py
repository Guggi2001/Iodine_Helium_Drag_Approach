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

from scipy.signal import savgol_filter, butter, sosfiltfilt, filtfilt


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



def plot_test_smoothing(dictionary, atom, z_component = False,  polyorder = 3):
    v_list = []
    resid_list = []
    t_list = []
    wl = []
    if z_component:
        y = dictionary[f"v{atom}_z"]
    else:
        y = dictionary[f"v{atom}"]

    for i in range(500, 8500, 2000):
        v_s, resid, dt = sg_smooth_v(dictionary["t"], y, window_length=i, polyorder = polyorder)
        v_list.append(v_s)
        wl.append(i)
        resid_list.append(resid)
        t_list.append(dt)
        t_full = np.asarray(dictionary["t"], float)
        v_full = np.asarray(y, float)
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))

        # Left: overlay raw v1 with smoothed signals
        axs[0].plot(t_full, v_full, color='k', lw=0.8, label=f'v{atom} (raw)')
        axs[0].plot(t_full, v_s, 'b--', lw=1.2, alpha=0.9, label=f'v{atom}_s wl={i}')
        axs[0].set_xlabel('t (ps)')
        axs[0].set_ylabel('velocity (Å/ps)')
        axs[0].set_title(f'v{atom}: raw vs smoothed')
        axs[0].legend(fontsize='small')
        axs[0].grid(True)

        # Right: residuals (v1 - v1_s) for each smoothing; include RMS in legend
        rms = np.sqrt(np.mean(resid ** 2))
        axs[1].plot(t_full, resid, lw=1.0, alpha=0.9, label=f'wl={i}  RMS={rms:.3e}')
        axs[1].axhline(0.0, color='k', lw=0.8, alpha=0.6)
        axs[1].set_xlabel('t (ps)')
        axs[1].set_ylabel('residual (Å/ps)')
        axs[1].set_title(f'Residuals: v{atom} - v{atom}_s')
        axs[1].legend(fontsize='small')
        axs[1].grid(True)

        plt.tight_layout()
        plt.show()


    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    t_full = np.asarray(dictionary["t"], float)
    axs[0].plot(t_full, v_full, color='k', lw=0.7, alpha=0.95, label='v1 (raw)')
    for i in range(len(v_list)):
        axs[0].plot(t_full, v_list[i], ls = '--', lw=1.2, label=f'v{atom}_s wl={wl[i]}')
        axs[0].set_xlabel('t (ps)')
        axs[0].set_ylabel('velocity (Å/ps)')
        axs[0].set_title('v1: smoothed')
        axs[0].legend(fontsize='small')
        axs[0].grid(True)

        rms = np.sqrt(np.mean(resid_list[i] ** 2))
        axs[1].plot(t_full, resid_list[i], lw=1.0, alpha=0.9, label=f'wl={wl[i]}  RMS={rms:.3e}')
        axs[1].axhline(0.0, color='k', lw=0.8, alpha=0.6)
        axs[1].set_xlabel('t (ps)')
        axs[1].set_ylabel('residual (Å/ps)')
        axs[1].set_title('Residuals: v1 - v1_s')
        axs[1].legend(fontsize='small')
        axs[1].grid(True)
    plt.tight_layout()
    plt.show()

#plot_test_smoothing(dict18, 2, z_component = False, polyorder = 3)



from scipy.signal import detrend as scipy_detrend, find_peaks


def dominant_periods_fft(
        t, y, t_min=None, t_max=None, detrend='linear',
        n_peaks=3, fmin=0.05, fmax=None,
        peak_min_prominence_ratio=0.05, peak_min_distance_hz=0.05
):
    t, y = np.asarray(t, float), np.asarray(y, float)

    # 1. Time window
    mask = np.ones_like(t, dtype=bool)
    if t_min is not None: mask &= (t >= t_min)

    if t_max is not None: mask &= (t <= t_max)
    tt, yy = t[mask], y[mask]

    # 2. Uniform spacing check
    dt = np.mean(np.diff(tt))
    if not np.allclose(np.diff(tt), dt, rtol=1e-3, atol=1e-12):
        raise ValueError("t is not uniformly spaced enough for FFT.")

    # 3. Detrending
    if detrend == 'mean':
        yy = yy - np.mean(yy)
    elif detrend == 'linear':
        yy = scipy_detrend(yy, type='linear')
    elif detrend != 'none':
        raise ValueError("detrend must be 'none', 'mean', or 'linear'.")

    # 4. Windowing and FFT
    yyw = yy * np.hanning(len(yy))
    Y = np.fft.rfft(yyw)
    freqs = np.fft.rfftfreq(len(yyw), d=dt)
    power = np.abs(Y) ** 2

    # 5. Frequency band filtering
    fmax = fmax if fmax is not None else (0.5 / dt)
    band = (freqs >= fmin) & (freqs <= fmax)
    freqs_b, power_b = freqs[band], power[band]

    if len(freqs_b) < 5:
        return [], freqs_b, power_b

    # 6. Peak finding
    df = freqs_b[1] - freqs_b[0]
    min_dist_bins = max(1, int(np.ceil(peak_min_distance_hz / df)))
    prom_abs = peak_min_prominence_ratio * np.max(power_b)

    peak_idx, props = find_peaks(power_b, prominence=prom_abs, distance=min_dist_bins)

    if peak_idx.size == 0:
        idx = int(np.argmax(power_b))
        return [(float(freqs_b[idx]), float(1.0 / freqs_b[idx]), float(power_b[idx]))], freqs_b, power_b

    # 7. Sorting and output
    prominences = props.get("prominences", np.zeros_like(peak_idx, dtype=float))
    order = np.lexsort((power_b[peak_idx], prominences))[::-1]
    peak_idx_sorted = peak_idx[order]

    peaks = [(float(freqs_b[i]), float(1.0 / freqs_b[i]), float(power_b[i])) for i in peak_idx_sorted[:n_peaks]]

    return peaks, freqs_b, power_b



def plot_spectrum(freqs, power, fmax=5.0, logy=False, t_min = None, t_max = None, title="FFT Power Spectrum"):
    title += f" (fmin={freqs[0]:.2f} 1/ps, fmax={fmax:.2f} 1/ps)"
    if t_max != None and t_min != None:
        title += f"\n (t={t_min:.1f}-{t_max:.1f} ps)"
    m = freqs <= fmax
    plt.figure(figsize=(8,4))
    plt.plot(freqs[m], power[m], lw=1.0)
    if logy:
        plt.yscale("log")
    plt.xlabel("frequency (1/ps)")
    plt.ylabel("power (arb.)")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_spectrum_stack(t, y, t_min_list = (2.0, 3.0, 3.5, 4.0), t_max = 8.5, detrend="linear",
    n_peaks=3, fmin=0.05, fmax_plot=5.0, logy=False, figsize=(9, 10)):
    """
    Make a 4x1 plot: power spectra for different t_min values stacked vertically.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)

    n = len(t_min_list)
    fig, axes = plt.subplots(n, 1, figsize=figsize, sharex=True)

    # If only one axis, wrap to list for uniform handling
    if n == 1:
        axes = [axes]

    for ax, t_min in zip(axes, t_min_list):
        peaks, freqs_b, power_b = dominant_periods_fft(
            t, y,
            t_min=t_min,
            t_max=t_max,
            detrend=detrend,
            n_peaks=n_peaks,
            fmin=fmin,
            fmax=None
        )

        m = freqs_b <= fmax_plot
        ax.plot(freqs_b[m], power_b[m], lw=1.0)

        if logy:
            ax.set_yscale("log")

        # Mark peaks
        for f, T, _p in peaks:
            if f <= fmax_plot:
                ax.axvline(f, ls="--", lw=1.0, alpha=0.8)
                ax.text(
                    f, 0.95, f"T≈{T:.2f} ps",
                    transform=ax.get_xaxis_transform(),
                    rotation=90, va="top", ha="right", fontsize=9
                )

        ax.set_title(f"Power spectrum (t_min={t_min} ps) | peaks: "
                     + ", ".join([f"{T:.2f} ps" for _, T, _ in peaks]))
        ax.grid(True)

    axes[-1].set_xlabel("frequency (1/ps)")
    fig.supylabel("power (arb.)")

    plt.tight_layout()
    plt.show()

# t_min = 2
# t_max = 9
# f_min = 0.05
# peaks, freqs_b, power_b = dominant_periods_fft(dict9["t"], dict9["v2"], t_min= t_min, t_max = t_max, detrend='linear', n_peaks=3, fmin=f_min)
# print(peaks)  # list of (freq, period, power)
# plot_spectrum(freqs_b, power_b, t_min= t_min, t_max = t_max, fmax=5,  logy=False)
# plot_spectrum_stack(dict9["t"], dict9["v2"], t_max = t_max, logy = False)

def oscillation_summary(data, t_star, t_out=9.0, fmin=0.05, fmax_plot=5.0, y_log = False, n_peaks=3, title=""):
    t = np.asarray(data["t"], float)
    v = np.asarray(data["v2"], float)

    peaks, freqs_b, power_b = dominant_periods_fft(
        t, v, t_min=t_star, t_max=t_out,
        detrend="linear", n_peaks=n_peaks, fmin=fmin
    )

    # plot (zoomed)
    m = freqs_b <= fmax_plot
    plt.figure(figsize=(7,4))
    plt.plot(freqs_b[m], power_b[m], lw=1.0)
    if y_log:
        plt.yscale("log")
    for f, T, _ in peaks:
        if f <= fmax_plot:
            plt.axvline(f, ls="--", lw=1.0)
            plt.text(f, np.max(power_b[m])*0.8, f"T≈{T:.2f} ps", rotation=90, va="top", fontsize=9)
    plt.xlabel("frequency (1/ps)")
    plt.ylabel("power (arb.)")
    plt.title(title + f" | window [{t_star},{t_out}] ps")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return peaks

# peaks_9  = oscillation_summary(dict9,  t_star=2.67, t_out=9.0, fmax_plot= 3,  title="Velocity Spectrum: Top Iodine (9 Å)", y_log=False)
# peaks_18 = oscillation_summary(dict18, t_star=4.90, t_out=8.5, fmax_plot= 3, title="Velocity Spectrum: Top Iodine (18 Å)", y_log=False)
# print("9Å peaks:", peaks_9[:3])
# print("18Å peaks:", peaks_18[:3])


def extract_smoothed_region(t,v, case_distance, t_start, t_end,  window_length, polyorder=3, enable_plot=False):
    """
    Smooths a specific time range and optionally plots the result with boundaries.
    """
    t = np.asarray(t)
    v = np.asarray(v)

    # 1. Identify indices
    mask = (t >= t_start) & (t <= t_end)
    indices = np.where(mask)[0]

    if len(indices) == 0:
        raise ValueError(f"No data found in range {t_start} to {t_end}")

    idx_start, idx_end = indices[0], indices[-1]

    # 2. Extract and Smooth
    t_segment = t[idx_start: idx_end + 1]
    v_segment = v[idx_start: idx_end + 1]
    v_s_segment, resid, dt = sg_smooth_v(t_segment, v_segment, window_length, polyorder)


    # 4. Optional Diagnostic Plot
    if enable_plot:
        plt.figure(figsize=(10, 4))
        plt.plot(t_segment, v_segment, color='gray', alpha=0.4, label='Original velocity')
        plt.plot(t_segment, v_s_segment, color='#0072BD', lw=1.5, label='Smoothed velocity')

        plt.title(f'{case_distance} Å Smoothing Verification: {t_start}s to {t_end}s')
        plt.xlabel('t / ps')
        plt.ylabel('v / Å/ps')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

    return t_segment, v_s_segment

def smooth_range_only(t, v, case_distance, t_start, t_end, window_length, polyorder=3, enable_plot=False):
    """
    Smooths a specific time range and optionally plots the result with boundaries.
    """
    t = np.asarray(t)
    v = np.asarray(v)

    # 1. Identify indices
    mask = (t >= t_start) & (t <= t_end)
    indices = np.where(mask)[0]

    if len(indices) == 0:
        raise ValueError(f"No data found in range {t_start} to {t_end}")

    idx_start, idx_end = indices[0], indices[-1]

    # 2. Extract and Smooth
    t_segment = t[idx_start: idx_end + 1]
    v_segment = v[idx_start: idx_end + 1]
    v_s_segment, resid, dt = sg_smooth_v(t_segment, v_segment, window_length, polyorder)

    # 3. Stitch
    v_full_smoothed = np.concatenate([
        v[:idx_start],
        v_s_segment,
        v[idx_end + 1:]
    ])

    # 4. Optional Diagnostic Plot
    if enable_plot:
        plt.figure(figsize=(10, 4))
        plt.plot(t, v, color='gray', alpha=0.4, label='Original (Full)')
        plt.plot(t, v_full_smoothed, color='#0072BD', lw=1.5, label='Stitched Result')

        # Mark boundaries
        plt.axvline(t_start, color='red', linestyle='--', alpha=0.7, label='Start Smoothing')
        plt.axvline(t_end, color='red', linestyle='--', alpha=0.7, label='End Smoothing')

        plt.title(f'{case_distance} Å Smoothing Verification: {t_start}s to {t_end}s')
        plt.xlabel('t / ps')
        plt.ylabel('v / Å/ps')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

    return v_full_smoothed

#v_part_smoothed = smooth_range_only(dict9["t"], dict9["v2"], 9,  t_start=2.67, t_end=8.5, window_length=2501, polyorder=1, enable_plot=True)
v_part_smoothed = smooth_range_only(dict18["t"], dict18["v2"], 18,  t_start=4.543, t_end=8, window_length=3401, polyorder=1, enable_plot=True)

t_seg, v_s_seg = extract_smoothed_region(dict9["t"], dict9["v2"], 9,  t_start=2.67, t_end=8.5, window_length=4901, polyorder=1, enable_plot=True)
t_seg, v_s_seg = extract_smoothed_region(dict18["t"], dict18["v2"], 18,  t_start=4.543, t_end=8, window_length=3401, polyorder=1, enable_plot=True)

t18 = dict18["t"]
R18 = dict18["R"]
v18 = dict18["v2"]
v18_1 = dict18['v1']
mask = (t18 >= 4.543) & (t18 <= 8)
t_18_w = t18[mask]
R_18_w = R18[mask]
v_18_w = v18[mask]
v_18_1_w = v18_1[mask]

a = 3
from scipy.integrate import cumulative_trapezoid
def reconstruct_R_from_v(v_w, t_w, R0):
    """
    Reconstruct R from smoothed velocity data using Savitzky-Golay filter.
    """
    return (R0 + cumulative_trapezoid(v_w, t_w, initial=0))


plt.figure(figsize=(7, 4))
plt.plot(t_seg, v_18_w,  alpha=0.4, label='Given v Iodine 2')
plt.plot(t_seg, v_18_1_w, alpha=0.4, label='Given v Iodine 1')
plt.plot(t_seg, v_s_seg, label='Smoothed v (SG wl=3401)')
plt.title('Given v and Reconstructed v')
plt.xlabel('Time (ps)')
plt.ylabel('Velocity (Å/ps)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()


wls = [1001, 2501,  3401]
plt.figure(figsize=(7, 4))
plt.plot(t_seg, v_18_w,  alpha=0.4, label='Original signal')
for wl in wls:
    v_s_seg, _, _ = sg_smooth_v(t_seg, v_18_w, window_length=wl, polyorder=1)
    plt.plot(t_seg, v_s_seg, lw=2.0, ls = '--', label='SG filtering (wl={})'.format(wl))
plt.title('Given v and Reconstructed v')
plt.xlabel('Time (ps)')
plt.ylabel('Velocity (Å/ps)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()


plt.figure(figsize=(7, 4))
plt.plot(t_seg, R_18_w, label='Given R')
plt.plot(t_seg, reconstruct_R_from_v(2*v_s_seg, t_seg, R0=R_18_w[0]), label='Reconstructed R')
plt.title('Given R and Reconstructed R in relevant time window for 9Å case')
plt.xlabel('Time (ps)')
plt.ylabel('Distance Å')
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()
