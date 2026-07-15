"""Plot the trusted experimental fragment curves (ihe_ked/IHe_KED_curves_n0..n4.csv).

Four panels, all five fragments (I+, I+He, I+He2, I+He3, I+He4) overlaid:

    2-D detector projection  vs velocity   |  2-D detector projection  vs energy
    3-D reconstruction P(v)  vs velocity   |  3-D reconstruction P(E)  vs energy

Each panel uses the signal column that is a proper density in its own x-axis
measure (signal_2d_Pv / signal_2d_PE / signal_3d_Pv / signal_3d_PE) -- see
ihe_ked/COLUMNS.md. Reminders encoded there and honored here:
  - against the 2-D columns the axis is the PROJECTED in-plane speed
    (v*sin(theta)) / apparent energy; against the 3-D columns it is the true
    speed / energy of the ion;
  - signal_2d_PE is elevated toward low E (finite center-fill plateau: the
    VMI image center is filled by fast ions projected inward) -- that rise
    is not slow-ion signal;
  - 2-D columns cover the full detector range; the 3-D columns carry NaN
    cuts (E < 0.01 eV, n=0 additionally E < 0.4 eV: Abel-center spike);
  - columns are independently envelope-normalized: compare shapes, never
    amplitudes across columns.

The dashed markers on the 3-D panels are the reference <E> per fragment
(IHe_KED_reference.csv), placed at <E> resp. v = sqrt(2<E>/m(n)).

Usage:  python experimental_reference_fragments.py
Writes: ../ihe_ked/experimental_reference_fragments.png
"""

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "ihe_ked")
OUT_PNG = os.path.join(DATA_DIR, "experimental_reference_fragments.png")

EKF = 1.03642697e-8          # u*(m/s)^2 -> eV
M_I, M_HE = 126.90, 4.0026
N_FRAGMENTS = range(5)
SMOOTH_BINS = 9              # display smoothing only (curves ship unsmoothed)

LABELS = {0: "I$^+$", 1: "I$^+$He"}
LABELS.update({n: f"I$^+$He$_{{{n}}}$" for n in range(2, 5)})


def load_curves(n):
    path = os.path.join(DATA_DIR, f"IHe_KED_curves_n{n}.csv")
    cols = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            for key, val in row.items():
                cols.setdefault(key, []).append(
                    float(val) if val not in ("", "NaN", "nan") else np.nan)
    return {k: np.asarray(v) for k, v in cols.items()}


def load_reference_means():
    path = os.path.join(DATA_DIR, "IHe_KED_reference.csv")
    means = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            means[int(row["n"])] = float(row["meanKE_eV"])
    return means


def movmean(y, k):
    out = np.full(len(y), np.nan)
    for i in range(len(y)):
        seg = y[max(0, i - k // 2):i + k // 2 + 1]
        if np.isfinite(seg).any():
            out[i] = np.nanmean(seg)
    out[~np.isfinite(y)] = np.nan   # do not extend into the NaN-cut regions
    return out


def main():
    means = load_reference_means()
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.6))
    panels = [
        (axs[0, 0], "v_mps", "signal_2d_Pv", "2-D detector projection vs velocity",
         "projected speed v$\\cdot$sin$\\theta$ / m s$^{-1}$"),
        (axs[0, 1], "E_eV", "signal_2d_PE", "2-D detector projection vs energy",
         "apparent energy / eV"),
        (axs[1, 0], "v_mps", "signal_3d_Pv", "3-D reconstruction P(v) vs velocity",
         "speed / m s$^{-1}$"),
        (axs[1, 1], "E_eV", "signal_3d_PE", "3-D reconstruction P(E) vs energy",
         "kinetic energy / eV"),
    ]
    cmap = plt.get_cmap("viridis")
    colors = [cmap(0.85 * n / (len(N_FRAGMENTS) - 1)) for n in N_FRAGMENTS]

    for n, color in zip(N_FRAGMENTS, colors):
        d = load_curves(n)
        mass_n = M_I + M_HE * n
        v_mean = np.sqrt(2.0 * means[n] / (mass_n * EKF))
        for ax, xcol, ycol, _, _ in panels:
            x, y = d[xcol], d[ycol]
            ax.plot(x, y, "-", color=color, lw=0.7, alpha=0.25)
            ax.plot(x, movmean(y, SMOOTH_BINS), "-", color=color, lw=1.6,
                    label=LABELS[n] if ax is axs[0, 0] else None)
        for ax, xpos in ((axs[1, 0], v_mean), (axs[1, 1], means[n])):
            ax.axvline(xpos, color=color, ls="--", lw=0.9, alpha=0.8)

    for ax, _, _, title, xlabel in panels:
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("signal, envelope-normalized")
        ax.set_ylim(0, 1.35)
        ax.set_xlim(0, 3100 if "m s" in xlabel else 6.3)
    axs[0, 0].legend(fontsize=9, loc="upper right")
    axs[1, 1].text(0.98, 0.86, "dashed: reference $\\langle E\\rangle$",
                   transform=axs[1, 1].transAxes, fontsize=8, ha="right")
    axs[1, 0].text(0.98, 0.86,
                   "dashed: $v(\\langle E\\rangle)$, not $\\langle v\\rangle$",
                   transform=axs[1, 0].transAxes, fontsize=8, ha="right")
    axs[0, 1].text(0.98, 0.86, "low-E elevation = center fill\n(projected fast ions), not slow ions",
                   transform=axs[0, 1].transAxes, fontsize=8, ha="right")
    fig.suptitle("Experimental reference fragments I$^+$He$_n$, n = 0..4 "
                 "(fragment scan 17.10.24, 600 mW series)", y=0.99)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=170)
    print(f"wrote {os.path.normpath(OUT_PNG)}")


if __name__ == "__main__":
    main()
