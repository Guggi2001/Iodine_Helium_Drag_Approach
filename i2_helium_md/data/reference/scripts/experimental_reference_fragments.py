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

A second figure holds only the two velocity panels (2-D projection | 3-D
P(v)), same drawing. Both are drawn at the thesis text width with the shared
``thesis_style``; the fragments are told apart by ``FRAGMENT_COLORS``, samples of
viridis.

Usage:  python data/reference/scripts/experimental_reference_fragments.py
Writes: data/figures/thesis_rest/experimental_reference_fragments.png (2 x 2)
        data/figures/thesis_rest/experimental_reference_fragments_velocity.png (1 x 2)
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.post_processing.thesis_style import (  # noqa: E402
    FIG_RC,
    FONT,
    TEXTWIDTH_IN,
    panel_letter,
    style_axis,
)

DATA_DIR: Path = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_rest"
OUT_PATH: Path = OUT_DIR / "experimental_reference_fragments.png"
OUT_PATH_VELOCITY: Path = OUT_DIR / "experimental_reference_fragments_velocity.png"

EKF = 1.03642697e-8          # u*(m/s)^2 -> eV
M_I, M_HE = 126.90, 4.0026
N_FRAGMENTS = range(5)
SMOOTH_BINS = 9              # display smoothing only (curves ship unsmoothed)
FRAGMENT_COLORS = tuple(plt.get_cmap("viridis")(np.linspace(0.0, 0.85, len(N_FRAGMENTS))))

FIGSIZE_OVERVIEW_IN: tuple[float, float] = (TEXTWIDTH_IN, 5.0)
FIGSIZE_VELOCITY_IN: tuple[float, float] = (TEXTWIDTH_IN, 2.9)
DPI: int = 200
RAW_LINE: dict = dict(lw=0.4, alpha=0.3)
SMOOTH_LINE: dict = dict(lw=1.2)
MEAN_LINE: dict = dict(ls="--", lw=0.8)
NOTE_BOX: dict = dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5)
NOTE_LOC: tuple[float, float] = (0.97, 0.95)   # axes fraction, top-right anchor
YLABEL: str = "Normalized signal"
Y_MAX: float = 1.35
V_MAX_MPS: float = 3100.0
E_MAX_EV: float = 6.3
Y_TICK_STEP: float = 0.25
V_TICK_STEP_MPS: float = 1000.0

LABELS = {0: "I$^+$", 1: "I$^+$He"}
LABELS.update({n: f"I$^+$He$_{{{n}}}$" for n in range(2, 5)})


def load_curves(n: int) -> dict[str, np.ndarray]:
    """Read the columns of ``IHe_KED_curves_n{n}.csv``; empty and NaN cells become NaN."""
    cols: dict[str, list[float]] = {}
    with open(DATA_DIR / f"IHe_KED_curves_n{n}.csv", newline="") as f:
        for row in csv.DictReader(f):
            for key, val in row.items():
                cols.setdefault(key, []).append(
                    float(val) if val not in ("", "NaN", "nan") else np.nan)
    return {k: np.asarray(v) for k, v in cols.items()}


def load_reference_means() -> dict[int, float]:
    """Reference mean kinetic energy per fragment n [eV]."""
    means = {}
    with open(DATA_DIR / "IHe_KED_reference.csv", newline="") as f:
        for row in csv.DictReader(f):
            means[int(row["n"])] = float(row["meanKE_eV"])
    return means


def movmean(y: np.ndarray, k: int) -> np.ndarray:
    """Centered moving mean over ``k`` bins that does not fill the NaN cuts of ``y``."""
    out = np.full(len(y), np.nan)
    for i in range(len(y)):
        seg = y[max(0, i - k // 2):i + k // 2 + 1]
        if np.isfinite(seg).any():
            out[i] = np.nanmean(seg)
    out[~np.isfinite(y)] = np.nan   # do not extend into the NaN-cut regions
    return out


# panel key -> (x column, signal column, title, x label, dashed-marker kind, note)
PANELS = {
    "2d_v": ("v_mps", "signal_2d_Pv", "2-D projection $P(v)$",
             r"Projected speed $v\,\sin\theta$ / m s$^{-1}$", None, None),
    "2d_E": ("E_eV", "signal_2d_PE", "2-D projection $P(E)$",
             "Apparent energy / eV", None,
             "Low-$E$ rise: center fill\n(projected fast ions)"),
    "3d_v": ("v_mps", "signal_3d_Pv", "3-D reconstruction $P(v)$",
             r"Speed $v$ / m s$^{-1}$", "v",
             r"Dashed: $v(\langle E\rangle)$, not $\langle v\rangle$"),
    "3d_E": ("E_eV", "signal_3d_PE", "3-D reconstruction $P(E)$",
             r"Kinetic energy $E$ / eV", "E", r"Dashed: reference $\langle E\rangle$"),
}


def draw_panel(
    ax: plt.Axes,
    key: str,
    curves: dict[int, dict[str, np.ndarray]],
    means: dict[int, float],
    letter: str,
) -> None:
    """All fragments on one panel: raw (faint) + display-smoothed curves and,
    on the 3-D panels, the reference <E> marker (at <E> or v(<E>))."""
    xcol, ycol, title, xlabel, marker, note = PANELS[key]
    for n, color in zip(N_FRAGMENTS, FRAGMENT_COLORS):
        x, y = curves[n][xcol], curves[n][ycol]
        ax.plot(x, y, color=color, **RAW_LINE)
        ax.plot(x, movmean(y, SMOOTH_BINS), color=color, label=LABELS[n], **SMOOTH_LINE)
        if marker == "v":
            ax.axvline(np.sqrt(2.0 * means[n] / ((M_I + M_HE * n) * EKF)),
                       color=color, **MEAN_LINE)
        elif marker == "E":
            ax.axvline(means[n], color=color, **MEAN_LINE)
    ax.set_title(title, fontsize=FONT["title_fontsize"])
    ax.set(xlabel=xlabel, ylabel=YLABEL, ylim=(0.0, Y_MAX),
           xlim=(0.0, V_MAX_MPS if xcol == "v_mps" else E_MAX_EV))
    ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_STEP))
    if xcol == "v_mps":
        ax.xaxis.set_major_locator(MultipleLocator(V_TICK_STEP_MPS))
        ax.xaxis.set_minor_locator(MultipleLocator(V_TICK_STEP_MPS / 2))
    if note:
        ax.text(*NOTE_LOC, note, transform=ax.transAxes, ha="right", va="top",
                fontsize=FONT["annot_fontsize"], bbox=NOTE_BOX)
    style_axis(ax)
    panel_letter(ax, letter)


def build_figure(
    panel_keys: tuple[str, ...],
    shape: tuple[int, int],
    figsize_in: tuple[float, float],
    curves: dict[int, dict[str, np.ndarray]],
    means: dict[int, float],
) -> Figure:
    """Draw ``panel_keys`` on a ``shape`` grid with one legend above the panels.

    Parameters
    ----------
    panel_keys
        Keys of ``PANELS``, row by row.
    shape
        Rows and columns of the grid.
    figsize_in
        Figure size in inches.
    curves
        Columns of each fragment n, as from :func:`load_curves`.
    means
        Reference mean kinetic energy per fragment [eV].

    Returns
    -------
    Figure
        The figure; the caller saves and closes it.
    """
    fig, axes = plt.subplots(*shape, figsize=figsize_in, constrained_layout=True)
    for i, (ax, key) in enumerate(zip(np.ravel(axes), panel_keys)):
        draw_panel(ax, key, curves, means, chr(ord("a") + i))
    fig.legend(*np.ravel(axes)[0].get_legend_handles_labels(), loc="outside upper center",
               ncols=len(N_FRAGMENTS), frameon=False, fontsize=FONT["legend_fontsize"])
    return fig


def main() -> None:
    """Load the references and write both figures to ``OUT_DIR``."""
    means = load_reference_means()
    curves = {n: load_curves(n) for n in N_FRAGMENTS}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    figures = (
        (("2d_v", "2d_E", "3d_v", "3d_E"), (2, 2), FIGSIZE_OVERVIEW_IN, OUT_PATH),
        (("2d_v", "3d_v"), (1, 2), FIGSIZE_VELOCITY_IN, OUT_PATH_VELOCITY),
    )
    with plt.rc_context(FIG_RC):
        for panel_keys, shape, figsize_in, path in figures:
            fig = build_figure(panel_keys, shape, figsize_in, curves, means)
            fig.savefig(path, dpi=DPI)
            plt.close(fig)
            print(f"wrote {path}")


if __name__ == "__main__":
    main()
