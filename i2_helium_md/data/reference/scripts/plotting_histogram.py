"""Abundance figures: measured I+He_n size distribution, without and with bare I+.

``I+He_n_Abundance``: bare I+ (n = 0) is excluded and the remaining bins are
reweighted to 100 %, linear axis. ``I+He_n_Abundance_with_bare``: all bins,
n = 0 included, in percent of all ions, logarithmic axis.

Drawn with the shared ``thesis_style`` at half the thesis text width, so the font
sizes are the printed sizes in a half-width subfigure; bars as the experiment
bars of the results histograms (``tier2_b2_thesis_figure``).

Reads ``data/reference/integrated_i_he_abundance.csv``.

Usage:  python data/reference/scripts/plotting_histogram.py
Writes: data/figures/thesis_rest/I+He_n_Abundance.png
        data/figures/thesis_rest/I+He_n_Abundance_with_bare.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess import load_he_abundance_reference  # noqa: E402
from i2_helium_md.postprocess.abundance_loader import HeAbundanceReference  # noqa: E402
from scripts.post_processing.thesis_style import (  # noqa: E402
    EXP_COLOR,
    FIG_RC,
    TEXTWIDTH_IN,
    style_axis,
)

ABUNDANCE_CSV: Path = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_rest"
OUT_PATH: Path = OUT_DIR / "I+He_n_Abundance.png"
OUT_PATH_WITH_BARE: Path = OUT_DIR / "I+He_n_Abundance_with_bare.png"

FIGSIZE_IN: tuple[float, float] = (TEXTWIDTH_IN / 2, 2.0)
DPI: int = 200
BAR_WIDTH: float = 0.8
BAR_EDGE_WIDTH: float = 0.4
X_MARGIN: float = 0.7                 # free n range beside the first and last bin
XLABEL: str = r"Shell size $n$"


def reweighted_percent(ref: HeAbundanceReference) -> tuple[np.ndarray, np.ndarray]:
    """Abundance of the n >= 1 bins, reweighted to sum to 100 %.

    Parameters
    ----------
    ref
        The abundance reference.

    Returns
    -------
    n, percent : np.ndarray, shape (K,)
        He count and abundance in percent of the n >= 1 ions.

    Raises
    ------
    ValueError
        If the n >= 1 bins carry no abundance.
    """
    solvated = ref.n != 0
    total = ref.ion_fraction[solvated].sum()
    if total <= 0.0:
        raise ValueError("the n >= 1 bins carry no abundance; cannot reweight")
    return ref.n[solvated], ref.ion_fraction[solvated] / total * 100.0


def _bar_figure(
    n: np.ndarray, percent: np.ndarray, ylabel: str, *, log_y: bool = False
) -> Figure:
    """Draw the abundance bars for the bins ``n``; the caller saves and closes it.

    Parameters
    ----------
    n, percent : np.ndarray, shape (K,)
        He count (ascending, contiguous) and abundance in percent.
    ylabel
        Label of the abundance axis.
    log_y
        Logarithmic abundance axis. Its lower limit is the decade below the
        smallest bin, so every bar is visible.

    Returns
    -------
    Figure
        The bar chart.

    Raises
    ------
    ValueError
        If ``log_y`` is set and a bin is not positive.
    """
    fig, ax = plt.subplots(figsize=FIGSIZE_IN, constrained_layout=True)
    ax.bar(n, percent, width=BAR_WIDTH, color=EXP_COLOR, edgecolor="black",
           linewidth=BAR_EDGE_WIDTH)
    ax.set(xlim=(n[0] - X_MARGIN, n[-1] + X_MARGIN), xlabel=XLABEL, ylabel=ylabel)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    if log_y:
        if np.any(percent <= 0.0):
            raise ValueError("a logarithmic abundance axis needs positive bins")
        ax.set_yscale("log")
        ax.set_ylim(bottom=10.0 ** np.floor(np.log10(percent.min())))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    style_axis(ax)
    return fig


def build_figure(ref: HeAbundanceReference) -> Figure:
    """Build the linear n >= 1 abundance figure; the caller saves and closes it."""
    return _bar_figure(*reweighted_percent(ref), r"Fraction ($n \geq 1$) / %")


def build_figure_with_bare(ref: HeAbundanceReference) -> Figure:
    """Build the logarithmic figure of all bins, n = 0 included, in percent of all ions.

    The caller saves and closes it.
    """
    return _bar_figure(ref.n, ref.ion_fraction * 100.0, "Fraction / %", log_y=True)


def main() -> None:
    """Load the reference and write both figures to ``OUT_DIR``."""
    ref = load_he_abundance_reference(ABUNDANCE_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(FIG_RC):
        for build, path in ((build_figure, OUT_PATH), (build_figure_with_bare, OUT_PATH_WITH_BARE)):
            fig = build(ref)
            fig.savefig(path, dpi=DPI)
            plt.close(fig)
            print(f"wrote {path}")


if __name__ == "__main__":
    main()
