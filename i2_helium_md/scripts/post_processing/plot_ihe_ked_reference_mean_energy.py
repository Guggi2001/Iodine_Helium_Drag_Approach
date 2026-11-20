"""Experiment-only mean-KE figure: the ihe_ked per-fragment reference alone.

Thesis figure for the *experimental reference* chapter: the
``Detected mean KE vs ihe_ked reference`` panel of
``plot_detection_summary.py`` with the **simulation overlay removed** --
no detected-ensemble points, no chi^2 annotation, hence no run directory
is needed. Only ``data/reference/ihe_ked/IHe_KED_reference.csv`` is read.

Layout: one full-width panel over all n on a logarithmic energy axis
(``YSCALE``), as panel (a) of the B2 thesis figure ``b2_2_mean_ke``
(``tier2_b2_thesis_figure.py``). Error model as the detection-summary
panel -- per-point ``sqrt(stat^2 + sys^2)`` error bars, the two
*correlated* fractional scale bands (calibration, condition) drawn as
multiplicative envelopes, never folded into the per-point bars -- and the
same n = 1 median anchor (I77) marker.

Two versions are written: ``ihe_ked_reference_mean_energy`` from ``n = 1``, as the
detection-summary panel, and ``ihe_ked_reference_mean_energy_with_bare`` from
``n = 0``, with the bare I+ point.

Units: energies in eV, ``n`` = number of attached He atoms.

Edit the ``USER SETTINGS`` block below and run the script::

    python scripts/post_processing/plot_ihe_ked_reference_mean_energy.py
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

from i2_helium_md.postprocess import load_ihe_ked_reference  # noqa: E402

from summary_sections import DETECTED_KE_N_MIN  # noqa: E402
from thesis_style import EXP_COLOR, FIG_RC, FONT, TEXTWIDTH_IN, style_axis  # noqa: E402

# =============================================================================
# USER SETTINGS -- edit these and run the script (e.g. from PyCharm)
# =============================================================================
# Frozen per-fragment reference table (one row per n, n = 0..17).
IHE_KED_REFERENCE_CSV: Path = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked" / "IHe_KED_reference.csv"
)

OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_rest"

# PNG file stem -> lowest fragment shown. DETECTED_KE_N_MIN reproduces the
# detection-summary panel, whose n = 0 bare bin is excluded by a *scoring*
# convention; 0 adds the bare I+ point -- measured data, legitimate in an
# experiment-only figure.
VARIANTS: dict[str, int] = {
    "ihe_ked_reference_mean_energy": DETECTED_KE_N_MIN,
    "ihe_ked_reference_mean_energy_with_bare": 0,
}

# Draw the n = 1 median-KE marker (I77 operative anchor) alongside <E>.
SHOW_MEDIAN_ANCHOR: bool = True

# Energy axis scale: "log" (as b2_2_mean_ke panel (a)) or "linear".
YSCALE: str = "log"

# Figure title; None = untitled (the thesis caption carries the title).
TITLE: str | None = None

# Thesis text width (thesis_style), so the fonts print at their set size.
FIGSIZE: tuple[float, float] = (TEXTWIDTH_IN, 3.0)
DPI: int = 300
# =============================================================================


def build_reference_mean_energy_figure(
    ked_ref,
    *,
    n_min: int = DETECTED_KE_N_MIN,
    show_median_anchor: bool = True,
    yscale: str = YSCALE,
    title: str | None = None,
    figsize: tuple[float, float] = FIGSIZE,
) -> plt.Figure:
    """Experiment-only single-panel mean-KE figure of the ihe_ked reference.

    Parameters
    ----------
    ked_ref : IHeKedReference
        Loaded reference table (energies eV, fractional bands in [0, 1]).
    n_min
        Lowest fragment n rendered (0 includes the bare I+ point).
    show_median_anchor
        Draw the n = 1 median-KE marker (I77 operative anchor).
    yscale
        ``"log"`` (plain-number tick labels, as b2_2_mean_ke panel (a)) or
        ``"linear"``.
    title
        Figure suptitle; ``None`` leaves the figure untitled.
    figsize
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
        One axis, all n >= ``n_min``.

    Raises
    ------
    ValueError
        If no fragment row satisfies ``n >= n_min``, or ``yscale`` is not
        ``"log"`` / ``"linear"``.
    """
    if yscale not in ("log", "linear"):
        raise ValueError(f"yscale must be 'log' or 'linear', got {yscale!r}")
    mean = ked_ref.mean_KE_eV
    selected = ked_ref.n >= n_min
    if not np.any(selected):
        raise ValueError(
            f"no reference rows with n >= {n_min} "
            f"(available n = {ked_ref.n.min()}..{ked_ref.n.max()})"
        )

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    for frac, label, alpha in (
        (ked_ref.calib_syst_frac, "Calibration band (correlated)", 0.20),
        (ked_ref.condition_syst_frac, "Condition band (correlated)", 0.12),
    ):
        ax.fill_between(
            ked_ref.n[selected],
            mean[selected] * (1.0 - frac[selected]),
            mean[selected] * (1.0 + frac[selected]),
            color=EXP_COLOR, alpha=alpha, linewidth=0, label=label,
        )
    ax.errorbar(
        ked_ref.n[selected], mean[selected],
        yerr=ked_ref.point_err_eV[selected], fmt="o",
        color=EXP_COLOR, markersize=4, capsize=2,
        label=r"Experiment $\bar{E}$ (stat $\oplus$ sys)")
    anchor_hit = selected & (ked_ref.n == 1)
    if show_median_anchor and np.any(anchor_hit):
        ax.plot(ked_ref.n[anchor_hit], ked_ref.median_KE_eV[anchor_hit],
                "D", markersize=6, markerfacecolor="none",
                markeredgecolor=EXP_COLOR, markeredgewidth=1.3,
                label="Experiment median ($n = 1$)")
    ax.set(yscale=yscale, xlabel=r"Shell size $n$", ylabel="Mean kinetic energy / eV")
    if yscale == "log":
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    # n is an integer count -- suppress matplotlib's fractional ticks.
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    style_axis(ax)
    ax.legend(frameon=False, fontsize=FONT["legend_fontsize"], loc="upper right")
    if title is not None:
        fig.suptitle(title, fontsize=FONT["title_fontsize"])
    return fig


def main() -> None:
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"reference: {ked_ref.source_path}")
    with plt.rc_context(FIG_RC):                 # thesis_style font (Palatino)
        for stem, n_min in VARIANTS.items():
            fig = build_reference_mean_energy_figure(
                ked_ref,
                n_min=n_min,
                show_median_anchor=SHOW_MEDIAN_ANCHOR,
                yscale=YSCALE,
                title=TITLE,
                figsize=FIGSIZE,
            )
            out_path = OUT_DIR / f"{stem}.png"
            fig.savefig(out_path, dpi=DPI)
            plt.close(fig)
            print(f"wrote {out_path} (n = {n_min}..{int(ked_ref.n.max())})")


if __name__ == "__main__":
    main()
