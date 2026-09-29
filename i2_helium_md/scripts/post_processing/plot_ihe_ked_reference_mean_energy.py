"""Experiment-only mean-KE figure: the ihe_ked per-fragment reference alone.

Thesis figure for the *experimental reference* chapter: the
``Detected mean KE vs ihe_ked reference`` panel of
``plot_detection_summary.py`` with the **simulation overlay removed** --
no detected-ensemble points, no chi^2 annotation, hence no run directory
is needed. Only ``data/reference/ihe_ked/IHe_KED_reference.csv`` is read.

Recipe parity (rule 1 -- this script re-uses the plot-tuning constants of
``summary_sections.py`` rather than restating them): identical split
layout (left ``n = N_MIN..split-1``, right ``n = split..n_max``, each on
its own linear scale), identical error model -- per-point
``sqrt(stat^2 + sys^2)`` error bars, the two *correlated* fractional
scale bands (calibration, condition) drawn as multiplicative envelopes,
never folded into the per-point bars -- and the same n = 1 median anchor
(I77) marker.

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
from matplotlib.ticker import MaxNLocator  # noqa: E402

from i2_helium_md.postprocess import load_ihe_ked_reference  # noqa: E402

from summary_sections import (  # noqa: E402
    DETECTED_KE_N_MIN,
    IHE_KED_MEAN_ENERGY_SPLIT_N,
)

# =============================================================================
# USER SETTINGS -- edit these and run the script (e.g. from PyCharm)
# =============================================================================
# Frozen per-fragment reference table (one row per n, n = 0..17).
IHE_KED_REFERENCE_CSV: Path = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked" / "IHe_KED_reference.csv"
)

# Output directory and basename; PDF (vector, for LaTeX) + PNG are written.
OUT_DIR: Path = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
OUT_BASENAME: str = "ihe_ked_reference_mean_energy"

# Lowest fragment shown. 1 reproduces the detection-summary panel, whose
# n = 0 bare bin is excluded by a *scoring* convention (DETECTED_KE_N_MIN).
# Set to 0 to also show the bare I+ reference point -- it is measured data
# and legitimate in an experiment-only figure.
N_MIN: int = DETECTED_KE_N_MIN

# Draw the n = 1 median-KE marker (I77 operative anchor) alongside <E>.
SHOW_MEDIAN_ANCHOR: bool = True

# Figure title; set to None for an untitled figure (thesis captions often
# carry the title instead).
TITLE: str | None = r"Experimental I$^+$He$_n$ mean kinetic energy (ihe_ked reference)"

FIGSIZE: tuple[float, float] = (12.5, 5.0)
DPI: int = 300
# =============================================================================


def build_reference_mean_energy_figure(
    ked_ref,
    *,
    n_min: int = DETECTED_KE_N_MIN,
    split: int = IHE_KED_MEAN_ENERGY_SPLIT_N,
    show_median_anchor: bool = True,
    title: str | None = None,
    figsize: tuple[float, float] = (12.5, 5.0),
) -> plt.Figure:
    """Experiment-only split-panel mean-KE figure of the ihe_ked reference.

    Parameters
    ----------
    ked_ref : IHeKedReference
        Loaded reference table (energies eV, fractional bands in [0, 1]).
    n_min
        Lowest fragment n rendered (0 includes the bare I+ point).
    split
        Panel split: left ``n_min <= n < split``, right ``n >= split``.
    show_median_anchor
        Draw the n = 1 median-KE marker (I77 operative anchor).
    title
        Figure suptitle; ``None`` leaves the figure untitled.
    figsize
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
        Two axes; each panel is linear in y. Empty panels are switched off.

    Raises
    ------
    ValueError
        If no fragment row satisfies ``n >= n_min``.
    """
    mean = ked_ref.mean_KE_eV
    selected = ked_ref.n >= n_min
    if not np.any(selected):
        raise ValueError(
            f"no reference rows with n >= {n_min} "
            f"(available n = {ked_ref.n.min()}..{ked_ref.n.max()})"
        )
    n_max = int(ked_ref.n.max())

    fig, axes = plt.subplots(1, 2, figsize=figsize, constrained_layout=True)
    panel_specs = (
        (axes[0], ked_ref.n < split, f"n = {n_min}–{split - 1}", True),
        (axes[1], ked_ref.n >= split, f"n = {split}–{n_max}", False),
    )
    for ax, subset, panel_title, show_legend in panel_specs:
        subset = subset & selected
        if not np.any(subset):
            ax.set_axis_off()
            continue
        for frac, label, alpha in (
            (ked_ref.calib_syst_frac, "calibration band (correlated)", 0.20),
            (ked_ref.condition_syst_frac, "condition band (correlated)", 0.12),
        ):
            ax.fill_between(
                ked_ref.n[subset],
                mean[subset] * (1.0 - frac[subset]),
                mean[subset] * (1.0 + frac[subset]),
                color="tab:blue", alpha=alpha, linewidth=0, label=label,
            )
        ax.errorbar(
            ked_ref.n[subset], mean[subset],
            yerr=ked_ref.point_err_eV[subset], fmt="o",
            color="tab:blue", markersize=4, capsize=2,
            label=r"experiment $\langle E\rangle$ (stat $\oplus$ sys)")
        anchor_hit = subset & (ked_ref.n == 1)
        if show_median_anchor and np.any(anchor_hit):
            ax.plot(ked_ref.n[anchor_hit], ked_ref.median_KE_eV[anchor_hit],
                    "D", markersize=7, markerfacecolor="none",
                    markeredgecolor="tab:purple", markeredgewidth=1.5,
                    label="n = 1 median anchor (I77, operative)")
        ax.set(title=panel_title, xlabel="n (attached He atoms)",
               ylabel="mean kinetic energy / eV")
        # n is an integer count -- suppress matplotlib's fractional ticks
        # (visible on the narrow left panel, where only a few n are shown).
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if show_legend:
            ax.legend(frameon=False, fontsize=8)
    if title is not None:
        fig.suptitle(title, fontsize=11)
    return fig


def main() -> None:
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    fig = build_reference_mean_energy_figure(
        ked_ref,
        n_min=N_MIN,
        show_median_anchor=SHOW_MEDIAN_ANCHOR,
        title=TITLE,
        figsize=FIGSIZE,
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    written = []
    for suffix in (".pdf", ".png"):
        out_path = OUT_DIR / f"{OUT_BASENAME}{suffix}"
        fig.savefig(out_path, dpi=DPI)
        written.append(out_path)
    plt.close(fig)
    print(f"reference: {ked_ref.source_path}")
    print(f"fragments shown: n = {N_MIN}..{int(ked_ref.n.max())}")
    for out_path in written:
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
