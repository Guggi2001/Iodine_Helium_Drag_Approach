"""TDDFT reference figure: iodine speeds |v|(t) of the 9 A and 18 A trajectories.

(a) 9 A: both atoms and the exit-speed lower bound (mean speed of the bottom atom
over ``EXIT_WINDOW_PS``). (b) 18 A: both atoms, the smoothed reference speed of
the bottom atom, the approximate bubble exit and the violent onset threshold.
Atom 1 of the CSVs is the bottom iodine, atom 2 the top one.

Drawn at the thesis text width with the shared ``thesis_style``; lines as in the
microscopic-calibration figures (``tier0_drag_comparison``).

Reads ``data/reference/TDDFT/{9A,18A}_All_Data.csv`` and
``data/reference/drag/18A/velocity_smoothed/cleaned_data_long.csv``.

Usage:  python data/reference/scripts/visualize_trajectory_velocity.py
Writes: data/figures/thesis_rest/TDDFT_Reference.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess import (  # noqa: E402
    load_hedft_trajectory,
    load_smoothed_speed_reference,
)
from i2_helium_md.postprocess.hedft_loader import (  # noqa: E402
    HedftTrajectory,
    SmoothedSpeedReference,
)
from scripts.post_processing.thesis_style import (  # noqa: E402
    FIG_RC,
    FONT,
    TEXTWIDTH_IN,
    panel_letter,
    style_axis,
)

REFERENCE_DIR: Path = PROJECT_ROOT / "data" / "reference"
HEDFT_9A_PATH: Path = REFERENCE_DIR / "TDDFT" / "9A_All_Data.csv"
HEDFT_18A_PATH: Path = REFERENCE_DIR / "TDDFT" / "18A_All_Data.csv"
SMOOTHED_18A_PATH: Path = (
    REFERENCE_DIR / "drag" / "18A" / "velocity_smoothed" / "cleaned_data_long.csv"
)
OUT_PATH: Path = PROJECT_ROOT / "data" / "figures" / "thesis_rest" / "TDDFT_Reference.png"

EXIT_WINDOW_PS: tuple[float, float] = (13.0, 14.0)  # 9 A exit-speed window
BUBBLE_EXIT_18A_PS: float = 8.5                     # approximate bubble exit, 18 A
VIOLENT_ONSET_18A_PS: float = 4.54                  # start of the 18 A fit window

FIGSIZE_IN: tuple[float, float] = (TEXTWIDTH_IN, 2.9)
DPI: int = 200
T_XLIM_PS: tuple[float, float] = (0.0, 15.0)
T_TICK_STEP_PS: float = 5.0
V_YLIM_9A_APS: tuple[float, float] = (0.0, 6.5)     # headroom for the panel letter
V_YLIM_18A_APS: tuple[float, float] = (0.0, 4.0)
BOTTOM_LINE: dict = dict(color="black", ls="-", lw=1.0, alpha=0.8)
TOP_LINE: dict = dict(color="dimgray", ls=":", lw=1.0, alpha=0.8)
SMOOTHED_LINE: dict = dict(color="tab:purple", lw=1.2)
EXIT_BOUND_LINE: dict = dict(color="green", ls="--", lw=1.2)
BUBBLE_EXIT_LINE: dict = dict(color="black", ls="--", lw=0.8)
VIOLENT_ONSET_LINE: dict = dict(color="black", ls=":", lw=1.0)
LEGEND_NCOLS: int = 3
XLABEL: str = r"Time $t$ / ps"
YLABEL: str = r"Speed $|v|$ / Å ps$^{-1}$"


def exit_speed_bound_Aps(hedft: HedftTrajectory) -> float:
    """Mean bottom-atom speed over ``EXIT_WINDOW_PS`` (inclusive).

    Parameters
    ----------
    hedft
        The 9 A trajectory.

    Returns
    -------
    float
        Exit-speed lower bound in angstrom/ps.

    Raises
    ------
    ValueError
        If no sample lies inside the window.
    """
    t_start_ps, t_end_ps = EXIT_WINDOW_PS
    in_window = (hedft.time_ps >= t_start_ps) & (hedft.time_ps <= t_end_ps)
    if not in_window.any():
        raise ValueError(f"no TDDFT sample in the exit window {EXIT_WINDOW_PS} ps")
    return float(hedft.v1_magnitude_Aps[in_window].mean())


def _plot_atom_speeds(ax: plt.Axes, hedft: HedftTrajectory) -> None:
    """Draw |v|(t) of the bottom and top iodine on ``ax``."""
    ax.plot(hedft.time_ps, hedft.v1_magnitude_Aps, **BOTTOM_LINE, label=r"HeDFT $|v_1|$")
    ax.plot(hedft.time_ps, hedft.v2_magnitude_Aps, **TOP_LINE, label=r"HeDFT $|v_2|$")


def _finish_axis(ax: plt.Axes, title: str, ylim: tuple[float, float], letter: str) -> None:
    """Set the title, labels and limits of ``ax`` and apply the shared style."""
    ax.set_title(title, fontsize=FONT["title_fontsize"])
    ax.set(xlim=T_XLIM_PS, ylim=ylim, xlabel=XLABEL, ylabel=YLABEL)
    ax.xaxis.set_major_locator(MultipleLocator(T_TICK_STEP_PS))
    style_axis(ax)
    panel_letter(ax, letter)


def build_figure(
    hedft_9A: HedftTrajectory,
    hedft_18A: HedftTrajectory,
    smoothed_18A: SmoothedSpeedReference,
) -> Figure:
    """Build the two-panel TDDFT reference figure.

    Parameters
    ----------
    hedft_9A, hedft_18A
        The 9 A and 18 A trajectories.
    smoothed_18A
        Smoothed 18 A reference speed of the bottom atom.

    Returns
    -------
    Figure
        The figure; the caller saves and closes it.
    """
    fig, (ax_9A, ax_18A) = plt.subplots(1, 2, figsize=FIGSIZE_IN, constrained_layout=True)

    _plot_atom_speeds(ax_9A, hedft_9A)
    ax_9A.axhline(exit_speed_bound_Aps(hedft_9A), **EXIT_BOUND_LINE, label="Lower bound exit velocity")
    _finish_axis(ax_9A, "9 Å", V_YLIM_9A_APS, "a")

    _plot_atom_speeds(ax_18A, hedft_18A)
    ax_18A.plot(smoothed_18A.time_ps, smoothed_18A.speed_Aps, **SMOOTHED_LINE,
                label=r"Smoothed $|v_1|$")
    ax_18A.axvline(BUBBLE_EXIT_18A_PS, **BUBBLE_EXIT_LINE, label="Approximate bubble exit")
    ax_18A.axvline(VIOLENT_ONSET_18A_PS, **VIOLENT_ONSET_LINE, label="Violent onset threshold")
    _finish_axis(ax_18A, "18 Å", V_YLIM_18A_APS, "b")

    entries: dict[str, object] = {}
    for ax in (ax_9A, ax_18A):
        for handle, label in zip(*ax.get_legend_handles_labels()):
            entries.setdefault(label, handle)
    fig.legend(entries.values(), entries.keys(), loc="outside upper center", ncols=LEGEND_NCOLS,
               frameon=False, fontsize=FONT["legend_fontsize"])
    return fig


def main() -> None:
    """Load the references and write ``OUT_PATH``."""
    hedft_9A = load_hedft_trajectory(HEDFT_9A_PATH)
    hedft_18A = load_hedft_trajectory(HEDFT_18A_PATH)
    smoothed_18A = load_smoothed_speed_reference(SMOOTHED_18A_PATH)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(FIG_RC):
        fig = build_figure(hedft_9A, hedft_18A, smoothed_18A)
        fig.savefig(OUT_PATH, dpi=DPI)
        plt.close(fig)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
