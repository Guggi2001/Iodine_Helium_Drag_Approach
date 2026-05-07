"""Plot the HeDFT comparison figures for one MD run.

Reproduces the matplotlib equivalent of the legacy MATLAB scripts

    legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
        simulation_image.m                    (full output, 3 plots)
        simulation_image_only_trajectories.m  (subset, 2 plots)

How to use this file
--------------------
Open this file, edit the values in USER SETTINGS, then run it from your
editor or with:

    python scripts/plot_hedft_comparison.py

Set ``SHOW_VMI_TILE = False`` to mirror the trajectories-only output of
``simulation_image_only_trajectories.m`` (Figure 1 + Figure 2 with only
the velocity-vs-time tile, no detector-velocity histograms).

Output is interactive: ``plt.show()`` opens two windows. No files are
written. Save manually from the matplotlib UI if needed.
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]


# =============================================================================
# USER SETTINGS
# =============================================================================
# Inputs
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
HEDFT_PATH = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
VMI_HE_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_iplus_he.csv"
VMI_GAS_PATH = PROJECT_ROOT / "data" / "reference" / "vmi_iplus_gas.csv"

# If False: skip the bottom (velocity-distribution) tile of Figure 2.
# Equivalent to running the legacy ``simulation_image_only_trajectories.m``.
SHOW_VMI_TILE = True

# Subset stride for the individual MD trajectories overlaid on
# Figure 1 (per molecule) and Figure 2 top tile (per atom). With the
# default 2000-molecule production run this gives ~40 distance lines
# and ~266 velocity lines, both readable. Lower numbers = denser plot.
MD_DISTANCE_STRIDE = 50
MD_VELOCITY_STRIDE = 15

# Mass selection (amu) for the simulation curves on the bottom tile.
MASS_I_HE_AMU = 131.0      # I + 1 He
MASS_I_HE2_AMU = 135.0     # I + 2 He

# Histogram resolution for the simulation traces.
HIST_NUM_BINS = 120
HIST_V_MAX_APS = 28.0


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    HedftTrajectory,
    VmiReference,
    compute_final_velocity_histogram,
    load_hedft_trajectory,
    load_vmi_reference,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


# Color and style constants ----------------------------------------------------
_HEDFT_COLOR = "black"
_MD_DISTANCE_COLOR = (1.0, 0.2, 0.6, 0.1)   # magenta with alpha=0.1
_MD_VELOCITY_COLOR = (0.2, 0.2, 0.6, 0.1)   # dark blue with alpha=0.1


def main() -> int:
    ion = RunDirectory(RUN_DIR).load_ion()
    hedft = load_hedft_trajectory(HEDFT_PATH)

    print(
        f"Loaded ion checkpoint: N={ion.num_molecules}, "
        f"time=[{ion.time_ps[0]:.3f}, {ion.time_ps[-1]:.3f}] ps "
        f"({ion.time_ps.size} samples)"
    )
    print(
        f"Loaded HeDFT reference: r0={hedft.droplet_radius_A:.1f} A, "
        f"time=[{hedft.time_ps[0]:.3f}, {hedft.time_ps[-1]:.3f}] ps "
        f"({hedft.time_ps.size} samples)"
    )

    _build_distance_figure(ion, hedft)
    _build_velocity_figure(ion, hedft, show_vmi_tile=SHOW_VMI_TILE)

    plt.show()
    return 0


# ===========================================================================
# Figure 1 -- distance trajectories
# ===========================================================================
def _build_distance_figure(
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7.0, 4.5))

    # Per-molecule MD I-I distance, shape (N, T_md).
    n = ion.num_molecules
    dx = ion.positions_x[:n] - ion.positions_x[n:]
    dy = ion.positions_y[:n] - ion.positions_y[n:]
    dz = ion.positions_z[:n] - ion.positions_z[n:]
    dR_per_mol = np.sqrt(dx * dx + dy * dy + dz * dz)

    # Subset to keep the figure legible.
    indices = range(0, n, max(1, MD_DISTANCE_STRIDE))
    for plot_idx, mol_idx in enumerate(indices):
        ax.plot(
            ion.time_ps,
            dR_per_mol[mol_idx],
            color=_MD_DISTANCE_COLOR,
            linewidth=1.0,
            label="MD trajectories" if plot_idx == 0 else None,
        )

    ax.plot(
        hedft.time_ps,
        hedft.distance_A,
        color=_HEDFT_COLOR,
        linewidth=1.5,
        label="HeDFT",
    )

    ax.set_xlim(0.0, 6.0)
    ax.set_ylim(8.0, 40.0)
    ax.set_xlabel("t / ps")
    ax.set_ylabel(r"$R_1 - R_2$ / $\mathrm{\AA}$")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


# ===========================================================================
# Figure 2 -- velocity panels
# ===========================================================================
def _build_velocity_figure(
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
    *,
    show_vmi_tile: bool,
) -> plt.Figure:
    if show_vmi_tile:
        fig = plt.figure(figsize=(9.5, 7.5), constrained_layout=True)
        gs = GridSpec(2, 1, figure=fig)
        ax_top = fig.add_subplot(gs[0, 0])
        ax_bot = fig.add_subplot(gs[1, 0])
    else:
        fig, ax_top = plt.subplots(figsize=(9.5, 4.0), constrained_layout=True)
        ax_bot = None

    _draw_velocity_top_tile(ax_top, ion, hedft)
    _add_subplot_label(ax_top, "a")

    if ax_bot is not None:
        vmi_he = load_vmi_reference(VMI_HE_PATH)
        vmi_gas = load_vmi_reference(VMI_GAS_PATH)
        hist_he = compute_final_velocity_histogram(
            ion,
            mass_amu=MASS_I_HE_AMU,
            num_bins=HIST_NUM_BINS,
            v_max_Aps=HIST_V_MAX_APS,
        )
        hist_he2 = compute_final_velocity_histogram(
            ion,
            mass_amu=MASS_I_HE2_AMU,
            num_bins=HIST_NUM_BINS,
            v_max_Aps=HIST_V_MAX_APS,
        )
        _draw_velocity_distribution_tile(
            ax_bot,
            vmi_gas=vmi_gas,
            vmi_he=vmi_he,
            sim_he=hist_he,
            sim_he2=hist_he2,
        )
        _add_subplot_label(ax_bot, "b")

    return fig


def _draw_velocity_top_tile(
    ax: plt.Axes,
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
) -> None:
    """Per-atom MD speed trajectories + HeDFT |v| + mean MD speed."""
    n = ion.num_molecules
    speed_per_atom = np.sqrt(
        ion.velocities_x**2 + ion.velocities_y**2 + ion.velocities_z**2
    )                             # shape (2N, T_md)
    mean_speed = np.mean(speed_per_atom, axis=0)

    indices = range(0, 2 * n, max(1, MD_VELOCITY_STRIDE))
    for plot_idx, atom_idx in enumerate(indices):
        ax.plot(
            ion.time_ps,
            speed_per_atom[atom_idx],
            color=_MD_VELOCITY_COLOR,
            linewidth=1.0,
            label="MD velocity" if plot_idx == 0 else None,
        )

    ax.plot(
        hedft.time_ps,
        hedft.v2_magnitude_Aps,
        linestyle=":",
        linewidth=2.0,
        color=_HEDFT_COLOR,
        label="HeDFT",
    )
    ax.plot(
        ion.time_ps,
        mean_speed,
        linestyle="--",
        color=_HEDFT_COLOR,
        linewidth=1.5,
        label="mean MD velocity",
    )

    ax.set_xlim(0.0, 12.0)
    ax.set_xlabel("t / ps")
    ax.set_ylabel(r"v / $\mathrm{\AA}/\mathrm{ps}$")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _draw_velocity_distribution_tile(
    ax: plt.Axes,
    *,
    vmi_gas: VmiReference,
    vmi_he: VmiReference,
    sim_he,
    sim_he2,
) -> None:
    """Bottom tile: experimental I+ gas, I+He droplet, simulation I+He, I+He2."""
    # 5-element discrete colormap that approximates colorcet('L08','N',5).
    palette = plt.colormaps["viridis"](np.linspace(0.05, 0.85, 5))
    c_gas, c_he, c_sim_he, c_sim_he2, _ = palette

    # Normalisation for the experimental traces, mirroring the mask-then-
    # max idiom in scripts/post_processing_comparison/compare.py:23-30.
    mask_gas = vmi_gas.velocity_Aps > 4.0
    max_gas = float(vmi_gas.signal_arb[mask_gas].max())
    max_he = float(vmi_he.signal_arb.max())

    ax.plot(
        vmi_gas.velocity_Aps,
        vmi_gas.signal_arb / max_gas,
        color=c_gas,
        linewidth=2.0,
        label=r"$I_2$:$I^+$",
    )
    ax.plot(
        vmi_he.velocity_Aps,
        vmi_he.signal_arb / max_he,
        linestyle=":",
        color=c_he,
        linewidth=2.0,
        label=r"$I_2 He_N$:$I^+ He$",
    )
    ax.plot(
        sim_he.bin_centers_Aps,
        sim_he.density / sim_he.density.max(),
        linestyle="--",
        color=c_sim_he,
        linewidth=2.0,
        label=r"simulation $I^+ He$",
    )
    ax.plot(
        sim_he2.bin_centers_Aps,
        sim_he2.density / max(sim_he2.density.max(), 1e-30),
        linestyle="-.",
        color=c_sim_he2,
        linewidth=2.0,
        label=r"simulation $I^+ He_2$",
    )

    ax.set_xlim(0.0, HIST_V_MAX_APS)
    ax.set_ylim(0.0, 1.1)
    ax.set_xlabel(r"v / $\mathrm{\AA}/\mathrm{ps}$")
    ax.set_ylabel("signal / arb. units")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _add_subplot_label(ax: plt.Axes, letter: str) -> None:
    """Top-left subplot letter, matching MATLAB add_letter_norm size."""
    ax.text(
        0.02,
        0.95,
        letter,
        transform=ax.transAxes,
        fontsize=20,
        fontweight="bold",
        va="top",
        ha="left",
    )


if __name__ == "__main__":
    raise SystemExit(main())
