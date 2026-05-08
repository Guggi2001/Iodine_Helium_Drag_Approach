"""Plot the HeDFT trajectory-comparison figures for one MD run.

Reproduces the trajectory-only part of the legacy MATLAB workflow:

    legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
        simulation_image_only_trajectories.m

How to use this file
--------------------
Open this file, edit the values in USER SETTINGS, then run it from your
editor or with:

    python scripts/post_processing/plot_hedft_comparison.py

This script intentionally does not plot the experimental VMI velocity
distribution. Use ``plot_experimental_comparison.py`` for that figure,
because it comes from the realistic experimental-condition run.

Output is interactive: ``plt.show()`` opens two windows. No files are
written. Save manually from the matplotlib UI if needed.
"""

from __future__ import annotations

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]


# =============================================================================
# USER SETTINGS
# =============================================================================
# Inputs
RUN_DIR = PROJECT_ROOT / "data" / "runs" / "18A_hedft_comparison"
HEDFT_PATH = PROJECT_ROOT / "data" / "reference" / "18A_All_Data.csv"

# Subset stride for the individual MD trajectories overlaid on Figure 1.
# With the default 2000-molecule production run this gives ~40 distance
# lines. Lower numbers = denser plot.
MD_DISTANCE_STRIDE = 50

# Figure 2: MATLAB sampled about 15 molecules and plotted both atoms,
# so cap the overlay near 30 individual velocity traces.
MD_VELOCITY_MAX_TRACES = 30


# =============================================================================
# IMPORT SETUP
# =============================================================================
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    HedftTrajectory,
    load_hedft_trajectory,
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
        f"Loaded HeDFT-comparison ion checkpoint: N={ion.num_molecules}, "
        f"time=[{ion.time_ps[0]:.3f}, {ion.time_ps[-1]:.3f}] ps "
        f"({ion.time_ps.size} samples)"
    )
    print(
        f"Loaded HeDFT reference: r0={hedft.droplet_radius_A:.1f} A, "
        f"time=[{hedft.time_ps[0]:.3f}, {hedft.time_ps[-1]:.3f}] ps "
        f"({hedft.time_ps.size} samples)"
    )

    _build_distance_figure(ion, hedft)
    _build_velocity_figure(ion, hedft)

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
# Figure 2 -- velocity trajectories
# ===========================================================================
def _build_velocity_figure(
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
    *,
    max_velocity_traces: int = MD_VELOCITY_MAX_TRACES,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(9.5, 4.0), constrained_layout=True)
    _draw_velocity_tile(
        ax,
        ion,
        hedft,
        max_traces=max_velocity_traces,
    )
    return fig


def _draw_velocity_tile(
    ax: plt.Axes,
    ion: IonCheckpoint,
    hedft: HedftTrajectory,
    *,
    max_traces: int,
) -> None:
    """Per-atom MD speed trajectories + HeDFT |v| + mean MD speed."""
    n = ion.num_molecules
    speed_per_atom = np.sqrt(
        ion.velocities_x**2 + ion.velocities_y**2 + ion.velocities_z**2
    )                             # shape (2N, T_md)
    mean_speed = np.mean(speed_per_atom, axis=0)

    indices = _velocity_trace_indices(n, max_traces=max_traces)
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


def _velocity_trace_indices(
    num_molecules: int,
    *,
    max_traces: int,
) -> list[int]:
    """Select paired atoms from evenly spaced molecules, capped by trace count."""
    if max_traces < 1:
        return []

    num_molecule_samples = min(num_molecules, max(1, (max_traces + 1) // 2))
    molecule_indices = np.linspace(
        0,
        num_molecules - 1,
        num=num_molecule_samples,
        dtype=int,
    )

    indices: list[int] = []
    for mol_idx in molecule_indices:
        indices.append(int(mol_idx))
        if len(indices) < max_traces:
            indices.append(int(mol_idx + num_molecules))
    return indices


if __name__ == "__main__":
    raise SystemExit(main())
