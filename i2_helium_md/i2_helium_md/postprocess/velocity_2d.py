"""2-D Cartesian velocity-density histogram for the lab-frame ion run.

Reproduces the (vx, vy) heat-map panel of
``legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse.m``
(velocity_bins = -22:0.5:22, ``hist2`` of ``vx_total, vy_total``). Uses a
plain rebin of the existing 3-D simulated final velocities -- no new physics,
no Abel inversion.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..physics.constants import U as U_KG
from ..simulation.checkpoint import IonCheckpoint


_AxisName = Literal["x", "y", "z"]


@dataclass(frozen=True)
class Velocity2DHistogram:
    """2-D histogram of two final-velocity components.

    Attributes
    ----------
    counts : np.ndarray, shape (n_bins, n_bins)
        Atom count per (axis_a, axis_b) bin.
    bin_centers_a_Aps, bin_centers_b_Aps : np.ndarray, shape (n_bins,)
        Bin centers along each axis in angstrom/ps.
    bin_edges_a_Aps, bin_edges_b_Aps : np.ndarray, shape (n_bins + 1,)
    axis_a, axis_b : str
        Axis labels, one of ``"x"``, ``"y"``, ``"z"``.
    mass_amu : float
    num_atoms_used : int
    """

    counts: np.ndarray
    bin_centers_a_Aps: np.ndarray
    bin_edges_a_Aps: np.ndarray
    bin_centers_b_Aps: np.ndarray
    bin_edges_b_Aps: np.ndarray
    axis_a: str
    axis_b: str
    mass_amu: float
    num_atoms_used: int


def velocity_density_2d(
    ion: IonCheckpoint,
    *,
    axes: tuple[_AxisName, _AxisName] = ("x", "y"),
    n_bins: int = 200,
    v_max_Aps: float = 22.0,
    mass_amu: float | None = None,
    mass_tolerance_amu: float = 0.5,
    require_outside: bool = True,
) -> Velocity2DHistogram:
    """Histogram two components of the final velocity into a square 2-D grid.

    Mirrors the legacy ``velocity_bins = -22:0.5:22`` choice by default.
    """
    if n_bins < 1:
        raise ValueError(f"n_bins must be >= 1, got {n_bins}")
    if v_max_Aps <= 0.0:
        raise ValueError(f"v_max_Aps must be > 0, got {v_max_Aps}")
    if axes[0] == axes[1]:
        raise ValueError(f"axes must differ, got {axes!r}")

    component = {
        "x": ion.velocities_final_x,
        "y": ion.velocities_final_y,
        "z": ion.velocities_final_z,
    }
    for ax in axes:
        if ax not in component:
            raise ValueError(f"axis {ax!r} must be one of 'x','y','z'")

    masses_amu = np.round(np.asarray(ion.mass_final_kg) / U_KG)
    if mass_amu is None:
        mass_mask = np.ones(masses_amu.shape, dtype=bool)
        mass_used = float("nan")
    else:
        mass_mask = np.abs(masses_amu - mass_amu) <= mass_tolerance_amu
        mass_used = float(mass_amu)

    if require_outside:
        outside = np.concatenate(
            [ion.b_ion_outside, ion.b_ion_outside]
        ).astype(bool)
        sel = mass_mask & outside
    else:
        sel = mass_mask

    a = np.asarray(component[axes[0]])[sel]
    b = np.asarray(component[axes[1]])[sel]

    edges = np.linspace(-v_max_Aps, v_max_Aps, n_bins + 1)
    counts, _, _ = np.histogram2d(a, b, bins=(edges, edges))
    centers = 0.5 * (edges[:-1] + edges[1:])

    return Velocity2DHistogram(
        counts=counts.astype(float),
        bin_centers_a_Aps=centers,
        bin_edges_a_Aps=edges,
        bin_centers_b_Aps=centers,
        bin_edges_b_Aps=edges,
        axis_a=axes[0],
        axis_b=axes[1],
        mass_amu=mass_used,
        num_atoms_used=int(sel.sum()),
    )
