"""Time-resolved radial distribution heatmap for an ion run.

Reproduces the time-evolution radial-distribution panel of
``legacy_matlab_repository/post_process_compare_radial_distributions.m``
(``histogram_data_radius(i, :)`` looped over time indices).

Bins ``|r_atom|`` (lab-frame distance from the origin) into a
``(n_time_slices, n_r_bins)`` grid. Time slices are spaced uniformly
across the stored time axis, *not* per stored simulation step, so the
heatmap remains affordable even on long runs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..simulation.checkpoint import IonCheckpoint, NeutralCheckpoint


@dataclass(frozen=True)
class RadialEvolution:
    """Time-evolved radial distribution (heatmap)."""

    counts: np.ndarray  # shape (n_time_slices, n_r_bins)
    time_centers_ps: np.ndarray
    time_indices: np.ndarray
    r_centers_A: np.ndarray
    r_edges_A: np.ndarray


def radial_distribution_evolution(
    ckpt: IonCheckpoint | NeutralCheckpoint,
    *,
    n_time_slices: int = 50,
    n_r_bins: int = 100,
    r_max_A: float | None = None,
) -> RadialEvolution:
    """Bin per-atom ``|r|`` over a sub-sampled time axis.

    The same primitive is used for both stages -- it only reads
    ``positions_x/y/z`` and ``time_ps`` which exist on both
    checkpoint types.

    Parameters
    ----------
    ckpt
        Ion or neutral checkpoint with full trajectory arrays.
    n_time_slices
        Target number of time slices in the heatmap. The actual number
        is ``min(n_time_slices, num_steps)``.
    n_r_bins
        Number of radial bins from 0 to ``r_max_A``.
    r_max_A
        Upper edge of the radial axis. If ``None``, uses ``1.05 *``
        the largest ``|r|`` seen across the chosen time slices.
    """
    if n_time_slices < 1:
        raise ValueError(f"n_time_slices must be >= 1, got {n_time_slices}")
    if n_r_bins < 1:
        raise ValueError(f"n_r_bins must be >= 1, got {n_r_bins}")

    time_ps = np.asarray(ckpt.time_ps)
    num_steps = time_ps.size
    if num_steps == 0:
        raise ValueError("checkpoint has no time samples")

    n_slices = min(int(n_time_slices), num_steps)
    indices = np.unique(
        np.linspace(0, num_steps - 1, n_slices).round().astype(int)
    )
    n_slices = indices.size

    px = np.asarray(ckpt.positions_x)[:, indices]
    py = np.asarray(ckpt.positions_y)[:, indices]
    pz = np.asarray(ckpt.positions_z)[:, indices]
    radial = np.sqrt(px * px + py * py + pz * pz)  # shape (2N, n_slices)

    if r_max_A is None:
        if not np.any(np.isfinite(radial)):
            r_max_A = 1.0
        else:
            r_max_A = float(np.nanmax(radial)) * 1.05
            if r_max_A <= 0.0:
                r_max_A = 1.0

    edges = np.linspace(0.0, r_max_A, n_r_bins + 1)
    counts = np.zeros((n_slices, n_r_bins), dtype=float)
    for k in range(n_slices):
        c, _ = np.histogram(radial[:, k], bins=edges)
        counts[k, :] = c

    return RadialEvolution(
        counts=counts,
        time_centers_ps=time_ps[indices],
        time_indices=indices,
        r_centers_A=0.5 * (edges[:-1] + edges[1:]),
        r_edges_A=edges,
    )
