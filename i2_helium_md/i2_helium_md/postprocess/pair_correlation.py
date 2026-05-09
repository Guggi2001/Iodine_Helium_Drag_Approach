"""Inter-particle distance histogram and angular pair covariance.

Final-state pair-correlation diagnostics for an ion run.

* :func:`interparticle_distance_histogram` mirrors the
  ``histogram_data_interatomic_distance`` block of
  ``legacy_matlab_repository/post_process_compare_radial_distributions.m``.
* :func:`angular_pair_covariance` mirrors the
  ``simulated_angular_covariance`` block of
  ``legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v4.m``,
  which bins ``atan2(vx, vy)`` for each ion-pair (atom_a, atom_b) into a
  2-D histogram with the diagonal (i.e. matched-bin self-correlation)
  removed.

Atom layout convention (same as the rest of the package): indices
``[0, N)`` are the first atom of each molecule, indices ``[N, 2 N)``
the second atom, so the per-molecule pair is ``(atom i, atom i + N)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..physics.constants import U as U_KG
from ..simulation.checkpoint import IonCheckpoint


@dataclass(frozen=True)
class DistanceHistogram:
    """1-D histogram of final per-molecule I-I separation."""

    bin_centers_A: np.ndarray
    bin_edges_A: np.ndarray
    counts: np.ndarray
    num_pairs_used: int


@dataclass(frozen=True)
class CovarianceMatrix:
    """2-D angular pair covariance ``(theta_a, theta_b)``.

    The diagonal of ``counts`` is set to 0 so that the overall
    histogram structure is dominated by the inter-atom angular
    correlation, not by trivial self-binning -- this matches the
    ``cov_angular - diag(...)`` step in the legacy MATLAB script.
    """

    counts: np.ndarray
    theta_centers_rad: np.ndarray
    theta_edges_rad: np.ndarray
    num_pairs_used: int


def interparticle_distance_histogram(
    ion: IonCheckpoint,
    *,
    num_bins: int = 100,
    max_distance_A: float | None = None,
) -> DistanceHistogram:
    """Bin the final per-molecule I-I separation ``|r_a - r_b|``.

    Uses ``positions_final_*`` (the asymptotic position recorded by the
    ion stage) so the result is independent of the chosen end-time of
    the trajectory.
    """
    if num_bins < 1:
        raise ValueError(f"num_bins must be >= 1, got {num_bins}")

    n = ion.num_molecules
    dx = np.asarray(ion.positions_final_x[:n] - ion.positions_final_x[n:])
    dy = np.asarray(ion.positions_final_y[:n] - ion.positions_final_y[n:])
    dz = np.asarray(ion.positions_final_z[:n] - ion.positions_final_z[n:])
    distance = np.sqrt(dx * dx + dy * dy + dz * dz)

    if max_distance_A is None:
        if distance.size == 0 or not np.any(np.isfinite(distance)):
            max_distance_A = 1.0
        else:
            max_distance_A = float(np.nanmax(distance)) * 1.05
            if max_distance_A <= 0.0:
                max_distance_A = 1.0

    edges = np.linspace(0.0, max_distance_A, num_bins + 1)
    counts, _ = np.histogram(distance, bins=edges)
    return DistanceHistogram(
        bin_centers_A=0.5 * (edges[:-1] + edges[1:]),
        bin_edges_A=edges,
        counts=counts.astype(int),
        num_pairs_used=int(distance.size),
    )


def angular_pair_covariance(
    ion: IonCheckpoint,
    *,
    n_theta_bins: int = 50,
    mass_amu: float | None = None,
    mass_tolerance_amu: float = 0.5,
    require_outside: bool = True,
    remove_diagonal: bool = True,
) -> CovarianceMatrix:
    """2-D angular histogram ``(theta_a, theta_b)`` of ion-pair velocities.

    For each molecule both atoms must pass the mass + outside selection;
    a single failing atom drops the molecule out of the pair sample.

    ``theta = arctan2(vx, vy) + pi`` is wrapped into ``[0, 2 pi)`` so the
    result has the same bin layout as the 1-D phi histogram in
    :mod:`postprocess.energy_balance`.
    """
    if n_theta_bins < 1:
        raise ValueError(f"n_theta_bins must be >= 1, got {n_theta_bins}")

    n = ion.num_molecules
    masses_amu = np.round(np.asarray(ion.mass_final_kg) / U_KG)
    if mass_amu is None:
        mass_mask_atom = np.ones(masses_amu.shape, dtype=bool)
    else:
        mass_mask_atom = np.abs(masses_amu - mass_amu) <= mass_tolerance_amu
    if require_outside:
        outside_atom = np.concatenate(
            [ion.b_ion_outside, ion.b_ion_outside]
        ).astype(bool)
        atom_ok = mass_mask_atom & outside_atom
    else:
        atom_ok = mass_mask_atom

    pair_ok = atom_ok[:n] & atom_ok[n:]

    vx_a = np.asarray(ion.velocities_final_x[:n])[pair_ok]
    vy_a = np.asarray(ion.velocities_final_y[:n])[pair_ok]
    vx_b = np.asarray(ion.velocities_final_x[n:])[pair_ok]
    vy_b = np.asarray(ion.velocities_final_y[n:])[pair_ok]

    theta_a = np.mod(np.arctan2(vx_a, vy_a) + np.pi, 2.0 * np.pi)
    theta_b = np.mod(np.arctan2(vx_b, vy_b) + np.pi, 2.0 * np.pi)

    edges = np.linspace(0.0, 2.0 * np.pi, n_theta_bins + 1)
    counts, _, _ = np.histogram2d(theta_a, theta_b, bins=(edges, edges))
    if remove_diagonal:
        np.fill_diagonal(counts, 0.0)
    return CovarianceMatrix(
        counts=counts.astype(float),
        theta_centers_rad=0.5 * (edges[:-1] + edges[1:]),
        theta_edges_rad=edges,
        num_pairs_used=int(pair_ok.sum()),
    )
