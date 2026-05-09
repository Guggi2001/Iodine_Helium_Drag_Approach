"""Tests for i2_helium_md/postprocess/pair_correlation.py."""

from __future__ import annotations

import numpy as np

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.pair_correlation import (
    CovarianceMatrix,
    DistanceHistogram,
    angular_pair_covariance,
    interparticle_distance_histogram,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint


def _make_ion(
    *,
    pos_final_x: np.ndarray,
    pos_final_y: np.ndarray,
    pos_final_z: np.ndarray,
    vx_final: np.ndarray,
    vy_final: np.ndarray,
    vz_final: np.ndarray,
    masses_amu: np.ndarray,
    num_steps: int = 4,
    b_outside: np.ndarray | None = None,
) -> IonCheckpoint:
    n2 = pos_final_x.size
    n = n2 // 2
    if b_outside is None:
        b_outside = np.ones(n, dtype=bool)
    masses_kg = masses_amu.astype(float) * U_KG
    diag_zero = np.zeros((2 * n, num_steps))
    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, num_steps),
        positions_x=np.zeros((2 * n, num_steps)),
        positions_y=np.zeros((2 * n, num_steps)),
        positions_z=np.zeros((2 * n, num_steps)),
        velocities_x=np.zeros((2 * n, num_steps)),
        velocities_y=np.zeros((2 * n, num_steps)),
        velocities_z=np.zeros((2 * n, num_steps)),
        positions_final_x=pos_final_x.astype(float),
        positions_final_y=pos_final_y.astype(float),
        positions_final_z=pos_final_z.astype(float),
        velocities_final_x=vx_final.astype(float),
        velocities_final_y=vy_final.astype(float),
        velocities_final_z=vz_final.astype(float),
        mass_kg=masses_kg.copy(),
        mass_final_kg=masses_kg,
        mass_history_kg=np.broadcast_to(
            masses_kg[:, None], (2 * n, num_steps)
        ).copy(),
        droplet_radii_angstrom=np.full(2 * n, 30.0),
        E_kin_eV=diag_zero,
        E_pot_eV=diag_zero,
        E_dissip_eV=diag_zero,
        E_mass_attach_defect_eV=diag_zero,
        b_ion_outside=b_outside,
        relative_loss_per_ps=diag_zero,
        number_of_collisions=np.zeros((2 * n, num_steps), dtype=int),
        temperature_diagnostic=np.full((num_steps, 3), np.nan),
    )


class TestInterparticleDistanceHistogram:
    def test_known_distance_lands_in_correct_bin(self):
        # Two molecules. Molecule 0: |r1 - r2| = 5.0. Molecule 1: |r1 - r2| = 3.0.
        ion = _make_ion(
            pos_final_x=np.array([5.0, 1.5, 0.0, -1.5]),
            pos_final_y=np.zeros(4),
            pos_final_z=np.zeros(4),
            vx_final=np.zeros(4), vy_final=np.zeros(4), vz_final=np.zeros(4),
            masses_amu=np.full(4, 131.0),
        )
        hist = interparticle_distance_histogram(
            ion, num_bins=10, max_distance_A=10.0,
        )
        assert isinstance(hist, DistanceHistogram)
        # Bin width = 1.0 -> distances 5.0 and 3.0 land in bins index 5 and 3.
        assert hist.counts[5] == 1
        assert hist.counts[3] == 1
        assert hist.counts.sum() == 2
        assert hist.num_pairs_used == 2

    def test_default_max_distance_uses_data_extent(self):
        ion = _make_ion(
            pos_final_x=np.array([0.0, 0.0, -2.0, 2.0]),
            pos_final_y=np.zeros(4),
            pos_final_z=np.zeros(4),
            vx_final=np.zeros(4), vy_final=np.zeros(4), vz_final=np.zeros(4),
            masses_amu=np.full(4, 131.0),
        )
        hist = interparticle_distance_histogram(ion, num_bins=10)
        assert hist.counts.sum() == 2
        assert hist.bin_edges_A[-1] >= 2.0


class TestAngularPairCovariance:
    def test_diagonal_removed(self):
        # Construct molecules whose pair angles align (theta_a = theta_b),
        # so all counts would land on the diagonal pre-removal.
        n = 8
        # I1 atoms at velocities (1, 0); I2 atoms at velocities (1, 0).
        vx = np.concatenate([np.ones(n), np.ones(n)])
        vy = np.concatenate([np.zeros(n), np.zeros(n)])
        vz = np.zeros(2 * n)
        ion = _make_ion(
            pos_final_x=np.zeros(2 * n),
            pos_final_y=np.zeros(2 * n),
            pos_final_z=np.zeros(2 * n),
            vx_final=vx, vy_final=vy, vz_final=vz,
            masses_amu=np.full(2 * n, 131.0),
        )
        cov = angular_pair_covariance(
            ion, n_theta_bins=12, mass_amu=131.0,
        )
        assert isinstance(cov, CovarianceMatrix)
        assert cov.counts.shape == (12, 12)
        np.testing.assert_array_equal(np.diag(cov.counts), 0.0)

    def test_offdiagonal_filled_when_theta_differs(self):
        # Two molecules. Molecule 0: I1 at +x, I2 at +y -> theta_a != theta_b.
        # Molecule 1: identical to Molecule 0.
        vx = np.array([1.0, 1.0, 0.0, 0.0])
        vy = np.array([0.0, 0.0, 1.0, 1.0])
        vz = np.zeros(4)
        ion = _make_ion(
            pos_final_x=np.zeros(4),
            pos_final_y=np.zeros(4),
            pos_final_z=np.zeros(4),
            vx_final=vx, vy_final=vy, vz_final=vz,
            masses_amu=np.full(4, 131.0),
        )
        cov = angular_pair_covariance(
            ion, n_theta_bins=8, mass_amu=131.0, remove_diagonal=False,
        )
        assert cov.counts.sum() == 2

    def test_outside_filter(self):
        n = 2
        vx = np.array([1.0, 1.0, 1.0, 1.0])
        vy = np.array([0.0, 0.0, 0.5, 0.5])
        vz = np.zeros(4)
        ion = _make_ion(
            pos_final_x=np.zeros(4),
            pos_final_y=np.zeros(4),
            pos_final_z=np.zeros(4),
            vx_final=vx, vy_final=vy, vz_final=vz,
            masses_amu=np.full(4, 131.0),
            b_outside=np.array([True, False]),
        )
        cov = angular_pair_covariance(
            ion, n_theta_bins=4, mass_amu=131.0,
            remove_diagonal=False,
        )
        assert cov.num_pairs_used == 1
