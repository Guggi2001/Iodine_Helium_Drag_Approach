"""Tests for i2_helium_md/postprocess/time_resolved.py."""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.time_resolved import (
    RadialEvolution,
    radial_distribution_evolution,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint


def _make_ion(*, positions: np.ndarray) -> IonCheckpoint:
    """positions has shape (3, 2N, T) for x, y, z."""
    n2, t_steps = positions.shape[1], positions.shape[2]
    n = n2 // 2
    diag_zero = np.zeros((2 * n, t_steps))
    masses_kg = np.full(2 * n, 127.0) * U_KG
    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 5.0, t_steps),
        positions_x=positions[0],
        positions_y=positions[1],
        positions_z=positions[2],
        velocities_x=np.zeros((2 * n, t_steps)),
        velocities_y=np.zeros((2 * n, t_steps)),
        velocities_z=np.zeros((2 * n, t_steps)),
        positions_final_x=positions[0, :, -1].copy(),
        positions_final_y=positions[1, :, -1].copy(),
        positions_final_z=positions[2, :, -1].copy(),
        velocities_final_x=np.zeros(2 * n),
        velocities_final_y=np.zeros(2 * n),
        velocities_final_z=np.zeros(2 * n),
        mass_kg=masses_kg.copy(),
        mass_final_kg=masses_kg,
        mass_history_kg=np.broadcast_to(
            masses_kg[:, None], (2 * n, t_steps)
        ).copy(),
        droplet_radii_angstrom=np.full(2 * n, 30.0),
        E_kin_eV=diag_zero,
        E_pot_eV=diag_zero,
        E_dissip_eV=diag_zero,
        E_mass_attach_defect_eV=diag_zero,
        b_ion_outside=np.ones(n, dtype=bool),
        relative_loss_per_ps=diag_zero,
        number_of_collisions=np.zeros((2 * n, t_steps), dtype=int),
        temperature_diagnostic=np.full((t_steps, 3), np.nan),
    )


class TestRadialDistributionEvolution:
    def test_shape_and_counts(self):
        n2, t_steps = 8, 20
        positions = np.zeros((3, n2, t_steps))
        # Linearly expanding: |r| = t for every atom on the x-axis.
        positions[0] = np.broadcast_to(
            np.linspace(0.0, 5.0, t_steps), (n2, t_steps)
        )
        ion = _make_ion(positions=positions)

        ev = radial_distribution_evolution(
            ion, n_time_slices=5, n_r_bins=10, r_max_A=6.0,
        )
        assert isinstance(ev, RadialEvolution)
        assert ev.counts.shape[1] == 10
        assert ev.counts.shape[0] <= 5
        # Each slice has all 8 atoms binned (assuming r_max>=5.0 covers them).
        np.testing.assert_array_equal(ev.counts.sum(axis=1), n2)

    def test_each_atom_appears_once_per_slice(self):
        # All atoms at fixed |r| = 2.0 -> all counts in the bin covering 2.0.
        n2, t_steps = 4, 10
        positions = np.zeros((3, n2, t_steps))
        positions[0] = 2.0  # constant
        ion = _make_ion(positions=positions)

        ev = radial_distribution_evolution(
            ion, n_time_slices=4, n_r_bins=5, r_max_A=5.0,
        )
        # Bin width 1.0 -> 2.0 falls in bin index 2.
        np.testing.assert_array_equal(ev.counts[:, 2], n2)
        np.testing.assert_array_equal(ev.counts[:, [0, 1, 3, 4]], 0)

    def test_bad_args(self):
        n2, t_steps = 2, 4
        positions = np.zeros((3, n2, t_steps))
        ion = _make_ion(positions=positions)
        with pytest.raises(ValueError, match="n_time_slices"):
            radial_distribution_evolution(ion, n_time_slices=0)
        with pytest.raises(ValueError, match="n_r_bins"):
            radial_distribution_evolution(ion, n_r_bins=0)
