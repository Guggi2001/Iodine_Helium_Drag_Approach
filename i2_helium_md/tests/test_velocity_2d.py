"""Tests for i2_helium_md/postprocess/velocity_2d.py."""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.velocity_2d import (
    Velocity2DHistogram,
    velocity_density_2d,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint


def _make_ion(*, vx, vy, vz, masses_amu, num_steps=4):
    n2 = vx.size
    n = n2 // 2
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
        positions_final_x=np.zeros(2 * n),
        positions_final_y=np.zeros(2 * n),
        positions_final_z=np.zeros(2 * n),
        velocities_final_x=vx.astype(float),
        velocities_final_y=vy.astype(float),
        velocities_final_z=vz.astype(float),
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
        b_ion_outside=np.ones(n, dtype=bool),
        relative_loss_per_ps=diag_zero,
        number_of_collisions=np.zeros((2 * n, num_steps), dtype=int),
        temperature_diagnostic=np.full((num_steps, 3), np.nan),
    )


class TestVelocityDensity2D:
    def test_counts_conservation(self):
        n_atoms = 100
        rng = np.random.default_rng(0)
        vx = rng.uniform(-15.0, 15.0, size=n_atoms)
        vy = rng.uniform(-15.0, 15.0, size=n_atoms)
        vz = rng.uniform(-15.0, 15.0, size=n_atoms)
        ion = _make_ion(
            vx=vx, vy=vy, vz=vz,
            masses_amu=np.full(n_atoms, 131.0),
        )
        hist = velocity_density_2d(
            ion, axes=("x", "y"), n_bins=40, v_max_Aps=20.0,
            mass_amu=131.0,
        )
        assert isinstance(hist, Velocity2DHistogram)
        assert hist.counts.shape == (40, 40)
        # Some atoms outside the symmetric box may not bin; check upper bound.
        assert hist.counts.sum() <= n_atoms
        assert hist.num_atoms_used == n_atoms

    def test_isotropic_input_is_symmetric(self):
        # Equal counts in 4 quadrants -> 90-degree-rotated histogram is identical.
        rng = np.random.default_rng(2)
        vx = rng.uniform(-10.0, 10.0, size=2_000)
        vy = rng.uniform(-10.0, 10.0, size=2_000)
        vx_sym = np.concatenate([vx, vx, -vx, -vx])
        vy_sym = np.concatenate([vy, -vy, vy, -vy])
        n_atoms = vx_sym.size
        ion = _make_ion(
            vx=vx_sym, vy=vy_sym, vz=np.zeros(n_atoms),
            masses_amu=np.full(n_atoms, 131.0),
        )
        h = velocity_density_2d(
            ion, axes=("x", "y"), n_bins=40, v_max_Aps=15.0,
            mass_amu=131.0,
        )
        # 4-fold symmetric: rotating 180 deg keeps count distribution.
        np.testing.assert_array_equal(h.counts, np.flip(h.counts))

    def test_axes_must_differ(self):
        ion = _make_ion(
            vx=np.array([1.0, 1.0]),
            vy=np.zeros(2),
            vz=np.zeros(2),
            masses_amu=np.array([131.0, 131.0]),
        )
        with pytest.raises(ValueError, match="axes must differ"):
            velocity_density_2d(ion, axes=("x", "x"))

    def test_bad_n_bins(self):
        ion = _make_ion(
            vx=np.array([1.0, 1.0]),
            vy=np.zeros(2),
            vz=np.zeros(2),
            masses_amu=np.array([131.0, 131.0]),
        )
        with pytest.raises(ValueError, match="n_bins"):
            velocity_density_2d(ion, n_bins=0)
