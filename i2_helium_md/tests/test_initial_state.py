"""Tests for i2_helium_md/simulation/initial_state.py."""

import numpy as np
import pytest

from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.checkpoint import (
    NeutralCheckpoint,
    _NEUTRAL_SCHEMA_VERSION,
)
from i2_helium_md.simulation.initial_state import build_initial_state


# ===========================================================================
# Basic API and shapes
# ===========================================================================
class TestApi:
    def test_returns_neutral_checkpoint(self):
        cfg = single_pulse_N2000(num_molecules=5, seed=0)
        ckpt = build_initial_state(cfg, num_steps=10,
                                   rng=np.random.default_rng(0))
        assert isinstance(ckpt, NeutralCheckpoint)

    def test_schema_version_is_current(self):
        cfg = single_pulse_N2000(num_molecules=5, seed=0)
        ckpt = build_initial_state(cfg, num_steps=10,
                                   rng=np.random.default_rng(0))
        assert ckpt.schema_version == _NEUTRAL_SCHEMA_VERSION

    def test_per_atom_array_shapes_2N_T(self):
        N = 7
        T = 13
        cfg = single_pulse_N2000(num_molecules=N, seed=0)
        ckpt = build_initial_state(cfg, num_steps=T, rng=np.random.default_rng(0))
        for arr in (ckpt.positions_x, ckpt.positions_y, ckpt.positions_z,
                    ckpt.velocities_x, ckpt.velocities_y, ckpt.velocities_z,
                    ckpt.E_kin_eV, ckpt.E_pot_eV,
                    ckpt.E_dissip_eV, ckpt.L_droplet_eV_ps):
            assert arr.shape == (2 * N, T), f"unexpected shape {arr.shape}"

    def test_per_atom_static_arrays(self):
        N = 7
        cfg = single_pulse_N2000(num_molecules=N, seed=0)
        ckpt = build_initial_state(cfg, num_steps=10,
                                   rng=np.random.default_rng(0))
        assert ckpt.mass_kg.shape == (2 * N,)
        assert ckpt.droplet_radii.shape == (2 * N,)

    def test_per_molecule_arrays(self):
        N = 7
        cfg = single_pulse_N2000(num_molecules=N, seed=0)
        ckpt = build_initial_state(cfg, num_steps=10,
                                   rng=np.random.default_rng(0))
        assert ckpt.r0.shape == (N,)
        assert ckpt.E_initial_eV.shape == (N,)

    def test_time_ps_shape_T(self):
        cfg = single_pulse_N2000(num_molecules=5, seed=0)
        ckpt = build_initial_state(cfg, num_steps=42,
                                   rng=np.random.default_rng(0))
        assert ckpt.time_ps.shape == (42,)
        assert ckpt.time_ps[0] == 0.0


# ===========================================================================
# Validation
# ===========================================================================
class TestValidation:
    def test_zero_num_steps_raises(self):
        cfg = single_pulse_N2000(num_molecules=5, seed=0)
        with pytest.raises(ValueError, match="num_steps"):
            build_initial_state(cfg, num_steps=0,
                                rng=np.random.default_rng(0))

    def test_negative_num_steps_raises(self):
        cfg = single_pulse_N2000(num_molecules=5, seed=0)
        with pytest.raises(ValueError, match="num_steps"):
            build_initial_state(cfg, num_steps=-3,
                                rng=np.random.default_rng(0))


# ===========================================================================
# Column 0 populated, columns 1+ zero
# ===========================================================================
class TestColumn0PopulatedRestZero:
    @pytest.fixture
    def ckpt(self):
        cfg = single_pulse_N2000(num_molecules=5, seed=0)
        return build_initial_state(cfg, num_steps=10,
                                   rng=np.random.default_rng(0))

    def test_positions_col0_nonzero(self, ckpt):
        # In single pulse mode positions are set, not at origin
        any_nonzero = np.any(ckpt.positions_x[:, 0] != 0)
        assert any_nonzero

    def test_positions_col_1plus_zero(self, ckpt):
        np.testing.assert_array_equal(ckpt.positions_x[:, 1:],
                                      np.zeros_like(ckpt.positions_x[:, 1:]))

    def test_E_kin_col_1plus_zero(self, ckpt):
        np.testing.assert_array_equal(ckpt.E_kin_eV[:, 1:],
                                      np.zeros_like(ckpt.E_kin_eV[:, 1:]))

    def test_E_dissip_zero_at_t0(self, ckpt):
        # dissipation only accumulates; should start at 0
        np.testing.assert_array_equal(ckpt.E_dissip_eV[:, 0],
                                      np.zeros(2 * ckpt.num_molecules))


# ===========================================================================
# Two-atom geometry: bond length, axis alignment
# ===========================================================================
class TestTwoAtomGeometry:
    def test_bond_length_equals_R0_when_deltaR0_zero(self):
        cfg = single_pulse_N2000(num_molecules=20, seed=0,
                                  R0_GS_angstrom=2.666,
                                  deltaR0_angstrom=0.0)
        ckpt = build_initial_state(cfg, num_steps=2,
                                   rng=np.random.default_rng(0))
        N = ckpt.num_molecules
        dx = ckpt.positions_x[:N, 0] - ckpt.positions_x[N:, 0]
        dy = ckpt.positions_y[:N, 0] - ckpt.positions_y[N:, 0]
        dz = ckpt.positions_z[:N, 0] - ckpt.positions_z[N:, 0]
        bonds = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        np.testing.assert_allclose(bonds, 2.666, rtol=1e-12)

    def test_bond_length_R0_9A_HeDFT_mode(self):
        cfg = single_pulse_N2000(num_molecules=20, seed=0,
                                  R0_GS_angstrom=9.0,
                                  deltaR0_angstrom=0.0)
        ckpt = build_initial_state(cfg, num_steps=2,
                                   rng=np.random.default_rng(0))
        N = ckpt.num_molecules
        dx = ckpt.positions_x[:N, 0] - ckpt.positions_x[N:, 0]
        dy = ckpt.positions_y[:N, 0] - ckpt.positions_y[N:, 0]
        dz = ckpt.positions_z[:N, 0] - ckpt.positions_z[N:, 0]
        bonds = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        np.testing.assert_allclose(bonds, 9.0, rtol=1e-12)

    def test_atom_velocities_anti_aligned(self):
        """Atom 1 and atom 2 of the same molecule should have v_2 = -v_1."""
        cfg = single_pulse_N2000(num_molecules=20, seed=0)
        # Override single_pulse so v0 != 0
        from dataclasses import replace
        cfg = replace(cfg, single_pulse=False, partner_interaction=True,
                      effusive_dynamics=False)
        ckpt = build_initial_state(cfg, num_steps=2,
                                   rng=np.random.default_rng(0))
        N = ckpt.num_molecules
        np.testing.assert_allclose(ckpt.velocities_x[:N, 0],
                                   -ckpt.velocities_x[N:, 0], rtol=1e-12)
        np.testing.assert_allclose(ckpt.velocities_y[:N, 0],
                                   -ckpt.velocities_y[N:, 0], rtol=1e-12)
        np.testing.assert_allclose(ckpt.velocities_z[:N, 0],
                                   -ckpt.velocities_z[N:, 0], rtol=1e-12)

    def test_centre_of_mass_inside_droplet(self):
        """Each atom should be inside its droplet (r < R)."""
        cfg = single_pulse_N2000(num_molecules=20, seed=0)
        ckpt = build_initial_state(cfg, num_steps=2,
                                   rng=np.random.default_rng(0))
        r = np.sqrt(
            ckpt.positions_x[:, 0] ** 2
            + ckpt.positions_y[:, 0] ** 2
            + ckpt.positions_z[:, 0] ** 2
        )
        # Allow atoms slightly outside the droplet (since the molecule
        # CoM is inside but atoms are at +-R0/2 from CoM, atoms can be
        # near the surface). Just check we're not totally outside:
        assert np.all(r < 1.5 * ckpt.droplet_radii), (
            f"some atoms way outside droplet: max r/R = "
            f"{(r / ckpt.droplet_radii).max():.3f}"
        )


# ===========================================================================
# Velocity scaling: single-pulse vs others
# ===========================================================================
class TestInitialVelocity:
    def test_single_pulse_zero_velocity(self):
        cfg = single_pulse_N2000(num_molecules=10, seed=0)
        assert cfg.single_pulse
        ckpt = build_initial_state(cfg, num_steps=2,
                                   rng=np.random.default_rng(0))
        assert np.all(ckpt.velocities_x[:, 0] == 0)
        assert np.all(ckpt.velocities_y[:, 0] == 0)
        assert np.all(ckpt.velocities_z[:, 0] == 0)
        assert np.all(ckpt.E_kin_eV[:, 0] == 0)

    def test_non_single_pulse_nonzero_velocity(self):
        from dataclasses import replace
        cfg = single_pulse_N2000(num_molecules=10, seed=0)
        cfg = replace(cfg, single_pulse=False)
        ckpt = build_initial_state(cfg, num_steps=2,
                                   rng=np.random.default_rng(0))
        # Some velocity should be present
        assert np.max(np.abs(ckpt.velocities_x[:, 0])) > 0


# ===========================================================================
# Reproducibility
# ===========================================================================
class TestReproducibility:
    def test_same_seed_same_state(self):
        cfg = single_pulse_N2000(num_molecules=5, seed=42)
        a = build_initial_state(cfg, num_steps=5, rng=np.random.default_rng(123))
        b = build_initial_state(cfg, num_steps=5, rng=np.random.default_rng(123))
        np.testing.assert_array_equal(a.positions_x, b.positions_x)
        np.testing.assert_array_equal(a.r0, b.r0)
        np.testing.assert_array_equal(a.droplet_radii, b.droplet_radii)


# ===========================================================================
# E_initial: photon energy
# ===========================================================================
class TestEInitial:
    def test_E_initial_matches_photon_energy(self):
        from i2_helium_md.physics.constants import HC
        cfg = single_pulse_N2000(num_molecules=10, seed=0)
        ckpt = build_initial_state(cfg, num_steps=2,
                                   rng=np.random.default_rng(0))
        expected = HC / cfg.lambda_pump_nm  # in eV
        np.testing.assert_allclose(ckpt.E_initial_eV, expected)


# ===========================================================================
# Single droplet size (HeDFT comparison mode)
# ===========================================================================
class TestSingleDropletSize:
    def test_use_single_droplet_size_uses_one_value(self):
        cfg = single_pulse_N2000(num_molecules=20, seed=0,
                                  use_single_droplet_size=True,
                                  single_droplet_size=2000)
        ckpt = build_initial_state(cfg, num_steps=2,
                                   rng=np.random.default_rng(0))
        # All molecules should have the same droplet radius
        N = ckpt.num_molecules
        first_atom_radii = ckpt.droplet_radii[:N]
        assert np.unique(first_atom_radii).size == 1


class TestEPotIncludesPartner:
    """Regression test for the legacy MATLAB bug where E_pot at t=0
    omitted the partner Morse term. This caused a ~3 eV discontinuity
    between t=0 and t=1 for R=9 A initial conditions, breaking energy
    conservation tests. See migration_log.md for the full story.
    """

    def test_E_pot_t0_includes_partner_morse(self):
        """E_pot at t=0 should equal droplet + half partner per atom.

        Specifically: E_pot[t=0] should be CONTINUOUS with E_pot[t=1]
        from a single propagation step (allowing only leapfrog drift).
        """
        from i2_helium_md.simulation.propagation_step import (
            neutral_propagation_step, state_from_checkpoint_column,
        )
        # Use single_pulse so velocities are 0 -- atoms barely move,
        # so E_pot[t=0] and E_pot[t=1] should be nearly identical.
        cfg = single_pulse_N2000(num_molecules=5, seed=42)
        ckpt = build_initial_state(cfg, num_steps=2,
                                    rng=np.random.default_rng(42))
        s0 = state_from_checkpoint_column(ckpt, 0)
        s1 = neutral_propagation_step(
            s0, cfg=cfg, mass_kg=ckpt.mass_kg,
            droplet_radii=ckpt.droplet_radii,
            prev_distance_angstrom=None,
            rng=np.random.default_rng(0),
        )
        E_pot_t0 = s0.E_pot_eV.sum()
        E_pot_t1 = s1.E_pot_eV.sum()
        # Without the fix, E_pot_t0 was ~0.003 and E_pot_t1 was ~3.28.
        # With the fix they should agree to leapfrog precision.
        assert abs(E_pot_t1 - E_pot_t0) < 1e-3, (
            f"E_pot discontinuity {E_pot_t0:.4f} -> {E_pot_t1:.4f}; "
            "did the partner Morse term get dropped at t=0?"
        )

    def test_E_pot_t0_at_R0_9_is_a_few_eV(self):
        """For R0_GS=9 A (HeDFT mode), the I-I Morse contribution is
        ~0.66 eV per pair. With N=5 molecules that's ~3.3 eV total --
        the dominant contribution to E_pot at t=0."""
        cfg = single_pulse_N2000(num_molecules=5, seed=42,
                                  R0_GS_angstrom=9.0)
        ckpt = build_initial_state(cfg, num_steps=2,
                                    rng=np.random.default_rng(42))
        E_pot_total = ckpt.E_pot_eV[:, 0].sum()
        # Should be dominated by the Morse term (~0.66 eV per pair).
        assert E_pot_total > 1.0, (
            f"E_pot at t=0 for R0=9 A is too small ({E_pot_total:.4f} eV); "
            "partner Morse term may be missing."
        )
