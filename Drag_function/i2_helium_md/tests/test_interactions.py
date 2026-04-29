"""Tests for i2_helium_md/physics/interactions.py."""

import numpy as np
import pytest

from i2_helium_md import SimConfig, single_pulse_N2000
from i2_helium_md.physics.constants import MASS_I_AMU, U
from i2_helium_md.physics.interactions import (
    _pair_geometry,
    _split_pair_coordinates,
    atom_interaction_potential,
    ion_interaction_potential,
    partner_interaction_ion,
    partner_interaction_neutral,
)
from i2_helium_md.physics.constants import EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2
from i2_helium_md.physics.potentials import I2_X_STATE, morse_X


# ---------------------------------------------------------------------------
# Conversion factor (derivation cross-check)
# ---------------------------------------------------------------------------
class TestConversionFactor:
    def test_factor_value(self):
        """Force-to-acceleration factor for kg masses: ≈ 1.602176634e-23.

        Derivation:
            a [A/ps^2] = F [eV/A] / m [kg] * EV * 1e-4
                       = F [eV/A] / m [kg] * 1.602176634e-23
        """
        assert EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2 == pytest.approx(
            1.602176634e-23, rel=1e-9,
        )

    def test_factor_matches_alternative_form(self):
        """Cross-check vs the equivalent 'eV/(A*u)' formulation: u * 9648.5."""
        from i2_helium_md.physics.constants import U
        alt = U * 9648.533
        assert EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2 == pytest.approx(alt, rel=1e-3)


# ---------------------------------------------------------------------------
# Pair geometry helpers
# ---------------------------------------------------------------------------
class TestPairGeometry:
    def test_split_pair_coordinates(self):
        x = np.array([1.0, 2.0, 3.0, 10.0, 20.0, 30.0])
        y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        z = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        r1x, r1y, r1z, r2x, r2y, r2z = _split_pair_coordinates(x, y, z)
        np.testing.assert_array_equal(r1x, [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(r2x, [10.0, 20.0, 30.0])

    def test_split_rejects_odd_length(self):
        with pytest.raises(ValueError):
            _split_pair_coordinates(np.array([1., 2., 3.]),
                                    np.array([1., 2., 3.]),
                                    np.array([1., 2., 3.]))

    def test_pair_geometry_simple_case(self):
        # two molecules: pair 0 along x (atom1 at 2.666, atom2 at 0),
        #                pair 1 along y (atom1 at 0, atom2 at -3)
        x = np.array([2.666, 0.0,    0.0,  0.0])
        y = np.array([0.0,   0.0,    0.0,  3.0])
        z = np.zeros(4)
        dr, dr_unit = _pair_geometry(x, y, z)

        np.testing.assert_allclose(dr, [2.666, 3.0])
        np.testing.assert_allclose(dr_unit[0], [1.0, 0.0, 0.0])
        np.testing.assert_allclose(dr_unit[1], [0.0, -1.0, 0.0])


# ---------------------------------------------------------------------------
# atom_interaction_potential == morse_X (delegation)
# ---------------------------------------------------------------------------
class TestAtomInteractionPotential:
    def test_delegates_to_morse_X(self):
        cfg = single_pulse_N2000()
        r = np.linspace(2.0, 10.0, 50)
        np.testing.assert_array_equal(
            atom_interaction_potential(r, cfg),
            morse_X(r, cfg),
        )


# ---------------------------------------------------------------------------
# ion_interaction_potential
# ---------------------------------------------------------------------------
class TestIonInteractionPotential:
    def test_pure_coulomb_doubly_charged(self):
        """q1 = q2 = 1: should give exactly E_coulomb_scale * 14.4 / r."""
        cfg = single_pulse_N2000()
        r = np.array([2.0, 5.0, 10.0])
        q1 = np.array([1, 1, 1])
        q2 = np.array([1, 1, 1])
        U = ion_interaction_potential(r, q1, q2, cfg)
        expected = cfg.E_coulomb_scale * 14.39964548 / r
        np.testing.assert_allclose(U, expected, rtol=1e-12)

    def test_coulomb_scale_knob(self):
        """E_coulomb_scale rescales the result linearly."""
        r = np.array([3.0])
        q = np.array([1])
        U_full = ion_interaction_potential(r, q, q, single_pulse_N2000(E_coulomb_scale=1.0))
        U_half = ion_interaction_potential(r, q, q, single_pulse_N2000(E_coulomb_scale=0.5))
        np.testing.assert_allclose(U_half * 2.0, U_full, rtol=1e-12)

    def test_zero_charge_gives_zero_coulomb(self):
        """q1 = q2 = 0 gives U = 0 (when single_charge is off)."""
        cfg = single_pulse_N2000()
        r = np.array([3.0, 5.0])
        z = np.zeros(2, dtype=int)
        U = ion_interaction_potential(r, z, z, cfg)
        np.testing.assert_allclose(U, 0.0)

    def test_single_charge_allowed_requires_state_ids(self):
        cfg = single_pulse_N2000(single_charge_ionization_allowed=True)
        r = np.array([2.666])
        with pytest.raises(ValueError):
            ion_interaction_potential(r, np.array([1]), np.array([0]), cfg)

    def test_single_charge_mode_adds_morse_on_singly_ionized(self):
        """For (q1,q2)=(1,0), U = Coulomb(1*0=0) + Morse = Morse only."""
        cfg = single_pulse_N2000(single_charge_ionization_allowed=True)
        r = np.array([2.666])
        q1 = np.array([1])
        q2 = np.array([0])
        state_ids = np.array([0])
        U = ion_interaction_potential(r, q1, q2, cfg, state_ids=state_ids)
        # Coulomb term = 1*0*14.4/2.666 = 0
        # Morse term at R_e gives 0 (minimum), then shifted by IP_0 + IP_rel[0] = 9.36
        assert U[0] == pytest.approx(9.36, rel=1e-10)


# ---------------------------------------------------------------------------
# partner_interaction_neutral
# ---------------------------------------------------------------------------
class TestPartnerInteractionNeutral:
    def _simple_setup(self, N=1, R=2.666):
        """Create N molecules aligned along x, separated by R."""
        # atom 1 of each molecule at +R/2, atom 2 at -R/2
        x = np.concatenate([np.full(N, R / 2), np.full(N, -R / 2)])
        y = np.zeros(2 * N)
        z = np.zeros(2 * N)
        mass = np.full(2 * N, MASS_I_AMU * U)  # kg
        return x, y, z, mass

    def test_shapes(self):
        x, y, z, mass = self._simple_setup(N=5, R=2.666)
        cfg = single_pulse_N2000()
        ax, ay, az, E_pot = partner_interaction_neutral(x, y, z, mass, cfg)
        assert ax.shape == (10,)
        assert ay.shape == (10,)
        assert az.shape == (10,)
        assert E_pot.shape == (5,)

    def test_at_equilibrium_force_is_zero(self):
        """At R=R_e without Xdip, Morse has a minimum -> force is near-zero.

        Note: we use finite differences (h=1e-4) to match MATLAB, so the
        measured force at a minimum is *not* exactly zero but proportional
        to h * U''(R_e). For I2 X state U''(R_e) ~ 11 eV/A^2, giving
        |a| ~ 0.04 A/ps^2 at h=1e-4. We check it's small, not exactly zero.
        """
        R_e = I2_X_STATE.R_e
        x, y, z, mass = self._simple_setup(N=1, R=R_e)
        cfg = single_pulse_N2000(Xdip_active=False)
        ax, ay, az, E_pot = partner_interaction_neutral(x, y, z, mass, cfg)
        # |ax| should be at most a few tenths of A/ps^2 from FD truncation.
        # Compare to the "real" force at R_e+-0.1 which is ~60-100 A/ps^2
        # -> FD residual is at least 100x smaller, so <1 A/ps^2 is a solid check.
        assert abs(ax[0]) < 1.0
        assert abs(ay[0]) < 1e-12
        assert abs(az[0]) < 1e-12

    def test_force_repulsive_near_R_e(self):
        """At R < R_e, Morse is repulsive -> atom 1 accelerates in +x (away from atom 2)."""
        R_e = I2_X_STATE.R_e
        x, y, z, mass = self._simple_setup(N=1, R=R_e - 0.1)
        cfg = single_pulse_N2000(Xdip_active=False)
        ax, ay, az, E_pot = partner_interaction_neutral(x, y, z, mass, cfg)
        # atom 1 (at +R/2) should feel +x force
        assert ax[0] > 0
        # atom 2 (at -R/2) should feel equal and opposite force
        assert ax[1] == pytest.approx(-ax[0], rel=1e-10)

    def test_force_attractive_just_beyond_R_e(self):
        """At R > R_e (just), Morse is attractive -> atom 1 pulled toward atom 2."""
        R_e = I2_X_STATE.R_e
        x, y, z, mass = self._simple_setup(N=1, R=R_e + 0.1)
        cfg = single_pulse_N2000(Xdip_active=False)
        ax, ay, az, _ = partner_interaction_neutral(x, y, z, mass, cfg)
        assert ax[0] < 0   # atom 1 pulled in -x
        assert ax[1] == pytest.approx(-ax[0], rel=1e-10)

    def test_newtons_third_law(self):
        """For any configuration, atom 1 and 2 get equal and opposite accelerations
        (because masses are equal here).
        """
        x, y, z, mass = self._simple_setup(N=3, R=4.0)
        cfg = single_pulse_N2000()
        ax, ay, az, _ = partner_interaction_neutral(x, y, z, mass, cfg)
        # atoms 0,1,2 are atom-1 of each pair; atoms 3,4,5 are atom-2
        np.testing.assert_allclose(ax[:3], -ax[3:], rtol=1e-10)
        np.testing.assert_allclose(ay[:3], -ay[3:], rtol=1e-10)
        np.testing.assert_allclose(az[:3], -az[3:], rtol=1e-10)


# ---------------------------------------------------------------------------
# partner_interaction_ion
# ---------------------------------------------------------------------------
class TestPartnerInteractionIon:
    def _simple_setup(self, N=1, R=2.666):
        x = np.concatenate([np.full(N, R / 2), np.full(N, -R / 2)])
        y = np.zeros(2 * N)
        z = np.zeros(2 * N)
        mass = np.full(2 * N, MASS_I_AMU * U)
        charge = np.ones(2 * N, dtype=int)
        return x, y, z, mass, charge

    def test_coulomb_repulsion_pushes_ions_apart(self):
        R = 3.0
        x, y, z, mass, charge = self._simple_setup(N=1, R=R)
        cfg = single_pulse_N2000()
        ax, ay, az, E_pot = partner_interaction_ion(x, y, z, mass, charge, cfg)
        # atom 1 at +R/2 should feel +x force (pushed away from atom 2)
        assert ax[0] > 0
        # symmetry
        assert ax[1] == pytest.approx(-ax[0], rel=1e-10)
        # no transverse components
        assert abs(ay[0]) < 1e-10
        assert abs(az[0]) < 1e-10

    def test_E_pot_per_atom_halved(self):
        """MATLAB convention: energy per atom is pair energy / 2."""
        R = 2.666
        x, y, z, mass, charge = self._simple_setup(N=1, R=R)
        cfg = single_pulse_N2000()
        _, _, _, E_pot = partner_interaction_ion(x, y, z, mass, charge, cfg)
        # Pair Coulomb U = 14.4/R; returned E_pot is U/2 per atom.
        expected_half = (cfg.E_coulomb_scale * 14.39964548 / R) / 2.0
        np.testing.assert_allclose(E_pot, [expected_half, expected_half], rtol=1e-10)

    def test_coulomb_scaling_affects_force(self):
        R = 3.0
        x, y, z, mass, charge = self._simple_setup(N=1, R=R)
        cfg_full = single_pulse_N2000(E_coulomb_scale=1.0)
        cfg_half = single_pulse_N2000(E_coulomb_scale=0.5)
        ax_full, _, _, _ = partner_interaction_ion(x, y, z, mass, charge, cfg_full)
        ax_half, _, _, _ = partner_interaction_ion(x, y, z, mass, charge, cfg_half)
        # Halving the Coulomb scale should halve the force
        np.testing.assert_allclose(ax_half * 2, ax_full, rtol=1e-10)

    def test_shapes(self):
        x, y, z, mass, charge = self._simple_setup(N=4, R=3.0)
        cfg = single_pulse_N2000()
        ax, ay, az, E_pot = partner_interaction_ion(x, y, z, mass, charge, cfg)
        assert ax.shape == (8,)
        assert E_pot.shape == (8,)   # per-atom layout


# ---------------------------------------------------------------------------
# Regression / MATLAB comparison
# ---------------------------------------------------------------------------
class TestAgainstMatlabConventions:
    """Checks that preserve MATLAB behaviour byte-for-byte."""

    def test_finite_difference_step_matches(self):
        """MATLAB uses h=1e-4 in the force derivative; we use the same default.

        This is a documentation test rather than a code test -- if someone
        changes the default h, they should consciously update this.
        """
        # Inspect the FD function's signature
        import inspect
        from i2_helium_md.physics import interactions
        sig = inspect.signature(interactions._force_from_potential_fd)
        assert sig.parameters["h"].default == 1e-4

    def test_sign_convention(self):
        """F = (U(r) - U(r+h))/h  (gives +F for repulsive short-range)."""
        from i2_helium_md.physics.interactions import _force_from_potential_fd

        def U(r):
            # simple repulsive 1/r potential
            return 1.0 / r

        F = _force_from_potential_fd(U, np.array([1.0]))
        # at r=1, dU/dr = -1, so F = -dU/dr = +1
        # force is +1 (repulsive) -- matches MATLAB sign convention
        assert F[0] > 0
