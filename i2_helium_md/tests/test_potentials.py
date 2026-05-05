"""Tests for i2_helium_md/physics/potentials.py."""

import numpy as np
import pytest

from i2_helium_md import SimConfig, single_pulse_N2000
from i2_helium_md.physics.constants import EV, U
from i2_helium_md.physics.potentials import (
    I2_X_STATE,
    I2PLUS_IP_0,
    I2PLUS_IP_REL,
    I2PLUS_STATES,
    MorseParams,
    _gaussian,
    _morse_a,
    droplet_force,
    droplet_potential,
    morse_I2plus,
    morse_I2plus_state_select,
    morse_X,
)


# ---------------------------------------------------------------------------
# _gaussian  -- MATLAB g()
# ---------------------------------------------------------------------------
class TestGaussian:
    def test_peak_at_mu_equals_one(self):
        """g(mu, sig, x=mu) == 1 (unit height)."""
        assert _gaussian(mu=9.0, sig=0.3, x=np.array(9.0)) == pytest.approx(1.0)

    def test_symmetry(self):
        """g is even about mu."""
        left = _gaussian(9.0, 0.3, np.array(8.5))
        right = _gaussian(9.0, 0.3, np.array(9.5))
        assert left == pytest.approx(right)

    def test_width_one_sigma(self):
        """At x = mu + sig, value equals exp(-1/2)."""
        val = _gaussian(9.0, 0.3, np.array(9.3))
        assert val == pytest.approx(np.exp(-0.5))

    def test_vectorized(self):
        x = np.linspace(7.0, 11.0, 101)
        out = _gaussian(9.0, 0.3, x)
        assert out.shape == x.shape
        # max should be close to 1 and at the grid point closest to 9.0
        assert out.max() == pytest.approx(1.0, abs=1e-3)


# ---------------------------------------------------------------------------
# _morse_a  -- Morse width parameter from spectroscopic constants
# ---------------------------------------------------------------------------
class TestMorseA:
    def test_matlab_formula(self):
        """Recompute a for the X state, compare to a direct MATLAB-formula calc."""
        c_m_per_s = 299_792_458.0
        mu_kg = 127.0 / 2.0 * U
        D_e_J = 1.556 * EV
        expected = (
            214.5 * (c_m_per_s * 1e2) * 2.0 * np.pi
            / np.sqrt(2.0 * D_e_J / mu_kg)
            * 1e-10
        )
        got = _morse_a(D_e_eV=1.556, omega_e_cm=214.5)
        assert got == pytest.approx(expected, rel=1e-12)

    def test_a_is_reasonable(self):
        """a should be of order 1 - 3 1/Angstrom for I2."""
        a = I2_X_STATE.a
        assert 1.0 < a < 3.0


# ---------------------------------------------------------------------------
# droplet_potential
# ---------------------------------------------------------------------------
class TestDropletPotential:
    def test_zero_well_inside(self):
        """Deep inside droplet (r << 0), potential tends to 0."""
        V = droplet_potential(
            r=np.array(-100.0), steepness=14.2, binding_energy=0.027
        )
        assert V == pytest.approx(0.0, abs=1e-15)

    def test_full_depth_outside(self):
        """Far outside (r >> 0), potential tends to the asymptote."""
        V = droplet_potential(
            r=np.array(100.0), steepness=14.2, binding_energy=0.027
        )
        assert V == pytest.approx(0.027, rel=1e-10)

    def test_half_depth_at_offset(self):
        """At r = offset, erf(0) = 0 -> potential is half the asymptote."""
        V = droplet_potential(
            r=np.array(5.0), steepness=14.2, binding_energy=0.027, offset=5.0
        )
        assert V == pytest.approx(0.027 / 2.0)

    def test_matlab_direct_port(self):
        """Direct substitution: MATLAB `((erf((x-c)/s)+1)/2)*E`."""
        from scipy.special import erf
        r = np.linspace(-50, 50, 500)
        s, E, c = 14.2, 0.027, 0.0
        expected = ((erf((r - c) / s) + 1.0) / 2.0) * E
        got = droplet_potential(r, s, E, c)
        np.testing.assert_allclose(got, expected, rtol=1e-12)

    def test_force_matches_finite_difference(self):
        """Analytical force should match MATLAB's finite-difference form."""
        r = np.linspace(-30, 30, 300)
        s, E, c = 14.2, 0.027, 0.0
        h = 1e-6
        fd = (
            droplet_potential(r + h, s, E, c) - droplet_potential(r, s, E, c)
        ) / h
        analytic = droplet_force(r, s, E, c)
        # Finite-difference is O(h), our analytic derivative is exact.
        np.testing.assert_allclose(analytic, fd, atol=1e-5)

    def test_force_positive_at_surface(self):
        """Force dU/dr is positive near r=0 (potential rises outward)."""
        assert droplet_force(np.array(0.0), 14.2, 0.027) > 0.0


# ---------------------------------------------------------------------------
# morse_X
# ---------------------------------------------------------------------------
class TestMorseX:
    def test_minimum_at_R_e(self):
        """Morse has its minimum at R_e (no Xdip contamination)."""
        cfg = single_pulse_N2000(Xdip_active=False)
        R_e = I2_X_STATE.R_e
        r = np.linspace(R_e - 0.5, R_e + 0.5, 200)
        U = morse_X(r, cfg)
        # argmin should coincide with the grid point closest to R_e
        assert abs(r[np.argmin(U)] - R_e) < 0.01

    def test_well_depth_no_xdip(self):
        """At R_e, U = 0 (no Xdip)."""
        cfg = single_pulse_N2000(Xdip_active=False)
        U = morse_X(np.array(I2_X_STATE.R_e), cfg)
        assert U == pytest.approx(0.0, abs=1e-12)

    def test_asymptote_no_xdip(self):
        """At r >> R_e, U -> D_e."""
        cfg = single_pulse_N2000(Xdip_active=False)
        U = morse_X(np.array(50.0), cfg)
        assert U == pytest.approx(I2_X_STATE.D_e, rel=1e-8)

    def test_xdip_lowers_potential_at_9A(self):
        """With Xdip on, the potential at r=9A is lower by ~0.9 eV."""
        cfg_on = single_pulse_N2000(Xdip_active=True)
        cfg_off = single_pulse_N2000(Xdip_active=False)
        r = np.array(9.0)
        diff = morse_X(r, cfg_off) - morse_X(r, cfg_on)
        # peak of the gaussian is 1 at r=mu=9, so difference = 0.9 eV exactly
        assert diff == pytest.approx(0.9, rel=1e-10)

    def test_xdip_has_no_effect_far_from_9A(self):
        """At r=3A, gaussian at center 9A, sigma 0.3 is negligibly small."""
        cfg_on = single_pulse_N2000(Xdip_active=True)
        cfg_off = single_pulse_N2000(Xdip_active=False)
        r = np.array(3.0)
        diff = abs(morse_X(r, cfg_off) - morse_X(r, cfg_on))
        # at (9-3)/0.3 = 20 sigma, gaussian is ~ exp(-200) ~ 0
        assert diff < 1e-80


# ---------------------------------------------------------------------------
# morse_I2plus
# ---------------------------------------------------------------------------
class TestMorseI2plus:
    def test_minimum_at_R_e(self):
        for state in I2PLUS_STATES:
            r = np.linspace(state.R_e - 0.5, state.R_e + 0.5, 200)
            U = morse_I2plus(r, state)
            assert abs(r[np.argmin(U)] - state.R_e) < 0.01, (
                f"Morse minimum mispositioned for state {state}"
            )

    def test_asymptotic_depth(self):
        for state in I2PLUS_STATES:
            U = morse_I2plus(np.array(50.0), state)
            assert U == pytest.approx(state.D_e, rel=1e-6)


# ---------------------------------------------------------------------------
# morse_I2plus_state_select
# ---------------------------------------------------------------------------
class TestMorseI2plusStateSelect:
    def test_single_molecule_each_state(self):
        """Check one molecule per state: the reference at r=2.666 equals IP_0 + IP_rel."""
        r = np.full(4, 2.666)
        state_ids = np.array([0, 1, 2, 3])
        U = morse_I2plus_state_select(r, state_ids)
        # at r = 2.666, V_morse(r) - V_morse(2.666) = 0, so U = IP_0 + IP_rel
        expected = I2PLUS_IP_0 + np.array(I2PLUS_IP_REL)
        np.testing.assert_allclose(U, expected, rtol=1e-10)

    def test_vectorized_matches_scalar(self):
        """Running the vectorized function on N molecules equals N scalar calls."""
        rng = np.random.default_rng(0)
        N = 50
        r = rng.uniform(2.0, 6.0, N)
        state_ids = rng.integers(0, 4, N)

        vec = morse_I2plus_state_select(r, state_ids)

        # recompute scalar
        scalar = np.zeros(N)
        for i in range(N):
            s = I2PLUS_STATES[state_ids[i]]
            scalar[i] = (
                morse_I2plus(np.array(r[i]), s)
                - morse_I2plus(np.array(2.666), s)
                + I2PLUS_IP_0
                + I2PLUS_IP_REL[state_ids[i]]
            )
        np.testing.assert_allclose(vec, scalar, rtol=1e-12)

    def test_rejects_bad_state_id(self):
        r = np.array([2.666, 2.666])
        bad = np.array([0, 99])
        with pytest.raises(ValueError):
            morse_I2plus_state_select(r, bad)

    def test_rejects_shape_mismatch(self):
        with pytest.raises(ValueError):
            morse_I2plus_state_select(
                np.array([2.0, 3.0]), np.array([0, 1, 2])
            )


# ---------------------------------------------------------------------------
# Regression check vs MATLAB baseline values
# ---------------------------------------------------------------------------
class TestAgainstMatlabBaseline:
    """Hardcoded values cross-checked with the MATLAB formulas by hand.

    If any of these fail after future refactoring, something broke physics.
    """

    def test_morse_X_at_R_e_without_xdip_is_zero(self):
        cfg = SimConfig(Xdip_active=False)
        assert morse_X(np.array(2.666), cfg) == pytest.approx(0.0, abs=1e-12)

    def test_morse_X_at_infinity_without_xdip_is_D_e(self):
        cfg = SimConfig(Xdip_active=False)
        assert morse_X(np.array(1e6), cfg) == pytest.approx(1.556, rel=1e-6)

    def test_droplet_potential_at_origin_half_depth(self):
        V = droplet_potential(np.array(0.0), steepness=14.2, binding_energy=0.1)
        assert V == pytest.approx(0.05, rel=1e-10)
