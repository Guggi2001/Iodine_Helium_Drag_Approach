"""Tests for i2_helium_md/physics/baoab.py (Slice 2: BAOAB ion-stage stepper).

The killer test (:meth:`TestAnchor.test_zero_gamma_matches_velocity_verlet`) is
the safety proof for the whole slice: with drag off (``gamma == 0``) and noise
off, the BAOAB step must recover the frozen baseline ``velocity_verlet_step`` on
the *same* ``acc_fn`` to round-off. That simultaneously proves (a) the BAOAB
scheme is a correct integrator and (b) the ``_kick``/``_drift`` extraction left
``velocity_verlet_step`` behaviourally unchanged.

All quantities here are mechanical and mass-in-amu (the module's locked unit
contract, ``SLICE2_GOALS_baoab_ion_stepper.md`` §2): velocities in A/ps, masses
in amu, gamma in amu/ps, dissipated energy in amu*A^2/ps^2. No kg, no eV.
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.baoab import make_ion_baoab_step
from i2_helium_md.physics.drag import LINEAR_CUBIC, DragCoefficients, drag_gamma
from i2_helium_md.physics.leapfrog import velocity_verlet_step


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _zero_gamma(v, depth):
    """gamma_fn that is identically zero (drag off)."""
    return np.zeros_like(np.asarray(v, dtype=float))


def _const_gamma(value):
    """gamma_fn returning a constant amu/ps regardless of v / depth."""

    def gamma_fn(v, depth):
        return np.full_like(np.asarray(v, dtype=float), float(value))

    return gamma_fn


def _zero_accel(p):
    """Conservative acceleration that is identically zero (free particle)."""
    z = np.zeros_like(p[0])
    return (z, z.copy(), z.copy()), np.zeros(p[0].shape[0] // 2)


def _harmonic_accel(k):
    """Position-dependent (isotropic harmonic) acceleration, a = -k*pos.

    Position-dependent so the two per-step ``acc_fn`` evaluations (at x0 and x1)
    differ -- the anchor test must exercise a non-trivial force.
    """

    def acc_fn(p):
        x, y, z = p
        E_pot = 0.5 * k * (x ** 2 + y ** 2 + z ** 2)
        N = x.shape[0] // 2
        E_pot_per_pair = E_pot[:N] + E_pot[N:]
        return (-k * x, -k * y, -k * z), E_pot_per_pair

    return acc_fn


class _CallCounter:
    """Wrap an ``acc_fn`` and count how many times it is evaluated.

    Lets the anchor test assert that BAOAB evaluates the conservative force
    *eval-for-eval* like ``velocity_verlet_step`` (exactly twice per step),
    not merely that outputs agree.
    """

    def __init__(self, fn):
        self._fn = fn
        self.calls = 0

    def __call__(self, p):
        self.calls += 1
        return self._fn(p)


def _sample_state(n_pairs=2):
    """A reproducible multi-atom (x, y, z), (vx, vy, vz) state, length 2N each."""
    rng = np.random.default_rng(12345)
    n = 2 * n_pairs
    pos = (rng.uniform(-2.0, 2.0, n), rng.uniform(-2.0, 2.0, n), rng.uniform(-2.0, 2.0, n))
    vel = (rng.uniform(-3.0, 3.0, n), rng.uniform(-3.0, 3.0, n), rng.uniform(-3.0, 3.0, n))
    return pos, vel


# ===========================================================================
# Anchor (killer) test
# ===========================================================================
class TestAnchor:
    def test_zero_gamma_matches_velocity_verlet(self):
        """gamma=0, noise off: BAOAB step == baseline velocity_verlet_step.

        Run many steps with a non-trivial position-dependent force and confirm
        positions, velocities, and E_pot agree to round-off.
        """
        pos0, vel0 = _sample_state()
        n = 8
        m = np.full(2 * 2, 203.0)            # amu, value irrelevant at gamma=0
        droplet_radii = np.full(2 * 2, 1.0e6)  # huge -> depth deep, irrelevant
        dt = 0.01
        # Independent counters over fresh (identical, pure) force instances:
        # the BAOAB and Verlet loops run sequentially below.
        baoab_acc = _CallCounter(_harmonic_accel(k=5.0))
        verlet_acc = _CallCounter(_harmonic_accel(k=5.0))

        baoab_step = make_ion_baoab_step(m, droplet_radii, baoab_acc, _zero_gamma)

        p, v = pos0, vel0
        for _ in range(n):
            p, v, E_pot_b, dE = baoab_step(p, v, dt)

        pp, vv = pos0, vel0
        for _ in range(n):
            pp, vv, E_pot_v = velocity_verlet_step(pp, vv, verlet_acc, dt)

        for a, b in zip(p, pp):
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-14)
        for a, b in zip(v, vv):
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(E_pot_b, E_pot_v, rtol=1e-12, atol=1e-14)
        # No drag -> no dissipation.
        np.testing.assert_allclose(dE, 0.0, atol=1e-14)

        # Eval-for-eval parity (load-bearing, not decorative): BAOAB evaluates
        # acc_fn exactly twice per step (entry + post-step position), the same
        # count velocity_verlet_step makes. This is the *present no-cache*
        # contract; when the planned final-B -> next-first-B caching lands
        # (~n+1 evals), this assertion is the intended trip-wire and updates
        # with it.
        assert baoab_acc.calls == 2 * n
        assert baoab_acc.calls == verlet_acc.calls


# ===========================================================================
# Analytic decay (constant linear gamma, no conservative force)
# ===========================================================================
class TestAnalyticDecay:
    def test_constant_gamma_exponential_decay(self):
        """With a=0 conservative force and constant gamma, v(t)=v0*exp(-gamma t/m).

        The BAOAB O-step is exact for constant gamma, so this holds to round-off
        over many steps (not just dt-order).
        """
        pos0, vel0 = _sample_state()
        m_val = 203.0
        m = np.full(4, m_val)
        droplet_radii = np.full(4, 1.0e6)
        gamma = 7.0  # amu/ps
        dt = 0.01
        n = 50

        step = make_ion_baoab_step(m, droplet_radii, _zero_accel, _const_gamma(gamma))
        p, v = pos0, vel0
        for _ in range(n):
            p, v, _, _ = step(p, v, dt)

        decay = np.exp(-gamma * (n * dt) / m_val)
        for vi, v0i in zip(v, vel0):
            np.testing.assert_allclose(vi, decay * v0i, rtol=1e-12, atol=1e-14)

    def test_mass_enters_only_through_o_step_exponent(self):
        """Doubling m halves the damping rate gamma/m and nothing else.

        With no conservative force, one step gives vel = exp(-gamma dt/m) v0;
        the conservative ``acc_fn`` output is never scaled by m inside the step.
        """
        pos0, vel0 = _sample_state()
        gamma = 5.0
        dt = 0.02
        for m_val in (100.0, 200.0):
            m = np.full(4, m_val)
            droplet_radii = np.full(4, 1.0e6)
            step = make_ion_baoab_step(m, droplet_radii, _zero_accel, _const_gamma(gamma))
            _, v, _, _ = step(pos0, vel0, dt)
            decay = np.exp(-gamma * dt / m_val)
            for vi, v0i in zip(v, vel0):
                np.testing.assert_allclose(vi, decay * v0i, rtol=1e-12, atol=1e-14)


# ===========================================================================
# Dissipated-energy bookkeeping
# ===========================================================================
class TestDissipatedEnergy:
    def test_dissipation_identity_and_units(self):
        """Returned dE == 0.5*m*(|v_in|^2 - |v_out|^2), per atom, in amu*A^2/ps^2.

        Use a=0 conservative force so v_in^O == v0 (the first B is a no-op) and
        v_out is the post-step velocity, making the identity checkable from the
        outside.
        """
        pos0, vel0 = _sample_state()
        m_val = 203.0
        m = np.full(4, m_val)
        droplet_radii = np.full(4, 1.0e6)
        gamma = 9.0
        dt = 0.03

        step = make_ion_baoab_step(m, droplet_radii, _zero_accel, _const_gamma(gamma))
        _, v_out, _, dE = step(pos0, vel0, dt)

        speed_in_sq = vel0[0] ** 2 + vel0[1] ** 2 + vel0[2] ** 2
        speed_out_sq = v_out[0] ** 2 + v_out[1] ** 2 + v_out[2] ** 2
        expected = 0.5 * m_val * (speed_in_sq - speed_out_sq)

        assert dE.shape == (4,)
        # The identity itself is the units guard: ``expected`` is in raw
        # amu*A^2/ps^2, so an accidental eV conversion (~1e-4 of this) would
        # break the rtol=1e-12 match below.
        np.testing.assert_allclose(dE, expected, rtol=1e-12, atol=1e-14)
        # Magnitude sanity: amu*A^2/ps^2 lands at O(0.1-10) here; the eV
        # equivalent would be ~1e-4 of this. Deterministic (fixed seed).
        assert np.all(dE > 0.1)

    def test_dissipativity_never_increases_speed(self):
        """For gamma>=0 and no conservative force, the O-step never raises |v|."""
        pos0, vel0 = _sample_state()
        m = np.full(4, 203.0)
        droplet_radii = np.full(4, 1.0e6)
        step = make_ion_baoab_step(m, droplet_radii, _zero_accel, _const_gamma(6.0))
        _, v_out, _, dE = step(pos0, vel0, 0.05)

        speed_in = np.sqrt(vel0[0] ** 2 + vel0[1] ** 2 + vel0[2] ** 2)
        speed_out = np.sqrt(v_out[0] ** 2 + v_out[1] ** 2 + v_out[2] ** 2)
        assert np.all(speed_out <= speed_in + 1e-14)
        assert np.all(dE >= 0.0)


# ===========================================================================
# Gate-off limit (real drag law, deep vacuum)
# ===========================================================================
class TestGateOff:
    def test_vacuum_is_identity_with_real_drag_law(self):
        """Deep in vacuum (depth >> 0) the erf gate -> 0, so drag -> 0.

        Wire the real ``drag_gamma`` (linear_cubic, a,b>0) and place atoms far
        outside the droplet: the O-step must be the identity and dissipate
        nothing.
        """
        coeffs = DragCoefficients(
            form=LINEAR_CUBIC,
            coefficients={"a": 3.0, "b": 0.5},
            extraction_mass_model="constant",
            extraction_mass_amu=203.0,
        )
        steepness = 14.2

        def gamma_fn(v, depth):
            return drag_gamma(v, depth, coeffs, steepness)

        m = np.full(4, 203.0)
        droplet_radii = np.full(4, 10.0)
        # Atoms ~1000 A out -> depth ~ +990 A >> steepness -> gate ~ 0.
        pos0 = (np.full(4, 1000.0), np.zeros(4), np.zeros(4))
        vel0 = (np.full(4, 4.0), np.full(4, -2.0), np.zeros(4))

        step = make_ion_baoab_step(m, droplet_radii, _zero_accel, gamma_fn)
        _, v_out, _, dE = step(pos0, vel0, 0.01)

        for vi, v0i in zip(v_out, vel0):
            np.testing.assert_allclose(vi, v0i, rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(dE, 0.0, atol=1e-14)


# ===========================================================================
# Noise dormant at Tier 0
# ===========================================================================
class TestNoiseDormant:
    def test_T_eff_zero_draws_no_rng_and_is_deterministic(self):
        """T_eff=0: passing an rng changes nothing and consumes no draws."""
        pos0, vel0 = _sample_state()
        m = np.full(4, 203.0)
        droplet_radii = np.full(4, 1.0e6)
        gamma_fn = _const_gamma(5.0)
        dt = 0.01

        # rng=None path
        step_none = make_ion_baoab_step(m, droplet_radii, _zero_accel, gamma_fn)
        p1, v1, _, dE1 = step_none(pos0, vel0, dt)

        # rng provided, T_eff still 0 -> must not be touched
        rng = np.random.default_rng(999)
        state_before = rng.bit_generator.state
        step_rng = make_ion_baoab_step(
            m, droplet_radii, _zero_accel, gamma_fn, T_eff=0.0, rng=rng
        )
        p2, v2, _, dE2 = step_rng(pos0, vel0, dt)
        state_after = rng.bit_generator.state

        # The dormancy proof is the *state* equality (no draw consumed), not
        # the output equality below -- outputs are identical at T_eff=0
        # regardless of whether an rng was touched.
        assert state_before == state_after
        for a, b in zip(v1, v2):
            np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(dE1, dE2)

    def test_active_noise_is_not_yet_implemented(self):
        """T_eff>0 must fail loudly: active Langevin noise is Slice >=3."""
        m = np.full(4, 203.0)
        droplet_radii = np.full(4, 1.0e6)
        with pytest.raises(NotImplementedError):
            make_ion_baoab_step(
                m, droplet_radii, _zero_accel, _const_gamma(5.0),
                T_eff=1.0, rng=np.random.default_rng(0),
            )
