"""Tests for i2_helium_md/physics/pickup.py (Tier-2 Phase-B Slice P).

The Poisson He-pickup gain channel: the per-step Bernoulli attach probability, the
Langmuir-capped rate, the momentum-conserving capture reset (delegated to
``mass_jump.capture``), and the S1 heat (delegated to
``internal_energy_budget.dE_int_pickup_eV``). Stochastic entry points take an injected
RNG; deterministic reset/heat composition is checked against the real neighbour modules,
statistical behaviour against Poisson/Bernoulli moment oracles with a seeded RNG.
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import MASS_HE_AMU, N_STAR
from i2_helium_md.physics.internal_energy_budget import dE_int_pickup_eV
from i2_helium_md.physics.mass_jump import capture
from i2_helium_md.physics.pickup import (
    PickupResult,
    attach_probability,
    lambda_attach,
    pickup_step,
    pickup_step_components,
)

V = np.array([3.0, -4.0, 12.0])       # |v| = 13 A/ps
M = 202.953908                        # a representative complex mass [amu]


class _ConstRNG:
    """Injected-RNG stub returning a constant uniform draw (scalar or array)."""

    def __init__(self, value: float):
        self.value = float(value)

    def random(self, size=None):
        if size is None:
            return self.value
        return np.full(size, self.value, dtype=float)


class _ArrayRNG:
    """Injected-RNG stub returning a preset draw array (size must match)."""

    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)

    def random(self, size=None):
        if size is None:
            return float(self.values.reshape(-1)[0])
        expected = int(np.prod(size))
        assert expected == self.values.size, f"size {size} != {self.values.size}"
        return self.values.reshape(size)


_ALWAYS = _ConstRNG(0.0)                 # draw 0 < any P>0 -> fire
_NEVER = _ConstRNG(0.9999999999)         # draw ~1, P<1 -> no fire


# ---------------------------------------------------------------------------
# attach_probability
# ---------------------------------------------------------------------------
class TestAttachProbability:
    def test_matches_closed_form(self):
        assert attach_probability(0.5, 0.02) == pytest.approx(1.0 - np.exp(-0.01))

    def test_rare_event_limit(self):
        # P_attach -> lambda*dt for small lambda*dt
        assert attach_probability(1e-4, 0.01) == pytest.approx(1e-6, rel=1e-3)

    def test_vectorized(self):
        lam = np.array([0.1, 1.0, 10.0])
        out = attach_probability(lam, 0.05)
        assert isinstance(out, np.ndarray)
        np.testing.assert_allclose(out, 1.0 - np.exp(-lam * 0.05))

    def test_monotone_and_bounded(self):
        lam = np.array([0.0, 0.1, 1.0, 10.0, 1e6])
        out = attach_probability(lam, 0.01)
        assert np.all(np.diff(out) >= 0.0)          # non-decreasing in lambda
        assert out[0] == 0.0                        # lambda=0 -> never fires
        assert np.all((0.0 <= out) & (out <= 1.0))

    def test_zero_dt_never_fires(self):
        assert attach_probability(5.0, 0.0) == 0.0


# ---------------------------------------------------------------------------
# lambda_attach
# ---------------------------------------------------------------------------
class TestLambdaAttach:
    def test_langmuir_p1_linear(self):
        lam = lambda_attach(0.5, 5, lambda0=0.8, n_star=21, p=1.0, cap="langmuir")
        assert lam == pytest.approx(0.8 * 0.5 * (1.0 - 5.0 / 21.0))

    def test_cap_none_is_density_only(self):
        lam = lambda_attach(0.5, 5, lambda0=0.8, cap="none")
        assert lam == pytest.approx(0.8 * 0.5)   # occupancy factor inert = 1

    def test_langmuir_and_none_agree_at_empty_shell(self):
        n0 = 0
        lam_l = lambda_attach(0.7, n0, lambda0=1.0, cap="langmuir")
        lam_n = lambda_attach(0.7, n0, lambda0=1.0, cap="none")
        assert lam_l == pytest.approx(lam_n)

    def test_cap_zero_at_full_shell(self):
        assert lambda_attach(1.0, N_STAR, lambda0=1.0, n_star=N_STAR) == pytest.approx(0.0)

    def test_cap_clamps_above_full_shell(self):
        # (1 - n/n*)_+ floors at 0 for n > n* (no negative rate)
        assert lambda_attach(1.0, N_STAR + 5, lambda0=1.0, n_star=N_STAR) == 0.0

    def test_p_greater_than_one_sharpens(self):
        base = 1.0 - 10.0 / 21.0
        lam = lambda_attach(1.0, 10, lambda0=1.0, n_star=21, p=2.0)
        assert lam == pytest.approx(base ** 2)

    def test_vectorized(self):
        rho = np.array([0.2, 0.5, 1.0])
        n = np.array([0, 10, 20])
        out = lambda_attach(rho, n, lambda0=1.0, n_star=21)
        assert isinstance(out, np.ndarray)
        np.testing.assert_allclose(out, rho * (1.0 - n / 21.0))

    def test_monotone_non_increasing_in_n(self):
        n = np.arange(0, 22)
        lam = lambda_attach(1.0, n, lambda0=1.0, n_star=21)
        assert np.all(np.diff(lam) <= 0.0)          # more He -> smaller pickup rate

    def test_zero_gate_gives_zero_rate(self):
        # rho_He/rho_bulk = 0 (ion outside the droplet) -> no pickup, regardless of n
        assert lambda_attach(0.0, 3, lambda0=1.0) == 0.0

    @pytest.mark.parametrize("form", ["sweeping", "dwell_time"])
    def test_unbuilt_rate_form_raises(self, form):
        with pytest.raises(NotImplementedError, match="rule-2"):
            lambda_attach(0.5, 3, lambda0=1.0, pickup_rate_form=form)

    def test_bad_cap_raises(self):
        with pytest.raises(ValueError, match="occupancy cap"):
            lambda_attach(0.5, 3, lambda0=1.0, cap="sigmoid")

    def test_negative_exponent_rejected(self):
        # p < 0 would diverge the Langmuir factor to inf at n >= n* (silent fire-always)
        with pytest.raises(ValueError, match="p must be >= 0"):
            lambda_attach(1.0, 3, lambda0=1.0, p=-1.0)

    def test_non_positive_n_star_rejected(self):
        with pytest.raises(ValueError, match="n_star must be > 0"):
            lambda_attach(1.0, 0, lambda0=1.0, n_star=0)

    def test_bad_knobs_ignored_when_cap_none(self):
        # p / n_star are unused under cap="none" -- must NOT raise
        assert lambda_attach(0.7, 3, lambda0=1.0, cap="none", p=-1.0, n_star=0) == \
            pytest.approx(0.7)

    def test_negative_rho_ratio_rejected(self):
        # Post-review fix (2026-07-02, Phase-B review): a negative density ratio gives
        # lambda < 0 -> P_attach < 0, silently disabling the channel -- the same class
        # as the lambda0 < 0 guard (helium_density can never produce one; this is
        # defense-in-depth against a bad caller).
        with pytest.raises(ValueError, match="rho_ratio must be >= 0"):
            lambda_attach(-0.1, 3, lambda0=1.0)

    def test_negative_rho_ratio_rejected_in_array(self):
        with pytest.raises(ValueError, match="rho_ratio must be >= 0"):
            lambda_attach(np.array([0.5, -0.1, 1.0]), np.array([3, 3, 3]), lambda0=1.0)

    def test_p_zero_makes_langmuir_cap_inert(self):
        # Characterization lock (post-review, 2026-07-02): p = 0 gives 0**0 == 1 even
        # at n >= n*, so the "cap" is structurally inert -- equivalent to cap="none"
        # with a langmuir label. Deliberately allowed (only p < 0 is rejected); locked
        # so any future tightening is a conscious choice, not silent drift.
        assert lambda_attach(1.0, N_STAR, lambda0=1.0, n_star=N_STAR, p=0.0) == \
            pytest.approx(1.0)
        assert lambda_attach(1.0, N_STAR + 5, lambda0=1.0, n_star=N_STAR, p=0.0) == \
            pytest.approx(1.0)


# ---------------------------------------------------------------------------
# pickup_step -- forced fire / no-fire
# ---------------------------------------------------------------------------
class TestPickupStepFire:
    def _fire(self, n=5, f_ret=0.3, kappa=1.0, picture="statistical_mixture"):
        return pickup_step(
            rng=_ALWAYS, n=n, v=V, m_amu=M, rho_ratio=1.0, lambda0=1.0,
            f_ret=f_ret, picture=picture, kappa=kappa, dt_ps=0.01,
        )

    def test_returns_pickupresult(self):
        assert isinstance(self._fire(), PickupResult)

    def test_fire_increments_n_and_mass(self):
        res = self._fire(n=5)
        assert res.fired is True
        assert res.n_plus == 6
        assert res.m_plus_amu == pytest.approx(M + MASS_HE_AMU)

    def test_fire_reset_matches_capture(self):
        res = self._fire()
        cap = capture(V, M)
        np.testing.assert_allclose(res.v_plus, cap.v_plus, atol=1e-13)
        assert res.dE_mass_transfer == pytest.approx(cap.dE_mass_transfer, rel=0, abs=1e-12)
        assert res.dE_mass_transfer > 0.0

    def test_fire_s1_heat_matches_U(self):
        # n = pre-pickup; picture/kappa threaded to dE_int_pickup_eV = +f_ret*D0(n+1)
        res = self._fire(n=5, f_ret=0.3, kappa=1.0, picture="statistical_mixture")
        expected = dE_int_pickup_eV(5, f_ret=0.3, picture="statistical_mixture", kappa=1.0)
        assert res.dE_int_eV == pytest.approx(expected)
        assert res.dE_int_eV > 0.0

    def test_s1_heat_picture_kappa_threaded(self):
        r_mix = self._fire(picture="statistical_mixture")
        r_x2 = self._fire(picture="x2_only")
        assert r_mix.dE_int_eV != r_x2.dE_int_eV     # picture actually threaded


class TestPickupStepNoFire:
    def test_no_fire_is_identity(self):
        res = pickup_step(
            rng=_NEVER, n=5, v=V, m_amu=M, rho_ratio=1.0, lambda0=1.0,
            f_ret=0.3, kappa=1.0, dt_ps=0.01,
        )
        assert res.fired is False
        assert res.n_plus == 5
        assert res.m_plus_amu == pytest.approx(M)
        np.testing.assert_array_equal(res.v_plus, V)
        assert res.dE_int_eV == 0.0
        assert res.dE_mass_transfer == 0.0

    def test_no_fire_velocity_is_independent_copy(self):
        v = V.copy()
        res = pickup_step(
            rng=_NEVER, n=5, v=v, m_amu=M, rho_ratio=1.0, lambda0=1.0,
            f_ret=0.3, kappa=1.0, dt_ps=0.01,
        )
        res.v_plus[0] = 999.0
        assert v[0] == V[0]

    def test_closed_gate_never_fires_even_with_forcing_rng(self):
        # rho_He/rho_bulk = 0 (ion exits the droplet) -> P_attach = 0 -> no pickup,
        # even though _ALWAYS draws 0.0. This is the physical exit-termination property.
        res = pickup_step(
            rng=_ALWAYS, n=5, v=V, m_amu=M, rho_ratio=0.0, lambda0=1.0,
            f_ret=0.3, kappa=1.0, dt_ps=0.01,
        )
        assert res.fired is False
        assert res.n_plus == 5

    def test_full_shell_never_fires_under_langmuir(self):
        # n = n* -> Langmuir factor 0 -> P_attach = 0, even with the forcing RNG
        res = pickup_step(
            rng=_ALWAYS, n=N_STAR, v=V, m_amu=M, rho_ratio=1.0, lambda0=1.0,
            f_ret=0.3, kappa=1.0, dt_ps=0.01, n_star=N_STAR, cap="langmuir",
        )
        assert res.fired is False

    def test_full_shell_still_fires_under_cap_none(self):
        # cap="none" removes the saturation -> a full shell can still capture
        res = pickup_step(
            rng=_ALWAYS, n=N_STAR, v=V, m_amu=M, rho_ratio=1.0, lambda0=1.0,
            f_ret=0.3, kappa=1.0, dt_ps=0.01, n_star=N_STAR, cap="none",
        )
        assert res.fired is True
        assert res.n_plus == N_STAR + 1


class TestPickupStepUnbuiltArms:
    def test_thermal_capture_raises(self):
        with pytest.raises(NotImplementedError, match="thermal"):
            pickup_step(
                rng=_ALWAYS, n=5, v=V, m_amu=M, rho_ratio=1.0, lambda0=1.0,
                f_ret=0.3, kappa=1.0, dt_ps=0.01, he_capture_velocity="thermal",
            )

    def test_sweeping_rate_form_raises(self):
        with pytest.raises(NotImplementedError, match="rule-2"):
            pickup_step(
                rng=_ALWAYS, n=5, v=V, m_amu=M, rho_ratio=1.0, lambda0=1.0,
                f_ret=0.3, kappa=1.0, dt_ps=0.01, pickup_rate_form="sweeping",
            )


# ---------------------------------------------------------------------------
# Statistical: Poisson / Bernoulli moments (seeded RNG, sample-size tolerance)
# ---------------------------------------------------------------------------
class TestPoissonMoments:
    def test_fire_count_matches_lambda_t(self):
        rng = np.random.default_rng(20260701)
        lam, dt, steps = 0.5, 0.02, 40000
        n_fires = 0
        for _ in range(steps):
            res = pickup_step(
                rng=rng, n=0, v=V, m_amu=M, rho_ratio=1.0, lambda0=lam,
                f_ret=0.3, kappa=1.0, dt_ps=dt, cap="none",   # constant rate: n fixed at 0
            )
            n_fires += int(res.fired)
        p = 1.0 - np.exp(-lam * dt)
        mean = steps * p                       # ~ lambda*t = 400
        std = np.sqrt(steps * p * (1.0 - p))
        assert abs(n_fires - mean) < 5.0 * std   # ~5-sigma band


class TestEnsembleIndependence:
    def test_fraction_fired_matches_probability(self):
        rng = np.random.default_rng(7)
        n_ions = 200_000
        lam, dt = 0.8, 0.05
        zeros = np.zeros(n_ions)
        n = np.zeros(n_ions, dtype=int)
        _, _, _, _, _, _, _, fired = pickup_step_components(
            rng=rng, n=n, vx=zeros + 1.0, vy=zeros, vz=zeros,
            m_amu=zeros + M, rho_ratio=zeros + 1.0, lambda0=lam,
            f_ret=0.3, kappa=1.0, dt_ps=dt, cap="none",
        )
        p = 1.0 - np.exp(-lam * dt)
        frac = fired.mean()
        std = np.sqrt(p * (1.0 - p) / n_ions)
        assert abs(frac - p) < 5.0 * std

    def test_disjoint_halves_have_consistent_rates(self):
        rng = np.random.default_rng(11)
        n_ions = 100_000
        lam, dt = 1.0, 0.03
        zeros = np.zeros(n_ions)
        _, _, _, _, _, _, _, fired = pickup_step_components(
            rng=rng, n=np.zeros(n_ions, dtype=int), vx=zeros + 1.0, vy=zeros,
            vz=zeros, m_amu=zeros + M, rho_ratio=zeros + 1.0, lambda0=lam,
            f_ret=0.3, kappa=1.0, dt_ps=dt, cap="none",
        )
        half = n_ions // 2
        r1, r2 = fired[:half].mean(), fired[half:].mean()
        # 5-sigma band on the difference of two independent Bernoulli means
        # (post-review fix 2026-07-02: the former 0.02 was ~17 sigma, nearly vacuous).
        p = 1.0 - np.exp(-lam * dt)
        se_diff = np.sqrt(2.0 * p * (1.0 - p) / half)
        assert abs(r1 - r2) < 5.0 * se_diff


# ---------------------------------------------------------------------------
# RNG-consumption parity with a REAL generator (post-review lock, 2026-07-02)
# ---------------------------------------------------------------------------
class TestRNGConsumptionParity:
    def test_scalar_and_components_consume_identical_stream(self):
        # Locks the Q5 / Slice-X draw-order contract with a real PCG64 generator (the
        # stub-based TestComponentsMatchScalar checks outputs but sidesteps actual RNG
        # consumption): M scalar pickup_step calls must consume the stream exactly like
        # one pickup_step_components call -- including the suppressed ion (P_attach = 0
        # at n = n*), which still consumes its draw. A refactor that skips the draw
        # when P == 0 would break this test while leaving every stub test green.
        n = np.array([0, N_STAR, 5, 10])          # ion 1: full shell -> P_attach = 0
        vx = np.array([1.0, 3.0, -2.0, 0.5])
        vy = np.array([0.0, -4.0, 1.0, 2.0])
        vz = np.array([0.0, 12.0, 2.0, -1.0])
        m = np.array([120.0, M, 300.0, 180.0])
        rho = np.array([1.0, 1.0, 0.5, 0.9])
        kw = dict(lambda0=14.0, f_ret=0.3, kappa=1.0, dt_ps=0.05)  # P ~ 0.5 in-gate

        rng_c = np.random.default_rng(20260702)
        n_plus, m_plus, vx_p, vy_p, vz_p, dE_int, dE_mt, fired = pickup_step_components(
            rng=rng_c, n=n, vx=vx, vy=vy, vz=vz, m_amu=m, rho_ratio=rho, **kw,
        )
        assert not fired[1]                        # the full-shell ion cannot fire

        rng_s = np.random.default_rng(20260702)
        for i in range(4):
            res = pickup_step(
                rng=rng_s, n=int(n[i]), v=np.array([vx[i], vy[i], vz[i]]),
                m_amu=float(m[i]), rho_ratio=float(rho[i]), **kw,
            )
            assert res.fired == bool(fired[i])
            assert res.n_plus == n_plus[i]
            assert res.m_plus_amu == pytest.approx(m_plus[i])
            assert res.dE_int_eV == pytest.approx(dE_int[i])
            assert res.dE_mass_transfer == pytest.approx(dE_mt[i])
        # Identical post-call stream state: the very next uniform must agree bitwise.
        assert rng_s.random() == rng_c.random()


# ---------------------------------------------------------------------------
# Components form vs the scalar single-ion oracle
# ---------------------------------------------------------------------------
class TestComponentsMatchScalar:
    def test_all_fire_matches_scalar(self):
        n = np.array([0, 5, 10])
        vx = np.array([1.0, 3.0, -2.0])
        vy = np.array([0.0, -4.0, 1.0])
        vz = np.array([0.0, 12.0, 2.0])
        m = np.array([120.0, M, 300.0])
        rho = np.array([1.0, 0.5, 0.9])
        out = pickup_step_components(
            rng=_ConstRNG(0.0), n=n, vx=vx, vy=vy, vz=vz, m_amu=m, rho_ratio=rho,
            lambda0=1.0, f_ret=0.3, kappa=1.0, dt_ps=0.01,
        )
        n_plus, m_plus, vx_p, vy_p, vz_p, dE_int, dE_mt, fired = out
        assert fired.all()
        for i in range(3):
            res = pickup_step(
                rng=_ConstRNG(0.0), n=int(n[i]), v=np.array([vx[i], vy[i], vz[i]]),
                m_amu=float(m[i]), rho_ratio=float(rho[i]), lambda0=1.0,
                f_ret=0.3, kappa=1.0, dt_ps=0.01,
            )
            assert n_plus[i] == res.n_plus
            assert m_plus[i] == pytest.approx(res.m_plus_amu)
            np.testing.assert_allclose(
                [vx_p[i], vy_p[i], vz_p[i]], res.v_plus, atol=1e-12
            )
            assert dE_int[i] == pytest.approx(res.dE_int_eV)
            assert dE_mt[i] == pytest.approx(res.dE_mass_transfer, rel=0, abs=1e-12)

    def test_mixed_fire_pattern_leaves_nonfired_untouched(self):
        n = np.array([2, 4, 6])
        vx = np.array([1.0, 2.0, 3.0])
        vy = np.zeros(3)
        vz = np.zeros(3)
        m = np.array([120.0, 150.0, 200.0])
        rho = np.ones(3)
        # per-ion draws: ions 0 and 2 fire (draw 0), ion 1 does not (draw ~1)
        rng = _ArrayRNG([0.0, 0.9999999999, 0.0])
        out = pickup_step_components(
            rng=rng, n=n, vx=vx, vy=vy, vz=vz, m_amu=m, rho_ratio=rho,
            lambda0=1.0, f_ret=0.3, kappa=1.0, dt_ps=0.01,
        )
        n_plus, m_plus, vx_p, _, _, dE_int, dE_mt, fired = out
        np.testing.assert_array_equal(fired, [True, False, True])
        # non-fired ion 1 untouched
        assert n_plus[1] == 4
        assert m_plus[1] == pytest.approx(150.0)
        assert vx_p[1] == pytest.approx(2.0)
        assert dE_int[1] == 0.0
        assert dE_mt[1] == 0.0
        # fired ions changed
        assert n_plus[0] == 3 and n_plus[2] == 7
        assert m_plus[0] == pytest.approx(120.0 + MASS_HE_AMU)
        assert dE_int[0] > 0.0 and dE_mt[0] > 0.0

    def test_none_fire_is_identity(self):
        n = np.array([2, 4])
        vx = np.array([1.0, 2.0])
        out = pickup_step_components(
            rng=_ConstRNG(0.9999999999), n=n, vx=vx, vy=np.zeros(2), vz=np.zeros(2),
            m_amu=np.array([120.0, 150.0]), rho_ratio=np.ones(2), lambda0=1.0,
            f_ret=0.3, kappa=1.0, dt_ps=0.01,
        )
        n_plus, m_plus, vx_p, _, _, dE_int, dE_mt, fired = out
        assert not fired.any()
        np.testing.assert_array_equal(n_plus, n)
        np.testing.assert_array_equal(vx_p, vx)
        np.testing.assert_array_equal(dE_int, np.zeros(2))
        np.testing.assert_array_equal(dE_mt, np.zeros(2))
