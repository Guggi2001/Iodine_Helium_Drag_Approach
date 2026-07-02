"""Tests for i2_helium_md/physics/evaporation.py (Tier-2 Phase-B Slice Q).

The energy-gated, RRK-rate-limited He-evaporation loss channel. Deterministic
oracles: the mode-count ``effective_dof(n)``, the saturating RRK rate ``k`` with the
self-bound gate + the ``n=1`` direct-dissociation branch, the signed gate margin, and
the cold-shed reset (delegated to the *real* ``mass_jump.cold_shed``) + the K1 drain
(the *real* ``internal_energy_budget.dE_int_shed_eV``). Stochastic entry points take an
injected RNG; the scalar step is the single-ion oracle for the vectorized components form.

Physics direction (MASS §4 / R9)
--------------------------------
Shedding happens in the band ``D_0(n) < E_int < Sigma(n)``: the upper bound is the
self-bound gate (``is_self_bound`` == ``E_int < Sigma(n)``, True = self-bound = shed
enabled), the lower bound is RRK-bracket positivity. ``n=1`` is the degenerate diatomic
(the band collapses), so it is the *direct-dissociation* special case with the inverted
gate ``E_int > D_0(1)``, ``k = nu``.
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import MASS_HE_AMU, N_STAR, NU_EVAP_PER_PS
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum
from i2_helium_md.physics.internal_energy_budget import dE_int_shed_eV
from i2_helium_md.physics.mass_jump import cold_shed
from i2_helium_md.physics.evaporation import (
    EvaporationResult,
    effective_dof,
    rrk_rate,
    shed_probability,
    is_self_bound,
    gate_margin_eV,
    evaporation_step,
    evaporation_step_components,
)

V = np.array([3.0, -4.0, 12.0])       # |v| = 13 A/ps
M = 202.953908                        # a representative I+He_n complex mass [amu]
KAPPA = 1.0
PIC = "statistical_mixture"


class _ConstRNG:
    """Injected-RNG stub returning a constant uniform draw (scalar or array)."""

    def __init__(self, value: float):
        self.value = float(value)
        self.n_calls = 0
        self.n_draws = 0

    def random(self, size=None):
        self.n_calls += 1
        if size is None:
            self.n_draws += 1
            return self.value
        self.n_draws += int(np.prod(size))
        return np.full(size, self.value, dtype=float)


class _ArrayRNG:
    """Injected-RNG stub returning a preset draw array (size must match)."""

    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)
        self.n_calls = 0

    def random(self, size=None):
        self.n_calls += 1
        if size is None:
            return float(self.values.reshape(-1)[0])
        expected = int(np.prod(size))
        assert expected == self.values.size, f"size {size} != {self.values.size}"
        return self.values.reshape(size)


_NEVER = _ConstRNG(0.9999999999)         # draw ~1, P<1 -> no fire


def _sigma(n):
    return ladder_cumsum(n, picture=PIC, kappa=KAPPA)


def _d0(n):
    return d0_of_n(n, picture=PIC, kappa=KAPPA)


# ---------------------------------------------------------------------------
# effective_dof
# ---------------------------------------------------------------------------
class TestEffectiveDof:
    def test_n2_is_four_linear(self):
        # I+He2 = ion + 2 He = 3-atom *linear* complex: 3N-5 = 4 (CALIBRATION row 10).
        assert effective_dof(2) == 4

    def test_n3_is_six(self):
        assert effective_dof(3) == 6         # nonlinear 3n-3

    def test_nstar_is_sixty(self):
        assert effective_dof(21) == 60       # 3*21-3

    def test_vectorized(self):
        out = effective_dof(np.array([2, 3, 21]))
        assert isinstance(out, np.ndarray)
        np.testing.assert_array_equal(out, [4, 6, 60])

    def test_raises_on_n1_direct_branch(self):
        # n=1 is the direct-dissociation branch (diatomic; 3n-3=0 degenerate).
        with pytest.raises(ValueError, match="n >= 2"):
            effective_dof(1)

    def test_raises_on_n0(self):
        with pytest.raises(ValueError, match="n >= 2"):
            effective_dof(0)

    def test_raises_on_array_with_bad_element(self):
        with pytest.raises(ValueError, match="n >= 2"):
            effective_dof(np.array([2, 1, 3]))

    def test_fractional_n_rejected_scalar(self):
        # Post-review fix (2026-07-02): a fractional n silently truncated the scalar
        # path (int(4.5) = 4) while the array path kept 4.5 -- now both fail loud,
        # mirroring the ladder_cumsum integer-occupancy contract.
        with pytest.raises(ValueError, match="integer occupancy"):
            effective_dof(2.5)

    def test_fractional_n_rejected_array(self):
        with pytest.raises(ValueError, match="integer occupancy"):
            effective_dof(np.array([2.5, 3.0]))

    def test_integer_valued_float_accepted(self):
        # Integer-valued floats cast cleanly (the ladder_cumsum contract).
        assert effective_dof(3.0) == 6


# ---------------------------------------------------------------------------
# shed_probability
# ---------------------------------------------------------------------------
class TestShedProbability:
    def test_matches_closed_form(self):
        assert shed_probability(0.5, 0.02) == pytest.approx(1.0 - np.exp(-0.01))

    def test_rare_event_limit(self):
        assert shed_probability(1e-4, 0.01) == pytest.approx(1e-6, rel=1e-3)

    def test_vectorized(self):
        k = np.array([0.1, 1.0, 10.0])
        out = shed_probability(k, 0.05)
        assert isinstance(out, np.ndarray)
        np.testing.assert_allclose(out, 1.0 - np.exp(-k * 0.05))

    def test_monotone_and_bounded(self):
        k = np.array([0.0, 0.1, 1.0, 10.0, 1e6])
        out = shed_probability(k, 0.01)
        assert np.all(np.diff(out) >= 0.0)
        assert out[0] == 0.0
        assert np.all((0.0 <= out) & (out <= 1.0))

    def test_zero_dt_never_fires(self):
        assert shed_probability(5.0, 0.0) == 0.0


# ---------------------------------------------------------------------------
# rrk_rate  (n >= 2 saturating RRK; n = 1 direct)
# ---------------------------------------------------------------------------
class TestRRKRate:
    def test_zero_below_threshold(self):
        # E_int <= D_0(n): cannot afford the top rung -> k = 0 exactly.
        n = 5
        E = _d0(n)                                   # exactly at threshold
        assert rrk_rate(E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0
        assert rrk_rate(0.5 * E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0

    def test_positive_in_band_and_bounded(self):
        n = 10
        E = 0.5 * (_d0(n) + _sigma(n))               # midpoint of (D_0, Sigma)
        k = rrk_rate(E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)
        assert 0.0 < k < NU_EVAP_PER_PS

    def test_monotone_increasing_in_E_within_band(self):
        n = 10
        lo, hi = _d0(n), _sigma(n)
        Es = np.linspace(lo + 1e-4 * (hi - lo), hi - 1e-4 * (hi - lo), 25)
        ks = rrk_rate(Es, n * np.ones_like(Es, dtype=int),
                      nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)
        assert np.all(np.diff(ks) > 0.0)

    def test_self_bound_gate_suppresses(self):
        # E_int >= Sigma(n): net self-unbound -> shedding suppressed -> k = 0.
        n = 8
        assert rrk_rate(_sigma(n), n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0
        assert rrk_rate(2.0 * _sigma(n), n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0

    def test_no_avalanche_saturation(self):
        # Even as E_int -> just below Sigma, k stays < nu (bounded, no gate-open burst).
        n = 15
        E = _sigma(n) * (1.0 - 1e-9)
        k = rrk_rate(E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)
        assert 0.0 <= k < NU_EVAP_PER_PS

    def test_n1_direct_fires_when_hot(self):
        # n=1: direct dissociation k = nu once E_int > D_0(1) (the INVERTED gate).
        E = 2.0 * _d0(1)
        assert rrk_rate(E, 1, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == pytest.approx(
            NU_EVAP_PER_PS
        )

    def test_n1_direct_zero_when_cold(self):
        E = 0.5 * _d0(1)
        assert rrk_rate(E, 1, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0

    def test_vectorized_over_E_and_n(self):
        n = np.array([1, 5, 8])
        E = np.array([2.0 * _d0(1), 0.5 * (_d0(5) + _sigma(5)), 2.0 * _sigma(8)])
        k = rrk_rate(E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)
        assert isinstance(k, np.ndarray)
        assert k[0] == pytest.approx(NU_EVAP_PER_PS)     # n=1 hot -> nu
        assert 0.0 < k[1] < NU_EVAP_PER_PS               # n=5 in band
        assert k[2] == 0.0                               # n=8 suppressed

    def test_dof_override_applies_to_n_ge2(self):
        n = 10
        E = 0.5 * (_d0(n) + _sigma(n))
        base = 1.0 - _d0(n) / E
        k = rrk_rate(E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA, evap_rrk_dof=4.0)
        assert k == pytest.approx(NU_EVAP_PER_PS * base ** (4.0 - 1.0))

    def test_dof_override_leaves_n1_direct(self):
        E = 2.0 * _d0(1)
        k = rrk_rate(E, 1, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA, evap_rrk_dof=5.0)
        assert k == pytest.approx(NU_EVAP_PER_PS)         # n=1 stays direct

    def test_dof_override_below_one_rejected(self):
        # Module-level defense (mirrors pickup's p<0 guard); s<1 diverges the rate.
        with pytest.raises(ValueError, match="s >= 1"):
            rrk_rate(0.05, 10, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA, evap_rrk_dof=0.5)

    def test_gate_onset_override_shifts_gate(self):
        # A forced (larger) threshold keeps a would-be-suppressed E in the shed band.
        n = 8
        E = _sigma(n) * 1.5                               # suppressed under Sigma(n)
        assert rrk_rate(E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0
        k = rrk_rate(E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                     gate_onset_eV=2.0 * E)
        assert 0.0 < k < NU_EVAP_PER_PS


# ---------------------------------------------------------------------------
# gate diagnostics
# ---------------------------------------------------------------------------
class TestGate:
    def test_gate_margin_value(self):
        n = 7
        E = 0.05
        assert gate_margin_eV(E, n, picture=PIC, kappa=KAPPA) == pytest.approx(E - _sigma(n))

    def test_gate_margin_sign(self):
        n = 7
        assert gate_margin_eV(0.5 * _sigma(n), n, picture=PIC, kappa=KAPPA) < 0.0   # self-bound
        assert gate_margin_eV(2.0 * _sigma(n), n, picture=PIC, kappa=KAPPA) > 0.0   # self-unbound

    def test_is_self_bound_direction(self):
        n = 7
        assert is_self_bound(0.5 * _sigma(n), n, picture=PIC, kappa=KAPPA) is True   # E < Sigma
        assert is_self_bound(2.0 * _sigma(n), n, picture=PIC, kappa=KAPPA) is False  # E > Sigma

    def test_is_self_bound_flips_at_threshold(self):
        n = 7
        s = _sigma(n)
        assert is_self_bound(s * (1.0 - 1e-9), n, picture=PIC, kappa=KAPPA) is True
        assert is_self_bound(s * (1.0 + 1e-9), n, picture=PIC, kappa=KAPPA) is False

    def test_gate_onset_override(self):
        n = 7
        E = 2.0 * _sigma(n)                               # self-unbound under Sigma(n)
        assert is_self_bound(E, n, picture=PIC, kappa=KAPPA) is False
        assert is_self_bound(E, n, picture=PIC, kappa=KAPPA, gate_onset_eV=3.0 * E) is True

    def test_vectorized(self):
        n = np.array([5, 7])
        E = np.array([0.01, 0.5])
        g = gate_margin_eV(E, n, picture=PIC, kappa=KAPPA)
        assert isinstance(g, np.ndarray)
        np.testing.assert_allclose(g, E - _sigma(n))
        sb = is_self_bound(E, n, picture=PIC, kappa=KAPPA)
        assert sb.dtype == bool


# ---------------------------------------------------------------------------
# evaporation_step (scalar single-ion oracle)
# ---------------------------------------------------------------------------
def _step(rng, E, n, v=V, m=M, **kw):
    return evaporation_step(
        rng=rng, E_int_eV=E, n=n, v=v, m_amu=m,
        nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA, dt_ps=0.01, **kw,
    )


class TestEvaporationStep:
    def test_fire_path_matches_cold_shed_and_K1(self):
        n = 10
        E = 0.5 * (_d0(n) + _sigma(n))
        res = _step(_ConstRNG(0.0), E, n)               # fire-forcing
        assert isinstance(res, EvaporationResult)
        assert res.fired is True
        assert res.n_plus == n - 1
        expected = cold_shed(V, M, m_he_amu=MASS_HE_AMU)
        assert res.m_plus_amu == pytest.approx(expected.m_plus_amu)
        np.testing.assert_allclose(res.v_plus, expected.v_plus)
        assert res.dE_mass_transfer == pytest.approx(expected.dE_mass_transfer)
        assert res.dE_mass_transfer < 0.0                # cold-shed defect is negative
        assert res.dE_int_eV == pytest.approx(
            float(dE_int_shed_eV(n, picture=PIC, kappa=KAPPA))
        )
        assert res.dE_int_eV < 0.0                       # K1 drain -D_0(n)

    def test_no_fire_is_identity(self):
        n = 10
        E = 0.5 * (_d0(n) + _sigma(n))
        res = _step(_NEVER, E, n)
        assert res.fired is False
        assert res.n_plus == n
        assert res.m_plus_amu == pytest.approx(M)
        np.testing.assert_allclose(res.v_plus, V)
        assert res.dE_int_eV == 0.0
        assert res.dE_mass_transfer == 0.0

    def test_suppressed_never_fires(self):
        # Self-unbound (E > Sigma): no fire even with a fire-forcing RNG.
        n = 8
        res = _step(_ConstRNG(0.0), 2.0 * _sigma(n), n)
        assert res.fired is False
        assert res.n_plus == n

    def test_below_threshold_never_fires(self):
        n = 8
        res = _step(_ConstRNG(0.0), 0.5 * _d0(n), n)
        assert res.fired is False

    def test_n1_direct_fires_when_hot(self):
        res = _step(_ConstRNG(0.0), 2.0 * _d0(1), 1)
        assert res.fired is True
        assert res.n_plus == 0

    def test_unconditional_draw_even_when_suppressed(self):
        # Decision Q5: one draw consumed every step (k=0 under suppression -> no fire).
        rng = _ConstRNG(0.0)
        _step(rng, 2.0 * _sigma(8), 8)                  # suppressed
        assert rng.n_draws == 1


# ---------------------------------------------------------------------------
# evaporation_step_components (vectorized ensemble form)
# ---------------------------------------------------------------------------
def _step_components(rng, E, n, vx, vy, vz, m, **kw):
    return evaporation_step_components(
        rng=rng, E_int_eV=E, n=n, vx=vx, vy=vy, vz=vz, m_amu=m,
        nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA, dt_ps=0.01, **kw,
    )


class TestEvaporationStepComponents:
    def _ensemble(self):
        n = np.array([10, 10, 8])                        # three in-band ions (n=10,10,8)
        E = np.array([
            0.5 * (_d0(10) + _sigma(10)),                # in band
            0.5 * (_d0(10) + _sigma(10)),                # in band
            0.5 * (_d0(8) + _sigma(8)),                  # in band
        ])
        vx = np.array([3.0, 1.0, -2.0])
        vy = np.array([-4.0, 0.0, 5.0])
        vz = np.array([12.0, -1.0, 0.0])
        m = np.array([M, M + 4.0026, M - 4.0026])        # per-ion masses (fires diverge)
        return n, E, vx, vy, vz, m

    def test_all_fire_matches_scalar_oracle_ion_by_ion(self):
        n, E, vx, vy, vz, m = self._ensemble()
        (n_p, m_p, vxp, vyp, vzp, dEi, dEmt, fired) = _step_components(
            _ConstRNG(0.0), E, n, vx, vy, vz, m,
        )
        assert np.all(fired)
        for i in range(3):
            ref = evaporation_step(
                rng=_ConstRNG(0.0), E_int_eV=E[i], n=int(n[i]),
                v=np.array([vx[i], vy[i], vz[i]]), m_amu=m[i],
                nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA, dt_ps=0.01,
            )
            assert n_p[i] == ref.n_plus
            assert m_p[i] == pytest.approx(ref.m_plus_amu)
            np.testing.assert_allclose([vxp[i], vyp[i], vzp[i]], ref.v_plus)
            assert dEi[i] == pytest.approx(ref.dE_int_eV)
            assert dEmt[i] == pytest.approx(ref.dE_mass_transfer)

    def test_non_fired_ions_untouched(self):
        n, E, vx, vy, vz, m = self._ensemble()
        # Mixed draws: fire ion 0, not ions 1,2.
        (n_p, m_p, vxp, vyp, vzp, dEi, dEmt, fired) = _step_components(
            _ArrayRNG([0.0, 0.999999, 0.999999]), E, n, vx, vy, vz, m,
        )
        assert fired[0] and not fired[1] and not fired[2]
        # non-fired keep everything
        np.testing.assert_array_equal([n_p[1], n_p[2]], [n[1], n[2]])
        np.testing.assert_allclose([m_p[1], m_p[2]], [m[1], m[2]])
        np.testing.assert_allclose([vxp[1], vxp[2]], [vx[1], vx[2]])
        assert dEi[1] == 0.0 and dEmt[2] == 0.0

    def test_suppressed_ion_never_fires(self):
        n = np.array([8])
        E = np.array([2.0 * _sigma(8)])                  # self-unbound
        (n_p, *_rest, fired) = _step_components(
            _ConstRNG(0.0), E, n, np.array([3.0]), np.array([4.0]), np.array([0.0]),
            np.array([M]),
        )
        assert not fired[0]
        assert n_p[0] == 8

    def test_single_vectorized_draw_per_step(self):
        n, E, vx, vy, vz, m = self._ensemble()
        rng = _ConstRNG(0.0)
        _step_components(rng, E, n, vx, vy, vz, m)
        assert rng.n_calls == 1                          # one vectorized draw
        assert rng.n_draws == 3

    def test_rng_consumption_parity_with_scalar(self):
        # Decision Q5: scalar and components consume the RNG identically (one draw/ion/step),
        # so the scalar stays the exact ion-by-ion oracle even for suppressed ions.
        n, E, vx, vy, vz, m = self._ensemble()
        rng_c = _ConstRNG(0.0)
        _step_components(rng_c, E, n, vx, vy, vz, m)
        scalar_draws = 0
        for i in range(3):
            rng_s = _ConstRNG(0.0)
            evaporation_step(
                rng=rng_s, E_int_eV=E[i], n=int(n[i]),
                v=np.array([vx[i], vy[i], vz[i]]), m_amu=m[i],
                nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA, dt_ps=0.01,
            )
            scalar_draws += rng_s.n_draws
        assert rng_c.n_draws == scalar_draws == 3


# ---------------------------------------------------------------------------
# statistical / moment checks (seeded RNG)
# ---------------------------------------------------------------------------
class TestBernoulliMoments:
    def test_fire_fraction_matches_P_shed(self):
        # Over many seeded ions at a fixed (E, n), the fired fraction ~ P_shed.
        n = 10
        E = 0.5 * (_d0(n) + _sigma(n))
        k = rrk_rate(E, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)
        p = shed_probability(k, 0.01)
        M_ions = 200_000
        rng = np.random.default_rng(20260701)
        na = np.full(M_ions, n)
        Ea = np.full(M_ions, E)
        z = np.zeros(M_ions)
        (*_ignore, fired) = _step_components(rng, Ea, na, z + 3.0, z - 4.0, z + 12.0,
                                             np.full(M_ions, M))
        frac = fired.mean()
        se = np.sqrt(p * (1.0 - p) / M_ions)
        assert abs(frac - p) < 5.0 * se


# ---------------------------------------------------------------------------
# Review + hardening locks (Slice Q review pass)
# ---------------------------------------------------------------------------
class TestHardening:
    def test_effective_dof_array_dtype(self):
        out = effective_dof(np.array([2, 3, 5, 21]))
        assert out.dtype.kind in ("i", "u")
        np.testing.assert_array_equal(out, [4, 6, 12, 60])

    def test_shed_probability_saturates_to_one(self):
        # Large k*dt underflows exp -> exactly 1.0 (fire-every-step limit), like
        # pickup.attach_probability. Bounded, non-decreasing.
        assert shed_probability(1e6, 1.0) == 1.0

    def test_rrk_rate_strictly_below_nu_across_band(self):
        # No avalanche: k < nu for every finite E in the self-bound band, all n>=2.
        for n in (2, 3, 7, 15, 21):
            lo, hi = _d0(n), _sigma(n)
            Es = np.linspace(lo + 1e-6 * (hi - lo), hi - 1e-9 * (hi - lo), 50)
            ks = rrk_rate(Es, n * np.ones_like(Es, dtype=int),
                          nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)
            assert np.all((0.0 <= ks) & (ks < NU_EVAP_PER_PS))

    def test_gate_onset_override_does_not_touch_n1(self):
        # The override shifts the n>=2 self-bound gate only; n=1 keeps its D_0(1) direct
        # onset regardless of gate_onset_eV.
        E = 2.0 * _d0(1)
        k_default = rrk_rate(E, 1, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)
        k_override = rrk_rate(E, 1, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                              gate_onset_eV=1e-9)
        assert k_default == k_override == pytest.approx(NU_EVAP_PER_PS)

    def test_no_fire_returns_defensive_copy(self):
        # A no-fire result must not alias the caller's velocity buffer.
        v = np.array([1.0, 2.0, 3.0])
        res = _step(_NEVER, 0.5 * (_d0(10) + _sigma(10)), 10, v=v)
        v[0] = 999.0
        assert res.v_plus[0] == 1.0

    def test_components_threads_dof_and_gate_override(self):
        # evap_rrk_dof + gate_onset_eV reach the components form identically to the scalar.
        n = np.array([8, 8])
        E = np.array([1.5 * _sigma(8), 1.5 * _sigma(8)])   # suppressed under Sigma(8)
        vx, vy, vz = np.array([3.0, 3.0]), np.array([-4.0, -4.0]), np.array([12.0, 12.0])
        m = np.array([M, M])
        (_np, _mp, vxp, *_rest, fired) = _step_components(
            _ConstRNG(0.0), E, n, vx, vy, vz, m,
            evap_rrk_dof=4.0, gate_onset_eV=2.0 * float(E[0]),
        )
        assert np.all(fired)                       # override opens the gate
        ref = evaporation_step(
            rng=_ConstRNG(0.0), E_int_eV=float(E[0]), n=8,
            v=np.array([3.0, -4.0, 12.0]), m_amu=M,
            nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA, dt_ps=0.01,
            evap_rrk_dof=4.0, gate_onset_eV=2.0 * float(E[0]),
        )
        assert ref.fired is True
        np.testing.assert_allclose(vxp[0], ref.v_plus[0])

    def test_cold_shed_components_per_ion_mass_matches_scalar(self):
        # Lock the mass_jump generalization: per-ion-mass cold_shed_velocity_components
        # agrees ion-by-ion with the scalar cold_shed (the single-atom oracle).
        from i2_helium_md.physics.mass_jump import (
            cold_shed,
            cold_shed_velocity_components,
        )
        vx = np.array([3.0, 1.0, -2.0])
        vy = np.array([-4.0, 0.0, 5.0])
        vz = np.array([12.0, -1.0, 0.0])
        m = np.array([202.0, 150.0, 130.0])        # per-ion masses
        vxp, vyp, vzp, mp, dE = cold_shed_velocity_components(vx, vy, vz, m)
        for i in range(3):
            ref = cold_shed(np.array([vx[i], vy[i], vz[i]]), m[i])
            np.testing.assert_allclose([vxp[i], vyp[i], vzp[i]], ref.v_plus)
            assert mp[i] == pytest.approx(ref.m_plus_amu)
            assert dE[i] == pytest.approx(ref.dE_mass_transfer)

    def test_cold_shed_components_rejects_underweight_ion(self):
        # Array-aware shed guard: any ion with m <= m_He fails loud (m+ would be <= 0).
        from i2_helium_md.physics.mass_jump import cold_shed_velocity_components
        with pytest.raises(ValueError, match="m\\+ = m - m_He"):
            cold_shed_velocity_components(
                np.array([1.0, 2.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]),
                np.array([200.0, 2.0]),         # second ion lighter than He
            )


# ---------------------------------------------------------------------------
# Extended review pass (deeper edge / parity / robustness locks)
# ---------------------------------------------------------------------------
class TestReviewExtensions:
    def test_rng_parity_with_a_suppressed_ion(self):
        # The load-bearing justification for the unconditional draw (Q5): a *suppressed*
        # ion in the components form still consumes its draw, so the components RNG stream
        # matches the scalar loop ion-for-ion even when some ions never fire.
        n = np.array([10, 8])                            # ion 0 in-band, ion 1 suppressed
        E = np.array([0.5 * (_d0(10) + _sigma(10)), 2.0 * _sigma(8)])
        vx, vy, vz = np.array([3.0, 3.0]), np.array([-4.0, -4.0]), np.array([12.0, 12.0])
        m = np.array([M, M])
        rng_c = _ConstRNG(0.0)
        (n_p, *_rest, fired) = _step_components(rng_c, E, n, vx, vy, vz, m)
        assert fired[0] and not fired[1]                 # ion 1 suppressed -> no fire
        assert n_p[1] == 8                               # untouched
        assert rng_c.n_draws == 2                        # but still drew for the suppressed ion

    def test_nan_E_gives_zero_rate_no_crash(self):
        # Deliberate non-guard (matches pickup/drag): a non-finite E_int is an upstream
        # bug, but must not silently mis-fire or raise -- it yields k = 0 (no fire).
        assert rrk_rate(np.nan, 8, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0

    def test_inf_E_shell_suppressed_n1_direct(self):
        # E -> inf is far above Sigma(n): the shell (n>=2) is suppressed; the n=1 diatomic
        # still direct-dissociates (E > D_0(1)).
        assert rrk_rate(np.inf, 8, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0
        assert rrk_rate(np.inf, 1, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == pytest.approx(
            NU_EVAP_PER_PS
        )

    def test_n0_has_no_rung(self):
        # n = 0 (fully stripped) is neither the n=1 direct branch nor an n>=2 shell -> k=0.
        assert rrk_rate(0.05, 0, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA) == 0.0
        k = rrk_rate(np.array([0.05, 0.05]), np.array([0, 10]),
                     nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)
        assert k[0] == 0.0 and k[1] > 0.0

    def test_negative_n_fails_loud(self):
        # Sigma(n) via ladder_cumsum rejects n < 0 (fail-loud, not a silent 0).
        with pytest.raises(ValueError):
            rrk_rate(0.05, -1, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA)

    def test_dof_override_s1_is_flat_nu_in_band(self):
        # s = 1 (the guard's lower edge): the bracket exponent is 0, so k = nu everywhere
        # in the self-bound band -- but still 0 below threshold and while suppressed.
        n = 10
        E_band = 0.5 * (_d0(n) + _sigma(n))
        assert rrk_rate(E_band, n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                        evap_rrk_dof=1.0) == pytest.approx(NU_EVAP_PER_PS)
        assert rrk_rate(0.5 * _d0(n), n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                        evap_rrk_dof=1.0) == 0.0                        # below threshold
        assert rrk_rate(2.0 * _sigma(n), n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                        evap_rrk_dof=1.0) == 0.0                        # suppressed

    def test_dof_override_respects_gate(self):
        # An override changes the in-band rate but never overrides the gate itself.
        n = 8
        assert rrk_rate(2.0 * _sigma(n), n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                        evap_rrk_dof=4.0) == 0.0                        # still suppressed
        assert rrk_rate(0.5 * _d0(n), n, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                        evap_rrk_dof=4.0) == 0.0                        # still below threshold

    def test_gate_diagnostics_vectorized_with_override(self):
        n = np.array([7, 7])
        E = np.array([0.5, 0.5])
        # Under Sigma(7) both are self-unbound (E large); a big override makes them bound.
        assert not np.any(is_self_bound(E, n, picture=PIC, kappa=KAPPA))
        sb = is_self_bound(E, n, picture=PIC, kappa=KAPPA, gate_onset_eV=10.0)
        assert np.all(sb)
        g = gate_margin_eV(E, n, picture=PIC, kappa=KAPPA, gate_onset_eV=10.0)
        np.testing.assert_allclose(g, E - 10.0)

    def test_fire_books_two_distinct_channels(self):
        # K1 drain (eV, into E_int) and the cold-shed KE defect (mechanical, into
        # E_mass_transfer) are distinct injections -- both present on a fire, different
        # quantities, no double count.
        n = 10
        E = 0.5 * (_d0(n) + _sigma(n))
        res = _step(_ConstRNG(0.0), E, n)
        assert res.dE_int_eV < 0.0
        assert res.dE_mass_transfer < 0.0
        assert res.dE_int_eV != res.dE_mass_transfer      # eV vs amu*A^2/ps^2 -- not the same

    def test_negative_nu_rejected(self):
        # Post-review fix (2026-07-02, Phase-B review): nu < 0 gives k < 0 ->
        # P_shed < 0, so evaporation silently never fires (draws live in [0, 1)) --
        # the same silent-shut-off class as the pickup lambda0 < 0 guard. Module-level
        # defense mirroring the config-load check (config.check_biphasic_config).
        with pytest.raises(ValueError, match="nu must be >= 0"):
            rrk_rate(0.05, 10, nu=-1.0, picture=PIC, kappa=KAPPA)

    def test_negative_nu_rejected_through_step(self):
        with pytest.raises(ValueError, match="nu must be >= 0"):
            evaporation_step(
                rng=_ConstRNG(0.0), E_int_eV=0.05, n=10, v=V, m_amu=M,
                nu=-1.0, picture=PIC, kappa=KAPPA, dt_ps=0.01,
            )

    def test_negative_n_fails_loud_with_gate_override(self):
        # Post-review fix (2026-07-02): the gate_onset_eV override used to bypass the
        # ladder_cumsum n-validation entirely (silent k = 0 for n = -1); both threshold
        # sources now validate n identically.
        with pytest.raises(ValueError, match="n >= 0"):
            rrk_rate(0.05, -1, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                     gate_onset_eV=1.0)

    def test_fractional_n_fails_loud_with_gate_override(self):
        # Message wording locked public-facing (Slice-G re-review fix 2026-07-02):
        # the context is "gate threshold", not the private helper name
        # _gate_threshold_eV a caller of rrk_rate has never seen.
        with pytest.raises(ValueError, match="gate threshold requires integer occupancy"):
            rrk_rate(0.05, 2.5, nu=NU_EVAP_PER_PS, picture=PIC, kappa=KAPPA,
                     gate_onset_eV=1.0)

    def test_gate_diagnostics_validate_n_with_override(self):
        # The bypass also covered is_self_bound / gate_margin_eV (they share
        # _gate_threshold_eV); locked here for all three entry points.
        with pytest.raises(ValueError, match="n >= 0"):
            is_self_bound(0.05, -1, picture=PIC, kappa=KAPPA, gate_onset_eV=1.0)
        with pytest.raises(ValueError, match="gate threshold requires integer occupancy"):
            gate_margin_eV(0.05, 2.5, picture=PIC, kappa=KAPPA, gate_onset_eV=1.0)

    def test_cold_shed_components_scalar_and_uniform_array_agree(self):
        # Lock the mass_jump generalization: the byte-identical scalar-mass path and the
        # new per-ion path agree when the per-ion masses are all equal to the scalar.
        from i2_helium_md.physics.mass_jump import cold_shed_velocity_components
        vx = np.array([3.0, 1.0, -2.0])
        vy = np.array([-4.0, 0.0, 5.0])
        vz = np.array([12.0, -1.0, 0.0])
        m_scalar = 200.0
        a = cold_shed_velocity_components(vx, vy, vz, m_scalar)
        b = cold_shed_velocity_components(vx, vy, vz, np.full(vx.shape, m_scalar))
        for j in range(3):                                # vx, vy, vz
            np.testing.assert_allclose(a[j], b[j])
        np.testing.assert_allclose(a[3], b[3])            # m_plus scalar == array value
        np.testing.assert_allclose(a[4], b[4])            # per-atom defect
