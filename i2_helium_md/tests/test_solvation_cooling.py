r"""Tests for i2_helium_md/physics/solvation_cooling.py (Tier-2 Phase-A Slice K).

Slice K is pure, closed-form energetics (no integrator state, no RNG, no mass and
no velocity), so every test is an oracle / identity comparison against the golden
values in ``docs/drag_port/Tier2/TIER2_PHASE_A_IMPLEMENTATION_PLAN.md`` (§2.2 K, §4).
It consumes Slice L's ``ladder_cumsum`` for ``Sigma(n)``; one independence test
stubs that call to prove K depends on L only through the cumulative gate.

Encoded form (plan §2.2 K)
--------------------------
    |S(N)| = |S| * Sigma(N) / Sigma(n*)          (collective solvation, occupancy-resolved)
    E_inf(N) = -|S(N)|                            (cooling asymptote)
    E_solv.struct(N) = -Sigma(N)  +  -(|S(N)|-Sigma(N))  +  E_int(N)
                       \__pair__/    \__E_elec <= 0 __/      (E_int owned by Slice U)
    dE_solv.struct/dt|K2 = -(E_solv.struct - E_inf)/tau     [eV/ps]

Sourced anchor: |S| = 0.308 eV (DFT first-shell solvation, MASS K2) = 2484 cm^-1
total, 118 cm^-1/atom at n* = 21. Because |S| = 0.308 > Sigma(n*) (mix 0.17-0.19 /
X2 0.25-0.28 eV), the electrostriction marginal is non-positive and |S(n*)| = |S|
exactly for any picture/kappa (the ratio Sigma(n*)/Sigma(n*) = 1).
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import (
    EV_PER_WAVENUMBER,
    N_STAR,
    S_ABS_EV,
)
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum
from i2_helium_md.physics.solvation_cooling import (
    e_bind_pair_eV,
    e_electrostriction_eV,
    e_infinity_eV,
    newton_cool_step,
    s_collective_eV,
)

# kappa sweep mirrors the Slice-L band/decoupling range (plan kappa in [0.3, 5]).
KAPPAS = (0.3, 0.5, 1.0, 2.0, 5.0)
PICTURES = ("statistical_mixture", "x2_only")


# ---------------------------------------------------------------------------
# |S(N)| collective solvation: occupancy-resolved, pinned to |S| at n*.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_s_collective_at_nstar_equals_S_abs(picture, kappa):
    # |S(n*)| = |S| * Sigma(n*)/Sigma(n*) = |S| exactly, for any picture/kappa.
    assert s_collective_eV(N_STAR, picture=picture, kappa=kappa) == pytest.approx(
        S_ABS_EV, abs=1e-9
    )


def test_s_collective_per_atom_marginal_is_118_wavenumber():
    # |S|/n* = 0.308/21 eV = 118 cm^-1 (170 K) at n* = 21 (3 figs; plan §4).
    per_atom_eV = s_collective_eV(N_STAR, kappa=1.0) / N_STAR
    per_atom_wavenumber = per_atom_eV / EV_PER_WAVENUMBER
    assert per_atom_wavenumber == pytest.approx(118.3, abs=0.5)


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_s_collective_zero_at_empty_shell(picture, kappa):
    assert s_collective_eV(0, picture=picture, kappa=kappa) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("kappa", KAPPAS)
def test_s_collective_monotone_increasing_in_N(kappa):
    n = np.arange(0, N_STAR + 1)
    s = s_collective_eV(n, kappa=kappa)
    assert np.all(np.diff(s) >= 0.0)


# ---------------------------------------------------------------------------
# E_inf(N) = -|S(N)|: cooling asymptote. Monotone decreasing, -> 0 as N -> 0.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("kappa", KAPPAS)
def test_e_infinity_is_negative_s_collective(kappa):
    n = np.arange(0, N_STAR + 1)
    assert np.allclose(
        e_infinity_eV(n, kappa=kappa), -s_collective_eV(n, kappa=kappa)
    )


def test_e_infinity_reaches_zero_at_empty_shell():
    # OQ6: occupancy-resolved -> the total strip stays dynamically reachable; a
    # *fixed* full-shell asymptote would mechanically halt shedding.
    assert e_infinity_eV(0, kappa=1.0) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("kappa", KAPPAS)
def test_e_infinity_monotone_decreasing_in_N(kappa):
    n = np.arange(0, N_STAR + 1)
    assert np.all(np.diff(e_infinity_eV(n, kappa=kappa)) <= 0.0)


# ---------------------------------------------------------------------------
# Binding split: E_inf = E_bind_pair + E_elec, with E_elec <= 0 everywhere.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_e_bind_pair_is_negative_sigma(picture, kappa):
    n = np.arange(0, N_STAR + 1)
    assert np.allclose(
        e_bind_pair_eV(n, picture=picture, kappa=kappa),
        -ladder_cumsum(n, picture=picture, kappa=kappa),
    )


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_electrostriction_nonpositive_everywhere(picture, kappa):
    n = np.arange(0, N_STAR + 1)
    assert np.all(e_electrostriction_eV(n, picture=picture, kappa=kappa) <= 0.0)


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_electrostriction_equals_split_definition(picture, kappa):
    n = np.arange(1, N_STAR + 1)
    expected = -(
        s_collective_eV(n, picture=picture, kappa=kappa)
        - ladder_cumsum(n, picture=picture, kappa=kappa)
    )
    assert np.allclose(
        e_electrostriction_eV(n, picture=picture, kappa=kappa), expected
    )


def test_electrostriction_zero_at_empty_shell():
    assert e_electrostriction_eV(0, kappa=1.0) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_split_closes_to_e_infinity(picture, kappa):
    # The central split identity: E_inf(N) = E_bind_pair(N) + E_elec(N).
    n = np.arange(0, N_STAR + 1)
    lhs = e_infinity_eV(n, picture=picture, kappa=kappa)
    rhs = e_bind_pair_eV(n, picture=picture, kappa=kappa) + e_electrostriction_eV(
        n, picture=picture, kappa=kappa
    )
    assert np.allclose(lhs, rhs)


# ---------------------------------------------------------------------------
# newton_cool_step: exact exponential relaxation, fixed point, dt-robustness.
# ---------------------------------------------------------------------------
def test_newton_cool_step_exact_decay_factor():
    N, tau, dt = 19, 6.55, 1.0
    e_inf = e_infinity_eV(N, kappa=1.0)
    e0 = 0.0   # a hot E_solv.struct above the asymptote
    e1 = newton_cool_step(e0, N, tau_ps=tau, dt_ps=dt, kappa=1.0)
    assert (e1 - e_inf) / (e0 - e_inf) == pytest.approx(np.exp(-dt / tau))


def test_newton_cool_step_fixed_point_at_asymptote():
    N, tau, dt = 19, 6.55, 3.0
    e_inf = e_infinity_eV(N, kappa=1.0)
    assert newton_cool_step(e_inf, N, tau_ps=tau, dt_ps=dt, kappa=1.0) == pytest.approx(
        e_inf, abs=1e-12
    )


def test_newton_cool_step_dt_robust_one_big_vs_many_small():
    N, tau, total = 19, 6.55, 10.0
    e0 = 0.05
    one_step = newton_cool_step(e0, N, tau_ps=tau, dt_ps=total, kappa=1.0)
    e = e0
    for _ in range(100):
        e = newton_cool_step(e, N, tau_ps=tau, dt_ps=total / 100.0, kappa=1.0)
    assert e == pytest.approx(one_step, abs=1e-12)


def test_newton_cool_step_relaxes_toward_asymptote():
    N, tau, dt = 19, 6.55, 1.0
    e_inf = e_infinity_eV(N, kappa=1.0)
    e0 = 0.05
    e1 = newton_cool_step(e0, N, tau_ps=tau, dt_ps=dt, kappa=1.0)
    # one step moves strictly toward (but not past) the asymptote
    assert e_inf < e1 < e0


def test_newton_cool_step_rejects_nonpositive_tau():
    with pytest.raises(ValueError, match="tau"):
        newton_cool_step(0.0, 19, tau_ps=0.0, dt_ps=1.0, kappa=1.0)


def test_newton_cool_step_rejects_negative_dt():
    # A negative dt gives exp(+|dt|/tau) > 1 -> silent anti-cooling (heats away
    # from E_inf). A backward cooling step is unphysical; fail loud (principle 4).
    with pytest.raises(ValueError, match="dt"):
        newton_cool_step(0.0, 19, tau_ps=6.55, dt_ps=-1.0, kappa=1.0)


def test_newton_cool_step_zero_dt_is_noop():
    e0 = 0.05
    assert newton_cool_step(e0, 19, tau_ps=6.55, dt_ps=0.0, kappa=1.0) == pytest.approx(e0)


def test_newton_cool_step_converges_to_asymptote_over_long_time():
    N, tau = 19, 6.55
    e_inf = e_infinity_eV(N, kappa=1.0)
    e = 0.05
    for _ in range(2000):
        e = newton_cool_step(e, N, tau_ps=tau, dt_ps=0.5, kappa=1.0)
    assert e == pytest.approx(e_inf, abs=1e-9)


def test_newton_cool_step_targets_picture_correct_asymptote():
    # The step relaxes toward the *configured-picture* asymptote, not a fixed one:
    # X2 binds deeper than the mixture, so its asymptote is more negative.
    N, tau, dt = 19, 6.55, 1000.0   # ~fully relaxed
    e_mix = newton_cool_step(0.0, N, tau_ps=tau, dt_ps=dt, picture="statistical_mixture", kappa=1.0)
    e_x2 = newton_cool_step(0.0, N, tau_ps=tau, dt_ps=dt, picture="x2_only", kappa=1.0)
    assert e_x2 < e_mix < 0.0
    assert e_mix == pytest.approx(e_infinity_eV(N, picture="statistical_mixture", kappa=1.0), abs=1e-6)
    assert e_x2 == pytest.approx(e_infinity_eV(N, picture="x2_only", kappa=1.0), abs=1e-6)


def test_newton_cool_step_array_N_scalar_E_returns_ndarray():
    N = np.array([10, 15, 21])
    out = newton_cool_step(0.0, N, tau_ps=6.55, dt_ps=1.0, kappa=1.0)
    assert isinstance(out, np.ndarray)
    assert out.shape == N.shape


def test_newton_cool_step_scalar_returns_float():
    out = newton_cool_step(0.0, 19, tau_ps=6.55, dt_ps=1.0, kappa=1.0)
    assert isinstance(out, float)


# ---------------------------------------------------------------------------
# newton_cool_step: rho_ratio spatial-gate factor (the cooling_spatial_gate arm).
# The `density_scaled` K2 arm passes rho_ratio in [0,1] from the shared erf gate;
# tau_eff = tau/rho_ratio, i.e. decay = exp(-dt*rho_ratio/tau). rho_ratio=1 recovers
# the ungated closed form byte-for-byte (the locked `none` arm); rho_ratio=0 disables
# cooling (no droplet bath to radiate into outside the bubble).
# ---------------------------------------------------------------------------
def test_newton_cool_step_rho_ratio_one_equals_ungated():
    # dt*1.0 is exact float, so passing rho_ratio=1.0 must be byte-identical to
    # omitting it -- the guarantee the default `none` arm relies on.
    N, tau, dt = 19, 6.55, 1.0
    e0 = 0.05
    ungated = newton_cool_step(e0, N, tau_ps=tau, dt_ps=dt, kappa=1.0)
    gated1 = newton_cool_step(e0, N, tau_ps=tau, dt_ps=dt, kappa=1.0, rho_ratio=1.0)
    assert gated1 == ungated


def test_newton_cool_step_rho_ratio_zero_is_noop():
    # rho_ratio=0 -> decay=exp(0)=1 -> E unchanged: cooling off outside the droplet.
    e0 = 0.05
    out = newton_cool_step(e0, 19, tau_ps=6.55, dt_ps=1.0, kappa=1.0, rho_ratio=0.0)
    assert out == pytest.approx(e0, abs=1e-12)


def test_newton_cool_step_rho_ratio_scales_decay():
    # The gate enters as tau_eff = tau/rho_ratio, i.e. decay = exp(-dt*rho_ratio/tau).
    N, tau, dt, rho = 19, 6.55, 1.0, 0.3
    e_inf = e_infinity_eV(N, kappa=1.0)
    e0 = 0.0
    e1 = newton_cool_step(e0, N, tau_ps=tau, dt_ps=dt, kappa=1.0, rho_ratio=rho)
    assert (e1 - e_inf) / (e0 - e_inf) == pytest.approx(np.exp(-dt * rho / tau))


def test_newton_cool_step_rejects_negative_rho_ratio():
    # A negative factor gives exp(+|.|) > 1 -> anti-cooling (heats away from E_inf),
    # the same failure class as a negative dt; fail loud (principle 4).
    with pytest.raises(ValueError, match="rho_ratio"):
        newton_cool_step(0.0, 19, tau_ps=6.55, dt_ps=1.0, kappa=1.0, rho_ratio=-0.1)


def test_newton_cool_step_array_rho_ratio_returns_ndarray():
    # A per-atom (2N,) rho_ratio broadcasts against the E array and returns an ndarray
    # (the driver passes one gate value per ion row) even with scalar E/N.
    rho = np.array([0.0, 0.5, 1.0])
    out = newton_cool_step(0.0, 19, tau_ps=6.55, dt_ps=1.0, kappa=1.0, rho_ratio=rho)
    assert isinstance(out, np.ndarray)
    assert out.shape == rho.shape


# ---------------------------------------------------------------------------
# Cold-shed neutrality on the pair + E_int sub-sum (K1 cap K2 non-overlap),
# modulo the A8 marginal-electrostriction bath booking (the collective
# electrostriction term |S| > Sigma(n*) is the bath booking, excluded here).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
@pytest.mark.parametrize("N", (5, 14, 21))
def test_cold_shed_neutral_on_pair_plus_e_int_subsum(picture, kappa, N):
    # A shed N -> N-1 drains the outermost rung D_0(N) from E_int, while the pair
    # binding gains +D_0(N) (Sigma(N) - Sigma(N-1) = D_0(N)). The discrete-trackable
    # sub-sum (pair + E_int) is therefore unchanged to machine precision.
    e_int = 0.123                                  # arbitrary reservoir occupancy
    d0_N = d0_of_n(N, picture=picture, kappa=kappa)
    before = e_bind_pair_eV(N, picture=picture, kappa=kappa) + e_int
    after = e_bind_pair_eV(N - 1, picture=picture, kappa=kappa) + (e_int - d0_N)
    assert after - before == pytest.approx(0.0, abs=1e-15)


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_e_int_eq_is_zero_tautology(picture, kappa):
    # SPLIT-CONSISTENCY (tautology), not an independent anchor: with E_int defined
    # as the residual E_solv.struct - E_bind_pair - E_elec and E_inf = E_bind_pair +
    # E_elec, the residual is identically 0 at equilibrium (E_solv.struct = E_inf)
    # by construction (MASS line 881).
    N = 19
    e_solv_struct_eq = e_infinity_eV(N, picture=picture, kappa=kappa)
    residual = (
        e_solv_struct_eq
        - e_bind_pair_eV(N, picture=picture, kappa=kappa)
        - e_electrostriction_eV(N, picture=picture, kappa=kappa)
    )
    assert residual == pytest.approx(0.0, abs=1e-15)


# ---------------------------------------------------------------------------
# Independence: K depends on L only through ladder_cumsum (Sigma).
# ---------------------------------------------------------------------------
def test_s_collective_depends_on_ladder_only_through_cumsum(monkeypatch):
    import i2_helium_md.physics.solvation_cooling as sc

    # Stub Sigma(n) = n: then |S(N)| = |S| * N / n* exactly, independent of the
    # real ladder shape -- proving K composes L only via ladder_cumsum.
    monkeypatch.setattr(
        sc, "ladder_cumsum", lambda n, *, picture, kappa, ladder=None: float(n)
    )
    s = sc.s_collective_eV(10, picture="statistical_mixture", kappa=1.0)
    assert s == pytest.approx(S_ABS_EV * 10 / N_STAR)


# ---------------------------------------------------------------------------
# Vectorization: scalar-in -> float, array-in -> ndarray (mirrors Slice L).
# ---------------------------------------------------------------------------
def test_scalar_in_returns_float():
    for fn in (s_collective_eV, e_infinity_eV, e_bind_pair_eV, e_electrostriction_eV):
        out = fn(14, kappa=1.0)
        assert isinstance(out, float)


def test_array_in_returns_ndarray_same_shape():
    n = np.arange(0, N_STAR + 1)
    for fn in (s_collective_eV, e_infinity_eV, e_bind_pair_eV, e_electrostriction_eV):
        out = fn(n, kappa=1.0)
        assert isinstance(out, np.ndarray)
        assert out.shape == n.shape


def test_s_collective_custom_s_abs_overrides_default():
    # The s_abs_eV knob is honoured (config field surfaces this override).
    assert s_collective_eV(N_STAR, kappa=1.0, s_abs_eV=0.5) == pytest.approx(0.5, abs=1e-9)


def test_e_infinity_custom_s_abs_overrides_default():
    assert e_infinity_eV(N_STAR, kappa=1.0, s_abs_eV=0.5) == pytest.approx(-0.5, abs=1e-9)


# ---------------------------------------------------------------------------
# Structural robustness: |S| > Sigma(n*) for ALL kappa > 0 (not just the
# sourced [0.3,5] band), because Sigma(n*) <= n* * D_0(1) caps below |S|. This
# is the precondition that makes E_elec <= 0 a structural guarantee, not a
# coincidence of the sourced kappa range.
# ---------------------------------------------------------------------------
WIDE_KAPPAS = (0.05, 0.1, 1.0, 10.0, 100.0)


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", WIDE_KAPPAS)
def test_sigma_at_nstar_below_S_abs_for_all_kappa(picture, kappa):
    # |S| = 0.308 eV exceeds the deepest possible full-shell sum 21*D_0(1)
    # (X2 0.278 / mix 0.194), so |S| > Sigma(n*) at every kappa.
    assert ladder_cumsum(N_STAR, picture=picture, kappa=kappa) < S_ABS_EV


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", WIDE_KAPPAS)
def test_electrostriction_nonpositive_for_all_kappa(picture, kappa):
    n = np.arange(0, N_STAR + 1)
    assert np.all(e_electrostriction_eV(n, picture=picture, kappa=kappa) <= 0.0)


# ---------------------------------------------------------------------------
# Fail-loud propagation: K inherits Slice L's ladder_cumsum guards (negative /
# genuinely-fractional occupancy) through every consumer -- no silent garbage.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "fn", [s_collective_eV, e_infinity_eV, e_bind_pair_eV, e_electrostriction_eV]
)
@pytest.mark.parametrize("bad_N", [-1, 2.5])
def test_consumers_propagate_ladder_failloud(fn, bad_N):
    with pytest.raises(ValueError):
        fn(bad_N, kappa=1.0)


def test_newton_cool_step_propagates_ladder_failloud():
    with pytest.raises(ValueError):
        newton_cool_step(0.0, -1, tau_ps=6.55, dt_ps=1.0, kappa=1.0)


# ---------------------------------------------------------------------------
# Slice T2 (§I.10): ladder= injection -- tabulated == form_u on the fed rungs.
# ---------------------------------------------------------------------------
class TestLadderInjectionSliceT2:
    KAPPA = 1.0

    def _form_u_ladder(self, length=32):
        from i2_helium_md.physics.dissociation_ladder import tabulated_ladder

        rungs = np.atleast_1d(
            d0_of_n(np.arange(1, length + 1), kappa=self.KAPPA)
        )
        return tabulated_ladder(rungs)

    def test_binding_split_form_u_fed_table_bit_identical(self):
        lad = self._form_u_ladder()
        N = np.arange(0, 22)
        for fn in (s_collective_eV, e_infinity_eV, e_bind_pair_eV,
                   e_electrostriction_eV):
            np.testing.assert_array_equal(
                np.asarray(fn(N, kappa=self.KAPPA, ladder=lad)),
                np.asarray(fn(N, kappa=self.KAPPA)),
            )

    def test_newton_cool_step_form_u_fed_table_bit_identical(self):
        lad = self._form_u_ladder()
        out_u = newton_cool_step(0.1, 21, tau_ps=6.55, dt_ps=0.01, kappa=self.KAPPA)
        out_t = newton_cool_step(
            0.1, 21, tau_ps=6.55, dt_ps=0.01, kappa=self.KAPPA, ladder=lad
        )
        assert out_t == out_u

    def test_distinct_table_moves_binding_pair(self):
        # Liveness: e_bind_pair = -Sigma(N) tracks the injected table exactly.
        from i2_helium_md.physics.dissociation_ladder import tabulated_ladder

        lad = tabulated_ladder(tuple(0.01 for _ in range(21)))
        assert e_bind_pair_eV(3, kappa=self.KAPPA, ladder=lad) == pytest.approx(
            -0.03, abs=1e-15
        )
