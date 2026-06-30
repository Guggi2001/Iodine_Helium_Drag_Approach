r"""Tests for i2_helium_md/physics/internal_energy_budget.py (Tier-2 Phase-A Slice U).

Slice U is the pure, stateless per-channel ``E_int`` bookkeeping layer (no reservoir
state, no RNG, no integrator, no mass and no velocity). Every test is an oracle /
identity comparison against the MASS §6 budget rules and the golden values in
``docs/drag_port/Tier2/TIER2_PHASE_A_IMPLEMENTATION_PLAN.md`` (§2.2 U, §3, §4). It
consumes Slice L's rungs (``d0_of_n``/``ladder_cumsum``) and Slice K's binding split
(``e_bind_pair_eV``/``e_electrostriction_eV``); independence tests stub those calls to
prove U composes its neighbours only through the documented surface.

Encoded form (plan §2.2 U; MASS §6)
-----------------------------------
    S2 onset:        E_int(0)      = f_int * E_avail^ion
    S1 pickup:       dE_int        = +f_ret * D_0(n+1)   (bath gets (1-f_ret)*D_0(n+1))
    K1 shed:         dE_int        = -D_0(n)
    reconstruction:  E_int         = E_solv.struct - E_bind^pair(N) - E_elec(N)   (post-t* only)
    f_int floor:     f_int^floor   = Sigma(n*) / E_avail^ion

Scenario-keyed budgets (MASS §6 S2, CALIBRATION row 15): 0.80 eV validation (d=9 A),
2.70 eV production (R_e). f_int floor bands (plan §4, MASS row 14): @0.80 eV X2
0.31-0.35 / mix 0.21-0.24; @2.70 eV X2 0.09-0.10 / mix ~0.065 (advisory, not enforced).
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics import internal_energy_budget as ieb
from i2_helium_md.physics.constants import N_STAR
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum
from i2_helium_md.physics.internal_energy_budget import (
    dE_int_pickup_eV,
    dE_int_shed_eV,
    e_int_onset_eV,
    f_int_floor,
    pickup_bath_release_eV,
    reconstruct_e_int_eV,
)
from i2_helium_md.physics.solvation_cooling import (
    e_bind_pair_eV,
    e_electrostriction_eV,
    e_infinity_eV,
)

KAPPAS = (0.3, 0.5, 1.0, 2.0, 5.0)
PICTURES = ("statistical_mixture", "x2_only")
E_AVAIL_VALIDATION = 0.80  # eV, d = 9 A
E_AVAIL_PRODUCTION = 2.70  # eV, R_e

# Scenario-keyed f_int floor bands (plan §4 / MASS row 14): (lo, hi) per picture, on
# round(floor, 2). Mirrors the Slice-L band-membership convention -- the rounded
# kappa-endpoints ARE the MASS-stated band edges (e.g. kappa=5 mix Sigma(21)=0.19306
# -> /0.80 = 0.2413 -> 0.24). The 2.70 eV mix band is MASS's "~0.065" at 2 decimals.
FLOOR_BANDS = {
    (E_AVAIL_VALIDATION, "statistical_mixture"): (0.21, 0.24),
    (E_AVAIL_VALIDATION, "x2_only"): (0.31, 0.35),
    (E_AVAIL_PRODUCTION, "statistical_mixture"): (0.06, 0.07),
    (E_AVAIL_PRODUCTION, "x2_only"): (0.09, 0.10),
}


# ---------------------------------------------------------------------------
# S2 onset deposit: E_int(0) = f_int * E_avail.
# ---------------------------------------------------------------------------
def test_onset_is_f_int_times_e_avail():
    assert e_int_onset_eV(f_int=0.3, e_avail_eV=0.80) == pytest.approx(0.24, abs=1e-12)


def test_onset_zero_f_int_is_zero():
    assert e_int_onset_eV(f_int=0.0, e_avail_eV=2.70) == pytest.approx(0.0, abs=1e-15)


def test_onset_full_f_int_is_e_avail():
    assert e_int_onset_eV(f_int=1.0, e_avail_eV=2.70) == pytest.approx(2.70, abs=1e-12)


def test_onset_does_not_enforce_floor():
    # The self-unbound floor is advisory (MASS S2 / CALIBRATION row 14, "not a
    # constraint"): a sub-floor f_int is accepted, not rejected. e_int_onset_eV is
    # pure arithmetic and must not consult f_int_floor.
    floor = f_int_floor(e_avail_eV=0.80, picture="x2_only", kappa=1.0)
    sub_floor = 0.5 * floor
    assert e_int_onset_eV(f_int=sub_floor, e_avail_eV=0.80) == pytest.approx(
        sub_floor * 0.80, abs=1e-12
    )


# ---------------------------------------------------------------------------
# S1 pickup: +f_ret * D_0(n+1); the bath gets the (1-f_ret) remainder; split closes.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
@pytest.mark.parametrize("n", [1, 5, 14, 20])
def test_pickup_is_f_ret_times_next_rung(picture, kappa, n):
    f_ret = 0.4
    expected = f_ret * d0_of_n(n + 1, picture=picture, kappa=kappa)
    assert dE_int_pickup_eV(n, f_ret=f_ret, picture=picture, kappa=kappa) == pytest.approx(
        expected, rel=1e-12
    )


def test_pickup_heats_e_int_positive():
    # Pickup is a source: f_ret in (0, 1] gives a strictly positive increment.
    assert dE_int_pickup_eV(10, f_ret=0.5, kappa=1.0) > 0.0


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
@pytest.mark.parametrize("n", [1, 5, 14, 20])
def test_pickup_bath_release_is_remainder(picture, kappa, n):
    f_ret = 0.4
    expected = (1.0 - f_ret) * d0_of_n(n + 1, picture=picture, kappa=kappa)
    assert pickup_bath_release_eV(
        n, f_ret=f_ret, picture=picture, kappa=kappa
    ) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
@pytest.mark.parametrize("n", [1, 5, 14, 20])
def test_pickup_split_closes_to_full_rung(picture, kappa, n):
    # S1 split closes exactly: f_ret*D_0 + (1-f_ret)*D_0 = D_0 (machine precision).
    f_ret = 0.37
    retained = dE_int_pickup_eV(n, f_ret=f_ret, picture=picture, kappa=kappa)
    released = pickup_bath_release_eV(n, f_ret=f_ret, picture=picture, kappa=kappa)
    full = d0_of_n(n + 1, picture=picture, kappa=kappa)
    assert retained + released == pytest.approx(full, abs=1e-15)


# ---------------------------------------------------------------------------
# K1 shed: -D_0(n).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
@pytest.mark.parametrize("n", [1, 5, 14, 21])
def test_shed_drains_one_rung(picture, kappa, n):
    expected = -d0_of_n(n, picture=picture, kappa=kappa)
    assert dE_int_shed_eV(n, picture=picture, kappa=kappa) == pytest.approx(
        expected, rel=1e-12
    )


def test_shed_is_negative():
    assert dE_int_shed_eV(10, kappa=1.0) < 0.0


# ---------------------------------------------------------------------------
# Jump-consistency (MASS §6): at a shed the K1 E_int drain exactly offsets the
# pair-binding gain, so the (pair + E_int) sub-sum is neutral (cold-shed identity).
# This is the variable-swap check expressed on the Slice-K split.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
@pytest.mark.parametrize("n", [2, 7, 21])
def test_shed_pair_plus_eint_subsum_neutral(picture, kappa, n):
    # dE_int_shed(n) = -D_0(n); the pair binding E_bind^pair = -Sigma rises by +D_0(n)
    # going n -> n-1. Their sum is 0 to machine precision (MASS no-double-count).
    d_eint = dE_int_shed_eV(n, picture=picture, kappa=kappa)
    d_pair = e_bind_pair_eV(n - 1, picture=picture, kappa=kappa) - e_bind_pair_eV(
        n, picture=picture, kappa=kappa
    )
    assert d_eint + d_pair == pytest.approx(0.0, abs=1e-15)


# ---------------------------------------------------------------------------
# Post-t* reconstruction: E_int = E_solv.struct - E_bind^pair(N) - E_elec(N).
# Guard: raises pre-t*. Equilibrium micro-case recovers E_int^eq = 0.
# ---------------------------------------------------------------------------
def test_reconstruct_raises_pre_crossing():
    with pytest.raises(ValueError, match="post_crossing"):
        reconstruct_e_int_eV(-0.2, 21, kappa=1.0, post_crossing=False)


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
@pytest.mark.parametrize("N", [3, 10, 21])
def test_reconstruct_matches_split(picture, kappa, N):
    e_solv = -0.15  # arbitrary post-t* E_solv.struct
    expected = (
        e_solv
        - e_bind_pair_eV(N, picture=picture, kappa=kappa)
        - e_electrostriction_eV(N, picture=picture, kappa=kappa)
    )
    assert reconstruct_e_int_eV(
        e_solv, N, picture=picture, kappa=kappa, post_crossing=True
    ) == pytest.approx(expected, abs=1e-15)


@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
@pytest.mark.parametrize("N", [3, 10, 21])
def test_reconstruct_equilibrium_is_zero(picture, kappa, N):
    # At equilibrium E_solv.struct = E_inf(N) = E_bind^pair + E_elec, so the
    # reconstructed residual E_int^eq is identically 0 (MASS line 881; a split-
    # consistency tautology, asserted here as a build check, not an independent anchor).
    e_solv_eq = e_infinity_eV(N, picture=picture, kappa=kappa)
    assert reconstruct_e_int_eV(
        e_solv_eq, N, picture=picture, kappa=kappa, post_crossing=True
    ) == pytest.approx(0.0, abs=1e-15)


# ---------------------------------------------------------------------------
# f_int floor = Sigma(n*) / E_avail: scenario-keyed bands, n* = N_STAR.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_f_int_floor_is_gate_over_e_avail(picture, kappa):
    expected = ladder_cumsum(N_STAR, picture=picture, kappa=kappa) / 0.80
    assert f_int_floor(e_avail_eV=0.80, picture=picture, kappa=kappa) == pytest.approx(
        expected, rel=1e-12
    )


@pytest.mark.parametrize("e_avail", [E_AVAIL_VALIDATION, E_AVAIL_PRODUCTION])
@pytest.mark.parametrize("picture", PICTURES)
@pytest.mark.parametrize("kappa", KAPPAS)
def test_f_int_floor_in_scenario_band(e_avail, picture, kappa):
    lo, hi = FLOOR_BANDS[(e_avail, picture)]
    floor = f_int_floor(e_avail_eV=e_avail, picture=picture, kappa=kappa)
    assert lo <= round(floor, 2) <= hi


def test_f_int_floor_rejects_nonpositive_e_avail():
    with pytest.raises(ValueError):
        f_int_floor(e_avail_eV=0.0, picture="statistical_mixture", kappa=1.0)


# ---------------------------------------------------------------------------
# Vectorisation: scalar-in -> float, array-in -> ndarray (mirrors L/K).
# ---------------------------------------------------------------------------
def test_pickup_vectorised_shape_parity():
    ns = np.array([1, 5, 14, 20])
    out = dE_int_pickup_eV(ns, f_ret=0.4, kappa=1.0)
    assert isinstance(out, np.ndarray) and out.shape == ns.shape
    for i, n in enumerate(ns):
        assert out[i] == pytest.approx(
            dE_int_pickup_eV(int(n), f_ret=0.4, kappa=1.0), rel=1e-12
        )


def test_shed_vectorised_shape_parity():
    ns = np.array([1, 5, 14, 21])
    out = dE_int_shed_eV(ns, kappa=1.0)
    assert isinstance(out, np.ndarray) and out.shape == ns.shape


def test_scalar_returns_float():
    assert isinstance(dE_int_pickup_eV(5, f_ret=0.4, kappa=1.0), float)
    assert isinstance(dE_int_shed_eV(5, kappa=1.0), float)
    assert isinstance(e_int_onset_eV(f_int=0.3, e_avail_eV=0.80), float)
    assert isinstance(f_int_floor(e_avail_eV=0.80, kappa=1.0), float)


# ---------------------------------------------------------------------------
# Fail-loud propagation: the reconstruction rides K's binding split -> L's
# ``ladder_cumsum``, whose negative/fractional-N guard must surface through U
# (the d0-based pickup/shed wrappers inherit L's *smooth, unguarded* d0_of_n by
# design -- only the cumulative/indexed path guards).
# ---------------------------------------------------------------------------
def test_reconstruct_propagates_negative_N_guard():
    with pytest.raises(ValueError):
        reconstruct_e_int_eV(-0.1, -1, kappa=1.0, post_crossing=True)


def test_reconstruct_propagates_fractional_N_guard():
    with pytest.raises(ValueError):
        reconstruct_e_int_eV(-0.1, 2.5, kappa=1.0, post_crossing=True)


# ---------------------------------------------------------------------------
# Independence: U composes L (rungs) only through d0_of_n. Stub it and confirm.
# ---------------------------------------------------------------------------
def test_pickup_uses_ladder_only_through_d0_of_n(monkeypatch):
    monkeypatch.setattr(ieb, "d0_of_n", lambda n, *, picture, kappa: 1.0)
    # With every rung stubbed to 1.0 eV, pickup = f_ret * 1.0 and shed = -1.0.
    assert dE_int_pickup_eV(5, f_ret=0.4, kappa=1.0) == pytest.approx(0.4, abs=1e-15)
    assert dE_int_shed_eV(5, kappa=1.0) == pytest.approx(-1.0, abs=1e-15)


# ===========================================================================
# Review + test-hardening pass (2026-06-30): coverage/regression locks on
# already-correct behavior (independent numerical re-verification confirmed every
# §6 oracle to machine precision; no physics bug surfaced).
# ===========================================================================

# --- Vectorised onset (de-speculates the array branch; rule-2) --------------
def test_onset_vectorised_over_f_int():
    f = np.array([0.0, 0.1, 0.5, 1.0])
    out = e_int_onset_eV(f_int=f, e_avail_eV=0.80)
    assert isinstance(out, np.ndarray) and out.shape == f.shape
    assert out == pytest.approx(f * 0.80, abs=1e-15)


def test_onset_vectorised_over_e_avail_at_fixed_f_int():
    # Phase-F-style sweep: scalar f_int, array of budgets -> array of onsets.
    ea = np.array([0.80, 2.70])
    out = e_int_onset_eV(f_int=0.2, e_avail_eV=ea)
    assert isinstance(out, np.ndarray) and out.shape == ea.shape
    assert out == pytest.approx(0.2 * ea, abs=1e-15)


# --- Module is config-agnostic: out-of-[0,1] f_int accepted (non-guard lock) -
def test_onset_module_does_not_guard_f_int_range():
    # Deliberate non-guard (mirrors Slice K): the module is config-agnostic -- the
    # [0,1] bound is a config-load concern (check_internal_energy_budget_config), not
    # the pure helper's. A super-unity or negative f_int is computed, not rejected.
    assert e_int_onset_eV(f_int=1.5, e_avail_eV=0.80) == pytest.approx(1.2, abs=1e-12)
    assert e_int_onset_eV(f_int=-0.1, e_avail_eV=0.80) == pytest.approx(-0.08, abs=1e-12)


# --- Vectorised / float-discipline on reconstruction ------------------------
@pytest.mark.parametrize("picture", PICTURES)
def test_reconstruct_vectorised_over_N(picture):
    Ns = np.array([3, 10, 21])
    e_solv = -0.15
    out = reconstruct_e_int_eV(e_solv, Ns, picture=picture, kappa=1.0, post_crossing=True)
    assert isinstance(out, np.ndarray) and out.shape == Ns.shape
    for i, N in enumerate(Ns):
        assert out[i] == pytest.approx(
            reconstruct_e_int_eV(
                e_solv, int(N), picture=picture, kappa=1.0, post_crossing=True
            ),
            abs=1e-15,
        )


def test_reconstruct_scalar_returns_float():
    assert isinstance(
        reconstruct_e_int_eV(-0.1, 10, kappa=1.0, post_crossing=True), float
    )


def test_reconstruct_composes_K_split_only(monkeypatch):
    # Prove reconstruct = E_solv - e_bind_pair - e_electrostriction and touches K only
    # through those two terms (not e.g. e_infinity directly). Stub both to constants.
    monkeypatch.setattr(ieb, "e_bind_pair_eV", lambda N, *, picture, kappa: -0.2)
    monkeypatch.setattr(
        ieb, "e_electrostriction_eV", lambda N, *, picture, kappa, s_abs_eV: -0.05
    )
    # E_int = E_solv - (-0.2) - (-0.05) = E_solv + 0.25.
    assert reconstruct_e_int_eV(-0.1, 10, kappa=1.0, post_crossing=True) == pytest.approx(
        0.15, abs=1e-15
    )


# --- f_int_floor negative-budget guard (only zero was covered) --------------
def test_f_int_floor_rejects_negative_e_avail():
    with pytest.raises(ValueError, match="e_avail_eV"):
        f_int_floor(e_avail_eV=-0.80, picture="x2_only", kappa=1.0)


# --- Return-dtype discipline on the vectorised helpers -----------------------
def test_bath_release_vectorised_shape_parity():
    ns = np.array([1, 5, 14, 20])
    out = pickup_bath_release_eV(ns, f_ret=0.4, kappa=1.0)
    assert isinstance(out, np.ndarray) and out.shape == ns.shape


def test_pickup_shed_array_return_is_ndarray():
    ns = np.array([2, 7, 19])
    assert isinstance(dE_int_pickup_eV(ns, f_ret=0.4, kappa=1.0), np.ndarray)
    assert isinstance(dE_int_shed_eV(ns, kappa=1.0), np.ndarray)


# --- Picture symmetry: cooling_relaxed orders between mix and x2 (carried L) --
@pytest.mark.parametrize("n", [1, 10, 20])
def test_pickup_cooling_relaxed_orders_between_pictures(n):
    # The S1 increment inherits L's picture ordering mix < cooling_relaxed < x2 for
    # the picture-set rung scale (a provisional-blend ordering lock, not a value oracle).
    f_ret, kappa = 0.5, 1.0
    mix = dE_int_pickup_eV(n, f_ret=f_ret, picture="statistical_mixture", kappa=kappa)
    cr = dE_int_pickup_eV(n, f_ret=f_ret, picture="cooling_relaxed", kappa=kappa)
    x2 = dE_int_pickup_eV(n, f_ret=f_ret, picture="x2_only", kappa=kappa)
    assert mix < cr < x2
