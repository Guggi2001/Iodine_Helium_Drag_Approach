"""Tests for i2_helium_md/physics/dissociation_ladder.py (Tier-2 Phase-A Slice L).

Slice L is pure, closed-form energetics, so every test is an oracle comparison
against the golden values in ``docs/drag_port/Tier2/TIER2_PHASE_A_IMPLEMENTATION_PLAN.md``
(§2/§4) -- no reference-data files. The Form-U ladder turns the sourced first rung
``D_0(1)`` (picture-keyed), the bulk floor ``D_floor``, and the steepness ``kappa``
into the per-rung dissociation cost ``D_0(n)`` and its cumulative self-bound gate
threshold ``Sigma(n) = sum_{i<=n} D_0(i)``.

Oracle conversions (EV_PER_WAVENUMBER = 1/8065.543937 eV/cm^-1)
--------------------------------------------------------------
* x2_only first rung           106.9 cm^-1 = 0.013254 eV
* statistical_mixture          74.4  cm^-1 = 0.009224 eV
* cooling_relaxed (provisional) 90.65 cm^-1 = 0.011240 eV  (average of the two;
  a rule-2 declared-but-unread stub -- ordering only, value pinned at Phase F)
* floor                        4.97  cm^-1 = 6.162e-4 eV
* cliff centre                 n* + 1/2 = 21.5

The integrated gate threshold ``Sigma(21)`` bands (plan §2/§4, stated to two
decimals) are the **decoupling** oracle: mixture 0.17-0.19 eV / X2 0.25-0.28 eV,
varying only ~11% over kappa in [0.3, 5] -- i.e. the gate threshold is
~kappa-independent and pinnable ahead of the kappa fit. The bands are asserted at
the plan's stated two-decimal precision (``round(value, 2)`` in the band), since
the endpoints the plan rounds *are* the computed values: the kappa=5 mixture value
is 0.19306 -> 0.19, the kappa=0.3 X2 value is 0.24973 -> 0.25.
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import (
    D0_1_COOLING_RELAXED_WAVENUMBER,
    D0_1_MIX_WAVENUMBER,
    D0_1_X2_WAVENUMBER,
    D_FLOOR_WAVENUMBER,
    EV_PER_WAVENUMBER,
    N_STAR,
)
from i2_helium_md.physics.dissociation_ladder import (
    d0_of_n,
    first_rung_d0_eV,
    gate_threshold,
    ladder_cumsum,
    sigma,
    tabulated_ladder,
)

# --- Sourced oracles (cm^-1 -> eV) -----------------------------------------
D0_1_X2_EV = 106.9 * EV_PER_WAVENUMBER          # 0.013254 eV
D0_1_MIX_EV = 74.4 * EV_PER_WAVENUMBER          # 0.009224 eV
D0_1_COOLING_RELAXED_EV = 90.65 * EV_PER_WAVENUMBER  # 0.011240 eV (provisional avg)
D_FLOOR_EV = 4.97 * EV_PER_WAVENUMBER           # 6.162e-4 eV

# Integrated gate-threshold bands (plan §2/§4), stated to two decimals. Asserted
# via round(value, 2) membership -- the plan's stated endpoints are the rounded
# computed values (mixture 0.19306 -> 0.19; X2 0.24973 -> 0.25).
MIX_SIGMA21_BAND_EV = (0.17, 0.19)
X2_SIGMA21_BAND_EV = (0.25, 0.28)

# E_bind drag cross-check (plan §4): 0.1168 eV = 942 cm^-1 is an *upper-bound*
# sanity, never a rung-sum target. A representative first-shell sum is the same
# order of magnitude (within the 835-1571 cm^-1 representative band).
E_BIND_DRAG_WAVENUMBER = 942.0

# kappa sweep used for the band / decoupling asserts (plan kappa in [0.3, 5]).
KAPPAS = (0.3, 0.5, 1.0, 2.0, 5.0)


# ---------------------------------------------------------------------------
# First rung: picture-keyed, kappa-independent (the Form-U normalisation pins
# D_0(1) == sourced first rung exactly for every kappa).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("kappa", KAPPAS)
def test_first_rung_x2_matches_oracle_4figs(kappa):
    assert d0_of_n(1, picture="x2_only", kappa=kappa) == pytest.approx(
        D0_1_X2_EV, abs=5e-7
    )


@pytest.mark.parametrize("kappa", KAPPAS)
def test_first_rung_mixture_matches_oracle_4figs(kappa):
    assert d0_of_n(1, picture="statistical_mixture", kappa=kappa) == pytest.approx(
        D0_1_MIX_EV, abs=5e-7
    )


def test_default_picture_is_statistical_mixture():
    assert d0_of_n(1, kappa=1.0) == pytest.approx(D0_1_MIX_EV, abs=5e-7)


def test_first_rung_resolver_matches_module_d0_of_n():
    for picture in ("statistical_mixture", "x2_only", "cooling_relaxed"):
        assert first_rung_d0_eV(picture) == pytest.approx(
            d0_of_n(1, picture=picture, kappa=1.0), abs=1e-12
        )


# ---------------------------------------------------------------------------
# cooling_relaxed: provisional-average stub, asserted by ORDERING only
# (statistical_mixture < cooling_relaxed < x2_only), no pinned 4-fig value.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n", [1, 10, 21, 22])
def test_cooling_relaxed_strictly_between_pictures(n):
    mix = d0_of_n(n, picture="statistical_mixture", kappa=1.0)
    relaxed = d0_of_n(n, picture="cooling_relaxed", kappa=1.0)
    x2 = d0_of_n(n, picture="x2_only", kappa=1.0)
    assert mix < relaxed < x2


def test_cooling_relaxed_first_rung_is_the_average():
    # Provisional blend = arithmetic mean of mixture and X2 (decision 2026-06-30).
    assert first_rung_d0_eV("cooling_relaxed") == pytest.approx(
        0.5 * (D0_1_MIX_EV + D0_1_X2_EV), abs=1e-12
    )


# ---------------------------------------------------------------------------
# Integrated gate threshold Sigma(21): band membership + ~kappa-independence.
# ---------------------------------------------------------------------------
def _sigma21_ev(picture, kappa):
    return float(ladder_cumsum(N_STAR, picture=picture, kappa=kappa))


@pytest.mark.parametrize("kappa", KAPPAS)
def test_sigma21_mixture_in_band(kappa):
    lo, hi = MIX_SIGMA21_BAND_EV
    assert lo <= round(_sigma21_ev("statistical_mixture", kappa), 2) <= hi


@pytest.mark.parametrize("kappa", KAPPAS)
def test_sigma21_x2_in_band(kappa):
    lo, hi = X2_SIGMA21_BAND_EV
    assert lo <= round(_sigma21_ev("x2_only", kappa), 2) <= hi


@pytest.mark.parametrize("picture", ["statistical_mixture", "x2_only"])
def test_sigma21_kappa_decoupling(picture):
    # The gate threshold varies only ~11% across kappa in [0.3, 5]: it is
    # picture-set and pinnable ahead of the kappa fit.
    vals = [_sigma21_ev(picture, k) for k in KAPPAS]
    spread = (max(vals) - min(vals)) / np.mean(vals)
    assert spread <= 0.12


def test_sigma21_mixture_meets_drag_binding_order_of_magnitude():
    # Non-binding upper-bound cross-check vs E_bind_drag (942 cm^-1); rungs are
    # NOT calibrated to it. Representative first-shell sum band: 835-1571 cm^-1.
    sigma21_wavenumber = _sigma21_ev("statistical_mixture", 1.0) / EV_PER_WAVENUMBER
    assert 835.0 <= sigma21_wavenumber <= 1571.0
    assert sigma21_wavenumber > E_BIND_DRAG_WAVENUMBER  # |S| order, above the bound


# ---------------------------------------------------------------------------
# Cliff geometry: full-depth in-shell, collapse to floor just past n*.
# ---------------------------------------------------------------------------
def test_cliff_full_depth_in_shell_large_kappa():
    # For a sharp cliff D_0(n*) recovers (nearly) the full first-rung depth.
    d1 = first_rung_d0_eV("statistical_mixture")
    assert d0_of_n(N_STAR, picture="statistical_mixture", kappa=8.0) == pytest.approx(
        d1, rel=0.02
    )


def test_cliff_collapses_to_floor_just_past_n_star_large_kappa():
    d1 = first_rung_d0_eV("statistical_mixture")
    val = d0_of_n(N_STAR + 1, picture="statistical_mixture", kappa=8.0)
    # Lands within 5% of the floor's distance below the first rung.
    assert val < D_FLOOR_EV + 0.05 * (d1 - D_FLOOR_EV)


def test_sigma_is_centered_at_n_star_plus_half():
    # sigma(n* + 1/2) = 1/2 exactly for any kappa (the cliff centre).
    for kappa in KAPPAS:
        assert sigma(N_STAR + 0.5, kappa=kappa) == pytest.approx(0.5, abs=1e-12)


# ---------------------------------------------------------------------------
# Monotonicity + bounds.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", ["statistical_mixture", "x2_only", "cooling_relaxed"])
@pytest.mark.parametrize("kappa", KAPPAS)
def test_monotone_non_increasing_and_bounded(picture, kappa):
    n = np.arange(1, 31)
    d = d0_of_n(n, picture=picture, kappa=kappa)
    assert np.all(np.diff(d) <= 1e-12)            # non-increasing
    d1 = first_rung_d0_eV(picture)
    assert np.all(d >= D_FLOOR_EV - 1e-12)
    assert np.all(d <= d1 + 1e-12)


# ---------------------------------------------------------------------------
# Vectorisation: scalar-in -> float, array-in -> ndarray, elementwise parity.
# ---------------------------------------------------------------------------
def test_d0_of_n_scalar_returns_float():
    out = d0_of_n(5, picture="statistical_mixture", kappa=1.0)
    assert isinstance(out, float)


def test_d0_of_n_vectorized_matches_scalar():
    n = np.array([1, 5, 19, 21, 25])
    arr = d0_of_n(n, picture="statistical_mixture", kappa=1.3)
    assert arr.shape == n.shape
    for i, ni in enumerate(n):
        assert arr[i] == pytest.approx(
            d0_of_n(int(ni), picture="statistical_mixture", kappa=1.3), abs=1e-12
        )


def test_ladder_cumsum_scalar_and_vector_parity():
    assert isinstance(ladder_cumsum(10, picture="x2_only", kappa=1.0), float)
    ns = np.array([0, 1, 2, 21])
    cum = ladder_cumsum(ns, picture="x2_only", kappa=1.0)
    assert cum.shape == ns.shape
    assert cum[0] == pytest.approx(0.0, abs=1e-12)          # Sigma(0) = 0
    assert cum[1] == pytest.approx(
        d0_of_n(1, picture="x2_only", kappa=1.0), abs=1e-12
    )
    assert cum[2] == pytest.approx(
        d0_of_n(1, picture="x2_only", kappa=1.0)
        + d0_of_n(2, picture="x2_only", kappa=1.0),
        abs=1e-12,
    )


def test_gate_threshold_is_ladder_cumsum_alias():
    for n in (0, 1, 5, 21):
        assert gate_threshold(n, picture="statistical_mixture", kappa=1.1) == (
            ladder_cumsum(n, picture="statistical_mixture", kappa=1.1)
        )


# ---------------------------------------------------------------------------
# Tabulated fallback: round-trips a hand-built ladder.
# ---------------------------------------------------------------------------
def test_tabulated_ladder_round_trips_rungs_and_cumsum():
    rungs_eV = [0.0100, 0.0080, 0.0050, 0.0010]
    lad = tabulated_ladder(rungs_eV)
    assert lad.d0_of_n(1) == pytest.approx(0.0100, abs=1e-12)
    assert lad.d0_of_n(3) == pytest.approx(0.0050, abs=1e-12)
    assert lad.ladder_cumsum(0) == pytest.approx(0.0, abs=1e-12)
    assert lad.ladder_cumsum(2) == pytest.approx(0.0180, abs=1e-12)
    assert lad.ladder_cumsum(4) == pytest.approx(0.0240, abs=1e-12)


def test_tabulated_ladder_vectorized():
    rungs_eV = [0.01, 0.008, 0.005]
    lad = tabulated_ladder(rungs_eV)
    out = lad.d0_of_n(np.array([1, 2, 3]))
    assert np.allclose(out, rungs_eV)


# ---------------------------------------------------------------------------
# Fail-loud guards.
# ---------------------------------------------------------------------------
def test_unknown_picture_rejected_d0_of_n():
    with pytest.raises(ValueError, match="picture"):
        d0_of_n(1, picture="nonsense", kappa=1.0)


def test_unknown_picture_rejected_resolver():
    with pytest.raises(ValueError, match="picture"):
        first_rung_d0_eV("nonsense")


def test_tabulated_ladder_rejects_index_out_of_range():
    lad = tabulated_ladder([0.01, 0.008])
    with pytest.raises((ValueError, IndexError)):
        lad.d0_of_n(3)


def test_tabulated_ladder_rejects_empty():
    with pytest.raises(ValueError):
        tabulated_ladder([])


def test_tabulated_ladder_rejects_fractional_n():
    # Review fix (2026-07-02): mirror the Form-U ladder_cumsum fractional-n
    # guard -- a genuinely fractional occupancy must raise a loud ValueError,
    # not numpy's cryptic "arrays used as indices must be ... integer"
    # IndexError that the bare table lookup produced.
    lad = tabulated_ladder([0.01, 0.008, 0.005])
    with pytest.raises(ValueError, match="integer"):
        lad.d0_of_n(2.5)
    with pytest.raises(ValueError, match="integer"):
        lad.ladder_cumsum(2.5)
    with pytest.raises(ValueError, match="integer"):
        lad.d0_of_n(np.array([1.0, 2.5]))
    with pytest.raises(ValueError, match="integer"):
        lad.ladder_cumsum(np.array([0.0, 1.5]))


def test_tabulated_ladder_accepts_integer_valued_floats():
    # Same contract as the Form-U ladder_cumsum: integer-valued floats are
    # accepted and cast (2.0 is occupancy 2, not a fractional n).
    lad = tabulated_ladder([0.01, 0.008, 0.005])
    assert lad.d0_of_n(2.0) == pytest.approx(0.008, abs=1e-15)
    assert lad.ladder_cumsum(3.0) == pytest.approx(0.023, abs=1e-15)
    out = lad.ladder_cumsum(np.array([0.0, 2.0]))
    assert np.allclose(out, [0.0, 0.018], atol=1e-15)


# ---------------------------------------------------------------------------
# Extended invariants and robustness (review pass).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("picture", ["statistical_mixture", "x2_only", "cooling_relaxed"])
@pytest.mark.parametrize("kappa", KAPPAS)
def test_cumsum_difference_equals_rung(picture, kappa):
    # Sigma(n) - Sigma(n-1) == D_0(n) for every n (cumulative consistency).
    ns = np.arange(1, 26)
    cum = ladder_cumsum(ns, picture=picture, kappa=kappa)
    cum_prev = ladder_cumsum(ns - 1, picture=picture, kappa=kappa)
    rungs = d0_of_n(ns, picture=picture, kappa=kappa)
    assert np.allclose(cum - cum_prev, rungs, atol=1e-12)


def test_ladder_cumsum_strictly_increasing():
    # Every rung exceeds the floor (> 0), so Sigma(n) is strictly increasing.
    ns = np.arange(0, 26)
    cum = ladder_cumsum(ns, picture="x2_only", kappa=1.0)
    assert np.all(np.diff(cum) > 0)


@pytest.mark.parametrize("kappa", KAPPAS)
def test_sigma_monotone_non_decreasing_in_n(kappa):
    # sigma is monotone non-decreasing; at a sharp cliff it saturates exactly to
    # 0/1 in the tails (consecutive equal values), so non-decreasing -- not
    # strictly increasing -- is the true property.
    s = sigma(np.linspace(1.0, 30.0, 200), kappa=kappa)
    assert np.all(np.diff(s) >= 0)
    # strictly increasing across the cliff window where it is not saturated
    s_cliff = sigma(np.linspace(20.0, 23.0, 50), kappa=kappa)
    assert np.all(np.diff(s_cliff) > 0)


def test_d0_of_n_continuous_between_integer_rungs():
    # Fractional n interpolates monotonically between adjacent integer rungs.
    lo = d0_of_n(20, picture="statistical_mixture", kappa=1.0)
    mid = d0_of_n(20.5, picture="statistical_mixture", kappa=1.0)
    hi = d0_of_n(21, picture="statistical_mixture", kappa=1.0)
    assert hi <= mid <= lo


@pytest.mark.parametrize("kappa", KAPPAS)
def test_cooling_relaxed_ordering_across_kappa(kappa):
    for n in (1, 15, 21, 22):
        mix = d0_of_n(n, picture="statistical_mixture", kappa=kappa)
        rel = d0_of_n(n, picture="cooling_relaxed", kappa=kappa)
        x2 = d0_of_n(n, picture="x2_only", kappa=kappa)
        assert mix < rel < x2


def test_d0_finite_and_no_overflow_warning_large_kappa():
    # The sigmoid must not overflow to inf intermediates for a sharp cliff.
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("error")          # promote RuntimeWarning -> error
        d = d0_of_n(np.arange(1, 31), picture="x2_only", kappa=25.0)
        s1 = sigma(1.0, kappa=60.0)
    assert np.all(np.isfinite(d))
    assert 0.0 <= s1 <= 1.0


def test_d0_of_n_defensive_floor_guard(monkeypatch):
    # A picture whose first rung sits below the floor must raise (the inverted
    # ladder guard), not silently produce a non-dissipative rung.
    import i2_helium_md.physics.dissociation_ladder as ladder

    monkeypatch.setitem(ladder._FIRST_RUNG_EV, "broken", 0.5 * ladder.D_FLOOR_EV)
    with pytest.raises(ValueError, match="floor"):
        ladder.d0_of_n(1, picture="broken", kappa=1.0)


def test_ladder_cumsum_rejects_fractional_n():
    with pytest.raises(ValueError, match="integer"):
        ladder_cumsum(2.5, picture="x2_only", kappa=1.0)


def test_ladder_cumsum_cache_is_bit_identical_across_call_orders():
    # Slice-G review fix regression: ladder_cumsum serves Sigma from a cached,
    # padded prefix table (_sigma_prefix_table). Prefix sums do not depend on
    # later rungs, so the same n must return the bit-identical value regardless
    # of which table height was built first (small-then-large, large-then-small)
    # and must equal an uncached fresh construction of the same prefix.
    import i2_helium_md.physics.dissociation_ladder as ladder

    ladder._sigma_prefix_table.cache_clear()
    small_first = ladder_cumsum(5, picture="x2_only", kappa=2.0)
    large_after = ladder_cumsum(500, picture="x2_only", kappa=2.0)

    ladder._sigma_prefix_table.cache_clear()
    large_first = ladder_cumsum(500, picture="x2_only", kappa=2.0)
    small_after = ladder_cumsum(5, picture="x2_only", kappa=2.0)

    assert small_first == small_after   # exact float equality, no tolerance
    assert large_after == large_first

    # Uncached oracle: the pre-cache construction of the same prefix.
    rungs = np.atleast_1d(d0_of_n(np.arange(1, 6), picture="x2_only", kappa=2.0))
    assert small_first == float(np.concatenate(([0.0], np.cumsum(rungs)))[5])


def test_ladder_cumsum_cached_table_not_mutable_via_result():
    # The cached table is shared; results are copies (fancy indexing), so
    # mutating a returned array must not poison later calls.
    n = np.array([3, 4, 5])
    first = ladder_cumsum(n, picture="x2_only", kappa=2.0)
    first += 99.0
    second = ladder_cumsum(n, picture="x2_only", kappa=2.0)
    assert not np.array_equal(first, second)
    assert np.all(second < 99.0)
