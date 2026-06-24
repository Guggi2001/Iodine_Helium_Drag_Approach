"""Tests for i2_helium_md/physics/shell_schedule.py (Tier-1a Slice S).

Slice S is pure, closed-form kinematics, so every test is an oracle comparison
against the golden values in ``TIER1A_IMPLEMENTATION_PLAN.md`` Sec.10 -- no
reference-data files. The schedule turns the anchors ``(t*,21),(10,19),(14,14)``
into the 7 cold-shed events and the piecewise-linear ``n_bar(t)`` evaluator.

Oracle note on the Sec.10 absolute masses
------------------------------------------
The Sec.10 pre-shed mass column (210.955, 206.952, 202.954, 198.951, 194.949,
190.946, 186.939) is *internally* rounded: the ``n=19`` row is pinned to the
config ``m_eff = 202.953908`` (which uses the precise iodine mass 126.9045) while
the formula labels I+ as ``126.90``. The two differ by ~0.0045 amu, so the
interior rows (n=16..19) sit ~0.004 amu off the pure ``126.90 + n*4.0026`` formula
this module uses. That residual is documented I-mass rounding, **not** a code bug:
the **kick factors and the telescoping product are ratios and cancel it exactly**,
so they are the real oracle (asserted tight); the absolute masses are asserted
against the formula (tight) and against the Sec.10 table (loose, ~5e-3 amu).
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU
from i2_helium_md.physics.shell_schedule import (
    ANCHOR_N_END,
    ANCHOR_N_MID,
    ANCHOR_N_START,
    ANCHOR_T_END_PS,
    ANCHOR_T_MID_PS,
    NUM_SHED_EVENTS,
    build_shell_schedule,
    complex_mass_amu,
)


# Sweep values used across the suite (TIER1A run matrix: t* in {0.5, 5, 9} ps).
T_STARS = (0.5, 5.0, 9.0)

# Sec.10 golden tables (n_before order 21->...->15).
SEG2_TIMES_PS = (10.4, 11.2, 12.0, 12.8, 13.6)            # t*-independent
ORACLE_KICKS = (1.0193, 1.0197, 1.0201, 1.0205, 1.0210, 1.0214, 1.0219)
ORACLE_PRESHED_MASSES = (
    210.955, 206.952, 202.954, 198.951, 194.949, 190.946, 186.939,
)
TELESCOPING = 1.1532


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------
class TestStructure:
    @pytest.mark.parametrize("t_star", T_STARS)
    def test_event_count_and_shell_sequence(self, t_star):
        sched = build_shell_schedule(t_star)
        assert len(sched.events) == NUM_SHED_EVENTS == 7
        # n strictly decreases by exactly one per shed, 21 -> 14.
        assert [e.n_before for e in sched.events] == list(range(21, 14, -1))
        assert all(e.n_after == e.n_before - 1 for e in sched.events)

    @pytest.mark.parametrize("t_star", T_STARS)
    def test_times_monotone_and_in_window(self, t_star):
        sched = build_shell_schedule(t_star)
        times = [e.time_ps for e in sched.events]
        assert all(b > a for a, b in zip(times, times[1:]))          # strictly increasing
        assert all(t_star < t <= ANCHOR_T_END_PS for t in times)     # in (t*, 14]


# ---------------------------------------------------------------------------
# Event timing -- analytic identities (tight)
# ---------------------------------------------------------------------------
class TestTiming:
    @pytest.mark.parametrize("t_star", (0.5, 9.0))
    def test_segment2_times_are_t_star_independent(self, t_star):
        sched = build_shell_schedule(t_star)
        seg2 = [e.time_ps for e in sched.events if e.crossing < ANCHOR_N_MID]
        # Five segment-2 events at {10.4..13.6} ps regardless of t*. Tight:
        # closed-form solve of a linear segment, only float round-off separates.
        np.testing.assert_allclose(seg2, SEG2_TIMES_PS, rtol=1e-12, atol=1e-12)

    @pytest.mark.parametrize("t_star", T_STARS)
    def test_segment1_times_match_closed_form(self, t_star):
        sched = build_shell_schedule(t_star)
        seg1 = [e.time_ps for e in sched.events if e.crossing > ANCHOR_N_MID]
        # The two segment-1 events at t* + {0.25, 0.75}*(10 - t*); both < 10.
        expected = [t_star + f * (ANCHOR_T_MID_PS - t_star) for f in (0.25, 0.75)]
        np.testing.assert_allclose(seg1, expected, rtol=1e-12, atol=1e-12)
        assert all(t < ANCHOR_T_MID_PS for t in seg1)


# ---------------------------------------------------------------------------
# Masses and kicks vs Sec.10 oracle
# ---------------------------------------------------------------------------
class TestMassesAndKicks:
    def test_masses_match_formula_exactly(self):
        sched = build_shell_schedule(0.5)
        for e in sched.events:
            # The module's own rule, machine precision.
            assert e.mass_before_amu == pytest.approx(
                MASS_I_ION_AMU + e.n_before * MASS_HE_AMU, rel=1e-12
            )
            assert e.mass_after_amu == pytest.approx(
                MASS_I_ION_AMU + e.n_after * MASS_HE_AMU, rel=1e-12
            )

    def test_preshed_masses_match_oracle_within_rounding(self):
        sched = build_shell_schedule(0.5)
        got = [e.mass_before_amu for e in sched.events]
        # Loose (5e-3 amu): the Sec.10 table pins n=19 to the precise-iodine m_eff
        # while the formula rounds I+ to 126.90 (see module docstring); the ~0.0045
        # amu interior offset is that documented rounding, not a code error.
        np.testing.assert_allclose(got, ORACLE_PRESHED_MASSES, atol=5e-3)

    def test_kick_factors_match_oracle_to_four_figures(self):
        sched = build_shell_schedule(0.5)
        got = [e.kick_factor for e in sched.events]
        # Tight: kicks are ratios m/(m-m_He) and cancel the I-mass rounding, so
        # they reproduce the Sec.10 column to its 4 printed decimals.
        np.testing.assert_allclose(got, ORACLE_KICKS, atol=5e-5)

    def test_kick_is_momentum_conserving_ratio(self):
        sched = build_shell_schedule(0.5)
        for e in sched.events:
            assert e.kick_factor == pytest.approx(
                e.mass_before_amu / (e.mass_before_amu - MASS_HE_AMU), rel=1e-12
            )


# ---------------------------------------------------------------------------
# Telescoping invariant -- the verdict-robustness check
# ---------------------------------------------------------------------------
class TestTelescoping:
    @pytest.mark.parametrize("t_star", T_STARS)
    @pytest.mark.parametrize("frac", (0.3, 0.5, 0.7))
    def test_kick_product_equals_endpoint_mass_ratio(self, t_star, frac):
        sched = build_shell_schedule(t_star, crossing_fraction=frac)
        product = float(np.prod([e.kick_factor for e in sched.events]))
        ratio = complex_mass_amu(ANCHOR_N_START) / complex_mass_amu(ANCHOR_N_END)
        # Algebraic identity (telescoping product == endpoint ratio): tight, and
        # invariant to t* and to the half-integer tie-break (only times move).
        assert product == pytest.approx(ratio, rel=1e-12)
        assert product == pytest.approx(TELESCOPING, abs=5e-5)
        assert len(sched.events) == NUM_SHED_EVENTS  # tie-break does not change count


# ---------------------------------------------------------------------------
# n_bar(t) evaluator
# ---------------------------------------------------------------------------
class TestNBarEvaluator:
    @pytest.mark.parametrize("t_star", T_STARS)
    def test_anchor_waypoints(self, t_star):
        sched = build_shell_schedule(t_star)
        assert sched.n_bar(t_star) == pytest.approx(ANCHOR_N_START, abs=1e-12)
        assert sched.n_bar(0.0) == pytest.approx(ANCHOR_N_START, abs=1e-12)
        assert sched.n_bar(ANCHOR_T_MID_PS) == pytest.approx(ANCHOR_N_MID, abs=1e-12)
        assert sched.n_bar(ANCHOR_T_END_PS) == pytest.approx(ANCHOR_N_END, abs=1e-12)

    @pytest.mark.parametrize("t_star", T_STARS)
    def test_flat_tail_past_14ps(self, t_star):
        sched = build_shell_schedule(t_star)
        for t in (14.0, 15.0, 20.0, 100.0):
            assert sched.n_bar(t) == pytest.approx(ANCHOR_N_END, abs=1e-12)

    def test_continuity_at_segment_joins(self):
        sched = build_shell_schedule(0.5)
        eps = 1e-9
        for t_join in (sched.t_star_ps, ANCHOR_T_MID_PS, ANCHOR_T_END_PS):
            assert sched.n_bar(t_join - eps) == pytest.approx(
                sched.n_bar(t_join + eps), abs=1e-6
            )

    def test_monotone_non_increasing(self):
        sched = build_shell_schedule(0.5)
        t = np.linspace(0.0, 20.0, 2001)
        n = sched.n_bar(t)
        assert np.all(np.diff(n) <= 1e-12)

    def test_vectorized_shape_preserved(self):
        sched = build_shell_schedule(0.5)
        # scalar in -> python float out
        assert isinstance(sched.n_bar(5.0), float)
        # array in -> same-shape ndarray out
        t = np.array([0.0, 5.0, 10.0, 12.0, 14.0, 18.0])
        out = sched.n_bar(t)
        assert isinstance(out, np.ndarray) and out.shape == t.shape


# ---------------------------------------------------------------------------
# Integer shell count n(t) -- the physical staircase
# ---------------------------------------------------------------------------
class TestIntegerShellCount:
    @pytest.mark.parametrize("t_star", T_STARS)
    def test_endpoints(self, t_star):
        sched = build_shell_schedule(t_star)
        # 21 before any shed (incl. the held onset region), 14 after the last.
        assert sched.n_of_t(0.0) == ANCHOR_N_START
        assert sched.n_of_t(t_star) == ANCHOR_N_START          # nothing sheds at/before t*
        for t in (14.0, 16.0, 20.0):
            assert sched.n_of_t(t) == ANCHOR_N_END

    @pytest.mark.parametrize("t_star", T_STARS)
    def test_steps_down_by_one_at_each_event(self, t_star):
        sched = build_shell_schedule(t_star)
        eps = 1e-9
        for e in sched.events:
            # Just before -> n_before; at/just after the fire time -> n_after.
            assert sched.n_of_t(e.time_ps - eps) == e.n_before
            assert sched.n_of_t(e.time_ps) == e.n_after
            assert sched.n_of_t(e.time_ps + eps) == e.n_after

    def test_always_integer_and_in_range(self):
        sched = build_shell_schedule(0.5)
        t = np.linspace(0.0, 20.0, 4001)
        n = sched.n_of_t(t)
        assert np.issubdtype(n.dtype, np.integer)
        assert n.min() == ANCHOR_N_END and n.max() == ANCHOR_N_START
        assert np.all(np.diff(n) <= 0)                          # monotone non-increasing

    def test_scalar_returns_python_int_array_preserves_shape(self):
        sched = build_shell_schedule(0.5)
        assert isinstance(sched.n_of_t(5.0), int)
        t = np.array([0.0, 11.0, 14.0, 20.0])
        out = sched.n_of_t(t)
        assert isinstance(out, np.ndarray) and out.shape == t.shape

    def test_count_brackets_n_bar(self):
        # The integer count never drifts more than half a shell from the loss
        # curve it discretizes (a sanity tie between the two quantities).
        sched = build_shell_schedule(0.5)
        t = np.linspace(0.0, 20.0, 4001)
        assert np.all(np.abs(sched.n_of_t(t) - sched.n_bar(t)) <= 0.5 + 1e-9)


# ---------------------------------------------------------------------------
# Fail-loud
# ---------------------------------------------------------------------------
class TestValidation:
    @pytest.mark.parametrize("bad", (-0.1, 10.0, 11.0, 14.0))
    def test_t_star_out_of_range_raises(self, bad):
        with pytest.raises(ValueError, match="t_star_ps"):
            build_shell_schedule(bad)

    @pytest.mark.parametrize("bad", (0.0, 1.0, -0.5, 1.5))
    def test_crossing_fraction_out_of_range_raises(self, bad):
        with pytest.raises(ValueError, match="crossing_fraction"):
            build_shell_schedule(0.5, crossing_fraction=bad)
