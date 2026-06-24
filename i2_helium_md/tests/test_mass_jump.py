"""Tests for i2_helium_md/physics/mass_jump.py (Tier-1a Slice M, SQ2).

Slice M is the pure, stateless **cold-shed mass-jump operator**: the
momentum-conserving velocity reset a single He shed applies, the post-jump mass
``m+`` the integrator (SQ3) needs, and the exact reduced-mass energy increment the
four-term ledger (Slice B) needs. No SimConfig, no integrator, no drag -- every
test here is an oracle / analytic-identity comparison (no reference-data files),
mirroring ``test_shell_schedule.py``.

Locked physics (``TIER1A_IMPLEMENTATION_PLAN.md`` Sec.2):

* reset           ``v+ = m/(m - m_He) * v-``         (He leaves at rest; m*v invariant)
* post-jump mass  ``m+ = m - m_He``
* ledger term     ``dE_mass_transfer = -0.5 * (m*m_He)/(m - m_He) * |v-|^2``

Units are mechanical amu throughout: ``m`` in amu, ``v`` in A/ps, energy in
``amu*A^2/ps^2`` (matching ``baoab.py``'s ``dE_dissip`` and ``shell_schedule.py``);
no eV conversion lives in Slice M.
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU
from i2_helium_md.physics.mass_jump import (
    ShedResult,
    apply_shed,
    cold_shed,
    kick_factor,
)
from i2_helium_md.physics.shell_schedule import (
    ANCHOR_N_END,
    ANCHOR_N_START,
    build_shell_schedule,
    complex_mass_amu,
)

# Sec.10 golden tables (n_before order 21 -> 15) -- the shared Tier-1a oracle.
ORACLE_KICKS = (1.0193, 1.0197, 1.0201, 1.0205, 1.0210, 1.0214, 1.0219)
TELESCOPING = 1.1532

# A generic non-axis-aligned velocity [A/ps] used across the reset/energy tests.
V_MINUS = np.array([3.0, -4.0, 12.0])           # |v| = 13 A/ps exactly
M_EFF_AMU = 202.953908                           # config m_eff (= 19-He complex)


# ---------------------------------------------------------------------------
# Momentum-conserving reset
# ---------------------------------------------------------------------------
class TestReset:
    @pytest.mark.parametrize("n_before", range(ANCHOR_N_START, ANCHOR_N_END, -1))
    def test_momentum_conserved_to_machine_precision(self, n_before):
        m = complex_mass_amu(n_before)
        res = cold_shed(V_MINUS, m)
        # (m - m_He) * v+ == m * v- elementwise (He shed at rest).
        np.testing.assert_allclose(
            res.m_plus_amu * res.v_plus, m * V_MINUS, rtol=0.0, atol=1e-12
        )

    def test_post_jump_mass_is_one_he_lighter(self):
        m = complex_mass_amu(21)
        res = cold_shed(V_MINUS, m)
        assert res.m_plus_amu == pytest.approx(m - MASS_HE_AMU, abs=1e-12)

    def test_direction_preserved_pure_scalar_scaling(self):
        res = cold_shed(V_MINUS, complex_mass_amu(21))
        # v+ is parallel to v- (zero cross product) and strictly longer (kick > 1).
        assert np.allclose(np.cross(res.v_plus, V_MINUS), 0.0, atol=1e-12)
        assert np.linalg.norm(res.v_plus) > np.linalg.norm(V_MINUS)

    def test_returns_independent_array(self):
        v = V_MINUS.copy()
        res = cold_shed(v, complex_mass_amu(21))
        res.v_plus[0] = 999.0          # mutating the result must not touch the input
        assert v[0] == V_MINUS[0]


# ---------------------------------------------------------------------------
# Kick factor vs the Sec.10 oracle AND the delivered Slice S
# ---------------------------------------------------------------------------
class TestKickFactor:
    def test_kick_matches_section10_oracle_4_figures(self):
        kicks = [kick_factor(complex_mass_amu(n))
                 for n in range(ANCHOR_N_START, ANCHOR_N_END, -1)]
        # 4 printed figures in the Sec.10 table.
        np.testing.assert_allclose(kicks, ORACLE_KICKS, atol=5e-5)

    @pytest.mark.parametrize("t_star", (0.5, 5.0, 9.0))
    def test_kick_agrees_with_shell_schedule_events(self, t_star):
        # Cross-check against Slice S without depending on its internals: the
        # operator's kick must equal each ShedEvent.kick_factor exactly.
        sched = build_shell_schedule(t_star)
        for ev in sched.events:
            k = kick_factor(ev.mass_before_amu)
            assert k == pytest.approx(ev.kick_factor, rel=0.0, abs=1e-12)
            # And the reset reproduces that kick on a real velocity.
            res = cold_shed(V_MINUS, ev.mass_before_amu)
            np.testing.assert_allclose(res.v_plus, ev.kick_factor * V_MINUS, atol=1e-12)

    def test_telescoping_product_is_endpoint_mass_ratio(self):
        # Chain the sheds 21 -> 14 on one velocity; cumulative speed boost telescopes
        # to m(21)/m(14) = 1.1532, independent of the per-event path.
        v = V_MINUS.copy()
        for n_before in range(ANCHOR_N_START, ANCHOR_N_END, -1):
            v = cold_shed(v, complex_mass_amu(n_before)).v_plus
        boost = np.linalg.norm(v) / np.linalg.norm(V_MINUS)
        ratio = complex_mass_amu(ANCHOR_N_START) / complex_mass_amu(ANCHOR_N_END)
        assert boost == pytest.approx(ratio, rel=1e-12)
        assert boost == pytest.approx(TELESCOPING, abs=5e-4)


# ---------------------------------------------------------------------------
# Energy increment (exact reduced-mass form, not heavy-ion)
# ---------------------------------------------------------------------------
class TestEnergyIncrement:
    @pytest.mark.parametrize("n_before", range(ANCHOR_N_START, ANCHOR_N_END, -1))
    def test_increment_is_negative_of_exact_ke_rise(self, n_before):
        m = complex_mass_amu(n_before)
        res = cold_shed(V_MINUS, m)
        speed_sq = float(V_MINUS @ V_MINUS)
        # Closed-form KE rise of the reset (Sec.2 derivation).
        ke_rise = 0.5 * (m * MASS_HE_AMU) / (m - MASS_HE_AMU) * speed_sq
        assert res.dE_mass_transfer == pytest.approx(-ke_rise, rel=1e-12)
        # Cross-check against the direct before/after kinetic-energy difference.
        ke_before = 0.5 * m * speed_sq
        ke_after = 0.5 * res.m_plus_amu * float(res.v_plus @ res.v_plus)
        assert (ke_after - ke_before) == pytest.approx(ke_rise, rel=1e-12)

    def test_increment_is_non_positive(self):
        for n in range(ANCHOR_N_START, ANCHOR_N_END, -1):
            assert cold_shed(V_MINUS, complex_mass_amu(n)).dE_mass_transfer <= 0.0

    def test_reduced_mass_form_diverges_from_heavy_ion_by_about_3pct_at_n1(self):
        # n=1 complex (I+ + 1 He); the exact reduced-mass coefficient
        # m*m_He/(m - m_He) exceeds the heavy-ion m_He by ~3% -- the Sec.2 flag
        # that the heavy-ion 0.5*m_He*v^2 approximation is NOT used here.
        m = MASS_I_ION_AMU + MASS_HE_AMU          # n = 1
        speed_sq = float(V_MINUS @ V_MINUS)
        exact_rise = -cold_shed(V_MINUS, m).dE_mass_transfer
        heavy_ion_rise = 0.5 * MASS_HE_AMU * speed_sq
        ratio = exact_rise / heavy_ion_rise
        assert 1.02 < ratio < 1.04                # ~3% (1.0316), documented divergence

    def test_heavy_ion_limit_recovered_for_large_mass(self):
        # As m -> inf the exact coefficient -> m_He (heavy-ion limit).
        speed_sq = float(V_MINUS @ V_MINUS)
        exact_rise = -cold_shed(V_MINUS, 1.0e6).dE_mass_transfer
        heavy_ion_rise = 0.5 * MASS_HE_AMU * speed_sq
        assert exact_rise == pytest.approx(heavy_ion_rise, rel=1e-4)


# ---------------------------------------------------------------------------
# Mode wrapper: fixed (null) vs anchored_discrete (reset)
# ---------------------------------------------------------------------------
class TestModeWrapper:
    def test_fixed_mode_is_a_no_op(self):
        res = apply_shed(V_MINUS, M_EFF_AMU, mode="fixed")
        np.testing.assert_array_equal(res.v_plus, V_MINUS)   # no reset
        assert res.m_plus_amu == M_EFF_AMU                   # mass held (m_eff)
        assert res.dE_mass_transfer == 0.0                   # zero defect

    def test_anchored_discrete_matches_cold_shed(self):
        m = complex_mass_amu(21)
        a = apply_shed(V_MINUS, m, mode="anchored_discrete")
        b = cold_shed(V_MINUS, m)
        np.testing.assert_array_equal(a.v_plus, b.v_plus)
        assert a.m_plus_amu == b.m_plus_amu
        assert a.dE_mass_transfer == b.dE_mass_transfer

    def test_unknown_mode_raises(self):
        with pytest.raises(ValueError, match="mode"):
            apply_shed(V_MINUS, M_EFF_AMU, mode="scenario_A_accretion")


# ---------------------------------------------------------------------------
# Fail-loud guards (CLAUDE.md rule 4)
# ---------------------------------------------------------------------------
class TestGuards:
    def test_mass_not_greater_than_he_raises(self):
        with pytest.raises(ValueError, match="m_he"):
            cold_shed(V_MINUS, MASS_HE_AMU)            # m+ would be 0
        with pytest.raises(ValueError, match="m_he"):
            cold_shed(V_MINUS, 0.5 * MASS_HE_AMU)      # m+ would be negative

    def test_non_positive_he_mass_raises(self):
        with pytest.raises(ValueError):
            cold_shed(V_MINUS, complex_mass_amu(21), m_he_amu=0.0)

    def test_non_finite_inputs_raise(self):
        with pytest.raises(ValueError):
            cold_shed(np.array([1.0, np.nan, 0.0]), complex_mass_amu(21))
        with pytest.raises(ValueError):
            cold_shed(V_MINUS, np.inf)

    def test_returns_shedresult_type(self):
        assert isinstance(cold_shed(V_MINUS, complex_mass_amu(21)), ShedResult)
