"""Tier-1a variable-mass integrator wiring unit tests.

The driver-facing ``anchored_discrete`` pre-step now uses continuous-velocity
shedding: the shell schedule fires at the same times, the mass drops by one He,
and the velocity components are copied through unchanged. The vectorized
cold-shed helper remains covered here only as a diagnostic bound.
"""

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.physics.constants import MASS_HE_AMU, U
from i2_helium_md.physics.mass_jump import (
    cold_shed,
    cold_shed_velocity_components,
    continuous_velocity_shed,
    continuous_velocity_shed_components,
)
from i2_helium_md.physics.shell_schedule import build_shell_schedule, complex_mass_amu
from i2_helium_md.simulation.ion_propagation_step import IonStepState, shed_step


def _state(time_ps, *, vx, vy=0.0, vz=0.0, mass_amu, n_atoms=2):
    """Minimal IonStepState with uniform velocity/mass across ``n_atoms``."""
    z = np.zeros(n_atoms)
    return IonStepState(
        x=z.copy(), y=z.copy(), z=z.copy(),
        vx=np.full(n_atoms, vx), vy=np.full(n_atoms, vy), vz=np.full(n_atoms, vz),
        mass_kg=np.full(n_atoms, mass_amu * U),
        E_kin_eV=z.copy(), E_pot_eV=z.copy(), E_dissip_eV=z.copy(),
        E_mass_transfer_eV=z.copy(),
        number_of_collisions=np.zeros(n_atoms, dtype=int),
        time_ps=time_ps,
    )


# ---------------------------------------------------------------------------
# cold_shed_velocity_components vs the scalar cold_shed diagnostic oracle
# ---------------------------------------------------------------------------
class TestVectorizedColdShed:
    def test_matches_scalar_oracle_per_atom(self):
        m = complex_mass_amu(21)
        vx = np.array([1.0, -2.0])
        vy = np.array([0.5, 0.0])
        vz = np.array([0.0, 3.0])
        vxp, vyp, vzp, m_plus, dE = cold_shed_velocity_components(vx, vy, vz, m)
        for i in range(2):
            oracle = cold_shed([vx[i], vy[i], vz[i]], m)
            assert m_plus == pytest.approx(oracle.m_plus_amu)
            np.testing.assert_allclose([vxp[i], vyp[i], vzp[i]], oracle.v_plus)
            assert dE[i] == pytest.approx(oracle.dE_mass_transfer)

    def test_defect_non_positive_and_mass_drops_one_he(self):
        m = complex_mass_amu(17)
        vxp, vyp, vzp, m_plus, dE = cold_shed_velocity_components(
            np.array([2.0]), np.array([0.0]), np.array([0.0]), m,
        )
        assert m_plus == pytest.approx(complex_mass_amu(16))
        assert np.all(dE <= 0.0)

    def test_guard_rejects_mass_below_he(self):
        with pytest.raises(ValueError, match="m_he|pre-shed"):
            cold_shed_velocity_components(
                np.array([1.0]), np.array([0.0]), np.array([0.0]), 2.0,
            )


# ---------------------------------------------------------------------------
# continuous_velocity_shed_components vs the scalar production oracle
# ---------------------------------------------------------------------------
class TestVectorizedContinuousVelocityShed:
    def test_matches_scalar_oracle_per_atom(self):
        m = complex_mass_amu(21)
        vx = np.array([1.0, -2.0])
        vy = np.array([0.5, 0.0])
        vz = np.array([0.0, 3.0])
        vxp, vyp, vzp, m_plus, dE = continuous_velocity_shed_components(vx, vy, vz, m)
        for i in range(2):
            oracle = continuous_velocity_shed([vx[i], vy[i], vz[i]], m)
            assert m_plus == pytest.approx(oracle.m_plus_amu)
            np.testing.assert_array_equal([vxp[i], vyp[i], vzp[i]], oracle.v_plus)
            assert dE[i] == pytest.approx(oracle.dE_mass_transfer)

    def test_velocity_unchanged_positive_transfer_and_mass_drops_one_he(self):
        m = complex_mass_amu(17)
        vx = np.array([2.0])
        vy = np.array([0.0])
        vz = np.array([0.0])
        vxp, vyp, vzp, m_plus, dE = continuous_velocity_shed_components(vx, vy, vz, m)
        assert m_plus == pytest.approx(complex_mass_amu(16))
        np.testing.assert_array_equal(vxp, vx)
        np.testing.assert_array_equal(vyp, vy)
        np.testing.assert_array_equal(vzp, vz)
        assert np.all(dE > 0.0)

    def test_batch_vectorized_energy_uses_removed_count(self):
        m = complex_mass_amu(21)
        vx = np.array([2.0])
        vy = np.array([0.0])
        vz = np.array([0.0])
        vxp, vyp, vzp, m_plus, dE = continuous_velocity_shed_components(
            vx, vy, vz, m, n_removed=21,
        )

        np.testing.assert_array_equal(vxp, vx)
        np.testing.assert_array_equal(vyp, vy)
        np.testing.assert_array_equal(vzp, vz)
        assert m_plus == pytest.approx(complex_mass_amu(0))
        assert dE[0] == pytest.approx(0.5 * 21 * MASS_HE_AMU * 4.0)

    def test_guard_rejects_mass_below_he(self):
        with pytest.raises(ValueError, match="m_he|pre-shed"):
            continuous_velocity_shed_components(
                np.array([1.0]), np.array([0.0]), np.array([0.0]), 2.0,
            )


# ---------------------------------------------------------------------------
# shed_step -- the driver-facing continuous-velocity pre-step
# ---------------------------------------------------------------------------
class TestShedStep:
    def test_no_fire_before_window_returns_unchanged(self):
        sched = build_shell_schedule(5.0)
        st = _state(0.0, vx=1.0, mass_amu=complex_mass_amu(21))
        new, idx = shed_step(st, sched, 0, dt=0.01)
        assert idx == 0
        assert new is st  # identity: nothing happened

    def test_fires_when_event_in_window(self):
        sched = build_shell_schedule(5.0)
        e0 = sched.events[0]
        st = _state(e0.time_ps - 0.005, vx=1.0, mass_amu=e0.mass_before_amu)
        new, idx = shed_step(st, sched, 0, dt=0.01)
        assert idx == 1
        # SQ3: post-jump mass = pre-shed mass - one He
        assert new.mass_kg[0] == pytest.approx(e0.mass_after_amu * U)
        # Continuous-velocity shed: no event-local speed kick.
        assert new.vx[0] == pytest.approx(1.0)
        assert new.vy[0] == pytest.approx(0.0)
        assert new.vz[0] == pytest.approx(0.0)
        # Co-moving He kinetic energy is booked positive (in eV).
        assert np.all(new.E_mass_transfer_eV > 0.0)

    def test_at_most_one_shed_per_call_under_dense_window(self):
        # A window wide enough to straddle two seg-2 events (0.8 ps apart) must
        # still fire only the next pending one (the <=1/step rule).
        sched = build_shell_schedule(5.0)
        assert sched.events[2].time_ps < 10.5 < sched.events[3].time_ps < 11.4
        st = _state(10.3, vx=1.0, mass_amu=sched.events[2].mass_before_amu)
        new, idx = shed_step(st, sched, 2, dt=1.0)
        assert idx == 3  # exactly one fired despite two events in (10.3, 11.3]

    def test_exhausted_schedule_is_noop(self):
        sched = build_shell_schedule(5.0)
        st = _state(20.0, vx=1.0, mass_amu=complex_mass_amu(14))
        new, idx = shed_step(st, sched, len(sched.events), dt=0.01)
        assert idx == len(sched.events)
        assert new is st


def _run_schedule(sched, dt, *, vx0, mass0, t_end=14.5):
    """Drive shed_step over a uniform time grid; return (n_fires, final_state)."""
    st = _state(0.0, vx=vx0, mass_amu=mass0)
    idx = 0
    fires = 0
    while st.time_ps < t_end:
        new, nidx = shed_step(st, sched, idx, dt)
        if nidx != idx:
            fires += 1
        st = replace(new, time_ps=new.time_ps + dt)
        idx = nidx
    return fires, st


class TestScheduleDriven:
    @pytest.mark.parametrize("dt", [0.001, 0.01, 0.1])
    def test_total_shed_count_invariant_to_dt(self, dt):
        # Jump-step measure-zero: refining dt leaves the event count fixed at 7.
        sched = build_shell_schedule(5.0)
        fires, st = _run_schedule(sched, dt, vx0=1.0, mass0=complex_mass_amu(21))
        assert fires == 7
        assert st.mass_kg[0] == pytest.approx(complex_mass_amu(14) * U)

    def test_no_telescoping_speed_boost_on_driver_path(self):
        # Force-free (no integration): composing the 7 production sheds leaves
        # velocity unchanged; the old telescoping boost is cold-shed-bound only.
        sched = build_shell_schedule(0.5)
        fires, st = _run_schedule(sched, 0.001, vx0=2.0, mass0=complex_mass_amu(21))
        assert fires == 7
        ratio = complex_mass_amu(21) / complex_mass_amu(14)
        assert ratio == pytest.approx(1.1532, abs=1e-4)
        assert st.vx[0] == pytest.approx(2.0)
