"""Tests for the Tier-2 Phase-C biphasic generative driver (Slice G).

Covers, in build order:

* **Plumbing** -- the ``IonStepState.n_shell`` genuine-state field, the v7 column
  read seam, the ``write_ion_state_to_checkpoint_column`` scenario branch (biphasic
  stores ``state.n_shell``; fixed/anchored keep the ``rint``-from-mass derivation
  byte-identical), and ``_check_drag_scope`` accepting ``biphasic``.
* **biphasic_step** -- the per-ion composition: K2 cooling drain (E_int -> E_dissip),
  the evaporation-then-pickup one-event order, the S1/K1 E_int increments, the
  E_mass_transfer eV conversion, and the ``m == m_I+ + n*m_He`` invariant.
* **Closure** -- the 5-term invariant to machine precision on a force-free run and on
  isolated pickup / shed events (the driver arm folds ``e_bind_pair(n)`` into E_pot).
* **Regression** -- ``fixed`` dispatch stays bit-identical (biphasic never entered).

The accounting oracle is ``scratchpad/sliceG_accounting.md`` (the worked 5-term
per-channel closure).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.config import SimConfig
from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU, U
from i2_helium_md.physics.shell_schedule import ANCHOR_N_START, complex_mass_amu
from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.ion_initial_state import build_initial_ion_state
from i2_helium_md.simulation.ion_propagation_step import (
    IonStepState,
    biphasic_step,
    ion_state_from_checkpoint_column,
    write_ion_state_to_checkpoint_column,
    _amu_ang2_ps2_to_eV,
    _check_drag_scope,
    _E_kin_eV,
)
from i2_helium_md.simulation.neutral import run_neutral_propagation

# Shared "minimal valid biphasic SimConfig" helper -- single source in the
# config-guard suite (CLAUDE.md principle 1); the local wrapper only flips the
# pickup channel off by default for these unit tests.
from tests.test_biphasic_config import _biphasic_cfg as _base_biphasic_cfg


@pytest.fixture(scope="module")
def fixed_ion_ckpt():
    """A tiny finished-init fixed-scenario ion checkpoint (v7, n_shell populated)."""
    cfg = single_pulse_N2000(num_molecules=6, seed=7)
    cfg = replace(cfg, t_max_neutral=0.1, dt_neutral=0.01, dt_ion=0.01)
    neutral = run_neutral_propagation(cfg, verbose=False)
    ion = build_initial_ion_state(cfg, neutral, num_steps_ion=8)
    return cfg, ion


# ===========================================================================
# Plumbing: n_shell genuine state + column seams + scope flip
# ===========================================================================
class TestNShellColumnSeams:
    def test_read_seam_populates_n_shell(self, fixed_ion_ckpt):
        _, ion_shared = fixed_ion_ckpt
        ion = replace(ion_shared, n_shell=ion_shared.n_shell.copy())
        two_N = ion.n_shell.shape[0]
        ion.n_shell[:, 0] = np.arange(two_N, dtype=float)  # known distinct values
        state = ion_state_from_checkpoint_column(ion, 0)
        np.testing.assert_array_equal(state.n_shell, ion.n_shell[:, 0])

    def test_writer_fixed_derives_n_shell_via_rint(self, fixed_ion_ckpt):
        _, ion_shared = fixed_ion_ckpt
        ion = replace(ion_shared, n_shell=ion_shared.n_shell.copy())
        two_N = ion.n_shell.shape[0]
        # A state at n=19 mass, but carrying a deliberately-wrong n_shell array.
        m19 = complex_mass_amu(19) * U
        state = ion_state_from_checkpoint_column(ion, 0)
        state = replace(
            state,
            mass_kg=np.full(two_N, m19),
            n_shell=np.full(two_N, 999.0),  # must be ignored off the biphasic path
        )
        write_ion_state_to_checkpoint_column(state, ion, 1, mass_scenario="fixed")
        expected = np.rint((m19 / U - MASS_I_ION_AMU) / MASS_HE_AMU)
        np.testing.assert_array_equal(ion.n_shell[:, 1], np.full(two_N, expected))

    def test_writer_biphasic_none_n_shell_raises(self, fixed_ion_ckpt):
        # Slice-G review fix: a biphasic write without genuine n_shell state must
        # fail loudly, not silently fall back to the rint-from-mass derivation
        # (which would mask exactly the channel/reset drift n_shell exists to expose).
        _, ion_shared = fixed_ion_ckpt
        ion = replace(ion_shared, n_shell=ion_shared.n_shell.copy())
        state = ion_state_from_checkpoint_column(ion, 0)
        state = replace(state, n_shell=None)
        with pytest.raises(ValueError, match="n_shell"):
            write_ion_state_to_checkpoint_column(
                state, ion, 1, mass_scenario="biphasic",
            )

    def test_writer_biphasic_stores_state_n_shell(self, fixed_ion_ckpt):
        _, ion_shared = fixed_ion_ckpt
        ion = replace(ion_shared, n_shell=ion_shared.n_shell.copy())
        two_N = ion.n_shell.shape[0]
        # Mass says n=21 but n_shell state says something else -> biphasic stores state.
        genuine_n = np.arange(two_N, dtype=float) + 3.0
        state = ion_state_from_checkpoint_column(ion, 0)
        state = replace(
            state,
            mass_kg=np.full(two_N, complex_mass_amu(ANCHOR_N_START) * U),
            n_shell=genuine_n.copy(),
        )
        write_ion_state_to_checkpoint_column(state, ion, 2, mass_scenario="biphasic")
        np.testing.assert_array_equal(ion.n_shell[:, 2], genuine_n)


class TestBiphasicInitialState:
    """build_initial_ion_state biphasic branch: n0=21, S2 onset, E_pot binding fold."""

    @pytest.fixture(scope="class")
    def biphasic_init(self):
        from i2_helium_md.physics.internal_energy_budget import e_int_onset_eV
        from i2_helium_md.physics.solvation_cooling import e_bind_pair_eV

        cfg = single_pulse_N2000(num_molecules=5, seed=11)
        cfg = replace(
            cfg,
            t_max_neutral=0.1, dt_neutral=0.01, dt_ion=0.01,
            mass_scenario="biphasic",
            internal_energy_partition_fraction=0.3,
            internal_energy_retained_fraction=0.2,
            pickup_rate_coefficient=0.5,
            coulomb_available_eV=0.80,
        )
        neutral = run_neutral_propagation(cfg, verbose=False)
        ion = build_initial_ion_state(cfg, neutral, num_steps_ion=6)
        # A fixed build from the SAME neutral: identical positions, no binding fold,
        # so its E_pot is the pure MD conservative term at t0.
        ion_fixed = build_initial_ion_state(
            replace(cfg, mass_scenario="fixed"), neutral, num_steps_ion=6,
        )
        return cfg, ion, ion_fixed, e_int_onset_eV, e_bind_pair_eV

    def test_starts_at_full_first_shell(self, biphasic_init):
        _, ion, _, _, _ = biphasic_init
        expected_mass = complex_mass_amu(ANCHOR_N_START) * U
        np.testing.assert_allclose(ion.mass_kg, expected_mass)
        np.testing.assert_array_equal(ion.n_shell[:, 0], ANCHOR_N_START)

    def test_s2_onset_seeded(self, biphasic_init):
        cfg, ion, _, e_int_onset_eV, _ = biphasic_init
        onset = e_int_onset_eV(
            f_int=cfg.internal_energy_partition_fraction,
            e_avail_eV=cfg.coulomb_available_eV,
        )
        assert onset == pytest.approx(0.3 * 0.80)
        np.testing.assert_allclose(ion.E_int_eV[:, 0], onset)

    def test_e_pot_includes_binding_fold(self, biphasic_init):
        cfg, ion, ion_fixed, _, e_bind_pair_eV = biphasic_init
        e_bind = e_bind_pair_eV(
            ANCHOR_N_START,
            picture=cfg.ladder_electronic_picture,
            kappa=cfg.ladder_steepness,
        )
        assert e_bind < 0.0  # -Sigma(21) < 0
        # Positions are identical (same neutral); the fixed build carries only the
        # MD conservative E_pot, so biphasic E_pot = MD E_pot + e_bind_pair(21).
        np.testing.assert_allclose(
            ion.E_pot_eV[:, 0], ion_fixed.E_pot_eV[:, 0] + e_bind
        )


def _biphasic_state(n0=ANCHOR_N_START, e_int=0.24, speed=3.0, two_N=6):
    """A hand-built consistent biphasic IonStepState (m == m_I+ + n*m_He)."""
    n_arr = np.full(two_N, float(n0))
    m = complex_mass_amu(n0) * U
    return IonStepState(
        x=np.linspace(-2.0, 2.0, two_N),
        y=np.zeros(two_N),
        z=np.zeros(two_N),
        vx=np.full(two_N, speed),
        vy=np.zeros(two_N),
        vz=np.zeros(two_N),
        mass_kg=np.full(two_N, m),
        E_kin_eV=np.zeros(two_N),
        E_pot_eV=np.zeros(two_N),
        E_dissip_eV=np.zeros(two_N),
        E_mass_transfer_eV=np.zeros(two_N),
        E_int_eV=np.full(two_N, e_int),
        number_of_collisions=np.zeros(two_N, dtype=int),
        time_ps=0.0,
        n_shell=n_arr,
    )


def _biphasic_cfg(**overrides):
    overrides.setdefault("pickup_rate_coefficient", 0.0)  # pickup off by default here
    return _base_biphasic_cfg(**overrides)


class TestCoolingDrainClosure:
    """K2 cooling: E_int drains by exp(-dt/tau); the drain books to E_dissip (closes)."""

    def test_cooling_only_drain_is_exact(self):
        cfg = _biphasic_cfg(pickup_rate_coefficient=0.0)
        # E_int = 0.24 eV sits ABOVE Sigma(21) (mix ~0.17-0.19) -> self-unbound ->
        # evaporation suppressed; lambda0=0 -> no pickup. Only cooling acts.
        state = _biphasic_state(e_int=0.24)
        two_N = state.x.shape[0]
        droplet = np.full(two_N, 30.0)   # deep inside; irrelevant (no pickup)
        rng = np.random.default_rng(0)

        new = biphasic_step(
            state, rng=rng, cfg=cfg, droplet_radii=droplet, gate_steepness=14.2,
        )
        decay = np.exp(-cfg.dt_ion / cfg.internal_energy_cooling_tau_ps)
        drain = 0.24 * (1.0 - decay)
        np.testing.assert_allclose(new.E_int_eV, 0.24 * decay)
        np.testing.assert_allclose(new.E_dissip_eV, drain)
        # velocity / mass / n untouched (no event)
        np.testing.assert_array_equal(new.vx, state.vx)
        np.testing.assert_array_equal(new.mass_kg, state.mass_kg)
        np.testing.assert_array_equal(new.n_shell, state.n_shell)
        # per-channel closure: E_int down == E_dissip up
        np.testing.assert_allclose(
            (state.E_int_eV - new.E_int_eV), (new.E_dissip_eV - state.E_dissip_eV)
        )

    def test_e_int_stays_nonnegative_and_time_unchanged(self):
        cfg = _biphasic_cfg()
        state = _biphasic_state(e_int=0.24)
        two_N = state.x.shape[0]
        rng = np.random.default_rng(1)
        new = biphasic_step(
            state, rng=rng, cfg=cfg, droplet_radii=np.full(two_N, 30.0),
            gate_steepness=14.2,
        )
        assert np.all(new.E_int_eV >= 0.0)
        assert new.time_ps == state.time_ps  # pre-step seam does not advance time


class TestMassOccupancyInvariant:
    """m == m_I+ + n*m_He holds every step (the deferred Phase-B Q3 guard)."""

    def test_invariant_holds_over_many_steps_with_events(self):
        # Hot enough + pickup on so both channels fire across the run.
        cfg = _biphasic_cfg(pickup_rate_coefficient=5.0)
        state = _biphasic_state(n0=5, e_int=0.05, two_N=8)  # low E_int -> evap opens
        two_N = state.x.shape[0]
        droplet = np.full(two_N, 30.0)
        rng = np.random.default_rng(3)
        for _ in range(50):
            state = biphasic_step(
                state, rng=rng, cfg=cfg, droplet_radii=droplet, gate_steepness=14.2,
            )
            m_amu = state.mass_kg / U
            np.testing.assert_allclose(
                m_amu, MASS_I_ION_AMU + state.n_shell * MASS_HE_AMU, atol=1e-6
            )
            assert np.all(state.n_shell >= 0)


def _five_term_residual(state, new, cfg):
    """Per-ion 5-term residual of one biphasic_step (no BAOAB / conservative force).

    Sums dE_kin (mechanical) + the e_bind_pair(n) E_pot fold delta (the driver's
    booking, recomputed here from n_shell) + dE_dissip + dE_mass_transfer + dE_int.
    Zero to machine precision when the step's bookings close (MASS §6).
    """
    from i2_helium_md.physics.solvation_cooling import e_bind_pair_eV

    dE_kin = _E_kin_eV(new.mass_kg, new.vx, new.vy, new.vz) - _E_kin_eV(
        state.mass_kg, state.vx, state.vy, state.vz
    )
    e_bind = lambda ns: e_bind_pair_eV(  # noqa: E731
        ns, picture=cfg.ladder_electronic_picture, kappa=cfg.ladder_steepness
    )
    dE_pot_bind = np.asarray(e_bind(new.n_shell)) - np.asarray(e_bind(state.n_shell))
    dE_dis = new.E_dissip_eV - state.E_dissip_eV
    dE_mt = new.E_mass_transfer_eV - state.E_mass_transfer_eV
    dE_int = new.E_int_eV - state.E_int_eV
    return dE_kin + dE_pot_bind + dE_dis + dE_mt + dE_int


class TestIsolatedEventClosure:
    """Per-ion 5-term closure to machine precision on a single forced event.

    biphasic_step's whole effect (cooling + reset + K1/S1 + bath) must close when
    combined with the driver's e_bind_pair(n) E_pot fold and the mechanical
    E_kin<->E_mass_transfer pair -- independently of the (here-absent) BAOAB O-step
    and conservative force. So the summed per-ion residual is ~0.
    """

    def _residual(self, state, new, cfg):
        return _five_term_residual(state, new, cfg)

    def test_forced_shed_closes(self):
        from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum

        cfg = _biphasic_cfg(
            pickup_rate_coefficient=0.0,        # no pickup
            evap_rate_prefactor_per_ps=1.0e4,   # force P_shed ~ 1
            evap_rrk_dof=1.0,                   # s=1 -> bracket=1, k=nu in band
        )
        n0 = 4
        d0 = d0_of_n(n0, picture=cfg.ladder_electronic_picture, kappa=cfg.ladder_steepness)
        sig = ladder_cumsum(n0, picture=cfg.ladder_electronic_picture, kappa=cfg.ladder_steepness)
        e_int = 0.5 * (d0 + sig)             # safely inside the shedding band
        state = _biphasic_state(n0=n0, e_int=e_int, two_N=6)
        rng = np.random.default_rng(0)
        new = biphasic_step(
            state, rng=rng, cfg=cfg, droplet_radii=np.full(6, 30.0), gate_steepness=14.2,
        )
        assert np.all(new.n_shell == n0 - 1)          # every ion shed
        np.testing.assert_allclose(self._residual(state, new, cfg), 0.0, atol=1e-12)

    def test_forced_pickup_closes(self):
        cfg = _biphasic_cfg(
            pickup_rate_coefficient=1.0e6,      # force P_attach ~ 1
        )
        # E_int high -> evaporation self-suppressed (E_int > Sigma(n)); only pickup.
        state = _biphasic_state(n0=2, e_int=1.0, two_N=6)
        rng = np.random.default_rng(0)
        new = biphasic_step(
            state, rng=rng, cfg=cfg, droplet_radii=np.full(6, 30.0), gate_steepness=14.2,
        )
        assert np.all(new.n_shell == 3)               # every ion picked up
        np.testing.assert_allclose(self._residual(state, new, cfg), 0.0, atol=1e-12)


class TestReproducibilityAndEvents:
    def _run_steps(self, seed, n_steps=40):
        cfg = _biphasic_cfg(pickup_rate_coefficient=5.0)
        state = _biphasic_state(n0=5, e_int=0.05, two_N=8)
        rng = np.random.default_rng(seed)
        for _ in range(n_steps):
            state = biphasic_step(
                state, rng=rng, cfg=cfg, droplet_radii=np.full(8, 30.0), gate_steepness=14.2,
            )
        return state

    def test_same_seed_reproduces(self):
        a = self._run_steps(seed=123)
        b = self._run_steps(seed=123)
        np.testing.assert_array_equal(a.n_shell, b.n_shell)
        np.testing.assert_array_equal(a.vx, b.vx)
        np.testing.assert_array_equal(a.E_int_eV, b.E_int_eV)

    def test_different_seed_diverges_and_events_fire(self):
        a = self._run_steps(seed=1)
        b = self._run_steps(seed=2)
        # Events actually fired (n moved off the n0=5 start for at least one seed).
        assert not np.all(a.n_shell == 5) or not np.all(b.n_shell == 5)
        assert not np.array_equal(a.n_shell, b.n_shell)


def _tiny_neutral():
    """A tiny synthetic neutral end-state: two ion pairs ~2.6 A apart, at rest,
    deep inside a 30 A droplet (mass is overridden to n=21 by the biphasic
    builder, so its value here is irrelevant)."""
    from i2_helium_md.simulation.checkpoint import NeutralCheckpoint

    N, two_N = 2, 4
    x = np.array([1.3, -1.3, 1.3, -1.3])
    y = np.zeros(two_N)
    z = np.array([0.0, 0.0, 6.0, 6.0])
    col = lambda a: a.reshape(two_N, 1)  # noqa: E731
    zeros = np.zeros((two_N, 1))
    return NeutralCheckpoint(
        num_molecules=N, time_ps=np.zeros(1),
        positions_x=col(x), positions_y=col(y), positions_z=col(z),
        velocities_x=zeros.copy(), velocities_y=zeros.copy(), velocities_z=zeros.copy(),
        mass_kg=np.full(two_N, 202.95 * U), droplet_radii=np.full(two_N, 30.0),
        r0=np.zeros(N), E_kin_eV=zeros.copy(), E_pot_eV=zeros.copy(),
        E_initial_eV=np.zeros(N), E_dissip_eV=zeros.copy(),
        L_droplet_eV_ps=zeros.copy(),
    )


def _biphasic_driver_cfg(**overrides):
    """A driver-runnable biphasic cfg on the drag preset (bundle + §6.6 override)."""
    from i2_helium_md.presets import single_pulse_N2000_drag

    cfg = single_pulse_N2000_drag(
        num_molecules=2, ion_simulation_time=0.2, dt_ion=0.01, seed=0,
    )
    return replace(
        cfg,
        mass_scenario="biphasic",
        internal_energy_partition_fraction=0.3,
        internal_energy_retained_fraction=0.2,
        pickup_rate_coefficient=0.5,
        coulomb_available_eV=0.80,
        allow_inconsistent_mass_pairing=True,  # §6.6 mid-window; constant coeffs
        **overrides,
    )


class TestDriverSmokeAndClosure:
    """Full run_ion_propagation biphasic dispatch: completes, closes, invariants hold."""

    @pytest.fixture(scope="class")
    def biphasic_run(self):
        from i2_helium_md.simulation.ion import run_ion_propagation

        cfg = _biphasic_driver_cfg()
        ck = run_ion_propagation(cfg, _tiny_neutral())
        return cfg, ck

    def test_completes_v7_and_finite(self, biphasic_run):
        _, ck = biphasic_run
        assert ck.schema_version == 7
        for name in ("positions_x", "velocities_x", "E_kin_eV", "E_pot_eV",
                     "E_dissip_eV", "E_mass_transfer_eV", "E_int_eV", "n_shell"):
            assert np.all(np.isfinite(getattr(ck, name))), name

    def test_onset_seeded_and_e_int_nonnegative(self, biphasic_run):
        cfg, ck = biphasic_run
        np.testing.assert_allclose(ck.E_int_eV[:, 0], 0.3 * 0.80)
        assert np.all(ck.E_int_eV >= 0.0)

    def test_mass_occupancy_consistency_across_columns(self, biphasic_run):
        _, ck = biphasic_run
        m_amu = ck.mass_history_kg / U
        np.testing.assert_allclose(
            m_amu, MASS_I_ION_AMU + ck.n_shell * MASS_HE_AMU, atol=1e-6
        )

    def test_five_term_invariant_closes_to_verlet(self, biphasic_run):
        from i2_helium_md.postprocess.energy_balance import ion_ledger_closure

        _, ck = biphasic_run
        closure = ion_ledger_closure(ck)
        e0 = abs(closure.E_system_eV[0])
        # Cooling / events / drag are booked exactly; only the conservative Coulomb
        # explosion drifts at Verlet O(dt^2) level (cf. the Tier-0 smoke's 1e-2 bound).
        assert closure.max_abs_residual_eV / e0 < 1e-2, (
            f"5-term drift {closure.max_abs_residual_eV / e0 * 100:.4f}% too large"
        )


class TestDriverGuards:
    """Slice-G review fixes: the driver refuses the config combinations that would
    silently corrupt a biphasic run."""

    def test_biphasic_without_drag_bundle_refused(self):
        # Without this guard the run dispatches onto the hard-sphere collision
        # path: E_int frozen, the E_pot binding fold lost after column 0, and a
        # frozen n_shell stored while attachment grows the mass.
        from i2_helium_md.simulation.ion import run_ion_propagation

        cfg = _biphasic_driver_cfg(drag_coefficients=None)
        with pytest.raises(NotImplementedError, match="drag_coefficients"):
            run_ion_propagation(cfg, _tiny_neutral())

    def test_tabulated_density_profile_refused_lazily(self):
        # 'tabulated' is a valid config enum whose sourced data array is a
        # deferred rule-2 carry; biphasic_step must refuse it at point-of-use
        # instead of silently substituting the analytic erf gate.
        cfg = _biphasic_cfg(helium_density_profile="tabulated")
        state = _biphasic_state()
        with pytest.raises(NotImplementedError, match="helium_density_profile"):
            biphasic_step(
                state, rng=np.random.default_rng(0), cfg=cfg,
                droplet_radii=np.full(6, 30.0), gate_steepness=14.2,
            )

    def test_tabulated_dissociation_ladder_refused_lazily(self):
        # Review fix (2026-07-02): 'tabulated' is a valid config enum, but the
        # biphasic energetics compose the Form-U module functions directly (the
        # TabulatedLadder fallback is a module-level object with no config data
        # path); selecting it must refuse at point-of-use instead of silently
        # running the Form-U ladder -- the helium_density_profile contract.
        cfg = _biphasic_cfg(dissociation_ladder="tabulated")
        state = _biphasic_state()
        with pytest.raises(NotImplementedError, match="dissociation_ladder"):
            biphasic_step(
                state, rng=np.random.default_rng(0), cfg=cfg,
                droplet_radii=np.full(6, 30.0), gate_steepness=14.2,
            )

    def test_strided_biphasic_run_keeps_genuine_n_shell_columns(self):
        # Slice-G review fix companion (ion.py writer-call audit): a strided
        # biphasic run (max_bytes forces stride=2 here) must keep every stored
        # column m<->n consistent through the biphasic verbatim writer branch.
        # (The post-loop tail write itself is unreachable defensive code --
        # 1 + floor((n-1)/s) == ceil(n/s) identically -- but now passes the
        # required mass_scenario kwarg, so it stays correct if ever reached;
        # an omitted kwarg on any writer call path is a TypeError since the
        # kwarg became required.)
        from i2_helium_md.simulation.ion import run_ion_propagation

        cfg = _biphasic_driver_cfg()
        ck = run_ion_propagation(cfg, _tiny_neutral(), max_bytes=5_000)
        assert ck.n_shell.shape[1] < 21  # self-check: stride actually > 1
        m_amu = ck.mass_history_kg / U
        np.testing.assert_allclose(
            m_amu, MASS_I_ION_AMU + ck.n_shell * MASS_HE_AMU, atol=1e-6,
        )
        assert np.all(np.isfinite(ck.n_shell))


class TestScopeAcceptsBiphasic:
    def test_biphasic_passes_check_drag_scope(self):
        cfg = SimConfig(
            mass_scenario="biphasic",
            internal_energy_partition_fraction=0.3,
            internal_energy_retained_fraction=0.2,
            pickup_rate_coefficient=0.5,
        )
        mass = np.full(4, complex_mass_amu(ANCHOR_N_START) * U)
        _check_drag_scope(cfg, mass)  # must not raise (biphasic now accepted)

    def test_biphasic_with_noise_still_rejected(self):
        cfg = SimConfig(
            mass_scenario="biphasic",
            internal_energy_partition_fraction=0.3,
            internal_energy_retained_fraction=0.2,
            noise_form="fdt",
        )
        mass = np.full(4, complex_mass_amu(ANCHOR_N_START) * U)
        with pytest.raises(NotImplementedError, match="noise"):
            _check_drag_scope(cfg, mass)

    def test_unknown_scenario_still_rejected(self):
        # Restores the driver-level reject coverage removed with the old
        # test_scope_guard_rejects_mass_scenario: a typo'd scenario must fail
        # loudly at _check_drag_scope instead of falling through every dispatch
        # branch and silently running a fixed-mass drag trajectory.
        cfg = SimConfig(mass_scenario="biphsic")  # deliberate typo
        mass = np.full(4, complex_mass_amu(ANCHOR_N_START) * U)
        with pytest.raises(NotImplementedError, match="mass_scenario"):
            _check_drag_scope(cfg, mass)


# ===========================================================================
# Slice-G review extension (2026-07-02): plan-spec oracles not previously
# covered -- the both-fire one-event order, the frozen RNG draw-order
# composition lock, the gate opening at t_x, E_int >= 0 under events, and the
# fail-loud guards of biphasic_step (TIER2_PHASE_C_IMPLEMENTATION_PLAN.md §3
# Slice G "Oracle values / analytic limits").
# ===========================================================================
class TestBothFireShedThenPickup:
    """Both Bernoullis fire on one step => shed THEN pickup (frozen event order).

    nu = 1e4/ps and lambda0 = 1e6/ps make both per-step fire probabilities
    1 - exp(-k*dt) = 1 to within ~e-100, so the draw outcome is deterministic
    in practice. The frozen order is proven by the E_int bookings: shed at n0
    drains -D_0(n0) (K1), then pickup reads the POST-shed rung n_e = n0-1 and
    deposits +f_ret*D_0(n_e+1) = +f_ret*D_0(n0) (S1). The swapped order would
    book +f_ret*D_0(n0+1) - D_0(n0+1) instead -- numerically distinct rungs.
    """

    def test_both_fire_net_n_unchanged_and_books_close(self):
        from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum

        cfg = _biphasic_cfg(
            pickup_rate_coefficient=1.0e6,     # P_attach ~ 1 (certain fire)
            evap_rate_prefactor_per_ps=1.0e4,  # P_shed ~ 1 (certain fire)
            evap_rrk_dof=1.0,                  # s=1 -> k=nu inside the band
        )
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        n0 = 4
        d0 = d0_of_n(n0, picture=picture, kappa=kappa)
        sig = ladder_cumsum(n0, picture=picture, kappa=kappa)
        e_int0 = 0.5 * (d0 + sig)              # inside the shedding band
        state = _biphasic_state(n0=n0, e_int=e_int0, two_N=6)
        new = biphasic_step(
            state, rng=np.random.default_rng(0), cfg=cfg,
            droplet_radii=np.full(6, 30.0), gate_steepness=14.2,
        )

        # One shed AND one pickup applied -> net n unchanged, mass round-trips
        # (-m_He then +m_He; default assert_allclose rtol covers the float trip).
        np.testing.assert_array_equal(new.n_shell, n0)
        np.testing.assert_allclose(new.mass_kg, state.mass_kg)

        # E_int: cooling decay, then K1 -D_0(n0), then S1 +f_ret*D_0(n0).
        # atol 1e-12: the driver's E_int passes through the e_inf +/- cancellation
        # in newton_cool_step (noise ~ eps*|e_inf| ~ 1e-16 eV), nothing looser.
        decay = np.exp(-cfg.dt_ion / cfg.internal_energy_cooling_tau_ps)
        f_ret = cfg.internal_energy_retained_fraction
        np.testing.assert_allclose(
            new.E_int_eV, e_int0 * decay - d0 + f_ret * d0, atol=1e-12,
        )
        # E_dissip: the cooling drain plus the S1 bath remainder (1-f_ret)*D_0(n0).
        np.testing.assert_allclose(
            new.E_dissip_eV, e_int0 * (1.0 - decay) + (1.0 - f_ret) * d0, atol=1e-12,
        )
        # Full 5-term residual closes to machine precision (same bound as the
        # isolated-event closure tests).
        np.testing.assert_allclose(_five_term_residual(state, new, cfg), 0.0, atol=1e-12)


class TestFrozenDrawOrderComposition:
    """Lock the frozen evaporation-then-pickup RNG draw order (forbidden-list item).

    Replays biphasic_step manually: the same cooling arithmetic, then
    evaporation_step_components followed by pickup_step_components on a fresh
    rng with the same seed. Exact (bitwise) equality proves biphasic_step
    consumes the RNG stream in exactly this order and count -- one (2N,)
    uniform block per channel, evaporation first. A swapped order or an extra
    draw would de-synchronize the streams and fail the equality.
    """

    def test_manual_composition_reproduces_step(self):
        from i2_helium_md.physics.evaporation import evaporation_step_components
        from i2_helium_md.physics.helium_density import rho_he_ratio
        from i2_helium_md.physics.internal_energy_budget import pickup_bath_release_eV
        from i2_helium_md.physics.pickup import pickup_step_components
        from i2_helium_md.physics.solvation_cooling import e_infinity_eV, newton_cool_step

        # Mid-range probabilities (P_shed ~ 0.3, P_attach ~ 0.5) so which uniform
        # block feeds which channel actually decides the outcome per ion.
        cfg = _biphasic_cfg(
            pickup_rate_coefficient=100.0,
            evap_rate_prefactor_per_ps=36.0,
            evap_rrk_dof=1.0,
        )
        from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum

        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        n0 = 5
        d0 = d0_of_n(n0, picture=picture, kappa=kappa)
        sig = ladder_cumsum(n0, picture=picture, kappa=kappa)
        state = _biphasic_state(n0=n0, e_int=0.5 * (d0 + sig), two_N=8)
        droplet = np.full(8, 30.0)
        seed, steepness = 42, 14.2

        new = biphasic_step(
            state, rng=np.random.default_rng(seed), cfg=cfg,
            droplet_radii=droplet, gate_steepness=steepness,
        )

        # Manual replay on a fresh rng with the same seed, reusing the same
        # physics functions so the float path is identical.
        rng = np.random.default_rng(seed)
        n = np.asarray(state.n_shell, dtype=float)
        e_inf = np.asarray(e_infinity_eV(
            n, picture=picture, kappa=kappa, s_abs_eV=cfg.solv_struct_asymptote_eV,
        ))
        e_solv_cooled = np.asarray(newton_cool_step(
            e_inf + state.E_int_eV, n,
            tau_ps=cfg.internal_energy_cooling_tau_ps, dt_ps=cfg.dt_ion,
            picture=picture, kappa=kappa, s_abs_eV=cfg.solv_struct_asymptote_eV,
        ))
        E_int = e_solv_cooled - e_inf
        E_dissip = state.E_dissip_eV + (state.E_int_eV - E_int)
        depth = np.sqrt(state.x ** 2 + state.y ** 2 + state.z ** 2) - droplet
        rho = rho_he_ratio(depth, steepness=steepness)

        n_e, m_e, vx_e, vy_e, vz_e, dE_int_e, dE_mt_e, _ = evaporation_step_components(
            rng=rng, E_int_eV=E_int, n=n, vx=state.vx, vy=state.vy, vz=state.vz,
            m_amu=state.mass_kg / U, nu=cfg.evap_rate_prefactor_per_ps,
            picture=picture, kappa=kappa, dt_ps=cfg.dt_ion,
            evap_rrk_dof=cfg.evap_rrk_dof, gate_onset_eV=cfg.gate_onset_override_eV,
        )
        E_int = E_int + dE_int_e
        E_mt = state.E_mass_transfer_eV + _amu_ang2_ps2_to_eV(dE_mt_e)
        n_p, m_p, vx_p, vy_p, vz_p, dE_int_p, dE_mt_p, fired_p = pickup_step_components(
            rng=rng, n=n_e, vx=vx_e, vy=vy_e, vz=vz_e, m_amu=m_e, rho_ratio=rho,
            lambda0=cfg.pickup_rate_coefficient,
            f_ret=cfg.internal_energy_retained_fraction,
            picture=picture, kappa=kappa, dt_ps=cfg.dt_ion,
            p=cfg.pickup_occupancy_exponent, cap=cfg.pickup_occupancy_cap,
            pickup_rate_form=cfg.pickup_rate_form,
            he_capture_velocity=cfg.he_capture_velocity,
        )
        E_int = E_int + dE_int_p
        E_mt = E_mt + _amu_ang2_ps2_to_eV(dE_mt_p)
        bath = np.where(
            fired_p,
            pickup_bath_release_eV(
                n_e, f_ret=cfg.internal_energy_retained_fraction,
                picture=picture, kappa=kappa,
            ),
            0.0,
        )
        E_dissip = E_dissip + bath

        # The replay must be non-trivial: at least one ion shed and at least one
        # picked up under this seed (mixed outcomes, or the order lock is vacuous).
        assert np.any(np.asarray(n_e) < n0), "no evaporation fired; adjust seed/rates"
        assert np.any(fired_p), "no pickup fired; adjust seed/rates"

        np.testing.assert_array_equal(new.n_shell, np.asarray(n_p, dtype=float))
        np.testing.assert_array_equal(new.mass_kg, np.asarray(m_p, dtype=float) * U)
        np.testing.assert_array_equal(new.vx, np.asarray(vx_p, dtype=float))
        np.testing.assert_array_equal(new.vy, np.asarray(vy_p, dtype=float))
        np.testing.assert_array_equal(new.vz, np.asarray(vz_p, dtype=float))
        np.testing.assert_array_equal(new.E_int_eV, E_int)
        np.testing.assert_array_equal(new.E_dissip_eV, E_dissip)
        np.testing.assert_array_equal(new.E_mass_transfer_eV, E_mt)


class TestGateOpensAtCrossing:
    """The self-bound gate opens exactly as cooling drains E_int below Sigma(n).

    E_int(0) = Sigma(4)/decay^3.5 puts the crossing between cooling application
    3 and 4 *independently of tau*: after k steps E_int = E_int(0)*decay^k, so
    steps 1-3 sit above Sigma (self-unbound, k = 0 exactly -> no draw can fire)
    and step 4 lands at Sigma*decay^0.5 -- inside the band (> D_0(4)), where
    s=1 and nu=1e4/ps give P_shed = 1 - exp(-100), a certain fire to ~e-100.
    """

    def test_no_shed_above_sigma_then_fires_on_crossing(self):
        from i2_helium_md.physics.dissociation_ladder import ladder_cumsum

        cfg = _biphasic_cfg(
            pickup_rate_coefficient=0.0,       # evaporation-only
            evap_rate_prefactor_per_ps=1.0e4,
            evap_rrk_dof=1.0,
        )
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        n0 = 4
        sig = ladder_cumsum(n0, picture=picture, kappa=kappa)
        decay = np.exp(-cfg.dt_ion / cfg.internal_energy_cooling_tau_ps)
        state = _biphasic_state(n0=n0, e_int=sig / decay ** 3.5, two_N=6)
        rng = np.random.default_rng(0)
        droplet = np.full(6, 30.0)

        n_trace = []
        for _ in range(4):
            state = biphasic_step(
                state, rng=rng, cfg=cfg, droplet_radii=droplet, gate_steepness=14.2,
            )
            n_trace.append(state.n_shell.copy())

        # Steps 1-3: post-cooling E_int = sig*decay^(k-3.5) > Sigma(4) -> gate shut.
        for k in (0, 1, 2):
            np.testing.assert_array_equal(n_trace[k], n0)
        # Step 4: post-cooling E_int = sig*decay^0.5 < Sigma(4) -> every ion sheds.
        np.testing.assert_array_equal(n_trace[3], n0 - 1)


class TestEIntNonNegativeUnderEvents:
    def test_e_int_stays_nonnegative_across_event_steps(self):
        # Plan §3 Slice-G oracle "E_int finite + gate correct": with both channels
        # firing across 50 steps (same hot regime as the m<->n invariant test),
        # the reservoir never goes negative -- K1 can only fire while
        # E_int > D_0(n), so the drain cannot overshoot.
        cfg = _biphasic_cfg(pickup_rate_coefficient=5.0)
        state = _biphasic_state(n0=5, e_int=0.05, two_N=8)
        rng = np.random.default_rng(7)
        for _ in range(50):
            state = biphasic_step(
                state, rng=rng, cfg=cfg, droplet_radii=np.full(8, 30.0),
                gate_steepness=14.2,
            )
            assert np.all(state.E_int_eV >= 0.0)
            assert np.all(np.isfinite(state.E_int_eV))


class TestBiphasicStepGuards:
    def test_none_n_shell_raises(self):
        # Documented Raises contract: biphasic requires genuine occupancy state.
        cfg = _biphasic_cfg()
        state = replace(_biphasic_state(), n_shell=None)
        with pytest.raises(ValueError, match="n_shell"):
            biphasic_step(
                state, rng=np.random.default_rng(0), cfg=cfg,
                droplet_radii=np.full(6, 30.0), gate_steepness=14.2,
            )

    def test_m_n_drift_raises(self):
        # The deferred Phase-B Q3 guard has teeth: a pre-existing half-He
        # mass/n_shell inconsistency passes through the (event-free: lambda0=0,
        # E_int=1.0 eV > Sigma(5) so the gate is shut) channels untouched and
        # must trip the 1e-6 amu per-step assert.
        cfg = _biphasic_cfg()
        state = _biphasic_state(n0=5, e_int=1.0)
        state = replace(state, mass_kg=state.mass_kg + 0.5 * U)  # +0.5 amu drift
        with pytest.raises(AssertionError, match="consistency drift"):
            biphasic_step(
                state, rng=np.random.default_rng(0), cfg=cfg,
                droplet_radii=np.full(6, 30.0), gate_steepness=14.2,
            )
