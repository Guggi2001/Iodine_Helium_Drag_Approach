"""Tests for simulation/relaxation_stage.py (Tier-2 Phase E, Slice E2).

The post-ejection relaxation stage composes the delivered ``biphasic_step``
verbatim under a lambda_0 = 0 relaxation view + a zero-gamma BAOAB closure. The
correctness gate is the **5-term invariant** (as for the ion stage); the sharp
oracles exercise the decoupling fact (coulomb vs free_flight give the identical
shed sequence), the freeze termination, and the frozen two-draw RNG consumption.
Tiny synthetic checkpoints + one tiny real biphasic run; no figures.
"""

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.config import check_relaxation_config
from i2_helium_md.physics.constants import MASS_HE_AMU, MASS_I_ION_AMU, U
from i2_helium_md.physics.dissociation_ladder import d0_of_n, ladder_cumsum
from i2_helium_md.physics.evaporation import is_self_bound
from i2_helium_md.physics.solvation_cooling import e_bind_pair_eV
from i2_helium_md.postprocess.energy_balance import ion_ledger_closure
from i2_helium_md.postprocess.size_distribution import compute_terminal_shell_distribution
from i2_helium_md.simulation.checkpoint import IonCheckpoint, load_ion_checkpoint
from i2_helium_md.simulation.ion import run_ion_propagation
from i2_helium_md.simulation.relaxation_stage import (
    RELAXATION_STREAM_KEY,
    RelaxationResult,
    run_relaxation_stage,
)

from tests.test_biphasic_step import _biphasic_driver_cfg, _tiny_neutral


def _relax_cfg(**overrides):
    """A biphasic driver cfg with the relaxation stage enabled."""
    base = dict(
        relaxation_stage_enabled=True,
        relaxation_time_ps=1.0,
        relaxation_forces="coulomb",
    )
    base.update(overrides)
    return _biphasic_driver_cfg(**base)


def _seed_checkpoint(cfg, *, n_shell=21, E_int_eV=0.15, time_ps_start=0.0):
    """A tiny one-column v7 biphasic IonCheckpoint with m == m_I+ + n*m_He.

    Two molecules (2N=4), paired atoms (i, i+N) separated in z (mirrors
    ``_tiny_neutral`` so the per-pair Coulomb is finite). Positions are inside the
    30 A droplet. E_pot/E_kin are not self-consistent MD values (fine for the
    n_shell/E_int dynamics tests, which do not check the 5-term closure -- that
    uses the real-driver seed instead). ``time_ps_start`` sets the seed's absolute
    time (a production ion stage ends near ~20 ps; default 0 keeps the historic
    fixtures unchanged).
    """
    N = 2
    two_n = 2 * N
    T = 1
    m_amu = MASS_I_ION_AMU + n_shell * MASS_HE_AMU

    def col(val):
        return np.full((two_n, T), float(val))

    x = np.array([[1.3], [-1.3], [1.3], [-1.3]])
    y = np.zeros((two_n, T))
    z = np.array([[0.0], [0.0], [6.0], [6.0]])

    return IonCheckpoint(
        num_molecules=N,
        time_ps=np.full(T, float(time_ps_start)),
        positions_x=x, positions_y=y, positions_z=z,
        velocities_x=col(0.2), velocities_y=col(0.0), velocities_z=col(0.0),
        positions_final_x=np.zeros(two_n), positions_final_y=np.zeros(two_n),
        positions_final_z=np.zeros(two_n),
        velocities_final_x=np.zeros(two_n), velocities_final_y=np.zeros(two_n),
        velocities_final_z=np.zeros(two_n),
        mass_kg=np.full(two_n, m_amu * U), mass_final_kg=np.full(two_n, m_amu * U),
        mass_history_kg=np.full((two_n, T), m_amu * U),
        droplet_radii_angstrom=np.full(two_n, 30.0),
        E_kin_eV=col(0.01), E_pot_eV=col(0.0), E_dissip_eV=np.zeros((two_n, T)),
        E_mass_transfer_eV=np.zeros((two_n, T)), E_int_eV=col(E_int_eV),
        n_shell=col(n_shell),
        b_ion_outside=np.zeros(N, dtype=bool),
        relative_loss_per_ps=np.zeros((two_n, T)),
        number_of_collisions=np.zeros((two_n, T), dtype=int),
        temperature_diagnostic=np.full((T, 3), np.nan),
        mass_scenario="biphasic",
    )


# ---------------------------------------------------------------------------
# 5-term invariant (the correctness gate) -- on the real-driver seed
# ---------------------------------------------------------------------------
class TestInvariantClosure:
    @pytest.fixture(scope="class")
    def real_seed(self):
        cfg = _biphasic_driver_cfg()
        ion = run_ion_propagation(cfg, _tiny_neutral())
        return cfg, ion

    def _closes(self, ckpt, tol=1e-4):
        # Tolerance justification (review 2026-07-03): the measured closure
        # residual on this fixture is ~2e-5 relative (FP accumulation over the
        # ~300-step window); 1e-4 gives ~5x headroom while still catching a
        # single mis-booked D_0 rung (~9 meV ~ 4e-3 of E_system ~ 2.5 eV). The
        # earlier 1e-2 band (~25 meV absolute) was loose enough to hide
        # per-rung fold errors.
        closure = ion_ledger_closure(ckpt)
        e0 = abs(closure.E_system_eV[0])
        return closure.max_abs_residual_eV / e0 < tol

    def _relaxed(self, cfg, ion, forces):
        # evap_rrk_dof=2 + fast cooling so the window actually sheds -> the
        # invariant covers the K1 drain / e_bind fold / cold-shed bookings, not
        # only K2 cooling.
        cfg = replace(cfg, relaxation_stage_enabled=True, relaxation_time_ps=3.0,
                      relaxation_forces=forces, evap_rrk_dof=2.0,
                      internal_energy_cooling_tau_ps=1.0)
        return run_relaxation_stage(ion, cfg)

    def test_coulomb_arm_closes(self, real_seed):
        cfg, ion = real_seed
        result = self._relaxed(cfg, ion, "coulomb")
        assert np.any(result.terminal_n < result.checkpoint.n_shell[:, 0])  # sheds fired
        assert self._closes(result.checkpoint)

    def test_free_flight_arm_closes(self, real_seed):
        cfg, ion = real_seed
        result = self._relaxed(cfg, ion, "free_flight")
        assert np.any(result.terminal_n < result.checkpoint.n_shell[:, 0])  # sheds fired
        assert self._closes(result.checkpoint)


# ---------------------------------------------------------------------------
# Decoupling fact: coulomb and free_flight give the identical shed sequence
# ---------------------------------------------------------------------------
class TestDecoupling:
    def test_identical_shed_sequence(self):
        # evap_rrk_dof=2 -> a moderate RRK rate so sheds actually fire (the native
        # n=21 dof s=60 makes the per-step rate vanishingly small); the decoupling
        # claim is only meaningful when the shed sequence is non-trivial.
        cfg = _relax_cfg(relaxation_time_ps=3.0, internal_energy_cooling_tau_ps=1.0,
                         evap_rrk_dof=2.0)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.98 * sigma21)

        r_coul = run_relaxation_stage(ion, replace(cfg, relaxation_forces="coulomb"))
        r_ff = run_relaxation_stage(ion, replace(cfg, relaxation_forces="free_flight"))

        # The shed sequence is non-trivial (the test would be vacuous otherwise).
        assert np.any(r_coul.terminal_n < 21)
        # n and E_int trajectories are identical; only x/v/E_kin/E_pot differ.
        np.testing.assert_array_equal(
            r_coul.checkpoint.n_shell, r_ff.checkpoint.n_shell
        )
        np.testing.assert_allclose(
            r_coul.checkpoint.E_int_eV, r_ff.checkpoint.E_int_eV, rtol=0, atol=0
        )
        np.testing.assert_array_equal(r_coul.terminal_n, r_ff.terminal_n)


# ---------------------------------------------------------------------------
# Shed dynamics: monotone, no avalanche, terminal n <= seed n
# ---------------------------------------------------------------------------
class TestShedDynamics:
    def test_hot_input_sheds_monotonically(self):
        cfg = _relax_cfg(relaxation_time_ps=5.0, internal_energy_cooling_tau_ps=1.0,
                         evap_rrk_dof=2.0)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.98 * sigma21)

        result = run_relaxation_stage(ion, cfg)
        n = result.checkpoint.n_shell

        # Pickup off -> n is monotone non-increasing over the window, per ion.
        assert np.all(np.diff(n, axis=1) <= 0)
        # No avalanche: biphasic_step sheds at most one rung per step (the plan's
        # per-step bound; a multi-rung loop bug would still be monotone).
        assert np.all(np.diff(n, axis=1) >= -1.0)
        # Terminal n never exceeds the seed n* = 21 and stays >= 0.
        assert np.all(result.terminal_n <= 21)
        assert np.all(result.terminal_n >= 0)
        # A hot input actually sheds (some ion drops below 21).
        assert np.any(result.terminal_n < 21)
        # terminal_n is the last checkpoint column.
        np.testing.assert_array_equal(result.terminal_n, n[:, -1])

    def test_freeze_termination_reached(self):
        cfg = _relax_cfg(relaxation_time_ps=50.0, internal_energy_cooling_tau_ps=0.5,
                         evap_rrk_dof=2.0)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.98 * sigma21)

        result = run_relaxation_stage(ion, cfg)
        # Fast cooling freezes the whole ensemble well before 50 ps.
        assert bool(result.freeze_flags.all())
        assert result.time_relaxed_ps < 50.0


# ---------------------------------------------------------------------------
# Per-channel sharp oracles (review 2026-07-03): E_dissip == cumulative K2
# drain; free_flight E_pot == held MD + e_bind fold. The 5-term closure alone
# cannot catch channel *mis-attribution* (the sum is conserved), so these pin
# the individual bookings.
# ---------------------------------------------------------------------------
class TestEnergyChannels:
    def _hot_run(self, forces):
        cfg = _relax_cfg(relaxation_time_ps=3.0, internal_energy_cooling_tau_ps=1.0,
                         evap_rrk_dof=2.0, relaxation_forces=forces)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.98 * sigma21)
        return cfg, run_relaxation_stage(ion, cfg)

    def test_e_dissip_equals_cumulative_k2_drain(self):
        # With gamma = 0 and lambda_0 = 0 the ONLY E_dissip writer is the K2
        # cooling drain (plan E2 acceptance oracle). Reconstruct the drain from
        # the stored E_int / n_shell columns: per step, K2 acts first
        # (E_pre -> E_afterK2), then K1 subtracts D_0(n_pre) on a fire, so
        # E_afterK2 = E_int[t] + shed_t * D_0(n_pre) and
        # drain_t = E_pre - E_afterK2.
        cfg, result = self._hot_run("coulomb")
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        ckpt = result.checkpoint
        # The oracle needs every step stored (stride 1) -- assert the assumption.
        np.testing.assert_allclose(np.diff(ckpt.time_ps), cfg.dt_ion, rtol=0, atol=1e-12)

        n = ckpt.n_shell
        E_int = ckpt.E_int_eV
        shed = n[:, :-1] - n[:, 1:]                     # 0 or 1 per step
        d0_pre = np.asarray(d0_of_n(np.maximum(n[:, :-1], 1.0),
                                    picture=picture, kappa=kappa))
        e_after_k2 = E_int[:, 1:] + shed * d0_pre
        drain = E_int[:, :-1] - e_after_k2
        cum_drain = np.cumsum(drain, axis=1)

        dE_dissip = ckpt.E_dissip_eV[:, 1:] - ckpt.E_dissip_eV[:, [0]]
        # Non-vacuous: real drain accumulated and sheds fired.
        assert np.all(cum_drain[:, -1] > 0.0)
        assert np.any(shed > 0)
        # Tolerance: the reconstruction differs from the booked floats only by
        # re-association round-off (<= ~1e-14 eV over the window); 1e-9 eV is
        # ~6 orders below the smallest physical booking (per-step K2 drain ~meV,
        # D_0 ~ 9 meV), so any mis-attributed channel fails loudly.
        np.testing.assert_allclose(dE_dissip, cum_drain, rtol=0, atol=1e-9)

    def test_free_flight_e_pot_is_held_plus_fold(self):
        # Under free flight no conservative work is done: E_pot(t) must equal
        # the held MD potential (seed E_pot minus the seed fold) plus the
        # e_bind_pair fold at the current n -- i.e. a shed moves E_pot by
        # exactly +D_0(n). A constant-offset bug (e.g. forgetting to remove the
        # seed fold) cancels in the closure residual, so it is pinned here.
        cfg, result = self._hot_run("free_flight")
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        ckpt = result.checkpoint

        held_md = ckpt.E_pot_eV[:, 0] - e_bind_pair_eV(
            ckpt.n_shell[:, 0], picture=picture, kappa=kappa,
        )
        expected = held_md[:, None] + e_bind_pair_eV(
            ckpt.n_shell, picture=picture, kappa=kappa,
        )
        assert np.any(result.terminal_n < 21)           # sheds fired (non-vacuous)
        # Exact arithmetic chain (same fold call, same floats): tight tolerance.
        np.testing.assert_allclose(ckpt.E_pot_eV, expected, rtol=0, atol=1e-12)

    def _hot_run_gated(self, droplet_radius, *, gate="density_scaled"):
        # The _hot_run seed with the cooling_spatial_gate arm selectable and the ions
        # placed near a small droplet so rho_ratio < 1 actively attenuates the K2
        # drain during relaxation (the 30 A default keeps rho == 1 -> a vacuous twin).
        cfg = _relax_cfg(relaxation_time_ps=3.0, internal_energy_cooling_tau_ps=1.0,
                         evap_rrk_dof=2.0, relaxation_forces="coulomb",
                         cooling_spatial_gate=gate)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.98 * sigma21)
        ion = replace(
            ion,
            droplet_radii_angstrom=np.full(
                ion.droplet_radii_angstrom.shape, float(droplet_radius)
            ),
        )
        return cfg, run_relaxation_stage(ion, cfg)

    def test_e_dissip_oracle_holds_under_density_scaled_gate(self):
        # The E_dissip == cumulative-K2-drain reconstruction is gate-invariant: it
        # reads the ACTUAL E_int/n_shell trajectory, so a spatially-attenuated drain
        # still closes. Ions straddle a small (3 A) droplet -> rho_ratio in (0, 1).
        cfg, result = self._hot_run_gated(3.0)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        ckpt = result.checkpoint
        np.testing.assert_allclose(np.diff(ckpt.time_ps), cfg.dt_ion, rtol=0, atol=1e-12)

        n = ckpt.n_shell
        E_int = ckpt.E_int_eV
        shed = n[:, :-1] - n[:, 1:]
        d0_pre = np.asarray(d0_of_n(np.maximum(n[:, :-1], 1.0),
                                    picture=picture, kappa=kappa))
        e_after_k2 = E_int[:, 1:] + shed * d0_pre
        drain = E_int[:, :-1] - e_after_k2
        cum_drain = np.cumsum(drain, axis=1)
        dE_dissip = ckpt.E_dissip_eV[:, 1:] - ckpt.E_dissip_eV[:, [0]]
        assert np.all(cum_drain[:, -1] > 0.0)           # non-vacuous drain accrued
        np.testing.assert_allclose(dE_dissip, cum_drain, rtol=0, atol=1e-9)

    def test_density_scaled_first_step_drain_below_ungated(self):
        # First step, identical seed for both arms -> no trajectory divergence yet:
        # the gated per-ion K2 drain (rho_ratio < 1) is strictly below the ungated
        # `none` drain at the same small droplet. This proves the gate is active.
        _, gated = self._hot_run_gated(3.0, gate="density_scaled")
        _, none = self._hot_run_gated(3.0, gate="none")
        gated_step1 = gated.checkpoint.E_dissip_eV[:, 1]     # drain over the 1st step
        none_step1 = none.checkpoint.E_dissip_eV[:, 1]
        assert np.all(gated_step1 > 0.0)
        assert np.all(gated_step1 < none_step1)


# ---------------------------------------------------------------------------
# Gate + ladder rungs (review 2026-07-03): the self-bound gate suppresses, a
# gate-suppressed ion is NOT frozen (it un-freezes into the RRK band as K2
# cools it), and the n = 1 direct-dissociation rung reaches the bare ion.
# ---------------------------------------------------------------------------
class TestGateAndRungs:
    def test_self_bound_seed_gate_suppresses_then_sheds(self):
        cfg = _relax_cfg(relaxation_time_ps=10.0, internal_energy_cooling_tau_ps=0.5,
                         evap_rrk_dof=2.0)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        # E_int > Sigma(21): net self-unbound -> shedding suppressed at entry.
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=1.5 * sigma21)

        result = run_relaxation_stage(ion, cfg)
        ckpt = result.checkpoint
        np.testing.assert_allclose(np.diff(ckpt.time_ps), cfg.dt_ion, rtol=0, atol=1e-12)

        n = ckpt.n_shell
        E_int = ckpt.E_int_eV
        shed = n[:, :-1] - n[:, 1:]
        d0_pre = np.asarray(d0_of_n(np.maximum(n[:, :-1], 1.0),
                                    picture=picture, kappa=kappa))
        # The gate saw the post-K2 E_int (reconstructed exactly as in the
        # E_dissip oracle above).
        e_after_k2 = E_int[:, 1:] + shed * d0_pre
        n_pre = n[:, :-1]
        sb = np.asarray(is_self_bound(
            e_after_k2, n_pre, picture=picture, kappa=kappa,
            gate_onset_eV=cfg.gate_onset_override_eV,
        ))
        rung_ge2 = n_pre >= 2
        # Gate suppression: no n>=2 shed ever fires while net self-unbound.
        assert np.all(shed[rung_ge2 & ~sb] == 0)
        # Non-vacuous both ways: suppressed steps existed at entry, and the ion
        # was NOT frozen there (a freeze mask conflating the self-bound gate
        # with freeze would terminate at step 1 with terminal_n == 21).
        assert np.any(rung_ge2 & ~sb)
        assert np.any(result.terminal_n < 21)

    def test_n1_direct_dissociation_reaches_bare_ion(self):
        cfg = _relax_cfg(relaxation_time_ps=5.0,
                         internal_energy_cooling_tau_ps=1e6)  # cooling ~off: stay hot
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        d0_1 = float(d0_of_n(1, picture=picture, kappa=kappa))
        # n = 1 with E_int > D_0(1): the direct-dissociation rung, k = nu.
        ion = _seed_checkpoint(cfg, n_shell=1, E_int_eV=5.0 * d0_1)

        result = run_relaxation_stage(ion, cfg)
        n = result.checkpoint.n_shell

        # Every ion sheds its last He (P_miss ~ (1-nu*dt)^500 ~ 5e-6 per ion at
        # the seeded stream -- deterministic under seed) and freezes at n = 0.
        np.testing.assert_array_equal(result.terminal_n, np.zeros(4))
        assert bool(result.freeze_flags.all())
        assert result.time_relaxed_ps < 5.0             # freeze early-exit fired
        assert np.all(np.diff(n, axis=1) >= -1.0)       # one rung per step


# ---------------------------------------------------------------------------
# Cold input: zero sheds while both draws are consumed (RNG parity)
# ---------------------------------------------------------------------------
class TestColdInputAndRNG:
    def test_cold_input_zero_sheds_and_draw_parity(self):
        cfg = _relax_cfg(relaxation_time_ps=0.05)  # up to 5 steps; freezes at step 1
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        d0_21 = float(d0_of_n(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.5 * d0_21)  # < D_0 -> frozen

        rng = np.random.default_rng(123)
        result = run_relaxation_stage(ion, cfg, rng=rng)

        # Cold self-bound input sheds nothing.
        np.testing.assert_array_equal(result.terminal_n, np.full(4, 21.0))
        assert bool(result.freeze_flags.all())

        # Exactly one relaxation step ran (freeze after step 1), consuming the
        # frozen two-draw stream: evaporation random(2N) then the inert pickup
        # random(2N).
        expected = np.random.default_rng(123)
        expected.random(size=4)
        expected.random(size=4)
        assert rng.bit_generator.state == expected.bit_generator.state


# ---------------------------------------------------------------------------
# Artifact: a bona fide v7 checkpoint that round-trips
# ---------------------------------------------------------------------------
class TestArtifact:
    def test_checkpoint_round_trips(self, tmp_path):
        cfg = _relax_cfg(relaxation_time_ps=1.0, internal_energy_cooling_tau_ps=1.0)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.98 * sigma21)

        path = tmp_path / "relaxation.npz"
        result = run_relaxation_stage(ion, cfg, save_path=path)
        assert path.exists()

        loaded = load_ion_checkpoint(path)
        assert loaded.schema_version == 7
        assert loaded.mass_scenario == "biphasic"
        np.testing.assert_array_equal(loaded.n_shell, result.checkpoint.n_shell)
        # The 5-term ledger machinery applies to the artifact unchanged.
        assert np.isfinite(ion_ledger_closure(loaded).max_abs_residual_eV)


# ---------------------------------------------------------------------------
# E1 relaxed-input admission (the RelaxationResult -> ShellDistribution hook)
# ---------------------------------------------------------------------------
class TestE1Admission:
    def test_e1_accepts_relaxation_result(self):
        cfg = _relax_cfg(relaxation_time_ps=3.0, internal_energy_cooling_tau_ps=1.0)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.98 * sigma21)

        result = run_relaxation_stage(ion, cfg)
        assert isinstance(result, RelaxationResult)

        dist = compute_terminal_shell_distribution(result)
        assert dist.source == "relaxed"
        assert int(dist.counts.sum()) == 4  # 2N atoms
        # Histogram matches the terminal_n by hand.
        for n_val in np.rint(result.terminal_n).astype(int):
            assert dist.counts[n_val] >= 1


# ---------------------------------------------------------------------------
# dt honored + default-stream derivation (review 2026-07-03)
# ---------------------------------------------------------------------------
class TestDtAndSeeding:
    def test_relaxation_dt_honored(self):
        # A gate-suppressed, cooling-off seed: never sheds, never freezes
        # (E_int > D_0), so the stage runs the full window -- this also covers
        # the time-exhaustion (non-freeze) termination arm.
        cfg = _relax_cfg(relaxation_time_ps=0.1, relaxation_dt_ps=0.02,
                         internal_energy_cooling_tau_ps=1e6)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=1.2 * sigma21)

        result = run_relaxation_stage(ion, cfg)
        ckpt = result.checkpoint

        # relaxation_dt_ps (not dt_ion = 0.01) drives the step: 5 steps of 0.02.
        assert ckpt.time_ps.shape == (6,)
        np.testing.assert_allclose(np.diff(ckpt.time_ps), 0.02, rtol=0, atol=1e-12)
        assert not bool(result.freeze_flags.any())
        assert result.time_relaxed_ps == pytest.approx(0.1, abs=1e-12)
        np.testing.assert_array_equal(result.terminal_n, np.full(4, 21.0))

    def test_default_rng_derivation_pinned(self):
        # The default (rng=None) stream MUST be
        # default_rng(SeedSequence((cfg.seed, RELAXATION_STREAM_KEY))): seeding
        # from cfg.seed directly would replay the ion-stage stream (the exact
        # correlation the key exists to prevent) and no other test executes the
        # default path.
        cfg = _relax_cfg(relaxation_time_ps=3.0, internal_energy_cooling_tau_ps=1.0,
                         evap_rrk_dof=2.0, seed=777)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        sigma21 = float(ladder_cumsum(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.98 * sigma21)

        r_default = run_relaxation_stage(ion, cfg)
        r_explicit = run_relaxation_stage(
            ion, cfg,
            rng=np.random.default_rng(
                np.random.SeedSequence((777, RELAXATION_STREAM_KEY))
            ),
        )
        assert np.any(r_default.terminal_n < 21)        # draws mattered (non-vacuous)
        np.testing.assert_array_equal(
            r_default.checkpoint.n_shell, r_explicit.checkpoint.n_shell
        )
        np.testing.assert_allclose(
            r_default.checkpoint.E_int_eV, r_explicit.checkpoint.E_int_eV,
            rtol=0, atol=0,
        )


# ---------------------------------------------------------------------------
# Input contracts (review 2026-07-03): window-relative time_relaxed_ps; the
# seed checkpoint must be biphasic and its last column must be the true final
# ion state.
# ---------------------------------------------------------------------------
class TestInputContracts:
    def test_time_relaxed_is_window_relative(self):
        # A production seed enters at ~20 ps absolute; time_relaxed_ps is the
        # window-relative duration (here: one dt_ion step before the cold
        # freeze), while the checkpoint's time_ps axis stays absolute.
        cfg = _relax_cfg(relaxation_time_ps=0.05)
        picture, kappa = cfg.ladder_electronic_picture, cfg.ladder_steepness
        d0_21 = float(d0_of_n(21, picture=picture, kappa=kappa))
        ion = _seed_checkpoint(cfg, n_shell=21, E_int_eV=0.5 * d0_21,
                               time_ps_start=20.0)

        result = run_relaxation_stage(ion, cfg)
        assert result.time_relaxed_ps == pytest.approx(cfg.dt_ion, abs=1e-12)
        # Early-freeze detection works on a real (nonzero-time) seed.
        assert result.time_relaxed_ps < cfg.relaxation_time_ps
        # The checkpoint axis continues the ion stage's absolute time.
        assert result.checkpoint.time_ps[0] == pytest.approx(20.0, abs=1e-12)
        assert result.checkpoint.time_ps[-1] == pytest.approx(20.0 + cfg.dt_ion,
                                                              abs=1e-12)

    def test_non_biphasic_seed_checkpoint_rejected(self):
        cfg = _relax_cfg()
        ion = replace(_seed_checkpoint(cfg), mass_scenario="fixed")
        with pytest.raises(ValueError, match="mass_scenario"):
            run_relaxation_stage(ion, cfg)

    def test_stale_seed_column_rejected(self):
        # mass_history_kg[:, -1] != mass_final_kg means the ion run's stride
        # dropped the true final state (E_int/n_shell unrecoverable) -- the
        # stage must refuse rather than silently relax a stale seed.
        cfg = _relax_cfg()
        ion = _seed_checkpoint(cfg)
        ion = replace(ion, mass_final_kg=ion.mass_final_kg + MASS_HE_AMU * U)
        with pytest.raises(ValueError, match="final"):
            run_relaxation_stage(ion, cfg)


# ---------------------------------------------------------------------------
# Config validation (check_relaxation_config arms)
# ---------------------------------------------------------------------------
class TestConfigValidation:
    def test_disabled_is_noop(self):
        cfg = _biphasic_driver_cfg()  # relaxation disabled by default
        check_relaxation_config(cfg)  # no raise

    def test_non_biphasic_raises(self):
        from i2_helium_md.presets import single_pulse_N2000_drag

        cfg = single_pulse_N2000_drag(
            num_molecules=2, ion_simulation_time=0.2, dt_ion=0.01, seed=0,
        )
        cfg = replace(
            cfg, mass_scenario="fixed", relaxation_stage_enabled=True,
            relaxation_time_ps=1.0,
        )
        with pytest.raises(ValueError, match="biphasic"):
            check_relaxation_config(cfg)

    def test_missing_time_raises(self):
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=None,
        )
        with pytest.raises(ValueError, match="relaxation_time_ps"):
            check_relaxation_config(cfg)

    @pytest.mark.parametrize("bad_time", [0.0, -1.0])
    def test_zero_or_negative_time_raises(self, bad_time):
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=bad_time,
        )
        with pytest.raises(ValueError, match="relaxation_time_ps"):
            check_relaxation_config(cfg)

    @pytest.mark.parametrize("bad_dt", [0.0, -0.01])
    def test_nonpositive_dt_raises(self, bad_dt):
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=1.0,
            relaxation_dt_ps=bad_dt,
        )
        with pytest.raises(ValueError, match="relaxation_dt_ps"):
            check_relaxation_config(cfg)

    def test_nu_dt_guard_raises(self):
        # nu = 2.42; dt_relax = 0.05 -> nu*dt = 0.121 > 0.1.
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=1.0,
            relaxation_dt_ps=0.05,
        )
        with pytest.raises(ValueError, match="nu"):
            check_relaxation_config(cfg)

    def test_nu_dt_boundary_accepted(self):
        # The guard is nu*dt <= 0.1 inclusive: exactly 0.1 must pass (1.0 * 0.1
        # is exact in binary floats, so this pins the > vs >= arm).
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=1.0,
            relaxation_dt_ps=0.1, evap_rate_prefactor_per_ps=1.0,
        )
        check_relaxation_config(cfg)  # no raise

    def test_bad_forces_raises(self):
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=1.0,
            relaxation_forces="verlet",
        )
        with pytest.raises(ValueError, match="relaxation_forces"):
            check_relaxation_config(cfg)

    def test_density_scaled_rejects_runaway_relaxation_cap(self):
        # density_scaled voids the freeze early-exit, so a runaway cap must fail
        # loudly at config-load: 2e5 ps / 0.01 ps = 2e7 steps > the 1e7 ceiling.
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=2.0e5,
            relaxation_dt_ps=0.01, cooling_spatial_gate="density_scaled",
        )
        with pytest.raises(ValueError, match="density_scaled"):
            check_relaxation_config(cfg)

    def test_density_scaled_modest_cap_accepted(self):
        # The total-strip probe's 1000 ps / 0.01 = 1e5 steps is well under the
        # ceiling, so the guarded arm's own intended usage stays valid.
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=1000.0,
            relaxation_dt_ps=0.01, cooling_spatial_gate="density_scaled",
        )
        check_relaxation_config(cfg)  # no raise

    def test_ungated_cooling_allows_generous_cap(self):
        # 'none' keeps the guaranteed freeze early-exit, so the same generous cap
        # that trips the density_scaled guard stays allowed (guard is arm-specific).
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=2.0e5,
            relaxation_dt_ps=0.01, cooling_spatial_gate="none",
        )
        check_relaxation_config(cfg)  # no raise

    def test_run_requires_enabled(self):
        cfg = _biphasic_driver_cfg()  # disabled
        ion = _seed_checkpoint(cfg)
        with pytest.raises(ValueError, match="relaxation_stage_enabled"):
            run_relaxation_stage(ion, cfg)

    def test_stream_key_is_stable(self):
        # A guard against an accidental key change (would silently re-seed runs).
        assert RELAXATION_STREAM_KEY == 0xE2_2026
