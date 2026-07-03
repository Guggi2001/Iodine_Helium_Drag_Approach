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


def _seed_checkpoint(cfg, *, n_shell=21, E_int_eV=0.15):
    """A tiny one-column v7 biphasic IonCheckpoint with m == m_I+ + n*m_He.

    Two molecules (2N=4), paired atoms (i, i+N) separated in z (mirrors
    ``_tiny_neutral`` so the per-pair Coulomb is finite). Positions are inside the
    30 A droplet. E_pot/E_kin are not self-consistent MD values (fine for the
    n_shell/E_int dynamics tests, which do not check the 5-term closure -- that
    uses the real-driver seed instead).
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
        time_ps=np.zeros(T),
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

    def _closes(self, ckpt, tol=1e-2):
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

    def test_nu_dt_guard_raises(self):
        # nu = 2.42; dt_relax = 0.05 -> nu*dt = 0.121 > 0.1.
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=1.0,
            relaxation_dt_ps=0.05,
        )
        with pytest.raises(ValueError, match="nu"):
            check_relaxation_config(cfg)

    def test_bad_forces_raises(self):
        cfg = _biphasic_driver_cfg(
            relaxation_stage_enabled=True, relaxation_time_ps=1.0,
            relaxation_forces="verlet",
        )
        with pytest.raises(ValueError, match="relaxation_forces"):
            check_relaxation_config(cfg)

    def test_run_requires_enabled(self):
        cfg = _biphasic_driver_cfg()  # disabled
        ion = _seed_checkpoint(cfg)
        with pytest.raises(ValueError, match="relaxation_stage_enabled"):
            run_relaxation_stage(ion, cfg)

    def test_stream_key_is_stable(self):
        # A guard against an accidental key change (would silently re-seed runs).
        assert RELAXATION_STREAM_KEY == 0xE2_2026
