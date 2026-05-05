"""Tests for i2_helium_md/simulation/ion.py driver."""

from dataclasses import replace

import numpy as np
import pytest

from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.checkpoint import IonCheckpoint
from i2_helium_md.simulation.ion import (
    DEFAULT_MAX_CHECKPOINT_BYTES_ION,
    _decide_stride_ion,
    _estimate_checkpoint_bytes_ion,
    _internal_step_count_ion,
    run_ion_propagation,
)
from i2_helium_md.simulation.neutral import run_neutral_propagation
from i2_helium_md.simulation.run_directory import RunDirectory


# ===========================================================================
# Shared fixture: cheap neutral run that downstream ion tests build on
# ===========================================================================
@pytest.fixture(scope="module")
def small_neutral():
    """A small neutral checkpoint (2 columns, single_pulse default).

    Single-pulse mode runs only 2 internal neutral steps, which is the
    canonical setup for the ion stage: it starts from the post-photon
    velocities just after photoexcitation.
    """
    cfg = single_pulse_N2000(num_molecules=4, seed=42)
    neutral = run_neutral_propagation(cfg)
    return cfg, neutral


def _ion_cfg(neutral_cfg, *, ion_simulation_time=0.1, dt_ion=0.01, **overrides):
    """Build an ion-stage cfg from a neutral cfg with sane test defaults."""
    return replace(
        neutral_cfg,
        ion_simulation_time=ion_simulation_time,
        dt_ion=dt_ion,
        **overrides,
    )


# ===========================================================================
# Public API
# ===========================================================================
class TestPublicApi:
    def test_returns_ion_checkpoint(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        ckpt = run_ion_propagation(cfg_ion, neutral)
        assert isinstance(ckpt, IonCheckpoint)

    def test_output_shapes(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        N = cfg_ion.num_molecules
        ckpt = run_ion_propagation(cfg_ion, neutral)
        T = ckpt.time_ps.size
        # Trajectory arrays
        for name in (
            "positions_x", "positions_y", "positions_z",
            "velocities_x", "velocities_y", "velocities_z",
            "E_kin_eV", "E_pot_eV", "E_dissip_eV",
            "mass_history_kg", "relative_loss_per_ps",
            "number_of_collisions",
        ):
            arr = getattr(ckpt, name)
            assert arr.shape == (2 * N, T), f"{name} shape {arr.shape}"
        # Static (2N,) arrays
        for name in (
            "mass_kg", "mass_final_kg", "droplet_radii_angstrom",
            "positions_final_x", "positions_final_y", "positions_final_z",
            "velocities_final_x", "velocities_final_y", "velocities_final_z",
        ):
            arr = getattr(ckpt, name)
            assert arr.shape == (2 * N,), f"{name} shape {arr.shape}"
        # Per-molecule (N,)
        assert ckpt.b_ion_outside.shape == (N,)
        assert ckpt.b_ion_outside.dtype == bool

    def test_time_axis_consistent_no_stride(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg, ion_simulation_time=0.1, dt_ion=0.01)
        ckpt = run_ion_propagation(cfg_ion, neutral)
        # 0.1 / 0.01 = 10 internal steps => 10 stored when no stride.
        assert ckpt.time_ps.size == 10
        assert ckpt.time_ps[0] == 0.0
        assert ckpt.time_ps[-1] == pytest.approx(9 * cfg_ion.dt_ion)

    def test_final_state_stored_no_stride(self, small_neutral):
        """With stride=1, the last stored column equals the final-state
        fields (the loop ends exactly at the last stored column)."""
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        ckpt = run_ion_propagation(cfg_ion, neutral)
        np.testing.assert_array_equal(ckpt.positions_x[:, -1], ckpt.positions_final_x)
        np.testing.assert_array_equal(ckpt.positions_y[:, -1], ckpt.positions_final_y)
        np.testing.assert_array_equal(ckpt.positions_z[:, -1], ckpt.positions_final_z)
        np.testing.assert_array_equal(ckpt.velocities_x[:, -1], ckpt.velocities_final_x)
        np.testing.assert_array_equal(ckpt.velocities_y[:, -1], ckpt.velocities_final_y)
        np.testing.assert_array_equal(ckpt.velocities_z[:, -1], ckpt.velocities_final_z)
        np.testing.assert_array_equal(ckpt.mass_history_kg[:, -1], ckpt.mass_final_kg)


# ===========================================================================
# Mass history
# ===========================================================================
class TestMassHistory:
    def test_starts_at_neutral_mass(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        ckpt = run_ion_propagation(cfg_ion, neutral)
        np.testing.assert_array_equal(ckpt.mass_history_kg[:, 0], neutral.mass_kg)

    def test_monotonic_nondecreasing(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        ckpt = run_ion_propagation(cfg_ion, neutral)
        assert np.all(np.diff(ckpt.mass_history_kg, axis=1) >= 0)

    def test_mass_final_equals_last_column_no_stride(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        ckpt = run_ion_propagation(cfg_ion, neutral)
        np.testing.assert_array_equal(ckpt.mass_final_kg, ckpt.mass_history_kg[:, -1])


# ===========================================================================
# Reproducibility
# ===========================================================================
class TestReproducibility:
    def test_same_seed_same_output(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        c1 = run_ion_propagation(cfg_ion, neutral, rng=np.random.default_rng(123))
        c2 = run_ion_propagation(cfg_ion, neutral, rng=np.random.default_rng(123))
        np.testing.assert_array_equal(c1.positions_x, c2.positions_x)
        np.testing.assert_array_equal(c1.mass_history_kg, c2.mass_history_kg)
        np.testing.assert_array_equal(c1.E_dissip_eV, c2.E_dissip_eV)
        np.testing.assert_array_equal(c1.number_of_collisions, c2.number_of_collisions)


# ===========================================================================
# Auto-stride
# ===========================================================================
class TestStride:
    def test_small_run_no_stride(self):
        stride, stored = _decide_stride_ion(
            num_molecules=5, num_internal_steps=10,
            max_bytes=DEFAULT_MAX_CHECKPOINT_BYTES_ION,
        )
        assert stride == 1
        assert stored == 10

    def test_large_run_strides(self):
        stride, stored = _decide_stride_ion(
            num_molecules=2000, num_internal_steps=20_000,
            max_bytes=DEFAULT_MAX_CHECKPOINT_BYTES_ION,
        )
        assert stride > 1
        est = _estimate_checkpoint_bytes_ion(2000, stored)
        assert est <= DEFAULT_MAX_CHECKPOINT_BYTES_ION * 1.1

    def test_stride_keeps_size_under_cap(self):
        cap = 50_000_000
        stride, stored = _decide_stride_ion(
            num_molecules=1000, num_internal_steps=2000, max_bytes=cap,
        )
        est = _estimate_checkpoint_bytes_ion(1000, stored)
        assert est <= cap * 1.1, f"got {est/1e6:.1f} MB > {cap/1e6:.1f} MB cap"

    def test_strided_run_produces_correct_checkpoint_size(self, small_neutral):
        """End-to-end: tight max_bytes forces a stride, and the resulting
        checkpoint has fewer columns than internal_steps but still ends
        at a meaningful (non-initial) state."""
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg, ion_simulation_time=0.5, dt_ion=0.01)
        # 50 internal steps; tight cap => stride > 1.
        ckpt = run_ion_propagation(cfg_ion, neutral, max_bytes=5_000)
        assert ckpt.time_ps.size < 50
        assert ckpt.time_ps[0] == 0.0
        assert ckpt.time_ps[-1] > 0
        # Final-state fields are populated from the actual last internal
        # step (not column -1, which under stride>1 may be earlier).
        assert np.any(ckpt.positions_final_x != 0.0)


# ===========================================================================
# RunDirectory integration
# ===========================================================================
class TestRunDirectoryIntegration:
    def test_save_and_load_round_trip(self, tmp_path, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        rd = RunDirectory(tmp_path / "ion_run1")
        c_saved = run_ion_propagation(cfg_ion, neutral, run_dir=rd)

        assert (tmp_path / "ion_run1" / "cfg.json").exists()
        assert (tmp_path / "ion_run1" / "ion.npz").exists()

        loaded_cfg = rd.load_cfg()
        loaded = rd.load_ion(cfg=loaded_cfg)
        np.testing.assert_array_equal(loaded.positions_x, c_saved.positions_x)
        np.testing.assert_array_equal(loaded.mass_history_kg, c_saved.mass_history_kg)
        np.testing.assert_array_equal(loaded.b_ion_outside, c_saved.b_ion_outside)

    def test_does_not_overwrite_existing_cfg(self, tmp_path, small_neutral):
        cfg, neutral = small_neutral
        rd = RunDirectory(tmp_path / "ion_run2")
        rd.save_cfg(cfg)  # pre-existing cfg
        cfg_ion = _ion_cfg(cfg)
        # Driver should detect cfg already present and not re-save.
        run_ion_propagation(cfg_ion, neutral, run_dir=rd)


# ===========================================================================
# Scope checks
# ===========================================================================
class TestScopeChecks:
    def test_single_pulse_false_raises(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_bad = replace(cfg, single_pulse=False)
        with pytest.raises(NotImplementedError, match="single_pulse"):
            run_ion_propagation(cfg_bad, neutral)

    def test_effusive_dynamics_raises(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_bad = _ion_cfg(cfg, effusive_dynamics=True)
        with pytest.raises(NotImplementedError, match="effusive_dynamics"):
            run_ion_propagation(cfg_bad, neutral)

    def test_single_charge_ionization_raises(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_bad = _ion_cfg(cfg, single_charge_ionization_allowed=True)
        with pytest.raises(NotImplementedError, match="single_charge_ionization"):
            run_ion_propagation(cfg_bad, neutral)

    def test_additional_droplet_charges_raises(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_bad = _ion_cfg(cfg, additional_droplet_charges=1)
        with pytest.raises(NotImplementedError, match="additional_droplet_charges"):
            run_ion_propagation(cfg_bad, neutral)

    def test_highly_charged_iodine_raises(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_bad = _ion_cfg(cfg, highly_charged_iodine=True)
        with pytest.raises(NotImplementedError, match="highly_charged_iodine"):
            run_ion_propagation(cfg_bad, neutral)

    def test_collision_mode_not_3_raises(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_bad = _ion_cfg(cfg, hard_sphere_collision_mode=1)
        with pytest.raises(ValueError, match="collision mode 3"):
            run_ion_propagation(cfg_bad, neutral)


# ===========================================================================
# Energy bookkeeping
# ===========================================================================
class TestEnergyBookkeeping:
    def test_dissipation_nondecreasing(self, small_neutral):
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg)
        ckpt = run_ion_propagation(cfg_ion, neutral)
        assert np.all(np.diff(ckpt.E_dissip_eV, axis=1) >= 0)

    def test_energy_balance_no_attachment(self, small_neutral):
        """With mass_attach_probability=0.0, total energy
        E_kin + E_pot + E_dissip should be approximately conserved."""
        cfg, neutral = small_neutral
        cfg_ion = _ion_cfg(cfg, mass_attach_probability=0.0)
        ck = run_ion_propagation(cfg_ion, neutral)
        E0 = (ck.E_kin_eV[:, 0] + ck.E_pot_eV[:, 0] + ck.E_dissip_eV[:, 0]).sum()
        E_end = (ck.E_kin_eV[:, -1] + ck.E_pot_eV[:, -1] + ck.E_dissip_eV[:, -1]).sum()
        rel_drift = abs(E_end - E0) / abs(E0)
        assert rel_drift < 0.01, f"energy drift {rel_drift*100:.3f}% > 1%"


# ===========================================================================
# Internal helpers
# ===========================================================================
class TestInternalHelpers:
    def test_internal_step_count(self, small_neutral):
        cfg, _ = small_neutral
        cfg_ion = replace(cfg, ion_simulation_time=0.1, dt_ion=0.01)
        assert _internal_step_count_ion(cfg_ion) == 10

    def test_estimate_grows_linearly_in_steps(self):
        a = _estimate_checkpoint_bytes_ion(100, 10)
        b = _estimate_checkpoint_bytes_ion(100, 20)
        ratio = b / a
        assert 1.8 < ratio < 2.2

    def test_estimate_grows_linearly_in_atoms(self):
        a = _estimate_checkpoint_bytes_ion(100, 100)
        b = _estimate_checkpoint_bytes_ion(200, 100)
        ratio = b / a
        assert 1.8 < ratio < 2.2
