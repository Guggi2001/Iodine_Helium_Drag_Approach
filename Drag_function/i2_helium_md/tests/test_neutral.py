"""Tests for i2_helium_md/simulation/neutral.py driver."""

from dataclasses import replace
import math
from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.checkpoint import NeutralCheckpoint
from i2_helium_md.simulation.neutral import (
    DEFAULT_MAX_CHECKPOINT_BYTES,
    _decide_stride,
    _estimate_checkpoint_bytes,
    _internal_step_count,
    run_neutral_propagation,
)
from i2_helium_md.simulation.run_directory import RunDirectory


# ===========================================================================
# Public API: returns valid checkpoint
# ===========================================================================
class TestPublicApi:
    def test_returns_neutral_checkpoint(self):
        cfg = single_pulse_N2000(num_molecules=3, seed=42)
        ckpt = run_neutral_propagation(cfg)
        assert isinstance(ckpt, NeutralCheckpoint)

    def test_single_pulse_2_steps(self):
        cfg = single_pulse_N2000(num_molecules=3, seed=42)
        assert cfg.single_pulse
        ckpt = run_neutral_propagation(cfg)
        assert ckpt.time_ps.size == 2
        assert ckpt.time_ps[0] == 0.0
        assert ckpt.time_ps[-1] == pytest.approx(cfg.dt_neutral)

    def test_long_run_consistent_time(self):
        cfg = single_pulse_N2000(num_molecules=5, seed=42)
        cfg = replace(cfg, single_pulse=False, t_max_neutral=0.5)
        ckpt = run_neutral_propagation(cfg)
        # 0.5 / 0.01 = 50 steps
        assert ckpt.time_ps.size == 50
        # final time should be (50-1) * dt = 0.49 ps
        assert ckpt.time_ps[-1] == pytest.approx(49 * cfg.dt_neutral)


# ===========================================================================
# Reproducibility
# ===========================================================================
class TestReproducibility:
    def test_same_seed_same_output(self):
        cfg = single_pulse_N2000(num_molecules=3, seed=42)
        cfg = replace(cfg, single_pulse=False, t_max_neutral=0.05)
        ck1 = run_neutral_propagation(cfg, rng=np.random.default_rng(123))
        ck2 = run_neutral_propagation(cfg, rng=np.random.default_rng(123))
        np.testing.assert_array_equal(ck1.positions_x, ck2.positions_x)
        np.testing.assert_array_equal(ck1.E_dissip_eV, ck2.E_dissip_eV)


# ===========================================================================
# DFT pre-fill stub
# ===========================================================================
class TestDftPrefill:
    def test_custom_dft_start_raises(self):
        cfg = single_pulse_N2000(num_molecules=3, seed=42)
        cfg = replace(cfg, custom_DFT_start=True)
        with pytest.raises(NotImplementedError, match="custom_DFT_start"):
            run_neutral_propagation(cfg)


# ===========================================================================
# Auto-stride
# ===========================================================================
class TestStride:
    def test_small_run_no_stride(self):
        stride, stored = _decide_stride(
            num_molecules=5, num_internal_steps=10,
            max_bytes=DEFAULT_MAX_CHECKPOINT_BYTES,
        )
        assert stride == 1
        assert stored == 10

    def test_large_run_strides(self):
        stride, stored = _decide_stride(
            num_molecules=2000, num_internal_steps=20_000,
            max_bytes=DEFAULT_MAX_CHECKPOINT_BYTES,
        )
        assert stride > 1
        # stored * bytes_per_step must fit in budget
        est = _estimate_checkpoint_bytes(2000, stored)
        assert est <= DEFAULT_MAX_CHECKPOINT_BYTES * 1.1  # allow 10% slack

    def test_stride_keeps_size_under_cap(self):
        """For a run that exceeds the budget, stride must be chosen so
        the resulting checkpoint fits."""
        cap = 50_000_000  # 50 MB cap
        stride, stored = _decide_stride(
            num_molecules=1000, num_internal_steps=2000, max_bytes=cap,
        )
        est = _estimate_checkpoint_bytes(1000, stored)
        assert est <= cap * 1.1, f"got {est/1e6:.1f} MB > {cap/1e6:.1f} MB cap"

    def test_strided_run_produces_correct_checkpoint_size(self):
        """End-to-end: tight max_bytes forces a stride, and the
        resulting checkpoint has fewer steps than internal_steps."""
        cfg = single_pulse_N2000(num_molecules=10, seed=42, R0_GS_angstrom=2.666)
        cfg = replace(cfg, single_pulse=False, t_max_neutral=0.5)
        # internal: 50 steps. With cap = 5000 bytes, stride should kick in.
        ckpt = run_neutral_propagation(cfg, max_bytes=5_000)
        assert ckpt.time_ps.size < 50
        assert ckpt.time_ps[0] == 0.0
        # Last column should be the final state, not zero
        assert ckpt.time_ps[-1] > 0


# ===========================================================================
# RunDirectory integration
# ===========================================================================
class TestRunDirectoryIntegration:
    def test_save_and_load_round_trip(self, tmp_path):
        rd = RunDirectory(tmp_path / "run1")
        cfg = single_pulse_N2000(num_molecules=4, seed=42)
        cfg = replace(cfg, single_pulse=False, t_max_neutral=0.05)
        ck_saved = run_neutral_propagation(cfg, run_dir=rd)

        # Files exist
        assert (tmp_path / "run1" / "cfg.json").exists()
        assert (tmp_path / "run1" / "neutral.npz").exists()

        # Load back
        loaded_cfg = rd.load_cfg()
        loaded = rd.load_neutral(cfg=loaded_cfg)
        np.testing.assert_array_equal(loaded.positions_x, ck_saved.positions_x)
        np.testing.assert_array_equal(loaded.E_dissip_eV, ck_saved.E_dissip_eV)

    def test_does_not_overwrite_existing_cfg(self, tmp_path):
        rd = RunDirectory(tmp_path / "run2")
        cfg_orig = single_pulse_N2000(num_molecules=3, seed=42)
        rd.save_cfg(cfg_orig)
        # Run with the same cfg; should not raise about cfg already there
        cfg_run = single_pulse_N2000(num_molecules=3, seed=42)
        cfg_run = replace(cfg_run, single_pulse=False, t_max_neutral=0.02)
        # Use a freshly constructed cfg with same fields; driver should
        # detect cfg already saved and skip saving (has_cfg() returns True).
        run_neutral_propagation(cfg_run, run_dir=rd)


# ===========================================================================
# Energy bookkeeping
# ===========================================================================
class TestEnergyBookkeeping:
    def test_energy_balance_long_run(self):
        cfg = single_pulse_N2000(num_molecules=20, seed=42, R0_GS_angstrom=2.666)
        cfg = replace(cfg, single_pulse=False, t_max_neutral=0.5)
        ck = run_neutral_propagation(cfg)
        # E_kin + E_pot + E_dissip should be approximately conserved
        E0 = (ck.E_kin_eV[:, 0] + ck.E_pot_eV[:, 0] + ck.E_dissip_eV[:, 0]).sum()
        E_end = (ck.E_kin_eV[:, -1] + ck.E_pot_eV[:, -1] + ck.E_dissip_eV[:, -1]).sum()
        rel_drift = abs(E_end - E0) / abs(E0)
        assert rel_drift < 0.10, f"energy drift {rel_drift*100:.2f}% > 10%"


# ===========================================================================
# Internal helpers
# ===========================================================================
class TestInternalHelpers:
    def test_internal_step_count_single_pulse(self):
        cfg = single_pulse_N2000(num_molecules=3, seed=0)
        assert cfg.single_pulse
        assert _internal_step_count(cfg) == 2

    def test_internal_step_count_long(self):
        cfg = single_pulse_N2000(num_molecules=3, seed=0)
        cfg = replace(cfg, single_pulse=False, t_max_neutral=1.0,
                      dt_neutral=0.01)
        assert _internal_step_count(cfg) == 100

    def test_estimate_grows_linearly_in_steps(self):
        a = _estimate_checkpoint_bytes(100, 10)
        b = _estimate_checkpoint_bytes(100, 20)
        # Should roughly double (small constant overhead from static arrays)
        ratio = b / a
        assert 1.8 < ratio < 2.2

    def test_estimate_grows_linearly_in_atoms(self):
        a = _estimate_checkpoint_bytes(100, 100)
        b = _estimate_checkpoint_bytes(200, 100)
        ratio = b / a
        assert 1.8 < ratio < 2.2
