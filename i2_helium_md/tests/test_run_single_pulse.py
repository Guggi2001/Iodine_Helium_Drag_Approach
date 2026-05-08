"""Tests for scripts/run_single_pulse.py."""

from __future__ import annotations

import pytest

import scripts.run_single_pulse as script
from i2_helium_md.simulation.run_directory import RunDirectory


def test_smoke_settings_write_run_directory(tmp_path, monkeypatch):
    run_path = tmp_path / "smoke"
    monkeypatch.setattr(script, "INPUT_PRESET", "single_pulse_N2000")
    monkeypatch.setattr(script, "RUN_SIZE", "smoke")
    monkeypatch.setattr(script, "RUN_DIR", run_path)
    monkeypatch.setattr(script, "OVERWRITE_EXISTING_RUN", False)
    monkeypatch.setattr(script, "NUM_MOLECULES", 2)
    monkeypatch.setattr(script, "SEED", 123)
    monkeypatch.setattr(script, "ION_TIME_PS", 0.02)
    monkeypatch.setattr(script, "VERBOSE", False)

    assert script.main() == 0

    run = RunDirectory(run_path)
    assert run.has_cfg()
    assert run.has_neutral()
    assert run.has_ion()

    cfg = run.load_cfg()
    assert cfg.num_molecules == 2
    assert cfg.seed == 123
    assert cfg.ion_simulation_time == pytest.approx(0.02)

    neutral = run.load_neutral()
    ion = run.load_ion()
    assert neutral.num_molecules == 2
    assert ion.num_molecules == 2
    assert neutral.time_ps.size == 2
    assert ion.time_ps.size == 2


def test_existing_outputs_require_explicit_overwrite(tmp_path, monkeypatch):
    run_path = tmp_path / "existing"
    monkeypatch.setattr(script, "INPUT_PRESET", "single_pulse_N2000")
    monkeypatch.setattr(script, "RUN_SIZE", "smoke")
    monkeypatch.setattr(script, "RUN_DIR", run_path)
    monkeypatch.setattr(script, "OVERWRITE_EXISTING_RUN", False)
    monkeypatch.setattr(script, "NUM_MOLECULES", 2)
    monkeypatch.setattr(script, "SEED", 1)
    monkeypatch.setattr(script, "ION_TIME_PS", 0.02)
    monkeypatch.setattr(script, "VERBOSE", False)
    assert script.main() == 0

    monkeypatch.setattr(script, "SEED", 2)
    with pytest.raises(SystemExit, match="OVERWRITE_EXISTING_RUN"):
        script.main()


def test_overwrite_setting_allows_rerun(tmp_path, monkeypatch):
    run_path = tmp_path / "overwrite"
    monkeypatch.setattr(script, "INPUT_PRESET", "single_pulse_N2000")
    monkeypatch.setattr(script, "RUN_SIZE", "smoke")
    monkeypatch.setattr(script, "RUN_DIR", run_path)
    monkeypatch.setattr(script, "NUM_MOLECULES", 2)
    monkeypatch.setattr(script, "ION_TIME_PS", 0.02)
    monkeypatch.setattr(script, "VERBOSE", False)

    monkeypatch.setattr(script, "SEED", 1)
    monkeypatch.setattr(script, "OVERWRITE_EXISTING_RUN", False)
    assert script.main() == 0

    monkeypatch.setattr(script, "SEED", 2)
    monkeypatch.setattr(script, "OVERWRITE_EXISTING_RUN", True)
    assert script.main() == 0

    cfg = RunDirectory(run_path).load_cfg()
    assert cfg.seed == 2


def test_production_settings_keep_preset_defaults(monkeypatch):
    monkeypatch.setattr(script, "INPUT_PRESET", "single_pulse_N2000")
    monkeypatch.setattr(script, "RUN_SIZE", "production")
    monkeypatch.setattr(script, "PRODUCTION_NUM_MOLECULES", None)
    monkeypatch.setattr(script, "PRODUCTION_SEED", None)
    monkeypatch.setattr(script, "PRODUCTION_ION_TIME_PS", None)

    cfg = script.build_config()
    assert cfg.num_molecules == 2000
    assert cfg.seed is None
    assert cfg.ion_simulation_time == pytest.approx(20.0)


def test_production_can_use_droplet_distribution_preset(monkeypatch):
    monkeypatch.setattr(script, "INPUT_PRESET", "single_pulse_droplet_distribution")
    monkeypatch.setattr(script, "RUN_SIZE", "production")
    monkeypatch.setattr(script, "PRODUCTION_NUM_MOLECULES", None)
    monkeypatch.setattr(script, "PRODUCTION_SEED", None)
    monkeypatch.setattr(script, "PRODUCTION_ION_TIME_PS", None)

    cfg = script.build_config()
    assert cfg.num_molecules == 8000
    assert cfg.R0_GS_angstrom == pytest.approx(2.666)
    assert cfg.E_coulomb_scale == pytest.approx(0.8)
    assert cfg.single_initial_position is False
    assert cfg.use_single_droplet_size is False


def test_production_can_use_18_angst_preset(monkeypatch):
    monkeypatch.setattr(script, "INPUT_PRESET", "single_pulse_N2000_18Angst")
    monkeypatch.setattr(script, "RUN_SIZE", "production")
    monkeypatch.setattr(script, "PRODUCTION_NUM_MOLECULES", None)
    monkeypatch.setattr(script, "PRODUCTION_SEED", None)
    monkeypatch.setattr(script, "PRODUCTION_ION_TIME_PS", None)

    cfg = script.build_config()
    assert cfg.num_molecules == 200
    assert cfg.R0_GS_angstrom == pytest.approx(18.0)
    assert cfg.geometric_scattering_crosssection_Iplus == pytest.approx(200.0)
    assert cfg.binding_energy_I_ion_eV == pytest.approx(0.05)
    assert cfg.mass_attach_probability == pytest.approx(0.005)


def test_custom_uses_selected_preset_with_overrides(monkeypatch):
    monkeypatch.setattr(script, "INPUT_PRESET", "single_pulse_droplet_distribution")
    monkeypatch.setattr(script, "RUN_SIZE", "custom")
    monkeypatch.setattr(script, "NUM_MOLECULES", 12)
    monkeypatch.setattr(script, "SEED", 321)
    monkeypatch.setattr(script, "ION_TIME_PS", 0.5)

    cfg = script.build_config()
    assert cfg.num_molecules == 12
    assert cfg.seed == 321
    assert cfg.ion_simulation_time == pytest.approx(0.5)
    assert cfg.R0_GS_angstrom == pytest.approx(2.666)
    assert cfg.use_single_droplet_size is False


def test_unknown_run_size_raises(monkeypatch):
    monkeypatch.setattr(script, "INPUT_PRESET", "single_pulse_N2000")
    monkeypatch.setattr(script, "RUN_SIZE", "bad")
    with pytest.raises(ValueError, match="RUN_SIZE"):
        script.build_config()


def test_unknown_input_preset_raises(monkeypatch):
    monkeypatch.setattr(script, "INPUT_PRESET", "bad")
    monkeypatch.setattr(script, "RUN_SIZE", "production")
    with pytest.raises(ValueError, match="INPUT_PRESET"):
        script.build_config()
