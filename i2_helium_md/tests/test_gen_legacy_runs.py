"""Tests for scripts/gen_legacy_runs.py."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from i2_helium_md import presets
from i2_helium_md.presets import LEGACY_THESIS_RUNS, same_model
from i2_helium_md.simulation import battery
from i2_helium_md.simulation.battery import RunDirectoryMismatch
from i2_helium_md.simulation.thesis_runs import RUNS_ROOT
import scripts.gen_legacy_runs as script
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
    assert cfg.num_molecules == 2000
    assert cfg.R0_GS_angstrom == pytest.approx(18.0)
    assert cfg.geometric_scattering_crosssection_Iplus == pytest.approx(1600.0)
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


def test_overwriting_is_off_by_default():
    """The committed setting must never replace a finished run."""
    import ast
    from pathlib import Path

    source = Path(__file__).resolve().parents[1] / "scripts" / "gen_legacy_runs.py"
    values = {t.id: n.value for n in ast.parse(source.read_text(encoding="utf-8")).body
              if isinstance(n, ast.Assign) for t in n.targets if isinstance(t, ast.Name)}
    assert ast.literal_eval(values["OVERWRITE_EXISTING_RUN"]) is False


# --- the seeded thesis runs ----------------------------------------------------

THESIS_BASES = {
    "legacy_9A_N2000": (presets.single_pulse_N2000, 2000),
    "legacy_18A_N2000": (presets.single_pulse_N2000_18Angst, 2000),
    "legacy_droplet_N20000": (presets.single_pulse_droplet_distribution, 20000),
}


@pytest.fixture
def thesis(monkeypatch, tmp_path):
    """Thesis settings with output in ``tmp_path``; every MD call is recorded instead."""
    for name, value in dict(INPUT_PRESET="thesis_legacy", RUN_SIZE="production",
                            OUTPUT_ROOT=tmp_path, OVERWRITE_EXISTING_RUN=False,
                            NUM_MOLECULES=4, SEED=1, ION_TIME_PS=0.5).items():
        monkeypatch.setattr(script, name, value)
    for name in script.PRODUCTION_OVERRIDES:
        monkeypatch.setattr(script, name, None)
    calls = {}
    monkeypatch.setattr(script, "run_neutral_and_ion",
                        lambda cfg, run_dir: calls.__setitem__(run_dir, cfg))
    return calls


def test_default_settings_make_the_thesis_runs():
    assert script.INPUT_PRESET == "thesis_legacy"
    assert script.RUN_SIZE == "production"
    assert all(getattr(script, name) is None for name in script.PRODUCTION_OVERRIDES)


def test_thesis_runs_are_the_presets_with_a_size_and_a_seed():
    assert set(LEGACY_THESIS_RUNS) == set(THESIS_BASES)
    for name, (base, n) in THESIS_BASES.items():
        cfg = LEGACY_THESIS_RUNS[name]()
        assert set(same_model(cfg, base())) <= {"num_molecules", "seed"}
        assert cfg.num_molecules == n


def test_thesis_seeds_are_distinct():
    seeds = [build().seed for build in LEGACY_THESIS_RUNS.values()]
    assert None not in seeds and len(set(seeds)) == 3


def test_thesis_legacy_makes_all_three_under_their_names(thesis, tmp_path):
    assert script.main() == 0
    assert set(thesis) == {tmp_path / name for name in THESIS_BASES}
    for name in THESIS_BASES:
        assert same_model(thesis[tmp_path / name], LEGACY_THESIS_RUNS[name]()) == []


def test_one_thesis_run_by_name(thesis, tmp_path, monkeypatch):
    monkeypatch.setattr(script, "INPUT_PRESET", "legacy_18A_N2000")
    assert script.main() == 0
    assert set(thesis) == {tmp_path / "legacy_18A_N2000"}


def test_thesis_smoke_stays_apart_from_the_thesis_runs(thesis, tmp_path, monkeypatch):
    monkeypatch.setattr(script, "RUN_SIZE", "smoke")
    assert script.main() == 0
    assert set(thesis) == {tmp_path / "smoke" / name for name in THESIS_BASES}
    assert all((c.num_molecules, c.ion_simulation_time) == (4, 0.5) for c in thesis.values())


@pytest.mark.parametrize("name, value", [
    ("RUN_SIZE", "custom"), ("RUN_SIZE", "bad"), ("PRODUCTION_SEED", 3),
    ("PRODUCTION_CROSSSECTION", 1.0), ("OVERWRITE_EXISTING_RUN", True)])
def test_settings_that_contradict_the_name_are_refused(thesis, monkeypatch, name, value):
    monkeypatch.setattr(script, name, value)
    with pytest.raises(SystemExit, match=name):
        script.main()
    assert thesis == {}


def test_thesis_run_of_another_model_is_refused_and_untouched(tmp_path, monkeypatch):
    run_dir = tmp_path / "legacy_9A_N2000"
    cfg = dataclasses.replace(LEGACY_THESIS_RUNS["legacy_9A_N2000"](), seed=None,
                              num_molecules=2, ion_simulation_time=0.02)
    battery.run_neutral_and_ion(cfg, run_dir)
    before = {p.name: p.read_bytes() for p in run_dir.iterdir()}
    monkeypatch.setattr(script, "INPUT_PRESET", "legacy_9A_N2000")
    monkeypatch.setattr(script, "OUTPUT_ROOT", tmp_path)
    with pytest.raises(RunDirectoryMismatch, match="seed None"):
        script.main()
    assert {p.name: p.read_bytes() for p in run_dir.iterdir()} == before


def test_stored_thesis_runs_are_recognised_and_nothing_is_run(monkeypatch):
    if not all((RUNS_ROOT / name / "ion.npz").is_file() for name in THESIS_BASES):
        pytest.skip("legacy runs not on disk")
    monkeypatch.setattr(battery, "run_neutral_propagation",
                        lambda *a, **k: pytest.fail("would run MD"))
    assert script.OUTPUT_ROOT == RUNS_ROOT
    assert script.main() == 0


def test_manual_and_thesis_runs_store_the_same_columns(tmp_path, monkeypatch):
    """Same seed and preset: the manual path and the thesis path write equal arrays."""
    cfg = presets.single_pulse_N2000(num_molecules=3, seed=5, ion_simulation_time=0.1)
    battery.run_neutral_and_ion(cfg, tmp_path / "thesis")
    for name, value in dict(INPUT_PRESET="single_pulse_N2000", RUN_SIZE="custom",
                            RUN_DIR=tmp_path / "manual", NUM_MOLECULES=3, SEED=5,
                            ION_TIME_PS=0.1, VERBOSE=False).items():
        monkeypatch.setattr(script, name, value)
    assert script.main() == 0
    for stage in ("neutral.npz", "ion.npz"):
        with np.load(tmp_path / "thesis" / stage) as a, np.load(tmp_path / "manual" / stage) as b:
            assert a.files == b.files
            assert all(np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == "f")
                       for k in a.files)


def test_manual_and_thesis_checkpoint_budgets_are_equal():
    from i2_helium_md.simulation.ion import DEFAULT_MAX_CHECKPOINT_BYTES_ION
    from i2_helium_md.simulation.neutral import DEFAULT_MAX_CHECKPOINT_BYTES

    assert script.NEUTRAL_MAX_BYTES == DEFAULT_MAX_CHECKPOINT_BYTES
    assert script.ION_MAX_BYTES == DEFAULT_MAX_CHECKPOINT_BYTES_ION
