"""Tier-0 runs of the thesis figures: names, configs, and that readers run nothing."""

from __future__ import annotations

import ast
import dataclasses
import json
from pathlib import Path

import pytest

from i2_helium_md.presets import same_model
from i2_helium_md.simulation import battery
from i2_helium_md.simulation.battery import LOCK_NAME, RunDirectoryMismatch
from i2_helium_md.simulation.run_directory import RunDirectory
from i2_helium_md.simulation.thesis_runs import RUNS_ROOT
from scripts import gen_tier0_runs, tier0_runs
from scripts.extraction import method_b_g412_ridge as ridge
from scripts.tier0_common import run_dir_name

ROOT = Path(__file__).resolve().parents[1]
B2_9A = "tier0_9A_B2"
NAMES = {"tier0_18A_B2_stuck", "tier0_18A_B2", "tier0_9A_B2", "tier0_9A_B3", "tier0_9A_h405"}
#: Directory of the ridge script's B2 exit-bound run.
RIDGE_B2 = "9A_drag_srho3p5A18A_linear_quadratic_N50_g11m21B2"


def fake_run(path, cfg):
    """A complete run directory of ``cfg`` with empty checkpoints."""
    run = RunDirectory(path)
    run.save_cfg(cfg)
    empty_zip = b"PK\x05\x06" + bytes(18)
    for name in ("neutral.npz", "ion.npz"):
        (path / name).write_bytes(empty_zip)


def contents(path):
    return {p.name: p.read_bytes() for p in path.iterdir()}


def test_names():
    assert set(tier0_runs.tier0_runs()) == NAMES


@pytest.mark.parametrize("run", tier0_runs.CALIBRATION_RUNS, ids=lambda r: r.name)
def test_calibration_configs(run):
    cfg = run.cfg()
    assert cfg.drag_coefficients.coefficients == {"a": 18.007189923460135, "c": 1.6}
    assert cfg.m_eff_amu == 126.90 + 21 * 4.0026
    assert (cfg.num_molecules, cfg.seed) == (50, 20260604)
    assert cfg.ion_simulation_time == run.ion_time_ps
    assert cfg.binding_energy_I_ion_eV == (run.e_bind_eV or 0.11727328616561973)


def test_the_five_runs_are_five_models():
    cfgs = list(tier0_runs.tier0_runs().values())
    assert all(same_model(a, b) for i, a in enumerate(cfgs) for b in cfgs[:i])


def test_names_state_geometry_and_law():
    for name, cfg in tier0_runs.tier0_runs().items():
        _, case, law, *_ = name.split("_")
        assert cfg.R0_GS_angstrom == float(case[:-1])
        want = {"B2": {"a": 18.007189923460135, "c": 1.6}, "B3": {"a": 19.8, "c": 1.6}}
        if law == "h405":
            assert cfg.drag_coefficients.coefficients["v_c"] == 5.4
        else:
            assert cfg.drag_coefficients.coefficients == want[law]


@pytest.mark.parametrize("law", ["B2", "B3", "h405"])
def test_exit_bound_run_of_each_law(law):
    assert tier0_runs.bound_run_name(law) == f"tier0_9A_{law}"


def test_missing_run_names_the_driver(tmp_path):
    with pytest.raises(FileNotFoundError, match="gen_tier0_runs.py"):
        tier0_runs.finished_tier0_run(B2_9A, root=tmp_path)


def test_empty_directory_names_the_driver(tmp_path):
    (tmp_path / B2_9A).mkdir()
    with pytest.raises(FileNotFoundError, match="gen_tier0_runs.py"):
        tier0_runs.finished_tier0_run(B2_9A, root=tmp_path)


def test_unknown_name_lists_the_choices(tmp_path):
    with pytest.raises(KeyError, match="tier0_18A_B2_stuck"):
        tier0_runs.finished_tier0_run("nope", root=tmp_path)
    with pytest.raises(KeyError, match="h405"):
        tier0_runs.bound_run_name("B9")


def test_another_model_is_refused_and_untouched(tmp_path):
    cfg = dataclasses.replace(tier0_runs.tier0_runs()[B2_9A], num_molecules=7)
    fake_run(tmp_path / B2_9A, cfg)
    before = contents(tmp_path / B2_9A)
    with pytest.raises(RunDirectoryMismatch, match="num_molecules"):
        tier0_runs.finished_tier0_run(B2_9A, root=tmp_path)
    assert contents(tmp_path / B2_9A) == before


@pytest.mark.parametrize("spoil", ["partial", "claimed"])
def test_partial_or_claimed_run_is_refused(tmp_path, spoil):
    path = tmp_path / B2_9A
    fake_run(path, tier0_runs.tier0_runs()[B2_9A])
    if spoil == "partial":
        (path / "ion.npz").unlink()
    else:
        (path / LOCK_NAME).touch()
    with pytest.raises(RunDirectoryMismatch):
        tier0_runs.finished_tier0_run(B2_9A, root=tmp_path)


def test_matching_run_is_returned(tmp_path):
    fake_run(tmp_path / B2_9A, tier0_runs.tier0_runs()[B2_9A])
    assert tier0_runs.finished_tier0_run(B2_9A, root=tmp_path).load_cfg().num_molecules == 50


def test_stored_runs_are_the_ones_the_figures_expect():
    if not all((RUNS_ROOT / name / "ion.npz").is_file() for name in NAMES):
        pytest.skip("Tier-0 runs not on disk")
    for name in NAMES:
        tier0_runs.finished_tier0_run(name)


@pytest.mark.parametrize("script", ["tier0_calibration_figures", "tier0_thesis_figure"])
def test_figure_scripts_import_no_run_writer(script):
    path = ROOT / "scripts" / "post_processing" / f"{script}.py"
    names = {a.name for node in ast.walk(ast.parse(path.read_text(encoding="utf-8")))
             if isinstance(node, (ast.Import, ast.ImportFrom)) for a in node.names}
    assert not names & {"run_neutral_propagation", "run_neutral_and_ion"}


def test_tier0_generator_refuses_another_model(tmp_path, monkeypatch):
    """The default settings name a thesis directory: never replaced in place."""
    name = run_dir_name(gen_tier0_runs.CASE, gen_tier0_runs.LAW, gen_tier0_runs.N,
                        gen_tier0_runs.RUN_TAG)
    cfg = dataclasses.replace(tier0_runs.tier0_runs()[B2_9A], num_molecules=7)
    fake_run(tmp_path / name, cfg)
    before = contents(tmp_path / name)
    monkeypatch.setattr(gen_tier0_runs, "OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(gen_tier0_runs, "RUN_SET", "custom")
    with pytest.raises(RunDirectoryMismatch):
        gen_tier0_runs.main()
    assert contents(tmp_path / name) == before


def test_exit_bound_mode_refuses_another_model(tmp_path, monkeypatch):
    name = RIDGE_B2
    fake_run(tmp_path / name, dataclasses.replace(tier0_runs.tier0_runs()[B2_9A], num_molecules=7))
    before = contents(tmp_path / name)
    monkeypatch.setattr(ridge, "RUNS_ROOT", tmp_path)
    monkeypatch.setattr(ridge, "dump_json", lambda row, path: pytest.fail(f"wrote {path}"))
    with pytest.raises(RunDirectoryMismatch):
        ridge.run_bound_m21(["B2"])
    assert contents(tmp_path / name) == before
    assert json.loads(before["cfg.json"])["num_molecules"] == 7


def test_thesis_set_makes_the_five_runs(tmp_path, monkeypatch):
    calls = {}
    monkeypatch.setattr(gen_tier0_runs, "OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(gen_tier0_runs, "run_neutral_and_ion",
                        lambda cfg, run_dir: calls.__setitem__(run_dir.name, cfg))
    assert gen_tier0_runs.RUN_SET == "thesis"
    assert gen_tier0_runs.main() == 0
    assert set(calls) == NAMES
    for name, cfg in tier0_runs.tier0_runs().items():
        assert same_model(calls[name], cfg) == []


def test_custom_set_makes_the_one_run_of_the_settings(tmp_path, monkeypatch):
    calls = {}
    for name, value in dict(OUTPUT_ROOT=tmp_path, RUN_SET="custom", RUN_TAG="mine").items():
        monkeypatch.setattr(gen_tier0_runs, name, value)
    monkeypatch.setattr(gen_tier0_runs, "run_neutral_and_ion",
                        lambda cfg, run_dir: calls.__setitem__(run_dir.name, cfg))
    assert gen_tier0_runs.main() == 0
    assert set(calls) == {"9A_drag_srho3p5A18A_linear_quadratic_N50_mine"}


def test_unknown_set_is_refused(monkeypatch):
    monkeypatch.setattr(gen_tier0_runs, "RUN_SET", "tiny")
    with pytest.raises(SystemExit, match="thesis"):
        gen_tier0_runs.main()


def test_stored_thesis_set_is_recognised_and_nothing_is_run(monkeypatch):
    if not all((RUNS_ROOT / name / "ion.npz").is_file() for name in NAMES):
        pytest.skip("Tier-0 runs not on disk")
    monkeypatch.setattr(battery, "run_neutral_propagation",
                        lambda *a, **k: pytest.fail("would run MD"))
    assert gen_tier0_runs.OUTPUT_ROOT == RUNS_ROOT
    assert gen_tier0_runs.main() == 0


def test_every_thesis_config_runs(tmp_path):
    for name, cfg in tier0_runs.tier0_runs().items():
        battery.run_neutral_and_ion(
            dataclasses.replace(cfg, num_molecules=2, ion_simulation_time=0.05), tmp_path / name)
        assert (tmp_path / name / "ion.npz").is_file()
