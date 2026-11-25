"""Seeded legacy hard-sphere runs and the Tier-0 runs of the thesis figures."""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from i2_helium_md import presets
from i2_helium_md.presets import LEGACY_THESIS_RUNS, same_model
from i2_helium_md.simulation.battery import LOCK_NAME, RunDirectoryMismatch
from i2_helium_md.simulation.run_directory import RunDirectory
from scripts import run_legacy_and_tier0 as driver

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

BASES = {
    "legacy_9A_N2000": (presets.single_pulse_N2000, 2000),
    "legacy_18A_N2000": (presets.single_pulse_N2000_18Angst, 2000),
    "legacy_droplet_N20000": (presets.single_pulse_droplet_distribution, 20000),
}
TIER0_DIRS = {
    "18A_drag_srho3p5A18A_linear_quadratic_N50_B2m21stuck",
    "18A_drag_srho3p5A18A_linear_quadratic_N50_B2m21",
    "9A_drag_srho3p5A18A_linear_quadratic_N50_B2m21",
    "9A_drag_shared_pure_cubic_N50_g11fig_h405pp35m21",
    "9A_drag_srho3p5A18A_linear_quadratic_N50_g11m21B2",
    "9A_drag_srho3p5A18A_linear_quadratic_N50_g11m21B3",
}


@pytest.fixture(autouse=True)
def settings(monkeypatch, tmp_path):
    for name, value in dict(RUN_LEGACY=True, RUN_TIER0=True, RUN_SIZE="production",
                            OUTPUT_ROOT=tmp_path / "runs").items():
        monkeypatch.setattr(driver, name, value)


def legacy_smoke(seed):
    cfg = LEGACY_THESIS_RUNS["legacy_9A_N2000"]()
    return dataclasses.replace(cfg, num_molecules=4, ion_simulation_time=0.5, seed=seed)


def test_legacy_runs_are_the_presets_with_a_size_and_a_seed():
    assert set(LEGACY_THESIS_RUNS) == set(BASES)
    for name, (base, n) in BASES.items():
        cfg = LEGACY_THESIS_RUNS[name]()
        assert set(same_model(cfg, base())) <= {"num_molecules", "seed"}
        assert cfg.num_molecules == n


def test_legacy_seeds_are_distinct():
    seeds = [build().seed for build in LEGACY_THESIS_RUNS.values()]
    assert None not in seeds and len(set(seeds)) == 3


def test_one_seed_repeats_and_another_does_not(tmp_path):
    for name, seed in (("a", 1), ("b", 1), ("c", 2)):
        driver.run_one(legacy_smoke(seed), tmp_path / name)
    with np.load(tmp_path / "a/ion.npz") as a, np.load(tmp_path / "b/ion.npz") as b, \
            np.load(tmp_path / "c/ion.npz") as c:
        assert all(np.array_equal(a[k], b[k], equal_nan=a[k].dtype.kind == "f")
                   for k in a.files)
        assert not np.array_equal(a["velocities_x"], c["velocities_x"])


def test_finished_run_is_skipped(tmp_path):
    cfg, run_dir = legacy_smoke(1), tmp_path / "a"
    driver.run_one(cfg, run_dir)
    before = {p: p.stat().st_mtime_ns for p in run_dir.iterdir()}
    driver.run_one(cfg, run_dir)
    assert {p: p.stat().st_mtime_ns for p in run_dir.iterdir()} == before


def test_an_unseeded_directory_is_refused(tmp_path):
    run_dir = tmp_path / "a"
    driver.run_one(legacy_smoke(None), run_dir)
    assert RunDirectory(run_dir).load_cfg().seed is None
    with pytest.raises(RunDirectoryMismatch, match="seed None"):
        driver.run_one(legacy_smoke(1), run_dir)


def test_a_partial_directory_is_refused(tmp_path):
    run_dir = tmp_path / "a"
    driver.run_one(legacy_smoke(1), run_dir)
    (run_dir / "ion.npz").unlink()
    with pytest.raises(RunDirectoryMismatch, match="ion.npz"):
        driver.run_one(legacy_smoke(1), run_dir)


def test_a_claimed_directory_is_refused(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "a" / LOCK_NAME).touch()
    with pytest.raises(RunDirectoryMismatch, match=LOCK_NAME):
        driver.run_one(legacy_smoke(1), tmp_path / "a")


def test_tier0_runs_have_the_names_the_figures_read():
    assert set(driver.tier0_runs()) == TIER0_DIRS


def test_stored_tier0_runs_are_recognised_and_nothing_is_run(monkeypatch):
    root = driver.PROJECT_ROOT / "data" / "runs"
    if not all((root / name / "ion.npz").is_file() for name in TIER0_DIRS):
        pytest.skip("Tier-0 runs not on disk")

    def refuse(cfg, run_dir):
        raise AssertionError(f"would run {run_dir}")

    monkeypatch.setattr(driver, "run_tier0", refuse)
    monkeypatch.setattr(driver, "OUTPUT_ROOT", root)
    monkeypatch.setattr(driver, "RUN_LEGACY", False)
    assert driver.main() == 0


def test_smoke_runs_everything_apart_from_the_thesis_runs(monkeypatch, tmp_path):
    monkeypatch.setattr(driver, "RUN_SIZE", "smoke")
    assert driver.main() == 0
    root = tmp_path / "runs"
    assert [p.name for p in root.iterdir()] == ["smoke"]
    dirs = list((root / "smoke").iterdir())
    assert {d.name for d in dirs} == set(LEGACY_THESIS_RUNS) | TIER0_DIRS
    assert all(RunDirectory(d).load_cfg().num_molecules == 4 for d in dirs)


def test_unknown_run_size_is_refused(monkeypatch):
    monkeypatch.setattr(driver, "RUN_SIZE", "tiny")
    with pytest.raises(SystemExit, match="smoke"):
        driver.main()
