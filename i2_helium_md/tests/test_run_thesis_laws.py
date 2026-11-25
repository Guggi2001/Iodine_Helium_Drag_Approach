"""Thesis driver: settings are checked before anything runs; smoke runs end to end."""

from __future__ import annotations

import csv
import dataclasses

import numpy as np
import pytest

from i2_helium_md.presets import drag_production_B2, same_model
from i2_helium_md.simulation import battery
from i2_helium_md.simulation.battery import ARTIFACTS, run_member
from i2_helium_md.simulation.detection_stage import load_detection_result
from i2_helium_md.simulation.run_directory import RunDirectory
from scripts import run_thesis_laws as driver

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


@pytest.fixture(autouse=True)
def settings(monkeypatch, tmp_path):
    """Default settings on an empty output root."""
    for name, value in dict(LAWS=("B2", "B3", "h405"), MEMBERS="all", RUN_SIZE="production",
                            CONCURRENCY=1, OUTPUT_ROOT=tmp_path / "runs", TAG="",
                            OVERRIDES={}, SCORES_CSV=None).items():
        monkeypatch.setattr(driver, name, value)


def set_(monkeypatch, **values):
    for name, value in values.items():
        monkeypatch.setattr(driver, name, value)


def test_directory_names(monkeypatch):
    assert driver.member_dir("B2", "s1").name == "B2_N1000_s1"
    assert driver.pool_dir("B2").name == "B2_N6000_pooled"
    set_(monkeypatch, TAG="_rrk10")
    assert driver.member_dir("h405", "s6").name == "h405_rrk10_N1000_s6"
    assert driver.pool_dir("h405").name == "h405_rrk10_N6000_pooled"


def test_overrides_need_a_tag(monkeypatch, tmp_path):
    set_(monkeypatch, OVERRIDES={"evap_rrk_dof": 10.0})
    with pytest.raises(SystemExit, match="TAG"):
        driver.main()
    assert not (tmp_path / "runs").exists()


@pytest.mark.parametrize("field", ["seed", "num_molecules"])
def test_overrides_refuse_the_fields_the_directory_name_states(monkeypatch, field):
    set_(monkeypatch, OVERRIDES={field: 7}, TAG="_x")
    with pytest.raises(SystemExit, match=field):
        driver.main()


def test_unknown_law_lists_the_choices(monkeypatch):
    set_(monkeypatch, LAWS=("B2", "B4"))
    with pytest.raises(KeyError, match="h405"):
        driver.main()


def test_unknown_member_lists_the_choices(monkeypatch):
    set_(monkeypatch, MEMBERS=["s1", "r6"])
    with pytest.raises(KeyError, match="s6"):
        driver.main()


def test_unknown_run_size_is_refused(monkeypatch):
    set_(monkeypatch, RUN_SIZE="tiny")
    with pytest.raises(SystemExit, match="smoke"):
        driver.main()


def test_smoke_runs_end_to_end_without_pool_or_score(monkeypatch, tmp_path):
    set_(monkeypatch, RUN_SIZE="smoke", LAWS=("B2", "h405"), MEMBERS=["s1"])
    assert driver.main() == 0
    root = tmp_path / "runs"
    assert sorted(p.name for p in root.iterdir()) == ["smoke"]
    dirs = sorted((root / "smoke").iterdir())
    assert [d.name for d in dirs] == ["B2_N20_s1", "h405_N20_s1"]
    for run_dir in dirs:
        assert {p.name for p in run_dir.iterdir()} == set(ARTIFACTS)
        assert RunDirectory(run_dir).load_cfg().num_molecules == 20


@pytest.fixture
def short_production(monkeypatch):
    """Production flow with members cut to N = 20, 1 ps and no handover check. Nothing is
    detected that early, so the score is a stand-in (test_battery_score covers it)."""
    real = run_member
    monkeypatch.setattr(battery, "run_member",
                        lambda cfg, run_dir, *, smoke=False: real(cfg, run_dir, smoke=True))
    monkeypatch.setattr(driver, "_build", lambda law, member: drag_production_B2(
        member, num_molecules=20, ion_simulation_time=1.0))
    monkeypatch.setattr(driver, "score_pool", lambda pool, label: {
        **dict.fromkeys(driver.SUMMARY_COLUMNS, 0.0), "label": label, "container": pool.name})


def test_production_pools_scores_and_skips_on_the_second_call(
        monkeypatch, tmp_path, short_production):
    set_(monkeypatch, LAWS=("B2",), SCORES_CSV=tmp_path / "scores.csv")
    assert driver.main() == 0
    pool = driver.pool_dir("B2")
    assert load_detection_result(pool / "detection.npz").num_molecules == 120
    with open(tmp_path / "scores.csv", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert [(r["label"], r["container"]) for r in rows] == [("B2", "B2_N6000_pooled")]
    files = [p for d in (tmp_path / "runs").iterdir() for p in d.iterdir()]
    before = {p: p.stat().st_mtime_ns for p in files}
    assert driver.main() == 0
    assert {p: p.stat().st_mtime_ns for p in files} == before


def test_scoring_one_law_keeps_the_rows_of_the_others(monkeypatch, tmp_path, short_production):
    set_(monkeypatch, LAWS=("B2", "B3"))
    assert driver.main() == 0
    set_(monkeypatch, LAWS=("B3",))
    assert driver.main() == 0
    with open(driver.scores_csv(), newline="", encoding="utf-8") as fh:
        assert [r["label"] for r in csv.DictReader(fh)] == ["B2", "B3"]


def test_a_member_subset_is_not_pooled(monkeypatch, tmp_path, short_production):
    set_(monkeypatch, LAWS=("B2",), MEMBERS=["s1", "s2"])
    assert driver.main() == 0
    assert sorted(p.name for p in (tmp_path / "runs").iterdir()) == [
        "B2_N1000_s1", "B2_N1000_s2"]


def test_score_table_lands_next_to_the_runs_and_carries_the_tag(monkeypatch, tmp_path):
    assert driver.scores_csv() == tmp_path / "runs" / "thesis_scores.csv"
    set_(monkeypatch, TAG="_x")
    assert driver.scores_csv() == tmp_path / "runs" / "thesis_scores_x.csv"


def test_build_applies_the_overrides(monkeypatch):
    set_(monkeypatch, OVERRIDES={"evap_rrk_dof": 10.0}, TAG="_rrk10")
    assert same_model(driver._build("B2", "s3"), drag_production_B2("s3")) == ["evap_rrk_dof"]
    set_(monkeypatch, RUN_SIZE="smoke")
    cfg = driver._build("B2", "s3")
    assert (cfg.num_molecules, cfg.ion_simulation_time) == (20, 1.0)
    assert dataclasses.is_dataclass(cfg) and np.isclose(cfg.evap_rrk_dof, 10.0)
