"""Battery module: run directories are run, skipped or refused, never overwritten."""

from __future__ import annotations

import dataclasses
import os
import shutil
import warnings
from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.presets import (
    THESIS_PRESETS,
    drag_effective_B3,
    drag_production_B2,
    same_model,
)
from i2_helium_md.simulation import battery
from i2_helium_md.simulation.battery import (
    ARTIFACTS,
    LOCK_NAME,
    RunDirectoryMismatch,
    pool_battery,
    pool_members,
    run_battery,
    run_member,
    verify_handover_alignment,
)
from i2_helium_md.simulation.detection_stage import DetectionResult, load_detection_result
from i2_helium_md.simulation.run_directory import RunDirectory

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

RUNS = Path(__file__).resolve().parents[1] / "data" / "runs"


def smoke_cfg(member="s1", build=drag_production_B2, **overrides):
    return build(member, num_molecules=4, ion_simulation_time=0.5, **overrides)


def contents(run_dir):
    return {p.name: p.read_bytes() for p in run_dir.iterdir() if p.is_file()}


def mtimes(run_dir):
    return {p.name: p.stat().st_mtime_ns for p in run_dir.iterdir()}


@pytest.fixture(scope="module")
def finished(tmp_path_factory):
    """One finished smoke member; tests that change it work on a copy."""
    run_dir = tmp_path_factory.mktemp("finished") / "B2_N4_s1"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        run_member(smoke_cfg(), run_dir, smoke=True)
    return run_dir


@pytest.fixture
def member_dir(finished, tmp_path):
    return Path(shutil.copytree(finished, tmp_path / "member"))


# --- run_member --------------------------------------------------------------

def test_smoke_member_writes_every_artifact(finished):
    assert ARTIFACTS == ("cfg.json", "neutral.npz", "ion.npz", "stage1_events.npz",
                         "detection.npz")
    assert {p.name for p in finished.iterdir()} == set(ARTIFACTS)
    assert load_detection_result(finished / "detection.npz").num_molecules == 4
    assert same_model(RunDirectory(finished).load_cfg(), smoke_cfg()) == []


def test_complete_member_is_skipped(member_dir):
    before = mtimes(member_dir)
    assert run_member(smoke_cfg(), member_dir, smoke=True) == member_dir
    assert mtimes(member_dir) == before


def test_extra_entries_are_ignored(member_dir):
    (member_dir / "figures").mkdir()
    run_member(smoke_cfg(), member_dir, smoke=True)


def test_different_model_is_refused_and_untouched(member_dir):
    before = contents(member_dir)
    with pytest.raises(RunDirectoryMismatch, match="evap_rrk_dof"):
        run_member(smoke_cfg(evap_rrk_dof=10.0), member_dir, smoke=True)
    assert contents(member_dir) == before


@pytest.mark.parametrize("missing", ARTIFACTS)
def test_partial_directory_is_refused(member_dir, missing):
    (member_dir / missing).unlink()
    before = contents(member_dir)
    with pytest.raises(RunDirectoryMismatch, match=missing):
        run_member(smoke_cfg(), member_dir, smoke=True)
    assert contents(member_dir) == before


@pytest.mark.parametrize("truncated", [a for a in ARTIFACTS if a.endswith(".npz")])
@pytest.mark.parametrize("keep", [0, 300])
def test_a_file_cut_off_by_a_crash_is_refused(member_dir, truncated, keep):
    path = member_dir / truncated
    path.write_bytes(path.read_bytes()[:keep])
    with pytest.raises(RunDirectoryMismatch, match=f"{member_dir.name}.*{truncated}"):
        run_member(smoke_cfg(), member_dir, smoke=True)


def test_a_corrupt_cfg_is_refused_with_its_path(member_dir):
    (member_dir / "cfg.json").write_text('{"num_mol', encoding="utf-8")
    with pytest.raises(RunDirectoryMismatch, match=f"{member_dir.name}.*cfg.json"):
        run_member(smoke_cfg(), member_dir, smoke=True)


@pytest.mark.parametrize("complete", [True, False])
def test_a_directory_claimed_by_another_process_is_refused(member_dir, tmp_path, complete):
    run_dir = member_dir if complete else tmp_path / "new"
    run_dir.mkdir(exist_ok=True)
    (run_dir / LOCK_NAME).touch()
    before = contents(run_dir)
    with pytest.raises(RunDirectoryMismatch, match=LOCK_NAME):
        run_member(smoke_cfg(), run_dir, smoke=True)
    assert contents(run_dir) == before


def test_a_run_finished_while_waiting_for_the_claim_is_not_overwritten(member_dir):
    """Another process finished it between the first check and the claim."""
    with battery.claimed(member_dir, smoke_cfg(), ARTIFACTS) as start:
        assert not start


def test_a_failed_run_releases_its_claim(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(battery, "run_neutral_propagation", fail)
    with pytest.raises(OSError, match="disk full"):
        run_member(smoke_cfg(), tmp_path / "m", smoke=True)
    assert [p.name for p in (tmp_path / "m").iterdir()] == ["cfg.json"]


def test_a_list_override_is_the_same_model_on_the_second_call(tmp_path):
    cfg = smoke_cfg(tabulated_ladder_rungs_eV=list(
        drag_production_B2().tabulated_ladder_rungs_eV))
    run_member(cfg, tmp_path / "m", smoke=True)
    run_member(cfg, tmp_path / "m", smoke=True)


def test_historical_bundle_label_is_the_same_model(member_dir):
    b3 = smoke_cfg(build=drag_effective_B3)
    RunDirectory(member_dir).save_cfg(dataclasses.replace(
        b3, drag_bundle="LQ",
        drag_coefficients=dataclasses.replace(b3.drag_coefficients, bundle_id="LQ")))
    before = contents(member_dir)
    run_member(b3, member_dir, smoke=True)
    assert contents(member_dir) == before


def test_production_run_refuses_a_window_that_misses_the_handover(tmp_path):
    with pytest.raises(ValueError, match="500"):
        run_member(smoke_cfg(), tmp_path / "m")
    assert not (tmp_path / "m").exists()


def test_thesis_handover_alignment():
    assert verify_handover_alignment(drag_production_B2()) == (50001, 16)


def test_handover_alignment_refuses_another_stride():
    with pytest.raises(ValueError, match="stride"):
        verify_handover_alignment(drag_production_B2(num_molecules=500))


# --- run_battery -------------------------------------------------------------

def test_failure_stops_new_members_and_names_it(member_dir, tmp_path):
    first, third = tmp_path / "first", tmp_path / "third"
    cfgs = {first: smoke_cfg("s2"), member_dir: smoke_cfg(evap_rrk_dof=10.0),
            third: smoke_cfg("s3")}
    with pytest.raises(RuntimeError, match=rf"^failed: {member_dir.name}: ") as failure:
        run_battery(cfgs, concurrency=1, smoke=True)
    assert isinstance(failure.value.__cause__, RunDirectoryMismatch)
    assert all((first / name).is_file() for name in ARTIFACTS)
    assert not third.exists()


def test_parallel_battery_finishes_started_members_and_starts_no_more(member_dir, tmp_path):
    """One slot fails at once while the other still runs: queued members stay unstarted."""
    dirs = [tmp_path / f"m{k}" for k in range(3)]
    slow = drag_production_B2("s2", num_molecules=1000, ion_simulation_time=2.0)
    cfgs = {member_dir: smoke_cfg(evap_rrk_dof=10.0), dirs[0]: slow,
            dirs[1]: smoke_cfg("s3"), dirs[2]: smoke_cfg("s4")}
    with pytest.raises(RuntimeError, match=rf"^failed: {member_dir.name}: ") as failure:
        run_battery(cfgs, concurrency=2, smoke=True)
    assert isinstance(failure.value.__cause__, RunDirectoryMismatch)
    assert all((dirs[0] / name).is_file() for name in ARTIFACTS)
    assert not dirs[1].exists() and not dirs[2].exists()


def test_battery_returns_the_directories_in_order(member_dir, tmp_path):
    cfgs = {tmp_path / "new": smoke_cfg("s2"), member_dir: smoke_cfg()}
    assert run_battery(cfgs, concurrency=2, smoke=True) == list(cfgs)


class KillsItsWorker:
    """Stands in for a config; unpickling it ends the worker process."""

    def __reduce__(self):
        return os._exit, (1,)


def test_a_killed_worker_is_reported_as_such(tmp_path):
    cfgs = {tmp_path / "dead": KillsItsWorker(), tmp_path / "never": smoke_cfg()}
    with pytest.raises(RuntimeError, match="dead.*worker process died") as failure:
        run_battery(cfgs, concurrency=2, smoke=True)
    assert "move the incomplete directories away" in str(failure.value)


def test_battery_refuses_zero_workers(tmp_path):
    with pytest.raises(ValueError, match="concurrency"):
        run_battery({tmp_path / "a": smoke_cfg()}, concurrency=0, smoke=True)


# --- pooling -----------------------------------------------------------------

def differing_fields(a, b):
    return [f.name for f in dataclasses.fields(DetectionResult)
            if not np.array_equal(getattr(a, f.name), getattr(b, f.name))]


@pytest.fixture(scope="module")
def two_members(finished, tmp_path_factory):
    other = tmp_path_factory.mktemp("other") / "B2_N4_s2"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        run_member(smoke_cfg("s2"), other, smoke=True)
    return [finished, other]


@pytest.mark.parametrize("law", list(THESIS_PRESETS))
def test_thesis_pool_is_rebuilt_array_for_array(law, tmp_path):
    members = [RUNS / f"{law}_N1000_s{k}" for k in range(1, 7)]
    committed = RUNS / f"{law}_N6000_pooled"
    if not all((d / "detection.npz").is_file() for d in [*members, committed]):
        pytest.skip("thesis runs not on disk")
    out = pool_battery(members, tmp_path / "pool")
    assert differing_fields(load_detection_result(out / "detection.npz"),
                            load_detection_result(committed / "detection.npz")) == []
    assert same_model(RunDirectory(out).load_cfg(), RunDirectory(committed).load_cfg()) == []


def test_pool_keeps_fragment_pairs_together(two_members):
    members = [load_detection_result(d / "detection.npz") for d in two_members]
    pooled = pool_members(members)
    assert pooled.num_molecules == 8
    for molecule in range(8):
        member, j = members[molecule // 4], molecule % 4
        assert pooled.vx_detected[molecule] == member.vx_detected[j]
        assert pooled.vx_detected[molecule + 8] == member.vx_detected[j + 4]


def test_pool_readme_names_members_and_seeds(two_members, tmp_path):
    out = pool_battery(two_members, tmp_path / "pool")
    text = (out / "README_POOLED.txt").read_text(encoding="utf-8")
    for run_dir in two_members:
        assert f"{run_dir.name} (seed {RunDirectory(run_dir).load_cfg().seed})" in text


def test_pool_refuses_a_duplicate_seed(finished, tmp_path):
    with pytest.raises(ValueError, match="20260730"):
        pool_battery([finished, finished], tmp_path / "pool")
    assert not (tmp_path / "pool").exists()


def test_pool_refuses_an_unseeded_member(two_members, member_dir, tmp_path):
    run = RunDirectory(member_dir)
    run.save_cfg(dataclasses.replace(run.load_cfg(), seed=None))
    with pytest.raises(ValueError, match="None"):
        pool_battery([member_dir, two_members[1]], tmp_path / "pool")


def test_pool_names_the_member_without_detection(two_members, member_dir, tmp_path):
    (member_dir / "detection.npz").unlink()
    with pytest.raises(FileNotFoundError, match=member_dir.name):
        pool_battery([member_dir, two_members[1]], tmp_path / "pool")


def test_pool_refuses_members_of_different_models(two_members, member_dir, tmp_path):
    run = RunDirectory(member_dir)
    run.save_cfg(dataclasses.replace(run.load_cfg(), evap_rrk_dof=10.0))
    with pytest.raises(RunDirectoryMismatch, match="evap_rrk_dof"):
        pool_battery([member_dir, two_members[1]], tmp_path / "pool")


def test_pool_names_a_corrupt_cfg(two_members, member_dir, tmp_path):
    (member_dir / "cfg.json").write_text("{", encoding="utf-8")
    with pytest.raises(RunDirectoryMismatch, match=f"{member_dir.name}.*cfg.json"):
        pool_battery([member_dir, two_members[1]], tmp_path / "pool")


def test_pool_names_a_truncated_detection(two_members, member_dir, tmp_path):
    (member_dir / "detection.npz").write_bytes(b"PK")
    with pytest.raises(RunDirectoryMismatch, match=member_dir.name):
        pool_battery([member_dir, two_members[1]], tmp_path / "pool")


def test_pool_refuses_the_same_ions_under_two_seeds(two_members, member_dir, tmp_path):
    shutil.copyfile(two_members[1] / "detection.npz", member_dir / "detection.npz")
    with pytest.raises(RunDirectoryMismatch, match=f"{member_dir.name}.*{two_members[1].name}"):
        pool_battery([member_dir, two_members[1]], tmp_path / "pool")


def test_pool_refuses_a_detection_of_another_size(two_members, member_dir, tmp_path):
    for d, n in ((member_dir, 7), (other := Path(shutil.copytree(two_members[1], tmp_path / "o")), 7)):
        run = RunDirectory(d)
        run.save_cfg(dataclasses.replace(run.load_cfg(), num_molecules=n))
    with pytest.raises(RunDirectoryMismatch, match="num_molecules"):
        pool_battery([member_dir, other], tmp_path / "pool")


def test_pool_refuses_no_members(tmp_path):
    with pytest.raises(ValueError, match="no members"):
        pool_battery([], tmp_path / "pool")


def test_equal_pool_is_left_as_is(two_members, tmp_path):
    out = pool_battery(two_members, tmp_path / "pool")
    before = mtimes(out)
    assert pool_battery(two_members, out) == out
    assert mtimes(out) == before


def test_different_pool_is_refused(two_members, tmp_path):
    """Same members in another order: same model fields, other arrays."""
    out = pool_battery(two_members, tmp_path / "pool")
    before = contents(out)
    with pytest.raises(RunDirectoryMismatch, match="pool"):
        pool_battery(two_members[::-1], out)
    assert contents(out) == before


def test_pool_with_other_arrays_is_refused(two_members, tmp_path):
    out = pool_battery(two_members, tmp_path / "pool")
    shutil.copyfile(two_members[0] / "detection.npz", out / "detection.npz")
    with pytest.raises(RunDirectoryMismatch, match="num_molecules"):
        pool_battery(two_members, out)


def test_partial_pool_is_refused(two_members, tmp_path):
    out = pool_battery(two_members, tmp_path / "pool")
    (out / "detection.npz").unlink()
    with pytest.raises(RunDirectoryMismatch, match="detection.npz"):
        pool_battery(two_members, out)
