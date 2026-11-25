"""Thesis presets B2 / B3 / h405: equal to the stored thesis runs, field by field."""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from i2_helium_md.presets import (
    THESIS_MEMBER_SEEDS,
    THESIS_PRESETS,
    drag_effective_B3,
    drag_production_B2,
    same_model,
)
from i2_helium_md.simulation.run_directory import RunDirectory

RUNS = Path(__file__).resolve().parents[1] / "data" / "runs"
FIXTURES = Path(__file__).resolve().parent / "fixtures" / "thesis_laws"
MEMBERS = ["s1", "s2", "s3", "s4", "s5", "s6"]


@pytest.mark.parametrize("law", ["B2", "B3", "h405"])
@pytest.mark.parametrize("member", MEMBERS)
def test_preset_equals_the_stored_thesis_cfg(law, member):
    run = RUNS / f"{law}_N1000_{member}"
    if not (run / "cfg.json").is_file():
        pytest.skip("thesis run not on disk")
    assert same_model(THESIS_PRESETS[law](member), RunDirectory(run).load_cfg()) == []


@pytest.mark.parametrize("law", ["B2", "B3", "h405"])
def test_preset_equals_the_committed_law(law):
    """Needs no run on disk: ``fixtures/thesis_laws`` holds the stored s1 configs."""
    assert same_model(THESIS_PRESETS[law](), RunDirectory(FIXTURES / law).load_cfg()) == []


def test_member_seeds_are_the_thesis_seeds():
    assert list(THESIS_MEMBER_SEEDS.values()) == [
        20260730, 20260731, 20260732, 20260733, 20260734, 20260729]


def test_b2_and_b3_differ_only_in_their_law():
    assert set(same_model(drag_production_B2(), drag_effective_B3())) == {
        "binding_energy_I_ion_eV", "drag_coefficients", "internal_energy_partition_fraction"}


def test_members_differ_only_in_seed():
    assert same_model(drag_production_B2("s1"), drag_production_B2("s6")) == ["seed"]


def test_six_distinct_member_seeds():
    assert list(THESIS_MEMBER_SEEDS) == MEMBERS
    assert len(set(THESIS_MEMBER_SEEDS.values())) == 6


def test_unknown_member_lists_the_choices():
    with pytest.raises(KeyError, match="s6"):
        drag_production_B2("r6")


def test_overrides_apply_and_every_preset_validates():
    assert drag_production_B2(evap_rrk_dof=10.0).evap_rrk_dof == 10.0
    for build in THESIS_PRESETS.values():
        build().validate()


def test_same_model_sees_a_changed_coefficient():
    b2 = drag_production_B2()
    c = dataclasses.replace(b2.drag_coefficients, coefficients={"a": 18.0, "c": 1.6})
    assert "drag_coefficients" in same_model(b2, dataclasses.replace(b2, drag_coefficients=c))


def test_same_model_ignores_the_bundle_label():
    b3 = drag_effective_B3()
    relabelled = dataclasses.replace(
        b3, drag_bundle="LQ",
        drag_coefficients=dataclasses.replace(b3.drag_coefficients, bundle_id="LQ"))
    assert same_model(b3, relabelled) == []


def test_same_model_handles_a_run_without_drag_coefficients():
    b2 = drag_production_B2()
    bare = dataclasses.replace(b2, drag_coefficients=None)
    assert same_model(bare, bare) == []
    assert same_model(b2, bare) == ["drag_coefficients"]
