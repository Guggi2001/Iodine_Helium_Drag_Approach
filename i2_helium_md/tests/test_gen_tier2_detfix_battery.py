"""Focused tests for the detector-stage-fix production battery generator
(TIER2_DetectorStageFix.md §5.4 D6).

Config-construction and guard tests only — no run, no trajectory. They lock
the design: the production pipeline keys, handover at exactly 500.00 ps with
stride 16 so the detection seed is the final state, the seed table (E2-era
battery twins + r6), and run dirs outside the E2-era namespace. The twin
oracle reads the committed run dirs (gitignored ``data/runs``) and skips when
they are absent.
"""

from __future__ import annotations

import math
import warnings

import pytest

from scripts.gen_tier2_detfix_battery import (
    MAX_BYTES_ION,
    MEMBER_SEEDS,
    STORAGE_STRIDE,
    T_HANDOVER_PS,
    TWIN_DIFF_KEYS,
    build_member,
    member_run_dir_name,
    twin_run_dir,
    verify_handover_alignment,
    verify_member_vs_twin,
)
from scripts.gen_tier2atlas_g4finals import build_cell, finals_run_dir_name
from scripts.gen_tier2atlas_g4step2_battery import (
    H405, member_run_dir_name as e2_member_run_dir_name,
)
from scripts.gen_tier2atlas_linclone import cfg_field_diff
from i2_helium_md.simulation.ion import _decide_stride_ion, _internal_step_count_ion

MEMBERS = tuple(MEMBER_SEEDS)


def _build(member):
    # the documented biphasic mass<->coefficient pairing RuntimeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return build_member(member), build_cell(H405)


def test_seed_table():
    assert MEMBER_SEEDS == {"s1": 20260730, "s2": 20260731, "s3": 20260732,
                            "s4": 20260733, "s5": 20260734, "r6": 20260729}


@pytest.mark.parametrize("member", MEMBERS)
def test_member_differs_from_h405_cell_only_in_pipeline_keys_and_seed(member):
    cfg, base = _build(member)
    assert set(cfg_field_diff(base, cfg)) == (
        TWIN_DIFF_KEYS | {"detection_coulomb_closure", "seed"}
        if member != "r6" else TWIN_DIFF_KEYS | {"detection_coulomb_closure"}
    )
    assert cfg.relaxation_stage_enabled is False
    assert cfg.detection_coulomb_closure == "partner_aware"
    assert cfg.evaporation_shed_convention == "co_moving"
    assert cfg.seed == MEMBER_SEEDS[member]


@pytest.mark.parametrize("member", MEMBERS)
def test_handover_at_500_ps_with_stride_16(member):
    cfg, _ = _build(member)
    steps, stride = verify_handover_alignment(cfg)
    assert steps == 50001 and stride == STORAGE_STRIDE == 16
    assert math.isclose((steps - 1) * cfg.dt_ion, T_HANDOVER_PS, abs_tol=1e-9)
    assert _decide_stride_ion(cfg.num_molecules, steps, MAX_BYTES_ION)[0] == 16


@pytest.mark.parametrize("member", MEMBERS)
def test_run_dir_outside_e2_era_namespace(member):
    name = member_run_dir_name(member)
    assert "detfix" in name and "tier2atlas" not in name
    twin = (finals_run_dir_name("h405") if member == "r6"
            else e2_member_run_dir_name(member))
    assert twin_run_dir(member).name == twin != name


@pytest.mark.parametrize("member", MEMBERS)
def test_twin_oracle_against_committed_cfg(member):
    if not (twin_run_dir(member) / "cfg.json").exists():
        pytest.skip("E2-era twin run dir not present (gitignored data/runs)")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        assert set(verify_member_vs_twin(member)) == TWIN_DIFF_KEYS
