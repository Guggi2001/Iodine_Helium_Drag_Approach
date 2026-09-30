"""Focused tests for the detector-stage-fix validation generator
(TIER2_DetectorStageFix.md §3b).

Config-construction and guard tests only — no run, no trajectory. They lock
the design: the validation cfg differs from the finals cell in exactly the
two stage keys, Stage I ends at t_h with the seed column equal to the final
step, and the run dir never collides with the committed finals namespace.
Guards that read the committed run dirs (gitignored ``data/runs``) skip when
those are absent.
"""

from __future__ import annotations

import math
import warnings

import pytest

from scripts.gen_detstagefix_validation import (
    ION_TIME_PS,
    T_HANDOVER_PS,
    VALIDATION_DIFF_KEYS,
    VALIDATION_LABELS,
    build_validation_cfg,
    committed_run_dir,
    validation_run_dir_name,
    verify_handover_alignment,
    verify_validation_cfg,
)
from scripts.gen_tier2atlas_g4finals import build_cell, finals_run_dir_name
from scripts.gen_tier2atlas_g4finals import FINAL_MATRIX
from scripts.gen_tier2atlas_linclone import cfg_field_diff
from i2_helium_md.simulation.ion import (
    DEFAULT_MAX_CHECKPOINT_BYTES_ION, _decide_stride_ion, _internal_step_count_ion,
)

_SPEC = {c.label: c for c in FINAL_MATRIX}


def _build(label: str):
    # the documented biphasic mass<->coefficient pairing RuntimeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return build_validation_cfg(label), build_cell(_SPEC[label])


@pytest.mark.parametrize("label", VALIDATION_LABELS)
def test_validation_cfg_differs_only_in_stage_keys(label):
    cfg, base = _build(label)
    assert set(cfg_field_diff(base, cfg)) == VALIDATION_DIFF_KEYS
    assert cfg.relaxation_stage_enabled is False
    assert cfg.ion_simulation_time == ION_TIME_PS


@pytest.mark.parametrize("label", VALIDATION_LABELS)
def test_last_stored_column_is_final_step_at_t_h(label):
    cfg, _ = _build(label)
    steps = _internal_step_count_ion(cfg)
    assert math.isclose((steps - 1) * cfg.dt_ion, T_HANDOVER_PS, abs_tol=1e-9)
    stride, _ = _decide_stride_ion(
        cfg.num_molecules, steps, DEFAULT_MAX_CHECKPOINT_BYTES_ION)
    assert (steps - 1) % stride == 0


@pytest.mark.parametrize("label", VALIDATION_LABELS)
def test_run_dir_distinct_from_committed_finals(label):
    assert validation_run_dir_name(label) != finals_run_dir_name(label)


@pytest.mark.parametrize("label", VALIDATION_LABELS)
def test_guards_against_committed_run(label):
    if not (committed_run_dir(label) / "relaxation.npz").exists():
        pytest.skip("committed finals run dir not present (gitignored data/runs)")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        assert set(verify_validation_cfg(label)) == VALIDATION_DIFF_KEYS
        cfg = build_validation_cfg(label)
    steps, stride = verify_handover_alignment(cfg, label)
    assert (steps - 1) % stride == 0
