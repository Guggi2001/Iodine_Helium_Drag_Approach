"""Focused tests for the droplet-size sampler A/B generator at h405
(``gen_tier2_detfix_postpickup.py``).

Config-construction and guard tests only — no run, no trajectory. They lock
the design: the A/B flips exactly ``droplet_size_sampler_mode`` relative to
the production battery, keeps its seeds and pipeline, writes to its own run-dir
namespace, and the pooled-container target points at those run dirs. The
committed-cfg oracle reads ``data/runs`` (gitignored) and skips when absent.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from scripts.build_pooled_detection_container import TARGETS
from scripts.gen_tier2_detfix_battery import (
    MEMBER_SEEDS,
    build_member,
    member_run_dir_name as production_run_dir_name,
)
from scripts.gen_tier2_detfix_postpickup import (
    AB_DIFF_KEYS,
    POOL_MEMBERS,
    POSTPICKUP_MEAN_BAND,
    build_postpickup_member,
    member_run_dir_name,
    production_run_dir,
    verify_member_vs_production,
    verify_postpickup_sizes,
)
from scripts.gen_tier2atlas_linclone import cfg_field_diff
from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes

MEMBERS = tuple(MEMBER_SEEDS)


def _build(member):
    # the documented biphasic mass<->coefficient pairing RuntimeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return build_postpickup_member(member), build_member(member)


@pytest.mark.parametrize("member", MEMBERS)
def test_only_the_sampler_mode_differs_from_production(member):
    pp, prod = _build(member)
    assert cfg_field_diff(pp, prod) == sorted(AB_DIFF_KEYS)
    assert (prod.droplet_size_sampler_mode, pp.droplet_size_sampler_mode) == (
        "raw", "post_pickup")
    assert pp.seed == MEMBER_SEEDS[member]


@pytest.mark.parametrize("member", MEMBERS)
def test_run_dirs_are_distinct_from_production(member):
    name = member_run_dir_name(member)
    assert name != production_run_dir_name(member)
    assert f"h405pp{member}_th500" in name


def test_pool_excludes_r6_only():
    assert POOL_MEMBERS == ("s1", "s2", "s3", "s4", "s5")


def test_pooled_targets_pair_the_same_members():
    pp_members, pp_seeds, pp_container, _ = TARGETS["detfix_h405_postpickup"]
    raw_members, raw_seeds, raw_container, _ = TARGETS["detfix_h405_s1s5"]
    assert pp_members == [member_run_dir_name(m) for m in POOL_MEMBERS]
    assert raw_members == [production_run_dir_name(m) for m in POOL_MEMBERS]
    assert pp_seeds == raw_seeds == [MEMBER_SEEDS[m] for m in POOL_MEMBERS]
    containers = {pp_container, raw_container, TARGETS["detfix_h405"][2]}
    assert len(containers) == 3   # never overwrites the N = 6000 production pool


def test_postpickup_size_oracle_separates_the_arms():
    """The size oracle accepts post_pickup and rejects raw at the same seed
    (raw mean ≈ 12.7k He lies below the band)."""
    pp, prod = _build("s1")
    mean_pp = verify_postpickup_sizes(pp)
    assert POSTPICKUP_MEAN_BAND[0] <= mean_pp <= POSTPICKUP_MEAN_BAND[1]
    raw = sample_droplet_sizes(prod, mode="raw", rng=np.random.default_rng(prod.seed))
    assert raw.mean() < POSTPICKUP_MEAN_BAND[0]
    with pytest.raises(AssertionError):
        verify_postpickup_sizes(prod)


def test_unknown_member_refused():
    with pytest.raises(KeyError):
        member_run_dir_name("s9")


@pytest.mark.parametrize("member", MEMBERS)
def test_committed_production_cfg_oracle(member):
    if not (production_run_dir(member) / "cfg.json").exists():
        pytest.skip("production run dir not present")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        assert set(verify_member_vs_production(member)) == AB_DIFF_KEYS
