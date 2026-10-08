"""Focused tests for the droplet-size sampler A/B generator at h405 and on the
lq A5 cell, G5 (``gen_tier2_detfix_postpickup.py``).

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
    B1_E0_EV,
    B2_E0_EV,
    POOL_MEMBERS,
    POSTPICKUP_MEAN_BAND,
    a5_member_cell,
    build_lq_member,
    build_postpickup_member,
    lq_label,
    lq_run_dir,
    main as postpickup_main,
    member_run_dir_name,
    production_run_dir,
    verify_member_vs_production,
    verify_postpickup_sizes,
)
from scripts.gen_tier2_g413_ridge_battery import _E0_KEY
from scripts.gen_tier2_g414_tauline import build as a5_build
from scripts.gen_tier2atlas_g4finals import BUDGET_EV
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


# --- G5: the lq A5 cell (``--cell lq``) -------------------------------------

def _build_lq(member, e0=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return build_lq_member(member, e0), a5_build(a5_member_cell(member))


@pytest.mark.parametrize("e0", [None, B1_E0_EV[1]])
def test_lq_flips_sampler_and_only_replaced_e0(e0):
    pp, a5 = _build_lq("s1", e0)
    want = sorted(AB_DIFF_KEYS | ({_E0_KEY} if e0 is not None else set()))
    assert cfg_field_diff(pp, a5) == want
    assert pp.seed == MEMBER_SEEDS["s1"]
    if e0 is not None:
        assert getattr(pp, _E0_KEY) * BUDGET_EV == pytest.approx(e0, abs=1e-12)


def test_lq_labels_and_run_dirs_are_distinct():
    names = {lq_run_dir(m, e0).name for m in POOL_MEMBERS for e0 in (None, *B1_E0_EV)}
    assert len(names) == len(POOL_MEMBERS) * (1 + len(B1_E0_EV))
    assert lq_label("s1") == "g5lqpps1"
    assert lq_label("s1", 0.3455) == "g5lqppe3455s1"
    assert all("srho3p5A18A_linear_quadratic" in n for n in names)


def test_lq_pooled_targets_pair_the_same_members():
    pp_members, pp_seeds, pp_container, _ = TARGETS["g5_lqpp"]
    _, raw_seeds, raw_container, _ = TARGETS["g5_lq_s1s5"]
    assert pp_members == [lq_run_dir(m).name for m in POOL_MEMBERS]
    assert pp_seeds == raw_seeds == [MEMBER_SEEDS[m] for m in POOL_MEMBERS]
    assert len({pp_container, raw_container, TARGETS["g414_c200t64"][2]}) == 3


def test_b2_refused_until_e0_chosen():
    if B2_E0_EV is not None:
        pytest.skip("B2_E0_EV already set")
    with pytest.raises(SystemExit):
        postpickup_main(["--cell", "lq", "--b2", "--dry-run"])
    with pytest.raises(SystemExit):
        postpickup_main(["--b1", "--dry-run"])   # h405 has no E0 column
