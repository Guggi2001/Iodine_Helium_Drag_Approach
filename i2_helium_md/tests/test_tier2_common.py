"""Coverage for the Tier-2 campaign harness (Slice F1, scripts/tier2_common.py).

The builder extends the Phase-D bridge subset with the campaign knobs
(None-sentinel pass-through) and adds the knob-encoding run-tag / run-dir
helpers the run-matrix generator and scoreboard share. These tests stay on the
cheap config surface -- no propagation, no run dirs (that is the F2 smoke test).
"""

from __future__ import annotations

import pytest

from scripts.tier2_common import (
    BRIDGE_F_INT,
    BRIDGE_F_RET,
    BRIDGE_LAMBDA0_PER_PS,
    EXPERIMENTAL_RELAXATION_TIME_PS,
    PRODUCTION_BUDGET_EV,
    TIER2_BRIDGE_TAG,
    VALIDATION_BUDGET_EV,
    build_biphasic_cfg,
    tier2_bridge_run_dir_name,
    tier2_run_dir_name,
    tier2_run_tag,
)
from i2_helium_md.config import SimConfig
from i2_helium_md.physics.internal_energy_budget import f_int_floor


_COMMON = dict(num_molecules=2, ion_time_ps=0.02, dt_ion_ps=0.01, seed=123)


def _default_cfg() -> SimConfig:
    """A pristine SimConfig for reading the ride-the-default values."""
    return SimConfig()


# ---------------------------------------------------------------------------
# build_biphasic_cfg -- scenario + budget + generative knobs
# ---------------------------------------------------------------------------


def test_bridge_call_rides_config_defaults_and_stays_reproducible():
    """A call with only the bridge kwargs must be byte-identical to the bridge:
    picture/kappa/tau ride the config defaults; relaxation stays off; budget 0.80."""
    cfg = build_biphasic_cfg("9A", "shared_pure_cubic", **_COMMON)
    defaults = _default_cfg()

    cfg.validate()
    assert cfg.mass_scenario == "biphasic"
    assert cfg.coulomb_available_eV == VALIDATION_BUDGET_EV == 0.80
    assert cfg.pickup_rate_coefficient == pytest.approx(BRIDGE_LAMBDA0_PER_PS)
    assert cfg.internal_energy_partition_fraction == pytest.approx(BRIDGE_F_INT)
    assert cfg.internal_energy_retained_fraction == pytest.approx(BRIDGE_F_RET)
    assert cfg.allow_inconsistent_mass_pairing is True
    # None-sentinel: the campaign knobs were not passed -> config defaults.
    assert cfg.ladder_electronic_picture == defaults.ladder_electronic_picture
    assert cfg.ladder_steepness == defaults.ladder_steepness
    assert cfg.internal_energy_cooling_tau_ps == defaults.internal_energy_cooling_tau_ps
    # Relaxation stage stays off (default scope; bridge byte-reproducible).
    assert cfg.relaxation_stage_enabled is False
    assert cfg.relaxation_time_ps is None


def test_campaign_knobs_override_config_defaults():
    cfg = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        picture="x2_only",
        kappa=4.0,
        tau_ps=12.5,
        lambda0_per_ps=0.7,
        f_int=0.24,
        f_ret=0.05,
        **_COMMON,
    )
    assert cfg.ladder_electronic_picture == "x2_only"
    assert cfg.ladder_steepness == pytest.approx(4.0)
    assert cfg.internal_energy_cooling_tau_ps == pytest.approx(12.5)
    assert cfg.pickup_rate_coefficient == pytest.approx(0.7)
    assert cfg.internal_energy_partition_fraction == pytest.approx(0.24)
    assert cfg.internal_energy_retained_fraction == pytest.approx(0.05)


def test_production_budget_stamp():
    cfg = build_biphasic_cfg(
        "9A", "shared_pure_cubic", coulomb_available_eV=PRODUCTION_BUDGET_EV, **_COMMON
    )
    assert cfg.coulomb_available_eV == 2.70


def test_relaxation_time_enables_stage_and_stamps_cap():
    cfg = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        relaxation_time_ps=EXPERIMENTAL_RELAXATION_TIME_PS,
        **_COMMON,
    )
    cfg.validate()
    assert cfg.relaxation_stage_enabled is True
    assert cfg.relaxation_time_ps == pytest.approx(8.53e6)
    # forces ride the config default (coulomb) unless overridden.
    assert cfg.relaxation_forces == _default_cfg().relaxation_forces == "coulomb"


def test_relaxation_forces_override():
    cfg = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        relaxation_time_ps=5.0,
        relaxation_forces="free_flight",
        **_COMMON,
    )
    assert cfg.relaxation_forces == "free_flight"


def test_experimental_relaxation_time_is_8530_ns_in_ps():
    # Provenance guard: 8530 ns expressed in ps (config field unit).
    assert EXPERIMENTAL_RELAXATION_TIME_PS == pytest.approx(8530.0 * 1e3)


# ---------------------------------------------------------------------------
# run-tag / run-dir helpers
# ---------------------------------------------------------------------------


def test_run_tag_encodes_all_grid_coordinates():
    tag = tier2_run_tag(
        picture="statistical_mixture",
        kappa=2.0,
        lambda0_per_ps=0.9,
        f_int=0.24,
        f_ret=0.10,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    assert tag == "tier2_b080_mix_k2.00_l0.90_fi0.24_fr0.10_tau6.55"


def test_run_tag_budget_and_picture_abbreviations():
    assert tier2_run_tag(
        picture="x2_only", kappa=1.0, lambda0_per_ps=0.9, f_int=0.1, f_ret=0.1,
        tau_ps=6.5, budget_eV=2.70,
    ).startswith("tier2_b270_x2_")
    assert tier2_run_tag(
        picture="cooling_relaxed", kappa=1.0, lambda0_per_ps=0.9, f_int=0.1, f_ret=0.1,
        tau_ps=6.5, budget_eV=0.80,
    ).startswith("tier2_b080_cool_")


def test_run_tag_total_strip_suffix():
    base = dict(
        picture="statistical_mixture", kappa=1.0, lambda0_per_ps=0.9, f_int=0.1,
        f_ret=0.1, tau_ps=6.5, budget_eV=0.80,
    )
    assert not tier2_run_tag(**base).endswith("_totalstrip")
    assert tier2_run_tag(**base, total_strip=True).endswith("_totalstrip")


def test_run_tag_appends_and_omits_s_eff_suffix():
    """s_eff promotion (2026-07-06): a set evap_rrk_dof appends _sNN.NN to the
    campaign tag; None appends nothing (delivered campaign tags byte-identical);
    the s_eff suffix precedes the total_strip variant marker."""
    base = dict(
        picture="statistical_mixture", kappa=1.0, lambda0_per_ps=0.9, f_int=0.5,
        f_ret=0.1, tau_ps=6.55, budget_eV=0.80,
    )
    assert tier2_run_tag(**base, evap_rrk_dof=None) == tier2_run_tag(**base)
    assert tier2_run_tag(**base, evap_rrk_dof=8.0).endswith("_tau6.55_s8.00")
    assert tier2_run_tag(
        **base, evap_rrk_dof=8.0, total_strip=True
    ).endswith("_s8.00_totalstrip")


def test_run_dir_name_carries_s_eff():
    """Two campaign grid points differing only in s_eff must not collide."""
    common = dict(
        picture="statistical_mixture", kappa=1.0, lambda0_per_ps=0.9, f_int=0.5,
        f_ret=0.1, tau_ps=6.55, budget_eV=0.80,
    )
    dirs = {
        tier2_run_dir_name("9A", "shared_pure_cubic", 500, **common, evap_rrk_dof=s)
        for s in (None, 5.0, 8.0, 12.0)
    }
    assert len(dirs) == 4


def test_run_tag_rejects_unknown_picture():
    with pytest.raises(ValueError, match="picture"):
        tier2_run_tag(
            picture="bogus", kappa=1.0, lambda0_per_ps=0.9, f_int=0.1, f_ret=0.1,
            tau_ps=6.5, budget_eV=0.80,
        )


def test_run_tag_separates_tau_at_two_decimals():
    """Review finding #1: the Stage-2 tau sweep must not alias 6.55 with 6.5."""
    common = dict(
        picture="statistical_mixture", kappa=1.0, lambda0_per_ps=0.9, f_int=0.24,
        f_ret=0.1, budget_eV=0.80,
    )
    assert tier2_run_tag(**common, tau_ps=6.55) != tier2_run_tag(**common, tau_ps=6.50)


def test_run_tag_separates_lambda0():
    """Review finding #2: two campaigns differing only in lambda_0 must not collide."""
    common = dict(
        picture="statistical_mixture", kappa=1.0, f_int=0.24, f_ret=0.1,
        tau_ps=6.55, budget_eV=0.80,
    )
    assert tier2_run_tag(**common, lambda0_per_ps=0.7) != tier2_run_tag(**common, lambda0_per_ps=1.1)


def test_picture_tags_stay_in_lockstep_with_config_enum():
    """Review finding #7: the abbreviation map must cover exactly the config enum."""
    import typing

    from i2_helium_md.config import LadderElectronicPicture
    from scripts.tier2_common import _PICTURE_TAGS

    assert set(_PICTURE_TAGS) == set(typing.get_args(LadderElectronicPicture))


def test_run_dir_name_shares_tier0_convention():
    name = tier2_run_dir_name(
        "9A",
        "shared_pure_cubic",
        500,
        picture="statistical_mixture",
        kappa=2.0,
        lambda0_per_ps=0.9,
        f_int=0.24,
        f_ret=0.10,
        tau_ps=6.55,
        budget_eV=0.80,
    )
    assert name == (
        "9A_drag_shared_pure_cubic_N500_tier2_b080_mix_k2.00_l0.90_fi0.24_fr0.10_tau6.55"
    )


@pytest.mark.parametrize("budget_eV", [VALIDATION_BUDGET_EV, PRODUCTION_BUDGET_EV])
@pytest.mark.parametrize("picture", ["statistical_mixture", "x2_only", "cooling_relaxed"])
@pytest.mark.parametrize("kappa", [0.5, 1.0, 2.0, 4.0, 8.0])
def test_every_grid_point_builds_with_floor_f_int(budget_eV, picture, kappa):
    """Grid-wide guarantee: at the floor f_int pin, every documented grid point
    builds and validates at both budgets, and the floor lands inside (0, 1] so it
    never trips the config [0,1] bound."""
    floor = f_int_floor(e_avail_eV=budget_eV, picture=picture, kappa=kappa)
    assert 0.0 < floor <= 1.0
    cfg = build_biphasic_cfg(
        "9A",
        "shared_pure_cubic",
        picture=picture,
        kappa=kappa,
        f_int=floor,
        coulomb_available_eV=budget_eV,
        **_COMMON,
    )
    cfg.validate()
    assert cfg.mass_scenario == "biphasic"
    assert cfg.coulomb_available_eV == budget_eV
    assert cfg.internal_energy_partition_fraction == pytest.approx(floor)


def test_bridge_tag_reserved_and_unchanged():
    """Regression lock: F1 must not repurpose the reserved Phase-D bridge tag."""
    assert TIER2_BRIDGE_TAG == "tier2_bridge_biphasic"
    assert (
        tier2_bridge_run_dir_name("9A", "shared_pure_cubic", 50)
        == "9A_drag_shared_pure_cubic_N50_tier2_bridge_biphasic"
    )
    # The campaign tag namespace ("tier2_bNNN_...") never collides with it.
    assert not tier2_run_tag(
        picture="statistical_mixture", kappa=1.0, lambda0_per_ps=0.9, f_int=0.5,
        f_ret=0.1, tau_ps=6.55, budget_eV=0.80,
    ).startswith(TIER2_BRIDGE_TAG)


def test_grid_points_get_distinct_run_dirs():
    def name(picture, kappa, budget):
        return tier2_run_dir_name(
            "9A", "shared_pure_cubic", 500,
            picture=picture, kappa=kappa, lambda0_per_ps=0.9, f_int=0.24, f_ret=0.1,
            tau_ps=6.55, budget_eV=budget,
        )

    names = {
        name(p, k, b)
        for p in ("statistical_mixture", "x2_only", "cooling_relaxed")
        for k in (0.5, 1.0, 2.0, 4.0, 8.0)
        for b in (0.80, 2.70)
    }
    assert len(names) == 3 * 5 * 2  # every grid point is uniquely named
