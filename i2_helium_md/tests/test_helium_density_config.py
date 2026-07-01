"""Tier-2 Phase-B Slice rho -- SimConfig helium-density surface + guard.

Config-construction + guard tests only (no run, no trajectory), mirroring
``test_ladder_config.py``. The guard is a load-time fail-loud enum reject arm
(``Literal`` is not runtime-enforced -- this is the recovery for that, mirroring the
``mass_scenario`` / ``dissociation_ladder`` reject arms). The ``helium_density_profile``
field **repurposes** the former ``Optional[object]`` G4 placeholder; rho stays G2.
"""

import typing

import pytest

from i2_helium_md import SimConfig
from i2_helium_md.config import (
    HeliumDensityProfile,
    check_helium_density_config,
)


def test_default_config_helium_density_surface_is_inert_and_valid():
    cfg = SimConfig()
    assert cfg.helium_density_profile == "erf_complement"
    check_helium_density_config(cfg)   # does not raise
    cfg.validate()                     # full validate path does not raise


@pytest.mark.parametrize("profile", ["erf_complement", "tabulated"])
def test_all_valid_profiles_pass_guard(profile):
    cfg = SimConfig(helium_density_profile=profile)
    check_helium_density_config(cfg)   # no raise
    cfg.validate()


def test_unknown_profile_rejected_by_guard():
    cfg = SimConfig(helium_density_profile="gaussian_typo")
    with pytest.raises(ValueError, match="helium_density_profile"):
        check_helium_density_config(cfg)


def test_unknown_profile_rejected_via_validate():
    cfg = SimConfig(helium_density_profile="gaussian_typo")
    with pytest.raises(ValueError, match="helium_density_profile"):
        cfg.validate()


def test_helium_density_profile_literal_members():
    assert set(typing.get_args(HeliumDensityProfile)) == {
        "erf_complement",
        "tabulated",
    }
