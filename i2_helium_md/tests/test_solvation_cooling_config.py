"""Tier-2 Phase-A Slice K -- SimConfig solvation-cooling surface + config guard.

Config-construction + guard tests only (no run, no trajectory), mirroring
``test_ladder_config.py``. The guard is a load-time fail-loud check (no silent
clamp): it refuses a non-positive Newton-cooling tau and a non-positive collective
asymptote |S| (a non-positive |S| inverts the binding split). The cooling-time band
[2.6, 16.5] ps is a soft R8 prior and is intentionally *not* enforced (decision
2026-06-30) -- Phase F may probe its edges.
"""

import pytest

from i2_helium_md import SimConfig
from i2_helium_md.config import check_solvation_cooling_config
from i2_helium_md.physics.constants import N_STAR, S_ABS_EV
from i2_helium_md.physics.solvation_cooling import e_electrostriction_eV


def test_default_config_cooling_surface_is_valid():
    cfg = SimConfig()
    assert cfg.solv_struct_asymptote_eV == pytest.approx(S_ABS_EV)
    assert cfg.solv_struct_asymptote_eV > 0.0
    assert cfg.internal_energy_cooling_tau_ps > 0.0
    check_solvation_cooling_config(cfg)   # does not raise
    cfg.validate()                        # full validate path does not raise


def test_default_tau_is_in_R8_band():
    # geometric mid of [2.6, 16.5] ps (decision 2026-06-30); soft prior, not enforced.
    tau = SimConfig().internal_energy_cooling_tau_ps
    assert 2.6 <= tau <= 16.5
    assert tau == pytest.approx(6.55, abs=0.05)


@pytest.mark.parametrize("bad_tau", [0.0, -1.0, -6.55])
def test_nonpositive_tau_rejected_by_guard(bad_tau):
    cfg = SimConfig(internal_energy_cooling_tau_ps=bad_tau)
    with pytest.raises(ValueError, match="tau"):
        check_solvation_cooling_config(cfg)


def test_nonpositive_tau_rejected_via_validate():
    with pytest.raises(ValueError, match="tau"):
        SimConfig(internal_energy_cooling_tau_ps=-2.0).validate()


@pytest.mark.parametrize("good_tau", [2.6, 6.55, 16.5, 20.0])
def test_positive_tau_accepted_band_not_enforced(good_tau):
    # 20.0 is outside the R8 band but accepted: the band is a soft prior, not a guard.
    check_solvation_cooling_config(SimConfig(internal_energy_cooling_tau_ps=good_tau))


@pytest.mark.parametrize("bad_s", [0.0, -0.1])
def test_nonpositive_asymptote_rejected_by_guard(bad_s):
    cfg = SimConfig(solv_struct_asymptote_eV=bad_s)
    with pytest.raises(ValueError, match="solv_struct_asymptote"):
        check_solvation_cooling_config(cfg)


def test_nonpositive_asymptote_rejected_via_validate():
    with pytest.raises(ValueError, match="solv_struct_asymptote"):
        SimConfig(solv_struct_asymptote_eV=-1.0).validate()


def test_electrostriction_binding_property_matches_module():
    # Derived inspection property: E_elec(n*) using the cfg picture/kappa/|S|.
    cfg = SimConfig()
    expected = e_electrostriction_eV(
        N_STAR,
        picture=cfg.ladder_electronic_picture,
        kappa=cfg.ladder_steepness,
        s_abs_eV=cfg.solv_struct_asymptote_eV,
    )
    assert cfg.electrostriction_binding_eV == pytest.approx(expected)
    assert cfg.electrostriction_binding_eV < 0.0   # bound state, non-positive marginal
