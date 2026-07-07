"""Tier-2 drag-port -- SimConfig cooling-spatial-gate arm + guard.

Config-construction + guard tests only (no run, no trajectory), mirroring
``test_helium_density_config.py``. The guard is a load-time fail-loud enum reject
arm (``Literal`` is not runtime-enforced -- this is the recovery for that, mirroring
the ``helium_density_profile`` / ``mass_scenario`` reject arms). The
``cooling_spatial_gate`` field selects whether the K2 Newton-cooling drain is
spatially gated by the shared erf-complement He-density surface: ``none`` (default,
byte-identical to the locked ungated bath dissipation) or ``density_scaled`` (cooling
attenuates as ``rho_He/rho_bulk`` and switches off outside the droplet).
"""

import typing

import pytest

from i2_helium_md import SimConfig
from i2_helium_md.config import (
    CoolingSpatialGate,
    check_cooling_spatial_gate_config,
)


def test_default_config_cooling_gate_is_inert_and_valid():
    cfg = SimConfig()
    assert cfg.cooling_spatial_gate == "none"
    check_cooling_spatial_gate_config(cfg)   # does not raise
    cfg.validate()                           # full validate path does not raise


@pytest.mark.parametrize("gate", ["none", "density_scaled"])
def test_all_valid_gates_pass_guard(gate):
    cfg = SimConfig(cooling_spatial_gate=gate)
    check_cooling_spatial_gate_config(cfg)   # no raise
    cfg.validate()


def test_unknown_gate_rejected_by_guard():
    cfg = SimConfig(cooling_spatial_gate="density_typo")
    with pytest.raises(ValueError, match="cooling_spatial_gate"):
        check_cooling_spatial_gate_config(cfg)


def test_unknown_gate_rejected_via_validate():
    cfg = SimConfig(cooling_spatial_gate="density_typo")
    with pytest.raises(ValueError, match="cooling_spatial_gate"):
        cfg.validate()


def test_cooling_spatial_gate_literal_members():
    from i2_helium_md.config import _KNOWN_COOLING_SPATIAL_GATES

    members = set(typing.get_args(CoolingSpatialGate))
    assert members == {"none", "density_scaled"}
    assert members == set(_KNOWN_COOLING_SPATIAL_GATES)
