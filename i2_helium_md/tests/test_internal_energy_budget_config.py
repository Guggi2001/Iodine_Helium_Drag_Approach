"""Tier-2 Phase-A Slice U -- SimConfig E_int-budget surface + config-load guard.

Config-construction + guard tests only (no run, no trajectory), mirroring
``test_ladder_config.py`` / ``test_solvation_cooling_config.py``. Per OQ-U1
(decision 2026-06-30), ``internal_energy_partition_fraction`` (f_int) and
``internal_energy_retained_fraction`` (f_ret) both land as ``Optional[float] = None``
declared-but-unread rule-2 stubs (the *module* takes them as kwargs; the Phase-C
driver reads the config). Their committed defaults are Phase-F calibration outputs
(the f_int floor is scenario-split, so no single default serves both budgets). The
guard is a load-time fail-loud **bounded-when-set** check: it does nothing when a
field is None, and asserts ``0 <= f <= 1`` (hard cap; MASS S2 / CALIBRATION rows
13/14) when a value is set. The soft ~0.2 f_int ceiling and the self-unbound floor
are advisory and deliberately NOT enforced.
"""

import pytest

from i2_helium_md import SimConfig
from i2_helium_md.config import check_internal_energy_budget_config


def test_default_config_budget_surface_is_none_and_valid():
    cfg = SimConfig()
    assert cfg.internal_energy_partition_fraction is None
    assert cfg.internal_energy_retained_fraction is None
    check_internal_energy_budget_config(cfg)   # None -> no-op, does not raise
    cfg.validate()                             # full validate path does not raise


@pytest.mark.parametrize("f_int", [0.0, 0.065, 0.21, 0.5, 1.0])
def test_in_range_f_int_accepted(f_int):
    check_internal_energy_budget_config(SimConfig(internal_energy_partition_fraction=f_int))


@pytest.mark.parametrize("f_ret", [0.0, 0.1, 0.4, 1.0])
def test_in_range_f_ret_accepted(f_ret):
    check_internal_energy_budget_config(SimConfig(internal_energy_retained_fraction=f_ret))


@pytest.mark.parametrize("bad", [-0.01, 1.01, 2.0, -1.0])
def test_out_of_range_f_int_rejected(bad):
    cfg = SimConfig(internal_energy_partition_fraction=bad)
    with pytest.raises(ValueError, match="internal_energy_partition_fraction"):
        check_internal_energy_budget_config(cfg)


@pytest.mark.parametrize("bad", [-0.01, 1.01, 2.0, -1.0])
def test_out_of_range_f_ret_rejected(bad):
    cfg = SimConfig(internal_energy_retained_fraction=bad)
    with pytest.raises(ValueError, match="internal_energy_retained_fraction"):
        check_internal_energy_budget_config(cfg)


def test_out_of_range_f_int_rejected_via_validate():
    with pytest.raises(ValueError, match="internal_energy_partition_fraction"):
        SimConfig(internal_energy_partition_fraction=1.5).validate()


def test_out_of_range_f_ret_rejected_via_validate():
    with pytest.raises(ValueError, match="internal_energy_retained_fraction"):
        SimConfig(internal_energy_retained_fraction=-0.2).validate()


def test_sub_floor_f_int_accepted_floor_is_advisory():
    # The self-unbound floor (~0.065-0.35) is advisory, not a constraint (MASS S2 /
    # CALIBRATION row 14): a sub-floor but in-[0,1] f_int passes the guard.
    check_internal_energy_budget_config(SimConfig(internal_energy_partition_fraction=0.01))


def test_both_fields_set_and_valid_pass():
    cfg = SimConfig(
        internal_energy_partition_fraction=0.22,
        internal_energy_retained_fraction=0.4,
    )
    check_internal_energy_budget_config(cfg)
    cfg.validate()


def test_both_set_one_bad_names_the_offender():
    # f_int valid, f_ret out of range -> error must name f_ret (the offender), not f_int.
    cfg = SimConfig(
        internal_energy_partition_fraction=0.22,
        internal_energy_retained_fraction=1.5,
    )
    with pytest.raises(ValueError, match="internal_energy_retained_fraction"):
        check_internal_energy_budget_config(cfg)


def test_boundary_values_zero_and_one_accepted():
    # [0,1] is inclusive (hard cap); both endpoints pass.
    check_internal_energy_budget_config(
        SimConfig(
            internal_energy_partition_fraction=0.0,
            internal_energy_retained_fraction=1.0,
        )
    )
