"""Tests for the Tier-2 Phase-B pickup-channel config-load guard (Slice P).

``check_pickup_config`` is a load-time fail-loud enum reject arm over the three pickup
selectors (mirrors the ``mass_scenario`` / ``dissociation_ladder`` /
``helium_density_profile`` reject arms). It accepts the valid-but-unbuilt rule-2 arms
(``sweeping`` / ``dwell_time`` / ``thermal``) -- those raise lazily at point-of-use in
``physics/pickup.py``, not at config-load -- and rejects only genuine typos.
"""

from __future__ import annotations

import typing

import pytest

from i2_helium_md.config import (
    SimConfig,
    PickupRateForm,
    PickupOccupancyCap,
    HeCaptureVelocity,
    check_pickup_config,
)


class TestDefaults:
    def test_default_pickup_surface_is_inert(self):
        cfg = SimConfig()
        assert cfg.pickup_rate_form == "density_only"
        assert cfg.pickup_rate_coefficient == 0.0
        assert cfg.pickup_occupancy_cap == "langmuir"
        assert cfg.pickup_occupancy_exponent == 1.0
        assert cfg.he_capture_velocity == "at_rest"
        check_pickup_config(cfg)          # does not raise
        cfg.validate()                    # wired into validate()

    def test_retired_fields_are_gone(self):
        cfg = SimConfig()
        assert not hasattr(cfg, "mass_rate_form")
        assert not hasattr(cfg, "mass_rate_coefficient")
        assert not hasattr(cfg, "mass_relaxation_tau_ps")


class TestAcceptsValidMembers:
    @pytest.mark.parametrize("form", ["density_only", "sweeping", "dwell_time"])
    def test_rate_form_members_accepted(self, form):
        check_pickup_config(SimConfig(pickup_rate_form=form))   # incl. unbuilt rule-2 arms

    @pytest.mark.parametrize("cap", ["langmuir", "none"])
    def test_occupancy_cap_members_accepted(self, cap):
        check_pickup_config(SimConfig(pickup_occupancy_cap=cap))

    @pytest.mark.parametrize("vel", ["at_rest", "thermal"])
    def test_capture_velocity_members_accepted(self, vel):
        check_pickup_config(SimConfig(he_capture_velocity=vel))  # thermal accepted at load


class TestRejectsTypos:
    def test_bad_rate_form_rejected(self):
        with pytest.raises(ValueError, match="pickup_rate_form"):
            check_pickup_config(SimConfig(pickup_rate_form="densityonly_typo"))

    def test_bad_cap_rejected(self):
        with pytest.raises(ValueError, match="pickup_occupancy_cap"):
            check_pickup_config(SimConfig(pickup_occupancy_cap="sigmoid_typo"))

    def test_bad_capture_velocity_rejected(self):
        with pytest.raises(ValueError, match="he_capture_velocity"):
            check_pickup_config(SimConfig(he_capture_velocity="cold_typo"))

    def test_reject_propagates_through_validate(self):
        with pytest.raises(ValueError, match="pickup_rate_form"):
            SimConfig(pickup_rate_form="nope").validate()


class TestLiteralMembers:
    def test_literal_member_sets(self):
        assert set(typing.get_args(PickupRateForm)) == {
            "density_only", "sweeping", "dwell_time"
        }
        assert set(typing.get_args(PickupOccupancyCap)) == {"langmuir", "none"}
        assert set(typing.get_args(HeCaptureVelocity)) == {"at_rest", "thermal"}
