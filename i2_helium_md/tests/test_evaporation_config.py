"""Tests for the Tier-2 Phase-B evaporation-channel config-load guard (Slice Q).

``check_evaporation_config`` mirrors the sibling Phase-A/B guards:

* an ``evap_rrk_dof`` ``s >= 1`` load check (the divergent-rate regime is refused),
  applied **only when the override is set** (``None`` -> per-``n`` ``effective_dof``);
* a ``gate_onset_override_eV`` <-> ``allow_gate_onset_override`` **provenance refuse**
  (mirrors ``allow_unvalidated_binding_pairing``): a forced fixed gate threshold can
  never silently enter a production / Tier-2-lock run (R10 diagnostic-lever, not a knob).

The ``NU_EVAP_PER_PS = 2.42`` constants anchor exists and is the default for
``evap_rate_prefactor_per_ps``.
"""

from __future__ import annotations

import warnings

import pytest

from i2_helium_md.config import SimConfig, check_evaporation_config
from i2_helium_md.physics.constants import NU_EVAP_PER_PS


class TestConstantsAnchor:
    def test_nu_anchor_value(self):
        assert NU_EVAP_PER_PS == 2.42

    def test_default_prefactor_is_the_anchor(self):
        assert SimConfig().evap_rate_prefactor_per_ps == NU_EVAP_PER_PS


class TestDefaults:
    def test_default_evap_surface_is_inert(self):
        cfg = SimConfig()
        assert cfg.evap_rate_prefactor_per_ps == NU_EVAP_PER_PS
        assert cfg.evap_rrk_dof is None
        assert cfg.gate_onset_override_eV is None
        assert cfg.allow_gate_onset_override is False
        check_evaporation_config(cfg)     # does not raise
        cfg.validate()                    # wired into validate()


class TestDofGuard:
    def test_none_is_noop(self):
        check_evaporation_config(SimConfig(evap_rrk_dof=None))

    @pytest.mark.parametrize("s", [1.0, 1.5, 4.0, 60.0])
    def test_dof_at_or_above_one_accepted(self, s):
        check_evaporation_config(SimConfig(evap_rrk_dof=s))

    @pytest.mark.parametrize("s", [0.0, 0.5, -1.0])
    def test_dof_below_one_rejected(self, s):
        with pytest.raises(ValueError, match="evap_rrk_dof"):
            check_evaporation_config(SimConfig(evap_rrk_dof=s))

    def test_reject_propagates_through_validate(self):
        with pytest.raises(ValueError, match="evap_rrk_dof"):
            SimConfig(evap_rrk_dof=0.5).validate()


class TestGateOnsetProvenanceGuard:
    def test_override_none_is_parameter_free(self):
        check_evaporation_config(SimConfig(gate_onset_override_eV=None))

    def test_override_without_allow_refused(self):
        with pytest.raises(ValueError, match="allow_gate_onset_override"):
            check_evaporation_config(SimConfig(gate_onset_override_eV=0.5))

    def test_override_with_allow_warns_not_raises(self):
        cfg = SimConfig(gate_onset_override_eV=0.5, allow_gate_onset_override=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            check_evaporation_config(cfg)
        assert any(issubclass(w.category, RuntimeWarning) for w in caught)

    def test_refuse_propagates_through_validate(self):
        with pytest.raises(ValueError, match="allow_gate_onset_override"):
            SimConfig(gate_onset_override_eV=1.0).validate()
