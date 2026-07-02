"""Tests for the Tier-2 Phase-C biphasic driver config-load guard (Slice G).

``check_biphasic_config`` fires only when ``mass_scenario == "biphasic"`` (the Tier-2
generative production scenario). It **hard-requires** the two partition fractions the
onset / S1 heat cannot be computed without -- ``internal_energy_partition_fraction``
(f_int) and ``internal_energy_retained_fraction`` (f_ret) -- plus (Slice-G review
fixes) a **drag_coefficients bundle** (biphasic is a drag-path scenario; a bundle-less
config would silently dispatch onto the hard-sphere collision path) and a
**non-negative** ``pickup_rate_coefficient`` (a negative lambda_0 gives P_attach < 0,
silently disabling pickup). It **warns (does not raise)** when
``pickup_rate_coefficient`` (lambda_0) or ``evap_rate_prefactor_per_ps`` (nu) is
exactly ``0.0`` (the channel is structurally inert; the evaporation-only-from-n0 /
pickup-only-growth limits are legitimate diagnostic runs). f_int/f_ret are still
range-checked by ``check_internal_energy_budget_config``; this guard only enforces
their *presence* under ``biphasic``.

``_biphasic_cfg`` here is the shared "minimal valid biphasic SimConfig" helper; the
step-level suite (``test_biphasic_step.py``) imports it rather than keeping a drifting
copy (CLAUDE.md principle 1).
"""

from __future__ import annotations

import warnings

import pytest

from i2_helium_md.config import SimConfig, check_biphasic_config
from i2_helium_md.physics.drag import LINEAR_CUBIC, DragCoefficients


def _drag_bundle() -> DragCoefficients:
    """A constant-mass linear_cubic bundle satisfying the biphasic bundle guard.

    Values mirror ``test_drag_config._constant_coeffs`` (the Tier-0 extraction
    numbers); the §6.5 constant-vs-evolving pairing tension is bypassed by
    ``allow_inconsistent_mass_pairing=True`` in ``_biphasic_cfg`` (the §6.6
    mid-window defence, same as the production pairing).
    """
    return DragCoefficients(
        form=LINEAR_CUBIC,
        coefficients={"a": 13.86, "b": 2.58},
        extraction_mass_model="constant",
        extraction_mass_amu=202.953908,  # the frozen 9 A extraction mass (m_eff)
        effective_binding_energy_I_ion_eV=0.3,
    )


def _biphasic_cfg(**overrides):
    """The minimal valid biphasic SimConfig (shared with test_biphasic_step)."""
    base = dict(
        mass_scenario="biphasic",
        internal_energy_partition_fraction=0.3,
        internal_energy_retained_fraction=0.2,
        pickup_rate_coefficient=0.5,
        coulomb_available_eV=0.80,
        drag_coefficients=_drag_bundle(),
        allow_inconsistent_mass_pairing=True,  # §6.6 mid-window; constant coeffs
    )
    base.update(overrides)
    return SimConfig(**base)


class TestNonBiphasicIsNoop:
    @pytest.mark.parametrize("scenario", ["fixed", "anchored_discrete"])
    def test_non_biphasic_scenarios_skip_the_guard(self, scenario):
        # f_int/f_ret None is fine off the biphasic path.
        cfg = SimConfig(mass_scenario=scenario)
        check_biphasic_config(cfg)  # does not raise

    def test_default_config_is_noop(self):
        check_biphasic_config(SimConfig())  # fixed default


class TestRequiredPartitionFractions:
    def test_valid_biphasic_passes(self):
        check_biphasic_config(_biphasic_cfg())

    def test_missing_f_int_rejected(self):
        with pytest.raises(ValueError, match="internal_energy_partition_fraction"):
            check_biphasic_config(_biphasic_cfg(internal_energy_partition_fraction=None))

    def test_missing_f_ret_rejected(self):
        with pytest.raises(ValueError, match="internal_energy_retained_fraction"):
            check_biphasic_config(_biphasic_cfg(internal_energy_retained_fraction=None))

    def test_reject_propagates_through_validate(self):
        with pytest.raises(ValueError, match="internal_energy_partition_fraction"):
            _biphasic_cfg(internal_energy_partition_fraction=None).validate()


class TestRequiredDragBundle:
    """Slice-G review fix: biphasic is a drag-path scenario -- a bundle-less config
    would silently dispatch onto the hard-sphere collision path (E_int frozen, the
    E_pot binding fold lost after column 0), so it is refused at config-load."""

    def test_missing_drag_bundle_rejected(self):
        with pytest.raises(ValueError, match="drag_coefficients"):
            check_biphasic_config(_biphasic_cfg(drag_coefficients=None))

    def test_missing_drag_bundle_rejected_through_validate(self):
        with pytest.raises(ValueError, match="drag_coefficients"):
            _biphasic_cfg(drag_coefficients=None).validate()

    def test_non_biphasic_without_bundle_unaffected(self):
        # The collision path stays bundle-less as before.
        check_biphasic_config(SimConfig(mass_scenario="fixed"))


class TestNegativeLambda0Rejected:
    """Slice-G review fix: lambda_0 < 0 gives P_attach < 0 (channel silently never
    fires), so the sign typo is refused; 0.0 stays advisory (see TestPickupInertWarns)."""

    def test_negative_lambda0_rejected(self):
        with pytest.raises(ValueError, match="pickup_rate_coefficient"):
            check_biphasic_config(_biphasic_cfg(pickup_rate_coefficient=-0.5))


class TestNegativeNuRejected:
    """Post-review fix (2026-07-02, Phase-B review): nu < 0 gives k < 0 -> P_shed < 0,
    so the evaporation channel silently never fires -- the same silent-shut-off class
    as the lambda_0 sign guard above. The Slice-Q decision that nu carries no load-time
    bound covers the *value* (Sourced/pinned 2.42), not the sign;
    ``physics.evaporation.rrk_rate`` carries the mirroring module-level defense."""

    def test_negative_nu_rejected(self):
        with pytest.raises(ValueError, match="evap_rate_prefactor_per_ps"):
            check_biphasic_config(_biphasic_cfg(evap_rate_prefactor_per_ps=-2.42))

    def test_negative_nu_rejected_through_validate(self):
        with pytest.raises(ValueError, match="evap_rate_prefactor_per_ps"):
            _biphasic_cfg(evap_rate_prefactor_per_ps=-2.42).validate()

    def test_default_nu_passes(self):
        check_biphasic_config(_biphasic_cfg())  # NU_EVAP_PER_PS = 2.42 default


class TestPickupInertWarns:
    def test_lambda0_zero_warns_not_raises(self):
        cfg = _biphasic_cfg(pickup_rate_coefficient=0.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            check_biphasic_config(cfg)  # must not raise
        assert any(
            issubclass(w.category, UserWarning)
            and "pickup_rate_coefficient" in str(w.message)
            for w in caught
        )

    def test_lambda0_positive_does_not_warn(self):
        cfg = _biphasic_cfg(pickup_rate_coefficient=0.5)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            check_biphasic_config(cfg)
        assert not any(
            "pickup_rate_coefficient" in str(w.message) for w in caught
        )


class TestEvapInertWarns:
    """Slice-G re-review fix (2026-07-02): nu == 0.0 mirrors the lambda_0 == 0.0
    advisory -- the evaporation channel is structurally inert (k = nu * (...) = 0,
    P_shed = 0), so the run is pickup-only growth from n_0. Legitimate as a
    diagnostic, but nu defaults to the Sourced 2.42, so an explicit 0 deserves the
    same warn-not-raise treatment as the pickup side (no silent asymmetry)."""

    def test_nu_zero_warns_not_raises(self):
        cfg = _biphasic_cfg(evap_rate_prefactor_per_ps=0.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            check_biphasic_config(cfg)  # must not raise
        assert any(
            issubclass(w.category, UserWarning)
            and "evap_rate_prefactor_per_ps" in str(w.message)
            for w in caught
        )

    def test_default_nu_does_not_warn(self):
        cfg = _biphasic_cfg()  # NU_EVAP_PER_PS = 2.42 default
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            check_biphasic_config(cfg)
        assert not any(
            "evap_rate_prefactor_per_ps" in str(w.message) for w in caught
        )
