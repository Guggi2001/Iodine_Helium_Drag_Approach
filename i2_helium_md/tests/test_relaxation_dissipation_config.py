"""Config surface for the E2 Landau-gated drag arm (plan §I.11.2 item 2, arm (c)).

The relaxation stage's dissipation is a byte-inert enum ``relaxation_dissipation``:
``zero_gamma`` (default -- the delivered conservative closure) vs
``landau_gated_drag`` (gamma = 0 below the Landau cutoff ``v_limit``, the locked
pure-cubic drag above). Guarded at config-load like every other string-enum
surface (``_reject_unknown_enum``), with the T5/T7 no-silent-inert rule: the
``free_flight`` translation arm never reads a dissipation coefficient, so pairing
it with ``landau_gated_drag`` is refused rather than silently ignored.
"""

from dataclasses import replace

import pytest

from i2_helium_md.config import SimConfig, check_relaxation_config

from tests.test_relaxation_stage import _relax_cfg


def test_default_is_zero_gamma():
    assert _relax_cfg().relaxation_dissipation == "zero_gamma"


def test_field_default_backs_missing_json_key():
    # An old cfg.json predating the field omits the key; RunDirectory.load_cfg
    # reconstructs via SimConfig(**payload), so the dataclass default is the
    # back-compat guarantee (the internal_energy_partition_law precedent).
    field = SimConfig.__dataclass_fields__["relaxation_dissipation"]
    assert field.default == "zero_gamma"


def test_landau_coulomb_biphasic_accepted():
    cfg = _relax_cfg(relaxation_dissipation="landau_gated_drag",
                     relaxation_forces="coulomb")
    check_relaxation_config(cfg)  # no raise


def test_unknown_dissipation_rejected():
    cfg = _relax_cfg(relaxation_dissipation="bogus")
    with pytest.raises(ValueError, match="relaxation_dissipation"):
        check_relaxation_config(cfg)


def test_landau_free_flight_rejected():
    # No-silent-inert: free_flight has no drag O-step to read the coefficient.
    cfg = _relax_cfg(relaxation_dissipation="landau_gated_drag",
                     relaxation_forces="free_flight")
    with pytest.raises(ValueError, match="landau_gated_drag"):
        check_relaxation_config(cfg)


def test_disabled_stage_ignores_dissipation():
    # Stage off -> check_relaxation_config is a no-op (default-scope preserved),
    # even under an otherwise-rejected pairing.
    cfg = replace(
        _relax_cfg(relaxation_dissipation="landau_gated_drag",
                   relaxation_forces="free_flight"),
        relaxation_stage_enabled=False,
    )
    check_relaxation_config(cfg)  # no raise
