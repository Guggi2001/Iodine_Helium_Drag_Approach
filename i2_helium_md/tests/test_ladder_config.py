"""Tier-2 Phase-A Slice L -- SimConfig ladder surface + config-load guard.

Config-construction + guard tests only (no run, no trajectory), mirroring
``test_drag_config.py``. The ladder guard is a load-time fail-loud check (it does
not silently clamp): it rejects an unknown electronic picture or ladder selector
(``Literal`` is not runtime-enforced -- this is the recovery for that, mirroring
the ``mass_scenario`` / ``drag_form`` reject arms), and asserts the resolved first
rung ``D_0(1)`` exceeds the bulk-He floor ``D_floor`` for the configured picture.
"""

import typing

import pytest

from i2_helium_md import SimConfig
from i2_helium_md.config import (
    LadderElectronicPicture,
    LadderForm,
    check_ladder_config,
)
from i2_helium_md.physics.dissociation_ladder import _FIRST_RUNG_EV


def test_default_config_ladder_surface_is_inert_and_valid():
    cfg = SimConfig()
    assert cfg.dissociation_ladder == "form_u"
    assert cfg.ladder_electronic_picture == "statistical_mixture"
    assert cfg.ladder_steepness > 0.0
    check_ladder_config(cfg)        # does not raise
    cfg.validate()                  # full validate path does not raise


@pytest.mark.parametrize("picture", ["statistical_mixture", "x2_only", "cooling_relaxed"])
def test_all_valid_pictures_pass_guard(picture):
    cfg = SimConfig(ladder_electronic_picture=picture)
    check_ladder_config(cfg)        # D_0(1) > D_floor holds for every picture


def test_unknown_picture_rejected_by_guard():
    cfg = SimConfig(ladder_electronic_picture="nonsense")
    with pytest.raises(ValueError, match="picture"):
        check_ladder_config(cfg)


def test_unknown_picture_rejected_via_validate():
    cfg = SimConfig(ladder_electronic_picture="nonsense")
    with pytest.raises(ValueError, match="picture"):
        cfg.validate()


def test_unknown_ladder_form_rejected_by_guard():
    cfg = SimConfig(dissociation_ladder="quadratic_typo")
    with pytest.raises(ValueError, match="dissociation_ladder"):
        check_ladder_config(cfg)


def test_picture_literal_members_match_module_source_of_truth():
    # The config Literal must list exactly the module's picture keys (principle 1:
    # single source of truth for the picture keying).
    literal_members = set(typing.get_args(LadderElectronicPicture))
    assert literal_members == set(_FIRST_RUNG_EV)


def test_ladder_form_literal_members():
    assert set(typing.get_args(LadderForm)) == {"form_u", "tabulated"}


# --- ladder_steepness (kappa) positivity guard -----------------------------
@pytest.mark.parametrize("bad_kappa", [0.0, -1.0, -0.3])
def test_nonpositive_kappa_rejected(bad_kappa):
    # kappa <= 0 inverts / flattens the cliff and drives (1 - sigma(1)) -> 0
    # (inf rungs); must fail loud at config-load (CLAUDE.md principle 4).
    cfg = SimConfig(ladder_steepness=bad_kappa)
    with pytest.raises(ValueError, match="ladder_steepness"):
        check_ladder_config(cfg)


@pytest.mark.parametrize("good_kappa", [0.3, 1.0, 5.0])
def test_positive_kappa_accepted(good_kappa):
    check_ladder_config(SimConfig(ladder_steepness=good_kappa))   # no raise


def test_nonpositive_kappa_rejected_via_validate():
    with pytest.raises(ValueError, match="ladder_steepness"):
        SimConfig(ladder_steepness=-2.0).validate()


# --- Slice T2 (§I.10): tabulated_ladder_rungs_eV config surface -------------
# The config-load validation is single-sourced in
# physics.dissociation_ladder.resolve_ladder (called from check_ladder_config);
# these tests exercise it through the config guard.
import numpy as np

from i2_helium_md.physics.constants import N_STAR
from i2_helium_md.physics.dissociation_ladder import d0_of_n


def _form_u_rungs(length=N_STAR + 11):
    return tuple(
        np.atleast_1d(d0_of_n(np.arange(1, length + 1), kappa=1.0))
    )


class TestTabulatedLadderConfig:
    def test_default_rungs_none_and_inert(self):
        cfg = SimConfig()
        assert cfg.tabulated_ladder_rungs_eV is None
        cfg.validate()                  # inert-default regression

    def test_tabulated_with_valid_rungs_passes(self):
        cfg = SimConfig(
            dissociation_ladder="tabulated",
            tabulated_ladder_rungs_eV=_form_u_rungs(),
        )
        check_ladder_config(cfg)        # does not raise
        cfg.validate()

    def test_tabulated_without_rungs_rejected(self):
        cfg = SimConfig(dissociation_ladder="tabulated")
        with pytest.raises(ValueError, match="tabulated_ladder_rungs_eV"):
            check_ladder_config(cfg)

    def test_form_u_with_rungs_rejected(self):
        # Off-diagonal pairing refused loudly (mirrors the drag form/coefficient
        # cross-check style; a silently-ignored table is a stale-intent hazard).
        cfg = SimConfig(tabulated_ladder_rungs_eV=_form_u_rungs())
        with pytest.raises(ValueError, match="form_u"):
            check_ladder_config(cfg)

    def test_short_table_rejected(self):
        cfg = SimConfig(
            dissociation_ladder="tabulated",
            tabulated_ladder_rungs_eV=_form_u_rungs(length=N_STAR - 1),
        )
        with pytest.raises(ValueError, match=str(N_STAR)):
            check_ladder_config(cfg)

    def test_nonpositive_rung_rejected(self):
        rungs = list(_form_u_rungs())
        rungs[0] = 0.0
        cfg = SimConfig(
            dissociation_ladder="tabulated",
            tabulated_ladder_rungs_eV=tuple(rungs),
        )
        with pytest.raises(ValueError, match="positive"):
            check_ladder_config(cfg)

    def test_tabulated_without_rungs_rejected_via_validate(self):
        with pytest.raises(ValueError, match="tabulated_ladder_rungs_eV"):
            SimConfig(dissociation_ladder="tabulated").validate()
