"""The three thesis drag bundles and the ids finished runs keep in their cfg.json.

``data/runs`` is not committed, so the finished-run contract is pinned on a
minimal stored ``cfg.json`` carrying the historical ids verbatim. Every comparison
here is exact: all values are read back from JSON, nothing is computed.
"""

from __future__ import annotations

import dataclasses
import json

import pytest

from i2_helium_md import (
    DRAG_BUNDLE_CATALOG,
    SimConfig,
    load_capped_cc,
    load_drag_bundle,
)
from i2_helium_md.config import HISTORICAL_DRAG_BUNDLE_IDS, current_drag_bundle_id
from i2_helium_md.physics.drag import DragCoefficients
from i2_helium_md.simulation.run_directory import RunDirectory
from scripts.gen_tier2_detfix_postpickup import B1W_C, B1W_E_BIND_EV, B3_A
from scripts.gen_tier2_g413_ridge_battery import ridge_point
from scripts.tier2_common import cfg_diff_vs_reference, with_current_bundle_ids

_M_EFF = 202.953908
_OLD_LQ = "srho3p5A18A_linear_quadratic"
_OLD_CC = "shared_pure_cubic"


def _stored_coeffs(bundle_id):
    """As a finished run stores them (the B2 / h405 cfg.json payload shape)."""
    if bundle_id is None:
        return DragCoefficients(
            form="capped_cubic",
            coefficients={"b": 2.5153508541760052, "v_c": 5.4, "p_tail": -1.0},
            extraction_mass_model="constant", extraction_mass_amu=_M_EFF,
            extraction_method="trajectory_matching",
            effective_binding_energy_I_ion_eV=0.11675778353879479)
    return DragCoefficients(
        form="linear_quadratic", coefficients={"a": 18.007189923460135, "c": 1.6},
        extraction_mass_model="constant", extraction_mass_amu=_M_EFF,
        extraction_method="trajectory_matching",
        effective_binding_energy_I_ion_eV=0.11727328616561973,
        extraction_gate_steepness_angstrom=3.5, bundle_id=bundle_id)


def _finished_run(tmp_path, drag_bundle, bundle_id):
    """A run directory whose cfg.json carries ``drag_bundle`` / ``bundle_id`` verbatim."""
    coeffs = _stored_coeffs(bundle_id)
    gate = ({} if bundle_id is None else
            dict(drag_spatial_gate="erf_independent", drag_gate_steepness=3.5))
    cfg = SimConfig(
        drag_form=coeffs.form, drag_coefficients=coeffs, drag_bundle=drag_bundle,
        binding_energy_I_ion_eV=coeffs.effective_binding_energy_I_ion_eV, **gate)
    run = RunDirectory(tmp_path / "run")
    run.save_cfg(cfg)
    return run, cfg


class TestFinishedRuns:
    @pytest.mark.parametrize("drag_bundle, bundle_id",
                             [(_OLD_LQ, _OLD_LQ), (_OLD_CC, None)])
    def test_stored_ids_load_validate_and_stay_verbatim(self, tmp_path, drag_bundle, bundle_id):
        run, cfg = _finished_run(tmp_path, drag_bundle, bundle_id)
        loaded = run.load_cfg()
        loaded.validate()
        assert loaded == cfg
        assert loaded.drag_bundle == drag_bundle
        assert loaded.drag_coefficients.bundle_id == bundle_id

    def test_a_run_without_the_key_reads_as_the_unstamped_default(self, tmp_path):
        """Pre-G4 cfg.json files have no ``drag_bundle`` key at all."""
        run, _ = _finished_run(tmp_path, _OLD_CC, None)
        payload = json.loads(run.cfg_path.read_text(encoding="utf-8"))
        del payload["drag_bundle"]
        run.cfg_path.write_text(json.dumps(payload), encoding="utf-8")
        loaded = run.load_cfg()
        loaded.validate()
        assert loaded.drag_bundle == "CC"

    def test_unknown_id_is_still_refused(self):
        with pytest.raises(ValueError, match="drag_bundle"):
            SimConfig(drag_bundle="srho3p5A18A_cubic").validate()

    @pytest.mark.parametrize("old", [_OLD_LQ, _OLD_CC])
    def test_historical_ids_cannot_select_a_bundle(self, old):
        with pytest.raises(ValueError, match="drag bundle must be one of"):
            load_drag_bundle(old, expected_m_eff_amu=_M_EFF)


class TestCurrentIds:
    def test_mapping(self):
        assert HISTORICAL_DRAG_BUNDLE_IDS == {"LQ": _OLD_LQ, "CC": _OLD_CC}
        assert current_drag_bundle_id(_OLD_LQ) == "LQ"
        assert current_drag_bundle_id(_OLD_CC) == "CC"
        assert current_drag_bundle_id("LQ_eff") == "LQ_eff"
        assert current_drag_bundle_id(None) is None

    @pytest.mark.parametrize("old, new, bundle_id",
                             [(_OLD_LQ, "LQ", _OLD_LQ), (_OLD_CC, "CC", None)])
    def test_stored_ids_are_no_cfg_difference(self, tmp_path, old, new, bundle_id):
        run, cfg = _finished_run(tmp_path, old, bundle_id)
        built = with_current_bundle_ids(cfg)
        assert built.drag_bundle == new
        assert cfg_diff_vs_reference(built, run.cfg_path) == []

    def test_a_real_difference_still_shows(self, tmp_path):
        run, cfg = _finished_run(tmp_path, _OLD_LQ, _OLD_LQ)
        built = dataclasses.replace(with_current_bundle_ids(cfg), seed=7)
        assert cfg_diff_vs_reference(built, run.cfg_path) == ["seed"]

    def test_with_current_bundle_ids_leaves_every_other_field(self, tmp_path):
        _, cfg = _finished_run(tmp_path, _OLD_LQ, _OLD_LQ)
        built = with_current_bundle_ids(cfg)
        assert dataclasses.replace(built, drag_bundle=_OLD_LQ, drag_coefficients=None) \
            == dataclasses.replace(cfg, drag_coefficients=None)
        assert built.drag_coefficients.coefficients == cfg.drag_coefficients.coefficients


class TestThesisLaws:
    def test_catalog(self):
        assert set(DRAG_BUNDLE_CATALOG) == {"LQ", "LQ_eff", "CC"}

    def test_lq_is_b2(self):
        c = load_drag_bundle("LQ", expected_m_eff_amu=_M_EFF)
        assert c.coefficients == {"a": 18.007189923460135, "c": 1.6}
        assert c.effective_binding_energy_I_ion_eV == 0.11727328616561973
        assert c.extraction_gate_steepness_angstrom == 3.5
        assert c.bundle_id == "LQ"

    def test_lq_eff_is_the_b3_cell_of_the_generator(self):
        c = load_drag_bundle("LQ_eff", expected_m_eff_amu=_M_EFF)
        assert c.coefficients == {"a": B3_A, "c": B1W_C}
        assert c.effective_binding_energy_I_ion_eV == B1W_E_BIND_EV
        assert c.extraction_gate_steepness_angstrom == 3.5
        assert c.bundle_id == "LQ_eff"

    def test_cc_bundle_stays_the_uncapped_preset_law(self):
        c = load_drag_bundle("CC", expected_m_eff_amu=_M_EFF)
        assert (c.form, c.coefficients) == (
            "linear_cubic", {"a": 0.0, "b": 2.5153508541760052})
        assert c.bundle_id is None and c.extraction_gate_steepness_angstrom is None

    def test_capped_cc_is_the_stored_h405_law(self):
        assert load_capped_cc() == _stored_coeffs(None)


class TestCommittedLaws:
    def test_b2_law_is_the_lq_bundle(self):
        lq = load_drag_bundle("LQ", expected_m_eff_amu=_M_EFF)
        assert ridge_point(1.6) == {
            "a": lq.coefficients["a"], "c": 1.6,
            "e_bind_eV": lq.effective_binding_energy_I_ion_eV}

    def test_a5_law(self):
        assert ridge_point(2.0) == {
            "a": 16.71720991124528, "c": 2.0, "e_bind_eV": 0.11747528461884696}

    def test_other_points_are_refused(self):
        with pytest.raises(ValueError, match="no committed law"):
            ridge_point(1.8)
