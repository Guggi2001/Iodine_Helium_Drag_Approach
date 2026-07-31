"""Focused tests for scripts/gen_tier2atlas_linring.py (plan §6.1 ring).

Config-construction + oracle tests only — **no run, no trajectory** (the
gen_tier2atlas_lqbattery test pattern). The MD itself is the instrument;
these lock the frozen design: cell matrix, CRN seed sharing, the
free-form provenance posture, the pre-registered cfg-diff sets, and the
LR-P1 committed-artifact oracle.
"""

from __future__ import annotations

import warnings

import pytest

from scripts.gen_tier2atlas_linring import (
    EBIND_EV,
    H405_ANCHOR,
    N,
    RING_MATRIX,
    SEED,
    TWIN_ROWS,
    build_cell,
    expected_diff_keys,
    ring_run_dir_name,
    verify_corrected_geometry,
    verify_twin_preregistration,
)


class TestRingMatrix:
    def test_eight_cells_one_partner(self):
        assert len(RING_MATRIX) == 8
        lin = [c for c in RING_MATRIX if c.a is not None]
        capped = [c for c in RING_MATRIX if c.a is None]
        assert len(lin) == 7
        assert [c.label for c in capped] == ["h405p"]

    def test_partner_is_h405_pins(self):
        p = next(c for c in RING_MATRIX if c.label == "h405p")
        assert (p.v_c, p.eb_tag, p.tau_ps, p.e0_eV) == (5.5, "eb1168", 4.4, 0.405)

    def test_lin_cells_are_the_frozen_sub_plateau_core(self):
        # Plan §6.1 table, frozen 2026-07-31.
        got = {(c.label, c.a, c.eb_tag, c.tau_ps, c.e0_eV)
               for c in RING_MATRIX if c.a is not None}
        assert got == {
            ("lr1", 27.5, "eb0482", 4.8, 0.35),
            ("lr2", 27.5, "eb1168", 4.8, 0.35),
            ("lr3", 30.0, "eb0482", 4.8, 0.36),
            ("lr4", 30.0, "eb1168", 4.8, 0.36),
            ("lr5", 32.5, "eb0482", 4.8, 0.36),
            ("lr6", 35.0, "eb0482", 4.8, 0.37),
            ("lr7", 27.5, "eb0482", 4.8, 0.37),
        }

    def test_frozen_twin_rows_cover_exactly_the_lin_cells(self):
        assert set(TWIN_ROWS) == {c.label for c in RING_MATRIX
                                  if c.a is not None}
        assert len(H405_ANCHOR) == 7

    def test_namespace(self):
        for cell in RING_MATRIX:
            name = ring_run_dir_name(cell.label)
            assert "_tier2_" not in name and "tier2probe" not in name
            assert f"N{N}" in name and name.endswith(f"linr{cell.label}")


class TestCellBuild:
    @pytest.mark.parametrize("cell", RING_MATRIX, ids=lambda c: c.label)
    def test_cell_builds_validates_and_shares_the_crn_seed(self, cell):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            cfg = build_cell(cell)
            verify_corrected_geometry(cfg, cell)
        assert cfg.seed == SEED
        assert cfg.num_molecules == N
        assert cfg.mass_scenario == "biphasic"
        assert cfg.allow_inconsistent_mass_pairing is True

    def test_lin_cells_carry_free_form_posture(self):
        # extraction_method free_form, NO binding stamp, hatch set, well
        # from the frozen table — the honest-provenance adjudication.
        for cell in RING_MATRIX:
            if cell.a is None:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                cfg = build_cell(cell)
            coeffs = cfg.drag_coefficients
            assert cfg.drag_form == "pure_linear"
            assert coeffs.form == "pure_linear"
            assert dict(coeffs.coefficients) == {"a": cell.a}
            assert coeffs.extraction_method == "free_form"
            assert coeffs.effective_binding_energy_I_ion_eV is None
            assert cfg.allow_unvalidated_binding_pairing is True
            assert cfg.binding_energy_I_ion_eV == EBIND_EV[cell.eb_tag]

    def test_partner_carries_the_bundle_pairing(self):
        p = next(c for c in RING_MATRIX if c.label == "h405p")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            cfg = build_cell(p)
        coeffs = cfg.drag_coefficients
        assert cfg.drag_form == "capped_cubic"
        assert float(coeffs.coefficients["v_c"]) == 5.5
        assert float(coeffs.coefficients["p_tail"]) == -1.0
        assert cfg.binding_energy_I_ion_eV == (
            coeffs.effective_binding_energy_I_ion_eV
        )
        assert cfg.allow_unvalidated_binding_pairing is False

    def test_lin_binding_warning_fires_loudly(self):
        # The §6.5.1 refuse->warn arm must actually fire on a lin cell —
        # the hatch is documented, not silent.
        cell = next(c for c in RING_MATRIX if c.a is not None)
        with pytest.warns(RuntimeWarning, match="free_form"):
            build_cell(cell)


class TestPreregisteredDiff:
    def test_expected_diff_key_sets(self):
        for cell in RING_MATRIX:
            keys = expected_diff_keys(cell)
            assert "drag_coefficients" in keys
            assert "internal_energy_cooling_tau_ps" in keys
            if cell.a is None:
                assert "drag_form" not in keys        # capped like reference
                assert "binding_energy_I_ion_eV" not in keys
            else:
                assert "drag_form" in keys
                assert "allow_unvalidated_binding_pairing" in keys
                assert (("binding_energy_I_ion_eV" in keys)
                        == (cell.eb_tag == "eb0482"))


class TestCommittedArtifactOracle:
    def test_lr_p1_passes_against_the_committed_csvs(self):
        # Gate-on-committed-artifacts rule: the frozen rows must reproduce
        # string-exact from the repo-committed CSVs on every machine.
        verify_twin_preregistration()
