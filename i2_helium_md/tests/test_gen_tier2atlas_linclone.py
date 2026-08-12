"""Focused tests for the h405-clone MD battery (plan §6.5).

Config-construction, guard and pure-helper tests only — **no run, no
trajectory** (the `test_gen_tier2atlas_linring.py` pattern). These lock the
frozen §6.5 design: the clone cell and its three seeds, the free-form
provenance posture inherited from the ring, the **CRN guard** (whose whole
content is that `seed` and `num_molecules` must NOT differ from the
committed partner), the seed-only relation between members, the LC-P1
committed-artifact oracle, and the scorer's three-way band verdict.
"""

from __future__ import annotations

import dataclasses
import warnings

import pytest

from scripts.gen_tier2atlas_linclone import (
    CLONE_CELL,
    CLONE_TWIN_COLS,
    CLONE_TWIN_ROW,
    CRN_MEMBER,
    MEMBER_SEEDS,
    PARTNER_DIFF_KEYS,
    build_member,
    cfg_field_diff,
    clone_run_dir_name,
    verify_clone_preregistration,
    verify_crn_pairing,
    verify_member,
    verify_seed_only_vs_crn_member,
)
from scripts.gen_tier2atlas_linring import EBIND_EV, N


def _build(member: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return build_member(member)


class TestFrozenDesign:
    def test_clone_cell_is_the_frozen_6_4_row(self):
        assert (CLONE_CELL.a, CLONE_CELL.v_c, CLONE_CELL.eb_tag,
                CLONE_CELL.tau_ps, CLONE_CELL.e0_eV) == (
            42.5, None, "eb0482", 6.4, 0.31)

    def test_three_seeds_with_the_ring_seed_as_the_crn_member(self):
        assert MEMBER_SEEDS == {"s1": 20260731, "s2": 20260812, "s3": 20260813}
        # s1 IS the §6.3 ring seed — that identity is what buys the free
        # CRN pair against the committed h405p partner.
        assert MEMBER_SEEDS[CRN_MEMBER] == 20260731

    def test_frozen_twin_row_shape(self):
        assert len(CLONE_TWIN_COLS) == len(CLONE_TWIN_ROW)
        row = dict(zip(CLONE_TWIN_COLS, CLONE_TWIN_ROW))
        assert row["phi"] == "0.985" and row["gate"] == "1"

    def test_namespace(self):
        for member in MEMBER_SEEDS:
            name = clone_run_dir_name(member)
            assert "_tier2_" not in name and "tier2probe" not in name
            assert f"N{N}" in name and name.endswith(f"linclone{member}")


class TestMemberBuild:
    @pytest.mark.parametrize("member", sorted(MEMBER_SEEDS))
    def test_seed_applied_and_pins_kept(self, member):
        cfg = _build(member)
        assert cfg.seed == MEMBER_SEEDS[member]
        assert cfg.num_molecules == N
        assert cfg.mass_scenario == "biphasic"
        assert cfg.internal_energy_cooling_tau_ps == CLONE_CELL.tau_ps

    def test_free_form_posture_inherited_from_the_ring(self):
        cfg = _build(CRN_MEMBER)
        coeffs = cfg.drag_coefficients
        assert cfg.drag_form == "pure_linear"
        assert dict(coeffs.coefficients) == {"a": CLONE_CELL.a}
        assert coeffs.extraction_method == "free_form"
        assert coeffs.effective_binding_energy_I_ion_eV is None
        assert cfg.allow_unvalidated_binding_pairing is True
        assert cfg.binding_energy_I_ion_eV == EBIND_EV[CLONE_CELL.eb_tag]

    def test_members_differ_from_the_crn_member_in_exactly_the_seed(self):
        for member in MEMBER_SEEDS:
            cfg = _build(member)
            diff = cfg_field_diff(cfg, _build(CRN_MEMBER))
            assert diff == ([] if member == CRN_MEMBER else ["seed"])
            verify_seed_only_vs_crn_member(cfg, member)

    def test_a_drifted_pin_is_caught(self):
        cfg = dataclasses.replace(_build("s2"),
                                  internal_energy_cooling_tau_ps=4.8)
        with pytest.raises(AssertionError, match=r"exactly \['seed'\]"):
            verify_seed_only_vs_crn_member(cfg, "s2")

    def test_unknown_member_rejected(self):
        with pytest.raises(KeyError):
            build_member("s9")


class TestCrnGuard:
    """The pairing IS the guard (plan §6.5 oracle 2)."""

    @pytest.mark.parametrize("member", sorted(MEMBER_SEEDS))
    def test_diff_vs_partner_is_the_pre_registered_set(self, member):
        cfg = _build(member)
        diff = set(verify_crn_pairing(cfg, member))
        expected = (PARTNER_DIFF_KEYS if member == CRN_MEMBER
                    else PARTNER_DIFF_KEYS | {"seed"})
        assert diff == set(expected)
        assert "num_molecules" not in diff

    def test_crn_member_may_not_differ_in_seed(self):
        # s2's cfg presented as the CRN member: the seed now differs from the
        # committed partner, which is exactly what destroys the pairing.
        cfg = _build("s2")
        with pytest.raises(AssertionError, match="CRN guard FAILED"):
            verify_crn_pairing(cfg, CRN_MEMBER)

    def test_num_molecules_may_never_differ(self):
        cfg = dataclasses.replace(_build(CRN_MEMBER), num_molecules=250)
        with pytest.raises(AssertionError, match="CRN guard FAILED"):
            verify_crn_pairing(cfg, CRN_MEMBER)

    def test_full_member_verification_passes(self):
        for member in MEMBER_SEEDS:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                cfg, diff = verify_member(member)
            assert cfg.seed == MEMBER_SEEDS[member]
            assert diff


class TestLcP1Oracle:
    def test_frozen_clone_row_reproduces_from_the_committed_csv(self, capsys):
        verify_clone_preregistration()
        out = capsys.readouterr().out
        assert "LC-P1 clone oracle PASSED" in out
        assert "LR-P1 twin oracle PASSED" in out


class TestBandVerdict:
    """CL-P1's three-way outcome (plan §6.5) — the N-cut consequence."""

    def test_pass_fail_and_marginal(self):
        from scripts.post_processing.tier2atlas_linclone_table import (
            band_verdict,
        )
        band = (0.19, 0.30)
        # Comfortably inside / outside.
        assert band_verdict(0.24, band, 0.008) == "pass"
        assert band_verdict(0.15, band, 0.008) == "fail"
        # Within 1 SE of an edge, from either side: neither pass nor fail.
        assert band_verdict(0.195, band, 0.008) == "gate-marginal"
        assert band_verdict(0.185, band, 0.008) == "gate-marginal"
        # Zero/absent SE degenerates to the plain in-band test.
        assert band_verdict(0.19, band, 0.0) == "gate-marginal"
        assert band_verdict(0.25, band, float("nan")) == "pass"

    def test_twin_row_maps_the_frozen_strings(self):
        from scripts.post_processing.tier2atlas_linclone_table import twin_row
        twin = twin_row()
        assert twin["twin_KE1"] == 0.6153
        assert twin["twin_midHot"] == 1.0144
        assert twin["twin_nbar"] == 4.576
        assert twin["twin_n1"] == 0.1966
        assert twin["twin_trap"] == 0.0476
