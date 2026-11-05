"""Unit coverage for scripts/tier0_common.py (Method-B Tier-0 run wiring).

No I/O beyond reading the committed fit_parameters.json bundles and no
propagation: every test builds a ``SimConfig`` and inspects it. The expensive
neutral/ion runs are exercised by the scripts themselves, not here.
"""

from __future__ import annotations

import itertools

import pytest

from i2_helium_md import DRAG_BUNDLE_CATALOG, load_capped_cc, resolve_drag_bundle_dir
from scripts.tier0_common import (
    M_EFF_AMU,
    TIER0_WINDOWS_PS,
    build_drag_cfg,
    build_tier0_cfg,
    run_dir_name,
)

_CASES = ("9A", "18A")
_BUILD_KW = dict(num_molecules=2, ion_time_ps=0.02, dt_ion_ps=0.01, seed=123)


@pytest.mark.parametrize(
    "case,variant", list(itertools.product(_CASES, sorted(DRAG_BUNDLE_CATALOG)))
)
def test_every_catalog_variant_builds_and_validates(case, variant):
    """Each catalog bundle wires onto either geometry and validates.

    Pins the section-6.5.1 pairing identity: the cfg's binding equals the
    bundle's stamped effective binding, so ``validate()`` passes with no escape
    hatch, and the cfg form matches the bundle's stamped form.
    """
    cfg = build_drag_cfg(case, variant, **_BUILD_KW)
    cfg.validate()  # raises on any guard failure

    assert cfg.drag_coefficients is not None
    assert cfg.drag_form == cfg.drag_coefficients.form
    assert cfg.mass_scenario == "fixed"
    assert cfg.m_eff_amu == M_EFF_AMU
    assert cfg.drag_coefficients.extraction_mass_amu == M_EFF_AMU
    # The coupled-pair identity (binding wired straight from the bundle).
    assert (
        cfg.binding_energy_I_ion_eV
        == cfg.drag_coefficients.effective_binding_energy_I_ion_eV
    )


def test_default_build_matches_committed_bundle_coefficients():
    """With no overrides the cfg carries the committed coefficients verbatim."""
    from i2_helium_md import load_drag_coefficients

    coeffs = load_drag_coefficients(
        resolve_drag_bundle_dir("CC"),
        expected_m_eff_amu=M_EFF_AMU,
    )
    cfg = build_drag_cfg("18A", "CC", **_BUILD_KW)
    assert cfg.drag_coefficients.coefficients == coeffs.coefficients
    assert (
        cfg.drag_coefficients.effective_binding_energy_I_ion_eV
        == coeffs.effective_binding_energy_I_ion_eV
    )


def test_coeff_override_lands_in_cfg():
    """A coefficient override replaces just that key; the rest stay committed."""
    cfg = build_drag_cfg(
        "18A", "CC", coeff_overrides={"b": 9.99}, **_BUILD_KW
    )
    assert cfg.drag_form == "linear_cubic"
    assert cfg.drag_coefficients.coefficients["b"] == 9.99
    # a is untouched (still the committed a = 0 of CC).
    assert cfg.drag_coefficients.coefficients["a"] == 0.0
    cfg.validate()


def test_e_bind_override_ties_both_fields_and_validates():
    """An E_bind override is wired into both the cfg and the bundle identity."""
    cfg = build_drag_cfg(
        "9A", "CC", e_bind_override=0.222, **_BUILD_KW
    )
    assert cfg.binding_energy_I_ion_eV == 0.222
    assert cfg.drag_coefficients.effective_binding_energy_I_ion_eV == 0.222
    cfg.validate()  # identity holds, no escape hatch needed


_LQ = "LQ"
_M21_AMU = 126.90 + 21 * 4.0026   # I+ + 21 He, constants.py m(n) convention


def test_m_eff_override_sets_mass_and_stamp_and_validates():
    """G11.2: a full lq override at m = 21 sets cfg mass and the bundle stamp."""
    cfg = build_drag_cfg("9A", _LQ, coeff_overrides={"a": 18.0, "c": 1.6},
                         e_bind_override=0.117, m_eff_amu=_M21_AMU, **_BUILD_KW)
    assert cfg.m_eff_amu == cfg.mass_initial_amu == _M21_AMU
    assert cfg.drag_coefficients.extraction_mass_amu == _M21_AMU
    cfg.validate()   # the §6.5 mass guard accepts the restamped pair


def test_m_eff_override_at_reference_mass_is_identity():
    """The mass parameter at M_EFF_AMU builds the same cfg as no parameter."""
    kw = dict(coeff_overrides={"a": 18.0, "c": 1.6}, e_bind_override=0.117, **_BUILD_KW)
    assert (build_drag_cfg("9A", _LQ, m_eff_amu=M_EFF_AMU, **kw)
            == build_drag_cfg("9A", _LQ, **kw))


def test_m_eff_override_without_coeff_overrides_runs_committed_law():
    """The committed bundle law at m = 21: coefficients verbatim, run cfg restamped."""
    plain = build_drag_cfg("18A", _LQ, **_BUILD_KW)
    heavy = build_drag_cfg("18A", _LQ, m_eff_amu=_M21_AMU, **_BUILD_KW)
    assert heavy.drag_coefficients.coefficients == plain.drag_coefficients.coefficients
    assert heavy.binding_energy_I_ion_eV == plain.binding_energy_I_ion_eV
    assert heavy.m_eff_amu == heavy.drag_coefficients.extraction_mass_amu == _M21_AMU
    heavy.validate()


def test_m_eff_override_must_be_positive():
    with pytest.raises(ValueError, match="must be positive"):
        build_drag_cfg("9A", _LQ, m_eff_amu=0.0, **_BUILD_KW)


def test_unknown_override_key_for_form_raises():
    """An override key that is not a coefficient of the form is refused."""
    with pytest.raises(ValueError, match="not coefficients of form"):
        # 'c' is a linear_quadratic key, not valid for linear_cubic.
        build_drag_cfg(
            "18A", "CC", coeff_overrides={"c": 1.0}, **_BUILD_KW
        )


def test_unknown_case_and_variant_raise():
    with pytest.raises(ValueError, match="CASE must be one of"):
        build_drag_cfg("12A", "CC", **_BUILD_KW)
    with pytest.raises(ValueError, match="VARIANT must be one of"):
        build_drag_cfg("9A", "nope", **_BUILD_KW)


def test_tier0_lq_laws_are_the_bundles():
    for law in ("LQ", "LQ_eff"):
        assert build_tier0_cfg("18A", law, **_BUILD_KW) == build_drag_cfg("18A", law, **_BUILD_KW)


def test_tier0_cc_runs_capped_at_the_3p5_gate():
    cfg = build_tier0_cfg("9A", "CC", m_eff_amu=_M21_AMU, **_BUILD_KW)
    assert cfg.drag_form == "capped_cubic"
    assert cfg.drag_coefficients.coefficients == load_capped_cc().coefficients
    assert cfg.drag_coefficients.extraction_mass_amu == _M21_AMU
    assert (cfg.drag_spatial_gate, cfg.drag_gate_steepness) == ("erf_independent", 3.5)
    assert cfg.allow_inconsistent_gate_pairing


def test_tier0_ridge_point_is_the_committed_m21_fit():
    cfg = build_tier0_cfg("18A", "LQ", ridge_c=1.6, m_eff_amu=_M21_AMU, **_BUILD_KW)
    assert cfg.drag_coefficients.coefficients == {"a": 17.75706787001758, "c": 1.6}
    assert cfg.binding_energy_I_ion_eV == 0.11766211392031967
    override = build_tier0_cfg("18A", "LQ", ridge_c=1.6, e_bind_override=0.3,
                               m_eff_amu=_M21_AMU, **_BUILD_KW)
    assert override.binding_energy_I_ion_eV == 0.3


def test_tier0_ridge_needs_lq_and_an_existing_point():
    with pytest.raises(ValueError, match="ridge_c needs"):
        build_tier0_cfg("18A", "CC", ridge_c=1.6, **_BUILD_KW)
    with pytest.raises(FileNotFoundError):
        build_tier0_cfg("18A", "LQ", ridge_c=1.8, **_BUILD_KW)


def test_run_dir_name_keeps_the_stored_id_of_renamed_bundles():
    """Run names keep the id on disk, so generators still find finished runs."""
    assert (
        run_dir_name("9A", "CC", 50, run_tag="tuned")
        == "9A_drag_shared_pure_cubic_N50_tuned"
    )
    assert (run_dir_name("18A", "LQ", 50, run_tag="B2m21")
            == "18A_drag_srho3p5A18A_linear_quadratic_N50_B2m21")
    assert run_dir_name("18A", "LQ_eff", 50) == "18A_drag_LQ_eff_N50"


def test_tier0_windows_are_the_extraction_windows():
    """The values the per-case fit_parameters.json files stamped (spec section 3.1)."""
    assert TIER0_WINDOWS_PS == {"9A": (2.67, 14.081), "18A": (4.54, 14.768)}
