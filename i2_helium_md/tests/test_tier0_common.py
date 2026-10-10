"""Unit coverage for scripts/tier0_common.py (Method-B Tier-0 run wiring).

No I/O beyond reading the committed fit_parameters.json bundles and no
propagation: every test builds a ``SimConfig`` and inspects it. The expensive
neutral/ion runs are exercised by the scripts themselves, not here.
"""

from __future__ import annotations

import itertools

import pytest

from i2_helium_md.config import _UNSTAMPED_LEGACY_DRAG_BUNDLES
from scripts.tier0_common import (
    CATALOG,
    M_EFF_AMU,
    build_drag_cfg,
    resolve_bundle_dir,
    run_dir_name,
    window_source_dir,
)

_CASES = ("9A", "18A")
_BUILD_KW = dict(num_molecules=2, ion_time_ps=0.02, dt_ion_ps=0.01, seed=123)


@pytest.mark.parametrize(
    "case,variant", list(itertools.product(_CASES, sorted(CATALOG)))
)
def test_every_catalog_variant_builds_and_validates(case, variant):
    """Each catalog bundle wires onto either geometry and validates.

    Pins the section-6.5.1 pairing identity: the cfg's binding equals the
    bundle's stamped effective binding, so ``validate()`` passes with no escape
    hatch, and the cfg form matches the bundle's stamped form. A G4 catalog id
    whose bundle has not been produced yet (or whose fit trapped and wrote no
    loadable bundle) is skipped, not failed.
    """
    if (variant not in _UNSTAMPED_LEGACY_DRAG_BUNDLES
            and not (resolve_bundle_dir(case, variant) / "fit_parameters.json").is_file()):
        pytest.skip(f"bundle {variant!r} not produced (G4 refit pending/trapped)")
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
        resolve_bundle_dir("18A", "shared_pure_cubic"),
        expected_m_eff_amu=M_EFF_AMU,
    )
    cfg = build_drag_cfg("18A", "shared_pure_cubic", **_BUILD_KW)
    assert cfg.drag_coefficients.coefficients == coeffs.coefficients
    assert (
        cfg.drag_coefficients.effective_binding_energy_I_ion_eV
        == coeffs.effective_binding_energy_I_ion_eV
    )


def test_coeff_override_lands_in_cfg():
    """A coefficient override replaces just that key; the rest stay committed."""
    cfg = build_drag_cfg(
        "18A", "shared_pure_cubic", coeff_overrides={"b": 9.99}, **_BUILD_KW
    )
    assert cfg.drag_form == "linear_cubic"
    assert cfg.drag_coefficients.coefficients["b"] == 9.99
    # a is untouched (still the committed a = 0 of shared_pure_cubic).
    assert cfg.drag_coefficients.coefficients["a"] == 0.0
    cfg.validate()


def test_e_bind_override_ties_both_fields_and_validates():
    """An E_bind override is wired into both the cfg and the bundle identity."""
    cfg = build_drag_cfg(
        "9A", "shared_pure_cubic", e_bind_override=0.222, **_BUILD_KW
    )
    assert cfg.binding_energy_I_ion_eV == 0.222
    assert cfg.drag_coefficients.effective_binding_energy_I_ion_eV == 0.222
    cfg.validate()  # identity holds, no escape hatch needed


_LQ = "srho3p5A18A_linear_quadratic"
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
            "18A", "shared_pure_cubic", coeff_overrides={"c": 1.0}, **_BUILD_KW
        )


def test_unknown_case_and_variant_raise():
    with pytest.raises(ValueError, match="CASE must be one of"):
        resolve_bundle_dir("12A", "shared_pure_cubic")
    with pytest.raises(ValueError, match="VARIANT must be one of"):
        resolve_bundle_dir("9A", "nope")


def test_per_case_variant_resolves_to_own_geometry():
    """The geometry guard: a per-case key only ever points at its own case."""
    assert resolve_bundle_dir("9A", "percase_pl").parts[-3] == "9A"
    assert resolve_bundle_dir("18A", "percase_pl").parts[-3] == "18A"


def test_run_dir_name_and_window_source():
    assert run_dir_name("18A", "shared_pl", 50) == "18A_drag_shared_pl_N50"
    assert (
        run_dir_name("9A", "shared_pure_cubic", 50, run_tag="tuned")
        == "9A_drag_shared_pure_cubic_N50_tuned"
    )
    assert window_source_dir("9A").parts[-2:] == ("9A", "trajectory_matching")
