"""THESIS_FINAL_GOALS G4.11: cross-case drag-law test (no MD, no reference files).

1. ``fit_drag_form`` recovers exact coefficients for the two new forms
   (``quadratic``, ``linear_quadratic``) and leaves the G4.10 forms unchanged;
2. ``bulk_law_force`` is the integrator's ungated law (same arithmetic as
   ``physics.drag``) for each family's bundle form;
3. ``cross_case_scale`` / ``shape_residual_rms`` on synthetic forces: the
   exact law gives s = 1 and shape 0; a scaled law gives the scale;
4. the pre-registered decision rule.
"""

import numpy as np
import pytest

from i2_helium_md.extraction.force_balance import (
    bulk_law_force,
    cross_case_scale,
    fit_drag_form,
    shape_residual_rms,
)
from i2_helium_md.physics.drag import LINEAR_CUBIC, LINEAR_QUADRATIC, DragCoefficients
from scripts.extraction.method_a_g411_cross_case import decide

V = np.linspace(2.7, 4.9, 60)   # the joint 18 A + 9 A speed span [A/ps]


def _coeffs(form: str, **c) -> DragCoefficients:
    return DragCoefficients(
        form=form, coefficients=c, extraction_mass_model="constant",
        extraction_mass_amu=202.953908, extraction_method="trajectory_matching",
        effective_binding_energy_I_ion_eV=0.12,
    )


class TestFitDragFormQuadratic:
    def test_quadratic_exact_recovery(self):
        fit = fit_drag_form(V, 7.2 * V**2, "quadratic")
        assert fit.c == pytest.approx(7.2, rel=1e-12)   # exact data: rounding only
        assert fit.a == 0.0 and fit.b == 0.0
        assert fit.rmse_force == pytest.approx(0.0, abs=1e-10)

    def test_linear_quadratic_exact_recovery(self):
        fit = fit_drag_form(V, 4.0 * V + 6.5 * V**2, "linear_quadratic")
        assert fit.a == pytest.approx(4.0, rel=1e-9)    # 2-column lstsq conditioning
        assert fit.c == pytest.approx(6.5, rel=1e-9)
        assert fit.b == 0.0

    def test_cubic_forms_unchanged_and_c_zero(self):
        fit = fit_drag_form(V, 3.0 * V + 2.2 * V**3, "linear_cubic")
        assert fit.a == pytest.approx(3.0, rel=1e-9)
        assert fit.b == pytest.approx(2.2, rel=1e-9)
        assert fit.c == 0.0 and fit.c_se == 0.0

    def test_unknown_form_raises(self):
        with pytest.raises(ValueError, match="form must be one of"):
            fit_drag_form(V, V, "power_law")


class TestBulkLawForce:
    def test_linear_quadratic_law(self):
        law = _coeffs(LINEAR_QUADRATIC, a=4.0, c=6.5)
        np.testing.assert_allclose(bulk_law_force(V, law), 4.0 * V + 6.5 * V**2,
                                   rtol=1e-15)

    def test_cubic_law_is_ungated(self):
        law = _coeffs(LINEAR_CUBIC, a=0.0, b=2.2286)
        np.testing.assert_allclose(bulk_law_force(V, law), 2.2286 * V**3, rtol=1e-15)


class TestCrossCaseScale:
    def test_exact_law_gives_one_and_zero_shape(self):
        L = 2.2 * V**3
        s = cross_case_scale(L.copy(), L)
        assert s == pytest.approx(1.0, rel=1e-14)
        assert shape_residual_rms(L, L, s) == pytest.approx(0.0, abs=1e-14)

    def test_scaled_law_gives_scale(self):
        L = 23.15 * V
        assert cross_case_scale(1.81 * L, L) == pytest.approx(1.81, rel=1e-14)

    def test_wrong_slope_shows_in_shape_not_only_scale(self):
        # Cubic forces scored with a linear law: some scale, nonzero shape.
        F, L = 2.2 * V**3, 23.15 * V
        s = cross_case_scale(F, L)
        assert shape_residual_rms(F, L, s) > 0.05

    def test_validation(self):
        with pytest.raises(ValueError, match="equal shape"):
            cross_case_scale(V, V[:-1])
        with pytest.raises(ValueError, match="vanishes"):
            cross_case_scale(V, np.zeros_like(V))
        with pytest.raises(ValueError, match="finite"):
            cross_case_scale(np.r_[V[:-1], np.nan], V)


class TestDecisionRule:
    def test_only_cubic(self):
        d = decide({"linear": False, "quadratic": False,
                    "linear_quadratic": False, "cubic": True})
        assert d["outcome"] == "cubic_supported"

    def test_several_consistent(self):
        d = decide({"linear": False, "quadratic": True,
                    "linear_quadratic": False, "cubic": True})
        assert d["outcome"] == "not_discriminated"
        assert d["consistent_families"] == ["cubic", "quadratic"]

    def test_cubic_not_consistent(self):
        d = decide({"linear": True, "quadratic": True,
                    "linear_quadratic": True, "cubic": False})
        assert d["outcome"] == "open_decision"


class TestTier0Selectable:
    """All five 18 A-only s_rho = 3.5 A fits are Tier-0 VARIANT keys.

    Reads the committed bundles; ``build_drag_cfg`` must load each one, bring
    the 3.5 A gate it was fitted under, and pass the config guards.
    """

    @pytest.mark.parametrize("family, form", [
        ("linear", LINEAR_CUBIC),
        ("cubic", LINEAR_CUBIC),
        ("linear_cubic", LINEAR_CUBIC),
        ("quadratic", LINEAR_QUADRATIC),
        ("linear_quadratic", LINEAR_QUADRATIC),
    ])
    def test_variant_builds_with_its_gate(self, family, form):
        from scripts.tier0_common import build_drag_cfg

        cfg = build_drag_cfg("18A", f"srho3p5A18A_{family}", num_molecules=2,
                             ion_time_ps=1.0, dt_ion_ps=0.01, seed=1)
        assert cfg.drag_bundle == f"srho3p5A18A_{family}"
        assert cfg.drag_coefficients.form == form
        assert cfg.drag_spatial_gate == "erf_independent"
        assert cfg.drag_gate_steepness == 3.5

    def test_figure_label_names_the_bundle(self):
        from scripts.post_processing.tier0_drag_comparison import drag_bundle_label
        from scripts.tier0_common import build_drag_cfg

        cfg = build_drag_cfg("18A", "srho3p5A18A_quadratic", num_molecules=2,
                             ion_time_ps=1.0, dt_ion_ps=0.01, seed=1)
        label = drag_bundle_label(cfg)
        assert label.startswith("srho3p5A18A_quadratic: linear_quadratic a=0, c=7.224")
        assert "E_bind 0.1191 eV" in label and "gate 3.5 Å" in label
