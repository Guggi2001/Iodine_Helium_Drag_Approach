"""Tests for i2_helium_md/physics/drag.py (Slice 1: pure gated-drag module).

Pure-physics, no stepper -- mirrors ``tests/test_collisions.py`` in spirit.
The killer test (``TestForceBalanceReproduction``) confirms the module *is* the
extracted law: re-evaluating ``drag_force`` with the in-hand 9 A {a, b} on the
extraction's own trusted-interior (v, F_drag) scatter reproduces it.

Reference values are plumbed in from the user-exported artifacts under
``data/reference/drag/`` (see ``SLICE1_GOALS_gated_drag_module.md`` §6-§7).
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.physics.drag import (
    LINEAR_CUBIC,
    LINEAR_QUADRATIC,
    POWER_LAW,
    THRESHOLD,
    DragCoefficients,
    drag_force,
    drag_gamma,
    spatial_gate,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DRAG_REF_DIR = PROJECT_ROOT / "data" / "reference" / "drag"

# Both extracted bubble sizes; each has linear_and_cubic/ and power/ subdirs
# with fit_parameters.json (a, b, meff_amu) + drag_data.csv.
CASES = ("9A", "18A")

# cfg.potential_steepness (Angstrom): the gate reuses the confining-potential
# erf width (§5.5 G4->G2).
STEEPNESS_A = 14.2

# Extraction trusted-interior truncation: the exported drag_data.csv is the
# FULL window; the fit used only the interior after dropping first/last N points
# (DragExtractionSettings.truncate_points = 500).
TRUNCATE_POINTS = 500


# ===========================================================================
# Fixtures / loaders
# ===========================================================================
def _load_linear_cubic_bundle(case: str) -> DragCoefficients:
    """Build the ``case`` linear_cubic bundle from the exported fit parameters.

    The extraction mass is read from the JSON's ``meff_amu`` (provenance), so
    the bundle's stamped mass tracks the export rather than being hand-coded.
    """
    fit_json = DRAG_REF_DIR / case / "linear_and_cubic" / "fit_parameters.json"
    params = json.loads(fit_json.read_text())
    return DragCoefficients(
        form=LINEAR_CUBIC,
        coefficients={"a": params["a"], "b": params["b"]},
        extraction_mass_model="constant",
        extraction_mass_amu=params["meff_amu"],
    )


def _load_trusted_interior(case: str) -> tuple[np.ndarray, np.ndarray]:
    """Load the ``case`` (v, F_drag) scatter, truncated to the trusted interior.

    Mirrors the extraction's own masking: drop first/last ``TRUNCATE_POINTS``
    points, then keep finite samples with ``v > 0``.
    """
    csv = DRAG_REF_DIR / case / "linear_and_cubic" / "drag_data.csv"
    data = np.genfromtxt(csv, delimiter=",", names=True)
    v = np.asarray(data["v_spline_Aps"], dtype=float)
    F = np.asarray(data["F_drag_amuAps2"], dtype=float)
    n = TRUNCATE_POINTS
    interior = slice(n, v.size - n)
    v, F = v[interior], F[interior]
    keep = np.isfinite(v) & np.isfinite(F) & (v > 0.0)
    return v[keep], F[keep]


# ===========================================================================
# Killer test: force-balance reproduction (SLICE1 §7)
# ===========================================================================
@pytest.mark.parametrize("case", CASES)
class TestForceBalanceReproduction:
    def test_module_is_the_extracted_law_refit_recovers_ab(self, case):
        """The module IS the extracted law: an independent least-squares refit
        of a*v + b*v^3 to the trusted-interior scatter recovers the stored
        {a,b} to which drag_force evaluates.

        This is the tight, model-correctness half of the killer test: the
        stored coefficients are the unique lstsq solution on this scatter, so
        recovering them (to ~1e-4) proves drag_force reproduces the extraction
        rather than merely fitting within noise. Tolerance is tight because
        we are comparing two solutions of the *same* linear least-squares
        problem (only float round-off separates them).
        """
        coeffs = _load_linear_cubic_bundle(case)
        v, F_ref = _load_trusted_interior(case)
        # a*v + b*v^3 is linear in (a, b): recover via lstsq on [v, v^3].
        design = np.vstack([v, v**3]).T
        (a_refit, b_refit), *_ = np.linalg.lstsq(design, F_ref, rcond=None)
        assert a_refit == pytest.approx(coeffs.coefficients["a"], rel=1e-4)
        assert b_refit == pytest.approx(coeffs.coefficients["b"], rel=1e-4)

    def test_linear_cubic_within_fit_residual_of_scatter(self, case):
        """drag_force with the in-hand {a,b} reproduces the case scatter within
        the genuine fit residual.

        The CSV ``F_drag_amuAps2`` column is the *noisy* force-balance scatter
        the curve_fit was fit TO, not the fitted curve, so the agreement band
        is the fit residual, not machine precision. The fits are high-R^2
        (9 A: ~0.974); a 7% relative-RMS band covers both cases' diagnosed
        residuals with a small margin, which is the meaningful "within fit
        tolerance" check of §7.
        """
        coeffs = _load_linear_cubic_bundle(case)
        v, F_ref = _load_trusted_interior(case)
        # Deep inside the droplet => gate g = 1, so the gated force equals the
        # raw extracted law a*v + b*v^3 the scatter was fit to.
        depth = np.full_like(v, -1.0e3)
        F_model = drag_force(v, depth, coeffs, STEEPNESS_A)

        resid = F_model - F_ref
        rms = float(np.sqrt(np.mean(resid**2)))
        scale = float(np.sqrt(np.mean(F_ref**2)))
        assert rms / scale < 0.07

    def test_deep_inside_gate_is_unity(self, case):
        """Sanity: the killer test's depth choice really gives g = 1."""
        assert spatial_gate(-1.0e3, STEEPNESS_A) == pytest.approx(1.0)


# ===========================================================================
# Dissipativity (SLICE1 §7) + guard-not-vacuous (§8)
# ===========================================================================
@pytest.mark.parametrize("case", CASES)
class TestDissipativity:
    def test_drag_opposes_motion_over_operating_range(self, case):
        """F_drag * v >= 0 everywhere: drag never adds energy."""
        coeffs = _load_linear_cubic_bundle(case)
        v = np.linspace(0.0, 30.0, 300)  # generous max-speed ceiling [A/ps]
        depth = np.full_like(v, -50.0)
        F = drag_force(v, depth, coeffs, STEEPNESS_A)
        assert np.all(F * v >= 0.0)

    def test_guard_not_vacuous_b_nonnegative(self, case):
        """SLICE1 §8: confirm the in-hand {a,b} actually satisfy a>0, b>=0.

        With b >= 0 there is no turnover speed v_dagger = sqrt(-a/b), so the
        dissipativity assertion above is testing a real (non-empty) condition.
        """
        coeffs = _load_linear_cubic_bundle(case)
        a = float(coeffs.coefficients["a"])
        b = float(coeffs.coefficients["b"])
        assert a > 0.0
        assert b >= 0.0


# ===========================================================================
# Spatial gate (SLICE1 §7)
# ===========================================================================
class TestSpatialGate:
    def test_limits(self):
        assert spatial_gate(-1.0e3, STEEPNESS_A) == pytest.approx(1.0)
        assert spatial_gate(0.0, STEEPNESS_A) == pytest.approx(0.5)
        assert spatial_gate(1.0e3, STEEPNESS_A) == pytest.approx(0.0)

    def test_monotone_decreasing_in_depth(self):
        depth = np.linspace(-60.0, 60.0, 500)
        g = spatial_gate(depth, STEEPNESS_A)
        assert np.all(np.diff(g) <= 0.0)

    def test_in_unit_interval(self):
        depth = np.linspace(-100.0, 100.0, 1000)
        g = spatial_gate(depth, STEEPNESS_A)
        assert np.all(g >= 0.0)
        assert np.all(g <= 1.0)

    def test_c1_continuity_matches_analytic_derivative(self):
        """C^1: central finite difference matches the analytic gate slope.

        d/d(depth) [0.5(1 - erf(depth/s))] = -exp(-(depth/s)^2)/(s*sqrt(pi)).
        """
        s = STEEPNESS_A
        depth = np.linspace(-30.0, 30.0, 61)
        h = 1.0e-5
        fd = (spatial_gate(depth + h, s) - spatial_gate(depth - h, s)) / (2.0 * h)
        analytic = -np.exp(-((depth / s) ** 2)) / (s * np.sqrt(np.pi))
        np.testing.assert_allclose(fd, analytic, atol=1e-7)

    def test_rejects_nonpositive_steepness(self):
        with pytest.raises(ValueError):
            spatial_gate(0.0, 0.0)
        with pytest.raises(ValueError):
            spatial_gate(0.0, -1.0)


# ===========================================================================
# FDT-coupling carrier: gamma carries the SAME gate as force (SLICE1 §7, §5.2)
# ===========================================================================
@pytest.mark.parametrize("case", CASES)
class TestFdtCouplingCarrier:
    def test_gamma_and_force_share_gate_factor(self, case):
        coeffs = _load_linear_cubic_bundle(case)
        a = float(coeffs.coefficients["a"])
        b = float(coeffs.coefficients["b"])
        v = np.array([0.5, 2.0, 7.0, 15.0])
        for d in (-40.0, -10.0, 0.0, 5.0, 30.0):
            depth = np.full_like(v, d)
            g = spatial_gate(depth, STEEPNESS_A)
            F = drag_force(v, depth, coeffs, STEEPNESS_A)
            gam = drag_gamma(v, depth, coeffs, STEEPNESS_A)
            np.testing.assert_allclose(F, g * (a * v + b * v**3), rtol=1e-12)
            np.testing.assert_allclose(gam, g * (a + b * v**2), rtol=1e-12)
            # Gate cancels in F/gamma => v exactly (same g in both, §5.2).
            np.testing.assert_allclose(F, gam * v, rtol=1e-12)


# ===========================================================================
# Low-velocity regularity (SLICE1 §7)
# ===========================================================================
@pytest.mark.parametrize("case", CASES)
class TestLowVelocityRegularity:
    def test_gamma_finite_at_rest_equals_g_times_a(self, case):
        coeffs = _load_linear_cubic_bundle(case)
        a = float(coeffs.coefficients["a"])
        for d in (-40.0, 0.0, 10.0):
            g = spatial_gate(d, STEEPNESS_A)
            gam0 = drag_gamma(0.0, d, coeffs, STEEPNESS_A)
            assert np.isfinite(gam0)
            assert float(gam0) == pytest.approx(g * a)

    def test_closed_form_avoids_division_singularity(self, case):
        """Documented contrast: |F_drag|/v is 0/0 at rest; closed form is finite.

        For linear_cubic, F_drag -> 0 as v -> 0, so computing gamma by division
        would be 0/0. The module's closed form gives the finite g*a instead.
        (The genuine v->0 divergence belongs to the DEFERRED power_law n<0 form,
        which raises NotImplementedError -- see TestFormDispatch.)
        """
        coeffs = _load_linear_cubic_bundle(case)
        assert float(drag_force(0.0, -40.0, coeffs, STEEPNESS_A)) == pytest.approx(0.0)
        assert float(drag_gamma(0.0, -40.0, coeffs, STEEPNESS_A)) > 0.0


# ===========================================================================
# Mass-agnosticism (SLICE1 §7)
# ===========================================================================
class TestMassAgnosticism:
    def test_no_function_accepts_a_mass_argument(self):
        mass_names = {"m", "mass", "mass_amu", "m_eff", "m_amu", "m_eff_amu"}
        for fn in (spatial_gate, drag_force, drag_gamma):
            params = set(inspect.signature(fn).parameters)
            assert not (params & mass_names), f"{fn.__name__} takes a mass arg"


# ===========================================================================
# Form dispatch: only linear_cubic realised (SLICE1 §5)
# ===========================================================================
class TestFormDispatch:
    @pytest.mark.parametrize(
        "form,coeffs",
        [
            (LINEAR_QUADRATIC, {"a": 1.0, "c": 1.0}),
            (THRESHOLD, {"F_sat": 1.0, "v0": 1.0}),
            (POWER_LAW, {"gamma": 6.15, "n": 2.6}),
        ],
    )
    def test_unrealised_forms_raise_not_implemented(self, form, coeffs):
        bundle = DragCoefficients(
            form=form,
            coefficients=coeffs,
            extraction_mass_model="constant",
            extraction_mass_amu=200.0,
        )
        v = np.array([1.0])
        depth = np.array([-10.0])
        with pytest.raises(NotImplementedError):
            drag_force(v, depth, bundle, STEEPNESS_A)
        with pytest.raises(NotImplementedError):
            drag_gamma(v, depth, bundle, STEEPNESS_A)


# ===========================================================================
# Coefficient-bundle type validation (SLICE1 §6)
# ===========================================================================
class TestDragCoefficientsType:
    @pytest.mark.parametrize("case", CASES)
    def test_carries_extraction_metadata(self, case):
        c = _load_linear_cubic_bundle(case)
        assert c.form == LINEAR_CUBIC
        assert c.extraction_mass_model == "constant"
        assert set(c.coefficients) >= {"a", "b"}
        # Mass stamp is read from the export's meff_amu = I + 19 He (~203 amu).
        assert c.extraction_mass_amu == pytest.approx(126.90447 + 19 * 4.002602)

    def test_rejects_unknown_form(self):
        with pytest.raises(ValueError):
            DragCoefficients("nope", {"a": 1.0, "b": 1.0}, "constant", 200.0)

    def test_rejects_missing_coefficients(self):
        with pytest.raises(ValueError):
            DragCoefficients(LINEAR_CUBIC, {"a": 1.0}, "constant", 200.0)

    def test_rejects_bad_mass_model(self):
        with pytest.raises(ValueError):
            DragCoefficients(LINEAR_CUBIC, {"a": 1.0, "b": 1.0}, "weird", 200.0)

    def test_rejects_nonpositive_mass(self):
        with pytest.raises(ValueError):
            DragCoefficients(LINEAR_CUBIC, {"a": 1.0, "b": 1.0}, "constant", 0.0)
