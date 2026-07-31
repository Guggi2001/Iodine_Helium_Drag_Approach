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
    CAPPED_CUBIC,
    CAPPED_LINEAR_QUADRATIC,
    LINEAR_CUBIC,
    LINEAR_QUADRATIC,
    POWER_LAW,
    PURE_LINEAR,
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
        (The genuine v->0 divergence belongs to a hypothetical power_law n<1
        law, which the §3.3 config guard refuses -- the realised power_law is
        bounded to n >= 1; see TestFormPhaseFamilies.)
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
# Form dispatch: threshold reserved; lq/pl realised (SLICE1 §5, METHOD_B §10)
# ===========================================================================
def _bundle(form: str, coeffs: dict) -> DragCoefficients:
    return DragCoefficients(
        form=form,
        coefficients=coeffs,
        extraction_mass_model="constant",
        extraction_mass_amu=200.0,
    )


class TestFormDispatch:
    def test_threshold_raises_not_implemented(self):
        bundle = _bundle(THRESHOLD, {"F_sat": 1.0, "v0": 1.0})
        v = np.array([1.0])
        depth = np.array([-10.0])
        with pytest.raises(NotImplementedError):
            drag_force(v, depth, bundle, STEEPNESS_A)
        with pytest.raises(NotImplementedError):
            drag_gamma(v, depth, bundle, STEEPNESS_A)

    @pytest.mark.parametrize(
        "form,coeffs",
        [
            (LINEAR_QUADRATIC, {"a": 1.0, "c": 1.0}),
            (POWER_LAW, {"C": 6.15, "n": 2.6}),
        ],
    )
    def test_form_phase_families_are_realised(self, form, coeffs):
        bundle = _bundle(form, coeffs)
        v = np.array([1.0, 3.0])
        depth = np.array([-10.0, -10.0])
        assert np.all(np.isfinite(drag_force(v, depth, bundle, STEEPNESS_A)))
        assert np.all(np.isfinite(drag_gamma(v, depth, bundle, STEEPNESS_A)))


# ===========================================================================
# METHOD_B §10 form-phase families: closed forms, v=0 regularity, nesting
# ===========================================================================
class TestFormPhaseFamilies:
    """linear_quadratic / power_law governing equations (METHOD_B §10.3)."""

    V = np.linspace(0.2, 6.0, 25)         # away from v=0 (division is legal)
    DEPTHS = np.array([-40.0, -5.0, 0.0, 10.0])

    @pytest.mark.parametrize(
        "form,coeffs",
        [
            (LINEAR_QUADRATIC, {"a": 4.0, "c": 11.0}),
            (LINEAR_QUADRATIC, {"a": 0.0, "c": 11.0}),   # pure-quadratic
            (POWER_LAW, {"C": 10.36, "n": 2.056}),
            (POWER_LAW, {"C": 3.0, "n": 1.0}),           # n = 1 boundary
        ],
    )
    def test_gamma_closed_form_equals_force_over_v(self, form, coeffs):
        # gamma(v) = |F_drag(v)| / v away from rest -- the defining identity
        # of the force-coefficient convention; analytical, tight tolerance.
        bundle = _bundle(form, coeffs)
        for d in self.DEPTHS:
            F = drag_force(self.V, d, bundle, STEEPNESS_A)
            gam = drag_gamma(self.V, d, bundle, STEEPNESS_A)
            np.testing.assert_allclose(gam, F / self.V, rtol=1e-12)

    def test_linear_quadratic_explicit_values(self):
        # F = g*(a*v + c*v^2), gamma = g*(a + c*v) at unit gate (deep inside).
        bundle = _bundle(LINEAR_QUADRATIC, {"a": 2.0, "c": 3.0})
        d = -400.0  # gate -> 1 to round-off
        v = 2.0
        assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            2.0 * v + 3.0 * v**2, rel=1e-12
        )
        assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            2.0 + 3.0 * v, rel=1e-12
        )

    def test_power_law_explicit_values(self):
        # F = g*C*v^n, gamma = g*C*v^(n-1) at unit gate (deep inside).
        bundle = _bundle(POWER_LAW, {"C": 5.0, "n": 2.5})
        d = -400.0
        v = 3.0
        assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            5.0 * v**2.5, rel=1e-12
        )
        assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            5.0 * v**1.5, rel=1e-12
        )

    def test_linear_quadratic_gamma_finite_at_rest(self):
        # gamma -> g*a at v = 0 (closed form; division would be 0/0).
        bundle = _bundle(LINEAR_QUADRATIC, {"a": 4.0, "c": 11.0})
        for d in (-40.0, 0.0, 10.0):
            g = spatial_gate(d, STEEPNESS_A)
            gam0 = drag_gamma(0.0, d, bundle, STEEPNESS_A)
            assert np.isfinite(gam0)
            assert float(gam0) == pytest.approx(g * 4.0)

    def test_power_law_gamma_at_rest_n_above_one_is_zero(self):
        bundle = _bundle(POWER_LAW, {"C": 10.0, "n": 2.5})
        gam0 = drag_gamma(0.0, -40.0, bundle, STEEPNESS_A)
        assert float(gam0) == 0.0

    def test_power_law_gamma_at_rest_n_equal_one_is_g_times_C(self):
        # numpy 0.0**0.0 == 1.0 realizes the n = 1 limit exactly.
        bundle = _bundle(POWER_LAW, {"C": 10.0, "n": 1.0})
        for d in (-40.0, 0.0):
            g = spatial_gate(d, STEEPNESS_A)
            gam0 = drag_gamma(0.0, d, bundle, STEEPNESS_A)
            assert float(gam0) == pytest.approx(g * 10.0)

    @pytest.mark.parametrize(
        "form,coeffs",
        [
            (LINEAR_QUADRATIC, {"a": 4.0, "c": 11.0}),
            (POWER_LAW, {"C": 10.36, "n": 2.056}),
        ],
    )
    def test_dissipative_and_gate_shared(self, form, coeffs):
        # Drag opposes motion (F >= 0 magnitude convention) and gamma carries
        # the SAME gate factor as the force (hard FDT coupling carrier, §5.2).
        bundle = _bundle(form, coeffs)
        for d in self.DEPTHS:
            F = drag_force(self.V, d, bundle, STEEPNESS_A)
            gam = drag_gamma(self.V, d, bundle, STEEPNESS_A)
            assert np.all(F >= 0.0)
            assert np.all(gam >= 0.0)
            np.testing.assert_allclose(F, gam * self.V, rtol=1e-12)


class TestNestingIdentities:
    """§10.3 cross-check obligations: power_law nests both incumbents EXACTLY.

    power_law(n=2, C=c) == linear_quadratic(a=0, c)   (pure-quadratic)
    power_law(n=3, C=b) == linear_cubic(a=0, b)       (pure-cubic)

    Analytical identities -- tight tolerance (rtol 1e-12; the only float
    difference is x*x*x vs x**3.0 evaluation order).
    """

    V = np.linspace(0.0, 6.0, 31)          # includes v = 0
    DEPTHS = np.array([-40.0, -5.0, 0.0, 10.0])

    def test_power_law_n2_equals_pure_quadratic(self):
        c = 11.016050300970692              # the locked c0 scale (§10.4.1)
        pl = _bundle(POWER_LAW, {"C": c, "n": 2.0})
        lq = _bundle(LINEAR_QUADRATIC, {"a": 0.0, "c": c})
        for d in self.DEPTHS:
            np.testing.assert_allclose(
                drag_force(self.V, d, pl, STEEPNESS_A),
                drag_force(self.V, d, lq, STEEPNESS_A),
                rtol=1e-12, atol=0.0,
            )
            np.testing.assert_allclose(
                drag_gamma(self.V, d, pl, STEEPNESS_A),
                drag_gamma(self.V, d, lq, STEEPNESS_A),
                rtol=1e-12, atol=0.0,
            )

    def test_power_law_n3_equals_pure_cubic(self):
        b = 2.5153654   # the shared_pure_cubic production scale
        pl = _bundle(POWER_LAW, {"C": b, "n": 3.0})
        lc = _bundle(LINEAR_CUBIC, {"a": 0.0, "b": b})
        for d in self.DEPTHS:
            np.testing.assert_allclose(
                drag_force(self.V, d, pl, STEEPNESS_A),
                drag_force(self.V, d, lc, STEEPNESS_A),
                rtol=1e-12, atol=0.0,
            )
            np.testing.assert_allclose(
                drag_gamma(self.V, d, pl, STEEPNESS_A),
                drag_gamma(self.V, d, lc, STEEPNESS_A),
                rtol=1e-12, atol=0.0,
            )


# ===========================================================================
# capped_cubic (Tier-2 Addendum I §I.10 Slice T1): locked in-band pure cubic
# + high-v tail gamma = b*v_c^2*(v/v_c)^p_tail, p_tail in {0, -1}
# ===========================================================================
class TestCappedCubic:
    """§I.3/§I.10 tail family: gamma = g*b*v^2 (v <= v_c), g*b*v_c^2*(v/v_c)^p
    (v > v_c); F = gamma*v. Continuous at v_c; exactly the locked pure cubic
    in-band. Analytical pins are tight (rtol 1e-12); byte-identity checks are
    exact (``==``, no tolerance) because the in-band branch must be the SAME
    arithmetic as ``linear_cubic(a=0)`` -- the Tier-0 lock rides on it.
    """

    B = 2.0          # amu*ps/A^2 (arbitrary positive test value)
    V_C = 3.0        # A/ps
    DEPTHS = np.array([-400.0, -5.0, 0.0, 10.0])
    # speeds straddling the cap, including rest and the cap itself
    V = np.array([0.0, 0.5, 1.5, 3.0, 3.5, 4.5, 6.0])

    def _capped(self, p_tail, v_c=V_C, b=B):
        return _bundle(
            CAPPED_CUBIC, {"b": b, "v_c": v_c, "p_tail": p_tail}
        )

    # --- closed-form pins per branch (§I.10 oracle 1) ---
    def test_in_band_explicit_values(self):
        # v < v_c: F = g*b*v^3, gamma = g*b*v^2 at unit gate (deep inside).
        bundle = self._capped(p_tail=0.0)
        d, v = -400.0, 2.0
        assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            self.B * v**3, rel=1e-12
        )
        assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            self.B * v**2, rel=1e-12
        )

    def test_tail_p0_explicit_values(self):
        # p_tail = 0 (Stokes-like tail): gamma = g*b*v_c^2 constant above the
        # cap; F = g*b*v_c^2*v grows linearly.
        bundle = self._capped(p_tail=0.0)
        d, v = -400.0, 6.0
        assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            self.B * self.V_C**2, rel=1e-12
        )
        assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            self.B * self.V_C**2 * v, rel=1e-12
        )

    def test_tail_pm1_explicit_values(self):
        # p_tail = -1 (saturated/plastic tail): F = g*b*v_c^3 CONSTANT above
        # the cap; gamma = g*b*v_c^3/v falls off as 1/v.
        bundle = self._capped(p_tail=-1.0)
        d = -400.0
        F_sat = self.B * self.V_C**3
        for v in (4.0, 6.0):
            assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
                F_sat, rel=1e-12
            )
            assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
                F_sat / v, rel=1e-12
            )

    def test_tail_pm2_explicit_values(self):
        # p_tail = -2 (§3.5h softening tail): F = g*b*v_c^4/v falls off as
        # 1/v above the cap — an explicit exponent-law check the generic
        # identities cannot catch (an off-by-one in p_tail would survive
        # F = gamma*v and continuity-at-the-cap).
        bundle = self._capped(p_tail=-2.0)
        d = -400.0
        for v in (4.0, 6.0):
            expected_F = self.B * self.V_C**4 / v
            assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
                expected_F, rel=1e-12
            )
            assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
                expected_F / v, rel=1e-12
            )

    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    def test_continuous_at_the_cap(self, p_tail):
        # Both branches evaluate to g*b*v_c^2 at v = v_c; approach from both
        # sides agrees to first order (eps*rel band, analytical continuity).
        bundle = self._capped(p_tail=p_tail)
        d = -5.0
        g = float(spatial_gate(d, STEEPNESS_A))
        expected = g * self.B * self.V_C**2
        eps = 1e-9
        for v in (self.V_C - eps, self.V_C, self.V_C + eps):
            assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
                expected, rel=1e-6
            )

    # --- convention identities (mirror TestFormPhaseFamilies) ---
    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    def test_gamma_closed_form_equals_force_over_v(self, p_tail):
        bundle = self._capped(p_tail=p_tail)
        v = self.V[self.V > 0.0]  # away from rest (division is legal)
        for d in self.DEPTHS:
            F = drag_force(v, d, bundle, STEEPNESS_A)
            gam = drag_gamma(v, d, bundle, STEEPNESS_A)
            np.testing.assert_allclose(gam, F / v, rtol=1e-12)

    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    def test_dissipative_and_gate_shared(self, p_tail):
        # F >= 0 (magnitude convention) and gamma carries the SAME gate factor
        # as the force (hard FDT coupling carrier, §5.2).
        bundle = self._capped(p_tail=p_tail)
        for d in self.DEPTHS:
            F = drag_force(self.V, d, bundle, STEEPNESS_A)
            gam = drag_gamma(self.V, d, bundle, STEEPNESS_A)
            assert np.all(F >= 0.0)
            assert np.all(gam >= 0.0)
            np.testing.assert_allclose(F, gam * self.V, rtol=1e-12)

    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    def test_gamma_at_rest_is_zero(self, p_tail):
        # v = 0 sits on the in-band pure-cubic branch: gamma -> 0 regularly
        # (the p_tail = -1 tail formula is never evaluated at rest).
        bundle = self._capped(p_tail=p_tail)
        gam0 = drag_gamma(0.0, -40.0, bundle, STEEPNESS_A)
        assert float(gam0) == 0.0

    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    def test_no_floating_point_warnings_across_the_cap(self, p_tail):
        # Mixed-speed arrays (rest + in-band + tail) and v_c = inf must never
        # touch a singular intermediate (0**-1, inf/inf): errstate-raise makes
        # any such leak a hard failure.
        for v_c in (self.V_C, np.inf):
            bundle = self._capped(p_tail=p_tail, v_c=v_c)
            with np.errstate(all="raise"):
                F = drag_force(self.V, -5.0, bundle, STEEPNESS_A)
                gam = drag_gamma(self.V, -5.0, bundle, STEEPNESS_A)
            assert np.all(np.isfinite(F))
            assert np.all(np.isfinite(gam))

    # --- §I.10 oracle 2: v_c = inf (or >= v_max) byte-identity with the
    # --- locked pure cubic. Exact ``==``: same arithmetic, not "close".
    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    @pytest.mark.parametrize("v_c", [np.inf, 6.0])  # 6.0 == max(V) >= v_max
    def test_vc_at_or_above_vmax_is_byte_identical_to_pure_cubic(
        self, p_tail, v_c
    ):
        b = 2.5153509  # the locked shared_pure_cubic scale (§I.3)
        capped = self._capped(p_tail=p_tail, v_c=v_c, b=b)
        pure = _bundle(LINEAR_CUBIC, {"a": 0.0, "b": b})
        for d in self.DEPTHS:
            F_capped = drag_force(self.V, d, capped, STEEPNESS_A)
            F_pure = drag_force(self.V, d, pure, STEEPNESS_A)
            assert np.array_equal(F_capped, F_pure)
            gam_capped = drag_gamma(self.V, d, capped, STEEPNESS_A)
            gam_pure = drag_gamma(self.V, d, pure, STEEPNESS_A)
            assert np.array_equal(gam_capped, gam_pure)

    def test_in_band_speeds_byte_identical_below_finite_cap(self):
        # Even with a finite in-range cap, every v <= v_c is the SAME
        # arithmetic as the locked pure cubic (the Tier-0 in-band lock).
        b = 2.5153509
        capped = self._capped(p_tail=-1.0, v_c=self.V_C, b=b)
        pure = _bundle(LINEAR_CUBIC, {"a": 0.0, "b": b})
        v = self.V[self.V <= self.V_C]
        for d in self.DEPTHS:
            assert np.array_equal(
                drag_force(v, d, capped, STEEPNESS_A),
                drag_force(v, d, pure, STEEPNESS_A),
            )
            assert np.array_equal(
                drag_gamma(v, d, capped, STEEPNESS_A),
                drag_gamma(v, d, pure, STEEPNESS_A),
            )

    def test_missing_coefficients_rejected(self):
        with pytest.raises(ValueError, match="missing"):
            _bundle(CAPPED_CUBIC, {"b": 2.0})  # no v_c, no p_tail


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

    def test_extraction_method_defaults_to_force_balance(self):
        # Legacy construction signature unchanged -> Method-A provenance with
        # no jointly-validated binding (None, handled by the §6.5.1 guard).
        c = DragCoefficients(LINEAR_CUBIC, {"a": 1.0, "b": 1.0}, "constant", 200.0)
        assert c.extraction_method == "force_balance"
        assert c.effective_binding_energy_I_ion_eV is None

    def test_trajectory_matching_carries_binding(self):
        c = DragCoefficients(
            LINEAR_CUBIC,
            {"a": 1.0, "b": 1.0},
            "constant",
            200.0,
            extraction_method="trajectory_matching",
            effective_binding_energy_I_ion_eV=0.2,
        )
        assert c.extraction_method == "trajectory_matching"
        assert c.effective_binding_energy_I_ion_eV == 0.2

    def test_rejects_unknown_extraction_method(self):
        with pytest.raises(ValueError, match="extraction_method"):
            DragCoefficients(
                LINEAR_CUBIC,
                {"a": 1.0, "b": 1.0},
                "constant",
                200.0,
                extraction_method="hand_tuned",
            )

    def test_rejects_nonpositive_binding(self):
        with pytest.raises(ValueError, match="effective_binding_energy"):
            DragCoefficients(
                LINEAR_CUBIC,
                {"a": 1.0, "b": 1.0},
                "constant",
                200.0,
                effective_binding_energy_I_ion_eV=0.0,
            )


class TestSpatialGateSingleSource:
    """Slice-rho single-source regression: ``spatial_gate`` routes through the
    shared ``_gates._erf_complement`` helper (no re-inlined erf formula), so the
    drag gate and the Tier-2 ``rho_He/rho_bulk`` density gate stay one formula in
    one place. Locks the Tier-0 no-behaviour-change refactor from the drag side.
    """

    def test_spatial_gate_equals_erf_complement_on_grid(self):
        from i2_helium_md.physics._gates import _erf_complement

        depth = np.linspace(-60.0, 60.0, 241)
        np.testing.assert_array_equal(
            np.asarray(spatial_gate(depth, STEEPNESS_A)),
            np.asarray(_erf_complement(depth, STEEPNESS_A)),
        )

    def test_spatial_gate_still_fails_loud_on_nonpositive_steepness(self):
        # The guard moved into _erf_complement; spatial_gate must still raise.
        with pytest.raises(ValueError):
            spatial_gate(0.0, 0.0)
        with pytest.raises(ValueError):
            spatial_gate(0.0, -1.0)


class TestCappedLinearQuadratic:
    """Atlas §6.6 counterfactual form: gamma = g*(a + c*v) (v <= v_c),
    g*(a + c*v_c)*(v/v_c)^p (v > v_c); F = gamma*v. Mirrors the
    ``capped_cubic`` tail conventions exactly (shared tail scaffold), and the
    in-band branch must be the SAME arithmetic as ``linear_quadratic`` --
    the byte-identity the twin/MD comparison rides on. Analytical pins tight
    (rtol 1e-12); byte-identity exact (``==``).
    """

    A = 9.805022771384936e-05  # amu/ps (the shared lq artifact value)
    C = 12.792201727123915     # amu/A  (the shared lq artifact value)
    V_C = 9.0                  # A/ps   (spot-check case-A cap)
    DEPTHS = np.array([-400.0, -5.0, 0.0, 10.0])
    V = np.array([0.0, 0.5, 2.0, 9.0, 9.5, 10.5, 12.0])

    def _capped(self, p_tail, v_c=V_C, a=A, c=C):
        return _bundle(
            CAPPED_LINEAR_QUADRATIC,
            {"a": a, "c": c, "v_c": v_c, "p_tail": p_tail},
        )

    # --- closed-form pins per branch ---
    def test_in_band_explicit_values(self):
        bundle = self._capped(p_tail=-1.0)
        d, v = -400.0, 2.0
        assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            self.A * v + self.C * v**2, rel=1e-12
        )
        assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            self.A + self.C * v, rel=1e-12
        )

    def test_tail_pm1_constant_force(self):
        # p_tail = -1: F = g*(a + c*v_c)*v_c CONSTANT above the cap.
        bundle = self._capped(p_tail=-1.0)
        d = -400.0
        F_sat = (self.A + self.C * self.V_C) * self.V_C
        for v in (10.5, 12.0):
            assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
                F_sat, rel=1e-12
            )
            assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
                F_sat / v, rel=1e-12
            )

    def test_in_band_byte_identity_with_linear_quadratic(self):
        # v <= v_c must be the SAME arithmetic as the delivered lq branch
        # (exact ==, no tolerance): the Tier-0 lq artifact semantics carry.
        capped = self._capped(p_tail=-1.0)
        lq = _bundle(LINEAR_QUADRATIC, {"a": self.A, "c": self.C})
        v_in = self.V[self.V <= self.V_C]
        for d in self.DEPTHS:
            np.testing.assert_array_equal(
                np.asarray(drag_force(v_in, d, capped, STEEPNESS_A)),
                np.asarray(drag_force(v_in, d, lq, STEEPNESS_A)),
            )
            np.testing.assert_array_equal(
                np.asarray(drag_gamma(v_in, d, capped, STEEPNESS_A)),
                np.asarray(drag_gamma(v_in, d, lq, STEEPNESS_A)),
            )

    def test_v_c_inf_byte_identity_with_linear_quadratic(self):
        # v_c = inf never enters the tail branch: full byte-identity.
        capped = self._capped(p_tail=-1.0, v_c=np.inf)
        lq = _bundle(LINEAR_QUADRATIC, {"a": self.A, "c": self.C})
        for d in self.DEPTHS:
            np.testing.assert_array_equal(
                np.asarray(drag_force(self.V, d, capped, STEEPNESS_A)),
                np.asarray(drag_force(self.V, d, lq, STEEPNESS_A)),
            )

    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    def test_continuous_at_the_cap(self, p_tail):
        bundle = self._capped(p_tail=p_tail)
        d = -5.0
        g = float(spatial_gate(d, STEEPNESS_A))
        expected = g * (self.A + self.C * self.V_C)
        eps = 1e-9
        for v in (self.V_C - eps, self.V_C, self.V_C + eps):
            assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
                expected, rel=1e-6
            )

    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    def test_dissipative_gate_shared_and_gamma_is_force_over_v(self, p_tail):
        bundle = self._capped(p_tail=p_tail)
        for d in self.DEPTHS:
            F = drag_force(self.V, d, bundle, STEEPNESS_A)
            gam = drag_gamma(self.V, d, bundle, STEEPNESS_A)
            assert np.all(F >= 0.0)
            assert np.all(gam >= 0.0)
            np.testing.assert_allclose(F, gam * self.V, rtol=1e-12)

    @pytest.mark.parametrize("p_tail", [0.0, -1.0, -2.0, -3.0])
    def test_no_floating_point_warnings_across_the_cap(self, p_tail):
        # Mixed arrays (rest + in-band + tail) and v_c = inf: no singular
        # intermediate may be touched (errstate-raise makes leaks fatal).
        for v_c in (self.V_C, np.inf):
            bundle = self._capped(p_tail=p_tail, v_c=v_c)
            with np.errstate(all="raise"):
                drag_force(self.V, -5.0, bundle, STEEPNESS_A)
                drag_gamma(self.V, -5.0, bundle, STEEPNESS_A)

    def test_gamma_at_rest_is_regular_a(self):
        # v = 0 sits in-band: gamma -> g*a (the lq rest limit), never the tail.
        bundle = self._capped(p_tail=-1.0)
        gam0 = float(drag_gamma(0.0, -400.0, bundle, STEEPNESS_A))
        assert gam0 == pytest.approx(self.A, rel=1e-12)


class TestPureLinear:
    """Free-form linear ring form (TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN §1):
    gamma = g*a constant in v, F = g*a*v. The family is exactly the a-corner
    of ``linear_cubic`` (b = 0) -- the identity is locked EXACTLY (``==``),
    because the MD ring's twin comparison rides on the form being nothing
    but the constant-gamma law. Analytical pins tight (rtol 1e-12).
    """

    A = 27.5                    # amu/ps (a ring sub-plateau value)
    DEPTHS = np.array([-400.0, -40.0, -5.0, 0.0, 10.0])
    V = np.linspace(0.0, 14.0, 29)   # includes v = 0 and fast-class speeds

    def _lin(self, a=A):
        return _bundle(PURE_LINEAR, {"a": a})

    def test_explicit_values_at_unit_gate(self):
        bundle = self._lin()
        d, v = -400.0, 3.0
        assert float(drag_force(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            self.A * v, rel=1e-12
        )
        assert float(drag_gamma(v, d, bundle, STEEPNESS_A)) == pytest.approx(
            self.A, rel=1e-12
        )

    def test_gamma_constant_in_v_and_regular_at_rest(self):
        # The defining property of the family: one number sets the force at
        # every velocity. gamma(v) is exactly constant across the whole grid
        # (including v = 0 -- no 0/0; closed form).
        bundle = self._lin()
        for d in self.DEPTHS:
            g = spatial_gate(d, STEEPNESS_A)
            gam = np.asarray(drag_gamma(self.V, d, bundle, STEEPNESS_A))
            np.testing.assert_array_equal(gam, np.full_like(self.V, g * self.A))

    def test_exact_identity_with_linear_cubic_b0(self):
        # pure_linear(a) == linear_cubic(a, b = 0) bitwise, force AND gamma:
        # the a-corner nesting identity (the §10.3 obligation pattern).
        lin = self._lin()
        lc = _bundle(LINEAR_CUBIC, {"a": self.A, "b": 0.0})
        for d in self.DEPTHS:
            np.testing.assert_array_equal(
                np.asarray(drag_force(self.V, d, lin, STEEPNESS_A)),
                np.asarray(drag_force(self.V, d, lc, STEEPNESS_A)),
            )
            np.testing.assert_array_equal(
                np.asarray(drag_gamma(self.V, d, lin, STEEPNESS_A)),
                np.asarray(drag_gamma(self.V, d, lc, STEEPNESS_A)),
            )

    def test_dissipative_and_gate_shared(self):
        # F >= 0 magnitude convention; gamma carries the SAME gate as the
        # force (FDT coupling carrier, §5.2): F == gamma * v on the grid.
        bundle = self._lin()
        for d in self.DEPTHS:
            F = np.asarray(drag_force(self.V, d, bundle, STEEPNESS_A))
            gam = np.asarray(drag_gamma(self.V, d, bundle, STEEPNESS_A))
            assert np.all(F >= 0.0)
            assert np.all(gam >= 0.0)
            np.testing.assert_allclose(F, gam * self.V, rtol=1e-12)

    def test_requires_coefficient_a(self):
        with pytest.raises(ValueError, match="missing"):
            _bundle(PURE_LINEAR, {})

    def test_free_form_extraction_method_accepted(self):
        # The ring bundles carry the honest provenance vocabulary.
        bundle = DragCoefficients(
            form=PURE_LINEAR,
            coefficients={"a": self.A},
            extraction_mass_model="constant",
            extraction_mass_amu=202.953908,
            extraction_method="free_form",
            effective_binding_energy_I_ion_eV=None,
        )
        assert bundle.extraction_method == "free_form"

    def test_free_form_forbids_binding_stamp(self):
        # A free-form parameter was never jointly calibrated with any
        # binding: a stamp would claim a nonexistent §6.5.1 validation.
        with pytest.raises(ValueError, match="free_form"):
            DragCoefficients(
                form=PURE_LINEAR,
                coefficients={"a": self.A},
                extraction_mass_model="constant",
                extraction_mass_amu=202.953908,
                extraction_method="free_form",
                effective_binding_energy_I_ion_eV=0.1168,
            )
