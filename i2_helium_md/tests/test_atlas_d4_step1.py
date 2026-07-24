"""Unit coverage for scripts/extraction/atlas_d4_step1_report.py.

Covers only the *new* arithmetic that script introduces: the two parked-form
gamma(v) formulas (Pade saturating cubic, subtractive gated cubic -- atlas
plan section 6.1 entry-gate arithmetic) and the gate = 1 stand-in used to
evaluate realised forms through :func:`i2_helium_md.physics.drag.drag_gamma`.
The form-table section is pure artifact echo (no arithmetic) and the fits
themselves are Tier-0 artifacts -- neither is re-tested here. No propagation,
no figures.

Tolerances: the parked-form checks are closed-form identities (exact or
float-tight 1e-12); the Pade high-v asymptote uses the analytic remainder
(v/v_s)^3 / (1 + (v/v_s)^3) at the probe point, not a loose band.
"""

from __future__ import annotations

import numpy as np
import pytest

from scripts.extraction.atlas_d4_step1_report import (
    build_gamma_rows,
    gamma_pade,
    gamma_realised,
    gamma_subtractive,
)
from i2_helium_md.physics.drag import CAPPED_CUBIC, LINEAR_CUBIC


B = 2.5154  # locked Tier-0 shared b [amu*ps/A^2] (value irrelevant to shape)


class TestGammaPade:
    def test_reduces_to_cubic_at_low_v(self):
        # (v/v_s)^3 = 1e-6 at v = 0.01*v_s -> ratio to b*v^2 within 1e-5
        v, v_s = 0.0725, 7.25
        assert gamma_pade(B, v_s, v) == pytest.approx(B * v**2, rel=1e-5)

    def test_half_cubic_exactly_at_v_s(self):
        v_s = 7.25
        assert gamma_pade(B, v_s, v_s) == pytest.approx(
            B * v_s**2 / 2.0, rel=1e-12
        )

    def test_force_asymptote_is_constant_b_vs_cubed(self):
        # F = gamma*v -> b*v_s^3 * x/(1+x), x = (v/v_s)^3; at v = 20*v_s the
        # analytic remainder is 8000/8001, so compare against that exactly.
        v_s = 7.25
        v = 20.0 * v_s
        force = float(gamma_pade(B, v_s, v)) * v
        assert force == pytest.approx(
            B * v_s**3 * 8000.0 / 8001.0, rel=1e-12
        )


class TestGammaSubtractive:
    def test_zero_at_and_below_v_f(self):
        v_f = 1.5
        assert gamma_subtractive(B, v_f, v_f) == 0.0
        assert gamma_subtractive(B, v_f, 0.5 * v_f) == 0.0

    def test_twin_window_suppression_at_v2(self):
        # the plan's quoted x0.44 at v = 2 with v_f = 1.5 is exactly
        # 1 - v_f^2/v^2 = 1 - 2.25/4 = 0.4375
        v, v_f = 2.0, 1.5
        ratio = float(gamma_subtractive(B, v_f, v)) / (B * v**2)
        assert ratio == pytest.approx(0.4375, rel=1e-12)

    def test_approaches_cubic_at_high_v(self):
        v, v_f = 50.0, 2.0
        ratio = float(gamma_subtractive(B, v_f, v)) / (B * v**2)
        assert ratio == pytest.approx(1.0 - v_f**2 / v**2, rel=1e-12)


class TestGateOneStandIn:
    def test_linear_cubic_at_deep_depth_is_ungated(self):
        # depth = -50 A with steepness 1 A puts erf at its exact 1.0 plateau,
        # so gamma must equal the closed form with no gate attenuation.
        v = np.array([2.0, 4.95, 10.5])
        got = gamma_realised(LINEAR_CUBIC, {"a": 0.0, "b": B}, v)
        np.testing.assert_allclose(got, B * v**2, rtol=1e-12)

    def test_capped_cubic_tail_constant_force(self):
        # p_tail = -1 above v_c = 7.25 -> gamma*v == b*v_c^3 identically
        v_c = 7.25
        v = np.array([8.0, 10.5])
        got = gamma_realised(
            CAPPED_CUBIC, {"b": B, "v_c": v_c, "p_tail": -1.0}, v
        )
        np.testing.assert_allclose(got * v, B * v_c**3, rtol=1e-12)


def test_build_gamma_rows_reads_committed_artifacts():
    """The row builder loads the committed bundles and yields finite gammas.

    Also pins the artifact-integrity anchor the whole report rides on: the
    shared pure-cubic b it reads must round to the locked 2.5154.
    """
    rows = build_gamma_rows()
    labels = [label for label, _ in rows]
    assert labels[0].startswith("capped_cubic")
    assert any("PARKED Pade" in x for x in labels)
    assert sum("PARKED subtractive" in x for x in labels) == 2
    for label, g in rows:
        assert np.all(np.isfinite(g)), label
        assert np.all(g >= 0.0), label
    # locked-b anchor: production row at v = 2 (in-band) is b*4 with gate 1
    b_read = float(rows[0][1][0]) / 4.0
    assert round(b_read, 4) == 2.5154
