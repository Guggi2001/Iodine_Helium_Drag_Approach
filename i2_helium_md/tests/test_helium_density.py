"""Tests for the helium density gate (Tier-2 Phase-B Slice rho).

Slice rho is the ``rho_He(depth)/rho_bulk in [0, 1]`` field that gates pickup
(``physics/helium_density.py``), plus the single-source erf-complement helper it
shares with ``drag.spatial_gate`` (``physics/_gates.py``).

Oracle (TIER2_PHASE_B_IMPLEMENTATION_PLAN.md §2.2 (rho) / §3.1 / §4):
* ratio = 0.5 at depth 0; -> 1 deep inside (depth << 0); -> 0 outside (depth >> 0);
* monotone non-increasing in depth;
* **parity / single-source regression:** ``rho_he_ratio`` equals ``drag.spatial_gate``
  to machine precision (both route through ``_erf_complement``);
* ``_erf_complement`` (and through it both callers) fail-loud on ``steepness <= 0``;
* the tabulated arm round-trips a hand-built ``[0, 1]`` profile and clamps the tails.

Pure / config-agnostic: depth + steepness only (no ``n``, no Langmuir cap, no ``gamma``,
no mass -- those live in Slice P / the drag law).
"""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics._gates import _erf_complement
from i2_helium_md.physics.drag import spatial_gate
from i2_helium_md.physics.helium_density import (
    TabulatedDensityProfile,
    rho_he_ratio,
    tabulated_density_profile,
)

# The sourced surface width (confining-potential steepness, 14.2 A; CALIBRATION
# row 5). rho reuses the *same* number the drag gate uses.
STEEPNESS_A = 14.2


class TestErfComplementHelper:
    """The extracted single-source ``0.5*(1 - erf(depth/steepness))``."""

    def test_half_at_surface(self):
        assert _erf_complement(0.0, STEEPNESS_A) == pytest.approx(0.5)

    def test_saturates_exactly_at_tails(self):
        # erf saturates to +-1 in float at large |depth| -- exact, no epsilon
        # (the Slice-K/L tail-saturation precedent).
        assert _erf_complement(-1.0e3, STEEPNESS_A) == 1.0
        assert _erf_complement(1.0e3, STEEPNESS_A) == 0.0

    def test_fail_loud_on_nonpositive_steepness(self):
        with pytest.raises(ValueError):
            _erf_complement(0.0, 0.0)
        with pytest.raises(ValueError):
            _erf_complement(0.0, -1.0)


class TestRhoHeRatioOracle:
    """The rho_He/rho_bulk gate values."""

    def test_half_at_surface(self):
        assert rho_he_ratio(0.0, steepness=STEEPNESS_A) == pytest.approx(0.5)

    def test_one_deep_inside_zero_outside(self):
        assert rho_he_ratio(-1.0e3, steepness=STEEPNESS_A) == 1.0
        assert rho_he_ratio(1.0e3, steepness=STEEPNESS_A) == 0.0

    def test_in_unit_interval(self):
        depth = np.linspace(-60.0, 60.0, 121)
        r = rho_he_ratio(depth, steepness=STEEPNESS_A)
        assert np.all(r >= 0.0) and np.all(r <= 1.0)

    def test_monotone_non_increasing_in_depth(self):
        depth = np.linspace(-60.0, 60.0, 121)
        r = rho_he_ratio(depth, steepness=STEEPNESS_A)
        assert np.all(np.diff(r) <= 0.0)

    def test_fail_loud_on_nonpositive_steepness(self):
        with pytest.raises(ValueError):
            rho_he_ratio(0.0, steepness=0.0)
        with pytest.raises(ValueError):
            rho_he_ratio(0.0, steepness=-14.2)


class TestVectorization:
    """Scalar-in -> float, array-in -> ndarray (shell_schedule idiom)."""

    def test_scalar_returns_float(self):
        out = rho_he_ratio(3.0, steepness=STEEPNESS_A)
        assert isinstance(out, float)

    def test_array_returns_ndarray_shape_parity(self):
        depth = np.linspace(-20.0, 20.0, 17)
        out = rho_he_ratio(depth, steepness=STEEPNESS_A)
        assert isinstance(out, np.ndarray)
        assert out.shape == depth.shape


class TestSpatialGateParity:
    """Single-source regression: rho and drag route through one helper."""

    def test_parity_to_machine_precision_on_grid(self):
        depth = np.linspace(-50.0, 50.0, 201)
        r = rho_he_ratio(depth, steepness=STEEPNESS_A)
        g = spatial_gate(depth, STEEPNESS_A)
        np.testing.assert_array_equal(np.asarray(r), np.asarray(g))

    def test_spatial_gate_routes_through__erf_complement(self):
        depth = np.linspace(-50.0, 50.0, 201)
        np.testing.assert_array_equal(
            np.asarray(spatial_gate(depth, STEEPNESS_A)),
            np.asarray(_erf_complement(depth, STEEPNESS_A)),
        )


class TestTabulatedProfile:
    """The declared sourced-profile fallback machinery (data deferred, rule-2)."""

    def test_round_trips_grid_nodes(self):
        # depth increasing -> ratio decreasing (1 inside, 0 outside).
        depth_grid = np.array([-30.0, -10.0, 0.0, 10.0, 30.0])
        ratio_grid = np.array([1.0, 0.8, 0.5, 0.2, 0.0])
        prof = tabulated_density_profile(depth_grid, ratio_grid)
        out = prof.ratio(depth_grid)
        np.testing.assert_allclose(out, ratio_grid, rtol=0, atol=0)

    def test_linear_interpolation_between_nodes(self):
        depth_grid = np.array([-10.0, 10.0])
        ratio_grid = np.array([1.0, 0.0])
        prof = tabulated_density_profile(depth_grid, ratio_grid)
        assert prof.ratio(0.0) == pytest.approx(0.5)

    def test_clamps_to_tail_values_outside_grid(self):
        depth_grid = np.array([-10.0, 0.0, 10.0])
        ratio_grid = np.array([1.0, 0.5, 0.0])
        prof = tabulated_density_profile(depth_grid, ratio_grid)
        assert prof.ratio(-1.0e3) == 1.0   # deep inside -> tail value
        assert prof.ratio(1.0e3) == 0.0    # far outside -> tail value

    def test_scalar_returns_float(self):
        prof = tabulated_density_profile(
            np.array([-10.0, 10.0]), np.array([1.0, 0.0])
        )
        assert isinstance(prof.ratio(0.0), float)

    def test_rejects_non_increasing_depth_grid(self):
        with pytest.raises(ValueError):
            tabulated_density_profile(
                np.array([10.0, -10.0]), np.array([0.0, 1.0])
            )

    def test_rejects_ratio_outside_unit_interval(self):
        with pytest.raises(ValueError):
            tabulated_density_profile(
                np.array([-10.0, 10.0]), np.array([1.2, 0.0])
            )

    def test_rejects_length_mismatch(self):
        with pytest.raises(ValueError):
            tabulated_density_profile(
                np.array([-10.0, 0.0, 10.0]), np.array([1.0, 0.0])
            )


class TestTabulatedProfileDirectConstruction:
    """Validation must fire at *construction*, not only in the builder -- a bad
    ``depth_grid`` passed straight to the frozen dataclass would otherwise make
    ``np.interp`` return silently-wrong values (worse than a crash; principle 4).
    """

    def test_direct_construction_rejects_non_increasing_depth(self):
        with pytest.raises(ValueError):
            TabulatedDensityProfile(depth_grid=(10.0, -10.0), ratio_grid=(0.0, 1.0))

    def test_direct_construction_rejects_ratio_outside_unit_interval(self):
        with pytest.raises(ValueError):
            TabulatedDensityProfile(depth_grid=(-10.0, 10.0), ratio_grid=(1.5, 0.0))

    def test_direct_construction_rejects_length_mismatch(self):
        with pytest.raises(ValueError):
            TabulatedDensityProfile(depth_grid=(-10.0, 0.0, 10.0), ratio_grid=(1.0, 0.0))

    def test_direct_construction_rejects_empty(self):
        with pytest.raises(ValueError):
            TabulatedDensityProfile(depth_grid=(), ratio_grid=())


class TestTabulatedProfileDeliberateAllowances:
    """Behaviours intentionally *not* rejected -- locked so a later 'tightening'
    is a conscious choice, not silent drift.
    """

    def test_array_input_returns_ndarray(self):
        prof = tabulated_density_profile(
            np.array([-10.0, 0.0, 10.0]), np.array([1.0, 0.5, 0.0])
        )
        out = prof.ratio(np.array([-5.0, 5.0]))
        assert isinstance(out, np.ndarray)
        assert out.shape == (2,)

    def test_non_monotone_ratio_grid_is_accepted(self):
        # A sourced TDDFT profile may show a near-surface density lip; only the
        # *depth* axis must be sorted (np.interp requirement), not the ratio.
        prof = tabulated_density_profile(
            np.array([-20.0, -10.0, 0.0, 10.0]),
            np.array([0.9, 1.0, 0.6, 0.0]),
        )
        assert prof.ratio(-10.0) == pytest.approx(1.0)

    def test_single_node_profile_is_constant(self):
        # A one-node table is degenerate but valid (matches tabulated_ladder's
        # >=1 contract); np.interp returns the lone value everywhere.
        prof = tabulated_density_profile(np.array([0.0]), np.array([0.5]))
        assert prof.ratio(-100.0) == 0.5
        assert prof.ratio(100.0) == 0.5


class TestErfComplementNonFiniteSteepness:
    """NaN / non-positive steepness is rejected by the single guard (both callers)."""

    def test_nan_steepness_rejected_in_helper(self):
        with pytest.raises(ValueError):
            _erf_complement(0.0, float("nan"))

    def test_nan_steepness_rejected_in_rho_ratio(self):
        with pytest.raises(ValueError):
            rho_he_ratio(0.0, steepness=float("nan"))


class TestRhoRatioInputForms:
    """Scalar-return discipline across python-scalar / 0-d / list / array inputs."""

    def test_python_int_returns_float(self):
        assert isinstance(rho_he_ratio(0, steepness=STEEPNESS_A), float)

    def test_zero_d_array_returns_float(self):
        assert isinstance(rho_he_ratio(np.float64(3.0), steepness=STEEPNESS_A), float)

    def test_list_input_returns_ndarray(self):
        out = rho_he_ratio([-10.0, 0.0, 10.0], steepness=STEEPNESS_A)
        assert isinstance(out, np.ndarray)
        assert out.shape == (3,)
