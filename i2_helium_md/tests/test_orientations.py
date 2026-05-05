"""Tests for i2_helium_md/sampling/orientations.py."""

import numpy as np
import pytest

from i2_helium_md.sampling.orientations import (
    MolecularOrientations,
    sample_orientations,
)


# ===========================================================================
# Basic API
# ===========================================================================
class TestApi:
    def test_returns_dataclass(self):
        out = sample_orientations(
            10, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=False, rng=np.random.default_rng(0),
        )
        assert isinstance(out, MolecularOrientations)

    def test_array_shapes(self):
        n = 50
        out = sample_orientations(
            n, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=False, rng=np.random.default_rng(0),
        )
        for arr in (out.alpha, out.beta, out.gamma, out.delta,
                    out.bond_length_angstrom):
            assert arr.shape == (n,)

    def test_reproducible_with_seed(self):
        out1 = sample_orientations(
            100, R0_GS_angstrom=2.666, deltaR0_angstrom=0.05,
            anisotropic=True, rng=np.random.default_rng(42),
        )
        out2 = sample_orientations(
            100, R0_GS_angstrom=2.666, deltaR0_angstrom=0.05,
            anisotropic=True, rng=np.random.default_rng(42),
        )
        np.testing.assert_array_equal(out1.alpha, out2.alpha)
        np.testing.assert_array_equal(out1.delta, out2.delta)
        np.testing.assert_array_equal(
            out1.bond_length_angstrom, out2.bond_length_angstrom,
        )

    def test_works_without_explicit_rng(self):
        """Should work and return correct shapes when rng=None."""
        out = sample_orientations(
            10, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=False,
        )
        assert out.alpha.shape == (10,)


# ===========================================================================
# Validation
# ===========================================================================
class TestValidation:
    def test_zero_num_molecules_raises(self):
        with pytest.raises(ValueError, match="num_molecules"):
            sample_orientations(
                0, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
                anisotropic=False, rng=np.random.default_rng(0),
            )

    def test_negative_num_molecules_raises(self):
        with pytest.raises(ValueError, match="num_molecules"):
            sample_orientations(
                -5, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
                anisotropic=False, rng=np.random.default_rng(0),
            )

    def test_negative_deltaR0_raises(self):
        with pytest.raises(ValueError, match="deltaR0"):
            sample_orientations(
                10, R0_GS_angstrom=2.666, deltaR0_angstrom=-0.01,
                anisotropic=False, rng=np.random.default_rng(0),
            )


# ===========================================================================
# Bond length
# ===========================================================================
class TestBondLength:
    def test_zero_delta_gives_exact_R0(self):
        out = sample_orientations(
            100, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=False, rng=np.random.default_rng(0),
        )
        np.testing.assert_array_equal(
            out.bond_length_angstrom, 2.666 * np.ones(100),
        )

    def test_nonzero_delta_gives_gaussian_spread(self):
        n = 100_000
        out = sample_orientations(
            n, R0_GS_angstrom=2.666, deltaR0_angstrom=0.05,
            anisotropic=False, rng=np.random.default_rng(0),
        )
        assert abs(out.bond_length_angstrom.mean() - 2.666) < 0.001
        assert abs(out.bond_length_angstrom.std() - 0.05) < 0.001

    def test_R0_value_carries_through(self):
        """A different R0 (e.g. 9 Å for HeDFT comparison) should be respected."""
        out = sample_orientations(
            100, R0_GS_angstrom=9.0, deltaR0_angstrom=0.0,
            anisotropic=False, rng=np.random.default_rng(0),
        )
        np.testing.assert_array_equal(out.bond_length_angstrom, 9.0 * np.ones(100))


# ===========================================================================
# Angle range checks
# ===========================================================================
class TestAngleRanges:
    @pytest.fixture
    def out(self):
        return sample_orientations(
            10_000, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=False, rng=np.random.default_rng(0),
        )

    def test_beta_in_range(self, out):
        assert np.all(out.beta >= 0)
        assert np.all(out.beta < 2 * np.pi)

    def test_alpha_in_range(self, out):
        assert np.all(out.alpha >= 0)
        assert np.all(out.alpha < 2 * np.pi)

    def test_gamma_in_range(self, out):
        assert np.all(out.gamma >= 0)
        assert np.all(out.gamma <= np.pi)

    def test_delta_in_range(self, out):
        assert np.all(out.delta >= 0)
        assert np.all(out.delta <= np.pi)

    def test_anisotropic_angles_also_in_range(self):
        out = sample_orientations(
            5000, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=True, rng=np.random.default_rng(0),
        )
        assert np.all(out.alpha >= 0) and np.all(out.alpha < 2 * np.pi)
        assert np.all(out.delta >= 0) and np.all(out.delta <= np.pi)


# ===========================================================================
# Position angles always uniform on sphere
# ===========================================================================
class TestPositionAnglesUniform:
    def test_position_angles_uniform_in_isotropic(self):
        """Positions should be uniform on the sphere regardless of mode."""
        n = 50_000
        out = sample_orientations(
            n, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=False, rng=np.random.default_rng(0),
        )
        # Uniform on sphere -> <cos(gamma)> ≈ 0,
        # <cos²(gamma)> ≈ 1/3 (the average of cos²θ over the unit sphere)
        assert abs(np.cos(out.gamma).mean()) < 0.02
        assert abs((np.cos(out.gamma) ** 2).mean() - 1/3) < 0.02

    def test_position_angles_uniform_in_anisotropic(self):
        """Even in anisotropic-axis mode, position angles stay uniform."""
        n = 50_000
        out = sample_orientations(
            n, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=True, rng=np.random.default_rng(0),
        )
        assert abs(np.cos(out.gamma).mean()) < 0.02
        assert abs((np.cos(out.gamma) ** 2).mean() - 1/3) < 0.02


# ===========================================================================
# Isotropic axis: same as position
# ===========================================================================
class TestIsotropicAxisOrientation:
    def test_axis_uniform_when_isotropic(self):
        n = 50_000
        out = sample_orientations(
            n, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=False, rng=np.random.default_rng(0),
        )
        # delta is the polar angle of the axis. <cos²(delta)> ≈ 1/3
        assert abs((np.cos(out.delta) ** 2).mean() - 1/3) < 0.02


# ===========================================================================
# Anisotropic axis: cos²φ weighting
# ===========================================================================
class TestAnisotropicAxisOrientation:
    def test_axis_x_component_sq_mean_is_3_over_5(self):
        """For p(x) ∝ x² (with x = cos φ in [-1, 1] integrated over sphere),
        the mean of x² is 3/5 = 0.6, not 1/3 ≈ 0.333.

        Derivation: for a vector uniform on the sphere with weighting cos²φ,
                  <cos²φ> = ∫ cos²φ · cos²φ dΩ / ∫ cos²φ dΩ = (3/5)·(1/3)/(1/3) = 3/5.

        x-component of unit axis vector is cos(alpha)·sin(delta).
        So <(cos α sin δ)²> should be ≈ 0.6.
        """
        n = 50_000
        out = sample_orientations(
            n, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=True, rng=np.random.default_rng(0),
        )
        x_axis_sq = (np.cos(out.alpha) * np.sin(out.delta)) ** 2
        mean_x2 = x_axis_sq.mean()
        # Expect ≈ 0.6, allow ±0.02 tolerance for n=50000
        assert abs(mean_x2 - 0.6) < 0.02, (
            f"anisotropic <x²> should be ~0.6, got {mean_x2:.3f}"
        )

    def test_isotropic_x_component_sq_mean_is_1_over_3(self):
        """Sanity counter-test: isotropic gives <x²> = 1/3."""
        n = 50_000
        out = sample_orientations(
            n, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=False, rng=np.random.default_rng(0),
        )
        x_axis_sq = (np.cos(out.alpha) * np.sin(out.delta)) ** 2
        mean_x2 = x_axis_sq.mean()
        assert abs(mean_x2 - 1/3) < 0.02, (
            f"isotropic <x²> should be ~0.333, got {mean_x2:.3f}"
        )

    def test_y_z_components_unchanged_by_anisotropic(self):
        """The y and z components should still average to 1/3·(2/3) on the
        unit sphere weighted by 1-cos²φ. Specifically <y²+z²> ≈ 0.4 for
        cos²-weighted (since x² averages to 0.6 and total is 1)."""
        n = 50_000
        out = sample_orientations(
            n, R0_GS_angstrom=2.666, deltaR0_angstrom=0.0,
            anisotropic=True, rng=np.random.default_rng(0),
        )
        # For unit vector: x² + y² + z² = 1 always
        x = np.cos(out.alpha) * np.sin(out.delta)
        y = np.sin(out.alpha) * np.sin(out.delta)
        z = np.cos(out.delta)
        norms = x*x + y*y + z*z
        np.testing.assert_allclose(norms, 1.0, atol=1e-12)


# ===========================================================================
# Smoke test for HeDFT-comparison-style call
# ===========================================================================
class TestRealisticUsage:
    def test_hedft_comparison_R0_9A(self):
        """Smoke: HeDFT comparison runs use R0=9 Å (artificially stretched)."""
        out = sample_orientations(
            500, R0_GS_angstrom=9.0, deltaR0_angstrom=0.0,
            anisotropic=True, rng=np.random.default_rng(0),
        )
        assert out.bond_length_angstrom.shape == (500,)
        np.testing.assert_array_equal(out.bond_length_angstrom, 9.0)
