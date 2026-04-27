"""Tests for i2_helium_md/sampling/radial_positions.py."""

import numpy as np
import pytest

from i2_helium_md import single_pulse_N2000
from i2_helium_md.physics.constants import EV, K_B
from i2_helium_md.sampling.radial_positions import (
    _radial_probability_density,
    _sample_radial_for_one_droplet,
    sample_radial_positions,
)


# ===========================================================================
# _radial_probability_density: the Boltzmann density itself
# ===========================================================================
class TestRadialProbabilityDensity:
    def test_zero_at_origin(self):
        """p(0) = 0 because of the r^2 volume element."""
        p = _radial_probability_density(
            r=np.array([0.0]),
            droplet_radius=30.0,
            T_K=0.4,
            steepness=14.3,
            binding_eV=0.025,
        )
        assert p[0] == 0.0

    def test_positive_inside_droplet(self):
        """For r inside the droplet, p > 0."""
        r = np.array([5.0, 10.0, 20.0])
        p = _radial_probability_density(
            r=r,
            droplet_radius=30.0,
            T_K=0.4,
            steepness=14.3,
            binding_eV=0.025,
        )
        assert (p > 0).all()

    def test_decays_far_outside_droplet(self):
        """For r >> R, the Boltzmann factor exp(-U/kT) suppresses p (U -> binding_eV).
        We expect p to be smaller than its peak, not strictly monotone in r.
        """
        R = 30.0
        # peak is somewhere near R
        p_peak = _radial_probability_density(
            r=np.array([R]),
            droplet_radius=R,
            T_K=0.4,
            steepness=14.3,
            binding_eV=0.025,
        )[0]
        # far outside, p should be much smaller
        p_far = _radial_probability_density(
            r=np.array([2 * R + 10.0]),
            droplet_radius=R,
            T_K=0.4,
            steepness=14.3,
            binding_eV=0.025,
        )[0]
        # because U(r >> R) ~ binding_eV, exp(-binding/kT) at T=0.4 K is ~0
        assert p_far < p_peak * 0.1

    def test_temperature_dependence(self):
        """At very high T, the Boltzmann factor is ~1 everywhere -> distribution
        is dominated by the r^2 Jacobian."""
        r = np.array([5.0, 25.0, 50.0])
        # very hot
        p_hot = _radial_probability_density(
            r=r, droplet_radius=30.0, T_K=1e6, steepness=14.3, binding_eV=0.025,
        )
        # ratios should look like r^2 ratios
        ratio_25_to_5 = p_hot[1] / p_hot[0]
        assert ratio_25_to_5 == pytest.approx((25 / 5) ** 2, rel=0.01)


# ===========================================================================
# _sample_radial_for_one_droplet
# ===========================================================================
class TestSingleDropletSampler:
    def test_returns_correct_count(self):
        rng = np.random.default_rng(42)
        r = _sample_radial_for_one_droplet(
            droplet_radius=30.0, n_samples=500, T_K=0.4,
            steepness=14.3, binding_eV=0.025, r_step=0.05, rng=rng,
        )
        assert r.shape == (500,)

    def test_samples_within_bounds(self):
        rng = np.random.default_rng(42)
        R = 30.0
        r = _sample_radial_for_one_droplet(
            droplet_radius=R, n_samples=1000, T_K=0.4,
            steepness=14.3, binding_eV=0.025, r_step=0.05, rng=rng,
        )
        # all samples must lie in [0, 2R)
        assert (r >= 0).all()
        assert (r < 2 * R).all()

    def test_seeded_reproducibility(self):
        kwargs = dict(
            droplet_radius=30.0, n_samples=200, T_K=0.4,
            steepness=14.3, binding_eV=0.025, r_step=0.05,
        )
        rng1 = np.random.default_rng(7)
        rng2 = np.random.default_rng(7)
        r1 = _sample_radial_for_one_droplet(**kwargs, rng=rng1)
        r2 = _sample_radial_for_one_droplet(**kwargs, rng=rng2)
        np.testing.assert_array_equal(r1, r2)

    def test_distribution_peaks_inside_droplet(self):
        """For cold-T sampling, samples concentrate inside the droplet (r < R)."""
        rng = np.random.default_rng(0)
        R = 30.0
        r = _sample_radial_for_one_droplet(
            droplet_radius=R, n_samples=10_000, T_K=0.4,
            steepness=14.3, binding_eV=0.025, r_step=0.05, rng=rng,
        )
        # at low T the molecule strongly prefers being inside; >90% of mass < R
        frac_inside = float((r < R).sum()) / r.size
        assert frac_inside > 0.9, f"only {frac_inside*100:.1f}% inside droplet"


# ===========================================================================
# sample_radial_positions: the public entry point
# ===========================================================================
class TestSampleRadialPositions:
    def test_returns_one_per_molecule(self):
        cfg = single_pulse_N2000(num_molecules=200, seed=42)
        droplet_radii = np.full(200, 30.0)
        r = sample_radial_positions(cfg, droplet_radii)
        assert r.shape == (200,)

    def test_handles_heterogeneous_radii(self):
        """Different droplet sizes should each get their own samples
        in the right slot."""
        cfg = single_pulse_N2000(num_molecules=300, seed=42)
        # Mix of three radii
        droplet_radii = np.concatenate([
            np.full(100, 30.0),
            np.full(100, 40.0),
            np.full(100, 50.0),
        ])
        r = sample_radial_positions(cfg, droplet_radii)
        assert r.shape == (300,)
        # Slot 0..99 was generated with R=30, so all should be < 60
        assert (r[:100] < 60).all()
        assert (r[100:200] < 80).all()
        assert (r[200:] < 100).all()

    def test_seeded_reproducibility(self):
        cfg = single_pulse_N2000(num_molecules=50, seed=99)
        droplet_radii = np.full(50, 30.0)
        r1 = sample_radial_positions(cfg, droplet_radii)
        r2 = sample_radial_positions(cfg, droplet_radii)
        np.testing.assert_array_equal(r1, r2)

    def test_concentrates_inside_droplet_at_low_T(self):
        """At T=0.4 K with binding energy, molecules are mostly inside."""
        cfg = single_pulse_N2000(num_molecules=2000, seed=42, T_particles_K=0.4)
        droplet_radii = np.full(2000, 30.0)
        r = sample_radial_positions(cfg, droplet_radii)
        # at low T, most molecules are inside the droplet
        frac_inside = float((r < 30.0).sum()) / 2000
        assert frac_inside > 0.9

    def test_shape_with_2D_input_droplet_radii(self):
        """Pass droplet_radii as a (N,) array; ravel/flatten should be transparent."""
        cfg = single_pulse_N2000(num_molecules=50, seed=42)
        droplet_radii = np.full(50, 30.0).reshape(50, 1)  # 2D
        r = sample_radial_positions(cfg, droplet_radii)
        assert r.shape == (50,)
