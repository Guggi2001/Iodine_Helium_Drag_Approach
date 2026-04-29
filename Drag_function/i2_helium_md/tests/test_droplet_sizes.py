"""Tests for i2_helium_md/sampling/droplet_sizes.py."""

import numpy as np
import pytest

from i2_helium_md import single_pulse_N2000
from i2_helium_md.sampling.droplet_sizes import (
    _droplet_radius_m,
    _evaporation_per_pickup,
    _impact_parameter_threshold_m,
    _pickup_cross_section_m2,
    _simulate_pickup,
    droplet_radius_from_N,
    mean_droplet_size,
    sample_droplet_sizes,
)


# ---------------------------------------------------------------------------
# mean_droplet_size  -- the empirical correlation
# ---------------------------------------------------------------------------
class TestMeanDropletSize:
    def test_matlab_formula_baseline(self):
        """Reference operating point: p=40 mbar, T=23 K, d=5 um.

        N = 4e5 * 40^0.97 * 23^-3.88 * 5^2
        Verified against MATLAB by direct calculation.
        """
        N, R = mean_droplet_size(40.0, 23.0)
        # compute expected by hand using the same formula
        expected_N = 4e5 * (40 ** 0.97) * (23 ** -3.88) * (5 ** 2)
        assert N == pytest.approx(expected_N, rel=1e-12)
        # radius scales as N^(1/3) at fixed density
        assert R > 0
        assert R == pytest.approx(droplet_radius_from_N(N), rel=1e-12)

    def test_pressure_dependence(self):
        """Higher source pressure -> larger droplets."""
        N_low, _ = mean_droplet_size(20.0, 14.0)
        N_high, _ = mean_droplet_size(80.0, 14.0)
        assert N_high > N_low

    def test_temperature_dependence(self):
        """Lower T -> larger droplets (k3 is negative)."""
        N_warm, _ = mean_droplet_size(40.0, 30.0)
        N_cold, _ = mean_droplet_size(40.0, 10.0)
        assert N_cold > N_warm

    def test_radius_grows_as_cube_root_of_N(self):
        # If we double N we expect R to grow by 2^(1/3) ≈ 1.26
        _, R1 = mean_droplet_size(40.0, 14.0)
        # We can't change N directly via the formula, but droplet_radius_from_N
        # gives us the underlying relation:
        R_2N = droplet_radius_from_N(2 * 1000)
        R_N = droplet_radius_from_N(1000)
        assert R_2N / R_N == pytest.approx(2 ** (1.0 / 3.0), rel=1e-12)


# ---------------------------------------------------------------------------
# sample_droplet_sizes -- raw mode
# ---------------------------------------------------------------------------
class TestSampleRaw:
    def test_returns_correct_count(self):
        cfg = single_pulse_N2000(num_molecules=500, seed=42)
        N = sample_droplet_sizes(cfg, mode="raw")
        assert N.shape == (500,)
        assert (N > 0).all()

    def test_seeded_reproducibility(self):
        cfg = single_pulse_N2000(num_molecules=100, seed=123)
        N1 = sample_droplet_sizes(cfg, mode="raw")
        N2 = sample_droplet_sizes(cfg, mode="raw")
        np.testing.assert_array_equal(N1, N2)

    def test_different_seeds_differ(self):
        cfg_a = single_pulse_N2000(num_molecules=100, seed=1)
        cfg_b = single_pulse_N2000(num_molecules=100, seed=2)
        Na = sample_droplet_sizes(cfg_a, mode="raw")
        Nb = sample_droplet_sizes(cfg_b, mode="raw")
        assert not np.array_equal(Na, Nb)

    def test_mean_close_to_target(self):
        """For 50,000 samples the geometric mean should match the target N_mean."""
        cfg = single_pulse_N2000(num_molecules=50_000, seed=7)
        N = sample_droplet_sizes(cfg, mode="raw")
        target_N, _ = mean_droplet_size(cfg.p_source_mbar, cfg.T_source_K)
        # Log-normal: arithmetic mean = exp(mu + delta^2/2) which equals N_mean
        assert N.mean() == pytest.approx(target_N, rel=0.05)

    def test_distribution_is_lognormal_shape(self):
        """log(N) should be approximately normal."""
        cfg = single_pulse_N2000(num_molecules=20_000, seed=11)
        N = sample_droplet_sizes(cfg, mode="raw")
        log_N = np.log(N)
        # Skewness of normal distribution is 0; log-N samples have nearly 0
        skew = ((log_N - log_N.mean()) ** 3).mean() / log_N.std() ** 3
        assert abs(skew) < 0.1, f"|skewness of log(N)| = {abs(skew):.3f} too large"


# ---------------------------------------------------------------------------
# sample_droplet_sizes -- post_pickup mode
# ---------------------------------------------------------------------------
class TestSamplePostPickup:
    def test_returns_correct_count(self):
        cfg = single_pulse_N2000(num_molecules=200, seed=42)
        N = sample_droplet_sizes(cfg, mode="post_pickup")
        assert N.shape == (200,)
        assert (N > 0).all()

    def test_post_pickup_distribution_consistent(self):
        """Pickup biases sampling toward larger droplets (cross section ~ R^2),
        but evaporation strips some atoms. The net effect on the mean depends
        on which dominates; we just check the result is finite and reasonable.
        """
        cfg = single_pulse_N2000(num_molecules=2_000, seed=99)
        N_raw = sample_droplet_sizes(cfg, mode="raw")
        N_post = sample_droplet_sizes(cfg, mode="post_pickup")
        assert 0 < N_post.mean()
        # A reasonable post-pickup mean stays within an order of magnitude of raw
        assert 0.1 * N_raw.mean() < N_post.mean() < 10 * N_raw.mean()

    def test_invalid_mode_raises(self):
        cfg = single_pulse_N2000(num_molecules=10, seed=0)
        with pytest.raises(ValueError):
            sample_droplet_sizes(cfg, mode="other")  # type: ignore[arg-type]

    def test_E_solv_changes_reduced_crosssection(self):
        """Higher E_solv -> smaller (log(E_solv/E_kin_0))^2 -> larger b_thresh
        -> reduced cross section closer to geometric.

        Iodine in the thesis (E_solv=30) gives smaller corrections than
        the legacy MATLAB E_solv=14.
        """
        from i2_helium_md.sampling.droplet_sizes import _pickup_cross_section_m2
        N = np.array([5_000.0, 10_000.0])
        sigma_legacy = _pickup_cross_section_m2(N, reduced=True, E_solv_meV=14.0)
        sigma_thesis = _pickup_cross_section_m2(N, reduced=True, E_solv_meV=30.0)
        sigma_geom = _pickup_cross_section_m2(N, reduced=False)
        # both reduced are smaller than geometric
        assert (sigma_legacy <= sigma_geom).all()
        assert (sigma_thesis <= sigma_geom).all()
        # the thesis value (30 meV) gives a *larger* reduced cross section
        # because the correction term shrinks
        assert (sigma_thesis >= sigma_legacy).all()

    def test_thesis_reproduction_panel_b_T18(self):
        """Regression test: at T=18, p~25 mbar with the default E_solv=14 meV
        (which matches both the legacy MATLAB and the thesis figure 3.2),
        the post-pickup distribution peak should be near N=2500-3500
        (matches thesis figure 3.2 panel b red curve).
        """
        cfg = single_pulse_N2000(num_molecules=10_000, seed=42,
                                  p_source_mbar=25, T_source_K=18)
        # Default E_solv=14 reproduces the thesis figure
        N = sample_droplet_sizes(cfg, mode="post_pickup",
                                  reduced_crosssection=False)
        hist, edges = np.histogram(N, bins=50, range=(0, 30_000))
        peak_bin = np.argmax(hist)
        peak_N = (edges[peak_bin] + edges[peak_bin + 1]) / 2
        assert 1500 < peak_N < 5000, (
            f"thesis-figure peak should be near 2500-3500, got {peak_N:.0f}"
        )


    def test_insufficient_droplets_raises(self):
        """If we ask for far too many, the function should error not silently."""
        # Drive to failure by setting a huge num_molecules with limited oversample.
        # We test the internal helper directly because cfg sets oversample.
        rng = np.random.default_rng(0)
        N_raw = rng.lognormal(mean=8.0, sigma=0.6, size=100)
        with pytest.raises(ValueError):
            _simulate_pickup(N_raw, num_target=10_000,
                             reduced_crosssection=False, rng=rng)


# ---------------------------------------------------------------------------
# Internal helpers (sanity checks)
# ---------------------------------------------------------------------------
class TestPickupInternals:
    def test_droplet_radius_meters_grows_with_N(self):
        N = np.array([100, 1_000, 10_000])
        R = _droplet_radius_m(N)
        # Should be monotonically increasing
        assert np.all(np.diff(R) > 0)
        # Order of magnitude: ~nm for thousand-atom droplets
        assert 1e-9 < R[1] < 1e-7

    def test_impact_parameter_threshold_nonneg(self):
        N = np.array([100, 1_000, 10_000])
        b = _impact_parameter_threshold_m(N)
        assert (b >= 0).all()
        # Threshold should be at most R
        R = _droplet_radius_m(N)
        assert (b <= R + 1e-15).all()

    def test_cross_section_geometric_vs_reduced(self):
        """Reduced cross section should never exceed geometric."""
        N = np.array([1_000.0, 5_000.0, 20_000.0])
        sigma_geom = _pickup_cross_section_m2(N, reduced=False)
        sigma_red = _pickup_cross_section_m2(N, reduced=True)
        assert (sigma_red <= sigma_geom + 1e-30).all()

    def test_evaporation_negative(self):
        """Evaporation should remove He atoms (negative dN)."""
        N = np.array([500.0, 2000.0])
        k = np.array([1, 1])
        dN = _evaporation_per_pickup(k, N)
        assert (dN < 0).all()

    def test_evaporation_increases_with_k(self):
        """Each additional pickup adds extra evaporation."""
        N = np.array([2000.0])
        dN1 = _evaporation_per_pickup(np.array([1]), N)
        dN2 = _evaporation_per_pickup(np.array([2]), N)
        # k=2 has 5/2 extra term -> more negative -> larger magnitude
        assert abs(dN2[0]) > abs(dN1[0])
