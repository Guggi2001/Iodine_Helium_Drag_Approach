"""Tests for i2_helium_md/sampling/droplet_sizes_diagnostics.py."""

import numpy as np
import pytest

# matplotlib must be in non-interactive mode for pytest
import matplotlib
matplotlib.use("Agg")

from i2_helium_md import single_pulse_N2000
from i2_helium_md.sampling.droplet_sizes_diagnostics import diagnose_pickup


class TestDiagnosePickup:
    def test_returns_dict_and_figure(self):
        cfg = single_pulse_N2000(num_molecules=200, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diag, fig = diagnose_pickup(cfg, plot=True, print_report=False)
        assert isinstance(diag, dict)
        assert fig is not None

    def test_returns_none_fig_when_plot_false(self):
        cfg = single_pulse_N2000(num_molecules=200, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diag, fig = diagnose_pickup(cfg, plot=False, print_report=False)
        assert fig is None

    def test_diagnostics_has_expected_keys(self):
        cfg = single_pulse_N2000(num_molecules=200, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diag, _ = diagnose_pickup(cfg, plot=False, print_report=False)
        for key in ("N_raw", "N_one_pickup", "N_no_pickup",
                    "per_round", "source_conditions"):
            assert key in diag, f"missing key: {key}"

    def test_per_round_has_required_fields(self):
        cfg = single_pulse_N2000(num_molecules=200, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diag, _ = diagnose_pickup(cfg, plot=False, print_report=False)
        # at least one round must have been recorded
        assert len(diag["per_round"]) >= 1
        for r in diag["per_round"]:
            for field in ("iter", "n_alive_before", "mean_p_pickup",
                          "mean_evap", "n_picked_up"):
                assert field in r, f"per_round entry missing {field}"

    def test_no_pickup_plus_one_pickup_le_raw_total(self):
        """Conservation: no_pickup + one_pickup droplets should not exceed raw total.

        (Some samples can be destroyed via evaporation, hence inequality not equality.)
        """
        cfg = single_pulse_N2000(num_molecules=200, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diag, _ = diagnose_pickup(cfg, plot=False, print_report=False)
        n_raw = diag["N_raw"].size
        n_zero = diag["N_no_pickup"].size
        n_one = diag["N_one_pickup"].size
        assert n_zero + n_one <= n_raw

    def test_source_conditions_round_trip(self):
        """The source_conditions dict should reflect what we passed in."""
        cfg = single_pulse_N2000(num_molecules=200, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diag, _ = diagnose_pickup(
            cfg, reduced_crosssection=True, E_solv_meV=30.0,
            plot=False, print_report=False,
        )
        src = diag["source_conditions"]
        assert src["p_source_mbar"] == 25
        assert src["T_source_K"] == 15
        assert src["reduced_crosssection"] is True
        assert src["E_solv_meV"] == 30.0

    def test_reduced_vs_normal_yield_differs(self):
        """Reduced cross section should give a different one-pickup yield
        than the normal cross section (by virtue of pruning small droplets)."""
        cfg = single_pulse_N2000(num_molecules=500, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diag_n, _ = diagnose_pickup(cfg, reduced_crosssection=False,
                                    plot=False, print_report=False)
        diag_r, _ = diagnose_pickup(cfg, reduced_crosssection=True,
                                    plot=False, print_report=False)
        assert diag_n["N_one_pickup"].size != diag_r["N_one_pickup"].size

    def test_print_report_does_not_crash(self, capsys):
        """Smoke test: print_report=True writes something to stdout."""
        cfg = single_pulse_N2000(num_molecules=200, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diagnose_pickup(cfg, plot=False, print_report=True)
        captured = capsys.readouterr()
        assert "Pickup-cell simulation diagnostics" in captured.out
        assert "round" in captured.out

    def test_per_round_includes_new_samples(self):
        """Each per_round entry should carry the post-pickup new_samples array
        (used by Panel 4 of the diagnostic figure)."""
        cfg = single_pulse_N2000(num_molecules=200, seed=42,
                                  p_source_mbar=25, T_source_K=15)
        diag, _ = diagnose_pickup(cfg, plot=False, print_report=False)
        assert len(diag["per_round"]) >= 1
        for r in diag["per_round"]:
            assert "new_samples" in r
            assert isinstance(r["new_samples"], np.ndarray)


class TestBayesHist:
    """Tests for the Bayesian-smoothed histogram helper."""

    def test_normalized_density(self):
        """Total integral should be ~1, modulo Laplace smoothing.
        With smoothing, sum(p_i) = sum((n_i+1)/(N+nbins+1)) = (N+nbins)/(N+nbins+1)
        which approaches 1 as N grows."""
        from i2_helium_md.sampling.droplet_sizes_diagnostics import _bayes_hist
        rng = np.random.default_rng(0)
        samples = rng.lognormal(np.log(10000), 0.5, size=5000)
        bins = np.linspace(1, samples.max(), 60)
        centers, h, sigma_h = _bayes_hist(samples, bins)
        # sum p_i ~ N/(N+nbins+1) which is close to 1 for large N
        binwidth = bins[1] - bins[0]
        total_p = np.sum(h * binwidth)
        assert 0.95 < total_p < 1.0

    def test_no_zero_bins(self):
        """Laplace smoothing should ensure no bin has p == 0."""
        from i2_helium_md.sampling.droplet_sizes_diagnostics import _bayes_hist
        # use a tight sample that doesn't fill all bins
        samples = np.array([5.0, 5.1, 5.2])
        bins = np.linspace(0, 10, 20)
        _, h, _ = _bayes_hist(samples, bins)
        assert np.all(h > 0), "Bayesian smoothing should keep all bins nonzero"

    def test_sigma_nonnegative(self):
        from i2_helium_md.sampling.droplet_sizes_diagnostics import _bayes_hist
        samples = np.linspace(0, 10, 100)
        bins = np.linspace(0, 10, 20)
        _, _, sigma_h = _bayes_hist(samples, bins)
        assert np.all(sigma_h >= 0)

    def test_centers_correct_length(self):
        from i2_helium_md.sampling.droplet_sizes_diagnostics import _bayes_hist
        samples = np.linspace(0, 10, 100)
        bins = np.linspace(0, 10, 20)
        centers, h, sigma_h = _bayes_hist(samples, bins)
        assert centers.size == bins.size - 1
        assert h.size == centers.size
        assert sigma_h.size == centers.size


class TestAnalyticalFormula:
    """Tests for conditional_size_distributions_analytical (Treiber's formula)."""

    def test_returns_two_arrays(self):
        from i2_helium_md.sampling.droplet_sizes import (
            conditional_size_distributions_analytical,
        )
        N_grid = np.arange(1, 100_001, dtype=float)
        p_n, p_r = conditional_size_distributions_analytical(
            N_grid, p_source_mbar=40.0, T_source_K=15.0,
        )
        assert p_n.shape == N_grid.shape
        assert p_r.shape == N_grid.shape

    def test_normal_normalized_to_one(self):
        """The normal-σ conditional integrates to 1."""
        from i2_helium_md.sampling.droplet_sizes import (
            conditional_size_distributions_analytical,
        )
        N_grid = np.arange(1, 100_001, dtype=float)
        p_n, _ = conditional_size_distributions_analytical(
            N_grid, p_source_mbar=40.0, T_source_K=15.0,
        )
        assert np.trapezoid(p_n, N_grid) == pytest.approx(1.0, rel=1e-6)

    def test_reduced_integrates_to_less_than_one(self):
        """In Treiber's normalisation, reduced-σ integrates to <1
        (it represents the fraction of one-pickup events that pass the
        kinetic-energy threshold)."""
        from i2_helium_md.sampling.droplet_sizes import (
            conditional_size_distributions_analytical,
        )
        N_grid = np.arange(1, 100_001, dtype=float)
        _, p_r = conditional_size_distributions_analytical(
            N_grid, p_source_mbar=40.0, T_source_K=15.0,
        )
        I = np.trapezoid(p_r, N_grid)
        assert 0 < I < 1.0

    def test_thesis_T18_normal_peak(self):
        """Regression: at T=18, p=40 mbar, normal σ peak should be near N=2500
        (matches thesis figure 3.2)."""
        from i2_helium_md.sampling.droplet_sizes import (
            conditional_size_distributions_analytical,
        )
        N_grid = np.arange(1, 100_001, dtype=float)
        p_n, _ = conditional_size_distributions_analytical(
            N_grid, p_source_mbar=40.0, T_source_K=18.0,
        )
        peak_N = N_grid[np.argmax(p_n)]
        assert 1500 < peak_N < 4000, f"peak should be near 2500, got {peak_N}"

    def test_thesis_T12_reduced_peak_above_normal(self):
        """Thesis-figure signature: at T=12, the reduced-σ peak (dashed) is
        TALLER than the normal-σ peak (solid)."""
        from i2_helium_md.sampling.droplet_sizes import (
            conditional_size_distributions_analytical,
        )
        N_grid = np.arange(1, 100_001, dtype=float)
        p_n, p_r = conditional_size_distributions_analytical(
            N_grid, p_source_mbar=40.0, T_source_K=12.0,
        )
        assert p_r.max() > p_n.max(), (
            f"thesis: at T=12, dashed peak ({p_r.max():.2e}) "
            f"should exceed solid peak ({p_n.max():.2e})"
        )

    def test_thesis_T18_reduced_peak_below_normal(self):
        """Mirror: at T=18, the reduced-σ peak (dashed) is MUCH shorter
        than the normal-σ peak (solid)."""
        from i2_helium_md.sampling.droplet_sizes import (
            conditional_size_distributions_analytical,
        )
        N_grid = np.arange(1, 100_001, dtype=float)
        p_n, p_r = conditional_size_distributions_analytical(
            N_grid, p_source_mbar=40.0, T_source_K=18.0,
        )
        assert p_r.max() < p_n.max() * 0.4, (
            f"thesis: at T=18, dashed peak ({p_r.max():.2e}) "
            f"should be much shorter than solid ({p_n.max():.2e})"
        )


class TestThesisFigureReproduction:
    def test_returns_figure(self):
        """plot_thesis_figure_3_2 should return a matplotlib Figure with 2 axes."""
        from i2_helium_md.sampling.droplet_sizes_diagnostics import (
            plot_thesis_figure_3_2,
        )
        fig = plot_thesis_figure_3_2()
        assert fig is not None
        assert len(fig.axes) == 2

    def test_custom_temperatures(self):
        """Should accept a custom temperature tuple."""
        from i2_helium_md.sampling.droplet_sizes_diagnostics import (
            plot_thesis_figure_3_2,
        )
        fig = plot_thesis_figure_3_2(temperatures_K=(10.0, 14.0, 20.0))
        assert fig is not None
