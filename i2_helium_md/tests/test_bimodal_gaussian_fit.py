"""Tests for the bimodal Gaussian fit added to velocity_distribution.py."""

from __future__ import annotations

import numpy as np

from i2_helium_md.postprocess.velocity_distribution import (
    BimodalGaussianFit,
    FinalVelocityHistogram,
    bimodal_gaussian_fit,
)


def _hist_from_density(v: np.ndarray, density: np.ndarray) -> FinalVelocityHistogram:
    edges = np.concatenate(([v[0] - 0.5 * (v[1] - v[0])], 0.5 * (v[1:] + v[:-1]),
                             [v[-1] + 0.5 * (v[-1] - v[-2])]))
    return FinalVelocityHistogram(
        bin_centers_Aps=v,
        bin_centers_mps=v * 100.0,
        bin_edges_Aps=edges,
        bin_edges_mps=edges * 100.0,
        counts=density.astype(int),
        density=density,
        mass_amu=131.0,
        num_atoms_used=int(density.sum()),
    )


class TestBimodalGaussianFit:
    def test_recovers_two_peaks(self):
        v = np.linspace(0.0, 30.0, 200)
        true_a1, true_mu1, true_s1 = 1.0, 8.0, 1.5
        true_a2, true_mu2, true_s2 = 0.6, 18.0, 2.0
        density = (
            true_a1 * np.exp(-((v - true_mu1) ** 2) / (2 * true_s1 ** 2))
            + true_a2 * np.exp(-((v - true_mu2) ** 2) / (2 * true_s2 ** 2))
        )
        # Add a tiny floor to keep histogram counts positive in the helper.
        hist = _hist_from_density(v, density + 1e-3)
        fit = bimodal_gaussian_fit(hist)
        assert isinstance(fit, BimodalGaussianFit)
        assert fit.success
        # Identify peaks regardless of the order returned.
        means = sorted([fit.mean_1_Aps, fit.mean_2_Aps])
        assert abs(means[0] - true_mu1) < 1.0
        assert abs(means[1] - true_mu2) < 1.0
        sigmas = sorted([fit.sigma_1_Aps, fit.sigma_2_Aps])
        assert abs(sigmas[0] - true_s1) < 0.5
        assert abs(sigmas[1] - true_s2) < 0.5

    def test_flat_input_fails_gracefully(self):
        v = np.linspace(0.0, 30.0, 50)
        density = np.zeros_like(v)
        density[5] = 1.0  # only one non-zero bin
        hist = _hist_from_density(v, density)
        fit = bimodal_gaussian_fit(hist)
        assert not fit.success
        assert np.isnan(fit.mean_1_Aps)
        assert np.isnan(fit.amplitude_1)
