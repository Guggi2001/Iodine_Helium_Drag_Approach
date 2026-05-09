"""Tests for i2_helium_md/postprocess/boltzmann_overlay.py."""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.postprocess.boltzmann_overlay import (
    BoltzmannCurve,
    boltzmann_population,
)


class TestBoltzmannPopulation:
    def test_density_normalises_to_unity(self):
        curve = boltzmann_population(
            droplet_radius_A=9.0,
            temperature_K=0.4,
            steepness_A=1.0,
            binding_energy_eV=0.05,
            n_points=400,
        )
        assert isinstance(curve, BoltzmannCurve)
        integral = float(np.trapezoid(curve.density, curve.r_grid_A))
        assert integral == pytest.approx(1.0, rel=1e-6)

    def test_unnormalised_starts_high_inside(self):
        # Inside the droplet (r << R) V ~ 0 -> exp(0) = 1 -> high.
        # Outside (r >> R) V -> binding_energy > 0 -> exp(-V/kT) -> small.
        curve = boltzmann_population(
            droplet_radius_A=9.0,
            temperature_K=0.4,
            steepness_A=1.0,
            binding_energy_eV=0.1,
        )
        # Find a sample well inside and well outside the droplet.
        idx_inside = int(np.searchsorted(curve.r_grid_A, 1.0))
        idx_outside = int(np.searchsorted(curve.r_grid_A, 17.0))
        assert curve.unnormalised[idx_inside] > curve.unnormalised[idx_outside]

    def test_explicit_grid_used(self):
        grid = np.linspace(0.0, 12.0, 50)
        curve = boltzmann_population(
            droplet_radius_A=9.0,
            temperature_K=1.0,
            steepness_A=0.5,
            binding_energy_eV=0.05,
            r_grid_A=grid,
        )
        np.testing.assert_array_equal(curve.r_grid_A, grid)
        assert curve.density.size == grid.size

    def test_bad_temperature(self):
        with pytest.raises(ValueError, match="temperature_K"):
            boltzmann_population(
                droplet_radius_A=9.0, temperature_K=0.0,
                steepness_A=1.0, binding_energy_eV=0.05,
            )

    def test_bad_radius(self):
        with pytest.raises(ValueError, match="droplet_radius_A"):
            boltzmann_population(
                droplet_radius_A=0.0, temperature_K=0.4,
                steepness_A=1.0, binding_energy_eV=0.05,
            )
