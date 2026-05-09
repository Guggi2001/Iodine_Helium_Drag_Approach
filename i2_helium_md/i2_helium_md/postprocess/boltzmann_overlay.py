"""Analytic Boltzmann reference curve over the existing droplet potential.

Reproduces the ``p_boltzmann = exp(-V/(k_B*T))`` overlay block of
``legacy_matlab_repository/post_process_compare_radial_distributions.m``.
The potential ``V(r)`` itself is the simulation's own
:func:`i2_helium_md.physics.potentials.droplet_potential`; this module
just normalises ``exp(-V / k_B T)`` onto a user-chosen radial grid so it
can be plotted on top of the initial-population histogram.

This is an analytic curve only -- no new physics -- and is consistent
with the Boltzmann sampling already performed in
``i2_helium_md/sampling``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..physics.constants import EV, K_B
from ..physics.potentials import droplet_potential


@dataclass(frozen=True)
class BoltzmannCurve:
    """Boltzmann population vs radial distance.

    ``density`` is normalised so that the trapezoidal integral over
    ``r_grid_A`` equals 1; ``unnormalised`` keeps the bare
    ``exp(-V / k_B T)`` values for diagnostics.
    """

    r_grid_A: np.ndarray
    density: np.ndarray
    unnormalised: np.ndarray
    droplet_radius_A: float
    temperature_K: float


def boltzmann_population(
    *,
    droplet_radius_A: float,
    temperature_K: float,
    steepness_A: float,
    binding_energy_eV: float,
    r_grid_A: np.ndarray | None = None,
    n_points: int = 400,
) -> BoltzmannCurve:
    """Compute ``exp(-V_droplet(r - R) / k_B T)`` and normalise on ``r_grid_A``.

    Parameters
    ----------
    droplet_radius_A
        Droplet radius R in angstrom; the potential argument is
        ``r - R`` so ``r = R`` is the surface.
    temperature_K
        Temperature in Kelvin used in the Boltzmann factor.
    steepness_A, binding_energy_eV
        Pass-through to :func:`droplet_potential` (``beta(1)`` and
        ``beta(2)`` in the legacy MATLAB code).
    r_grid_A
        Optional explicit radial grid. If ``None``, a uniform grid
        from 0 to ``2 * droplet_radius_A`` with ``n_points`` samples
        is used.
    n_points
        Number of grid points when ``r_grid_A`` is ``None``.
    """
    if temperature_K <= 0.0:
        raise ValueError(f"temperature_K must be > 0, got {temperature_K}")
    if droplet_radius_A <= 0.0:
        raise ValueError(
            f"droplet_radius_A must be > 0, got {droplet_radius_A}"
        )

    if r_grid_A is None:
        r_grid_A = np.linspace(0.0, 2.0 * droplet_radius_A, n_points)
    else:
        r_grid_A = np.asarray(r_grid_A, dtype=float)
        if r_grid_A.size < 2:
            raise ValueError("r_grid_A must have at least 2 samples")

    V_eV = droplet_potential(
        r_grid_A - droplet_radius_A,
        steepness=steepness_A,
        binding_energy=binding_energy_eV,
    )
    kT_eV = K_B * temperature_K / EV
    unnormalised = np.exp(-V_eV / kT_eV)

    integral = float(np.trapezoid(unnormalised, r_grid_A))
    if integral > 0.0:
        density = unnormalised / integral
    else:
        density = np.zeros_like(unnormalised)

    return BoltzmannCurve(
        r_grid_A=r_grid_A,
        density=density,
        unnormalised=unnormalised,
        droplet_radius_A=float(droplet_radius_A),
        temperature_K=float(temperature_K),
    )
