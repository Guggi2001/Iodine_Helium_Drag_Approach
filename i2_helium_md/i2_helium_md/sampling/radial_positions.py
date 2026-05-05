"""Sampling molecule positions inside helium droplets.

Ports ``generate_radial_samples_3d.m``.

Each molecule's distance from its droplet center is drawn from a thermal
Boltzmann distribution in the droplet solvation potential. The full
probability density (in 3D, after integrating out angles) is

.. math::

    p(r) \\propto r^2 \\cdot \\exp\\left(-\\frac{U_{\\text{drop}}(r - R)}{k_B T}\\right)

where the ``r^2`` factor is the spherical volume element (``4 pi r^2 dr``
with ``4 pi`` absorbed into the normalization). The ``r^2`` Jacobian
already accounts for the ``sin(theta) d\\theta d\\phi`` from the angular
integration.

We sample by rejection sampling: propose ``r ~ Uniform(0, 2R)``,
accept with probability ``p(r) / p_max``.

Conventions
-----------
* All inputs in Angstrom, temperatures in Kelvin.
* Energies internally in eV (then divided by ``k_B * T`` to be unitless).
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from ..config import SimConfig
from ..physics.constants import EV, K_B
from ..physics.potentials import droplet_potential


# ===========================================================================
# Public API
# ===========================================================================
def sample_radial_positions(
    cfg: SimConfig,
    droplet_radii: np.ndarray,
    *,
    E_max_meV: float = 200.0,
    r_step_angstrom: float = 0.01,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Sample one radial position per molecule, Boltzmann-weighted in the droplet.

    For each unique droplet radius in the ensemble, builds the radial
    probability distribution
    ``p(r) ~ r^2 * exp(-U_drop(r - R) / (k_B T))`` and rejection-samples
    the requested number of points from it.

    Parameters
    ----------
    cfg : SimConfig
        Simulation config. Uses ``T_particles_K``, ``potential_steepness_molecule``,
        ``binding_energy_molecule_meV``, and ``seed``.
    droplet_radii : np.ndarray, shape (N,)
        One droplet radius per molecule, in Angstrom.
    E_max_meV : float, optional
        Maximum energy considered for the Boltzmann normalization. Default
        200 meV (matches the MATLAB default; in practice the integral
        converges much earlier).
    r_step_angstrom : float, optional
        Grid spacing for the radial probability density evaluation, used
        only for finding the proposal envelope ``y_max``. Default 0.01 A.
    rng : np.random.Generator, optional
        Provide a numpy RNG for reproducibility. If None, build one from
        ``cfg.seed``.

    Returns
    -------
    r : np.ndarray, shape (N,)
        Radial distance from droplet center for each molecule, in Angstrom.

    Notes
    -----
    The MATLAB version had an inner loop with batches of 1000 proposals.
    We use a single oversampled vectorized batch sized from an estimated
    acceptance rate. This is cleaner and equally fast in practice.

    The ``E`` array used in MATLAB to "normalize p(E)" is unnecessary --
    the radial density gets re-normalized via the ``r``-axis integral
    anyway. We drop it.
    """
    if rng is None:
        rng = np.random.default_rng(cfg.seed)

    droplet_radii = np.asarray(droplet_radii, dtype=float).ravel()

    # Convert binding energy to eV (the unit of `droplet_potential`'s output).
    # Note: cfg.binding_energy_molecule_meV is meV; we want eV for the potential.
    binding_eV = cfg.binding_energy_molecule_meV / 1000.0
    steepness = cfg.potential_steepness_molecule
    T = cfg.T_particles_K

    # Each unique droplet radius gets its own p(r) and rejection batch.
    unique_radii, inverse = np.unique(droplet_radii, return_inverse=True)

    # We'll fill positions in the order droplet_radii appears.
    r_out = np.empty(droplet_radii.shape[0], dtype=float)

    for i, R in enumerate(unique_radii):
        n_needed = int((inverse == i).sum())
        if n_needed == 0:
            continue

        r_samples = _sample_radial_for_one_droplet(
            droplet_radius=float(R),
            n_samples=n_needed,
            T_K=T,
            steepness=steepness,
            binding_eV=binding_eV,
            r_step=r_step_angstrom,
            rng=rng,
        )
        r_out[inverse == i] = r_samples

    return r_out


# ===========================================================================
# Internals
# ===========================================================================
def _radial_probability_density(
    r: np.ndarray,
    droplet_radius: float,
    T_K: float,
    steepness: float,
    binding_eV: float,
) -> np.ndarray:
    """Unnormalized radial Boltzmann density.

    Returns ``r^2 * exp(-U/(k_B T))`` where ``U = droplet_potential(r - R)``.

    The ``r^2`` is the spherical volume element. The exponential is the
    Boltzmann weight in eV / (eV/K * K).
    """
    U_eV = droplet_potential(
        r - droplet_radius,
        steepness=steepness,
        binding_energy=binding_eV,
    )
    # k_B in eV/K
    k_B_eV_per_K = K_B / EV
    return r ** 2 * np.exp(-U_eV / (k_B_eV_per_K * T_K))


def _sample_radial_for_one_droplet(
    *,
    droplet_radius: float,
    n_samples: int,
    T_K: float,
    steepness: float,
    binding_eV: float,
    r_step: float,
    rng: np.random.Generator,
    safety_factor: float = 1.05,
) -> np.ndarray:
    """Rejection-sample n_samples radii inside one droplet.

    Strategy
    --------
    1. Build a fine grid r in [0, 2R].
    2. Evaluate p(r) on the grid; record p_max.
    3. Estimate acceptance rate via mean(p)/p_max.
    4. Propose ``ceil(n_samples / accept_rate * safety_factor)`` samples in
       one vectorized batch from Uniform(0, 2R).
    5. Accept where ``Uniform(0, p_max) < p(r_proposal)``.
    6. If we somehow under-shot (rare with safety_factor > 1), top up with
       another batch.
    """
    r_max = 2.0 * droplet_radius

    # 1. Grid for envelope estimation.
    r_grid = np.arange(0.0, r_max + r_step, r_step)
    p_grid = _radial_probability_density(
        r_grid, droplet_radius, T_K, steepness, binding_eV,
    )
    y_max = p_grid.max()
    if y_max <= 0.0:
        raise RuntimeError(
            f"radial probability density is identically zero "
            f"for droplet radius {droplet_radius} A; "
            "check temperature and binding energy"
        )

    # 2. Estimate acceptance rate from grid mean.
    accept_rate = max(p_grid.mean() / y_max, 1e-3)

    accepted: list[np.ndarray] = []
    n_have = 0
    while n_have < n_samples:
        n_propose = int(np.ceil(
            (n_samples - n_have) / accept_rate * safety_factor
        ))
        r_proposal = rng.uniform(0.0, r_max, size=n_propose)
        u_trial = rng.uniform(0.0, y_max, size=n_propose)
        p_eval = _radial_probability_density(
            r_proposal, droplet_radius, T_K, steepness, binding_eV,
        )
        accept_mask = u_trial < p_eval
        new_accepted = r_proposal[accept_mask]
        accepted.append(new_accepted)
        n_have += new_accepted.size

    return np.concatenate(accepted)[:n_samples]
