"""Sampling of molecular orientations and bond lengths.

This module produces the random angular degrees of freedom for each
molecule at the start of a simulation:

* **Position angles** ``(beta, gamma)`` -- the spherical angles for the
  molecule's centre-of-mass position inside the droplet. Always
  uniform on a sphere.
* **Orientation angles** ``(alpha, delta)`` -- the spherical angles for
  the molecule's interatomic axis. Distribution depends on the laser
  geometry:

  - **Anisotropic** (single-pulse mode): linear-polarisation cos²φ
    weighting. The probability that a molecule is excited is
    proportional to ``cos²(angle to polarisation axis)``, which we
    realise by rejection sampling. From: *Molecular reorientation
    during dissociative multiphoton ionisation*, PRA 1993.
  - **Isotropic**: uniform on a sphere.

* **Bond length** -- ``R0_GS + Gaussian(0, deltaR0)``, the equilibrium
  bond length plus zero-point fluctuation width. ``deltaR0=0`` reduces
  to a fixed bond length.

The samplers return only angles and lengths; the conversion to atomic
xyz coordinates is the responsibility of the simulation driver
(``simulation/neutral.py``) because it depends on the 2N array layout
convention used throughout the codebase.

Replaces the orientation-sampling block of
``vmi_sim_3d_neutral_propa_HeDFT_mimic.m`` (lines ~225-275).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ===========================================================================
# Result container
# ===========================================================================
@dataclass(frozen=True)
class MolecularOrientations:
    """Per-molecule angular degrees of freedom.

    All arrays have shape ``(num_molecules,)``. Angles are in radians.

    Attributes
    ----------
    beta : np.ndarray
        Azimuthal angle of the molecule's centre position, in [0, 2π).
    gamma : np.ndarray
        Polar angle of the molecule's centre position, in [0, π].
    alpha : np.ndarray
        Azimuthal angle of the molecule's interatomic axis, in [0, 2π).
    delta : np.ndarray
        Polar angle of the molecule's interatomic axis, in [0, π].
    bond_length_angstrom : np.ndarray
        Equilibrium bond length plus zero-point fluctuation, in Å.
    """

    beta: np.ndarray
    gamma: np.ndarray
    alpha: np.ndarray
    delta: np.ndarray
    bond_length_angstrom: np.ndarray


# ===========================================================================
# Public API
# ===========================================================================
def sample_orientations(
    num_molecules: int,
    *,
    R0_GS_angstrom: float,
    deltaR0_angstrom: float,
    anisotropic: bool,
    rng: np.random.Generator | None = None,
) -> MolecularOrientations:
    """Sample molecular orientations and bond lengths for the ensemble.

    Parameters
    ----------
    num_molecules : int
        Number of molecules to sample.
    R0_GS_angstrom : float
        Ground-state equilibrium bond length in Å. Bond lengths are
        drawn as ``R0_GS + N(0, deltaR0)``.
    deltaR0_angstrom : float
        Standard deviation of the bond length fluctuation, in Å.
        Use ``0.0`` for a fixed bond length.
    anisotropic : bool
        If True, the molecular axis orientation is sampled with the
        cos²φ weighting appropriate for dissociation by a linearly
        polarised pump pulse. If False, uniform on a sphere.
    rng : np.random.Generator, optional
        Reproducible RNG. If None, a fresh default RNG is constructed.

    Returns
    -------
    MolecularOrientations

    Raises
    ------
    ValueError
        If ``num_molecules <= 0`` or ``deltaR0_angstrom < 0``.
    """
    if num_molecules <= 0:
        raise ValueError(f"num_molecules must be positive, got {num_molecules}")
    if deltaR0_angstrom < 0:
        raise ValueError(
            f"deltaR0_angstrom must be non-negative, got {deltaR0_angstrom}"
        )

    if rng is None:
        rng = np.random.default_rng()

    # 1. Centre-of-mass position angles -- always uniform on the sphere.
    beta, gamma = _sample_uniform_sphere_angles(num_molecules, rng=rng)

    # 2. Interatomic axis orientation -- anisotropic or isotropic.
    if anisotropic:
        alpha, delta = _sample_anisotropic_axis_angles(num_molecules, rng=rng)
    else:
        alpha, delta = _sample_uniform_sphere_angles(num_molecules, rng=rng)

    # 3. Bond length with optional zero-point fluctuation.
    bond_length = (
        R0_GS_angstrom + deltaR0_angstrom * rng.standard_normal(num_molecules)
    )

    return MolecularOrientations(
        beta=beta,
        gamma=gamma,
        alpha=alpha,
        delta=delta,
        bond_length_angstrom=bond_length,
    )


# ===========================================================================
# Internal samplers
# ===========================================================================
def _sample_uniform_sphere_angles(
    n: int,
    *,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample (azimuth, polar) angle pairs uniformly on the unit sphere.

    Standard inverse-CDF sampling:

    * azimuth ~ Uniform(0, 2π)
    * cos(polar) ~ Uniform(-1, 1)  ->  polar = arccos(cos(polar))

    The ``cos(polar)`` step is what makes the distribution uniform on
    the sphere rather than uniform in (azimuth, polar), which would
    produce too many samples near the poles. This is the same fix
    that the MATLAB code applies.
    """
    azimuth = rng.uniform(0.0, 2.0 * np.pi, size=n)
    cos_polar = rng.uniform(-1.0, 1.0, size=n)
    polar = np.arccos(cos_polar)
    return azimuth, polar


def _sample_anisotropic_axis_angles(
    n: int,
    *,
    rng: np.random.Generator,
    max_oversample_factor: int = 20,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample axis angles ``(alpha, delta)`` weighted by cos²φ.

    Here φ is the angle between the molecular axis and the laser
    polarisation axis (the lab x-axis by convention). The probability
    of accepting a candidate ``(alpha, delta)`` proposed uniformly on
    the sphere is

        p = |cos(alpha) * sin(delta)|^2

    where ``cos(alpha) * sin(delta)`` is the x-component of the unit
    vector along the molecular axis.

    Implemented by rejection sampling with batched proposals to keep
    the number of RNG calls bounded.

    Parameters
    ----------
    n : int
        Number of accepted samples to return.
    rng : np.random.Generator
    max_oversample_factor : int, optional
        Safety guard. The mean acceptance rate of cos² rejection is
        1/3, so we aim for ``3*n`` proposals per round. If after
        ``max_oversample_factor`` rounds we still don't have enough
        accepted samples, raise. Default 20.

    Returns
    -------
    alpha, delta : np.ndarray
        Each of shape ``(n,)``, in radians.
    """
    accepted_alpha: list[np.ndarray] = []
    accepted_delta: list[np.ndarray] = []
    accepted_count = 0

    # Aim for slightly more than 3*n per round (acceptance rate is 1/3).
    batch_size = max(int(np.ceil(3.5 * n)), 32)

    for _ in range(max_oversample_factor):
        if accepted_count >= n:
            break
        a, d = _sample_uniform_sphere_angles(batch_size, rng=rng)
        cos_phi = np.cos(a) * np.sin(d)         # = x-component of axis unit vec
        p_accept = cos_phi ** 2                  # in [0, 1]
        u = rng.uniform(0.0, 1.0, size=batch_size)
        keep = u < p_accept
        accepted_alpha.append(a[keep])
        accepted_delta.append(d[keep])
        accepted_count += int(keep.sum())
    else:
        # else clause runs if loop exits without break
        raise RuntimeError(
            f"anisotropic angle sampling failed to produce {n} samples in "
            f"{max_oversample_factor} batches of {batch_size} -- "
            "this should not happen for cos² with ~33% acceptance; "
            "RNG or implementation issue?"
        )

    alpha = np.concatenate(accepted_alpha)[:n]
    delta = np.concatenate(accepted_delta)[:n]
    return alpha, delta
