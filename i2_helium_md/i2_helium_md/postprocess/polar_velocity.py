"""Polar (|v|, phi) histogram and cos^2 anisotropy fit.

Reproduces the simulation-side polar VMI panels of
``legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v3.m``
(lines 96-111: angular distribution + nlinfit cos^2 model + per-velocity beta(v)).

The legacy MATLAB code computed an anisotropy fit on top of an Abel-inverted
2-D experimental VMI image. We have the *full 3-D* simulated velocities, so
no Abel inversion is needed -- we bin the lab-frame final velocities directly
into a (|v|, phi) histogram and fit::

    f(phi) = a + b * cos(phi - phi0)^2

(equivalent to a Legendre-P2 expansion). The conventional anisotropy
parameter is recovered as ``beta = 2*b / (2*a + b)`` with the convention
used in photodissociation (beta in [-1, 2]; +2 = pure cos^2, -1 = pure sin^2).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

from ..physics.constants import U as U_KG
from ..simulation.checkpoint import IonCheckpoint


@dataclass(frozen=True)
class PolarHistogram:
    """2-D histogram of final ion velocities in polar (|v|, phi) coordinates.

    Attributes
    ----------
    counts : np.ndarray, shape (n_v_bins, n_phi_bins)
        Atom count in each (v, phi) bin.
    v_centers_Aps : np.ndarray, shape (n_v_bins,)
        Bin centers along the speed axis in angstrom/ps.
    v_edges_Aps : np.ndarray, shape (n_v_bins + 1,)
    phi_centers_rad : np.ndarray, shape (n_phi_bins,)
        Bin centers along the phi axis in radians, in [0, 2 pi).
    phi_edges_rad : np.ndarray, shape (n_phi_bins + 1,)
    mass_amu : float
    num_atoms_used : int
    """

    counts: np.ndarray
    v_centers_Aps: np.ndarray
    v_edges_Aps: np.ndarray
    phi_centers_rad: np.ndarray
    phi_edges_rad: np.ndarray
    mass_amu: float
    num_atoms_used: int


@dataclass(frozen=True)
class AnisotropyFit:
    """Result of a single cos^2 fit to a 1-D phi distribution.

    Model: ``f(phi) = a + b * cos(phi - phi0)**2``.

    ``beta`` is the conventional photodissociation anisotropy parameter,
    derived as ``beta = 2 * b / (2 * a + b)`` (range [-1, 2]). Returns
    ``np.nan`` parameters when the fit fails to converge or the input
    is empty.
    """

    a: float
    b: float
    phi0_rad: float
    beta: float
    residual: float
    success: bool


@dataclass(frozen=True)
class BetaCurve:
    """Anisotropy beta(v) computed by fitting cos^2 in each |v| bin."""

    v_centers_Aps: np.ndarray
    beta: np.ndarray
    beta_uncertainty: np.ndarray
    valid: np.ndarray  # bool mask: True where the fit converged


def polar_velocity_histogram(
    ion: IonCheckpoint,
    *,
    n_v_bins: int = 80,
    n_phi_bins: int = 72,
    v_max_Aps: float = 28.0,
    mass_amu: float | None = None,
    mass_tolerance_amu: float = 0.5,
    require_outside: bool = True,
) -> PolarHistogram:
    """Bin final ion velocities into (|v|, phi) polar coordinates.

    ``phi = arctan2(vy_final, vx_final) + pi`` so the result is
    consistent with :func:`postprocess.energy_balance.phi_histogram`.

    The 1-D radial histogram is recovered by ``counts.sum(axis=1)``;
    the 1-D phi histogram by ``counts.sum(axis=0)``.
    """
    if n_v_bins < 1:
        raise ValueError(f"n_v_bins must be >= 1, got {n_v_bins}")
    if n_phi_bins < 1:
        raise ValueError(f"n_phi_bins must be >= 1, got {n_phi_bins}")
    if v_max_Aps <= 0.0:
        raise ValueError(f"v_max_Aps must be > 0, got {v_max_Aps}")

    masses_amu = np.round(np.asarray(ion.mass_final_kg) / U_KG)
    if mass_amu is None:
        mass_mask = np.ones(masses_amu.shape, dtype=bool)
        mass_used = float("nan")
    else:
        mass_mask = np.abs(masses_amu - mass_amu) <= mass_tolerance_amu
        mass_used = float(mass_amu)

    if require_outside:
        outside = np.concatenate(
            [ion.b_ion_outside, ion.b_ion_outside]
        ).astype(bool)
        sel = mass_mask & outside
    else:
        sel = mass_mask

    vx = np.asarray(ion.velocities_final_x)[sel]
    vy = np.asarray(ion.velocities_final_y)[sel]
    vz = np.asarray(ion.velocities_final_z)[sel]
    speed = np.sqrt(vx * vx + vy * vy + vz * vz)
    phi = np.mod(np.arctan2(vy, vx) + np.pi, 2.0 * np.pi)

    v_edges = np.linspace(0.0, v_max_Aps, n_v_bins + 1)
    phi_edges = np.linspace(0.0, 2.0 * np.pi, n_phi_bins + 1)
    counts, _, _ = np.histogram2d(speed, phi, bins=(v_edges, phi_edges))

    return PolarHistogram(
        counts=counts.astype(float),
        v_centers_Aps=0.5 * (v_edges[:-1] + v_edges[1:]),
        v_edges_Aps=v_edges,
        phi_centers_rad=0.5 * (phi_edges[:-1] + phi_edges[1:]),
        phi_edges_rad=phi_edges,
        mass_amu=mass_used,
        num_atoms_used=int(sel.sum()),
    )


def _cos2_model(phi: np.ndarray, a: float, b: float, phi0: float) -> np.ndarray:
    return a + b * np.cos(phi - phi0) ** 2


def _fit_cos2(phi: np.ndarray, counts: np.ndarray) -> AnisotropyFit:
    """Single-shot cos^2 fit. Returns NaN params on failure or empty input."""
    if counts.size == 0 or not np.any(counts > 0):
        return AnisotropyFit(
            a=float("nan"), b=float("nan"), phi0_rad=float("nan"),
            beta=float("nan"), residual=float("nan"), success=False,
        )
    a0 = float(counts.min())
    b0 = float(counts.max() - counts.min())
    phi0_init = float(phi[int(np.argmax(counts))])
    try:
        popt, _pcov = curve_fit(
            _cos2_model, phi, counts,
            p0=(a0, b0, phi0_init),
            maxfev=5000,
        )
    except (RuntimeError, ValueError):
        return AnisotropyFit(
            a=float("nan"), b=float("nan"), phi0_rad=float("nan"),
            beta=float("nan"), residual=float("nan"), success=False,
        )
    a_fit, b_fit, phi0_fit = popt
    fit_curve = _cos2_model(phi, *popt)
    residual = float(np.sqrt(np.mean((fit_curve - counts) ** 2)))
    denom = 2.0 * a_fit + b_fit
    beta = float(2.0 * b_fit / denom) if denom != 0.0 else float("nan")
    return AnisotropyFit(
        a=float(a_fit), b=float(b_fit), phi0_rad=float(phi0_fit),
        beta=beta, residual=residual, success=True,
    )


def anisotropy_fit(polar: PolarHistogram) -> AnisotropyFit:
    """Fit a single cos^2 model to the velocity-summed phi distribution."""
    return _fit_cos2(polar.phi_centers_rad, polar.counts.sum(axis=0))


def beta_of_velocity(
    polar: PolarHistogram,
    *,
    min_counts_per_v_bin: int = 50,
) -> BetaCurve:
    """Per-|v|-bin cos^2 fit; bins with too few counts are flagged invalid."""
    n_v = polar.v_centers_Aps.size
    beta = np.full(n_v, np.nan, dtype=float)
    sigma = np.full(n_v, np.nan, dtype=float)
    valid = np.zeros(n_v, dtype=bool)
    for i in range(n_v):
        row = polar.counts[i]
        if row.sum() < min_counts_per_v_bin:
            continue
        fit = _fit_cos2(polar.phi_centers_rad, row)
        if fit.success:
            beta[i] = fit.beta
            sigma[i] = fit.residual
            valid[i] = True
    return BetaCurve(
        v_centers_Aps=polar.v_centers_Aps,
        beta=beta,
        beta_uncertainty=sigma,
        valid=valid,
    )
