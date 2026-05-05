"""Sampling helium-droplet sizes from log-normal source distributions.

Ports two MATLAB functions:

* ``get_dropletsize.m``         -> :func:`mean_droplet_size`
* ``generate_droplet_sizes.m``  -> :func:`sample_droplet_sizes`

The sampler can run in two modes, controlled by ``mode``:

* ``"raw"``         -- return raw log-normal samples (the simpler legacy file).
* ``"post_pickup"`` -- simulate I2 pickup + evaporation and return the
                       post-pickup distribution (default; matches the full
                       legacy `generate_droplet_sizes`).

Physics references
------------------
* Lackner dissertation (TU Graz)            -- nozzle correlation, Kornilov delta
* Kornilov (Mol. Phys. 2009)                -- log-normal width parameter
* Yang et al., 10.1007/978-3-030-94896-2_1  -- He-droplet chemical potential

Default source parameters reproduce the legacy `generate_droplet_sizes_simpler`
(p=40 mbar, T=23 K, d=5 um). For other operating conditions, override via
``cfg.p_source_mbar`` / ``cfg.T_source_K`` or pass values to the function.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from ..config import SimConfig


# ---------------------------------------------------------------------------
# Empirical constants from get_dropletsize.m
# ---------------------------------------------------------------------------
# <N> = k1 * p^k2 * T^k3 * d^k4
_K1: float = 4e5         # prefactor
_K2: float = 0.97        # pressure exponent
_K3: float = -3.88       # temperature exponent
_K4: float = 2.0         # nozzle-diameter exponent

_DEFAULT_NOZZLE_DIAMETER_UM: float = 5.0    # micrometers
_LOGNORMAL_DELTA: float = 0.625             # Kornilov 2009; constant for log-normal width

# Density used to convert N <-> R
_RHO_BULK_HE: float = 0.0218                # atoms / Angstrom^3
_RHO_DROPLET: float = 0.8 * _RHO_BULK_HE


# ---------------------------------------------------------------------------
# Pickup-simulation constants (from generate_droplet_sizes.m)
# ---------------------------------------------------------------------------
_HE_DENSITY_PER_M3: float = 2.18e28              # n_he, particles / m^3
_DROPLET_DENSITY_PER_M3: float = 0.8 * _HE_DENSITY_PER_M3

_MEAN_FREE_PATH_M: float = 4e-10                 # l, meters
_REL_ENERGY_LOSS_PER_COLLISION: float = 0.04     # epsilon

# Dopant-specific solvation energy and initial kinetic energy.
# IMPORTANT: There is an inconsistency between the thesis text and figure:
# - Thesis text states E_solv = 30 meV for iodine.
# - Thesis figure 3.2 was clearly produced with E_solv = 14 meV (the
#   reduced-σ distributions shift right by several thousand atoms, which
#   only happens with the stronger correction E_solv=14 gives).
# - The legacy MATLAB code uses E_solv = 14 meV.
#
# We default to 14 because that matches both the figure and the legacy
# MATLAB code, which is the more useful regression target. Pass
# `E_solv_meV=30.0` explicitly if you want the thesis-text value.
_DEFAULT_E_SOLV_meV: float = 14.0       # matches MATLAB and thesis figure 3.2
_E_KIN_THERMAL_meV: float = 25.85       # 3/2 k_B * 300 K, dopant thermal kinetic energy

_PICKUP_REGION_LENGTH_M: float = 2e-2            # a, meters
_PICKUP_GAS_TEMP_K: float = 275.0
_PICKUP_GAS_PRESSURE_MBAR: float = 0.5e-5        # pickup-cell pressure
_K_B: float = 1.381e-23                          # J/K (legacy local copy; matches constants.K_B)


# ===========================================================================
# Public functions
# ===========================================================================
def mean_droplet_size(
    p_source_mbar: float,
    T_source_K: float,
    nozzle_diameter_um: float = _DEFAULT_NOZZLE_DIAMETER_UM,
) -> tuple[float, float]:
    """Mean droplet size and radius for a given He source.

    Direct port of ``get_dropletsize.m``. Empirical correlation:

    .. math::

        \\langle N \\rangle = k_1 \\cdot p^{k_2} \\cdot T^{k_3} \\cdot d^{k_4}

    with ``k1=4e5, k2=0.97, k3=-3.88, k4=2``. The radius assumes a spherical
    droplet at 0.8 of bulk He density.

    Parameters
    ----------
    p_source_mbar : float
        Source pressure in mbar.
    T_source_K : float
        Source temperature in Kelvin.
    nozzle_diameter_um : float, optional
        Nozzle diameter in micrometers. Default 5 um.

    Returns
    -------
    N_mean : float
        Mean number of He atoms per droplet.
    R_mean : float
        Mean droplet radius in Angstrom.
    """
    N_mean = _K1 * (p_source_mbar ** _K2) * (T_source_K ** _K3) * (nozzle_diameter_um ** _K4)
    R_mean = (3.0 / (4.0 * np.pi * _RHO_DROPLET)) ** (1.0 / 3.0) * N_mean ** (1.0 / 3.0)
    return N_mean, R_mean


def droplet_radius_from_N(N: np.ndarray | float) -> np.ndarray | float:
    """Convert droplet size N to droplet radius in Angstrom.

    Uses the same density assumption as :func:`mean_droplet_size`
    (0.8 of bulk He density, ~0.0175 atoms/A^3).
    """
    return (3.0 / (4.0 * np.pi * _RHO_DROPLET)) ** (1.0 / 3.0) * np.asarray(N) ** (1.0 / 3.0)


def sample_droplet_sizes(
    cfg: SimConfig,
    *,
    mode: Literal["raw", "post_pickup"] = "post_pickup",
    reduced_crosssection: bool = False,
    nozzle_diameter_um: float = _DEFAULT_NOZZLE_DIAMETER_UM,
    E_solv_meV: float = _DEFAULT_E_SOLV_meV,
    rng: np.random.Generator | None = None,
    oversample_factor: int = 10,
    max_retries: int = 3,
) -> np.ndarray:
    """Sample helium-droplet sizes for ``cfg.num_molecules`` simulated molecules.

    Two modes:

    * ``mode="raw"``:
        Return ``cfg.num_molecules`` raw samples from a log-normal
        distribution with mean given by :func:`mean_droplet_size` and width
        delta=0.625 (Kornilov 2009).

    * ``mode="post_pickup"``:
        Generate ``oversample_factor * num_molecules`` raw samples, simulate
        I2 pickup with evaporation cooling, and return only droplets that
        picked up exactly one I2 molecule. This matches the legacy
        ``generate_droplet_sizes.m`` behaviour and is what experiments
        actually observe.

        For cold sources (large droplets), most droplets pick up multiple
        molecules and the single-pickup yield is low. The function
        automatically retries with a larger oversample factor up to
        ``max_retries`` times before raising.

    Parameters
    ----------
    cfg : SimConfig
        Used for ``num_molecules``, ``p_source_mbar``, ``T_source_K``, and
        ``seed``.
    mode : {"raw", "post_pickup"}, optional
        Which distribution to return. Default "post_pickup".
    reduced_crosssection : bool, optional
        If True (only relevant in "post_pickup" mode), use the
        kinetic-energy-dissipation-reduced pickup cross section instead of
        the geometric one. Default False, matching the legacy preset.
    nozzle_diameter_um : float, optional
        Override the default nozzle diameter (5 um).
    E_solv_meV : float, optional
        Dopant solvation energy in meV (used for the reduced cross
        section threshold). Default 14 meV (matches legacy MATLAB and
        thesis figure 3.2). The thesis *text* says 30 meV for iodine but
        that value does not reproduce the thesis *figure* -- the inconsistency
        is documented in `migration_log.md`. Pass `E_solv_meV=30.0`
        explicitly if you want to follow the thesis text instead.
    rng : np.random.Generator, optional
        Provide a numpy RNG for reproducibility. If None, build one from
        ``cfg.seed``.
    oversample_factor : int, optional
        How many raw samples to generate per requested molecule before the
        pickup simulation prunes them. Default 10. For cold sources where
        multi-pickup dominates, the function will automatically grow this
        up to ``max_retries`` times.
    max_retries : int, optional
        Maximum number of times to grow oversample_factor (each by 4x)
        when single-pickup yield is too low. Default 3 (-> up to 640x).

    Returns
    -------
    N : np.ndarray, shape (num_molecules,)
        Sampled droplet sizes (numbers of He atoms).

    Raises
    ------
    ValueError
        If after `max_retries` we still can't produce enough single-pickup
        droplets. Use ``mode="raw"`` or accept multi-pickup samples in that
        case.
    """
    if rng is None:
        rng = np.random.default_rng(cfg.seed)

    N_mean, _ = mean_droplet_size(
        cfg.p_source_mbar, cfg.T_source_K, nozzle_diameter_um,
    )

    # Log-normal parameter mu derived from N_mean and delta
    # (matches MATLAB: mu = log(N_mean) - delta^2/2)
    mu = np.log(N_mean) - _LOGNORMAL_DELTA ** 2 / 2.0
    delta = _LOGNORMAL_DELTA

    if mode == "raw":
        return rng.lognormal(mean=mu, sigma=delta, size=cfg.num_molecules)

    if mode == "post_pickup":
        factor = oversample_factor
        for attempt in range(max_retries + 1):
            oversample = max(cfg.num_molecules * factor, 10_000)
            N_raw = rng.lognormal(mean=mu, sigma=delta, size=oversample)
            try:
                return _simulate_pickup(
                    N_raw,
                    num_target=cfg.num_molecules,
                    reduced_crosssection=reduced_crosssection,
                    E_solv_meV=E_solv_meV,
                    rng=rng,
                )
            except ValueError:
                if attempt == max_retries:
                    raise
                factor *= 4   # try harder on the next round

    raise ValueError(f"unknown mode {mode!r}; expected 'raw' or 'post_pickup'")


# ===========================================================================
# Analytical Poisson convolution (Treiber's didactic formula)
# ===========================================================================
def conditional_size_distributions_analytical(
    N_grid: np.ndarray,
    *,
    p_source_mbar: float,
    T_source_K: float,
    nozzle_diameter_um: float = _DEFAULT_NOZZLE_DIAMETER_UM,
    k: int = 1,
    # Treiber-script parameter values:
    p_pickup_gas_mbar: float = 8.70e-5,
    T_pickup_gas_K: float = 293.0,
    E_solv_meV: float = 30.0,
    E_kin_thermal_meV: float = 38.78,
    n_he_per_m3: float = 2.18e28,        # Treiber uses bulk He density, NOT 0.8 x n_he
    mean_free_path_m: float = 4.0e-10,
    rel_energy_loss: float = 0.04,
    pickup_region_length_m: float = 2.0e-2,
) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form conditional droplet-size distributions from Treiber's script.

    Direct port of ``conditional_droplet_size_distribution_simplified.m`` --
    the didactic script that produced thesis figure 3.2. Returns BOTH the
    normal-σ and reduced-σ conditional distributions on ``N_grid``,
    using Treiber's normalization convention.

    Critical implementation details (verified by line-by-line comparison
    with the MATLAB script):

    * The droplet radius uses **bulk** helium density ``n_he = 2.18e28``,
      NOT the 0.8x droplet density used elsewhere in this codebase. This
      single 1.077x scaling on R has a 1.16x effect on σ and dramatically
      changes the threshold cutoff in the reduced-σ branch.
    * The reduced-σ conditional uses the **normal-σ** normalization
      (``p_k_normalization``) in its denominator, NOT a separate
      reduced-σ normalization. The reduced-σ distribution therefore does
      not integrate to 1; its integral is the fraction of one-pickup
      events that came from droplets above the kinetic-energy threshold.
      This is what the thesis figure shows -- dashed lines have **smaller**
      peaks than solid lines, not equal-area peaks.

    Use this function only for reproducing the thesis figure or for
    analytical comparison. It is a simpler physical model than our
    Monte-Carlo sampler -- no evaporation, no destroyed droplets.

    Parameters match the MATLAB script's variable names; defaults are
    the values used to produce the thesis figure exactly.

    Returns
    -------
    p_normal : np.ndarray
        Conditional p(N | k pickups) using the geometric cross section.
        Integrates to 1.
    p_reduced : np.ndarray
        Conditional p(N | k pickups) using the kinetic-energy-reduced
        cross section, normalised by the SAME factor as ``p_normal``
        (so its integral is < 1, as in the thesis figure).
    """
    from math import factorial

    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")

    N = np.asarray(N_grid, dtype=float)

    # ---- log-normal of source droplets ----
    N_mean, _ = mean_droplet_size(p_source_mbar, T_source_K, nozzle_diameter_um)
    mu = np.log(N_mean) - _LOGNORMAL_DELTA ** 2 / 2.0
    lognpdf = (
        1.0 / (N * _LOGNORMAL_DELTA * np.sqrt(2.0 * np.pi))
        * np.exp(-(np.log(N) - mu) ** 2 / (2.0 * _LOGNORMAL_DELTA ** 2))
    )

    # ---- droplet radius (Treiber: uses BULK He density) ----
    R = (3.0 * N / (4.0 * np.pi * n_he_per_m3)) ** (1.0 / 3.0)

    # ---- pickup-cell gas density ----
    # MATLAB: n_gas = p_gas[mbar] * 100[Pa/mbar] / (T_gas * kb)
    n_gas = p_pickup_gas_mbar * 100.0 / (T_pickup_gas_K * _K_B)
    a = pickup_region_length_m

    # ---- cross sections ----
    sigma_normal_m2 = np.pi * R ** 2

    E_kin_0 = E_kin_thermal_meV + E_solv_meV
    log_ratio_sq = np.log(E_solv_meV / E_kin_0) ** 2
    b_thresh_sq = R ** 2 - mean_free_path_m ** 2 / (4.0 * rel_energy_loss ** 2) * log_ratio_sq
    sigma_reduced_m2 = np.pi * np.maximum(b_thresh_sq, 0.0)

    # ---- Poisson PMF P(k | N) = lambda^k exp(-lambda) / k! ----
    lam_normal = a * n_gas * sigma_normal_m2
    lam_reduced = a * n_gas * sigma_reduced_m2

    # ---- normalised P(k | N) along N axis (lines 71, 75 of MATLAB) ----
    Z_normal = np.trapezoid(
        lam_normal ** k / factorial(k) * np.exp(-lam_normal), N
    )
    Z_reduced = np.trapezoid(
        lam_reduced ** k / factorial(k) * np.exp(-lam_reduced), N
    )
    p_k_normal = lam_normal ** k / factorial(k) * np.exp(-lam_normal) / Z_normal
    p_k_reduced = lam_reduced ** k / factorial(k) * np.exp(-lam_reduced) / Z_reduced

    # ---- conditional distributions, BOTH normalised by p_k_normalization ----
    # MATLAB lines 89-95 (the key normalisation choice):
    #   p_k_normalization = trapz(N, p_k(N) .* lognpdf(...))
    #   p_conditional = p_k_reduced(N) .* lognpdf(...) / p_k_normalization
    #   p_normal      = p_k(N)         .* lognpdf(...) / p_k_normalization
    p_k_normalization = np.trapezoid(p_k_normal * lognpdf, N)
    p_normal = p_k_normal * lognpdf / p_k_normalization
    p_reduced = p_k_reduced * lognpdf / p_k_normalization

    return p_normal, p_reduced


# ===========================================================================
# Pickup simulation (post_pickup mode internals)
# ===========================================================================
def _droplet_radius_m(N: np.ndarray) -> np.ndarray:
    """Droplet radius in meters (used inside the pickup-cell simulation)."""
    return (3.0 * N / (4.0 * np.pi * _DROPLET_DENSITY_PER_M3)) ** (1.0 / 3.0)


def _impact_parameter_threshold_m(
    N: np.ndarray,
    E_solv_meV: float = _DEFAULT_E_SOLV_meV,
) -> np.ndarray:
    """Threshold impact parameter beyond which kinetic energy can't dissipate
    enough for the dopant to be retained.

    Port of MATLAB ``b_thresh``.  Uses ``np.maximum(..., 0)`` only at the
    end (after the sqrt argument); for arguments that go negative this
    returns 0 (small droplets cannot capture the dopant).

    Parameters
    ----------
    N : np.ndarray
        Droplet sizes.
    E_solv_meV : float, optional
        Dopant solvation energy in meV. Default 30 (thesis value for I2).
    """
    R = _droplet_radius_m(N)
    E_kin_0 = _E_KIN_THERMAL_meV + E_solv_meV
    inside = (
        R ** 2
        - _MEAN_FREE_PATH_M ** 2
        / (4.0 * _REL_ENERGY_LOSS_PER_COLLISION ** 2)
        * np.log(E_solv_meV / E_kin_0) ** 2
    )
    return np.sqrt(np.maximum(inside, 0.0))


def _pickup_cross_section_m2(
    N: np.ndarray,
    reduced: bool,
    E_solv_meV: float = _DEFAULT_E_SOLV_meV,
) -> np.ndarray:
    """Geometric or reduced pickup cross section, in m^2."""
    if reduced:
        return np.pi * _impact_parameter_threshold_m(N, E_solv_meV) ** 2
    return np.pi * _droplet_radius_m(N) ** 2


def _evaporation_per_pickup(k: np.ndarray, N: np.ndarray) -> np.ndarray:
    """Number of He atoms lost when a droplet picks up its k-th molecule.

    Port of MATLAB ``num_evap``. The k>1 branch adds extra evaporation
    from the kinetic energy of the second-and-later dopants.

    `num_evap` returns negative numbers because the chemical potential ``mu``
    is negative and the evaporation reduces N -- matches MATLAB.

    Parameters
    ----------
    k : array
        Pickup count for each droplet (1-based, integer).
    N : array
        Current droplet size for each droplet.

    Returns
    -------
    delta_N : array
        Change in N due to evaporation (negative values).
    """
    # chemical potential in Kelvin (Yang et al., 2022)
    mu = -7.21 + 17.71 * N ** (-1.0 / 3.0) - 5.95 * N ** (-2.0 / 3.0)
    T_K = _PICKUP_GAS_TEMP_K
    base = T_K / mu
    extra = (k > 1) * (T_K * 5.0 / 2.0) / mu
    return base + extra


def _simulate_pickup(
    N_raw: np.ndarray,
    *,
    num_target: int,
    reduced_crosssection: bool,
    rng: np.random.Generator,
    E_solv_meV: float = _DEFAULT_E_SOLV_meV,
    return_diagnostics: bool = False,
) -> np.ndarray | tuple[np.ndarray, dict]:
    """Run the pickup-cell simulation on raw log-normal samples.

    Each iteration represents one "round trip" through a portion of the
    pickup region. Droplets that have picked up exactly one molecule are
    moved to the "completed" list; droplets that get destroyed (N<=0
    after evaporation) are dropped; the rest continue to the next round.

    Parameters
    ----------
    N_raw : np.ndarray
        Raw log-normal samples to feed into the pickup-cell simulation.
    num_target : int
        Number of one-pickup droplets to return.
    reduced_crosssection : bool
        Use kinetic-energy-dissipation-reduced cross section.
    rng : np.random.Generator
        Reproducible RNG.
    E_solv_meV : float, optional
        Dopant solvation energy in meV. Default 14.
    return_diagnostics : bool, optional
        If True, return a tuple ``(samples, diagnostics)`` where diagnostics
        is a dict with per-round statistics + raw + post-pickup arrays.
        Used by :func:`diagnose_pickup` for plotting and debugging.
        Default False (returns only the samples array, like before).

    Returns
    -------
    np.ndarray
        ``num_target`` one-pickup droplet sizes.
    diagnostics : dict (only if return_diagnostics=True)
        Keys:
        - ``"N_raw"`` -- input array
        - ``"N_one_pickup"`` -- all single-pickup droplets (not just first num_target)
        - ``"N_no_pickup"`` -- all droplets that never picked anything up
        - ``"per_round"`` -- list of dicts with keys "iter", "n_alive",
          "mean_p_pickup", "mean_evap", "n_picked_up", "n_destroyed".
    """
    n_gas = _PICKUP_GAS_PRESSURE_MBAR * 100.0 / (_PICKUP_GAS_TEMP_K * _K_B)
    a = _PICKUP_REGION_LENGTH_M

    # Working arrays (mutated each round)
    samples = np.asarray(N_raw, dtype=float).copy()
    total_pickup = np.zeros_like(samples)

    completed_samples: list[np.ndarray] = []
    completed_pickups: list[np.ndarray] = []

    per_round: list[dict] = []

    iter_id = 1
    max_pickup = 1
    max_iters = 100   # safety guard

    while samples.size > 0 and total_pickup.max(initial=0) < max_pickup:
        if iter_id > max_iters:
            break

        # 1. Probability of pickup this round (matches MATLAB)
        sigma = _pickup_cross_section_m2(
            samples, reduced=reduced_crosssection, E_solv_meV=E_solv_meV,
        )
        p_pickup = a * n_gas * sigma

        pickup_event = rng.random(samples.shape) < p_pickup
        total_pickup = total_pickup + pickup_event

        # 2. Evaporate He atoms from droplets that just picked up
        delta_N = _evaporation_per_pickup(total_pickup, samples) * pickup_event
        new_samples = samples + delta_N

        # Diagnostic capture (cheap; happens before mutation)
        if return_diagnostics:
            per_round.append({
                "iter": iter_id,
                "n_alive_before": int(samples.size),
                "mean_p_pickup": float(p_pickup.mean()),
                "mean_evap": float(delta_N[pickup_event].mean()) if pickup_event.any() else 0.0,
                "n_picked_up": int(pickup_event.sum()),
            })

        # 3. Drop droplets that got destroyed (N <= 0)
        alive = new_samples > 0
        n_destroyed = int((~alive).sum())
        if return_diagnostics and per_round:
            per_round[-1]["n_destroyed"] = n_destroyed
        new_samples = new_samples[alive]
        total_pickup = total_pickup[alive]

        # 4. Move droplets that did *not* pick up this round to "completed"
        no_pickup_this_round = total_pickup == iter_id - 1
        completed_samples.append(new_samples[no_pickup_this_round])
        completed_pickups.append(total_pickup[no_pickup_this_round])

        # 5. Carry the rest forward
        new_samples = new_samples[~no_pickup_this_round]
        total_pickup = total_pickup[~no_pickup_this_round]

        # Diagnostic capture of per-round new_samples distribution
        # (matches MATLAB's debug_pickup_plot: plotted AFTER destroyed and
        # no-pickup droplets are removed -- it's the droplets that just
        # picked up at least one molecule this round, after evaporation)
        if return_diagnostics:
            per_round[-1]["new_samples"] = new_samples.copy()

        samples = new_samples
        iter_id += 1

    # Catch any remaining samples that hit max_pickup during the last round.
    completed_samples.append(samples)
    completed_pickups.append(total_pickup)

    all_samples = np.concatenate(completed_samples)
    all_pickups = np.concatenate(completed_pickups)

    one_pickup = all_samples[all_pickups == 1]
    no_pickup = all_samples[all_pickups == 0]

    if one_pickup.size < num_target and not return_diagnostics:
        # Production mode: refuse to silently truncate.
        raise ValueError(
            f"only {one_pickup.size} droplets picked up exactly one molecule "
            f"out of {N_raw.size} initial samples (target {num_target}); "
            "increase oversampling factor or use mode='raw'."
        )

    # In diagnostics mode we may return fewer than num_target samples; the
    # caller is investigating *why* yield is low, not consuming the samples.
    result = one_pickup[: min(num_target, one_pickup.size)]

    if return_diagnostics:
        diagnostics = {
            "N_raw": np.asarray(N_raw),
            "N_one_pickup": one_pickup,
            "N_no_pickup": no_pickup,
            "per_round": per_round,
        }
        return result, diagnostics

    return result
