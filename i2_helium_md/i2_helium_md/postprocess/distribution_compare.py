"""Integer-support Wasserstein size-distribution comparison (Tier-2 Phase E, Slice E4).

Scores the E1 simulated I+He_n size distribution
(:class:`~i2_helium_md.postprocess.size_distribution.ShellDistribution`) against
the E3 experimental abundance reference
(:class:`~i2_helium_md.postprocess.abundance_loader.HeAbundanceReference`).

On unit-spaced integer support the 1-D Wasserstein-1 distance reduces to the
CDF-gap sum

    W1 = sum_n | F_sim(n) - F_ref(n) |,

with ``F`` the cumulative fraction over the **union** integer support
(zero-filled where a side is absent -- this is where the simulation's legal
``n = 21`` meets the reference's ``0..20`` support). Units: He atoms (unit rung
spacing). Pure numpy: a 1-D integer-support W1 is exact and trivial. This is a
style/precision choice, not a dependency constraint (scipy is already a hard
dependency); the ``scipy.stats.wasserstein_distance`` cross-check lives in tests
only.

``compare_size_distributions`` is the ``validation_histogram_metric`` config
dispatch -- the field's first physics-live reader. The unbuilt ``chi2`` / ``ks``
arms get point-of-use ``NotImplementedError`` refusals (the established
unbuilt-enum-arm convention).
"""

from __future__ import annotations

import numpy as np

from .abundance_loader import HeAbundanceReference
from .size_distribution import ShellDistribution

# Both producers normalize by their own sum (E1 counts/sum, E3 percent/sum),
# so a genuine fraction vector sums to 1 to a few ulp. The band is deliberately
# far looser than that -- it exists to catch misuse that is off by a *factor*
# (passing .counts or raw percent instead of the fraction), never to reject a
# legitimately rounded input.
_FRACTION_SUM_TOL: float = 1e-6


def _validated_distribution(n, fraction, *, side: str) -> tuple[np.ndarray, np.ndarray]:
    """Coerce and validate one side's (support, fraction) pair.

    Raises
    ------
    ValueError
        Empty support; support/fraction shape mismatch; duplicate support
        values (the union scatter is last-write-wins and would silently drop
        probability mass); non-finite or negative fractions; or a fraction
        vector that does not sum to 1 (e.g. ``.counts`` passed by mistake).
    """
    n_arr = np.asarray(n, dtype=int)
    f_arr = np.asarray(fraction, dtype=float)
    if n_arr.size == 0:
        raise ValueError("empty support: both distributions must be non-empty.")
    if n_arr.shape != f_arr.shape:
        raise ValueError(
            f"{side} support and fraction must share one shape; got "
            f"{n_arr.shape} and {f_arr.shape}."
        )
    if np.unique(n_arr).size != n_arr.size:
        raise ValueError(
            f"{side} support has duplicate n values: {n_arr.tolist()} -- "
            f"probability mass would be silently dropped."
        )
    if not np.all(np.isfinite(f_arr)):
        raise ValueError(
            f"{side} fraction has non-finite entries (NaN/inf): {f_arr}."
        )
    if np.any(f_arr < 0.0):
        raise ValueError(
            f"{side} fraction has negative entries: {f_arr[f_arr < 0.0]}."
        )
    total = float(f_arr.sum())
    if abs(total - 1.0) > _FRACTION_SUM_TOL:
        raise ValueError(
            f"{side} fraction must sum to 1 (got {total}); pass the normalized "
            f"fraction vector (ShellDistribution.fraction / "
            f"HeAbundanceReference.ion_fraction), not counts or percent."
        )
    return n_arr, f_arr


def wasserstein_integer_support(
    sim: ShellDistribution,
    ref: HeAbundanceReference,
) -> float:
    """1-D Wasserstein-1 distance between two integer-support distributions.

    Parameters
    ----------
    sim
        Simulated distribution -- read via ``sim.n_values`` (int support) and
        ``sim.fraction`` (sums to 1).
    ref
        Experimental reference -- read via ``ref.n`` (int support) and
        ``ref.ion_fraction`` (sums to 1).

    Returns
    -------
    float
        ``W1 = sum_n |F_sim(n) - F_ref(n)|`` over the contiguous union support,
        in He atoms (unit rung spacing).

    Raises
    ------
    ValueError
        If either support is empty or carries duplicates; either fraction
        vector is non-finite, negative, or does not sum to 1 (e.g. ``.counts``
        passed by mistake); or the two supports are disjoint (no shared ``n``
        -- a units/labelling bug, since both are indexed by shell count).
    """
    sim_n, sim_f = _validated_distribution(sim.n_values, sim.fraction, side="sim")
    ref_n, ref_f = _validated_distribution(ref.n, ref.ion_fraction, side="ref")

    if not np.intersect1d(sim_n, ref_n).size:
        raise ValueError(
            f"disjoint integer support: sim n in [{sim_n.min()}, {sim_n.max()}], "
            f"ref n in [{ref_n.min()}, {ref_n.max()}] share no shell count -- "
            f"likely a units/labelling bug."
        )

    # Contiguous union support guarantees unit spacing, so the CDF-gap sum is
    # exactly W1 (no per-bin spacing factor needed).
    n_lo = int(min(sim_n.min(), ref_n.min()))
    n_hi = int(max(sim_n.max(), ref_n.max()))
    support = np.arange(n_lo, n_hi + 1)

    f_sim = np.zeros(support.shape, dtype=float)
    f_ref = np.zeros(support.shape, dtype=float)
    f_sim[sim_n - n_lo] = sim_f
    f_ref[ref_n - n_lo] = ref_f

    cdf_sim = np.cumsum(f_sim)
    cdf_ref = np.cumsum(f_ref)
    return float(np.sum(np.abs(cdf_sim - cdf_ref)))


def compare_size_distributions(
    sim: ShellDistribution,
    ref: HeAbundanceReference,
    *,
    metric: str,
) -> float:
    """Score ``sim`` against ``ref`` with the config-selected metric.

    This is the ``SimConfig.validation_histogram_metric`` dispatch point: the
    caller passes ``metric = cfg.validation_histogram_metric``.

    Parameters
    ----------
    sim, ref
        The E1 simulated distribution and the E3 experimental reference.
    metric
        One of ``"wasserstein"`` (built), ``"chi2"`` / ``"ks"`` (declared enum
        arms, not yet implemented).

    Returns
    -------
    float
        The metric value (currently only ``wasserstein``).

    Raises
    ------
    NotImplementedError
        For the ``chi2`` / ``ks`` arms (point-of-use refusal).
    ValueError
        For any value outside the three declared enum arms.
    """
    if metric == "wasserstein":
        return wasserstein_integer_support(sim, ref)
    if metric == "chi2":
        raise NotImplementedError(
            "validation_histogram_metric='chi2' is not implemented; Tier 2 uses "
            "'wasserstein'."
        )
    if metric == "ks":
        raise NotImplementedError(
            "validation_histogram_metric='ks' is not implemented; Tier 2 uses "
            "'wasserstein'."
        )
    raise ValueError(
        f"unknown validation_histogram_metric {metric!r}; expected one of "
        f"'wasserstein', 'chi2', 'ks'."
    )
