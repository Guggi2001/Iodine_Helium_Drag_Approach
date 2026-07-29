"""Drag-shell state coupling s(n) (Tier-2 atlas §3.5i; design doc S1/S2).

The dimensionless per-ion state factor that lets the drag learn the ion has
stripped: ``gamma(v, d, n) = g(d) * s(n) * gamma_form(v)``
(``TIER2_DRAG_STATE_COUPLING_DESIGN.md`` §3). The moving object is the
snowball currently carrying ``n`` He; its ram/displacement coupling scales
with geometric cross-section, hence the two-parameter closure (OQ-B)::

    s(n)     = ( R_eff(n) / R_eff(n_ref) )**2          # dimensionless
    R_eff(n) = ( R_core**3 + 3*n / (4*pi*rho_shell) )**(1/3)   # Angstrom

Units: ``R_core**3`` [A^3] + ``n/rho_shell`` [A^3] -> ``R_eff`` [A] -> ``s``
dimensionless. The exponent 2 is fixed by geometry (area scaling); no free
exponent (OQ-B).

Architecture (design §5, BC-2): the drag module ``physics/drag.py`` stays
mass-agnostic and **state-blind** -- this module never touches it. The
*driver* computes ``s(n(t))`` from the per-ion shell state it already tracks
and multiplies the ``gamma_fn`` closure result (:func:`apply_state_factor`),
exactly as the spatial gate composes. Because ``s`` scales ``gamma`` itself,
the (Tier-3, stubbed) FDT noise amplitude ``sqrt(2*gamma*kB*T_eff)`` inherits
it automatically (OQ-E) and the drag power ``P = gamma*v**2`` feeds the
5-term invariant unchanged in structure.

Normalization (OQ-C): ``n_ref`` is **derived** from the coefficient bundle's
extraction-mass stamp -- the shell state the TDDFT calibration measured
(202.953908 amu = I+ + 19 He under the biphasic ``m(n)`` convention
``MASS_I_ION_AMU + n * MASS_HE_AMU``). Any other normalization is a stealth
rescale of the calibrated ``b``. :func:`derive_n_ref_amu` is the single
source for both the config-load guard and the driver seam (rule 1).
"""

from __future__ import annotations

import numpy as np

from .constants import MASS_HE_AMU, MASS_I_ION_AMU


#: Tolerance [shell counts] for the n_ref-vs-stamp integer consistency check.
#: The standing bundle stamp (202.953908 amu) derives to 19.0011 under the
#: constants.py biphasic mass convention (the stamp was written with
#: higher-precision atomic masses); 0.01 accepts that mass-table slack while
#: still refusing any genuinely non-integer (mis-stamped) bundle.
N_REF_INTEGER_TOL: float = 0.01


def derive_n_ref_amu(extraction_mass_amu: float) -> int:
    """Derive the s(n) normalization shell count from a bundle's mass stamp.

    ``n_ref = round((extraction_mass_amu - MASS_I_ION_AMU) / MASS_HE_AMU)``
    under the biphasic mass convention (constants.py: mass at shell count n is
    ``MASS_I_ION_AMU + n * MASS_HE_AMU``). Consistency-checked, not free
    (OQ-C): the derivation must land within :data:`N_REF_INTEGER_TOL` of an
    integer >= 1, else the bundle's stamp does not name a shell state and the
    coupling has no defensible normalization.

    Parameters
    ----------
    extraction_mass_amu : float
        The bundle's ``extraction_mass_amu`` provenance stamp [amu].

    Returns
    -------
    int
        The normalization shell count (19 for the standing shared bundle).

    Raises
    ------
    ValueError
        If the stamp is not positive, or does not derive to a near-integer
        shell count >= 1.
    """
    if not (extraction_mass_amu > 0):
        raise ValueError(
            f"extraction_mass_amu must be positive, got {extraction_mass_amu!r}"
        )
    derived = (float(extraction_mass_amu) - MASS_I_ION_AMU) / MASS_HE_AMU
    n_ref = int(round(derived))
    if abs(derived - n_ref) > N_REF_INTEGER_TOL or n_ref < 1:
        raise ValueError(
            f"extraction_mass_amu={extraction_mass_amu!r} derives to "
            f"{derived:.4f} He under m(n) = MASS_I_ION_AMU + n*MASS_HE_AMU -- "
            f"not within {N_REF_INTEGER_TOL} of an integer shell count >= 1, "
            "so s(n) has no defensible n_ref normalization (OQ-C: n_ref is "
            "derived from the stamp, never free)."
        )
    return n_ref


def effective_radius_angstrom(
    n, *, R_core_angstrom: float, rho_shell_per_A3: float,
) -> np.ndarray:
    """Effective collision radius ``R_eff(n)`` [Angstrom] of the dressed snowball.

    ``R_eff(n) = (R_core**3 + 3*n / (4*pi*rho_shell))**(1/3)``: the bare-core
    volume plus the attached-shell volume of ``n`` He at shell density
    ``rho_shell``, read as an equivalent sphere.

    Parameters
    ----------
    n : array_like
        Shell count(s), >= 0 (float-valued state accepted; the biphasic
        ``n_shell`` state is integer-valued but stored as float).
    R_core_angstrom : float
        Bare-core effective collision radius [A], > 0 (Bounded, prior 3.2).
    rho_shell_per_A3 : float
        Attached-shell He number density [A^-3], > 0 (Bounded, prior 0.030).

    Returns
    -------
    np.ndarray
        ``R_eff`` [Angstrom], shape of ``n``.

    Raises
    ------
    ValueError
        On non-positive/non-finite parameters or any negative ``n``.
    """
    if not (np.isfinite(R_core_angstrom) and R_core_angstrom > 0):
        raise ValueError(
            f"R_core_angstrom must be finite and > 0, got {R_core_angstrom!r}"
        )
    if not (np.isfinite(rho_shell_per_A3) and rho_shell_per_A3 > 0):
        raise ValueError(
            f"rho_shell_per_A3 must be finite and > 0, got {rho_shell_per_A3!r}"
        )
    n = np.asarray(n, dtype=float)
    if np.any(n < 0):
        raise ValueError("shell count n must be >= 0 everywhere")
    return (
        R_core_angstrom ** 3 + 3.0 * n / (4.0 * np.pi * rho_shell_per_A3)
    ) ** (1.0 / 3.0)


def shell_area_state_factor(
    n,
    *,
    R_core_angstrom: float,
    rho_shell_per_A3: float,
    n_ref: int,
) -> np.ndarray:
    """Dimensionless drag state factor ``s(n) = (R_eff(n)/R_eff(n_ref))**2``.

    ``s(n_ref) = 1`` exactly (the bundle's calibrated state is the identity);
    ``s < 1`` for a stripped ion, ``s > 1`` for an over-dressed one (n = 21
    early-window: ~1.06 at the priors -- deliberately uncapped, the closure
    is monotone in n by construction).

    Parameters
    ----------
    n : array_like
        Per-ion shell count(s), >= 0.
    R_core_angstrom, rho_shell_per_A3 : float
        The geometric closure parameters (see
        :func:`effective_radius_angstrom`).
    n_ref : int
        Normalization shell count, >= 1 -- **derived** from the bundle stamp
        via :func:`derive_n_ref_amu`, never free (OQ-C).

    Returns
    -------
    np.ndarray
        ``s(n)`` dimensionless, > 0, shape of ``n``.

    Raises
    ------
    ValueError
        Propagated from :func:`effective_radius_angstrom`, or on
        ``n_ref < 1``.
    """
    if int(n_ref) < 1:
        raise ValueError(f"n_ref must be an integer >= 1, got {n_ref!r}")
    r_eff = effective_radius_angstrom(
        n, R_core_angstrom=R_core_angstrom, rho_shell_per_A3=rho_shell_per_A3,
    )
    r_ref = effective_radius_angstrom(
        float(int(n_ref)),
        R_core_angstrom=R_core_angstrom, rho_shell_per_A3=rho_shell_per_A3,
    )
    return (r_eff / r_ref) ** 2


def apply_state_factor(gamma_fn, s):
    """Compose the per-ion state factor onto a ``gamma_fn`` closure (BC-2 seam).

    Returns a new ``GammaFn`` ``(speed, depth) -> s * gamma_fn(speed, depth)``
    with ``s`` frozen at the per-ion values of the current step (the driver
    rebuilds this wrapper every step from the live ``n_shell``, jump-then-O:
    the same post-event state the step's ``m(t)`` reads -- BC-1). The drag
    module itself is never touched; under ``drag_state_coupling='off'`` the
    driver passes its base closure verbatim, so the off path is structurally
    bit-identical.

    Parameters
    ----------
    gamma_fn : Callable[[speed, depth], np.ndarray]
        The base closure (``drag_gamma`` partial or the E2 Landau-gated arm).
    s : np.ndarray
        Per-ion state factors, shape broadcastable with the closure's output
        (``(2N,)`` in both drivers).

    Returns
    -------
    Callable[[speed, depth], np.ndarray]
        The scaled closure.
    """
    s = np.asarray(s, dtype=float)

    def _scaled(speed, depth):
        return s * gamma_fn(speed, depth)

    return _scaled
