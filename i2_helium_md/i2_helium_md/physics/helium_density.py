"""Helium density gate ``rho_He(depth)/rho_bulk in [0, 1]`` (Tier-2 Phase-B Slice rho).

The occupancy field that gates the pickup channel (Slice P): the local He density,
normalised to bulk, as a function of the ion's depth below the droplet surface. It
is the first, fully-independent Phase-B build unit -- pure, config-agnostic, and
mass-agnostic (``depth`` + ``steepness`` only; no ``n``, no Langmuir cap, no
``gamma``, no mass). Slice P *mocks* it; the generative driver (Phase C) threads it.

Encoded form (TIER2_PHASE_B_IMPLEMENTATION_PLAN.md §2.2 (rho))
-------------------------------------------------------------
With ``depth = r_atom - r_droplet`` (negative inside the droplet, positive
outside)::

    rho_He(depth) / rho_bulk = 0.5 * (1 - erf(depth / steepness))    # in [0, 1]

-> 1 deep inside (``depth << 0``), 0.5 at the nominal surface (``depth = 0``),
-> 0 outside (``depth >> 0``, ion exit -> pickup rate lambda -> 0 -> termination).

This reuses the **same** erf-complement form as the drag spatial gate, routed
through the single-source helper :func:`i2_helium_md.physics._gates._erf_complement`
so the density and drag surfaces are one formula in one place (CLAUDE.md rule 1).
rho is kept a *distinct* quantity feeding P -- its role (occupancy/capture gate) is
independent of the drag law, and it carries the declared sourced-profile fallback.

Profiles (config selector ``HeliumDensityProfile``)
---------------------------------------------------
* ``erf_complement`` (default) -- :func:`rho_he_ratio`, the analytic gate above;
* ``tabulated``               -- :func:`tabulated_density_profile`, the declared
  fallback for a *sourced* baseline/TDDFT ``rho_He(r)`` array (CALIBRATION_MAP row 8).
  The interpolation **machinery** is built and round-trip tested here; the concrete
  sourced data array is a deferred rule-2 carry (pinned when the profile exists).

Vectorised: scalar-in -> float-out, array-in -> ndarray-out (the ``shell_schedule`` /
``dissociation_ladder`` idiom).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._gates import _erf_complement


def rho_he_ratio(depth, *, steepness: float) -> float | np.ndarray:
    """``rho_He(depth)/rho_bulk in [0, 1]`` -- the erf-complement occupancy gate.

    Parameters
    ----------
    depth : array_like
        ``r_atom - r_droplet`` in Angstrom (negative inside, positive outside).
    steepness : float
        Width of the erf transition in Angstrom (``> 0``). **No default** --
        mirrors :func:`i2_helium_md.physics.drag.spatial_gate`; the Phase-C driver
        passes ``_drag_gate_steepness(cfg)`` (= ``cfg.potential_steepness`` under the
        production gate) so density and drag see the *same* surface.

    Returns
    -------
    float or np.ndarray
        Dimensionless ratio in ``[0, 1]``: ``float`` for scalar ``depth``,
        ``ndarray`` (shape of ``depth``) otherwise.

    Raises
    ------
    ValueError
        If ``steepness`` is not positive (propagated from ``_erf_complement``).
    """
    out = _erf_complement(depth, steepness)
    return float(out) if np.ndim(depth) == 0 else out


@dataclass(frozen=True)
class TabulatedDensityProfile:
    """A ``rho_He/rho_bulk`` profile from an explicit ``(depth, ratio)`` table.

    The declared sourced-profile fallback (not the default path). ``depth_grid`` is
    strictly increasing (Angstrom); ``ratio_grid`` is the matching ``[0, 1]`` density
    ratio, decreasing from ~1 deep inside to ~0 outside. Lookup is **linear
    interpolation on depth** with the tails **clamped to the endpoint ratios**
    (``ratio_grid[0]`` for ``depth < depth_grid[0]``, ``ratio_grid[-1]`` for
    ``depth > depth_grid[-1]``) -- matching the erf-complement asymptotes so both
    profile arms share one ``[0, 1]`` boundary contract.
    """

    depth_grid: tuple[float, ...]
    ratio_grid: tuple[float, ...]

    def __post_init__(self) -> None:
        """Validate the profile invariants at construction (principle 4).

        Placed here (not only in :func:`tabulated_density_profile`) so a direct
        dataclass construction cannot smuggle an unsorted / out-of-range table past
        ``np.interp`` -- which would return silently-wrong values rather than raise.
        Single source of validation for both the builder and any direct construction.
        """
        depth = np.asarray(self.depth_grid, dtype=float)
        ratio = np.asarray(self.ratio_grid, dtype=float)
        if depth.ndim != 1 or ratio.ndim != 1:
            raise ValueError("depth_grid and ratio_grid must be 1-D")
        if depth.shape != ratio.shape:
            raise ValueError(
                f"depth_grid and ratio_grid must have equal length; got "
                f"{depth.shape} vs {ratio.shape}"
            )
        if depth.size == 0:
            raise ValueError("tabulated density profile requires at least one node")
        # Post-review fix (2026-07-02): NaN passes the [0, 1] range check vacuously
        # (it fails both comparisons) and np.interp degenerates on infinite nodes --
        # both "smuggle past np.interp" cases this class exists to refuse. A sourced
        # TDDFT CSV with one NaN row is the realistic failure mode.
        if not (np.all(np.isfinite(depth)) and np.all(np.isfinite(ratio))):
            raise ValueError(
                "depth_grid and ratio_grid values must be finite (a NaN/inf node "
                "makes np.interp return silently-wrong ratios)"
            )
        if not np.all(np.diff(depth) > 0.0):
            raise ValueError("depth_grid must be strictly increasing")
        if np.any(ratio < 0.0) or np.any(ratio > 1.0):
            raise ValueError("ratio_grid values must lie in [0, 1]")

    def ratio(self, depth):
        """``rho_He/rho_bulk`` at ``depth`` [A] by clamped linear interpolation.

        Scalar-in -> ``float``; array-in -> ``ndarray``. ``np.interp`` clamps to the
        endpoint values outside the grid, which realises the tail-clamp contract.
        """
        xp = np.asarray(self.depth_grid, dtype=float)
        fp = np.asarray(self.ratio_grid, dtype=float)
        out = np.interp(np.asarray(depth, dtype=float), xp, fp)
        return float(out) if np.ndim(depth) == 0 else out


def tabulated_density_profile(depth_grid, ratio_grid) -> TabulatedDensityProfile:
    """Build a :class:`TabulatedDensityProfile` from a ``(depth, ratio)`` table.

    Thin coercion wrapper: accepts any ``array_like`` and constructs the frozen
    profile, whose ``__post_init__`` fail-loud (CLAUDE.md principle 4) on a length
    mismatch, a non-strictly-increasing ``depth_grid``, an empty table, a non-1-D
    input, a non-finite node, or a ``ratio`` value outside ``[0, 1]``. Round-trips
    the supplied nodes exactly and interpolates linearly between them.

    Parameters
    ----------
    depth_grid : array_like
        Strictly increasing depths [A].
    ratio_grid : array_like
        Matching ``rho_He/rho_bulk`` values in ``[0, 1]``.

    Raises
    ------
    ValueError
        On length mismatch, empty table, non-increasing ``depth_grid``, a non-1-D
        input, a non-finite node, or a ratio outside ``[0, 1]`` (all raised by
        ``__post_init__``).
    """
    return TabulatedDensityProfile(
        depth_grid=tuple(np.asarray(depth_grid, dtype=float).tolist()),
        ratio_grid=tuple(np.asarray(ratio_grid, dtype=float).tolist()),
    )
