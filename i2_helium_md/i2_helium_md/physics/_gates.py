"""Shared spatial-gate primitive (Tier-2 Phase-B Slice rho).

The single source of the erf-complement gate ``0.5 * (1 - erf(depth/steepness))``,
extracted here so the two callers with legitimately-distinct *roles* share **one
formula in one place** (CLAUDE.md rule 1):

* :func:`i2_helium_md.physics.drag.spatial_gate` -- the drag force-coefficient gate
  ``g(depth)`` (Tier-0, G2; turns drag off outside the droplet);
* :func:`i2_helium_md.physics.helium_density.rho_he_ratio` -- the ``rho_He/rho_bulk``
  occupancy gate that governs pickup (Tier-2 Phase-B, Slice rho).

This is a **neutral** home on purpose: ``helium_density`` must **not** import the
drag law (it is a density/occupancy quantity, independent of ``gamma``), and the
drag module must not import a Phase-B module. Both import this leaf instead. The
module is mass-agnostic and config-agnostic (``depth`` + ``steepness`` only; no
``n``, no mass, no ``gamma``).

The ``spatial_gate`` rewire onto this helper is a Tier-0 **no-behaviour-change**
refactor: the arithmetic and the ``steepness > 0`` guard are byte-identical to the
former in-lined form, locked by a parity regression test
(``tests/test_helium_density.py`` + ``tests/test_drag.py``).
"""

from __future__ import annotations

import numpy as np
from scipy.special import erf


def _erf_complement(depth, steepness: float) -> np.ndarray:
    """Erf-complement gate ``0.5 * (1 - erf(depth/steepness))`` in ``[0, 1]``.

    1 deep inside (``depth << 0``), 0.5 at the nominal surface (``depth = 0``),
    0 outside (``depth >> 0``). Smooth and ``C^1``. The erf saturates *exactly*
    to ``+-1`` in float at large ``|depth|``, so the tails are exact ``1.0`` /
    ``0.0`` (no epsilon). Broadcasts to the shape of ``depth``.

    Parameters
    ----------
    depth : array_like
        ``r_atom - r_droplet`` in Angstrom. Negative inside, positive outside.
    steepness : float
        Width of the erf transition in Angstrom (``> 0``).

    Returns
    -------
    np.ndarray
        Dimensionless gate factor in ``[0, 1]`` (0-d for scalar ``depth``,
        matching the historical ``spatial_gate`` return exactly).

    Raises
    ------
    ValueError
        If ``steepness`` is not positive.
    """
    if not (steepness > 0):
        raise ValueError(f"steepness must be positive, got {steepness!r}")
    depth = np.asarray(depth, dtype=float)
    return 0.5 * (1.0 - erf(depth / steepness))
