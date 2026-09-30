"""Exact asymptotic closure of a residual repulsive pair Coulomb interaction.

The detection stage hands over at a finite time ``t_h`` while the two I+
fragments of a molecule still repel each other (TIER2_DetectorStageFix.md
§2c). Outside helium the remaining flight is conservative (flat droplet well,
no drag), so the residual pair energy ends up as kinetic energy at infinity in
a way fixed entirely by the handover state. This module maps a handover state
to that asymptotic velocity in closed form (§2e):

* **two-body** (both fragments free): the centre of mass moves uniformly; the
  relative coordinate follows the repulsive Kepler hyperbola of the reduced
  mass. Energy and angular momentum fix the asymptotic relative speed, the
  Laplace-Runge-Lenz vector fixes its direction, and momentum conservation
  fixes each fragment's share.
* **fixed-centre** (partner held by the droplet): the ``m_partner -> inf``
  limit -- the free ion is repelled from a fixed point charge at the partner's
  position and takes the whole pair energy.

Units (MD units throughout, no eV): positions [A], velocities [A/ps],
masses [amu], coupling ``k`` [amu*A^3/ps^2] (pair potential ``k / r``, so
``k = E_c * r`` converted from eV*A). Pure functions; no config, no RNG.

Dimensional check. ``alpha = k / mu`` [A^3/ps^2]; ``2*alpha/r`` [A^2/ps^2]
adds to ``|v_rel|^2``; ``h = r x v`` [A^2/ps]; ``(v x h)/alpha`` and
``c = w_inf*|h|/alpha`` are dimensionless.
"""

from __future__ import annotations

import numpy as np

__all__ = ["two_body_coulomb_asymptote", "fixed_centre_coulomb_asymptote"]


def _repulsive_kepler_asymptote(r_rel, v_rel, alpha):
    """Asymptotic relative velocity of a repulsive Kepler orbit.

    Parameters
    ----------
    r_rel, v_rel : np.ndarray, shape (P, 3)
        Relative position [A] and velocity [A/ps].
    alpha : np.ndarray, shape (P,)
        ``k / mu`` [A^3/ps^2], strictly positive.

    Returns
    -------
    np.ndarray, shape (P, 3)
        ``w_inf * u_hat`` [A/ps] with ``w_inf^2 = |v|^2 + 2 alpha / r`` and
        ``u_hat = (-e + (w_inf/alpha) * (e x h)) / (1 + c^2)``, where
        ``e = -(v x h)/alpha - r_hat`` (LRL vector of the repulsive orbit,
        ``|e|^2 = 1 + c^2``) and ``c = w_inf |h| / alpha``. Written with
        ``e x h`` instead of ``c * (e x h_hat)`` so the radial orbit
        (``h = 0``) needs no division by ``|h|``: then ``e = -r_hat`` and
        ``u_hat = r_hat``.
    """
    r = np.linalg.norm(r_rel, axis=1)
    r_hat = r_rel / r[:, None]
    h = np.cross(r_rel, v_rel)
    w_inf = np.sqrt(np.sum(v_rel ** 2, axis=1) + 2.0 * alpha / r)
    e = -np.cross(v_rel, h) / alpha[:, None] - r_hat
    c2 = (w_inf * np.linalg.norm(h, axis=1) / alpha) ** 2
    u_hat = (-e + (w_inf / alpha)[:, None] * np.cross(e, h)) / (1.0 + c2)[:, None]
    return w_inf[:, None] * u_hat


def _check_pairs(name, *arrays):
    n = arrays[0].shape[0]
    for a in arrays:
        if a.ndim != 2 or a.shape != (n, 3):
            raise ValueError(f"{name}: vectors must share shape (P, 3); got {a.shape}")


def two_body_coulomb_asymptote(r1, v1, m1_amu, r2, v2, m2_amu, k):
    """Asymptotic velocities of two free point charges repelling as ``k / r``.

    Parameters
    ----------
    r1, v1, r2, v2 : array_like, shape (P, 3)
        Fragment positions [A] and velocities [A/ps] at handover.
    m1_amu, m2_amu : array_like, shape (P,)
        Fragment masses [amu], held fixed over the closure (handover masses).
    k : array_like, shape (P,)
        Pair coupling [amu*A^3/ps^2], ``>= 0``. ``k = 0`` pairs are returned
        unchanged.

    Returns
    -------
    v1_inf, v2_inf : np.ndarray, shape (P, 3)
        ``V + (m2/M) w_inf u_hat`` and ``V - (m1/M) w_inf u_hat`` [A/ps]
        (``V`` the centre-of-mass velocity). Total momentum is conserved and
        the kinetic-energy gain equals ``k / |r1 - r2|`` exactly.

    Raises
    ------
    ValueError
        On mismatched shapes, non-positive masses, negative ``k`` or
        coincident fragments.
    """
    r1, v1, r2, v2 = (np.asarray(a, dtype=float) for a in (r1, v1, r2, v2))
    _check_pairs("two_body_coulomb_asymptote", r1, v1, r2, v2)
    m1 = np.asarray(m1_amu, dtype=float)
    m2 = np.asarray(m2_amu, dtype=float)
    k = np.asarray(k, dtype=float)
    if np.any(m1 <= 0) or np.any(m2 <= 0):
        raise ValueError("two_body_coulomb_asymptote: masses must be > 0")
    if np.any(k < 0):
        raise ValueError("two_body_coulomb_asymptote: k must be >= 0 (repulsive)")
    r_rel = r1 - r2
    if np.any(np.linalg.norm(r_rel, axis=1) <= 0):
        raise ValueError("two_body_coulomb_asymptote: coincident fragments")

    M = m1 + m2
    mu = m1 * m2 / M
    V = (m1[:, None] * v1 + m2[:, None] * v2) / M[:, None]
    live = k > 0
    out1, out2 = v1.copy(), v2.copy()          # k = 0: verbatim, no round-off
    if np.any(live):
        w = _repulsive_kepler_asymptote(
            r_rel[live], (v1 - v2)[live], k[live] / mu[live])
        out1[live] = V[live] + (m2 / M)[live, None] * w
        out2[live] = V[live] - (m1 / M)[live, None] * w
    return out1, out2


def fixed_centre_coulomb_asymptote(r, v, m_amu, centre, k):
    """Asymptotic velocity of a free charge repelled by a fixed point charge.

    The ``m_partner -> inf`` limit of :func:`two_body_coulomb_asymptote`
    (``mu = m``, centre at rest): used when the partner is retained in the
    droplet, which absorbs its recoil (TIER2_DetectorStageFix.md §2e class B).

    Parameters
    ----------
    r, v : array_like, shape (P, 3)
        Free-ion position [A] and velocity [A/ps] at handover.
    m_amu : array_like, shape (P,)
        Free-ion mass [amu].
    centre : array_like, shape (P, 3)
        Fixed repelling charge position [A] (the partner at handover).
    k : array_like, shape (P,)
        Pair coupling [amu*A^3/ps^2], ``>= 0``; ``k = 0`` returns ``v``.

    Returns
    -------
    np.ndarray, shape (P, 3)
        Asymptotic velocity [A/ps]; the kinetic-energy gain is ``k / |r - centre|``.
    """
    r, v, centre = (np.asarray(a, dtype=float) for a in (r, v, centre))
    _check_pairs("fixed_centre_coulomb_asymptote", r, v, centre)
    m = np.asarray(m_amu, dtype=float)
    k = np.asarray(k, dtype=float)
    if np.any(m <= 0):
        raise ValueError("fixed_centre_coulomb_asymptote: masses must be > 0")
    if np.any(k < 0):
        raise ValueError("fixed_centre_coulomb_asymptote: k must be >= 0 (repulsive)")
    r_rel = r - centre
    if np.any(np.linalg.norm(r_rel, axis=1) <= 0):
        raise ValueError("fixed_centre_coulomb_asymptote: ion at the centre")
    out = v.copy()
    live = k > 0
    if np.any(live):
        out[live] = _repulsive_kepler_asymptote(r_rel[live], v[live], k[live] / m[live])
    return out
