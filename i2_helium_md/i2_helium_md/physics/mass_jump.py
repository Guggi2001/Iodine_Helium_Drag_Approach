"""Mass-shed operators for Tier-1a anchored validation.

Pure, stateless physics: a single He shed consumes a pre-shed velocity and mass
and emits the post-shed velocity, post-shed mass, and mechanical ledger term. It
performs no scheduling, integration, or drag evaluation; those live in
``shell_schedule.py`` and the ion driver.

The production Tier-1a ``anchored_discrete`` path uses a continuous-velocity
shed. The shed He leaves co-moving with the tracked complex at the instant of
bookkeeping::

    v+  = v-
    m+  = m - m_He
    dE_mass_transfer = +0.5 * m_He * |v-|^2

The tracked complex loses the kinetic energy carried away by the removed He, so
the positive ledger term keeps ``E_kin + E_pot + E_dissip + E_mass_transfer``
closed.

The older cold-shed reset is retained only as a diagnostic upper bound. It
removes one He atom at rest, so the remaining complex receives the scalar kick
``m/(m - m_He)`` and the ledger term is the negative of that reset's kinetic
energy rise.

Units / mass contract
---------------------
Mechanical amu throughout, matching ``baoab.py`` and ``shell_schedule.py``: ``m``
in amu, ``v`` in A/ps, energy in ``amu*A^2/ps^2``. There is no kg and no eV here --
the eV conversion of the ledger term is owned downstream (the ledger slice), the
same split ``baoab.py`` uses for ``dE_dissip``. This module takes ``m`` only as the
*shed* quantity; it is otherwise mass-agnostic and never sees the drag law.

Run modes
---------
:func:`apply_shed` selects between the two Tier-1a A/B modes by a **function
argument**, deliberately *not* the ``SimConfig.mass_scenario`` enum (that
config-surface change is the integrator-wiring slice's work):

* ``"fixed"`` -- the null: no reset, mass held, zero defect.
* ``"anchored_discrete"`` -- the continuous-velocity shed above.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .constants import MASS_HE_AMU

# The two Tier-1a A/B run modes, passed explicitly (NOT SimConfig.mass_scenario;
# the config-surface enum change is deferred to the integrator-wiring slice).
ShedMode = Literal["fixed", "anchored_discrete"]


@dataclass(frozen=True)
class ShedResult:
    """The outcome of one (possibly null) mass-shed applied to a velocity.

    Attributes
    ----------
    v_plus : np.ndarray
        Post-shed velocity [A/ps], same shape as the input ``v_minus``.
    m_plus_amu : float
        Post-shed complex mass [amu]: ``m - m_He`` under a real shed, the input
        mass under ``fixed``. This is the ``m+`` the post-jump O-step (SQ3) reads.
    dE_mass_transfer : float
        The ledger increment [amu*A^2/ps^2]. For the production continuous shed
        this is ``+0.5*m_He*|v-|^2``; for the cold-shed bound this is the negative
        of the reset's kinetic-energy rise.
    """

    v_plus: np.ndarray
    m_plus_amu: float
    dE_mass_transfer: float


def kick_factor(m_minus_amu: float, m_he_amu: float = MASS_HE_AMU) -> float:
    """Momentum-conserving speed kick ``m/(m - m_He)`` of one cold shed [dimensionless].

    The factor by which the reset scales the speed (``|v+| = kick*|v-|``); equals
    :attr:`~i2_helium_md.physics.shell_schedule.ShedEvent.kick_factor`. Raises
    :class:`ValueError` unless ``m_he_amu > 0`` and ``m_minus_amu > m_he_amu`` (so
    the post-shed mass is strictly positive).
    """
    _check_masses(m_minus_amu, m_he_amu)
    return m_minus_amu / (m_minus_amu - m_he_amu)


def _reduced_mass_defect_coeff(
    m_minus_amu: float, m_plus_amu: float, m_he_amu: float
) -> float:
    """Reduced-mass energy-defect coefficient ``-0.5*(m*m_He)/(m - m_He)`` [amu].

    The single source of the exact reduced-mass form (NOT the heavy-ion
    ``0.5*m_He`` approximation). Multiplying by ``|v-|^2`` gives the (negative)
    ledger increment ``dE_mass_transfer``. Shared by the scalar :func:`cold_shed`
    and the per-atom :func:`cold_shed_velocity_components` so the formula lives
    in exactly one place (CLAUDE.md rule 1).
    """
    return -0.5 * (m_minus_amu * m_he_amu) / m_plus_amu


def cold_shed(
    v_minus,
    m_minus_amu: float,
    *,
    m_he_amu: float = MASS_HE_AMU,
) -> ShedResult:
    """Apply one momentum-conserving cold shed to a velocity (SQ2).

    Parameters
    ----------
    v_minus : array_like
        Pre-shed velocity vector [A/ps]. ``|v-|^2`` is the full sum of squares of
        the components, so a 3-vector ``(vx, vy, vz)`` is the natural input; the
        kick is a scalar multiply, so any shape is accepted and returned unchanged.
    m_minus_amu : float
        Pre-shed complex mass [amu] (e.g. a Slice-S ``mass_before_amu``). Must be
        finite and ``> m_he_amu``.
    m_he_amu : float, optional
        Shed He mass [amu] (default :data:`~i2_helium_md.physics.constants.MASS_HE_AMU`).
        Must be ``> 0``.

    Returns
    -------
    ShedResult
        ``v_plus = m/(m - m_He) * v_minus``, ``m_plus_amu = m - m_He``, and
        ``dE_mass_transfer = -0.5*(m*m_He)/(m - m_He)*|v-|^2`` (``<= 0``).

    Raises
    ------
    ValueError
        If the masses violate ``m_he_amu > 0 < (m_minus_amu - m_he_amu)``, or if
        ``v_minus`` / ``m_minus_amu`` contain non-finite values.
    """
    _check_masses(m_minus_amu, m_he_amu)
    v = np.asarray(v_minus, dtype=float)
    if not np.all(np.isfinite(v)):
        raise ValueError(f"v_minus must be finite; got {v_minus!r}.")

    m_plus = m_minus_amu - m_he_amu
    kick = m_minus_amu / m_plus
    v_plus = kick * v

    speed_sq = float(v @ v) if v.ndim else float(v) ** 2
    # Exact KE rise of the reset (reduced-mass form, NOT heavy-ion 0.5*m_He*v^2);
    # booked negative so the four-term ledger closes by construction.
    dE_mass_transfer = _reduced_mass_defect_coeff(m_minus_amu, m_plus, m_he_amu) * speed_sq

    return ShedResult(v_plus=v_plus, m_plus_amu=m_plus, dE_mass_transfer=dE_mass_transfer)


def continuous_velocity_shed(
    v_minus,
    m_minus_amu: float,
    *,
    m_he_amu: float = MASS_HE_AMU,
) -> ShedResult:
    """Apply one production Tier-1a co-moving He shed to a velocity.

    The velocity is copied unchanged, the complex mass drops by one He, and the
    mass-transfer channel receives the positive kinetic energy carried away by
    the removed co-moving He.
    """
    _check_masses(m_minus_amu, m_he_amu)
    v = np.asarray(v_minus, dtype=float)
    if not np.all(np.isfinite(v)):
        raise ValueError(f"v_minus must be finite; got {v_minus!r}.")

    m_plus = m_minus_amu - m_he_amu
    speed_sq = float(np.sum(v ** 2))
    dE_mass_transfer = 0.5 * m_he_amu * speed_sq
    return ShedResult(
        v_plus=v.copy(),
        m_plus_amu=m_plus,
        dE_mass_transfer=dE_mass_transfer,
    )


def cold_shed_velocity_components(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    m_minus_amu: float,
    *,
    m_he_amu: float = MASS_HE_AMU,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
    """Per-atom vectorized cold shed over ``(2N,)`` velocity component arrays.

    The driver-facing form of :func:`cold_shed`: the same momentum-conserving
    reset applied independently to every atom's velocity, with a **uniform**
    pre-shed complex mass ``m_minus_amu`` (the anchored Tier-1a schedule sheds
    one He from every ion at the same instant, so the kick factor and ``m+`` are
    scalar across atoms; only the per-atom ``|v_i|^2`` -- and hence the defect --
    differs). Reuses :func:`kick_factor` and :func:`_reduced_mass_defect_coeff`,
    so the physics is identical to the scalar :func:`cold_shed` (the single-atom
    oracle in the tests); this routine only vectorizes its application.

    Parameters
    ----------
    vx, vy, vz : np.ndarray, shape (2N,)
        Pre-shed velocity components [A/ps].
    m_minus_amu : float
        Uniform pre-shed complex mass [amu]; must be finite and ``> m_he_amu``.
    m_he_amu : float, optional
        Shed He mass [amu] (default :data:`MASS_HE_AMU`).

    Returns
    -------
    (vx_plus, vy_plus, vz_plus, m_plus_amu, dE_mass_transfer) : tuple
        Kicked components (each ``(2N,)``), the scalar post-shed mass
        ``m+ = m - m_He`` [amu], and the per-atom ledger increment ``(2N,)``
        ``-0.5*(m*m_He)/(m - m_He)*|v_i|^2`` [amu*A^2/ps^2], ``<= 0``.

    Raises
    ------
    ValueError
        If the masses violate ``m_he_amu > 0 < (m_minus_amu - m_he_amu)``.
    """
    _check_masses(m_minus_amu, m_he_amu)
    m_plus = m_minus_amu - m_he_amu
    kick = m_minus_amu / m_plus
    vxf = np.asarray(vx, dtype=float)
    vyf = np.asarray(vy, dtype=float)
    vzf = np.asarray(vz, dtype=float)
    vx_plus = kick * vxf
    vy_plus = kick * vyf
    vz_plus = kick * vzf
    speed_sq = vxf ** 2 + vyf ** 2 + vzf ** 2
    dE_mass_transfer = _reduced_mass_defect_coeff(m_minus_amu, m_plus, m_he_amu) * speed_sq
    return vx_plus, vy_plus, vz_plus, m_plus, dE_mass_transfer


def continuous_velocity_shed_components(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    m_minus_amu: float,
    *,
    m_he_amu: float = MASS_HE_AMU,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
    """Per-atom vectorized production shed over velocity component arrays.

    Velocity components are copied unchanged. The scalar post-shed mass is
    ``m - m_He`` and the per-atom ledger increment is
    ``+0.5*m_He*|v_i|^2`` [amu*A^2/ps^2].
    """
    _check_masses(m_minus_amu, m_he_amu)
    vxf = np.asarray(vx, dtype=float)
    vyf = np.asarray(vy, dtype=float)
    vzf = np.asarray(vz, dtype=float)
    if not (np.all(np.isfinite(vxf)) and np.all(np.isfinite(vyf)) and np.all(np.isfinite(vzf))):
        raise ValueError("velocity components must be finite.")

    m_plus = m_minus_amu - m_he_amu
    speed_sq = vxf ** 2 + vyf ** 2 + vzf ** 2
    dE_mass_transfer = 0.5 * m_he_amu * speed_sq
    return vxf.copy(), vyf.copy(), vzf.copy(), m_plus, dE_mass_transfer


def apply_shed(
    v_minus,
    m_minus_amu: float,
    *,
    mode: ShedMode,
    m_he_amu: float = MASS_HE_AMU,
) -> ShedResult:
    """Apply a shed under the selected Tier-1a A/B mode.

    Parameters
    ----------
    v_minus : array_like
        Pre-shed velocity [A/ps].
    m_minus_amu : float
        Current complex mass [amu].
    mode : {"fixed", "anchored_discrete"}
        ``"fixed"`` -- the null: returns an unchanged velocity copy, the mass held,
        zero defect (no mass relationship is required). ``"anchored_discrete"`` --
        delegates to :func:`continuous_velocity_shed`.
    m_he_amu : float, optional
        Shed He mass [amu]; used only by ``"anchored_discrete"``.

    Returns
    -------
    ShedResult

    Raises
    ------
    ValueError
        If ``mode`` is unknown, or (under ``"fixed"``) ``v_minus`` / ``m_minus_amu``
        are non-finite, or (under ``"anchored_discrete"``) the production shed
        guards trip.
    """
    if mode == "anchored_discrete":
        return continuous_velocity_shed(v_minus, m_minus_amu, m_he_amu=m_he_amu)
    if mode == "fixed":
        v = np.asarray(v_minus, dtype=float)
        if not np.all(np.isfinite(v)) or not np.isfinite(m_minus_amu):
            raise ValueError(
                f"fixed-mode inputs must be finite; got v_minus={v_minus!r}, "
                f"m_minus_amu={m_minus_amu!r}."
            )
        return ShedResult(v_plus=v.copy(), m_plus_amu=float(m_minus_amu),
                          dE_mass_transfer=0.0)
    raise ValueError(
        f"unknown shed mode {mode!r}; expected 'fixed' or 'anchored_discrete'."
    )


def _check_masses(m_minus_amu: float, m_he_amu: float) -> None:
    """Fail loudly unless ``m_he_amu > 0`` and ``m_minus_amu - m_he_amu > 0``."""
    if not (np.isfinite(m_minus_amu) and np.isfinite(m_he_amu)):
        raise ValueError(
            f"masses must be finite; got m_minus_amu={m_minus_amu!r}, "
            f"m_he_amu={m_he_amu!r}."
        )
    if not (m_he_amu > 0.0):
        raise ValueError(f"m_he_amu must be > 0; got {m_he_amu!r}.")
    if not (m_minus_amu - m_he_amu > 0.0):
        raise ValueError(
            f"pre-shed mass must exceed the He mass so m+ = m - m_He > 0; got "
            f"m_minus_amu={m_minus_amu!r}, m_he_amu={m_he_amu!r}."
        )
