"""Mass-shed operators for Tier-1a anchored validation.

Pure, stateless physics: a He shed consumes a pre-shed velocity and mass
and emits the post-shed velocity, post-shed mass, and mechanical ledger term. It
performs no scheduling, integration, or drag evaluation; those live in
``shell_schedule.py`` and the ion driver.

The production Tier-1a ``anchored_discrete`` path uses a continuous-velocity
shed. Removed He leaves co-moving with the tracked complex at the instant of
bookkeeping; for the default one-He shed::

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


@dataclass(frozen=True)
class CaptureResult:
    """The outcome of one momentum-conserving He capture applied to a velocity.

    The +He counterpart of :class:`ShedResult`: the complex *gains* mass and a
    strictly-positive kinetic-energy defect (KE lost to the internal reservoir as the
    incoming He merges). The name matches the gain semantics (``m_plus_amu`` is larger,
    ``dE_mass_transfer > 0``), so a capture is never mistaken for a shed at the call site.

    Attributes
    ----------
    v_plus : np.ndarray
        Post-capture centre-of-mass velocity [A/ps], same shape as the input
        ``v_minus``: ``(m*v- + m_He*u_He)/(m + m_He)``.
    m_plus_amu : float
        Post-capture complex mass [amu]: ``m + m_He`` (a **gain**).
    dE_mass_transfer : float
        The ledger increment [amu*A^2/ps^2], ``+0.5*(m*m_He)/(m + m_He)*|v- - u_He|^2``
        (**> 0** for a non-zero relative velocity). This is the KE that leaves the
        mechanical channel and books into ``E_mass_transfer`` at Phase C -- the exact
        reduced-mass counterpart of the (negative) cold-shed defect.
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
    """Unsigned reduced-mass energy-defect coefficient ``0.5*(m*m_He)/m_plus`` [amu].

    The single source of the exact reduced-mass form (NOT the heavy-ion
    ``0.5*m_He`` approximation), returned as an **unsigned magnitude**. Multiplying
    by ``|v-|^2`` gives the magnitude of the ledger increment ``|dE_mass_transfer|``;
    the **sign is applied at the call site** -- ``cold_shed`` negates it (energy is
    booked negative so the shed ledger closes), ``capture`` keeps it positive (KE is
    absorbed into the merged complex). ``m_plus`` is the post-event complex mass:
    ``m - m_He`` for a shed, ``m + m_He`` for a capture, so the same reduced-mass
    coefficient serves both channels with one formula (CLAUDE.md rule 1). Shared by
    :func:`cold_shed`, :func:`cold_shed_velocity_components`, :func:`capture`, and
    :func:`capture_velocity_components`.
    """
    return 0.5 * (m_minus_amu * m_he_amu) / m_plus_amu


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
    # booked negative (the coeff is now an unsigned magnitude -- sign at the call
    # site) so the four-term ledger closes by construction.
    dE_mass_transfer = -_reduced_mass_defect_coeff(m_minus_amu, m_plus, m_he_amu) * speed_sq

    return ShedResult(v_plus=v_plus, m_plus_amu=m_plus, dE_mass_transfer=dE_mass_transfer)


def continuous_velocity_shed(
    v_minus,
    m_minus_amu: float,
    *,
    m_he_amu: float = MASS_HE_AMU,
    n_removed: int = 1,
) -> ShedResult:
    """Apply one production Tier-1a co-moving He shed to a velocity.

    The velocity is copied unchanged, the complex mass drops by ``n_removed`` He,
    and the mass-transfer channel receives the positive kinetic energy carried
    away by the removed co-moving He.
    """
    _check_masses(m_minus_amu, m_he_amu)
    n = _check_removed_count(n_removed, m_minus_amu, m_he_amu)
    v = np.asarray(v_minus, dtype=float)
    if not np.all(np.isfinite(v)):
        raise ValueError(f"v_minus must be finite; got {v_minus!r}.")

    m_plus = m_minus_amu - n * m_he_amu
    speed_sq = float(np.sum(v ** 2))
    dE_mass_transfer = 0.5 * n * m_he_amu * speed_sq
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
    reset applied independently to every atom's velocity. ``m_minus_amu`` may be a
    **uniform scalar** (the anchored Tier-1a schedule sheds one He from every ion at
    the same instant, so the kick factor and ``m+`` are scalar across atoms; only the
    per-atom ``|v_i|^2`` -- and hence the defect -- differs) **or a per-ion array**
    (the Tier-2 generative evaporation channel fires per ion independently, so masses
    diverge -- the symmetric counterpart of :func:`capture_velocity_components`). The
    arithmetic (``m+``, ``kick``, the reduced-mass coefficient) broadcasts either way,
    so the scalar path is byte-identical to before. Reuses
    :func:`_reduced_mass_defect_coeff`, so the physics is identical to the scalar
    :func:`cold_shed` (the single-atom oracle in the tests); this routine only
    vectorizes its application.

    Parameters
    ----------
    vx, vy, vz : np.ndarray, shape (2N,)
        Pre-shed velocity components [A/ps].
    m_minus_amu : float or np.ndarray
        Pre-shed complex mass(es) [amu]; a uniform scalar (Tier-1a) or a per-ion
        ``(2N,)`` array (Tier-2 evaporation). Finite and ``> m_he_amu`` (every element).
    m_he_amu : float, optional
        Shed He mass [amu] (default :data:`MASS_HE_AMU`).

    Returns
    -------
    (vx_plus, vy_plus, vz_plus, m_plus_amu, dE_mass_transfer) : tuple
        Kicked components (each ``(2N,)``), the post-shed mass ``m+ = m - m_He`` [amu]
        (scalar or ``(2N,)`` matching ``m_minus_amu``), and the per-atom ledger
        increment ``(2N,)`` ``-0.5*(m*m_He)/(m - m_He)*|v_i|^2`` [amu*A^2/ps^2], ``<= 0``.

    Raises
    ------
    ValueError
        If any mass violates ``m_he_amu > 0 < (m_minus_amu - m_he_amu)``.
    """
    _check_masses_shed(m_minus_amu, m_he_amu)
    m_plus = m_minus_amu - m_he_amu
    kick = m_minus_amu / m_plus
    vxf = np.asarray(vx, dtype=float)
    vyf = np.asarray(vy, dtype=float)
    vzf = np.asarray(vz, dtype=float)
    vx_plus = kick * vxf
    vy_plus = kick * vyf
    vz_plus = kick * vzf
    speed_sq = vxf ** 2 + vyf ** 2 + vzf ** 2
    # Sign at the call site (the coeff is an unsigned magnitude): shed books negative.
    dE_mass_transfer = -_reduced_mass_defect_coeff(m_minus_amu, m_plus, m_he_amu) * speed_sq
    return vx_plus, vy_plus, vz_plus, m_plus, dE_mass_transfer


def capture(
    v_minus,
    m_minus_amu,
    *,
    m_he_amu: float = MASS_HE_AMU,
    u_he: float = 0.0,
) -> CaptureResult:
    """Apply one momentum-conserving He capture to a velocity (Tier-2 pickup reset).

    The +He counterpart of :func:`cold_shed`: an incoming He atom (velocity ``u_he``,
    zero under the production ``at_rest`` arm) merges with the complex, the centre of
    mass conserves momentum, and the kinetic energy lost to the merge books into the
    positive mass-transfer defect. Shares :func:`_reduced_mass_defect_coeff` with the
    shed resets so the exact reduced-mass form lives in one place (CLAUDE.md rule 1);
    the sign is applied here (``+``).

    Parameters
    ----------
    v_minus : array_like
        Pre-capture complex velocity vector [A/ps]. ``|v- - u_He|^2`` is the full sum
        of squares over the components; any shape is accepted (the kick is a scalar
        multiply on ``v- - u_He``) and returned unchanged.
    m_minus_amu : float
        Pre-capture complex mass [amu]. Must be finite and ``> 0`` (a capture only
        *adds* mass, so -- unlike a shed -- no ``m > m_He`` precondition applies;
        :func:`_check_masses_gain` enforces positivity only).
    m_he_amu : float, optional
        Captured He mass [amu] (default :data:`MASS_HE_AMU`). Must be ``> 0``.
    u_he : float, optional
        Incoming He velocity [A/ps], broadcast against ``v_minus`` (default ``0.0`` --
        the production ``he_capture_velocity="at_rest"`` arm). A ``thermal`` He velocity
        is a Tier-3 rule-2 arm resolved by the caller, not here. NB a *scalar* ``u_he``
        broadcasts per-component, i.e. an He velocity vector ``(u, u, u)`` -- inert at
        ``u = 0`` but wrong for a thermal speed; the ``thermal`` arm must pass a
        per-component He velocity *vector*, not a scalar speed.

    Returns
    -------
    CaptureResult
        ``v_plus = (m*v- + m_He*u_He)/(m + m_He)``, ``m_plus_amu = m + m_He``, and
        ``dE_mass_transfer = +0.5*(m*m_He)/(m + m_He)*|v- - u_He|^2`` (``>= 0``).

    Raises
    ------
    ValueError
        If ``m_he_amu <= 0`` or ``m_minus_amu <= 0`` (or either non-finite), or if
        ``v_minus`` / ``u_he`` contain non-finite values.
    """
    _check_masses_gain(m_minus_amu, m_he_amu)
    v = np.asarray(v_minus, dtype=float)
    u = np.asarray(u_he, dtype=float)
    if not (np.all(np.isfinite(v)) and np.all(np.isfinite(u))):
        raise ValueError(f"v_minus and u_he must be finite; got {v_minus!r}, {u_he!r}.")

    m_plus = m_minus_amu + m_he_amu
    v_rel = v - u
    v_plus = (m_minus_amu * v + m_he_amu * u) / m_plus

    speed_sq = float(v_rel @ v_rel) if v_rel.ndim else float(v_rel) ** 2
    # Exact KE absorbed by the merge (reduced-mass form); booked positive (the coeff
    # is an unsigned magnitude -- sign at the call site) into the mass-transfer channel.
    dE_mass_transfer = _reduced_mass_defect_coeff(m_minus_amu, m_plus, m_he_amu) * speed_sq

    return CaptureResult(v_plus=v_plus, m_plus_amu=m_plus, dE_mass_transfer=dE_mass_transfer)


def capture_velocity_components(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    m_minus_amu,
    *,
    m_he_amu: float = MASS_HE_AMU,
    u_he: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-ion vectorized He capture over ``(M,)`` velocity component arrays.

    The ensemble-facing form of :func:`capture`: the same momentum-conserving reset
    applied independently per ion. Unlike :func:`cold_shed_velocity_components` (which
    takes a *uniform* scalar complex mass, since the anchored Tier-1a schedule sheds
    from every ion at one instant), pickup fires per ion independently, so
    ``m_minus_amu`` may be a **per-ion array**; scalar and array both broadcast.
    Reuses :func:`_reduced_mass_defect_coeff`, so the physics is identical to the
    scalar :func:`capture` (the single-ion oracle in the tests).

    Parameters
    ----------
    vx, vy, vz : np.ndarray, shape (M,)
        Pre-capture velocity components [A/ps].
    m_minus_amu : float or np.ndarray
        Pre-capture complex mass(es) [amu]; finite and ``> 0`` (scalar or ``(M,)``).
    m_he_amu : float, optional
        Captured He mass [amu] (default :data:`MASS_HE_AMU`).
    u_he : float, optional
        Incoming He velocity [A/ps] (default ``0.0``, the ``at_rest`` arm). NB a
        scalar here means the He vector ``(u, u, u)`` for *every* ion -- inert at
        ``u = 0``, but the Tier-3 ``thermal`` arm needs per-ion He velocity
        *components* (e.g. ``ux, uy, uz`` arrays); extend this signature then rather
        than threading a scalar speed.

    Returns
    -------
    (vx_plus, vy_plus, vz_plus, m_plus_amu, dE_mass_transfer) : tuple
        Merged components (each ``(M,)``), the post-capture mass ``m + m_He``
        (scalar or ``(M,)`` matching ``m_minus_amu``), and the per-ion defect ``(M,)``
        ``+0.5*(m*m_He)/(m + m_He)*|v_i - u_He|^2`` [amu*A^2/ps^2], ``>= 0``.

    Raises
    ------
    ValueError
        If any mass is non-positive/non-finite or a velocity component is non-finite.
    """
    _check_masses_gain(m_minus_amu, m_he_amu)
    vxf = np.asarray(vx, dtype=float)
    vyf = np.asarray(vy, dtype=float)
    vzf = np.asarray(vz, dtype=float)
    u = np.asarray(u_he, dtype=float)
    if not (
        np.all(np.isfinite(vxf)) and np.all(np.isfinite(vyf))
        and np.all(np.isfinite(vzf)) and np.all(np.isfinite(u))
    ):
        raise ValueError("velocity components and u_he must be finite.")

    m_plus = np.asarray(m_minus_amu, dtype=float) + m_he_amu
    coeff = _reduced_mass_defect_coeff(np.asarray(m_minus_amu, dtype=float), m_plus, m_he_amu)
    vrx = vxf - u
    vry = vyf - u
    vrz = vzf - u
    vx_plus = (np.asarray(m_minus_amu, dtype=float) * vxf + m_he_amu * u) / m_plus
    vy_plus = (np.asarray(m_minus_amu, dtype=float) * vyf + m_he_amu * u) / m_plus
    vz_plus = (np.asarray(m_minus_amu, dtype=float) * vzf + m_he_amu * u) / m_plus
    speed_sq = vrx ** 2 + vry ** 2 + vrz ** 2
    dE_mass_transfer = coeff * speed_sq
    return vx_plus, vy_plus, vz_plus, m_plus, dE_mass_transfer


def continuous_velocity_shed_components(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    m_minus_amu,
    *,
    m_he_amu: float = MASS_HE_AMU,
    n_removed: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
    """Per-atom vectorized co-moving shed over velocity component arrays.

    Velocity components are copied unchanged. The post-shed mass is
    ``m - n_removed*m_He`` and the per-atom ledger increment is
    ``+0.5*n_removed*m_He*|v_i|^2`` [amu*A^2/ps^2] — the exact lab-frame KE
    the co-moving He carries away (no reduced-mass correction: the He leaves
    at exactly ``v_i``).

    ``m_minus_amu`` may be a **uniform scalar** (the Tier-1a anchored schedule)
    or a **per-ion array** (the Tier-2 generative evaporation channel under the
    ``co_moving`` OQ-J arm — fires diverge, so masses differ across ions), the
    same contract as :func:`cold_shed_velocity_components`. The scalar path is
    byte-identical to the pre-generalization behaviour (the arithmetic
    broadcasts; only the mass check widened).
    """
    _check_masses_shed(m_minus_amu, m_he_amu)
    n = _check_removed_count(n_removed, m_minus_amu, m_he_amu)
    vxf = np.asarray(vx, dtype=float)
    vyf = np.asarray(vy, dtype=float)
    vzf = np.asarray(vz, dtype=float)
    if not (np.all(np.isfinite(vxf)) and np.all(np.isfinite(vyf)) and np.all(np.isfinite(vzf))):
        raise ValueError("velocity components must be finite.")

    m_plus = m_minus_amu - n * m_he_amu
    speed_sq = vxf ** 2 + vyf ** 2 + vzf ** 2
    dE_mass_transfer = 0.5 * n * m_he_amu * speed_sq
    return vxf.copy(), vyf.copy(), vzf.copy(), m_plus, dE_mass_transfer


def apply_shed(
    v_minus,
    m_minus_amu: float,
    *,
    mode: ShedMode,
    m_he_amu: float = MASS_HE_AMU,
    n_removed: int = 1,
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
        return continuous_velocity_shed(
            v_minus,
            m_minus_amu,
            m_he_amu=m_he_amu,
            n_removed=n_removed,
        )
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


def _check_masses_gain(m_minus_amu, m_he_amu: float) -> None:
    """Fail loudly unless ``m_he_amu > 0`` and every ``m_minus_amu > 0`` (a gain).

    The capture counterpart of :func:`_check_masses`. A capture *adds* mass
    (``m+ = m + m_He``), so the shed precondition ``m_minus - m_He > 0`` does not
    apply -- only positivity and finiteness are required (CLAUDE.md principle 4: fail
    loud for the right reason). ``m_minus_amu`` may be a scalar or a per-ion array.
    """
    m = np.asarray(m_minus_amu, dtype=float)
    if not (np.all(np.isfinite(m)) and np.isfinite(m_he_amu)):
        raise ValueError(
            f"masses must be finite; got m_minus_amu={m_minus_amu!r}, "
            f"m_he_amu={m_he_amu!r}."
        )
    if not (m_he_amu > 0.0):
        raise ValueError(f"m_he_amu must be > 0; got {m_he_amu!r}.")
    if not np.all(m > 0.0):
        raise ValueError(
            f"pre-capture complex mass must be > 0; got m_minus_amu={m_minus_amu!r}."
        )


def _check_masses_shed(m_minus_amu, m_he_amu: float) -> None:
    """Fail loudly unless ``m_he_amu > 0`` and every ``m_minus_amu - m_He > 0``.

    The array-aware form of :func:`_check_masses` (which is scalar-only). A cold shed
    *removes* mass, so the post-shed mass ``m+ = m - m_He`` must stay strictly positive
    for **every** ion (CLAUDE.md principle 4). ``m_minus_amu`` may be a scalar or a
    per-ion array; the scalar path reduces to exactly the :func:`_check_masses` checks,
    so :func:`cold_shed_velocity_components` stays byte-identical for uniform scalar mass
    (the Tier-1a contract) while also accepting the per-ion masses the Tier-2 generative
    evaporation channel produces (fires diverge, so masses differ across ions).
    """
    m = np.asarray(m_minus_amu, dtype=float)
    if not (np.all(np.isfinite(m)) and np.isfinite(m_he_amu)):
        raise ValueError(
            f"masses must be finite; got m_minus_amu={m_minus_amu!r}, "
            f"m_he_amu={m_he_amu!r}."
        )
    if not (m_he_amu > 0.0):
        raise ValueError(f"m_he_amu must be > 0; got {m_he_amu!r}.")
    if not np.all(m - m_he_amu > 0.0):
        raise ValueError(
            f"every pre-shed mass must exceed the He mass so m+ = m - m_He > 0; "
            f"got m_minus_amu={m_minus_amu!r}, m_he_amu={m_he_amu!r}."
        )


def _check_removed_count(n_removed: int, m_minus_amu, m_he_amu: float) -> int:
    """Return validated integer removed-He count.

    ``m_minus_amu`` may be a scalar or a per-ion array (the OQ-J co-moving
    generalization); every element must survive the removal. The scalar path
    is behaviour-identical to the original scalar-only check.
    """
    if isinstance(n_removed, (bool, np.bool_)) or not isinstance(n_removed, (int, np.integer)):
        raise ValueError(f"n_removed must be an integer; got {n_removed!r}.")
    n = int(n_removed)
    if n <= 0:
        raise ValueError(f"n_removed must be > 0; got {n_removed!r}.")
    if not np.all(np.asarray(m_minus_amu, dtype=float) - n * m_he_amu > 0.0):
        raise ValueError(
            f"n_removed={n} removes too much mass from m_minus_amu={m_minus_amu!r}."
        )
    return n
