r"""Poisson He-pickup channel for the I+He_n complex (Tier-2 Phase-B Slice P).

The **gain** channel of the biphasic mass mechanism: one independent Bernoulli draw
per ion per step and, on a fire, the momentum-conserving He **capture** reset plus the
S1 heat deposited into the internal reservoir. This module is a *stateful-per-event but
storage-free primitive* -- it takes the current ``(n, v, m)`` plus an **injected** RNG
and returns the post-event ``(n', m', v', dE_int, dE_mass_transfer, fired)``. It performs
**no** integration, no ``E_int`` state persistence, no schema I/O, and no 5-term closure
(those are the Phase-C schema + generative driver).

Encoded form (TIER2_PHASE_B_IMPLEMENTATION_PLAN.md §2.2 (P); MASS §4/§5)
------------------------------------------------------------------------
    P_attach(dt)   = 1 - exp(-lambda_attach * dt)
    lambda_attach  = lambda_0 * (rho_He/rho_bulk) * (1 - n/n*)_+^p     (langmuir cap)
    on fire:  n -> n+1,  m -> m + m_He,
              v+ = (m*v- + m_He*u_He)/(m + m_He),
              dE_int = +f_ret * D_0(n+1)                               (S1 heat, via U)

The capture reset (velocity/mass/reduced-mass KE defect) is delegated to
``mass_jump.capture`` / ``mass_jump.capture_velocity_components`` so the exact
reduced-mass form lives in one place (CLAUDE.md rule 1). The S1 heat is delegated to
``internal_energy_budget.dE_int_pickup_eV`` (``n`` = pre-pickup occupancy). The capture
KE defect (``dE_mass_transfer``, into ``E_mass_transfer``) and the S1 heat (into
``E_int``) are **two distinct injections** -- no double count; the closure combining them
is Phase C.

Units
-----
``lambda`` / ``lambda_0`` in ps^-1; ``dt`` in ps (``[lambda*dt] = 1``); ``rho_ratio`` /
occupancy factor dimensionless; ``v`` in A/ps; ``m`` in amu; ``dE_int`` in eV;
``dE_mass_transfer`` in amu*A^2/ps^2 (the mechanical convention of ``mass_jump``). The
module is **mass-agnostic in the drag sense** -- it never evaluates the friction
coefficient ``gamma``; ``m`` enters only as the captured/merged quantity in the reset.

RNG contract
------------
Every stochastic entry point takes an **injected** ``numpy.random.Generator`` (no
module-level global RNG, no implicit seeding). :func:`pickup_step` draws one scalar
uniform ``rng.random()``; :func:`pickup_step_components` draws one vectorized
``rng.random(size=n_ions)``. A fire is ``draw < P_attach``. The shed-then-pickup channel
order across the two channels is documented for one-event-per-step bookkeeping (A13) but
is *frozen* only at the Phase-C schema/driver, not here.

Rule-2 arms
-----------
Only ``pickup_rate_form="density_only"`` and ``he_capture_velocity="at_rest"`` are built.
``sweeping`` / ``dwell_time`` (rate form) and ``thermal`` (capture velocity) are valid
config enum members (a carrying config round-trips) but **unbuilt**: they raise
``NotImplementedError`` lazily here, at point-of-use, never at config-load.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import MASS_HE_AMU, N_STAR
from .internal_energy_budget import dE_int_pickup_eV
from .mass_jump import capture, capture_velocity_components

# Built vs. declared-but-unbuilt (rule-2) arms.
_BUILT_PICKUP_RATE_FORM = "density_only"
_BUILT_HE_CAPTURE_VELOCITY = "at_rest"
_KNOWN_OCCUPANCY_CAPS = ("langmuir", "none")


@dataclass(frozen=True)
class PickupResult:
    """The outcome of one (possibly null) pickup draw applied to a single ion.

    Symmetric to ``mass_jump.ShedResult`` / ``mass_jump.CaptureResult`` (a named record,
    not a bare tuple).

    Attributes
    ----------
    n_plus : int
        Post-step He occupancy: ``n + 1`` on a fire, ``n`` otherwise.
    m_plus_amu : float
        Post-step complex mass [amu]: ``m + m_He`` on a fire, ``m`` otherwise.
    v_plus : np.ndarray
        Post-step velocity [A/ps]: the capture centre-of-mass velocity on a fire, an
        unchanged copy otherwise.
    dE_int_eV : float
        The S1 internal-energy deposit [eV]: ``+f_ret*D_0(n+1)`` on a fire, ``0`` otherwise.
    dE_mass_transfer : float
        The capture KE defect [amu*A^2/ps^2] (``>= 0``) on a fire, ``0`` otherwise. Books
        into ``E_mass_transfer`` at Phase C.
    fired : bool
        Whether the Bernoulli draw fired this step.
    """

    n_plus: int
    m_plus_amu: float
    v_plus: np.ndarray
    dE_int_eV: float
    dE_mass_transfer: float
    fired: bool


def _resolve_u_he(he_capture_velocity: str) -> float:
    """Resolve the incoming-He velocity for the capture reset [A/ps].

    ``at_rest`` -> ``0.0`` (production). ``thermal`` is a Tier-3 rule-2 arm: it is a
    valid config member (config-load accepts it) but unbuilt, so it raises
    ``NotImplementedError`` here at point-of-use.
    """
    if he_capture_velocity == _BUILT_HE_CAPTURE_VELOCITY:
        return 0.0
    if he_capture_velocity == "thermal":
        raise NotImplementedError(
            "he_capture_velocity='thermal' is a declared-but-unbuilt rule-2 arm "
            "(Tier 3); only 'at_rest' (u_He=0) is built in Phase B."
        )
    raise ValueError(
        f"unknown he_capture_velocity {he_capture_velocity!r}; "
        f"expected 'at_rest' or 'thermal'."
    )


def lambda_attach(
    rho_ratio,
    n,
    *,
    lambda0: float,
    n_star: int = N_STAR,
    p: float = 1.0,
    cap: str = "langmuir",
    pickup_rate_form: str = "density_only",
):
    r"""Per-ion attachment rate ``lambda_0 * (rho_He/rho_bulk) * (1 - n/n*)_+^p`` [ps^-1].

    Parameters
    ----------
    rho_ratio : array_like
        The ``rho_He/rho_bulk`` density gate in ``[0, 1]`` (from ``helium_density``).
    n : array_like
        Current He occupancy (pre-pickup); dimensionless count.
    lambda0 : float
        Pickup rate coefficient ``lambda_0`` [ps^-1] (Sourced+Bounded prior).
    n_star : int, optional
        First-shell capacity ``n*`` (default :data:`N_STAR` = 21).
    p : float, optional
        Langmuir occupancy exponent (default ``1.0``; held fixed, not ``p=kappa``).
    cap : {"langmuir", "none"}, optional
        ``langmuir`` applies the saturation factor ``(1 - n/n*)_+^p``; ``none`` recovers
        the density-only limit (factor ``1``, cap inert).
    pickup_rate_form : {"density_only", ...}, optional
        Only ``density_only`` is built; ``sweeping`` / ``dwell_time`` raise
        ``NotImplementedError`` (rule-2 arms).

    Returns
    -------
    float or np.ndarray
        The rate, scalar-in -> float / array-in -> ndarray. ``>= 0`` for
        ``lambda_0, rho_ratio >= 0``; ``-> 0`` as ``n -> n*`` under the Langmuir cap.

    Raises
    ------
    NotImplementedError
        If ``pickup_rate_form`` is a valid-but-unbuilt rule-2 arm.
    ValueError
        If ``cap`` is not a recognised occupancy-cap selector.
    """
    if pickup_rate_form != _BUILT_PICKUP_RATE_FORM:
        raise NotImplementedError(
            f"pickup_rate_form={pickup_rate_form!r} is a declared-but-unbuilt rule-2 "
            f"arm; only 'density_only' is built in Phase B (MASS §5)."
        )
    if cap not in _KNOWN_OCCUPANCY_CAPS:
        raise ValueError(
            f"unknown pickup occupancy cap {cap!r}; expected one of {_KNOWN_OCCUPANCY_CAPS}."
        )

    rho = np.asarray(rho_ratio, dtype=float)
    if cap == "langmuir":
        # Fail loud on the two knobs that make the Langmuir factor produce a silent
        # inf/nan rate (CLAUDE.md principle 4): p < 0 gives 0**neg = inf at n >= n*,
        # and n_star <= 0 gives 0/0 = nan at n = 0. Both would fire every step.
        if p < 0.0:
            raise ValueError(
                f"pickup_occupancy_exponent p must be >= 0 (a negative exponent "
                f"diverges the Langmuir factor to inf at n >= n*); got {p!r}"
            )
        if not (n_star > 0):
            raise ValueError(
                f"n_star must be > 0 (a non-positive shell capacity makes the "
                f"occupancy factor nan/inf); got {n_star!r}"
            )
        occupancy = np.maximum(0.0, 1.0 - np.asarray(n, dtype=float) / n_star) ** p
    else:  # "none": density-only, cap inert (p / n_star unused)
        occupancy = 1.0
    out = lambda0 * rho * occupancy
    if np.ndim(rho_ratio) == 0 and np.ndim(n) == 0:
        return float(out)
    return out


def attach_probability(lambda_attach_per_ps, dt_ps: float):
    """Per-step Bernoulli attach probability ``P_attach = 1 - exp(-lambda*dt)``.

    The Poisson-to-Bernoulli reduction of one step. ``-> lambda*dt`` for small
    ``lambda*dt`` (the rare-event limit); non-decreasing in ``lambda`` and in ``[0, 1]``
    (it saturates to exactly ``1.0`` in float when ``lambda*dt`` is large enough that
    ``exp`` underflows -- the fire-every-step limit). Scalar-in -> float, array-in ->
    ndarray.
    """
    lam = np.asarray(lambda_attach_per_ps, dtype=float)
    out = 1.0 - np.exp(-lam * dt_ps)
    if np.ndim(lambda_attach_per_ps) == 0:
        return float(out)
    return out


def pickup_step(
    *,
    rng,
    n: int,
    v,
    m_amu: float,
    rho_ratio: float,
    lambda0: float,
    f_ret: float,
    picture: str = "statistical_mixture",
    kappa: float,
    dt_ps: float,
    n_star: int = N_STAR,
    p: float = 1.0,
    cap: str = "langmuir",
    pickup_rate_form: str = "density_only",
    he_capture_velocity: str = "at_rest",
    m_he_amu: float = MASS_HE_AMU,
) -> PickupResult:
    """Draw one Bernoulli pickup for a single ion; on fire, capture He + deposit S1 heat.

    The single-ion oracle for :func:`pickup_step_components`. Composes
    :func:`mass_jump.capture` (velocity/mass/KE-defect reset) and
    :func:`internal_energy_budget.dE_int_pickup_eV` (S1 heat, ``n`` = pre-pickup). Inputs
    are loose scalar kwargs -- the persistent ``(n, E_int, v, m)`` carrier is the Phase-C
    driver's ``IonStepState``, not a Phase-B container.

    Parameters
    ----------
    rng : numpy.random.Generator
        Injected RNG; one scalar ``rng.random()`` is drawn.
    n : int
        Current He occupancy (pre-pickup).
    v : array_like
        Current velocity vector [A/ps].
    m_amu : float
        Current complex mass [amu].
    rho_ratio : float
        ``rho_He/rho_bulk`` gate at the ion's depth.
    lambda0, f_ret, kappa, dt_ps : float
        Rate coefficient, S1 retained fraction, ladder steepness, integrator step.
    picture : str, optional
        Electronic picture threaded to ``dE_int_pickup_eV`` / the ladder.
    n_star, p, cap, pickup_rate_form, he_capture_velocity, m_he_amu
        Occupancy/cap/reset knobs (see :func:`lambda_attach`, :func:`_resolve_u_he`).

    Returns
    -------
    PickupResult
    """
    u_he = _resolve_u_he(he_capture_velocity)
    lam = lambda_attach(
        rho_ratio, n, lambda0=lambda0, n_star=n_star, p=p, cap=cap,
        pickup_rate_form=pickup_rate_form,
    )
    p_attach = attach_probability(lam, dt_ps)
    fired = bool(float(rng.random()) < p_attach)

    if not fired:
        return PickupResult(
            n_plus=int(n),
            m_plus_amu=float(m_amu),
            v_plus=np.asarray(v, dtype=float).copy(),
            dE_int_eV=0.0,
            dE_mass_transfer=0.0,
            fired=False,
        )

    cap_res = capture(v, m_amu, m_he_amu=m_he_amu, u_he=u_he)
    dE_int = dE_int_pickup_eV(n, f_ret=f_ret, picture=picture, kappa=kappa)
    return PickupResult(
        n_plus=int(n) + 1,
        m_plus_amu=cap_res.m_plus_amu,
        v_plus=cap_res.v_plus,
        dE_int_eV=float(dE_int),
        dE_mass_transfer=cap_res.dE_mass_transfer,
        fired=True,
    )


def pickup_step_components(
    *,
    rng,
    n,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    m_amu,
    rho_ratio,
    lambda0: float,
    f_ret: float,
    picture: str = "statistical_mixture",
    kappa: float,
    dt_ps: float,
    n_star: int = N_STAR,
    p: float = 1.0,
    cap: str = "langmuir",
    pickup_rate_form: str = "density_only",
    he_capture_velocity: str = "at_rest",
    m_he_amu: float = MASS_HE_AMU,
):
    """Vectorized per-ion pickup over ``(M,)`` ensemble arrays (one draw per step).

    The ensemble-facing form of :func:`pickup_step`: one vectorized Bernoulli draw
    ``rng.random(size=M)``, per-ion independence, and the same composition of
    :func:`mass_jump.capture_velocity_components` + :func:`internal_energy_budget.
    dE_int_pickup_eV`. Non-fired ions keep their velocity/mass/occupancy and receive zero
    energy increments.

    Parameters
    ----------
    rng : numpy.random.Generator
        Injected RNG; one ``rng.random(size=M)`` draw.
    n : np.ndarray, shape (M,)
        Per-ion He occupancy (pre-pickup); integer-valued.
    vx, vy, vz : np.ndarray, shape (M,)
        Per-ion velocity components [A/ps].
    m_amu : np.ndarray, shape (M,)
        Per-ion complex mass [amu].
    rho_ratio : np.ndarray, shape (M,)
        Per-ion ``rho_He/rho_bulk`` gate.
    (remaining kwargs)
        As :func:`pickup_step`.

    Returns
    -------
    (n_plus, m_plus_amu, vx_plus, vy_plus, vz_plus, dE_int_eV, dE_mass_transfer, fired)
        Each an ``(M,)`` array (``fired`` boolean, ``n_plus`` integer), mirroring the
        tuple discipline of ``mass_jump.cold_shed_velocity_components``.
    """
    u_he = _resolve_u_he(he_capture_velocity)
    n_arr = np.asarray(n)
    m_arr = np.asarray(m_amu, dtype=float)

    lam = lambda_attach(
        rho_ratio, n_arr, lambda0=lambda0, n_star=n_star, p=p, cap=cap,
        pickup_rate_form=pickup_rate_form,
    )
    p_attach = np.asarray(attach_probability(lam, dt_ps), dtype=float)
    draws = rng.random(size=p_attach.shape)
    fired = draws < p_attach

    # Would-be capture for every ion (cheap, vectorized); non-fired results discarded.
    vx_c, vy_c, vz_c, m_c, dE_mt_c = capture_velocity_components(
        vx, vy, vz, m_arr, m_he_amu=m_he_amu, u_he=u_he,
    )
    dE_int_c = np.asarray(
        dE_int_pickup_eV(n_arr, f_ret=f_ret, picture=picture, kappa=kappa), dtype=float
    )

    n_plus = np.where(fired, n_arr + 1, n_arr)
    m_plus = np.where(fired, m_c, m_arr)
    vx_plus = np.where(fired, vx_c, np.asarray(vx, dtype=float))
    vy_plus = np.where(fired, vy_c, np.asarray(vy, dtype=float))
    vz_plus = np.where(fired, vz_c, np.asarray(vz, dtype=float))
    dE_int = np.where(fired, dE_int_c, 0.0)
    dE_mt = np.where(fired, dE_mt_c, 0.0)
    return n_plus, m_plus, vx_plus, vy_plus, vz_plus, dE_int, dE_mt, fired
