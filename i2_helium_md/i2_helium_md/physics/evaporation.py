r"""Energy-gated, RRK-rate-limited He-evaporation channel for the I+He_n complex
(Tier-2 Phase-B Slice Q).

The **loss** channel of the biphasic mass mechanism: the parameter-free self-bound
gate, then the saturating RRK shed of the top rung. Like the pickup channel this is a
*stateful-per-event but storage-free primitive* -- it takes the current
``(E_int, n, v, m)`` plus an **injected** RNG and returns the post-event
``(n', m', v', dE_int, dE_mass_transfer, fired)``. It performs **no** integration, no
``E_int`` state persistence, no schema I/O, and no 5-term closure (those are the Phase-C
schema + generative driver).

Encoded form (TIER2_PHASE_B_IMPLEMENTATION_PLAN.md §2.2 (Q); MASS §4/§R9)
------------------------------------------------------------------------
Shedding happens in the band ``D_0(n) < E_int < Sigma(n)``:

* the **upper** bound is the parameter-free **self-bound gate** -- shedding is
  *suppressed* while ``E_int > Sigma(n) = sum_{i<=n} D_0(i)`` (the complex is net
  self-unbound; the droplet holds it together and an instantaneous gate would strip
  wholesale). Equivalently :func:`is_self_bound` (``E_int < Sigma(n)``) must hold.
* the **lower** bound is RRK-bracket positivity -- below the top rung
  (``E_int <= D_0(n)``) the rate is exactly zero (cannot afford to break the bond).

Within the band the top rung sheds with the saturating unimolecular rate::

    k(E_int, n) = nu * (1 - D_0(n)/E_int)^(s-1),   s = effective_dof(n),   n >= 2
    P_shed(dt)  = 1 - exp(-k * dt)

bounded to ``k in [0, nu)`` (no gate-open avalanche; at most ~nu*dt per step). The
**n = 1 diatomic is the degenerate boundary** (the band ``D_0(1) < E_int < Sigma(1) =
D_0(1)`` collapses; ``s = 3n-3 = 0``), so it is modelled as **direct dissociation**
``k = nu`` under the *inverted* gate ``E_int > D_0(1)`` -- the last atom leaves once the
complex is hotter than its single remaining bond (MASS A11; the boundary of validity of
the statistical picture). The fire path therefore keys on :func:`rrk_rate` (which encodes
both gates); :func:`is_self_bound` / :func:`gate_margin_eV` are the n>=2 diagnostics.

On a fire: ``n -> n-1``, ``m -> m - m_He``, the cold-shed velocity reset (delegated to the
*real* :func:`mass_jump.cold_shed`, He leaving cold), and the **K1** drain ``dE_int =
-D_0(n)`` (delegated to :func:`internal_energy_budget.dE_int_shed_eV`, ``n`` = pre-shed).

Units
-----
``E_int`` / ``D_0`` / ``Sigma`` in eV; ``nu`` / ``k`` in ps^-1; ``dt`` in ps
(``[k*dt] = 1``); ``v`` in A/ps; ``m`` in amu; ``dE_int`` in eV; ``dE_mass_transfer`` in
amu*A^2/ps^2 (the mechanical convention of ``mass_jump``). The module is **mass-agnostic
in the drag sense** -- it never evaluates the friction coefficient ``gamma``; ``m`` enters
only as the shed quantity in the reset.

RNG contract
------------
Every stochastic entry point takes an **injected** ``numpy.random.Generator``.
:func:`evaporation_step` draws one scalar ``rng.random()`` **unconditionally** every step
(``k = 0`` under suppression -> ``P_shed = 0`` -> the draw structurally never fires but is
still consumed), and :func:`evaporation_step_components` draws one vectorized
``rng.random(size=n_ions)``. This keeps the scalar and components forms' RNG consumption
identical, so the scalar step stays the exact ion-by-ion oracle for the components form
(and matches the pickup channel; the shed-then-pickup draw order is *frozen* only at the
Phase-C schema/driver, not here).

Diagnostic override
-------------------
``gate_onset_eV`` replaces the computed ``Sigma(n)`` self-bound threshold with a fixed
value (the n>=2 gate only; n=1 keeps its ``D_0(1)`` direct onset). This is the
falsification/sensitivity lever behind the loud ``allow_gate_onset_override`` config
provenance guard (R10 diagnostic-lever, not a knob); ``None`` is the parameter-free
production path.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import MASS_HE_AMU
from .dissociation_ladder import d0_of_n, ladder_cumsum
from .internal_energy_budget import dE_int_shed_eV
from .mass_jump import cold_shed, cold_shed_velocity_components


@dataclass(frozen=True)
class EvaporationResult:
    """The outcome of one (possibly null) evaporation draw applied to a single ion.

    Symmetric to ``pickup.PickupResult`` / ``mass_jump.ShedResult`` (a named record, not
    a bare tuple). The ``n_plus`` / ``m_plus_amu`` "post-event" names are kept though
    evaporation is a **loss** (``n_plus = n - 1``, ``m_plus_amu = m - m_He`` on a fire),
    so one post-event record shape serves both channels.

    Attributes
    ----------
    n_plus : int
        Post-step He occupancy: ``n - 1`` on a fire, ``n`` otherwise.
    m_plus_amu : float
        Post-step complex mass [amu]: ``m - m_He`` on a fire, ``m`` otherwise.
    v_plus : np.ndarray
        Post-step velocity [A/ps]: the cold-shed reset on a fire, an unchanged copy
        otherwise.
    dE_int_eV : float
        The K1 internal-energy drain [eV]: ``-D_0(n) < 0`` on a fire, ``0`` otherwise.
    dE_mass_transfer : float
        The cold-shed KE defect [amu*A^2/ps^2] (``<= 0``) on a fire, ``0`` otherwise.
        Books into ``E_mass_transfer`` at Phase C.
    fired : bool
        Whether the Bernoulli draw fired this step.
    """

    n_plus: int
    m_plus_amu: float
    v_plus: np.ndarray
    dE_int_eV: float
    dE_mass_transfer: float
    fired: bool


def effective_dof(n):
    r"""Classical-RRK effective mode count ``s`` for the I+He_n complex (``n >= 2``).

    ``s = 3N - 5 = 4`` for ``n = 2`` (the linear 3-atom complex I+He_2 = ion + 2 He;
    CALIBRATION_MAP row 10 "n=2 linear +1") and ``s = 3n - 3`` for ``n >= 3`` (the
    nonlinear full complex, MASS A11). ``s`` is always ``>= 4`` here, so the RRK bracket
    exponent ``s - 1 >= 3`` never diverges.

    ``n = 1`` is **not** in the domain: the diatomic has a single vibrational mode, so the
    statistical count degenerates (``3n - 3 = 0``) and evaporation is handled by the
    ``k = nu`` direct-dissociation branch instead (see the module docstring). Asking for
    ``effective_dof`` of a ``n < 2`` rung is therefore a caller error and raises loudly
    (CLAUDE.md principle 4).

    Scalar-in -> int, array-in -> ndarray of the same shape.
    """
    n_arr = np.asarray(n)
    if np.any(n_arr < 2):
        raise ValueError(
            f"effective_dof is defined for n >= 2 (n=1 uses the direct k=nu branch, "
            f"n<=0 has no rung); got {n!r}"
        )
    out = np.where(n_arr == 2, 4, 3 * n_arr - 3)
    return int(out) if np.ndim(n) == 0 else out


def shed_probability(k_per_ps, dt_ps: float):
    """Per-step Bernoulli shed probability ``P_shed = 1 - exp(-k*dt)``.

    The Poisson-to-Bernoulli reduction of one step (the shed-channel counterpart of
    ``pickup.attach_probability``; same closed form). ``-> k*dt`` for small ``k*dt`` (the
    rare-event limit), non-decreasing in ``k``, in ``[0, 1]``. Scalar-in -> float,
    array-in -> ndarray.
    """
    k = np.asarray(k_per_ps, dtype=float)
    out = 1.0 - np.exp(-k * dt_ps)
    if np.ndim(k_per_ps) == 0:
        return float(out)
    return out


def _gate_threshold_eV(n, *, picture: str, kappa: float, gate_onset_eV):
    """Resolve the self-bound gate threshold [eV]: ``Sigma(n)`` or a fixed override.

    ``gate_onset_eV=None`` -> the parameter-free integrated ladder ``Sigma(n)`` (Derived,
    MASS R9). A set float -> that fixed threshold (the diagnostic override; broadcast).
    """
    if gate_onset_eV is not None:
        return float(gate_onset_eV)
    return ladder_cumsum(n, picture=picture, kappa=kappa)


def gate_margin_eV(E_int_eV, n, *, picture: str = "statistical_mixture", kappa: float,
                   gate_onset_eV=None):
    r"""Signed self-unbound margin ``G = E_int - Sigma(n)`` [eV].

    The mechanism's natural readout: ``G < 0`` (``E_int < Sigma(n)``) is the self-bound
    regime where the shell sheds; ``G > 0`` is net self-unbound (shedding suppressed).
    Pure/cheap. Scalar-in -> float, array-in -> ndarray. See :func:`is_self_bound`.
    """
    E = np.asarray(E_int_eV, dtype=float)
    thresh = np.asarray(
        _gate_threshold_eV(n, picture=picture, kappa=kappa, gate_onset_eV=gate_onset_eV),
        dtype=float,
    )
    out = E - thresh
    if np.ndim(E_int_eV) == 0 and np.ndim(n) == 0:
        return float(out)
    return out


def is_self_bound(E_int_eV, n, *, picture: str = "statistical_mixture", kappa: float,
                  gate_onset_eV=None):
    r"""The self-bound gate: ``True`` iff ``E_int < Sigma(n)`` (== ``gate_margin_eV < 0``).

    **True = self-bound = shell evaporation enabled** (MASS R9 "once self-bound ... the
    top rung sheds"). Shedding is *suppressed* while ``False`` (net self-unbound). This is
    the shell (``n >= 2``) diagnostic; the ``n = 1`` diatomic uses the inverted
    direct-dissociation gate, so :func:`is_self_bound` is **not** the fire predictor at
    ``n = 1`` (that is :func:`rrk_rate`). Scalar-in -> bool, array-in -> bool ndarray.
    """
    margin = gate_margin_eV(
        E_int_eV, n, picture=picture, kappa=kappa, gate_onset_eV=gate_onset_eV
    )
    out = np.asarray(margin) < 0.0
    if np.ndim(E_int_eV) == 0 and np.ndim(n) == 0:
        return bool(out)
    return out


def rrk_rate(E_int_eV, n, *, nu: float, picture: str = "statistical_mixture",
             kappa: float, evap_rrk_dof=None, gate_onset_eV=None):
    r"""Energy-gated RRK shed rate ``k(E_int, n)`` [ps^-1] -- the full fire-path rate.

    Encodes **all** gating so a downstream draw fires on ``rng.random() < P_shed`` alone:

    * ``n >= 2``: ``k = nu * (1 - D_0(n)/E_int)^(s-1)`` inside the self-bound band
      (``D_0(n) < E_int < Sigma(n)``), and ``k = 0`` outside it (below the top rung, or
      net self-unbound). Bounded to ``[0, nu)`` -- no avalanche.
    * ``n == 1``: direct dissociation ``k = nu`` once ``E_int > D_0(1)`` (the inverted
      gate), else ``0``.
    * ``n <= 0``: ``0`` (no rung to shed).

    Parameters
    ----------
    E_int_eV, n : array_like
        Complex internal energy [eV] and He occupancy (integer). Vectorised together.
    nu : float
        RRK prefactor [ps^-1] (Sourced/pinned; :data:`constants.NU_EVAP_PER_PS`).
    picture : str, optional
        Electronic picture threaded to the ladder.
    kappa : float
        Ladder cliff steepness.
    evap_rrk_dof : float, optional
        Override for ``s`` applied to the ``n >= 2`` bracket only (``n = 1`` stays
        direct). ``None`` -> per-``n`` :func:`effective_dof`. Fail-loud on ``s < 1`` (the
        divergent-rate regime; a module-level defense mirroring the config-load guard and
        the pickup ``p < 0`` guard).
    gate_onset_eV : float, optional
        Diagnostic self-bound threshold override (n>=2 gate only); ``None`` -> ``Sigma(n)``.

    Returns
    -------
    float or np.ndarray
        ``k`` in ps^-1. Scalar-in -> float, array-in -> ndarray.
    """
    E = np.asarray(E_int_eV, dtype=float)
    n_arr = np.asarray(n)
    d0 = np.asarray(d0_of_n(n_arr, picture=picture, kappa=kappa), dtype=float)
    thresh = np.asarray(
        _gate_threshold_eV(n_arr, picture=picture, kappa=kappa, gate_onset_eV=gate_onset_eV),
        dtype=float,
    )

    # Effective dof for the n>=2 bracket. The override applies to n>=2 only; the per-n
    # form is evaluated on a clamped n (>=2) so the raising effective_dof never sees the
    # n<2 lanes -- those are handled by the direct / suppressed branches below and their
    # bracket value is discarded by the final np.where.
    if evap_rrk_dof is None:
        s = np.asarray(effective_dof(np.maximum(n_arr, 2)), dtype=float)
    else:
        if not (evap_rrk_dof >= 1.0):
            raise ValueError(
                f"evap_rrk_dof must satisfy s >= 1 (s < 1 diverges the RRK rate as "
                f"E_int -> D_0); got {evap_rrk_dof!r}"
            )
        s = float(evap_rrk_dof)

    # RRK bracket, guarded so no divide-by-zero / 0**0 / 0**neg is ever formed: E_safe is
    # never 0 (division safe), base is exactly 0 wherever E_int <= D_0(n) (or E_int <= 0)
    # giving k = 0 below threshold for any s, and the base>0 mask means base_safe**(s-1) is
    # only raised for a strictly-positive base. The exponent s-1 is >= 0 unconditionally:
    # per-n effective_dof(max(n,2)) >= 4, and the override is guarded s >= 1 above.
    E_pos = E > 0.0
    E_safe = np.where(E_pos, E, 1.0)
    base = np.where(E_pos, np.maximum(0.0, 1.0 - d0 / E_safe), 0.0)
    base_safe = np.where(base > 0.0, base, 1.0)
    k_rrk = nu * np.where(base > 0.0, base_safe ** (s - 1.0), 0.0)

    # Self-bound gate (n>=2): shed only while E_int < threshold (== is_self_bound).
    k_ge2 = np.where(E < thresh, k_rrk, 0.0)
    # n=1 direct dissociation: k = nu once E_int > D_0(1) (the inverted gate).
    k_direct = np.where(E > d0, nu, 0.0)

    k = np.where(n_arr == 1, k_direct, np.where(n_arr >= 2, k_ge2, 0.0))
    if np.ndim(E_int_eV) == 0 and np.ndim(n) == 0:
        return float(k)
    return k


def evaporation_step(
    *,
    rng,
    E_int_eV: float,
    n: int,
    v,
    m_amu: float,
    nu: float,
    picture: str = "statistical_mixture",
    kappa: float,
    dt_ps: float,
    evap_rrk_dof=None,
    gate_onset_eV=None,
    m_he_amu: float = MASS_HE_AMU,
) -> EvaporationResult:
    """Draw one Bernoulli evaporation for a single ion; on fire, cold-shed + drain K1.

    The single-ion oracle for :func:`evaporation_step_components`. Composes
    :func:`rrk_rate` (all gating), :func:`mass_jump.cold_shed` (velocity/mass/KE-defect
    reset), and :func:`internal_energy_budget.dE_int_shed_eV` (K1 drain, ``n`` = pre-shed).
    Inputs are loose scalar kwargs -- the persistent ``(n, E_int, v, m)`` carrier is the
    Phase-C driver's ``IonStepState``, not a Phase-B container. One uniform ``rng.random()``
    is drawn **unconditionally** (see the RNG contract in the module docstring).

    Returns
    -------
    EvaporationResult
    """
    k = rrk_rate(
        E_int_eV, n, nu=nu, picture=picture, kappa=kappa,
        evap_rrk_dof=evap_rrk_dof, gate_onset_eV=gate_onset_eV,
    )
    p_shed = shed_probability(k, dt_ps)
    fired = bool(float(rng.random()) < p_shed)

    if not fired:
        return EvaporationResult(
            n_plus=int(n),
            m_plus_amu=float(m_amu),
            v_plus=np.asarray(v, dtype=float).copy(),
            dE_int_eV=0.0,
            dE_mass_transfer=0.0,
            fired=False,
        )

    shed = cold_shed(v, m_amu, m_he_amu=m_he_amu)
    dE_int = dE_int_shed_eV(n, picture=picture, kappa=kappa)
    return EvaporationResult(
        n_plus=int(n) - 1,
        m_plus_amu=shed.m_plus_amu,
        v_plus=shed.v_plus,
        dE_int_eV=float(dE_int),
        dE_mass_transfer=shed.dE_mass_transfer,
        fired=True,
    )


def evaporation_step_components(
    *,
    rng,
    E_int_eV,
    n,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    m_amu,
    nu: float,
    picture: str = "statistical_mixture",
    kappa: float,
    dt_ps: float,
    evap_rrk_dof=None,
    gate_onset_eV=None,
    m_he_amu: float = MASS_HE_AMU,
):
    """Vectorized per-ion evaporation over ``(M,)`` ensemble arrays (one draw per step).

    The ensemble-facing form of :func:`evaporation_step`: one vectorized Bernoulli draw
    ``rng.random(size=M)``, per-ion independence, and the same composition of
    :func:`mass_jump.cold_shed_velocity_components` (**per-ion mass**, since fires diverge)
    + :func:`internal_energy_budget.dE_int_shed_eV`. Non-fired / suppressed ions keep their
    velocity/mass/occupancy and receive zero energy increments.

    Parameters
    ----------
    rng : numpy.random.Generator
        Injected RNG; one ``rng.random(size=M)`` draw.
    E_int_eV : np.ndarray, shape (M,)
        Per-ion complex internal energy [eV].
    n : np.ndarray, shape (M,)
        Per-ion He occupancy (pre-shed); integer-valued.
    vx, vy, vz : np.ndarray, shape (M,)
        Per-ion velocity components [A/ps].
    m_amu : np.ndarray, shape (M,)
        Per-ion complex mass [amu].
    (remaining kwargs)
        As :func:`evaporation_step`.

    Returns
    -------
    (n_plus, m_plus_amu, vx_plus, vy_plus, vz_plus, dE_int_eV, dE_mass_transfer, fired)
        Each an ``(M,)`` array (``fired`` boolean, ``n_plus`` integer), mirroring the
        tuple discipline of ``mass_jump.cold_shed_velocity_components`` /
        ``pickup.pickup_step_components``.
    """
    n_arr = np.asarray(n)
    m_arr = np.asarray(m_amu, dtype=float)

    k = rrk_rate(
        E_int_eV, n_arr, nu=nu, picture=picture, kappa=kappa,
        evap_rrk_dof=evap_rrk_dof, gate_onset_eV=gate_onset_eV,
    )
    p_shed = np.asarray(shed_probability(k, dt_ps), dtype=float)
    draws = rng.random(size=p_shed.shape)
    fired = draws < p_shed

    # Would-be cold shed for every ion (per-ion mass; non-fired results discarded).
    vx_s, vy_s, vz_s, m_s, dE_mt_s = cold_shed_velocity_components(
        vx, vy, vz, m_arr, m_he_amu=m_he_amu,
    )
    dE_int_s = np.asarray(dE_int_shed_eV(n_arr, picture=picture, kappa=kappa), dtype=float)

    n_plus = np.where(fired, n_arr - 1, n_arr)
    m_plus = np.where(fired, m_s, m_arr)
    vx_plus = np.where(fired, vx_s, np.asarray(vx, dtype=float))
    vy_plus = np.where(fired, vy_s, np.asarray(vy, dtype=float))
    vz_plus = np.where(fired, vz_s, np.asarray(vz, dtype=float))
    dE_int = np.where(fired, dE_int_s, 0.0)
    dE_mt = np.where(fired, dE_mt_s, 0.0)
    return n_plus, m_plus, vx_plus, vy_plus, vz_plus, dE_int, dE_mt, fired
