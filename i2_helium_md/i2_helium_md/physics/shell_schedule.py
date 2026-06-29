"""He-shell schedule generator for the Tier-1a anchored kinematic validation (Slice S).

Pure, stateless kinematics: turns the three 9 A TDDFT shell anchors
``(t*, 21), (10, 19), (14, 14)`` into the deterministic, monotone, sheds-only
integer schedule ``n_bar(t)`` and its **7 shed events**
(``TIER1A_IMPLEMENTATION_PLAN.md`` Sec.4 / Sec.10). No physics state, no
integrator coupling -- this is the first, fully-independent Tier-1a build unit.
Everything downstream (the mass-jump reset, the variable-mass BAOAB wiring, the
four-term ledger) *consumes* this output and lives in later slices.

Schedule (the locked piecewise-linear shell count, Sec.4)
---------------------------------------------------------
With ``t*`` the swept schedule onset (``n = 21`` held for ``t <= t*``)::

    n_bar(t) = 21                          ,  t <= t*
             = 21 - 2*(t - t*)/(10 - t*)   ,  t* < t <= 10
             = 19 - 1.25*(t - 10)          ,  10 < t <= 14
             = 14                          ,  t > 14

The 10 ps slope kink (segment loss rates ``2/(10-t*)`` then ``1.25`` /ps -- a real
~5x acceleration) is kept as a kink, not smoothed. A shed event fires at the
**half-integer downward crossing** ``n_bar = n - 1/2`` (the S4 discretization);
the crossing fraction is exposed as ``crossing_fraction`` so the verdict-robustness
of the telescoping invariant to the tie-break can be exercised.

``n_bar(t)`` vs ``n(t)`` -- continuous curve vs integer count
-------------------------------------------------------------
``n_bar(t)`` (:meth:`ShellSchedule.n_bar`) is the **continuous anchored loss
curve** and is fractional by construction -- it is the interpolant of the TDDFT
waypoints whose only role is to fix *when* sheds fire (the half-integer crossings
above). It is **not** the physical shell count. The **physical, integer shell
count** is ``n(t)`` (:meth:`ShellSchedule.n_of_t`): a piecewise-constant staircase
that starts at 21 and steps down by exactly one at each of the 7 shed times. ``n(t)``
-- never ``n_bar(t)`` -- is what the complex mass uses,
``m(t) = MASS_I_ION_AMU + n(t)*MASS_HE_AMU``.

Masses (Sec.2, Sec.10)
----------------------
The complex mass at shell count ``n`` is ``m(n) = MASS_I_ION_AMU + n*MASS_HE_AMU``
(both [amu]; ``MASS_I_ION_AMU = 126.90`` is the Tier-1a shell-schedule iodine-ion
reference, intentionally distinct from the MD ``MASS_I_AMU = 127.0`` -- see
``constants.py``). Each cold shed removes one He, so the momentum-conserving
velocity kick (applied later by the mass-jump slice, not here) is the ratio
``kick = m(n) / m(n-1) = m / (m - m_He)``. This module reports the kick as a
*property of the schedule*; it never takes or applies a velocity. Telescoping over
all sheds gives the force-free speed-boost ceiling
``prod(kick) = m(t*)/m(end) = m(21)/m(14) = 1.1532``, timing- and tie-break-independent.

Units-in-names follow the package convention (``*_ps``, ``*_amu``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import MASS_HE_AMU, MASS_I_ION_AMU


# ---------------------------------------------------------------------------
# Anchors (9 A TDDFT loss curve, author-confirmed; TIER1A Sec.2 / Sec.10).
# t* (schedule onset) is swept per run and is *not* an anchor constant.
# ---------------------------------------------------------------------------
ANCHOR_N_START: int = 21          # shell count held for t <= t*
ANCHOR_T_MID_PS: float = 10.0     # ps -- first anchored loss waypoint
ANCHOR_N_MID: int = 19            # shell count at ANCHOR_T_MID_PS
ANCHOR_T_END_PS: float = 14.0     # ps -- terminal anchored waypoint
ANCHOR_N_END: int = 14            # shell count at/after ANCHOR_T_END_PS (flat tail)

# Number of sheds = ANCHOR_N_START - ANCHOR_N_END (21 -> 14 = 7).
NUM_SHED_EVENTS: int = ANCHOR_N_START - ANCHOR_N_END


def complex_mass_amu(n) -> np.ndarray | float:
    """I+He_n complex mass [amu] at integer shell count ``n``.

    ``m(n) = MASS_I_ION_AMU + n*MASS_HE_AMU``. Accepts a scalar or an array;
    returns the same shape (scalar in -> float out).
    """
    m = MASS_I_ION_AMU + np.asarray(n, dtype=float) * MASS_HE_AMU
    return float(m) if np.ndim(n) == 0 else m


@dataclass(frozen=True)
class ShedEvent:
    """One cold-shed event in the anchored schedule.

    Attributes
    ----------
    index : int
        1-based ordinal within the schedule (1..7).
    crossing : float
        The ``n_bar`` value at which the shed fires (``n_before - crossing_fraction``;
        the S4 half-integer rule gives ``n_before - 0.5``).
    time_ps : float
        Fire time [ps], the analytic ``n_bar(t) = crossing`` solution.
    n_before, n_after : int
        Integer shell count immediately before / after the shed (``n_after = n_before - 1``).
    mass_before_amu, mass_after_amu : float
        Complex mass [amu] before / after the shed (``complex_mass_amu``).
    kick_factor : float
        Momentum-conserving speed kick the later mass-jump slice will apply,
        ``mass_before_amu / mass_after_amu = m / (m - m_He)``. Reported here as a
        schedule property; no velocity is touched in this module.
    """

    index: int
    crossing: float
    time_ps: float
    n_before: int
    n_after: int
    mass_before_amu: float
    mass_after_amu: float
    kick_factor: float


@dataclass(frozen=True)
class ShellSchedule:
    """Anchored He-shell schedule: the ``n_bar(t)`` evaluator plus the ordered sheds.

    Built by :func:`build_shell_schedule`; treat as immutable. ``events`` is the
    time-ordered tuple of :class:`ShedEvent` (length :data:`NUM_SHED_EVENTS`).
    """

    t_star_ps: float
    crossing_fraction: float
    events: tuple[ShedEvent, ...]

    def n_bar(self, t_ps) -> np.ndarray | float:
        """Continuous piecewise-linear shell count ``n_bar(t)`` [dimensionless].

        Vectorized: accepts a scalar or an ndarray of times [ps] and returns the
        same shape (scalar in -> float out). Flat at 21 for ``t <= t*`` and at 14
        for ``t > 14``; the two interior segments carry the anchored loss slopes.
        """
        t = np.asarray(t_ps, dtype=float)
        if _is_onset_strip_schedule(self):
            event = self.events[0]
            n = np.where(t < event.time_ps, float(event.n_before), float(event.n_after))
            return float(n) if np.ndim(t_ps) == 0 else n

        ts = self.t_star_ps
        slope1 = (ANCHOR_N_START - ANCHOR_N_MID) / (ANCHOR_T_MID_PS - ts)
        slope2 = (ANCHOR_N_MID - ANCHOR_N_END) / (ANCHOR_T_END_PS - ANCHOR_T_MID_PS)
        n = np.select(
            [t <= ts, t <= ANCHOR_T_MID_PS, t <= ANCHOR_T_END_PS],
            [
                float(ANCHOR_N_START),
                ANCHOR_N_START - slope1 * (t - ts),
                ANCHOR_N_MID - slope2 * (t - ANCHOR_T_MID_PS),
            ],
            default=float(ANCHOR_N_END),
        )
        return float(n) if np.ndim(t_ps) == 0 else n

    def n_of_t(self, t_ps) -> np.ndarray | int:
        """Physical **integer** shell count ``n(t)`` [dimensionless].

        The piecewise-constant staircase the complex mass actually uses
        (``m(t) = MASS_I_ION_AMU + n(t)*MASS_HE_AMU``): starts at
        :data:`ANCHOR_N_START` and steps down by one at each shed fire time. A shed
        at time ``t_e`` is taken to have already occurred at ``t == t_e`` (the count
        is decremented from that instant). Vectorized: scalar in -> ``int`` out;
        ndarray in -> integer ndarray of the same shape.

        Distinct from :meth:`n_bar`: ``n_bar`` is the continuous anchored loss curve
        whose half-integer crossings *define* these steps; this is the integer count
        itself, and is the quantity downstream slices (mass jump, integrator) consume.
        """
        t = np.asarray(t_ps, dtype=float)
        if _is_onset_strip_schedule(self):
            event = self.events[0]
            n = np.where(t < event.time_ps, event.n_before, event.n_after)
            return int(n) if np.ndim(t_ps) == 0 else n.astype(int)

        fire_times = np.array([e.time_ps for e in self.events])
        shed_count = np.count_nonzero(t[..., None] >= fire_times, axis=-1)
        n = ANCHOR_N_START - shed_count
        return int(n) if np.ndim(t_ps) == 0 else n.astype(int)


def _crossing_time_ps(crossing: float, t_star_ps: float) -> float:
    """Analytic ``t`` at which ``n_bar(t) == crossing`` (downward), no root-finding.

    ``crossing`` lies in the open interval ``(ANCHOR_N_END, ANCHOR_N_START)`` and
    is resolved on segment 1 if ``crossing >= ANCHOR_N_MID`` else segment 2 (the two
    expressions agree at the ``crossing == ANCHOR_N_MID`` join, ``t = 10`` ps).
    """
    if crossing >= ANCHOR_N_MID:  # segment 1: 21 - slope1*(t - t*)
        slope1 = (ANCHOR_N_START - ANCHOR_N_MID) / (ANCHOR_T_MID_PS - t_star_ps)
        return t_star_ps + (ANCHOR_N_START - crossing) / slope1
    # segment 2: 19 - slope2*(t - 10)
    slope2 = (ANCHOR_N_MID - ANCHOR_N_END) / (ANCHOR_T_END_PS - ANCHOR_T_MID_PS)
    return ANCHOR_T_MID_PS + (ANCHOR_N_MID - crossing) / slope2


def _is_onset_strip_schedule(schedule: ShellSchedule) -> bool:
    """The one-event crossing_fraction=1 sentinel marks diagnostic onset-strip mode."""
    return len(schedule.events) == 1 and schedule.crossing_fraction == 1.0


def _validate_onset_strip_n_final(n_final) -> int:
    message = f"n_final must be an integer in [0, {ANCHOR_N_START - 1}], got {n_final!r}."
    if isinstance(n_final, (bool, np.bool_)) or np.ndim(n_final) != 0:
        raise ValueError(message)

    try:
        n_value = float(n_final)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(message) from exc

    if not np.isfinite(n_value) or not n_value.is_integer():
        raise ValueError(message)

    n_after = int(n_value)
    if not (0 <= n_after < ANCHOR_N_START):
        raise ValueError(message)
    return n_after


def build_shell_schedule(
    t_star_ps: float,
    *,
    crossing_fraction: float = 0.5,
) -> ShellSchedule:
    """Build the anchored He-shell schedule and its 7 shed events.

    Parameters
    ----------
    t_star_ps : float
        Schedule onset [ps]; ``n = 21`` is held for ``t <= t_star_ps``. Must lie in
        ``[0, 10)`` so segment 1 has positive width (the Tier-1a sweep uses
        ``{0.5, 5, 9}``).
    crossing_fraction : float, optional
        Fraction below the integer at which a shed fires; the S4 half-integer rule
        is ``0.5`` (default). Must lie in ``(0, 1)``. The shed *count*, the shell-count
        sequence, and the telescoping kick product are invariant to this choice --
        only the fire *times* move -- which is exactly the verdict-robustness the
        Tier-1a telescoping unit test exercises.

    Returns
    -------
    ShellSchedule
        Frozen schedule with the ``n_bar`` evaluator and the time-ordered events.

    Raises
    ------
    ValueError
        If ``t_star_ps`` is not in ``[0, 10)`` or ``crossing_fraction`` not in ``(0, 1)``;
        or (defensively) if the assembled schedule is not the expected monotone
        7-event sheds-only sequence.
    """
    if not (0.0 <= t_star_ps < ANCHOR_T_MID_PS):
        raise ValueError(
            f"t_star_ps must be in [0, {ANCHOR_T_MID_PS}) ps so segment 1 has "
            f"positive width; got {t_star_ps}."
        )
    if not (0.0 < crossing_fraction < 1.0):
        raise ValueError(
            f"crossing_fraction must be in (0, 1); got {crossing_fraction}."
        )

    events: list[ShedEvent] = []
    # Sheds in firing order: n_before runs 21, 20, ..., 15 (n_after = n_before - 1).
    for i, n_before in enumerate(range(ANCHOR_N_START, ANCHOR_N_END, -1), start=1):
        n_after = n_before - 1
        crossing = n_before - crossing_fraction
        time_ps = _crossing_time_ps(crossing, t_star_ps)
        mass_before = complex_mass_amu(n_before)
        mass_after = complex_mass_amu(n_after)
        events.append(
            ShedEvent(
                index=i,
                crossing=crossing,
                time_ps=time_ps,
                n_before=n_before,
                n_after=n_after,
                mass_before_amu=mass_before,
                mass_after_amu=mass_after,
                kick_factor=mass_before / mass_after,
            )
        )

    # Fail-loud structural guards (TIER1A Sec.4 acceptance).
    if len(events) != NUM_SHED_EVENTS:
        raise ValueError(
            f"expected {NUM_SHED_EVENTS} shed events, assembled {len(events)}."
        )
    times = [e.time_ps for e in events]
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError(f"shed times not strictly increasing: {times}.")
    if not all(t_star_ps < t <= ANCHOR_T_END_PS for t in times):
        raise ValueError(
            f"shed times must lie in (t*={t_star_ps}, {ANCHOR_T_END_PS}] ps: {times}."
        )

    return ShellSchedule(
        t_star_ps=t_star_ps,
        crossing_fraction=crossing_fraction,
        events=tuple(events),
    )


def build_onset_strip_schedule(
    *,
    t_strip_ps: float = 0.5,
    n_final: int,
) -> ShellSchedule:
    """Build a one-event onset-violent stripping stress schedule.

    This is a diagnostic stress schedule, not a TDDFT-anchored physical loss curve.
    It strips directly from n=21 to ``n_final`` at ``t_strip_ps``.
    """
    if not np.isfinite(t_strip_ps) or t_strip_ps < 0.0:
        raise ValueError(f"t_strip_ps must be finite and >= 0; got {t_strip_ps!r}.")
    n_after = _validate_onset_strip_n_final(n_final)
    mass_before = complex_mass_amu(ANCHOR_N_START)
    mass_after = complex_mass_amu(n_after)
    event = ShedEvent(
        index=1,
        crossing=float(n_after),
        time_ps=float(t_strip_ps),
        n_before=ANCHOR_N_START,
        n_after=n_after,
        mass_before_amu=mass_before,
        mass_after_amu=mass_after,
        kick_factor=mass_before / mass_after,
    )
    return ShellSchedule(
        t_star_ps=float(t_strip_ps),
        crossing_fraction=1.0,
        events=(event,),
    )
