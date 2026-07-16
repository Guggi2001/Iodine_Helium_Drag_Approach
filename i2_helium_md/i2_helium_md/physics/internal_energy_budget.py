r"""Per-channel ``E_int`` budget bookkeeping for the I+He_n complex (Tier-2 Phase-A
Slice U).

Pure, stateless energetics: the closed-form ``E_int`` heating/draining increments and
the post-crossing reconstruction identity. **No reservoir state lives here** -- the
``E_int`` *state* is the Phase-C generative driver; this module only supplies the
per-event arithmetic that driver consumes. No RNG, no integrator, no mass and no
velocity (Phase A is mass-agnostic; the friction/force-coefficient ``gamma`` is never
touched here). The reduced-mass capture defect named in MASS S1 is **booked by Slice P
at Phase B**, not computed here.

This module *consumes* Slice L's rungs (``d0_of_n`` / ``ladder_cumsum``) and, for the
reconstruction, Slice K's binding split (``e_bind_pair_eV`` / ``e_electrostriction_eV``);
everything stateful (channels P/Q, Phase B) and composed (driver G, Phase C) lives later
and consumes this module.

Encoded form (TIER2_PHASE_A_IMPLEMENTATION_PLAN.md §2.2 (U); MASS §6)
--------------------------------------------------------------------
    S2 onset:        E_int(0)    = f_int * E_avail^ion
    S1 pickup:       dE_int      = +f_ret * D_0(n+1)     (bath gets (1-f_ret)*D_0(n+1))
    K1 shed:         dE_int      = -D_0(n)
    reconstruction:  E_int       = E_solv.struct - E_bind^pair(N) - E_elec(N)   (post-t* only)
    f_int floor:     f_int^floor = Sigma(n*) / E_avail^ion

**Index convention.** A pickup grows the complex ``n -> n+1`` and forms the bond at rung
``n+1``, so it releases ``D_0(n+1)`` (``n`` = pre-pickup occupancy); a shed breaks the
current rung ``n`` and consumes ``D_0(n)`` (``n`` = pre-shed occupancy). The S1 split
``f_ret*D_0 + (1-f_ret)*D_0 = D_0`` closes exactly.

**Reconstruction is post-t* only (A9).** Recovering ``E_int`` from ``E_solv.struct`` via
the static binding split is valid only after the gate-crossing time ``t_x``;
:func:`reconstruct_e_int_eV` therefore *raises* loudly when ``post_crossing`` is False so
no downstream slice reconstructs ``E_int`` inside the early window.

**The f_int floor is advisory, not a constraint** (MASS S2, CALIBRATION row 14): a
sub-floor ``f_int`` is physically self-unbound but Tier 2 may probe it -- so
:func:`e_int_onset_eV` performs no floor check, and :func:`f_int_floor` is a diagnostic
helper only.

Units follow the package convention: all energies in eV (``*_eV``); ``f_int``, ``f_ret``,
``kappa`` dimensionless; ``n``/``N`` dimensionless integer occupancy. Vectorised
scalar-in -> float / array-in -> ndarray, mirroring ``dissociation_ladder.py`` and
``solvation_cooling.py``.
"""

from __future__ import annotations

import numpy as np

from .constants import N_STAR, S_ABS_EV
from .dissociation_ladder import d0_of_n, ladder_cumsum
from .solvation_cooling import e_bind_pair_eV, e_electrostriction_eV


def e_int_onset_eV(*, f_int: float, e_avail_eV: float):
    """S2 Coulomb-onset internal-energy deposit ``E_int(0) = f_int * E_avail`` [eV].

    Pure arithmetic. Performs **no** floor check -- the self-unbound floor
    ``f_int^floor = Sigma(n*)/E_avail`` is an *advisory* bound (MASS S2 / CALIBRATION
    row 14, "not a constraint"); a sub-floor ``f_int`` is accepted. Use
    :func:`f_int_floor` for the diagnostic floor. Scalar-in -> float.
    """
    out = np.asarray(f_int) * np.asarray(e_avail_eV)
    if np.ndim(f_int) == 0 and np.ndim(e_avail_eV) == 0:
        return float(out)
    return out


def dE_int_pickup_eV(n, *, f_ret: float, picture: str = "statistical_mixture",
                     kappa: float, ladder=None):
    """S1 pickup heating ``dE_int = +f_ret * D_0(n+1)`` [eV] (``n`` = pre-pickup).

    The retained fraction of the bond energy released when the complex grows
    ``n -> n+1``; the remainder ``(1-f_ret)*D_0(n+1)`` goes to the bath
    (:func:`pickup_bath_release_eV`). Strictly positive for ``f_ret > 0``. Consumes
    Slice L's ``d0_of_n``. Scalar-in -> float, array-in -> ndarray.
    ``ladder`` (Slice T2, §I.10): when set, ``D_0`` comes from the injected
    :class:`~i2_helium_md.physics.dissociation_ladder.TabulatedLadder` and
    ``picture``/``kappa`` are ignored; ``None`` (default) is the byte-inert
    Form-U path (same convention across this module).
    """
    d0_next = d0_of_n(np.asarray(n) + 1, picture=picture, kappa=kappa, ladder=ladder)
    out = f_ret * np.asarray(d0_next)
    return float(out) if np.ndim(n) == 0 else out


def pickup_bath_release_eV(n, *, f_ret: float, picture: str = "statistical_mixture",
                           kappa: float, ladder=None):
    """S1 pickup bath release ``(1-f_ret) * D_0(n+1)`` [eV] (``n`` = pre-pickup).

    The complement of :func:`dE_int_pickup_eV`; together they close the S1 split to
    the full rung ``D_0(n+1)`` exactly. Scalar-in -> float, array-in -> ndarray.
    ``ladder``: injected table overrides Form-U (``picture``/``kappa`` ignored
    when set; see :func:`dE_int_pickup_eV`).
    """
    d0_next = d0_of_n(np.asarray(n) + 1, picture=picture, kappa=kappa, ladder=ladder)
    out = (1.0 - f_ret) * np.asarray(d0_next)
    return float(out) if np.ndim(n) == 0 else out


def dE_int_shed_eV(n, *, picture: str = "statistical_mixture", kappa: float,
                   ladder=None):
    """K1 evaporation drain ``dE_int = -D_0(n)`` [eV] (``n`` = pre-shed occupancy).

    One rung is spent breaking the bond as the complex sheds ``n -> n-1``. Strictly
    negative. Consumes Slice L's ``d0_of_n``. Scalar-in -> float, array-in -> ndarray.
    ``ladder``: injected table overrides Form-U (``picture``/``kappa`` ignored
    when set; see :func:`dE_int_pickup_eV`).
    """
    d0 = d0_of_n(n, picture=picture, kappa=kappa, ladder=ladder)
    return -d0 if np.ndim(n) == 0 else -np.asarray(d0)


def reconstruct_e_int_eV(E_solv_struct_eV, N, *, picture: str = "statistical_mixture",
                         kappa: float, s_abs_eV: float = S_ABS_EV,
                         post_crossing: bool, ladder=None):
    """A9 reconstruction ``E_int = E_solv.struct - E_bind^pair(N) - E_elec(N)`` [eV].

    Recovers ``E_int`` from the cooled master variable ``E_solv.struct`` using the
    *same* binding split Slice K cools toward, so ``E_int^eq = 0`` at equilibrium with
    no spurious offset (R12 eliminated; MASS line 881). Consumes Slice L (rungs, via K)
    and Slice K (``e_bind_pair_eV`` + ``e_electrostriction_eV``).

    **Valid only after the gate-crossing time ``t_x`` (A9).** ``post_crossing`` is a
    required keyword: a False (or omitted) value **raises** ``ValueError`` -- the static
    reconstruction errs in the early window, so pre-crossing use is refused loudly
    (CLAUDE.md principle 4). Scalar-in -> float, array-in -> ndarray.
    ``ladder``: injected table sets the binding split and ``picture``/``kappa``
    are ignored (see :func:`dE_int_pickup_eV`).
    """
    if not post_crossing:
        raise ValueError(
            "reconstruct_e_int_eV is valid only after the gate-crossing time t_x "
            "(A9); the static binding split errs in the early window. Pass "
            "post_crossing=True only when the reconstruction is known to be post-t_x."
        )
    e_bind = e_bind_pair_eV(N, picture=picture, kappa=kappa, ladder=ladder)
    e_elec = e_electrostriction_eV(N, picture=picture, kappa=kappa, s_abs_eV=s_abs_eV,
                                   ladder=ladder)
    out = np.asarray(E_solv_struct_eV) - np.asarray(e_bind) - np.asarray(e_elec)
    if np.ndim(E_solv_struct_eV) == 0 and np.ndim(N) == 0:
        return float(out)
    return out


def f_int_floor(*, e_avail_eV: float, picture: str = "statistical_mixture",
                kappa: float, ladder=None):
    """Self-unbound floor ``f_int^floor = Sigma(n*) / E_avail`` [dimensionless].

    The S2 lower bound below which ``E_int(0)`` cannot keep the full first shell
    net-bound (MASS S2). Diagnostic only -- it is *not* enforced by
    :func:`e_int_onset_eV` (CALIBRATION row 14: advisory). Scenario-keyed through
    ``e_avail_eV`` (0.80 eV validation / 2.70 eV production) and picture-set /
    near-kappa-independent through ``Sigma(n*)``. Consumes Slice L's ``ladder_cumsum``
    at ``n* = N_STAR``. Raises ``ValueError`` on a non-positive ``e_avail_eV``.
    ``ladder``: injected table sets ``Sigma(n*)`` and ``picture``/``kappa``
    are ignored (see :func:`dE_int_pickup_eV`).
    """
    if not (e_avail_eV > 0.0):
        raise ValueError(
            f"f_int_floor requires e_avail_eV > 0 (the per-ion Coulomb budget); "
            f"got {e_avail_eV!r}"
        )
    sigma_nstar = ladder_cumsum(N_STAR, picture=picture, kappa=kappa, ladder=ladder)
    return float(sigma_nstar / e_avail_eV)
