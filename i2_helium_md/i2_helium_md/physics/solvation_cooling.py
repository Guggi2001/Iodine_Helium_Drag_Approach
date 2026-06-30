r"""Newton cooling of the GAH25 solvation-structure energy ``E_solv.struct`` toward
its occupancy-resolved asymptote, with the pair / electrostriction binding split
(Tier-2 Phase-A Slice K).

Pure, stateless energetics: no integrator state, no RNG, no mass and no velocity
(Phase A is mass-agnostic; the friction/force-coefficient ``gamma`` is never touched
here). This module is named for the variable it relaxes -- ``E_solv.struct = E_bind +
E_int`` (MASS K2) -- not for ``E_int``, which is only one additive term *inside* it,
owned by Slice U (``internal_energy_budget.py``). It *consumes* Slice L's
``ladder_cumsum`` for the cumulative gate ``Sigma(n)``; everything composed (the
generative driver, Phase C) and stateful lives later and consumes this module.

Encoded form (TIER2_PHASE_A_IMPLEMENTATION_PLAN.md §2.2 (K); MASS K2)
--------------------------------------------------------------------
The collective first-shell solvation magnitude is occupancy-resolved by the gate
threshold, the cooling asymptote is its negative, and the binding splits into a
pair part (``-Sigma``) and a non-positive electrostriction marginal::

    |S(N)| = |S| * Sigma(N) / Sigma(n*)              (occupancy-resolved collective)
    E_inf(N) = -|S(N)|                                (cooling asymptote)
    E_solv.struct(N) = -Sigma(N)  +  -(|S(N)| - Sigma(N))  +  E_int(N)
                       \__pair__/    \___ E_elec <= 0 ___/      (E_int: Slice U)

``E_inf(N) -> 0`` as ``N -> 0`` (OQ6: the total strip stays dynamically reachable; a
*fixed* full-shell asymptote would mechanically halt shedding). Because the sourced
``|S| = 0.308 eV`` exceeds ``Sigma(n*)`` (mix 0.17-0.19 / X2 0.25-0.28 eV), the
electrostriction marginal is non-positive, and ``|S(n*)| = |S|`` exactly for every
picture/kappa (the ratio ``Sigma(n*)/Sigma(n*) = 1``).

The cooling itself is the closed-form exponential relaxation toward ``E_inf(N)``::

    dE_solv.struct/dt|K2 = -(E_solv.struct - E_inf(N)) / tau          [eV/ps]
    => E_solv.struct(t + dt) = E_inf + (E_solv.struct - E_inf) * exp(-dt/tau)

The exact ``exp(-dt/tau)`` step is ``dt``-robust (one big step == many small steps)
and needs no integrator state. Units follow the package convention: energies in eV,
rate in eV/ps, ``tau`` and ``dt`` in ps; ``N`` (occupancy) is a dimensionless integer.
Vectorised scalar-in -> float / array-in -> ndarray, mirroring ``shell_schedule.py``
and ``dissociation_ladder.py``.
"""

from __future__ import annotations

import numpy as np

from .constants import N_STAR, S_ABS_EV
from .dissociation_ladder import ladder_cumsum


def s_collective_eV(N, *, picture: str = "statistical_mixture", kappa: float,
                    s_abs_eV: float = S_ABS_EV):
    """Occupancy-resolved collective solvation magnitude ``|S(N)|`` [eV].

    ``|S(N)| = s_abs_eV * Sigma(N) / Sigma(n*)`` -- consumes Slice L's
    ``ladder_cumsum``. Pinned to ``s_abs_eV`` at ``N = n*`` (ratio 1) for any
    picture/kappa, and ``-> 0`` as ``N -> 0``. Scalar-in -> float, array-in ->
    ndarray.
    """
    sigma_N = ladder_cumsum(N, picture=picture, kappa=kappa)
    sigma_nstar = ladder_cumsum(N_STAR, picture=picture, kappa=kappa)
    out = s_abs_eV * np.asarray(sigma_N) / sigma_nstar
    return float(out) if np.ndim(N) == 0 else out


def e_infinity_eV(N, *, picture: str = "statistical_mixture", kappa: float,
                  s_abs_eV: float = S_ABS_EV):
    """Cooling asymptote ``E_inf(N) = -|S(N)|`` [eV] (occupancy-resolved).

    Monotone non-increasing in ``N``; ``-> 0`` as ``N -> 0`` (OQ6). Scalar-in ->
    float, array-in -> ndarray.
    """
    s = s_collective_eV(N, picture=picture, kappa=kappa, s_abs_eV=s_abs_eV)
    return -s if np.ndim(N) == 0 else -np.asarray(s)


def e_bind_pair_eV(N, *, picture: str = "statistical_mixture", kappa: float):
    """Pair binding ``E_bind^pair(N) = -Sigma(N)`` [eV] -- thin wrapper over L.

    The discrete-trackable binding term (the cumulative rung cost, negated).
    Scalar-in -> float, array-in -> ndarray.
    """
    sigma_N = ladder_cumsum(N, picture=picture, kappa=kappa)
    return -sigma_N if np.ndim(N) == 0 else -np.asarray(sigma_N)


def e_electrostriction_eV(N, *, picture: str = "statistical_mixture", kappa: float,
                          s_abs_eV: float = S_ABS_EV):
    """Electrostriction marginal ``E_elec(N) = -(|S(N)| - Sigma(N)) <= 0`` [eV].

    The collective (continuous, A8 bath-booked) part of the binding, non-positive
    everywhere because ``|S| > Sigma(n*)`` so ``|S(N)| >= Sigma(N)``. Together with
    :func:`e_bind_pair_eV` it closes the split: ``E_inf = E_bind_pair + E_elec``.
    Scalar-in -> float, array-in -> ndarray.
    """
    s = s_collective_eV(N, picture=picture, kappa=kappa, s_abs_eV=s_abs_eV)
    sigma_N = ladder_cumsum(N, picture=picture, kappa=kappa)
    out = -(np.asarray(s) - np.asarray(sigma_N))
    return float(out) if np.ndim(N) == 0 else out


def newton_cool_step(E_solv_struct_eV, N, *, tau_ps: float, dt_ps: float,
                     picture: str = "statistical_mixture", kappa: float,
                     s_abs_eV: float = S_ABS_EV):
    """One exact exponential Newton-cooling step toward ``E_inf(N)`` [eV].

    ``E_new = E_inf + (E_solv_struct - E_inf) * exp(-dt/tau)`` -- the closed-form
    solution of ``dE/dt = -(E - E_inf)/tau`` over ``dt``. ``dt``-robust (one big
    step equals many small steps) with no integrator state; fixed point at
    ``E_solv_struct = E_inf``.

    Raises ``ValueError`` on a non-positive ``tau_ps`` or a negative ``dt_ps``
    (fail-loud; CLAUDE.md principle 4 -- a non-positive relaxation time is
    unphysical and would invert the cooling, and a negative step gives
    ``exp(+|dt|/tau) > 1`` i.e. silent anti-cooling away from ``E_inf``). ``dt_ps =
    0`` is a valid no-op. Scalar-in -> float, array-in -> ndarray.
    """
    if not (tau_ps > 0.0):
        raise ValueError(
            f"newton_cool_step requires tau_ps > 0 (a non-positive relaxation "
            f"time is unphysical and inverts the cooling); got {tau_ps!r}"
        )
    if dt_ps < 0.0:
        raise ValueError(
            f"newton_cool_step requires dt_ps >= 0 (a negative step gives "
            f"exp(+|dt|/tau) > 1, i.e. anti-cooling away from E_inf); got {dt_ps!r}"
        )
    e_inf = e_infinity_eV(N, picture=picture, kappa=kappa, s_abs_eV=s_abs_eV)
    decay = np.exp(-dt_ps / tau_ps)
    out = np.asarray(e_inf) + (np.asarray(E_solv_struct_eV) - np.asarray(e_inf)) * decay
    if np.ndim(E_solv_struct_eV) == 0 and np.ndim(N) == 0:
        return float(out)
    return out
