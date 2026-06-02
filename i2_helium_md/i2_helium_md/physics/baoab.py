"""BAOAB operator-split ion-stage stepper (Slice 2).

Replaces the baseline :func:`i2_helium_md.physics.leapfrog.velocity_verlet_step`
**for the ion stage only** with a B-A-O-A-B operator-split integrator
(``DRAG_PORT_DESIGN_DECISIONS.md`` §4.6, specified in
``SLICE2_GOALS_baoab_ion_stepper.md``). The neutral stage keeps velocity-Verlet.

Scheme (one step ``dt``)::

    B   v <- v + (dt/2)*a_cons(x)        half-kick, conservative force only
    A   x <- x + (dt/2)*v                half-drift
    O   v <- e^(-gamma*dt/m)*v (+ noise) full Ornstein-Uhlenbeck: drag (+ dormant noise)
    A   x <- x + (dt/2)*v                half-drift
    B   v <- v + (dt/2)*a_cons(x)        half-kick, conservative force only

B and A reuse :func:`leapfrog._kick` / :func:`leapfrog._drift` -- the same
primitives velocity-Verlet uses -- so there is no duplicate integrator physics.
The conservative acceleration ``acc_fn`` (Coulomb + droplet) is the *existing*
ion acceleration, evaluated afresh twice per step (at the entry position and at
the post-step position; no cross-step force caching). O is the only new physics.

Units / mass contract (locked, ``SLICE2_GOALS_baoab_ion_stepper.md`` §2)
------------------------------------------------------------------------
This module is **pure-mechanical and mass-in-amu**. ``gamma`` is a force
coefficient in amu/ps (Slice 1 convention); ``m`` arrives in **amu**; the
damping exponent ``gamma*dt/m`` is dimensionless; dissipated energy is returned
in **amu*A^2/ps^2**. There is no kg and no eV conversion here -- Slice 4 owns the
kg<->amu plumbing and the eV conversion of the dissipated-energy accumulator.
This is the one place mass enters the drag model (Slice 1 is mass-free): once,
in the O-step exponent and the kinetic-energy bookkeeping.

Asymmetric gamma-freeze (intentional -- do not "fix" this)
----------------------------------------------------------
The nonlinear O-step is frozen to an O(dt) approximation by evaluating
``gamma`` once, at a *mixed* point:

* **velocity** frozen at the O-step *input* velocity ``v_in^O`` (the
  post-first-B velocity; the first A does not change velocity). This is the only
  choice available without an implicit solve and makes the returned dissipated
  energy *exact for the frozen gamma*.
* **depth (gate)** evaluated at the *current* O-step position (after the first
  half-drift): in B-A-O-A-B the position is already updated when O runs, so
  ``depth = r_atom - r_droplet`` is freshly and legitimately known.

Because Slice 1 guarantees ``gamma >= 0`` for ``linear_cubic`` (a, b > 0), the
damping factor ``e^(-gamma*dt/m)`` lies in ``(0, 1]`` and the O-step can never
*add* kinetic energy -- the dissipativity is exact.

Noise -- dormant at Tier 0
--------------------------
At ``T_eff = 0`` (the Tier-0 default) the Ornstein-Uhlenbeck step collapses to
pure multiplicative damping and **no random number is drawn**. The Langevin
fluctuation term (``DRAG_PORT_DESIGN_DECISIONS.md`` §1, §5.2) slots in at the
single marked site in the O-step without touching B, A, or the energy return.
Activating it (``T_eff > 0``) is Slice >=3 work (FDT amplitude, RNG draw-order
pinning, the second noise-energy accumulator); until then ``T_eff > 0`` raises
:class:`NotImplementedError` rather than silently running an unvalidated path.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from .leapfrog import AccelFn, Positions, Velocities, _drift, _kick

# A gamma_fn maps (speed, depth) -> friction force-coefficient gamma [amu/ps],
# each array of shape (2N,). In Tier-0 production this is Slice 1's
# ``drag_gamma`` closed over the coefficient bundle and gate steepness.
GammaFn = Callable[[np.ndarray, np.ndarray], np.ndarray]

# The step closure: (pos, vel, dt) -> (pos', vel', E_pot_per_pair, dE_dissip_per_atom).
BaoabStep = Callable[
    [Positions, Velocities, float],
    tuple[Positions, Velocities, np.ndarray, np.ndarray],
]


def make_ion_baoab_step(
    m: np.ndarray,
    droplet_radii: np.ndarray,
    acc_fn: AccelFn,
    gamma_fn: GammaFn,
    *,
    T_eff: float = 0.0,
    rng: np.random.Generator | None = None,
) -> BaoabStep:
    """Build a single BAOAB ion-stage step closure.

    Mirrors :func:`i2_helium_md.physics.leapfrog.make_ion_step`, including its
    call convention: ``dt`` is passed to the returned ``step`` **on each call**
    (not bound in the factory), sourced from ``SimConfig.dt_ion`` by the Slice-4
    driver (as ``ion_propagation_step.py`` does ``dt = cfg.dt_ion``). Slice 4
    rebuilds this closure every step because **mass** changes under future mass
    scenarios, exactly as ``ion_propagation_step.py`` rebuilds ``make_ion_step``;
    the factory is built to support that even though Tier-0 mass is fixed.

    Parameters
    ----------
    m : np.ndarray, shape (2N,)
        Per-atom physical mass in **amu** (not kg). Enters only the O-step
        damping exponent and the dissipated-energy bookkeeping.
    droplet_radii : np.ndarray, shape (2N,)
        Per-atom droplet radius in Angstrom. Used to form
        ``depth = r_atom - droplet_radii`` for the gate.
    acc_fn : callable
        Conservative ion acceleration ``(x, y, z) -> ((ax, ay, az), E_pot)`` with
        accelerations in A/ps^2 and ``E_pot`` per pair (length N) in eV --
        the existing ``leapfrog._ion_accel_fn`` closure. Evaluated twice per step.
    gamma_fn : callable
        Friction force-coefficient ``(speed, depth) -> gamma`` in amu/ps, each
        shape (2N,). Tier-0: Slice 1 ``drag_gamma`` closed over coeffs + steepness.
    T_eff : float, optional
        Effective temperature for the Langevin noise. Tier-0 default ``0.0``
        (noise dormant). ``> 0`` raises :class:`NotImplementedError` (Slice >=3).
    rng : np.random.Generator, optional
        Reserved for the future noise draw. Never consumed at ``T_eff = 0``.

    Returns
    -------
    step : callable
        ``step(pos, vel, dt) -> (new_pos, new_vel, E_pot_per_pair, dE_dissip)``
        where ``dt`` is the timestep in picoseconds (``SimConfig.dt_ion``),
        ``E_pot_per_pair`` has shape (N,) in eV (from the final B's ``acc_fn``
        evaluation, matching ``velocity_verlet_step``) and ``dE_dissip`` has
        shape (2N,) in **amu*A^2/ps^2**, per atom, ``>= 0``.

    Raises
    ------
    NotImplementedError
        If ``T_eff > 0`` -- active Langevin noise is Slice >=3, not yet realised.
    """
    if T_eff < 0.0:
        raise ValueError(f"T_eff must be non-negative, got {T_eff!r}")
    if T_eff > 0.0:
        raise NotImplementedError(
            "active Langevin noise (T_eff > 0) is Slice >=3: it requires the FDT "
            "amplitude, RNG draw-order pinning, and the separate noise-energy "
            "accumulator. Tier-0 runs deterministically at T_eff = 0."
        )

    def step(
        pos: Positions, vel: Velocities, dt: float
    ) -> tuple[Positions, Velocities, np.ndarray, np.ndarray]:
        half_dt = 0.5 * dt
        # --- B: half-kick with the conservative force at the entry position ---
        a0, _ = acc_fn(pos)
        v_half = _kick(vel, a0, half_dt)

        # --- A: half-drift ---
        pos_mid = _drift(pos, v_half, half_dt)

        # --- O: Ornstein-Uhlenbeck step (drag damping; dormant noise site) ---
        new_vel_o, dE_dissip = _o_step(pos_mid, v_half, m, droplet_radii, gamma_fn, dt)

        # --- A: half-drift ---
        new_pos = _drift(pos_mid, new_vel_o, half_dt)

        # --- B: half-kick with the conservative force at the new position ---
        a1, E_pot_end = acc_fn(new_pos)
        new_vel = _kick(new_vel_o, a1, half_dt)

        return new_pos, new_vel, E_pot_end, dE_dissip

    return step


def _o_step(
    pos: Positions,
    vel: Velocities,
    m: np.ndarray,
    droplet_radii: np.ndarray,
    gamma_fn: GammaFn,
    dt: float,
) -> tuple[Velocities, np.ndarray]:
    """The O-step: multiplicative drag damping (Tier-0, noise off).

    ``v_out = e^(-gamma*dt/m) * v``, with ``gamma`` frozen at the input speed
    ``|v|`` (velocity argument) and the current position's depth (gate argument)
    -- the deliberate asymmetric freeze documented in the module docstring.

    Returns the damped velocity and the per-atom dissipated energy
    ``dE = 0.5*m*(|v_in|^2 - |v_out|^2) = 0.5*m*|v_in|^2*(1 - e^(-2*gamma*dt/m))``
    in amu*A^2/ps^2 (``>= 0`` since the factor is in ``[0, 1]``).
    """
    vx, vy, vz = vel
    x, y, z = pos

    speed_in_sq = vx ** 2 + vy ** 2 + vz ** 2
    speed_in = np.sqrt(speed_in_sq)

    # Gate argument: depth at the current (post-first-drift) O-step position.
    r_atom = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    depth = r_atom - droplet_radii

    gamma = gamma_fn(speed_in, depth)            # amu/ps, shape (2N,)
    decay = np.exp(-gamma * dt / m)              # dimensionless, in (0, 1]

    # NOTE (noise site): when T_eff > 0 (Slice >=3) the Langevin term
    #   + sqrt((kB*T_eff/m)*(1 - decay**2)) * xi   (xi ~ N(0, 1), gate already in gamma)
    # is added here, and its injected energy is tracked in a *separate*
    # accumulator from dE_dissip below. Dormant at Tier 0 -> no RNG drawn.
    new_vel = (decay * vx, decay * vy, decay * vz)

    dE_dissip = 0.5 * m * speed_in_sq * (1.0 - decay ** 2)
    return new_vel, dE_dissip
