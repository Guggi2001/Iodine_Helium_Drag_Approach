"""One-timestep advance for the ion propagation loop (pure version).

Public API:

* :class:`IonStepState` -- a frozen dataclass holding everything
  needed to step (positions, velocities, current mass, cumulative
  diagnostics).
* :func:`ion_propagation_step` -- a pure function that takes an
  ``IonStepState`` plus the previous-step distance traveled per atom
  (for the Mode-3 collision sampler), returns a new ``IonStepState``.
  Does not mutate inputs.

Step sequence (mirrors ``vmi_sim_3d_ion_propa.m`` lines ~300-783 with
the unsupported branches stripped out):

1. Leapfrog integrate one ``dt`` -> candidate positions/velocities and
   the per-atom Coulomb potential ``E_pot_coulomb`` (already half-
   per-atom in 2N layout).
2. Compute per-atom depth into droplet (``r1 - droplet_radius``).
3. Compute per-atom velocity-dependent cross section if
   ``cfg.sigma_dependent_on_v``, else use the constant
   ``cfg.geometric_scattering_crosssection_Iplus``.
4. Sample hard-sphere collision events using **Mode 3** (probability
   per step = ``prev_distance * sigma_per_atom * rho_droplet``).
5. Apply collisions, replacing post-leapfrog velocities for colliders
   with elastically-scattered velocities. Track ``actual_dE = E0 - E1``
   per atom.
6. Apply mass attachment: each collider has independent probability
   ``cfg.mass_attach_probability`` to absorb a 4-amu helium atom.
7. Compute energy diagnostics (E_kin, E_pot per atom, cumulative
   E_dissip, cumulative number_of_collisions) and assemble the new
   state.

Out-of-scope branches (raise ValueError if cfg requests them):

- ``hard_sphere_collision_mode != 3``
- ``effusive_dynamics = True``
- ``single_charge_ionization_allowed = True``
- ``additional_droplet_charges > 0``

The legacy MATLAB also has a ``relative_energy_loss_ion`` alternative-
energy-loss-model branch which is not in our SimConfig (both
production input scripts leave it disabled) and is therefore not
implemented here.

The driver (``simulation/ion.py`` -- Step 11d) is responsible for
orchestrating calls to this function, tracking ``prev_distance``
between steps, and optionally storing only every K-th state to a
checkpoint.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from ..config import SimConfig, _MASS_COEFFICIENT_CONSISTENCY_TOL_AMU
from ..physics.collisions import (
    apply_collision,
    sample_collision_events,
    temperature_diagnostic_from_collision,
    velocity_dependent_cross_section,
)
from ..physics.baoab import BaoabStep
from ..physics.constants import EV, MASS_HE_AMU, MASS_I_ION_AMU, U
from ..physics.evaporation import evaporation_step_components
from ..physics.helium_density import rho_he_ratio
from ..physics.internal_energy_budget import pickup_bath_release_eV
from ..physics.mass_jump import continuous_velocity_shed_components
from ..physics.pickup import pickup_step_components
from ..physics.shell_schedule import ShellSchedule
from ..physics.solvation_cooling import e_infinity_eV, newton_cool_step
from ..physics.drag import REALIZED_FORMS
from ..physics.interactions import partner_interaction_ion
from ..physics.leapfrog import make_ion_step
from ..physics.potentials import droplet_potential

#: Tight tolerance [amu] for the per-step ``m == m_I+ + n*m_He`` consistency assert
#: (Slice G; the reduced-mass resets produce exact +/- m_He, so the residual is
#: float-rounding only -- distinct from the loose 8-amu drag mass-insensitivity band).
_M_N_CONSISTENCY_TOL_AMU: float = 1e-6


# ===========================================================================
# State container
# ===========================================================================
@dataclass(frozen=True)
class IonStepState:
    """The minimum-sufficient state to advance the ion propagation.

    All trajectory arrays have shape ``(2N,)`` where ``N`` is the
    number of molecules. Per-atom **dynamic** information includes
    ``mass_kg`` because helium attachment changes mass per atom over
    time -- this is the main difference from ``NeutralStepState``.
    Per-atom static information (``droplet_radii``, ``charge``) is
    NOT carried here; the driver passes it as separate arguments.

    Attributes
    ----------
    x, y, z : np.ndarray, shape (2N,)
        Cartesian positions in Angstrom.
    vx, vy, vz : np.ndarray, shape (2N,)
        Velocities in Angstrom/ps.
    mass_kg : np.ndarray, shape (2N,)
        Per-atom mass in kg. Changes via mass attachment.
    E_kin_eV : np.ndarray, shape (2N,)
        Kinetic energy in eV at this state's time.
    E_pot_eV : np.ndarray, shape (2N,)
        Per-atom potential energy in eV (ion-droplet + half partner Coulomb).
    E_dissip_eV : np.ndarray, shape (2N,)
        Cumulative energy dissipated per atom up to this state's time.
    E_mass_transfer_eV : np.ndarray, shape (2N,)
        Cumulative kinetic-energy defect from helium mass transfer, in
        eV -- attachment (collision path) OR a scheduled shed (Tier-1a
        anchored_discrete drag path). When mass changes by ``dm`` at the
        atom's current velocity, the recomputed E_kin shifts by
        ``1/2 * dm * v^2``; this field accumulates the negative of that
        increment so that
        ``E_kin + E_pot + E_dissip + E_mass_transfer`` is conserved
        (modulo Verlet drift). On the drag path the shed defect is the
        continuous-velocity co-moving-He form (``physics/mass_jump.py``); on the
        collision path it mirrors MATLAB ``E_mass_attach_defect`` at
        vmi_sim_3d_ion_propa.m:762 (renamed from ``E_mass_attach_defect_eV``
        at schema v6).
    E_int_eV : np.ndarray, shape (2N,)
        Per-atom internal-energy reservoir in eV (the Tier-2 fifth
        invariant term, checkpoint schema v7). Zero on the ``fixed`` /
        ``anchored_discrete`` paths (carried through untouched); it
        evolves only under the ``biphasic`` generative driver via the
        S1/S2/K1/K2 budget.
    number_of_collisions : np.ndarray, shape (2N,)
        Cumulative number of hard-sphere collisions per atom.
    time_ps : float
        Time at this state, in picoseconds.
    n_shell : np.ndarray, shape (2N,) or None
        Per-complex integer He-shell count (genuine state on the Tier-2
        ``biphasic`` path, where the pickup/evaporation channels advance ``n``
        independently of the mass reset so a channel/reset drift can be caught
        by the ``m == m_I+ + n*m_He`` per-step invariant). ``None`` on the
        ``fixed`` / ``anchored_discrete`` deterministic paths, where the
        checkpoint writer re-derives ``n`` from the realized mass via ``rint``
        (byte-identical to the pre-Slice-G behaviour). Defaults to ``None``.
    temperature_diagnostic : np.ndarray, shape (3,) or None
        Legacy MATLAB per-step temperature diagnostic
        ``[<T'/T>_actual, <T'/T>_from_mass_ratio, <theta_lab>_rad]``
        averaged over the atoms that collided in the transition that
        produced this state. ``None`` when this state was reconstructed
        from a checkpoint column (no per-step measurement available)
        and an all-NaN ``(3,)`` array when no atom collided in the
        step that produced this state. Mirrors ``diagnostic_array`` at
        ``vmi_sim_3d_ion_propa.m:683``.
    """
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    mass_kg: np.ndarray
    E_kin_eV: np.ndarray
    E_pot_eV: np.ndarray
    E_dissip_eV: np.ndarray
    E_mass_transfer_eV: np.ndarray
    E_int_eV: np.ndarray
    number_of_collisions: np.ndarray
    time_ps: float
    n_shell: np.ndarray | None = None
    temperature_diagnostic: np.ndarray | None = None


# ===========================================================================
# Public step function (pure)
# ===========================================================================
def ion_propagation_step(
    state: IonStepState,
    *,
    cfg: SimConfig,
    droplet_radii: np.ndarray,
    charge: np.ndarray,
    prev_distance_angstrom: np.ndarray | None,
    rng: np.random.Generator,
) -> IonStepState:
    """Advance the ion propagation by one ``dt``. Pure function.

    Parameters
    ----------
    state : IonStepState
        Current state (read only; not mutated).
    cfg : SimConfig
        Simulation config.
    droplet_radii : np.ndarray, shape (2N,)
        Per-atom droplet radius in Angstrom (constant across the run).
    charge : np.ndarray, shape (2N,)
        Per-atom charges (constant across the run; all 1.0 in our scope).
    prev_distance_angstrom : np.ndarray of shape (2N,) or None
        Distance traveled per atom during the previous step (used by
        the Mode-3 collision sampler). Pass ``None`` for the very
        first step of a run -- collisions are then disabled because
        there is no previous distance to use.
    rng : np.random.Generator
        Reproducible RNG for collision sampling and mass attachment.

    Returns
    -------
    IonStepState
        The new state at ``state.time_ps + dt_ion``.

    Raises
    ------
    ValueError
        If cfg requests out-of-scope features.
    """
    _check_scope(cfg)

    dt = cfg.dt_ion

    # 1. Leapfrog integration. The ion step closure must be rebuilt each
    #    iteration because mass changes due to attachment (Option C from
    #    the design discussion -- mass carried in state, closure rebuilt
    #    inside the pure step function).
    step_fn = make_ion_step(cfg, state.mass_kg, droplet_radii, charge)
    (x1, y1, z1), (vx1, vy1, vz1), E_pot_coulomb_per_pair = step_fn(
        (state.x, state.y, state.z),
        (state.vx, state.vy, state.vz),
        dt,
    )
    # NOTE: make_ion_step returns the partner Coulomb energy per pair
    # (shape (N,)), matching the neutral-branch convention. We split it
    # half-and-half between the two atoms of each molecule below to
    # populate the per-atom (2N,) E_pot array.

    # 2. Depth into droplet.
    depth = _depth(x1, y1, z1, droplet_radii)

    # 3. Pre-collision speed and energy. We compute these AFTER the
    #    leapfrog step so the collision uses the energy at the new
    #    position (matches MATLAB line 384: E0 = (v1*100)^2 * m / 2 / eV).
    v1_speed_sq = vx1 ** 2 + vy1 ** 2 + vz1 ** 2
    v1_speed = np.sqrt(v1_speed_sq)
    E0_eV = _E_kin_eV(state.mass_kg, vx1, vy1, vz1)

    # 4. Per-atom cross section (constant or v-dependent).
    if cfg.sigma_dependent_on_v:
        sigma_per_atom = velocity_dependent_cross_section(
            v1_speed,
            sigma_0_angstrom_sq=cfg.geometric_scattering_crosssection_Iplus,
            exponent=cfg.sigma_ion_exponent,
        )
    else:
        sigma_per_atom = cfg.geometric_scattering_crosssection_Iplus

    # 5. Mode-3 collision sampling. Uses the *previous* step's distance.
    n_atoms = state.x.shape[0]
    if prev_distance_angstrom is None:
        b_collision = np.zeros(n_atoms, dtype=bool)
    else:
        b_collision = sample_collision_events(
            distance_travelled_angstrom=prev_distance_angstrom,
            depth_angstrom=depth,
            E0_eV=E0_eV,
            sigma_angstrom_sq=sigma_per_atom,
            E_min_eV=cfg.E_min_eV,
            rng=rng,
        )

    # 6. Apply collisions to flagged atoms. Capture diagnostics so we
    #    can record the legacy MATLAB per-step temperature accumulator
    #    ([<T'/T>_actual, <T'/T>_mass, <theta>_rad]).
    masses_amu = state.mass_kg / U
    vx_after, vy_after, vz_after, dE_eV, coll_diag = apply_collision(
        vx=vx1, vy=vy1, vz=vz1,
        masses_amu=masses_amu,
        b_collision=b_collision,
        scatter_mass_amu=cfg.scatter_mass_ion_amu,
        neutral_scatter_angle_std_deg=cfg.ion_scatter_angle_std_deg,
        rng=rng,
        return_diagnostics=True,
    )
    temperature_diagnostic_step = temperature_diagnostic_from_collision(coll_diag)

    # 7. Mass attachment. Each collider has independent probability
    #    cfg.mass_attach_probability to gain 4 amu (one He atom).
    #    NOTE: MATLAB draws this random number for ALL atoms (line 727),
    #    not just colliders, then masks with b_collision. We mirror that
    #    so the rng stream is deterministic with the same call pattern.
    mass_attach_trial = rng.uniform(0.0, 1.0, size=n_atoms)
    b_mass_attach = (mass_attach_trial < cfg.mass_attach_probability) & b_collision
    new_mass_kg = state.mass_kg + b_mass_attach * 4.0 * U

    # 8. Energy diagnostics for the new state. Use the NEW (possibly
    #    increased) mass for E_kin -- this matches MATLAB line 761,
    #    which uses ``mass_i(:, t_id+1)`` (post-attachment mass).
    v_post_sq = vx_after ** 2 + vy_after ** 2 + vz_after ** 2
    E_kin_new_eV = _E_kin_eV(new_mass_kg, vx_after, vy_after, vz_after)

    # E_pot at NEW positions: ion-droplet binding + half partner Coulomb.
    E_pot_new_eV = _E_pot_per_atom(depth, E_pot_coulomb_per_pair, cfg)

    # 9. Cumulative bookkeeping.
    E_dissip_new = state.E_dissip_eV + dE_eV
    n_coll_new = state.number_of_collisions + b_collision.astype(state.number_of_collisions.dtype)

    # Kinetic-energy defect from mass attachment. MATLAB
    # vmi_sim_3d_ion_propa.m:762 computes
    #   E_defect[t+1] = E_defect[t] - (m_new - m_old) * (v*100)^2 / 2 / eV
    # using the post-collision, post-attachment velocity. The sign is
    # negative because attaching mass at velocity v adds ``1/2 dm v^2``
    # of fictitious E_kin; the defect compensates so that the system
    # invariant E_kin + E_pot + E_dissip + E_mass_transfer is
    # conserved up to Verlet drift.
    mass_diff_kg = new_mass_kg - state.mass_kg
    dE_defect_eV = -0.5 * mass_diff_kg * (v_post_sq * 100.0 ** 2) / EV
    E_mass_transfer_new = state.E_mass_transfer_eV + dE_defect_eV

    return IonStepState(
        x=x1, y=y1, z=z1,
        vx=vx_after, vy=vy_after, vz=vz_after,
        mass_kg=new_mass_kg,
        E_kin_eV=E_kin_new_eV,
        E_pot_eV=E_pot_new_eV,
        E_dissip_eV=E_dissip_new,
        E_mass_transfer_eV=E_mass_transfer_new,
        E_int_eV=state.E_int_eV,  # untouched on the collision/attach path
        number_of_collisions=n_coll_new,
        time_ps=state.time_ps + dt,
        n_shell=state.n_shell,  # carried through (None on the collision path)
        temperature_diagnostic=temperature_diagnostic_step,
    )


# ===========================================================================
# Drag-branch step function (pure) -- Slice 4, Tier 0
# ===========================================================================
def baoab_propagation_step(
    state: IonStepState,
    *,
    step: BaoabStep,
    cfg: SimConfig,
    droplet_radii: np.ndarray,
) -> IonStepState:
    """Advance the ion propagation by one ``dt`` via the drag (BAOAB) path. Pure.

    The Tier-0 drag-branch analog of :func:`ion_propagation_step`: **no** collision
    sampling, **no** mass attachment. The conservative B/A kicks and the
    dissipative O-step both live inside ``step`` (a BAOAB closure built by the
    driver via :func:`i2_helium_md.physics.baoab.make_ion_baoab_step`); this
    function does thin per-step accounting only -- depth, eV energies, the
    dissipation unit-conversion, and the Tier-0 checkpoint fills.

    Parameters
    ----------
    state : IonStepState
        Current state (read only; not mutated).
    step : BaoabStep
        Pre-built BAOAB closure
        ``(pos, vel, dt) -> (pos', vel', E_pot_per_pair, dE_dissip)`` where
        ``E_pot_per_pair`` has shape (N,) in eV and ``dE_dissip`` has shape (2N,)
        in **amu*A^2/ps^2** (per atom, ``>= 0``). The driver rebuilds it every
        step (mass-driven; Tier-1-ready though Tier-0 mass is fixed).
    cfg : SimConfig
        Simulation config. ``dt_ion`` is the timestep.
    droplet_radii : np.ndarray, shape (2N,)
        Per-atom droplet radius in Angstrom (constant across the run).

    Returns
    -------
    IonStepState
        The new state at ``state.time_ps + cfg.dt_ion``. At Tier 0: ``mass_kg``
        unchanged (fixed mass), ``E_mass_transfer_eV`` and
        ``number_of_collisions`` carried unchanged (both 0 under the drag branch),
        and ``temperature_diagnostic`` an all-NaN ``(3,)`` sentinel (no collisions
        to diagnose).
    """
    dt = cfg.dt_ion

    (x1, y1, z1), (vx1, vy1, vz1), E_pot_coulomb_per_pair, dE_dissip = step(
        (state.x, state.y, state.z),
        (state.vx, state.vy, state.vz),
        dt,
    )

    # Per-atom depth into droplet (shared geometry lift).
    depth = _depth(x1, y1, z1, droplet_radii)

    # E_kin at the new velocity. Mass is fixed at Tier 0, so the post-step mass
    # is the carried mass (no attachment).
    E_kin_new_eV = _E_kin_eV(state.mass_kg, vx1, vy1, vz1)
    E_pot_new_eV = _E_pot_per_atom(depth, E_pot_coulomb_per_pair, cfg)

    # Dissipated energy: amu*A^2/ps^2 -> eV via the shared baseline idiom, so
    # E_dissip_eV stays in a consistent eV with E_kin_eV / E_pot_eV.
    dE_dissip_eV = _amu_ang2_ps2_to_eV(dE_dissip)
    E_dissip_new = state.E_dissip_eV + dE_dissip_eV

    return IonStepState(
        x=x1, y=y1, z=z1,
        vx=vx1, vy=vy1, vz=vz1,
        mass_kg=state.mass_kg,                              # fixed mass (Tier 0)
        E_kin_eV=E_kin_new_eV,
        E_pot_eV=E_pot_new_eV,
        E_dissip_eV=E_dissip_new,
        E_mass_transfer_eV=state.E_mass_transfer_eV,  # stays 0 (no attach)
        E_int_eV=state.E_int_eV,  # carried through (0 on Tier-0; evolved by biphasic_step)
        number_of_collisions=state.number_of_collisions,       # stays 0 (no collisions)
        time_ps=state.time_ps + dt,
        n_shell=state.n_shell,  # carried through (None off biphasic; genuine state under it)
        temperature_diagnostic=np.full(3, np.nan, dtype=float),
    )


# ===========================================================================
# Tier-1a continuous-velocity shed pre-step (SQ2 + SQ3) -- variable-mass drag path
# ===========================================================================
def shed_step(
    state: IonStepState,
    schedule: ShellSchedule,
    next_shed_idx: int,
    dt: float,
    *,
    m_he_amu: float = MASS_HE_AMU,
) -> tuple[IonStepState, int]:
    """Apply at most one scheduled continuous-velocity shed before BAOAB.

    The Tier-1a ``anchored_discrete`` pre-step. If the next pending shed
    (``schedule.events[next_shed_idx]``) fires within this step's window --
    its analytic fire time is ``<= state.time_ps + dt`` -- copy **all** atoms'
    velocities unchanged, drop the scheduled number of He atoms from the
    (uniform) complex mass, book the per-atom co-moving-He kinetic energy into
    ``E_mass_transfer_eV``, and advance the pointer. **At most one shed per
    call** (the plan's <=1/step rule):
    if the schedule is dense relative to ``dt`` the surplus events fire on the
    following steps, so the total shed count is ``len(events)`` independent of
    ``dt`` (the jump-step measure-zero property).

    Applying the jump here, at the step seam, means the rebuilt
    :func:`~i2_helium_md.physics.baoab.make_ion_baoab_step` closure reads the
    post-shed mass ``m+`` for *both* the conservative kicks and the drag O-step
    (SQ3) -- the BAOAB O-step itself (SQ1) is reused unchanged. The order
    reduction vs an in-step "B/A -> jump -> O" placement is ``O(dt)`` in the
    first half-kick's mass and benign as ``dt -> 0`` (jumps are ``dt``-independent
    in count; plan §2).

    Parameters
    ----------
    state : IonStepState
        Current state (read only; not mutated). Its ``mass_kg`` is uniform across
        atoms under the anchored schedule.
    schedule : ShellSchedule
        The anchored He-shell schedule (Slice S), source of the fire times and the
        per-event pre-shed mass.
    next_shed_idx : int
        Index of the next not-yet-fired event in ``schedule.events``.
    dt : float
        Ion timestep [ps] (``cfg.dt_ion``).
    m_he_amu : float, optional
        Shed He mass [amu] (default :data:`MASS_HE_AMU`).

    Returns
    -------
    (new_state, new_next_shed_idx) : tuple[IonStepState, int]
        The post-shed state and the advanced pointer. When no shed fires this
        step, ``state`` and ``next_shed_idx`` are returned unchanged.
    """
    events = schedule.events
    if next_shed_idx >= len(events):
        return state, next_shed_idx
    event = events[next_shed_idx]
    if event.time_ps > state.time_ps + dt:
        return state, next_shed_idx

    n_removed = event.n_before - event.n_after
    vx_p, vy_p, vz_p, m_plus_amu, dE_amu = continuous_velocity_shed_components(
        state.vx,
        state.vy,
        state.vz,
        event.mass_before_amu,
        m_he_amu=m_he_amu,
        n_removed=n_removed,
    )
    # amu*A^2/ps^2 -> eV via the shared baseline idiom, so the booked defect
    # stays in a consistent eV with E_kin/E_pot/E_dissip.
    dE_eV = _amu_ang2_ps2_to_eV(dE_amu)
    new_state = replace(
        state,
        vx=vx_p, vy=vy_p, vz=vz_p,
        mass_kg=np.full_like(state.mass_kg, m_plus_amu * U),
        E_mass_transfer_eV=state.E_mass_transfer_eV + dE_eV,
    )
    return new_state, next_shed_idx + 1


# ===========================================================================
# Tier-2 biphasic generative pre-step (Phase-C Slice G) -- variable-mass drag path
# ===========================================================================
def biphasic_step(
    state: IonStepState,
    *,
    rng: np.random.Generator,
    cfg: SimConfig,
    droplet_radii: np.ndarray,
    gate_steepness: float,
) -> IonStepState:
    r"""Apply one biphasic generative pre-step (K2 cooling + <=1 mass event). Pure.

    The Tier-2 production analog of :func:`shed_step`: a per-ion **pre-step seam**
    that runs *before* the driver rebuilds the BAOAB closure at the post-jump mass
    ``m+``. It composes the accepted Phase-A/B channels into one step with the
    one-event-per-step / frozen-draw-order / 5-term-closing contract, and does **no**
    position or time advance (that is :func:`baoab_propagation_step`, run after).

    Per ion, in fixed order (MASS §6; ``scratchpad`` accounting oracle):

    1. **K2 cooling** -- ``E_int -> E_int*exp(-dt/tau)`` via
       :func:`~i2_helium_md.physics.solvation_cooling.newton_cool_step` on
       ``E_solv.struct = E_inf(n) + E_int`` (rule-1 reuse; at fixed ``n`` the
       asymptote cancels). The drain ``E_int*(1-decay) >= 0`` books to ``E_dissip``.
    2. **Evaporation draw** (frozen first), on the **post-cooling** ``E_int`` so the
       self-bound gate opens exactly as cooling drains ``E_int`` below ``Sigma(n)``:
       on a fire ``n->n-1``, ``m->m-m_He`` (cold-shed reset), ``E_int -= D_0(n)``
       (K1), the mechanical defect ``->`` ``E_mass_transfer``.
    3. **Pickup draw** (frozen second), on the **post-shed** ``(n, m, v)`` -- so a
       both-fire ion applies **shed then pickup** (net ``n`` unchanged, both bookings
       applied): on a fire ``n->n+1``, ``m->m+m_He`` (capture reset),
       ``E_int += f_ret*D_0(n+1)`` (S1), the bath remainder
       ``(1-f_ret)*D_0(n+1) -> E_dissip``, the capture defect ``-> E_mass_transfer``.

    The **binding-potential fold** ``e_bind_pair(n) -> E_pot`` that closes the ``+/-
    D_0`` bookings is applied by the **driver** after :func:`baoab_propagation_step`
    recomputes the MD ``E_pot`` (it is a pure function of the post-event ``n``), so it
    is *not* set here. ``E_kin`` / ``E_pot`` are likewise left to the BAOAB step.

    Parameters
    ----------
    state : IonStepState
        Current state (read only; not mutated). ``n_shell`` must be genuine ``(2N,)``
        state (raised if ``None``).
    rng : np.random.Generator
        Injected RNG. Draws exactly one ``rng.random(size=2N)`` for evaporation, then
        one for pickup (the frozen shed-then-pickup order).
    cfg : SimConfig
        Reads the Phase-A/B knob surface (tau, |S|, picture, kappa, nu, f_ret,
        lambda_0, p, cap, forms, dof/gate overrides) plus ``dt_ion``.
    droplet_radii : np.ndarray, shape (2N,)
        Per-atom droplet radius [A] (for the pickup density gate depth).
    gate_steepness : float
        The erf-gate steepness [A] for the pickup occupancy density (the driver
        passes the resolved drag-gate steepness so density and drag share a surface).

    Returns
    -------
    IonStepState
        Advanced ``vx/vy/vz``, ``mass_kg``, ``n_shell``, ``E_int_eV``, ``E_dissip_eV``,
        ``E_mass_transfer_eV`` at the **same** ``time_ps`` and positions.

    Raises
    ------
    ValueError
        If ``state.n_shell`` is ``None`` (biphasic requires genuine occupancy state).
    NotImplementedError
        If ``cfg.helium_density_profile`` selects the declared-but-unbuilt
        ``'tabulated'`` arm, or ``cfg.dissociation_ladder`` selects the
        ``'tabulated'`` fallback (the biphasic energetics compose the Form-U
        module functions directly; a :class:`TabulatedLadder` has no config
        data path). Both are lazy point-of-use refusals (rule-2 contract).
    AssertionError
        If the post-step ``m == m_I+ + n*m_He`` consistency drifts beyond
        :data:`_M_N_CONSISTENCY_TOL_AMU` (the deferred Phase-B Q3 guard).
    """
    if state.n_shell is None:
        raise ValueError(
            "biphasic_step requires genuine n_shell state on the IonStepState "
            "(got None); the biphasic driver reads it from the v7 checkpoint column."
        )

    # Ladder-form refusal at point-of-use: every energetics call below (K2
    # cooling, evaporation gate/drain, pickup heat/bath) composes the Form-U
    # module functions; the 'tabulated' fallback is a module-level
    # TabulatedLadder object with no config data path, so selecting it would
    # otherwise silently run the Form-U ladder -- the same lazy rule-2 contract
    # as helium_density_profile='tabulated' below.
    if cfg.dissociation_ladder != "form_u":
        raise NotImplementedError(
            f"dissociation_ladder={cfg.dissociation_ladder!r} is not wired into "
            "the biphasic driver; only the 'form_u' ladder is composed "
            "(physics/dissociation_ladder.py -- the TabulatedLadder fallback "
            "has no config data path yet)."
        )

    dt = cfg.dt_ion
    picture = cfg.ladder_electronic_picture
    kappa = cfg.ladder_steepness
    s_abs = cfg.solv_struct_asymptote_eV
    f_ret = cfg.internal_energy_retained_fraction
    n = np.asarray(state.n_shell, dtype=float)
    m_amu = state.mass_kg / U

    # 1. K2 cooling: reuse newton_cool_step on E_solv.struct = E_inf(n) + E_int. At
    #    fixed n the asymptote cancels, leaving E_int*exp(-dt/tau); the drain books
    #    to E_dissip (MASS §6 K2; rule-1 single source of the cooling form).
    e_inf = np.asarray(e_infinity_eV(n, picture=picture, kappa=kappa, s_abs_eV=s_abs))
    e_solv_cooled = np.asarray(
        newton_cool_step(
            e_inf + state.E_int_eV, n,
            tau_ps=cfg.internal_energy_cooling_tau_ps, dt_ps=dt,
            picture=picture, kappa=kappa, s_abs_eV=s_abs,
        )
    )
    E_int = e_solv_cooled - e_inf
    E_dissip = state.E_dissip_eV + (state.E_int_eV - E_int)   # cooling drain >= 0

    # Pickup occupancy density gate at the pre-step positions (shared drag surface).
    # Only the analytic erf-complement profile is built; 'tabulated' is a valid
    # config enum (declared rule-2 arm) whose sourced data array does not exist
    # yet, so selecting it fails lazily here at point-of-use -- the same contract
    # as pickup_rate_form='sweeping' / he_capture_velocity='thermal'.
    if cfg.helium_density_profile != "erf_complement":
        raise NotImplementedError(
            f"helium_density_profile={cfg.helium_density_profile!r} is a "
            "declared-but-unbuilt rule-2 arm for the biphasic pickup density "
            "gate; only 'erf_complement' is built (physics/helium_density.py)."
        )
    depth = _depth(state.x, state.y, state.z, droplet_radii)
    rho_ratio = rho_he_ratio(depth, steepness=gate_steepness)

    # 2. Evaporation draw FIRST (frozen order), on the post-cooling E_int.
    n_e, m_e, vx_e, vy_e, vz_e, dE_int_e, dE_mt_e, _fired_e = evaporation_step_components(
        rng=rng, E_int_eV=E_int, n=n, vx=state.vx, vy=state.vy, vz=state.vz,
        m_amu=m_amu, nu=cfg.evap_rate_prefactor_per_ps, picture=picture, kappa=kappa,
        dt_ps=dt, evap_rrk_dof=cfg.evap_rrk_dof,
        gate_onset_eV=cfg.gate_onset_override_eV,
    )
    E_int = E_int + dE_int_e                                  # K1 drain (-D_0(n)) on fire
    E_mass_transfer = state.E_mass_transfer_eV + _amu_ang2_ps2_to_eV(dE_mt_e)

    # 3. Pickup draw SECOND, on the post-shed (n, m, v) -> shed-then-pickup if both.
    n_p, m_p, vx_p, vy_p, vz_p, dE_int_p, dE_mt_p, fired_p = pickup_step_components(
        rng=rng, n=n_e, vx=vx_e, vy=vy_e, vz=vz_e, m_amu=m_e, rho_ratio=rho_ratio,
        lambda0=cfg.pickup_rate_coefficient, f_ret=f_ret, picture=picture, kappa=kappa,
        dt_ps=dt, p=cfg.pickup_occupancy_exponent, cap=cfg.pickup_occupancy_cap,
        pickup_rate_form=cfg.pickup_rate_form,
        he_capture_velocity=cfg.he_capture_velocity,
    )
    E_int = E_int + dE_int_p                                  # S1 heat (+f_ret*D_0) on fire
    E_mass_transfer = E_mass_transfer + _amu_ang2_ps2_to_eV(dE_mt_p)
    # S1 bath remainder (1-f_ret)*D_0(n+1) -> E_dissip (n = pre-pickup = post-shed n_e).
    bath = np.where(
        fired_p,
        pickup_bath_release_eV(n_e, f_ret=f_ret, picture=picture, kappa=kappa),
        0.0,
    )
    E_dissip = E_dissip + bath

    n_p = np.asarray(n_p, dtype=float)
    m_p_amu = np.asarray(m_p, dtype=float)

    # m <-> n consistency (deferred Phase-B Q3 guard, owned here): the two
    # bookkeepings must never drift (a cheap structural catch on a channel/reset bug).
    if not np.allclose(
        m_p_amu, MASS_I_ION_AMU + n_p * MASS_HE_AMU,
        rtol=0.0, atol=_M_N_CONSISTENCY_TOL_AMU,
    ):
        raise AssertionError(
            "biphasic_step m<->n consistency drift: mass and n_shell disagree beyond "
            f"{_M_N_CONSISTENCY_TOL_AMU} amu (max gap "
            f"{float(np.max(np.abs(m_p_amu - (MASS_I_ION_AMU + n_p * MASS_HE_AMU)))):.3e} "
            "amu). The pickup/evaporation resets and the integer n counter must stay "
            "in lockstep (capture/cold_shed are the single mass source)."
        )

    return replace(
        state,
        vx=np.asarray(vx_p, dtype=float),
        vy=np.asarray(vy_p, dtype=float),
        vz=np.asarray(vz_p, dtype=float),
        mass_kg=m_p_amu * U,
        n_shell=n_p,
        E_int_eV=E_int,
        E_dissip_eV=E_dissip,
        E_mass_transfer_eV=E_mass_transfer,
    )


# ===========================================================================
# Internal helpers
# ===========================================================================
def _depth(x, y, z, droplet_radii):
    """Per-atom depth into the droplet ``r_atom - r_droplet`` (negative inside).

    Pure geometry, shared by the collision and drag paths (CLAUDE.md rule 1).
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return r - droplet_radii


def _amu_ang2_ps2_to_eV(dE_amu_ang2_ps2):
    """Convert an energy increment from MD units [amu*A^2/ps^2] to eV.

    The baseline idiom (amu->kg via ``U``, A/ps->m/s via 100, J->eV via ``EV``),
    the single source for the drag-path energy bookings (CLAUDE.md rule 1):
    ``baoab_propagation_step`` (dE_dissip), ``shed_step`` (shed defect), and
    ``biphasic_step`` (channel dE_mass_transfer). The operation order
    ``((dE*U)*100^2)/EV`` is exactly the pre-refactor inline sequence, so the
    Tier-0/1a outputs stay byte-identical. Scalar or ndarray in -> same out.
    """
    return dE_amu_ang2_ps2 * U * (100.0 ** 2) / EV


def _E_kin_eV(mass_kg, vx, vy, vz):
    """Per-atom kinetic energy in eV from mass in kg and velocity in A/ps.

    ``0.5 * m * v^2`` with v in m/s (A/ps * 100), result in eV. The single source
    of the ion-stage E_kin idiom (matches ``ion_initial_state._compute_E_kin_per_atom``
    and the collision step); physics-free, so safe to share.
    """
    v_sq = vx ** 2 + vy ** 2 + vz ** 2
    return 0.5 * mass_kg * (v_sq * 100.0 ** 2) / EV


def _E_pot_per_atom(depth, E_pot_coulomb_per_pair, cfg):
    """Per-atom ion potential in eV: ion-droplet binding + half partner Coulomb.

    The conservative ion potential, identical for the collision and drag paths
    (both consume the same per-pair Coulomb from the ion acceleration). The
    per-pair Coulomb is split half-and-half between a molecule's two atoms,
    matching the neutral ``propagation_step`` convention.
    """
    E_droplet_eV = droplet_potential(
        depth,
        steepness=cfg.potential_steepness,
        binding_energy=cfg.binding_energy_I_ion_eV,
    )
    E_partner_per_atom = np.tile(E_pot_coulomb_per_pair, 2) / 2.0
    return E_droplet_eV + E_partner_per_atom


def _check_drag_scope(cfg: SimConfig, initial_mass_kg: np.ndarray) -> None:
    """Refuse to run the drag branch outside the deterministic mass envelope.

    The drag-branch analog of :func:`_check_scope`. ``_check_scope`` demands
    collision mode 3, which is irrelevant under drag; the drag path instead
    asserts the deterministic / realised-form envelope for ``fixed`` and the
    Tier-1a ``anchored_discrete`` schedule (``DRAG_PORT_DESIGN_DECISIONS.md``
    §6.4, form set widened by the METHOD_B §10 form phase) so an out-of-scope
    drag config fails at the driver, not deep in a half-implemented path.

    Distinct from ``config.check_drag_config``: that validates the config's
    *internal consistency* (form agreement, dissipativity, mass<->coefficient
    pairing); this validates the *Tier-0 runnability envelope* (config-valid is
    not the same as Tier-0-runnable). Mirrors how ``make_ion_baoab_step`` raises
    :class:`NotImplementedError` for the not-yet-realised tiers.

    Realized-mass trip-wire (§6.5). The final check reads the **realized** initial
    ion mass ``initial_mass_kg`` -- the mass the stepper will actually integrate,
    downstream of the ``build_initial_ion_state`` m_eff override -- *not* a config
    field. Reading the config would merely echo what the override just set and
    guard nothing. The drag law was extracted at ``m_eff``; integrating the
    ``fixed`` scenario at any other inertia (e.g. the inherited bare-I+ ~127 amu)
    silently applies the calibrated law at the wrong mass. A future refactor or a
    bypassing path that leaves the mass at ~127 amu trips this on the ~76-amu gap.
    Reuses the §6.5 ~8 amu mass-insensitivity band (no third tolerance).

    Parameters
    ----------
    cfg : SimConfig
    initial_mass_kg : np.ndarray, shape (2N,)
        The realized per-atom initial ion mass in kg (``ckpt.mass_kg`` after the
        Change-A override), i.e. what the BAOAB stepper integrates.
    """
    unsupported = []
    if cfg.noise_form != "none":
        unsupported.append(
            f"noise_form={cfg.noise_form!r} (active Langevin noise is Tier 3)"
        )
    if cfg.mass_scenario not in ("fixed", "anchored_discrete", "biphasic"):
        unsupported.append(
            f"mass_scenario={cfg.mass_scenario!r} (drag branch supports 'fixed', the "
            "Tier-1a 'anchored_discrete', and the Tier-2 generative 'biphasic')"
        )
    if cfg.drag_form not in REALIZED_FORMS:
        unsupported.append(
            f"drag_form={cfg.drag_form!r} (realised forms: {REALIZED_FORMS}; "
            "'threshold' is reserved and raises NotImplementedError in "
            "physics/drag.py)"
        )
    if cfg.effusive_dynamics:
        unsupported.append("effusive_dynamics")
    if cfg.single_charge_ionization_allowed:
        unsupported.append("single_charge_ionization_allowed")
    if cfg.additional_droplet_charges > 0:
        unsupported.append("additional_droplet_charges > 0")

    if unsupported:
        raise NotImplementedError(
            "drag-branch ion propagation supports the deterministic fixed, the "
            "Tier-1a anchored_discrete, and the Tier-2 biphasic mass scenarios, and "
            "does not support: "
            + ", ".join(unsupported)
            + ". Envelope = mass_scenario in ('fixed', 'anchored_discrete', "
            "'biphasic'), noise_form='none', a realised drag_form (METHOD_B §10 "
            "form phase)."
        )

    # Realized-mass trip-wire: under `fixed`, the integration mass must equal the
    # drag law's extraction mass m_eff within the §6.5 mass-insensitivity band.
    # Per-atom (max), since the override fills uniformly. Skipped for
    # `anchored_discrete`: its mass legitimately runs n=21 -> 14 (210.955 ->
    # 182.936 amu, off m_eff by more than the band by construction), defended by
    # the §6.6 mid-window argument rather than the m_eff band.
    if cfg.mass_scenario != "fixed":
        return
    max_mass_gap_amu = float(np.max(np.abs(initial_mass_kg / U - cfg.m_eff_amu)))
    if max_mass_gap_amu > _MASS_COEFFICIENT_CONSISTENCY_TOL_AMU:
        raise NotImplementedError(
            "drag-branch ion propagation must integrate at the drag-law "
            f"extraction mass m_eff_amu={cfg.m_eff_amu}, but the realized initial "
            f"ion mass is off by {max_mass_gap_amu:.3f} amu (> "
            f"{_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU} amu, §6.5). The fixed "
            "scenario must run at m_eff (DRAG_PORT_DESIGN_DECISIONS.md §6.5); "
            "the build_initial_ion_state m_eff override appears not to have "
            "fired."
        )


def _check_scope(cfg: SimConfig) -> None:
    """Refuse to step when cfg requests features not yet implemented."""
    if cfg.hard_sphere_collision_mode != 3:
        raise ValueError(
            "Only collision mode 3 is implemented; got "
            f"hard_sphere_collision_mode={cfg.hard_sphere_collision_mode}."
        )

    unsupported = []
    if cfg.effusive_dynamics:
        unsupported.append("effusive_dynamics")
    if cfg.single_charge_ionization_allowed:
        unsupported.append("single_charge_ionization_allowed")
    if cfg.additional_droplet_charges > 0:
        unsupported.append("additional_droplet_charges > 0")
    # Note: legacy MATLAB also has a `relative_energy_loss_ion` flag (a
    # different model from hard-sphere). It isn't in our SimConfig
    # because both production input scripts leave it disabled. If we
    # ever add the field, mirror the check pattern above.

    if unsupported:
        raise ValueError(
            "ion_propagation_step does not support: "
            + ", ".join(unsupported)
            + ". The two production input scripts leave all of these "
            "at their default-disabled values."
        )


# ===========================================================================
# Convenience helpers for checkpoint I/O
# ===========================================================================
def ion_state_from_checkpoint_column(ckpt, t_id: int) -> IonStepState:
    """Extract an ``IonStepState`` from column ``t_id`` of an IonCheckpoint.

    Used by the driver to bootstrap the inner loop after
    :func:`build_initial_ion_state`. Copies the underlying arrays so the
    state is independent of the checkpoint's storage. Mirrors the
    neutral helper in ``propagation_step.py``; the only ion-specific
    differences are that ``mass_kg`` lives in the state (sourced from
    ``mass_history_kg[:, t_id]`` because mass changes via attachment)
    and that ``number_of_collisions`` is tracked.
    """
    return IonStepState(
        x=ckpt.positions_x[:, t_id].copy(),
        y=ckpt.positions_y[:, t_id].copy(),
        z=ckpt.positions_z[:, t_id].copy(),
        vx=ckpt.velocities_x[:, t_id].copy(),
        vy=ckpt.velocities_y[:, t_id].copy(),
        vz=ckpt.velocities_z[:, t_id].copy(),
        mass_kg=ckpt.mass_history_kg[:, t_id].copy(),
        E_kin_eV=ckpt.E_kin_eV[:, t_id].copy(),
        E_pot_eV=ckpt.E_pot_eV[:, t_id].copy(),
        E_dissip_eV=ckpt.E_dissip_eV[:, t_id].copy(),
        E_mass_transfer_eV=ckpt.E_mass_transfer_eV[:, t_id].copy(),
        E_int_eV=ckpt.E_int_eV[:, t_id].copy(),
        number_of_collisions=ckpt.number_of_collisions[:, t_id].copy(),
        time_ps=float(ckpt.time_ps[t_id]),
        n_shell=ckpt.n_shell[:, t_id].copy(),
    )


def write_ion_state_to_checkpoint_column(
    state: IonStepState,
    ckpt,
    t_id: int,
    *,
    mass_scenario: str,
) -> None:
    """Write an ``IonStepState`` into column ``t_id`` of an IonCheckpoint.

    Parameters
    ----------
    state : IonStepState
        The state to store (all ``(2N,)`` per-atom fields).
    ckpt : IonCheckpoint
        Target checkpoint; column ``t_id`` of every trajectory array is written.
    t_id : int
        Column index.
    mass_scenario : str
        The run's ``cfg.mass_scenario``; selects the ``n_shell`` source (Slice G).
        Under ``"biphasic"`` the He-shell count is **genuine state** carried on
        ``state.n_shell`` (shape ``(2N,)``, advanced by the pickup/evaporation
        channels independently of the mass reset), so it is stored verbatim --
        and a ``None`` ``state.n_shell`` is refused loudly rather than silently
        re-derived. For every other scenario (``"fixed"``, ``"anchored_discrete"``,
        the collision path) it is re-derived from the realized mass via ``rint``
        -- byte-identical to the pre-Slice-G writer, so the deterministic paths
        are untouched. **Required** keyword (Slice-G review fix): an omitted
        scenario on a biphasic caller would silently revert the stored ``n_shell``
        to the mass-derived value, masking exactly the channel/reset drift the
        genuine state exists to expose.

    Raises
    ------
    ValueError
        If ``mass_scenario == "biphasic"`` and ``state.n_shell`` is ``None``.
    """
    if mass_scenario == "biphasic" and state.n_shell is None:
        raise ValueError(
            "write_ion_state_to_checkpoint_column: mass_scenario='biphasic' "
            "requires genuine n_shell state on the IonStepState (got None); "
            "silently falling back to the rint-from-mass derivation would mask "
            "channel/reset drift."
        )
    ckpt.positions_x[:, t_id] = state.x
    ckpt.positions_y[:, t_id] = state.y
    ckpt.positions_z[:, t_id] = state.z
    ckpt.velocities_x[:, t_id] = state.vx
    ckpt.velocities_y[:, t_id] = state.vy
    ckpt.velocities_z[:, t_id] = state.vz
    ckpt.mass_history_kg[:, t_id] = state.mass_kg
    ckpt.E_kin_eV[:, t_id] = state.E_kin_eV
    ckpt.E_pot_eV[:, t_id] = state.E_pot_eV
    ckpt.E_dissip_eV[:, t_id] = state.E_dissip_eV
    ckpt.E_mass_transfer_eV[:, t_id] = state.E_mass_transfer_eV
    ckpt.E_int_eV[:, t_id] = state.E_int_eV
    # Per-atom integer He-shell count. On the `biphasic` path it is genuine state
    # (channels advance n independently of the mass reset), so store it verbatim;
    # otherwise re-derive from the realized mass via the same rule the v5->v6 load
    # shim uses (so writer and shim agree): n = round((m/U - m_I+) / m_He).
    # Constant under `fixed`; the 21->14 staircase under `anchored_discrete`.
    if mass_scenario == "biphasic":
        ckpt.n_shell[:, t_id] = state.n_shell
    else:
        ckpt.n_shell[:, t_id] = np.rint(
            (state.mass_kg / U - MASS_I_ION_AMU) / MASS_HE_AMU
        )
    ckpt.number_of_collisions[:, t_id] = state.number_of_collisions
    ckpt.time_ps[t_id] = state.time_ps
