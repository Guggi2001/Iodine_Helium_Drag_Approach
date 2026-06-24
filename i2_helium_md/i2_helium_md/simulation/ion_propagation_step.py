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
from ..physics.mass_jump import continuous_velocity_shed_components
from ..physics.shell_schedule import ShellSchedule
from ..physics.drag import REALIZED_FORMS
from ..physics.interactions import partner_interaction_ion
from ..physics.leapfrog import make_ion_step
from ..physics.potentials import droplet_potential


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
    number_of_collisions : np.ndarray, shape (2N,)
        Cumulative number of hard-sphere collisions per atom.
    time_ps : float
        Time at this state, in picoseconds.
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
    number_of_collisions: np.ndarray
    time_ps: float
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
        number_of_collisions=n_coll_new,
        time_ps=state.time_ps + dt,
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

    # Dissipated energy: amu*A^2/ps^2 -> eV via the baseline idiom (amu->kg via U,
    # A/ps->m/s via 100, J->eV via EV) -- the same conversion path as the
    # mass-attach defect below in ion_propagation_step, so E_dissip_eV stays in a
    # consistent eV with E_kin_eV / E_pot_eV.
    dE_dissip_eV = dE_dissip * U * (100.0 ** 2) / EV
    E_dissip_new = state.E_dissip_eV + dE_dissip_eV

    return IonStepState(
        x=x1, y=y1, z=z1,
        vx=vx1, vy=vy1, vz=vz1,
        mass_kg=state.mass_kg,                              # fixed mass (Tier 0)
        E_kin_eV=E_kin_new_eV,
        E_pot_eV=E_pot_new_eV,
        E_dissip_eV=E_dissip_new,
        E_mass_transfer_eV=state.E_mass_transfer_eV,  # stays 0 (no attach)
        number_of_collisions=state.number_of_collisions,       # stays 0 (no collisions)
        time_ps=state.time_ps + dt,
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
    velocities unchanged, drop one He from the (uniform) complex mass, book the
    per-atom co-moving-He kinetic energy into ``E_mass_transfer_eV``, and advance
    the pointer. **At most one shed per call** (the plan's <=1/step rule):
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

    vx_p, vy_p, vz_p, m_plus_amu, dE_amu = continuous_velocity_shed_components(
        state.vx, state.vy, state.vz, event.mass_before_amu, m_he_amu=m_he_amu,
    )
    # amu*A^2/ps^2 -> eV via the baseline idiom (amu->kg via U, A/ps->m/s via 100,
    # J->eV via EV) -- the same path baoab_propagation_step uses for dE_dissip, so
    # the booked defect stays in a consistent eV with E_kin/E_pot/E_dissip.
    dE_eV = dE_amu * U * (100.0 ** 2) / EV
    new_state = replace(
        state,
        vx=vx_p, vy=vy_p, vz=vz_p,
        mass_kg=np.full_like(state.mass_kg, m_plus_amu * U),
        E_mass_transfer_eV=state.E_mass_transfer_eV + dE_eV,
    )
    return new_state, next_shed_idx + 1


# ===========================================================================
# Internal helpers
# ===========================================================================
def _depth(x, y, z, droplet_radii):
    """Per-atom depth into the droplet ``r_atom - r_droplet`` (negative inside).

    Pure geometry, shared by the collision and drag paths (CLAUDE.md rule 1).
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return r - droplet_radii


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
    """Refuse to run the drag branch outside the deterministic fixed-mass envelope.

    The drag-branch analog of :func:`_check_scope`. ``_check_scope`` demands
    collision mode 3, which is irrelevant under drag; the drag path instead
    asserts the deterministic / fixed-mass / realised-form envelope
    (``DRAG_PORT_DESIGN_DECISIONS.md`` §6.4, form set widened by the METHOD_B
    §10 form phase) so an out-of-scope drag config fails at the driver, not
    deep in a half-implemented path.

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
    if cfg.mass_scenario not in ("fixed", "anchored_discrete"):
        unsupported.append(
            f"mass_scenario={cfg.mass_scenario!r} (drag branch supports 'fixed' and "
            "the Tier-1a 'anchored_discrete'; 'biphasic' is the Tier-2 generative "
            "mechanism, not yet wired)"
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
            "drag-branch ion propagation is deterministic fixed-mass only "
            "and does not support: "
            + ", ".join(unsupported)
            + ". Envelope = mass_scenario='fixed', noise_form='none', a "
            "realised drag_form (METHOD_B §10 form phase)."
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
        number_of_collisions=ckpt.number_of_collisions[:, t_id].copy(),
        time_ps=float(ckpt.time_ps[t_id]),
    )


def write_ion_state_to_checkpoint_column(
    state: IonStepState,
    ckpt,
    t_id: int,
) -> None:
    """Write an ``IonStepState`` into column ``t_id`` of an IonCheckpoint."""
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
    # Per-atom integer He-shell count, derived from the realized mass via
    # the same rule the v5->v6 load shim uses (so writer and shim agree):
    # n = round((m/U - m_I+) / m_He). Constant under `fixed`; the 21->14
    # staircase under `anchored_discrete`.
    ckpt.n_shell[:, t_id] = np.rint(
        (state.mass_kg / U - MASS_I_ION_AMU) / MASS_HE_AMU
    )
    ckpt.number_of_collisions[:, t_id] = state.number_of_collisions
    ckpt.time_ps[t_id] = state.time_ps
