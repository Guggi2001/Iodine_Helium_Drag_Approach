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

from dataclasses import dataclass

import numpy as np

from ..config import SimConfig
from ..physics.collisions import (
    apply_collision,
    sample_collision_events,
    velocity_dependent_cross_section,
)
from ..physics.constants import EV, U
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
    number_of_collisions : np.ndarray, shape (2N,)
        Cumulative number of hard-sphere collisions per atom.
    time_ps : float
        Time at this state, in picoseconds.
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
    number_of_collisions: np.ndarray
    time_ps: float


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
    r1 = np.sqrt(x1 ** 2 + y1 ** 2 + z1 ** 2)
    depth = r1 - droplet_radii

    # 3. Pre-collision speed and energy. We compute these AFTER the
    #    leapfrog step so the collision uses the energy at the new
    #    position (matches MATLAB line 384: E0 = (v1*100)^2 * m / 2 / eV).
    v1_speed_sq = vx1 ** 2 + vy1 ** 2 + vz1 ** 2
    v1_speed = np.sqrt(v1_speed_sq)
    E0_eV = 0.5 * state.mass_kg * (v1_speed_sq * 100.0 ** 2) / EV

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

    # 6. Apply collisions to flagged atoms.
    masses_amu = state.mass_kg / U
    vx_after, vy_after, vz_after, dE_eV = apply_collision(
        vx=vx1, vy=vy1, vz=vz1,
        masses_amu=masses_amu,
        b_collision=b_collision,
        scatter_mass_amu=cfg.scatter_mass_ion_amu,
        neutral_scatter_angle_std_deg=cfg.ion_scatter_angle_std_deg,
        rng=rng,
    )

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
    E_kin_new_eV = 0.5 * new_mass_kg * (v_post_sq * 100.0 ** 2) / EV

    # E_pot_droplet at NEW positions. Uses ion binding energy.
    E_droplet_eV = droplet_potential(
        depth,
        steepness=cfg.potential_steepness,
        binding_energy=cfg.binding_energy_I_ion_eV,
    )
    # Split per-pair Coulomb energy half-and-half between the two atoms
    # of each molecule (matches the neutral propagation_step convention).
    E_partner_per_atom = np.tile(E_pot_coulomb_per_pair, 2) / 2.0
    E_pot_new_eV = E_droplet_eV + E_partner_per_atom

    # 9. Cumulative bookkeeping.
    E_dissip_new = state.E_dissip_eV + dE_eV
    n_coll_new = state.number_of_collisions + b_collision.astype(state.number_of_collisions.dtype)

    return IonStepState(
        x=x1, y=y1, z=z1,
        vx=vx_after, vy=vy_after, vz=vz_after,
        mass_kg=new_mass_kg,
        E_kin_eV=E_kin_new_eV,
        E_pot_eV=E_pot_new_eV,
        E_dissip_eV=E_dissip_new,
        number_of_collisions=n_coll_new,
        time_ps=state.time_ps + dt,
    )


# ===========================================================================
# Internal helpers
# ===========================================================================
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
    ckpt.number_of_collisions[:, t_id] = state.number_of_collisions
    ckpt.time_ps[t_id] = state.time_ps
