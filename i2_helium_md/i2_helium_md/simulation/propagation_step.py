"""One-timestep advance for the neutral propagation loop (pure version).

Public API:

* :class:`NeutralStepState` -- a frozen dataclass holding everything
  needed to step (positions, velocities, cumulative diagnostics).
* :func:`neutral_propagation_step` -- a pure function: takes a
  ``NeutralStepState`` plus the previous-step distance traveled
  per atom (for the Mode-3 collision sampler), returns a new
  ``NeutralStepState``. Does not mutate inputs.

The function ports the per-step body of the MATLAB ``while`` loop in
``vmi_sim_3d_neutral_propa_HeDFT_mimic.m`` (lines ~536-913) excluding
the ``attach_he`` mass-attachment branch (out of scope for the
neutral stage).

Step sequence:

1. Leapfrog integrate one ``dt`` -> candidate positions/velocities and
   the per-pair Morse potential ``E_pot_partner``.
2. Compute per-atom depth into droplet ``r1 - droplet_radius``.
3. Sample hard-sphere collision events using **Mode 3** (probability per
   step = ``prev_distance * sigma * rho_droplet``). Atoms outside the
   droplet (``depth >= 0``) and below the Landau cutoff ``E_min`` cannot
   collide. If ``prev_distance is None`` (first step), no collisions.
4. Apply collisions, replacing post-leapfrog velocities for colliders
   with elastically-scattered velocities. Track ``actual_dE = E0 - E1``
   per atom.
5. Compute energy diagnostics (E_kin, E_pot per atom, cumulative
   E_dissip, cumulative L_droplet) and assemble the new state.

The driver (``simulation/neutral.py``) is responsible for orchestrating
calls to this function, tracking the previous-step distance, and
optionally storing only every K-th state to a checkpoint.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import SimConfig
from ..physics.collisions import apply_collision, sample_collision_events
from ..physics.constants import EV, U
from ..physics.leapfrog import make_neutral_step
from ..physics.potentials import droplet_potential


# ===========================================================================
# State container
# ===========================================================================
@dataclass(frozen=True)
class NeutralStepState:
    """The minimum-sufficient state to advance the neutral propagation.

    All trajectory arrays have shape ``(2N,)`` where ``N`` is the
    number of molecules and the 2N layout convention is preserved.
    Per-atom static information (mass, droplet radius) is **not**
    carried in the state -- it doesn't change per step. The driver
    holds those arrays once and passes them in alongside the state.

    Attributes
    ----------
    x, y, z : np.ndarray, shape (2N,)
        Cartesian positions in Angstrom.
    vx, vy, vz : np.ndarray, shape (2N,)
        Velocities in Angstrom/ps.
    E_kin_eV : np.ndarray, shape (2N,)
        Kinetic energy in eV at this state's time.
    E_pot_eV : np.ndarray, shape (2N,)
        Per-atom potential energy in eV (droplet + half partner pair).
    E_dissip_eV : np.ndarray, shape (2N,)
        Cumulative energy dissipated per atom up to this state's time.
    L_droplet_eV_ps : np.ndarray, shape (2N,)
        Cumulative path length inside the droplet, per atom (Angstrom).
        Despite the field name, units are Angstrom (legacy MATLAB
        naming preserved for checkpoint compatibility).
    time_ps : float
        Time at this state, in picoseconds.
    """

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    E_kin_eV: np.ndarray
    E_pot_eV: np.ndarray
    E_dissip_eV: np.ndarray
    L_droplet_eV_ps: np.ndarray
    time_ps: float


# ===========================================================================
# Public API
# ===========================================================================
def neutral_propagation_step(
    state: NeutralStepState,
    *,
    cfg: SimConfig,
    mass_kg: np.ndarray,
    droplet_radii: np.ndarray,
    prev_distance_angstrom: np.ndarray | None,
    rng: np.random.Generator,
) -> NeutralStepState:
    """Advance the neutral propagation by one ``dt``. Pure function.

    Parameters
    ----------
    state : NeutralStepState
        Current state (read only; not mutated).
    cfg : SimConfig
        Simulation config.
    mass_kg : np.ndarray, shape (2N,)
        Per-atom mass (constant across the run for neutral).
    droplet_radii : np.ndarray, shape (2N,)
        Per-atom droplet radius in Angstrom (constant across the run).
    prev_distance_angstrom : np.ndarray of shape (2N,) or None
        Distance traveled per atom during the previous step (used by
        the Mode-3 collision sampler). Pass ``None`` for the very
        first step of a run -- collisions are then disabled because
        there is no previous distance to use.
    rng : np.random.Generator
        Reproducible RNG for collision sampling.

    Returns
    -------
    NeutralStepState
        The new state at ``state.time_ps + dt_neutral``.

    Raises
    ------
    ValueError
        If ``cfg.hard_sphere_collision_mode != 3``.
    """
    if cfg.hard_sphere_collision_mode != 3:
        raise ValueError(
            f"Only collision mode 3 is implemented; got "
            f"hard_sphere_collision_mode={cfg.hard_sphere_collision_mode}."
        )

    dt = cfg.dt_neutral

    # 1. Leapfrog integration.
    step_fn = make_neutral_step(cfg, mass_kg, droplet_radii)
    (x1, y1, z1), (vx1, vy1, vz1), E_pot_partner_per_pair = step_fn(
        (state.x, state.y, state.z),
        (state.vx, state.vy, state.vz),
        dt,
    )

    # 2. Depth into droplet.
    r1 = np.sqrt(x1 ** 2 + y1 ** 2 + z1 ** 2)
    depth = r1 - droplet_radii

    # 3. Pre-collision energy E0 (eV).
    v1_speed_sq = vx1 ** 2 + vy1 ** 2 + vz1 ** 2
    E0_eV = 0.5 * mass_kg * (v1_speed_sq * 100.0 ** 2) / EV

    # 4. Mode-3 collision sampling. Uses the *previous* step's distance.
    n_atoms = state.x.shape[0]
    if prev_distance_angstrom is None:
        b_collision = np.zeros(n_atoms, dtype=bool)
    else:
        b_collision = sample_collision_events(
            distance_travelled_angstrom=prev_distance_angstrom,
            depth_angstrom=depth,
            E0_eV=E0_eV,
            sigma_angstrom_sq=cfg.geometric_scattering_crosssection_I,
            E_min_eV=cfg.E_min_eV,
            rng=rng,
        )

    # 5. Apply collisions to flagged atoms.
    masses_amu = mass_kg / U
    vx_after, vy_after, vz_after, dE_eV = apply_collision(
        vx=vx1, vy=vy1, vz=vz1,
        masses_amu=masses_amu,
        b_collision=b_collision,
        scatter_mass_amu=cfg.scatter_mass_neutral_amu,
        neutral_scatter_angle_std_deg=cfg.neutral_scatter_angle_std_deg,
        rng=rng,
    )

    # 6. Energy diagnostics for the new state.
    v_post_sq = vx_after ** 2 + vy_after ** 2 + vz_after ** 2
    E_kin_new_eV = 0.5 * mass_kg * (v_post_sq * 100.0 ** 2) / EV

    E_droplet_eV = droplet_potential(
        r1 - droplet_radii,
        steepness=cfg.potential_steepness,
        binding_energy=cfg.binding_energy_I_atom_eV,
    )
    E_partner_per_atom = np.tile(E_pot_partner_per_pair, 2) / 2.0
    E_pot_new_eV = E_droplet_eV + E_partner_per_atom

    # 7. Cumulative bookkeeping.
    E_dissip_new = state.E_dissip_eV + dE_eV

    lx = x1 - state.x
    ly = y1 - state.y
    lz = z1 - state.z
    step_length = np.sqrt(lx ** 2 + ly ** 2 + lz ** 2)
    inside = (depth < 0).astype(float)
    L_droplet_new = state.L_droplet_eV_ps + inside * step_length

    return NeutralStepState(
        x=x1, y=y1, z=z1,
        vx=vx_after, vy=vy_after, vz=vz_after,
        E_kin_eV=E_kin_new_eV,
        E_pot_eV=E_pot_new_eV,
        E_dissip_eV=E_dissip_new,
        L_droplet_eV_ps=L_droplet_new,
        time_ps=state.time_ps + dt,
    )


# ===========================================================================
# Convenience helpers for checkpoint I/O
# ===========================================================================
def state_from_checkpoint_column(ckpt, t_id: int) -> NeutralStepState:
    """Extract a ``NeutralStepState`` from column ``t_id`` of a checkpoint.

    Used by the driver to bootstrap the inner loop after
    :func:`build_initial_state`. Copies the underlying arrays so the
    state is independent of the checkpoint's storage.
    """
    return NeutralStepState(
        x=ckpt.positions_x[:, t_id].copy(),
        y=ckpt.positions_y[:, t_id].copy(),
        z=ckpt.positions_z[:, t_id].copy(),
        vx=ckpt.velocities_x[:, t_id].copy(),
        vy=ckpt.velocities_y[:, t_id].copy(),
        vz=ckpt.velocities_z[:, t_id].copy(),
        E_kin_eV=ckpt.E_kin_eV[:, t_id].copy(),
        E_pot_eV=ckpt.E_pot_eV[:, t_id].copy(),
        E_dissip_eV=ckpt.E_dissip_eV[:, t_id].copy(),
        L_droplet_eV_ps=ckpt.L_droplet_eV_ps[:, t_id].copy(),
        time_ps=float(ckpt.time_ps[t_id]),
    )


def write_state_to_checkpoint_column(
    state: NeutralStepState,
    ckpt,
    t_id: int,
) -> None:
    """Write a ``NeutralStepState`` into column ``t_id`` of a checkpoint."""
    ckpt.positions_x[:, t_id] = state.x
    ckpt.positions_y[:, t_id] = state.y
    ckpt.positions_z[:, t_id] = state.z
    ckpt.velocities_x[:, t_id] = state.vx
    ckpt.velocities_y[:, t_id] = state.vy
    ckpt.velocities_z[:, t_id] = state.vz
    ckpt.E_kin_eV[:, t_id] = state.E_kin_eV
    ckpt.E_pot_eV[:, t_id] = state.E_pot_eV
    ckpt.E_dissip_eV[:, t_id] = state.E_dissip_eV
    ckpt.L_droplet_eV_ps[:, t_id] = state.L_droplet_eV_ps
    ckpt.time_ps[t_id] = state.time_ps
