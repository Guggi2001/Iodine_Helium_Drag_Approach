"""Build the initial physical state of the ion propagation.

This module is the integration point between the neutral checkpoint
(end-state of the neutral stage) and the ``IonCheckpoint`` data
structure. It owns the **2N array layout convention** (atom 1 of
molecule i at index i, atom 2 at index ``num_molecules + i``).

The function `build_initial_ion_state` produces a fully-allocated
`IonCheckpoint` with `num_steps_ion` columns, where only column 0
is populated with physical data. The driver
(``run_ion_propagation``) writes columns 1..num_steps_ion-1 by
time-stepping.

Replaces the initialization block of ``vmi_sim_3d_ion_propa.m``
(lines ~127-294) for the **single-pulse, no-effusive, no-additional-
charges** scope of the project.

Out of scope
------------
The legacy MATLAB code supports many features that are unused in our
two production input scripts (``single_pulse_N2000.m`` and
``single_pulse_droplet_distribution.m``). When the cfg has any of the
following set, this builder raises ``NotImplementedError`` rather
than silently producing wrong physics:

- ``effusive_dynamics = True``
- ``single_charge_ionization_allowed = True``
- ``additional_droplet_charges > 0``
- ``highly_charged_iodine = True``

The check fires here (at build time) so an unsupported run fails
loudly with a clear message before any expensive ion stepping.

Bug fixes vs. legacy MATLAB
---------------------------
The legacy ``vmi_sim_3d_ion_propa.m`` has two **t=0 bookkeeping bugs**
in its E_kin/E_pot recording (lines ~289-291):

1. ``E_kin_ion(:,1) = mass_i.*(vx² + vy²)²/2/eV`` -- **MISSING vz**
   AND ``(...)²`` squares the kinetic-energy expression so the formula
   computes ``m * v⁴ / 2 / eV`` instead of ``m * v² / 2 / eV``.

2. ``E_pot_ion(:,1) = droplet_potential(sqrt(x² + y²) - R)`` --
   **MISSING vz/z** in the radial coordinate, AND missing the partner
   Coulomb term that subsequent steps DO include.

These bugs are silent because ``E_pot_ion`` is only ever read for
diagnostic plotting at the end of the run (same pattern as the t=0
E_pot bug we found in the neutral stage in Step 10). We **fix both
in this Python port** per project principle #10 and add regression
tests that catch them.
"""

from __future__ import annotations

import numpy as np

from ..config import SimConfig
from ..physics.constants import EV, U
from ..physics.interactions import partner_interaction_ion
from ..physics.potentials import droplet_potential
from .checkpoint import IonCheckpoint, NeutralCheckpoint, _ION_SCHEMA_VERSION


# ===========================================================================
# Public API
# ===========================================================================
def build_initial_ion_state(
    cfg: SimConfig,
    neutral_ckpt: NeutralCheckpoint,
    *,
    num_steps_ion: int,
    start_id: int = -1,
    rng: np.random.Generator | None = None,
) -> IonCheckpoint:
    """Build the ion-stage t=0 state from a neutral checkpoint.

    The neutral stage runs from t=0 to t=t_neutral_max. At time
    ``start_id`` (column index in ``neutral_ckpt``) we "switch on"
    ionization: every atom gets charge +1, the ion-stage potentials
    take over, and a fresh ion trajectory starts. This builder
    produces the t=0 state of that ion run, with all later columns
    zero-allocated for the driver to fill.

    Parameters
    ----------
    cfg : SimConfig
        Simulation configuration.
    neutral_ckpt : NeutralCheckpoint
        End-state of the neutral stage. We read positions, velocities,
        masses, and droplet radii from column ``start_id``.
    num_steps_ion : int
        Number of timesteps to allocate in the ion trajectory arrays
        (= number of stored columns; may be less than the integrator's
        internal step count if the driver downsamples).
    start_id : int, optional
        Which column of ``neutral_ckpt`` to use as the ion start state.
        Default ``-1`` (last column = end of neutral stage), matching
        the production single-pulse use case. Other values exist for
        pump-probe experiments which are out of scope here.
    rng : np.random.Generator, optional
        Reproducible RNG. Currently unused (single-charge ionization
        is out of scope, and the ion stage has no rng-driven init for
        our scope). Reserved for future expansion.

    Returns
    -------
    IonCheckpoint
        With column 0 populated and columns 1..num_steps_ion-1 zero.

    Raises
    ------
    NotImplementedError
        If cfg requests unsupported features (see module docstring).
    ValueError
        If shapes don't match or num_steps_ion < 1.
    """
    # 1. Scope checks: refuse to run with unsupported features.
    _check_scope(cfg)

    if num_steps_ion < 1:
        raise ValueError(f"num_steps_ion must be >= 1, got {num_steps_ion}")

    # 2. Validate neutral_ckpt against cfg (shape consistency).
    if neutral_ckpt.num_molecules != cfg.num_molecules:
        raise ValueError(
            f"neutral_ckpt has num_molecules={neutral_ckpt.num_molecules} "
            f"but cfg has num_molecules={cfg.num_molecules}"
        )

    N = cfg.num_molecules
    two_N = 2 * N
    T = num_steps_ion

    # 3. Pick the start column.
    n_neutral_steps = neutral_ckpt.positions_x.shape[1]
    if start_id < 0:
        start_id = n_neutral_steps + start_id   # -1 -> last
    if not (0 <= start_id < n_neutral_steps):
        raise ValueError(
            f"start_id={start_id} out of range; "
            f"neutral_ckpt has {n_neutral_steps} columns"
        )

    x0 = neutral_ckpt.positions_x[:, start_id].copy()
    y0 = neutral_ckpt.positions_y[:, start_id].copy()
    z0 = neutral_ckpt.positions_z[:, start_id].copy()
    vx0 = neutral_ckpt.velocities_x[:, start_id].copy()
    vy0 = neutral_ckpt.velocities_y[:, start_id].copy()
    vz0 = neutral_ckpt.velocities_z[:, start_id].copy()

    # 4. Inherit static per-atom data from the neutral checkpoint.
    mass_kg_initial = neutral_ckpt.mass_kg.copy()
    droplet_radii_angstrom = neutral_ckpt.droplet_radii.copy()

    # 5. Charges: all +1, since single_charge_ionization_allowed=False
    #    (verified above by _check_scope).
    charge = np.ones(two_N, dtype=float)

    # 6. Compute t=0 energies. These are the "fixed" formulas (see
    #    module docstring): include z component, and include the
    #    partner Coulomb energy.
    E_kin_t0 = _compute_E_kin_per_atom(mass_kg_initial, vx0, vy0, vz0)

    E_drop_t0 = _compute_E_pot_droplet_per_atom(
        x0, y0, z0, droplet_radii_angstrom, cfg,
    )
    _, _, _, E_partner_per_atom_t0 = partner_interaction_ion(
        x0, y0, z0, mass_kg_initial, charge, cfg,
    )
    # partner_interaction_ion already returns per-atom half-pair energy
    E_pot_t0 = E_drop_t0 + E_partner_per_atom_t0

    # 7. Allocate trajectory arrays and fill column 0.
    positions_x = np.zeros((two_N, T))
    positions_y = np.zeros((two_N, T))
    positions_z = np.zeros((two_N, T))
    velocities_x = np.zeros((two_N, T))
    velocities_y = np.zeros((two_N, T))
    velocities_z = np.zeros((two_N, T))

    positions_x[:, 0] = x0
    positions_y[:, 0] = y0
    positions_z[:, 0] = z0
    velocities_x[:, 0] = vx0
    velocities_y[:, 0] = vy0
    velocities_z[:, 0] = vz0

    E_kin_eV = np.zeros((two_N, T))
    E_pot_eV = np.zeros((two_N, T))
    E_dissip_eV = np.zeros((two_N, T))
    E_mass_attach_defect_eV = np.zeros((two_N, T))
    relative_loss_per_ps = np.zeros((two_N, T))
    number_of_collisions = np.zeros((two_N, T), dtype=int)
    mass_history_kg = np.zeros((two_N, T))

    E_kin_eV[:, 0] = E_kin_t0
    E_pot_eV[:, 0] = E_pot_t0
    # E_dissip, E_mass_attach_defect, relative_loss, n_collisions all
    # start at 0 (already zeros).
    mass_history_kg[:, 0] = mass_kg_initial

    # 8. Static finals -- placeholder, the driver fills these at end.
    positions_final_x = np.zeros(two_N)
    positions_final_y = np.zeros(two_N)
    positions_final_z = np.zeros(two_N)
    velocities_final_x = np.zeros(two_N)
    velocities_final_y = np.zeros(two_N)
    velocities_final_z = np.zeros(two_N)
    mass_final_kg = mass_kg_initial.copy()  # driver may overwrite after attachment
    b_ion_outside = np.zeros(N, dtype=bool)

    # 9. Time axis: dt_ion * t_index, filled in by driver. Initialize to zeros.
    time_ps = np.zeros(T)

    return IonCheckpoint(
        num_molecules=N,
        time_ps=time_ps,
        positions_x=positions_x,
        positions_y=positions_y,
        positions_z=positions_z,
        velocities_x=velocities_x,
        velocities_y=velocities_y,
        velocities_z=velocities_z,
        positions_final_x=positions_final_x,
        positions_final_y=positions_final_y,
        positions_final_z=positions_final_z,
        velocities_final_x=velocities_final_x,
        velocities_final_y=velocities_final_y,
        velocities_final_z=velocities_final_z,
        mass_kg=mass_kg_initial,
        mass_final_kg=mass_final_kg,
        mass_history_kg=mass_history_kg,
        droplet_radii_angstrom=droplet_radii_angstrom,
        E_kin_eV=E_kin_eV,
        E_pot_eV=E_pot_eV,
        E_dissip_eV=E_dissip_eV,
        E_mass_attach_defect_eV=E_mass_attach_defect_eV,
        b_ion_outside=b_ion_outside,
        relative_loss_per_ps=relative_loss_per_ps,
        number_of_collisions=number_of_collisions,
        schema_version=_ION_SCHEMA_VERSION,
    )


# ===========================================================================
# Internal helpers
# ===========================================================================
def _check_scope(cfg: SimConfig) -> None:
    """Refuse to build an ion state for cfg flags we don't support."""
    unsupported = []
    if cfg.effusive_dynamics:
        unsupported.append("effusive_dynamics")
    if cfg.single_charge_ionization_allowed:
        unsupported.append("single_charge_ionization_allowed")
    if cfg.additional_droplet_charges > 0:
        unsupported.append(
            f"additional_droplet_charges={cfg.additional_droplet_charges} (must be 0)"
        )
    if cfg.highly_charged_iodine:
        unsupported.append("highly_charged_iodine")

    if unsupported:
        raise NotImplementedError(
            "ion stage does not yet support: "
            + ", ".join(unsupported)
            + ". The two production input scripts (single_pulse_N2000.m and "
            "single_pulse_droplet_distribution.m) leave all of these at "
            "their default disabled values. To re-enable, see the design "
            "discussion in migration_log.md."
        )


def _compute_E_kin_per_atom(
    mass_kg: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
) -> np.ndarray:
    """E_kin per atom in eV from velocities in Å/ps and mass in kg.

    Standard ``½ m v²`` with v in m/s (v_AA_ps × 100) and result in eV.
    Encapsulated as a helper so the conversion factor is in one place.
    """
    v_sq_m2_per_s2 = (vx ** 2 + vy ** 2 + vz ** 2) * (100.0 ** 2)
    return 0.5 * mass_kg * v_sq_m2_per_s2 / EV


def _compute_E_pot_droplet_per_atom(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    droplet_radii: np.ndarray,
    cfg: SimConfig,
) -> np.ndarray:
    """Ion-droplet potential per atom in eV.

    Uses the SAME ``droplet_potential`` shape function as neutral, but
    with ``cfg.binding_energy_I_ion_eV`` instead of the atom-droplet
    binding energy.
    """
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return droplet_potential(
        r - droplet_radii,
        steepness=cfg.potential_steepness,
        binding_energy=cfg.binding_energy_I_ion_eV,
    )
