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
from ..physics.constants import EV, MASS_HE_AMU, MASS_I_ION_AMU, U
from ..physics.dissociation_ladder import resolve_ladder
from ..physics.helium_density import rho_he_ratio
from ..physics.interactions import partner_interaction_ion
from ..physics.internal_energy_budget import e_int_onset_eV, sigma_partition_factor
from ..physics.potentials import droplet_potential
from ..physics.shell_schedule import ANCHOR_N_START, complex_mass_amu
from ..physics.solvation_cooling import e_bind_pair_eV
from ..sampling.ce_channels import (
    CE_CHANNEL_NONE,
    CE_CHANNEL_Q2,
    CE_CHANNEL_Q3,
    CE_CHANNEL_Q3_PARTNER,
    CE_CHANNEL_SINGLE,
    CE_CHANNEL_STREAM_KEY,
    sample_ce_channels,
)
from .checkpoint import (
    IonCheckpoint,
    NeutralCheckpoint,
    _ION_SCHEMA_VERSION,
    stage_stream_rng,
)


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
    if cfg.drag_coefficients is not None and cfg.mass_scenario == "fixed":
        # Drag fixed-scenario integrates at the drag law's extraction mass
        # m_eff, NOT the inherited neutral ~127 amu (bare I+). The linear_cubic
        # law was fit under m_eff (~203 amu, ~19 He), and
        # DRAG_PORT_DESIGN_DECISIONS.md §6.5 makes `fixed` the *only*
        # self-consistent pairing -- it must run at m_eff or the calibrated law
        # is applied at the wrong inertia. Uniform fill over all 2N ions. We
        # read `mass_initial_amu` (the field whose purpose is "ion-stage initial
        # mass", §2.8); at Tier 0 it equals `m_eff_amu`, but reading it keeps the
        # field semantics clean for Tier 1 (where the two diverge). Do NOT
        # "restore" the inherited neutral mass here -- the override is the fix.
        mass_kg_initial = np.full(two_N, cfg.mass_initial_amu * U)
    elif cfg.mass_scenario == "biphasic" or (
        cfg.drag_coefficients is not None and cfg.mass_scenario == "anchored_discrete"
    ):
        # Tier-1a anchored_discrete and Tier-2 generative biphasic both start at
        # the full first shell n0 = 21 (the 21->19->14 validation target): the
        # anchored schedule sheds it down to n=14, the generative pickup /
        # evaporation channels evolve n from there. NOT m_eff (= the n=19
        # mid-window mass) and NOT the inherited bare-I+ neutral mass -- the
        # onset shell count is the physical start. Biphasic is gated on the
        # scenario alone: it is a drag-path scenario by contract, enforced
        # upstream by config.check_biphasic_config + ion._check_scope_ion_driver
        # (a bundle-less biphasic cfg never reaches a propagation loop); the
        # biphasic column-0 seeds (S2 onset + E_pot binding fold) live below
        # with the other column-0 physics.
        if (
            cfg.mass_scenario == "biphasic"
            and cfg.initial_shell_model == "density_tied"
        ):
            # Slice T5 (plan §I.11; H.3b revived): dress each ion's t0 shell
            # by the local He availability at its birth position,
            # n_0i = round(n* * rho_hat(d_birth,i)), through the SAME
            # erf-complement surface the drag/pickup/cooling gates share --
            # zero new free parameters. A first-order occupancy statement:
            # under-dressed ions may re-fill via the live pickup channel.
            # density_tied is biphasic-only (config.check_initial_shell_config),
            # so the anchored_discrete path can never reach this branch.
            from .ion import drag_gate_steepness  # local: ion.py imports us

            depth_birth_angstrom = (
                np.sqrt(x0 ** 2 + y0 ** 2 + z0 ** 2)
                - neutral_ckpt.droplet_radii
            )
            n0_initial = np.rint(
                ANCHOR_N_START
                * rho_he_ratio(
                    depth_birth_angstrom, steepness=drag_gate_steepness(cfg)
                )
            )
            mass_kg_initial = complex_mass_amu(n0_initial) * U
        else:
            n0_initial = ANCHOR_N_START
            mass_kg_initial = np.full(
                two_N, complex_mass_amu(ANCHOR_N_START) * U
            )
    else:
        mass_kg_initial = neutral_ckpt.mass_kg.copy()
    droplet_radii_angstrom = neutral_ckpt.droplet_radii.copy()

    # 5. Charges: all +1, since single_charge_ionization_allowed=False
    #    (verified above by _check_scope).
    charge = np.ones(two_N, dtype=float)

    # 5b. Tier-2 (C) design (B): the per-molecule CE channel draw, on its own
    #     dedicated stream (SeedSequence((seed, CE_CHANNEL_STREAM_KEY))) --
    #     the passed driver rng (the ion-stage stream) is NEVER consumed here,
    #     so channels-off is byte-identical and channels-on leaves every
    #     pre-existing stream untouched (design §5; the s(n) S1 precedent).
    ce_draw = None
    if cfg.ce_channel_mode == "sampled":
        ce_draw = sample_ce_channels(
            N,
            stage_stream_rng(cfg.seed, CE_CHANNEL_STREAM_KEY),
            weights=cfg.ce_channel_weights,
            fraction_f=cfg.ce_fraction_f,
            sigma_eV=cfg.ce_channel_sigma_eV,
            single_ker_eV=cfg.ce_single_ker_eV,
        )
    pair_scale = None if ce_draw is None else ce_draw.pair_scale

    # 6. Compute t=0 energies. These are the "fixed" formulas (see
    #    module docstring): include z component, and include the
    #    partner Coulomb energy (per-pair CE-scaled under "sampled").
    E_kin_t0 = _compute_E_kin_per_atom(mass_kg_initial, vx0, vy0, vz0)

    E_drop_t0 = _compute_E_pot_droplet_per_atom(
        x0, y0, z0, droplet_radii_angstrom, cfg,
    )
    _, _, _, E_partner_per_atom_t0 = partner_interaction_ion(
        x0, y0, z0, mass_kg_initial, charge, cfg, pair_scale=pair_scale,
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
    E_mass_transfer_eV = np.zeros((two_N, T))
    # Tier-2 internal-energy reservoir (schema v7). All-zero at allocation;
    # it stays zero on the fixed / anchored_discrete paths and is filled by
    # the biphasic generative driver via the S1/S2/K1/K2 budget.
    E_int_eV = np.zeros((two_N, T))
    n_shell = np.zeros((two_N, T))
    relative_loss_per_ps = np.zeros((two_N, T))
    number_of_collisions = np.zeros((two_N, T), dtype=int)
    mass_history_kg = np.zeros((two_N, T))

    # Legacy MATLAB per-step temperature diagnostic. NaN means "no
    # collision in this stored step" (matches the "blank rows" we get
    # in MATLAB when an outer step has no colliders -- diagnostic_array
    # in vmi_sim_3d_ion_propa.m:683 is only appended when collisions
    # occurred). The driver overwrites rows where collisions happened.
    temperature_diagnostic = np.full((T, 3), np.nan, dtype=float)

    E_kin_eV[:, 0] = E_kin_t0
    E_pot_eV[:, 0] = E_pot_t0
    # E_dissip, E_mass_transfer, relative_loss, n_collisions all
    # start at 0 (already zeros).
    mass_history_kg[:, 0] = mass_kg_initial
    # Per-atom He-shell count at t=0, from the initial mass via the same
    # rule the driver writer and the v5->v6 load shim use.
    n_shell[:, 0] = np.rint(
        (mass_kg_initial / U - MASS_I_ION_AMU) / MASS_HE_AMU
    )
    # Biphasic column-0 seeds (grouped; Slice-G review fix -- one block owns the
    # scenario's t0 physics beyond the shared mass init above):
    # * E_pot binding fold: fold E_bind^pair(n) = -Sigma(n) into E_pot so the
    #   5-term invariant closes when the pickup/evaporation channels shift the
    #   rung ladder (a shed adds +D_0(n), a pickup removes -D_0(n+1) -- exactly
    #   the MASS §6 "E_pot += D_0 / -= D_0" bookings). Seeded at the t0 shell
    #   n_0 (21 under the `full` arm; per-ion dressed under T5 `density_tied`)
    #   so the offset is consistent from column 0; the electrostriction marginal stays
    #   the untracked "A8 -> bath" collective term. The driver re-applies the
    #   fold each step from the genuine n_shell state (solvation_cooling.
    #   e_bind_pair_eV; MASS §6).
    # * S2 onset: deposit E_int(0) = f_int * E_avail^ion once at t0 (MASS §6 S2).
    #   Only S1/K1/K2 modify E_int thereafter. f_int is required non-None under
    #   biphasic (config.check_biphasic_config), so it is set here by construction.
    if cfg.mass_scenario == "biphasic":
        # Slice T2 (§I.10): the t0 fold rides the same resolved ladder as the
        # driver/stages (fail-loud on a tabulated selector without a table).
        # Slice T5: the fold is seeded at the ion's actual t0 shell --
        # n0_initial is the scalar ANCHOR_N_START under the `full` arm
        # (byte-identical to the delivered call) or the per-ion dressed
        # array under `density_tied` (e_bind_pair_eV is vectorised).
        E_pot_eV[:, 0] += e_bind_pair_eV(
            n0_initial,
            picture=cfg.ladder_electronic_picture,
            kappa=cfg.ladder_steepness,
            ladder=resolve_ladder(
                cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV
            ),
        )
        # S2 onset E_int(0) = f_int * E_avail, then the Slice-T6 p-law dressing
        # coupling (Sigma(n0)/Sigma(n*))^p. Under the byte-inert "constant" arm
        # the factor is exactly 1.0 (the delivered per-ion onset, unchanged);
        # under "sigma_proportional" the onset rides the same per-ion t0 shell
        # n0_initial and resolved ladder as the E_pot binding fold above, so
        # under-dressed births get less onset. Inert at n0 = n* regardless of
        # law (factor 1) -- a no-op without the T5 dressing axis.
        # Under the (C) mixture (ce_channel_mode="sampled") the scalar
        # budget retires (design §3.1/§3.2): the per-ion onset becomes
        # f_int,c(m) * E_m with the T6 p-law factor preserved verbatim; the
        # emulated q3_partner rides the Q3 coupling.
        if ce_draw is not None:
            fc = cfg.ce_internal_energy_partition_fractions
            f_int_by_code = np.zeros(4, dtype=float)
            f_int_by_code[CE_CHANNEL_SINGLE] = fc[0]
            f_int_by_code[CE_CHANNEL_Q2] = fc[1]
            f_int_by_code[CE_CHANNEL_Q3] = fc[2]
            f_int_by_code[CE_CHANNEL_Q3_PARTNER] = fc[2]
            onset = e_int_onset_eV(
                f_int=f_int_by_code[ce_draw.channel],
                e_avail_eV=ce_draw.E_m_eV,
            )
        else:
            onset = e_int_onset_eV(
                f_int=cfg.internal_energy_partition_fraction,
                e_avail_eV=cfg.coulomb_available_eV,
            )
        E_int_eV[:, 0] = onset * sigma_partition_factor(
            n0_initial,
            law=cfg.internal_energy_partition_law,
            picture=cfg.ladder_electronic_picture,
            kappa=cfg.ladder_steepness,
            ladder=resolve_ladder(
                cfg.dissociation_ladder, cfg.tabulated_ladder_rungs_eV
            ),
        )

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

    # 10. Tier-2 (C) per-ion fields (schema v8, OQ-E scope): the channel
    #     assignment / KER stamp from the draw (or the exact no-channel
    #     sentinels), and a zero strip counter the drivers accumulate into.
    if ce_draw is not None:
        ce_channel = ce_draw.channel.copy()
        ce_E_m_eV = ce_draw.E_m_eV.copy()
    else:
        ce_channel = np.full(two_N, CE_CHANNEL_NONE, dtype=int)
        ce_E_m_eV = np.full(two_N, np.nan, dtype=float)
    ce_strip_count = np.zeros(two_N, dtype=int)

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
        E_mass_transfer_eV=E_mass_transfer_eV,
        E_int_eV=E_int_eV,
        n_shell=n_shell,
        b_ion_outside=b_ion_outside,
        relative_loss_per_ps=relative_loss_per_ps,
        number_of_collisions=number_of_collisions,
        temperature_diagnostic=temperature_diagnostic,
        ce_channel=ce_channel,
        ce_E_m_eV=ce_E_m_eV,
        ce_strip_count=ce_strip_count,
        mass_scenario=cfg.mass_scenario,
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
