"""Build the initial physical state of the neutral propagation.

This module is the integration point between the various sampling and
physics primitives (droplet sizes, radial positions, orientations,
potentials) and the `NeutralCheckpoint` data structure. It owns the
**2N array layout convention**: atom 1 of molecule i sits at index i,
and atom 2 sits at index ``num_molecules + i``.

The function `build_initial_state` produces a fully-allocated
`NeutralCheckpoint` with `num_steps` columns, where only column 0
is populated with physical data. The driver (`run_neutral_propagation`)
then writes columns 1..num_steps-1 by time-stepping.

Replaces the initialization block of
`vmi_sim_3d_neutral_propa_HeDFT_mimic.m` (lines ~284-483, excluding
the optional DFT-pre-fill block which is left unimplemented for now).
"""

from __future__ import annotations

import numpy as np

from ..config import SimConfig
from ..physics.constants import EV, HC, MASS_I_AMU, U, droplet_radius_bulk_angstrom
from ..physics.interactions import partner_interaction_neutral
from ..physics.potentials import droplet_potential
from ..sampling.droplet_sizes import sample_droplet_sizes
from ..sampling.orientations import sample_orientations
from ..sampling.radial_positions import sample_radial_positions
from .checkpoint import NeutralCheckpoint, _NEUTRAL_SCHEMA_VERSION


# ===========================================================================
# Public API
# ===========================================================================
def build_initial_state(
    cfg: SimConfig,
    *,
    num_steps: int,
    rng: np.random.Generator | None = None,
) -> NeutralCheckpoint:
    """Build the t=0 initial state and pre-allocate trajectory arrays.

    The returned `NeutralCheckpoint` has all per-step arrays sized
    `(2N, num_steps)` (or `(N, num_steps)` for `E_initial_eV` since
    that is per-molecule), with column 0 populated and the rest zero.
    The driver writes columns 1..num_steps-1 in place.

    Parameters
    ----------
    cfg : SimConfig
        Simulation configuration.
    num_steps : int
        Number of timesteps to allocate in the trajectory arrays
        (= number of stored columns; may be less than the integrator's
        internal step count if the driver downsamples).
    rng : np.random.Generator, optional
        Reproducible RNG. If None, built from `cfg.seed`.

    Returns
    -------
    NeutralCheckpoint
    """
    if rng is None:
        rng = np.random.default_rng(cfg.seed)
    if num_steps < 1:
        raise ValueError(f"num_steps must be >= 1, got {num_steps}")

    N = cfg.num_molecules

    # 1. Sample droplet sizes -> per-molecule droplet count.
    if cfg.use_single_droplet_size:
        droplet_counts = np.full(N, cfg.single_droplet_size, dtype=float)
    else:
        droplet_counts = sample_droplet_sizes(cfg, mode="post_pickup", rng=rng)

    # 2. Convert N -> droplet radius (Angstrom). Bulk-density formula.
    droplet_radii_per_molecule = droplet_radius_bulk_angstrom(droplet_counts)

    # 3. Radial position of molecule centre inside the droplet.
    r0 = sample_radial_positions(
        cfg,
        droplet_radii=droplet_radii_per_molecule,
        rng=rng,
    )

    # 4. Orientation angles + bond length.
    orient = sample_orientations(
        N,
        R0_GS_angstrom=cfg.R0_GS_angstrom,
        deltaR0_angstrom=cfg.deltaR0_angstrom,
        anisotropic=cfg.single_pulse,
        rng=rng,
    )

    # 5. Initial speed v0. MATLAB lines 92-115:
    #    E_initial = hc/lambda_pump * eV   (note: NB this is in J, not eV,
    #                                       despite the variable name)
    #    fwhm_v ~ FWHM of laser bandwidth -> velocity spread
    #    mean_v depends on whether the partner interaction is included
    #    (with partner: v from full E_initial; without: subtract E_diss)
    E_initial_J = (HC / cfg.lambda_pump_nm) * EV  # joules
    fwhm_E_eV = HC / cfg.lambda_pump_nm ** 2 * cfg.fwhm_lambda_nm
    # MATLAB:
    #   fwhm_v = (fwhm_E[eV] * eV[J/eV]) / sqrt(E_initial[J] * 127*u[kg]) / 100
    # Result in Å/ps (the /100 converts m/s -> Å/ps).
    mass_kg = MASS_I_AMU * U
    fwhm_v_a_per_ps = (
        0.5 * fwhm_E_eV * EV / np.sqrt(E_initial_J * mass_kg) / 100.0
    )
    if cfg.partner_interaction:
        mean_v_a_per_ps = np.sqrt(E_initial_J / mass_kg) / 100.0
    else:
        mean_v_a_per_ps = (
            np.sqrt((E_initial_J - cfg.E_diss_eV * EV) / mass_kg) / 100.0
        )
    if cfg.effusive_dynamics:
        fwhm_v_a_per_ps = 0.0

    v0 = rng.standard_normal(N) * fwhm_v_a_per_ps + mean_v_a_per_ps
    if cfg.single_pulse:
        # MATLAB: "if single_pulse, v0 = v0*0;" -- atoms start at rest
        v0 = np.zeros_like(v0)

    # 6. Assemble per-atom xyz and velocities (2N layout).
    #    Atom 1 sits at the molecule centre + (bond/2) * axis.
    #    Atom 2 sits at the molecule centre - (bond/2) * axis.
    #    MATLAB uses sin(delta+pi) = -sin(delta) for atom 2 in xy,
    #    and cos(delta+pi) = -cos(delta) for atom 2 in z.
    half_bond = orient.bond_length_angstrom / 2.0
    cb_sg = np.cos(orient.beta) * np.sin(orient.gamma)
    sb_sg = np.sin(orient.beta) * np.sin(orient.gamma)
    cg = np.cos(orient.gamma)
    ca_sd = np.cos(orient.alpha) * np.sin(orient.delta)
    sa_sd = np.sin(orient.alpha) * np.sin(orient.delta)
    cd = np.cos(orient.delta)

    x_centre = r0 * cb_sg
    y_centre = r0 * sb_sg
    z_centre = r0 * cg

    # Atom 1 (indices 0..N-1)
    x_atom1 = x_centre + ca_sd * half_bond
    y_atom1 = y_centre + sa_sd * half_bond
    z_atom1 = z_centre + cd * half_bond
    # Atom 2 (indices N..2N-1) -- offset by sin(delta+pi) = -sin(delta) and
    # cos(delta+pi) = -cos(delta) on each component, so it sits opposite.
    # MATLAB writes sin(delta+pi) explicitly, which equals -sin(delta).
    x_atom2 = x_centre - ca_sd * half_bond
    y_atom2 = y_centre - sa_sd * half_bond
    z_atom2 = z_centre - cd * half_bond

    # Velocities (along axis, with shared speed v0)
    vx_atom1 = v0 * ca_sd
    vy_atom1 = v0 * sa_sd
    vz_atom1 = v0 * cd
    vx_atom2 = -vx_atom1
    vy_atom2 = -vy_atom1
    vz_atom2 = -vz_atom1

    # 7. Allocate trajectory arrays and fill column 0.
    positions_x = np.zeros((2 * N, num_steps))
    positions_y = np.zeros((2 * N, num_steps))
    positions_z = np.zeros((2 * N, num_steps))
    velocities_x = np.zeros((2 * N, num_steps))
    velocities_y = np.zeros((2 * N, num_steps))
    velocities_z = np.zeros((2 * N, num_steps))

    positions_x[:N, 0] = x_atom1
    positions_x[N:, 0] = x_atom2
    positions_y[:N, 0] = y_atom1
    positions_y[N:, 0] = y_atom2
    positions_z[:N, 0] = z_atom1
    positions_z[N:, 0] = z_atom2

    velocities_x[:N, 0] = vx_atom1
    velocities_x[N:, 0] = vx_atom2
    velocities_y[:N, 0] = vy_atom1
    velocities_y[N:, 0] = vy_atom2
    velocities_z[:N, 0] = vz_atom1
    velocities_z[N:, 0] = vz_atom2

    # 8. Per-atom static arrays (replicated for atom 1 + atom 2).
    droplet_radii_2N = np.tile(droplet_radii_per_molecule, 2)
    mass_kg_2N = np.full(2 * N, mass_kg)

    # 9. Energy at t=0.
    # Per-atom kinetic energy (eV) = 0.5 * m * v^2 with units conversion.
    v_sq_2N = (
        velocities_x[:, 0] ** 2 + velocities_y[:, 0] ** 2 + velocities_z[:, 0] ** 2
    )
    # v in Å/ps -> m/s factor 100
    E_kin_eV_t0 = 0.5 * mass_kg_2N * (v_sq_2N * 100.0 ** 2) / EV

    # Per-atom droplet potential (eV).
    r_atom_2N = np.sqrt(
        positions_x[:, 0] ** 2 + positions_y[:, 0] ** 2 + positions_z[:, 0] ** 2
    )
    E_droplet_eV_t0 = droplet_potential(
        r_atom_2N - droplet_radii_2N,
        steepness=cfg.potential_steepness,
        binding_energy=cfg.binding_energy_I_atom_eV,
    )

    # Per-pair Morse potential energy (eV). Mirrors the per-step
    # convention in `propagation_step.py`: each atom carries half of
    # its pair's Morse energy, so summing per-atom E_pot across a
    # molecule recovers droplet[atom1] + droplet[atom2] + Morse[pair].
    #
    # Note: the legacy MATLAB code OMITS this Morse term at t=0 (line
    # 476 of vmi_sim_3d_neutral_propa_HeDFT_mimic.m) but INCLUDES it
    # for all subsequent timesteps (line 885). That makes E_pot
    # discontinuous between t=0 and t=1 by the full pair energy --
    # in our test case ~3 eV per molecule for R=9 A. We treat this as
    # a legacy bug and include the Morse term at t=0 here, per project
    # principle #10 ("don't preserve legacy approximations").
    _ax, _ay, _az, E_pot_partner_per_pair = partner_interaction_neutral(
        positions_x[:, 0], positions_y[:, 0], positions_z[:, 0],
        mass_kg_2N, cfg,
    )
    E_partner_per_atom_t0 = np.tile(E_pot_partner_per_pair, 2) / 2.0
    E_pot_eV_t0 = E_droplet_eV_t0 + E_partner_per_atom_t0

    # E_initial: per-molecule photon energy in eV
    # (The "per-molecule" part of the laser deposit: hc/lambda in eV.)
    E_initial_eV_per_mol = np.full(N, HC / cfg.lambda_pump_nm)

    E_kin_eV = np.zeros((2 * N, num_steps))
    E_kin_eV[:, 0] = E_kin_eV_t0
    E_pot_eV = np.zeros((2 * N, num_steps))
    E_pot_eV[:, 0] = E_pot_eV_t0
    E_dissip_eV = np.zeros((2 * N, num_steps))
    L_droplet_eV_ps = np.zeros((2 * N, num_steps))

    time_ps = np.zeros(num_steps)

    return NeutralCheckpoint(
        num_molecules=N,
        time_ps=time_ps,
        positions_x=positions_x,
        positions_y=positions_y,
        positions_z=positions_z,
        velocities_x=velocities_x,
        velocities_y=velocities_y,
        velocities_z=velocities_z,
        mass_kg=mass_kg_2N,
        droplet_radii=droplet_radii_2N,
        r0=r0,
        E_kin_eV=E_kin_eV,
        E_pot_eV=E_pot_eV,
        E_initial_eV=E_initial_eV_per_mol,
        E_dissip_eV=E_dissip_eV,
        L_droplet_eV_ps=L_droplet_eV_ps,
        schema_version=_NEUTRAL_SCHEMA_VERSION,
    )
