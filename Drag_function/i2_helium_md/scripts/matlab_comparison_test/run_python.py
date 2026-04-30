"""Run a deterministic neutral-propagation comparison run.

Setup:
    1 I2 molecule, atoms at rest, bond aligned along x-axis, 9 A apart,
    centered in a single He droplet of fixed radius. NO collisions
    (cross section forced to zero). Pure leapfrog dynamics over the
    Morse pair potential + droplet solvation potential.

Outputs:
    `python_trajectory.csv` with columns
        t_ps, x_atom1, y_atom1, z_atom1, x_atom2, y_atom2, z_atom2,
        vx_atom1, vy_atom1, vz_atom1, vx_atom2, vy_atom2, vz_atom2,
        E_kin_eV, E_pot_eV, E_total_eV
    one row per timestep.

Run with:
    python run_python.py

This script bypasses `build_initial_state` to use a fixed initial
condition (no RNG dependency on the sampling layer). It also bypasses
`run_neutral_propagation` because we want to disable collisions
deterministically -- which the driver currently doesn't support.
We call the leapfrog directly, just like MATLAB would.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

# Make i2_helium_md importable when this script is run from anywhere.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np

from i2_helium_md.physics.constants import EV, MASS_I_AMU, U
from i2_helium_md.physics.leapfrog import make_neutral_step
from i2_helium_md.physics.potentials import droplet_potential
from i2_helium_md.presets import single_pulse_N2000


# ---------------------------------------------------------------------------
# Setup parameters (must match run_matlab.m exactly)
# ---------------------------------------------------------------------------
NUM_STEPS_TO_RUN: int = 100
DT_PS: float = 0.01

# Initial geometry (atoms placed by hand; no RNG involvement).
# Both atoms at rest. Bond along +/- x. Center of mass at origin.
#
# Bond length: 3.0 Å -- chosen because:
#   * Morse equilibrium is at 2.666 Å, so 3.0 Å is on the repulsive
#     side -> atoms will oscillate around 2.666 (good for testing
#     leapfrog dynamics).
#   * Far from the Xdip perturbation centered at r=9 Å, so the
#     dynamics are simple Morse (no Gaussian dip artifact).
#   * Atoms remain deep inside the droplet, so the droplet
#     potential contributes essentially nothing -> Morse-only test.
BOND_LENGTH_A: float = 3.0      # Angstrom
DROPLET_RADIUS_A: float = 27.97 # Angstrom (= 2.22*2000^(1/3) MATLAB convention,
                                 #          ≈ 2.2173*2000^(1/3) modern)

# Output file
OUT_CSV = Path(__file__).parent / "python_trajectory.csv"


def main() -> None:
    # ------------------------------------------------------------------
    # Build cfg with collisions disabled and partner interaction enabled.
    # We reuse single_pulse_N2000 for sensible defaults (Morse params,
    # potential steepness, binding energy), then override the bits we
    # need.
    # ------------------------------------------------------------------
    cfg = single_pulse_N2000(num_molecules=1)
    cfg = replace(
        cfg,
        single_pulse=False,                          # we control v0 directly
        partner_interaction=True,
        Xdip_active=False,                           # disable Gaussian dip in Morse
        geometric_scattering_crosssection_I=0.0,    # disables collisions
        dt_neutral=DT_PS,
        seed=0,                                      # not used since no rng draws
    )

    # ------------------------------------------------------------------
    # Hand-built initial state (no rng involved).
    # 2N=2 atom layout: atom 1 at index 0, atom 2 at index 1.
    # ------------------------------------------------------------------
    N = 1   # one molecule
    half = BOND_LENGTH_A / 2.0
    x0 = np.array([+half, -half])  # atom 1 at +half, atom 2 at -half
    y0 = np.zeros(2)
    z0 = np.zeros(2)
    vx0 = np.zeros(2)
    vy0 = np.zeros(2)
    vz0 = np.zeros(2)
    mass_kg = np.full(2, MASS_I_AMU * U)
    droplet_radii = np.full(2, DROPLET_RADIUS_A)

    # ------------------------------------------------------------------
    # Build the leapfrog step closure once.
    # ------------------------------------------------------------------
    step_fn = make_neutral_step(cfg, mass_kg, droplet_radii)

    # ------------------------------------------------------------------
    # Output buffer: one row per step (including t=0).
    # Columns:
    # t, x1, y1, z1, x2, y2, z2, vx1, vy1, vz1, vx2, vy2, vz2,
    # E_kin, E_pot, E_total
    # ------------------------------------------------------------------
    rows = []

    # E_pot at t=0 = droplet + half partner Morse per atom
    def energies_at(x, y, z, vx, vy, vz, E_pot_partner_per_pair):
        """Compute E_kin, E_pot, E_total in eV for a given state."""
        v_sq = vx ** 2 + vy ** 2 + vz ** 2
        E_kin = 0.5 * mass_kg * (v_sq * 100.0 ** 2) / EV  # A/ps -> m/s factor 100
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        E_drop = droplet_potential(
            r - droplet_radii,
            steepness=cfg.potential_steepness,
            binding_energy=cfg.binding_energy_I_atom_eV,
        )
        E_partner_per_atom = np.tile(E_pot_partner_per_pair, 2) / 2.0
        E_pot = E_drop + E_partner_per_atom
        return E_kin.sum(), E_pot.sum(), (E_kin + E_pot).sum()

    # Compute initial Morse via one zero-dt step (or the partner_interaction
    # function directly). Easier: use the step function with dt=0 to get
    # the partner energy at the initial configuration. But step_fn doesn't
    # accept dt=0 cleanly; instead, call partner_interaction_neutral().
    from i2_helium_md.physics.interactions import partner_interaction_neutral
    _ax, _ay, _az, E_pot_partner_t0 = partner_interaction_neutral(
        x0, y0, z0, mass_kg, cfg,
    )
    Ek, Ep, Et = energies_at(x0, y0, z0, vx0, vy0, vz0, E_pot_partner_t0)
    rows.append((0.0, x0[0], y0[0], z0[0], x0[1], y0[1], z0[1],
                 vx0[0], vy0[0], vz0[0], vx0[1], vy0[1], vz0[1],
                 Ek, Ep, Et))

    # Step loop
    pos = (x0.copy(), y0.copy(), z0.copy())
    vel = (vx0.copy(), vy0.copy(), vz0.copy())

    for t_idx in range(1, NUM_STEPS_TO_RUN + 1):
        new_pos, new_vel, E_pot_partner = step_fn(pos, vel, DT_PS)
        x, y, z = new_pos
        vx, vy, vz = new_vel
        Ek, Ep, Et = energies_at(x, y, z, vx, vy, vz, E_pot_partner)
        t_ps = t_idx * DT_PS
        rows.append((t_ps, x[0], y[0], z[0], x[1], y[1], z[1],
                     vx[0], vy[0], vz[0], vx[1], vy[1], vz[1],
                     Ek, Ep, Et))
        pos, vel = new_pos, new_vel

    # ------------------------------------------------------------------
    # Write CSV
    # ------------------------------------------------------------------
    header = ("t_ps,"
              "x1_A,y1_A,z1_A,x2_A,y2_A,z2_A,"
              "vx1_Aps,vy1_Aps,vz1_Aps,vx2_Aps,vy2_Aps,vz2_Aps,"
              "E_kin_eV,E_pot_eV,E_total_eV")
    with open(OUT_CSV, "w") as f:
        f.write(header + "\n")
        for row in rows:
            f.write(",".join(f"{v:.16e}" for v in row) + "\n")

    print(f"Wrote {OUT_CSV} with {len(rows)} rows.")
    print(f"\nFinal state:")
    last = rows[-1]
    print(f"  t = {last[0]:.4f} ps")
    print(f"  atom 1 x = {last[1]:.6f} A   (started at {BOND_LENGTH_A/2:.4f})")
    print(f"  atom 2 x = {last[4]:.6f} A   (started at {-BOND_LENGTH_A/2:.4f})")
    print(f"  bond length = {last[1] - last[4]:.6f} A")
    print(f"  E_kin = {last[13]:.6f} eV")
    print(f"  E_pot = {last[14]:.6f} eV")
    print(f"  E_total = {last[15]:.6f} eV")
    print(f"\nE_total drift: "
          f"{abs(rows[-1][15] - rows[0][15]) / abs(rows[0][15]) * 100:.4f}%")


if __name__ == "__main__":
    main()
