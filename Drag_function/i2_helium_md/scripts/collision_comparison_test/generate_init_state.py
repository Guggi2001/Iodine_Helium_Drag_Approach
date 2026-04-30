"""Generate a deterministic initial state for the Python/MATLAB collision
comparison test. Writes init_state.csv with one row per atom.

Setup:
    NUM_ATOMS atoms placed deterministically inside a single He droplet.
    Each atom gets a fixed initial velocity (different per atom but
    reproducible from a hash of its index, no RNG involved). The atoms
    are non-interacting (no Morse, no droplet force) -- they just
    free-stream until they collide with He atoms in the droplet.

    This isolates the collision sampler and the elastic-scattering
    kinematics, the most stochastic part of the simulation.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
NUM_ATOMS = 1000             # for good statistics
DROPLET_RADIUS_A = 27.97     # Angstrom
INITIAL_SPEED_AA_PER_PS = 5.0  # atoms start with this speed (~ 0.4 eV)

OUT_CSV = Path(__file__).parent / "init_state.csv"


def main() -> None:
    # Deterministic placement: atoms on a coarse 3D grid inside the
    # droplet, with deterministic velocity directions from a hashed
    # (Halton-like) sequence. No RNG involved.

    # Place atoms on a cubic grid spanning a box smaller than the droplet.
    # We oversize the grid to ensure enough points fit after the spherical
    # filter, then take the first NUM_ATOMS.
    side = int(np.ceil((NUM_ATOMS * 2.5) ** (1.0 / 3.0)))
    box_half = DROPLET_RADIUS_A * 0.5  # keep deep inside the droplet
    coords_1d = np.linspace(-box_half, box_half, side)
    xx, yy, zz = np.meshgrid(coords_1d, coords_1d, coords_1d, indexing='ij')
    pos_all = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1)

    # Filter to those inside r = 0.7 R (well clear of surface)
    r = np.linalg.norm(pos_all, axis=1)
    inside = pos_all[r < 0.7 * DROPLET_RADIUS_A]

    # Take exactly NUM_ATOMS (truncate or warn if not enough)
    if inside.shape[0] < NUM_ATOMS:
        raise RuntimeError(
            f"only {inside.shape[0]} grid points fit inside; "
            "increase grid resolution"
        )
    pos = inside[:NUM_ATOMS]
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

    # Velocity directions from a deterministic Halton sequence so they're
    # reproducible across implementations
    def halton(i: int, base: int) -> float:
        """Halton sequence value for index i in given base."""
        f = 1.0
        result = 0.0
        while i > 0:
            f /= base
            result += f * (i % base)
            i //= base
        return result

    vx = np.zeros(NUM_ATOMS)
    vy = np.zeros(NUM_ATOMS)
    vz = np.zeros(NUM_ATOMS)
    for i in range(NUM_ATOMS):
        # Map two Halton values to spherical (theta, phi)
        u1 = halton(i + 1, 2)
        u2 = halton(i + 1, 3)
        theta = np.arccos(2 * u1 - 1)        # uniform on sphere
        phi = 2 * np.pi * u2
        vx[i] = INITIAL_SPEED_AA_PER_PS * np.sin(theta) * np.cos(phi)
        vy[i] = INITIAL_SPEED_AA_PER_PS * np.sin(theta) * np.sin(phi)
        vz[i] = INITIAL_SPEED_AA_PER_PS * np.cos(theta)

    # Write CSV
    with open(OUT_CSV, "w") as f:
        f.write("idx,x_A,y_A,z_A,vx_Aps,vy_Aps,vz_Aps\n")
        for i in range(NUM_ATOMS):
            f.write(f"{i},{x[i]:.16e},{y[i]:.16e},{z[i]:.16e},"
                    f"{vx[i]:.16e},{vy[i]:.16e},{vz[i]:.16e}\n")

    speed_check = np.sqrt(vx**2 + vy**2 + vz**2)
    print(f"Wrote {OUT_CSV} with {NUM_ATOMS} atoms.")
    print(f"  Position range: r in [{r[r < 0.7 * DROPLET_RADIUS_A].min():.2f}, "
          f"{r[r < 0.7 * DROPLET_RADIUS_A].max():.2f}] A "
          f"(droplet radius {DROPLET_RADIUS_A} A)")
    print(f"  Speed (all atoms): {speed_check.min():.4f}-{speed_check.max():.4f} A/ps "
          f"(target {INITIAL_SPEED_AA_PER_PS})")
    print(f"  Mean velocity vector: ({vx.mean():.4f}, {vy.mean():.4f}, {vz.mean():.4f}) "
          "A/ps  (should be ~0)")


if __name__ == "__main__":
    main()
