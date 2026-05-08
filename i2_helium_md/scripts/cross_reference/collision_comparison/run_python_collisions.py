"""Drive a long collision-only simulation in Python and write trajectories
+ collision events for statistical comparison with MATLAB.

Setup:
    Atoms loaded from `init_state.csv` (deterministic, no RNG).
    No Morse, no droplet force. Atoms free-stream between collisions.
    Mode-3 collision sampler is invoked each step. When an atom collides,
    its velocity is updated by the elastic-scattering kinematics in
    `apply_collision`. Collisions only happen if the atom is inside the
    droplet (always true here -- atoms are deep inside and won't reach
    the surface in 1 ps).

Run with:
    python run_python_collisions.py

Outputs:
    python_summary.csv  -- one row per timestep:
        t_ps, mean_E_kin_eV, var_E_kin_eV,
        mean_E_dissip_eV, var_E_dissip_eV,
        n_collisions_this_step, n_collisions_cumulative
    python_collision_events.csv -- one row per collision event:
        t_ps, atom_idx, dE_eV, E0_eV
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

# Make i2_helium_md importable
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np

from i2_helium_md.physics.collisions import apply_collision, sample_collision_events
from i2_helium_md.physics.constants import EV, MASS_I_AMU, U
from i2_helium_md.presets import single_pulse_N2000


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
NUM_STEPS = 200
DT_PS = 0.01
DROPLET_RADIUS_A = 27.97
RNG_SEED = 12345
SCRIPT_DIR = Path(__file__).parent
INIT_CSV = SCRIPT_DIR / "init_state.csv"
SUMMARY_CSV = SCRIPT_DIR / "python_summary.csv"
EVENTS_CSV = SCRIPT_DIR / "python_collision_events.csv"


def main() -> None:
    if not INIT_CSV.exists():
        raise FileNotFoundError(
            f"{INIT_CSV} not found. Run generate_init_state.py first."
        )

    # Load init state
    data = np.loadtxt(INIT_CSV, delimiter=",", skiprows=1)
    n_atoms = data.shape[0]
    x = data[:, 1].copy()
    y = data[:, 2].copy()
    z = data[:, 3].copy()
    vx = data[:, 4].copy()
    vy = data[:, 5].copy()
    vz = data[:, 6].copy()

    print(f"Loaded {n_atoms} atoms from {INIT_CSV}")

    # Build cfg with collisions enabled, partner & droplet disabled.
    cfg = single_pulse_N2000(num_molecules=n_atoms // 2)  # only used for parameters
    cfg = replace(
        cfg,
        single_pulse=False,
        partner_interaction=False,    # no Morse
        Xdip_active=False,
        # Use the production-default cross section etc.:
        # geometric_scattering_crosssection_I (default 30 A^2)
        # scatter_mass_neutral_amu (4.0 = He)
        # neutral_scatter_angle_std_deg (default 0)
        # E_min_eV from cfg property
        dt_neutral=DT_PS,
    )

    rng = np.random.default_rng(RNG_SEED)
    mass_kg = np.full(n_atoms, MASS_I_AMU * U)
    droplet_radii = np.full(n_atoms, DROPLET_RADIUS_A)
    masses_amu = mass_kg / U

    # Per-atom cumulative tracking
    E_dissip_cum = np.zeros(n_atoms)
    n_coll_cum = np.zeros(n_atoms, dtype=int)
    prev_distance = None  # first step has no previous distance

    # Output buffers
    summary_rows = []
    event_rows = []   # (t_ps, atom_idx, dE_eV, E0_eV)

    def compute_E_kin(vx, vy, vz):
        v_sq = vx ** 2 + vy ** 2 + vz ** 2
        return 0.5 * mass_kg * (v_sq * 100.0 ** 2) / EV  # per atom

    # Record initial state (t=0)
    E_kin_per_atom = compute_E_kin(vx, vy, vz)
    summary_rows.append((
        0.0,
        E_kin_per_atom.mean(),
        E_kin_per_atom.var(ddof=0),
        E_dissip_cum.mean(),
        E_dissip_cum.var(ddof=0),
        0,    # collisions this step
        0,    # cumulative
    ))

    # Step loop
    for step_idx in range(1, NUM_STEPS + 1):
        # 1. Free-stream the atoms (no forces, since partner & droplet off)
        x_new = x + vx * DT_PS
        y_new = y + vy * DT_PS
        z_new = z + vz * DT_PS

        # 2. Compute depth (atoms are deep inside; depth < 0)
        r_new = np.sqrt(x_new ** 2 + y_new ** 2 + z_new ** 2)
        depth = r_new - droplet_radii

        # 3. Pre-collision kinetic energy (after free-streaming)
        E0 = compute_E_kin(vx, vy, vz)   # velocity hasn't changed yet

        # 4. Sample collision events (Mode 3)
        if prev_distance is None:
            b_collision = np.zeros(n_atoms, dtype=bool)
        else:
            b_collision = sample_collision_events(
                distance_travelled_angstrom=prev_distance,
                depth_angstrom=depth,
                E0_eV=E0,
                sigma_angstrom_sq=cfg.geometric_scattering_crosssection_I,
                E_min_eV=cfg.E_min_eV,
                rng=rng,
            )

        # 5. Apply collisions
        vx_after, vy_after, vz_after, dE_eV = apply_collision(
            vx=vx, vy=vy, vz=vz,
            masses_amu=masses_amu,
            b_collision=b_collision,
            scatter_mass_amu=cfg.scatter_mass_neutral_amu,
            neutral_scatter_angle_std_deg=cfg.neutral_scatter_angle_std_deg,
            rng=rng,
        )

        # 6. Bookkeeping
        E_dissip_cum += dE_eV
        n_coll_cum += b_collision.astype(int)

        # Log per-collision events
        n_coll_this_step = int(b_collision.sum())
        if n_coll_this_step > 0:
            t_ps = step_idx * DT_PS
            for atom_idx in np.flatnonzero(b_collision):
                event_rows.append((
                    t_ps,
                    int(atom_idx),
                    float(dE_eV[atom_idx]),
                    float(E0[atom_idx]),
                ))

        # 7. Compute distance traveled this step (used as prev_distance next iter)
        prev_distance = np.sqrt(
            (x_new - x) ** 2 + (y_new - y) ** 2 + (z_new - z) ** 2
        )

        # 8. Advance to new state
        x, y, z = x_new, y_new, z_new
        vx, vy, vz = vx_after, vy_after, vz_after

        # Record summary
        E_kin_per_atom = compute_E_kin(vx, vy, vz)
        t_ps = step_idx * DT_PS
        summary_rows.append((
            t_ps,
            E_kin_per_atom.mean(),
            E_kin_per_atom.var(ddof=0),
            E_dissip_cum.mean(),
            E_dissip_cum.var(ddof=0),
            n_coll_this_step,
            int(n_coll_cum.sum()),
        ))

    # Write CSVs
    with open(SUMMARY_CSV, "w") as f:
        f.write("t_ps,mean_E_kin_eV,var_E_kin_eV,"
                "mean_E_dissip_eV,var_E_dissip_eV,"
                "n_coll_this_step,n_coll_cumulative\n")
        for row in summary_rows:
            f.write(",".join(
                (f"{v:.16e}" if not isinstance(v, int) else str(v))
                for v in row
            ) + "\n")

    with open(EVENTS_CSV, "w") as f:
        f.write("t_ps,atom_idx,dE_eV,E0_eV\n")
        for row in event_rows:
            f.write(f"{row[0]:.16e},{row[1]},{row[2]:.16e},{row[3]:.16e}\n")

    # Summary
    print(f"\nWrote {SUMMARY_CSV} with {len(summary_rows)} rows")
    print(f"Wrote {EVENTS_CSV} with {len(event_rows)} collision events")
    print(f"\nFinal state:")
    print(f"  Total collisions: {len(event_rows)}")
    print(f"  Avg collisions per atom: {len(event_rows) / n_atoms:.2f}")
    print(f"  Mean E_kin at end: {summary_rows[-1][1]:.6f} eV "
          f"(started at {summary_rows[0][1]:.6f})")
    print(f"  Mean E_dissip at end: {summary_rows[-1][3]:.6f} eV")
    print(f"  Energy balance check (mean E_kin start - end vs E_dissip end): "
          f"{summary_rows[0][1] - summary_rows[-1][1]:.6f} vs "
          f"{summary_rows[-1][3]:.6f} eV  (should agree)")


if __name__ == "__main__":
    main()
