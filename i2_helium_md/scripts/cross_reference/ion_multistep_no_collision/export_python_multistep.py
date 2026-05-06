"""Export the Python ion multi-step trajectory for cross-reference.

CLAUDE.md validation targets 3 and 4: several deterministic ion steps
with collisions disabled, plus energy bookkeeping.

This script uses the implemented modular functions
(``build_initial_ion_state`` and ``ion_propagation_step``) -- no
duplicated physics. All stochastic ion physics is disabled at the
configuration level:

* ``geometric_scattering_crosssection_Iplus = 0.0`` -> sigma = 0,
  so collision probability is identically zero,
* ``sigma_dependent_on_v = False``         -> use the constant 0,
* ``mass_attach_probability = 0.0``        -> no attachment,
* ``ion_scatter_angle_std_deg = 0.0``      -> elastic scattering kept
  unused since no collisions occur,
* ``additional_droplet_charges = 0``       -> no extra-charge model.

The RNG is consumed but its draws don't affect dynamics with the above
flags. We seed deterministically anyway for reproducibility.

Run:
    python scripts/cross_reference/ion_multistep_no_collision/export_python_multistep.py

Output: python_multistep.csv with one row per (step, atom).
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

# Make i2_helium_md importable when run directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from i2_helium_md.physics.constants import U
from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.checkpoint import (
    NeutralCheckpoint,
    _NEUTRAL_SCHEMA_VERSION,
)
from i2_helium_md.simulation.ion_initial_state import build_initial_ion_state
from i2_helium_md.simulation.ion_propagation_step import (
    ion_propagation_step,
    ion_state_from_checkpoint_column,
)


SCRIPT_DIR = Path(__file__).parent
INPUTS_JSON = SCRIPT_DIR / "inputs.json"
OUT_CSV = SCRIPT_DIR / "python_multistep.csv"


def _build_neutral_checkpoint(inputs: dict) -> NeutralCheckpoint:
    """Single-column NeutralCheckpoint built from hand-crafted state."""
    N = int(inputs["num_molecules"])
    two_N = 2 * N
    a0 = inputs["atom_0"]
    a1 = inputs["atom_1"]

    return NeutralCheckpoint(
        num_molecules=N,
        time_ps=np.array([float(inputs["time_ps_at_t0"])]),
        positions_x=np.array([[a0["x_A"]], [a1["x_A"]]]),
        positions_y=np.array([[a0["y_A"]], [a1["y_A"]]]),
        positions_z=np.array([[a0["z_A"]], [a1["z_A"]]]),
        velocities_x=np.array([[a0["vx_Aps"]], [a1["vx_Aps"]]]),
        velocities_y=np.array([[a0["vy_Aps"]], [a1["vy_Aps"]]]),
        velocities_z=np.array([[a0["vz_Aps"]], [a1["vz_Aps"]]]),
        mass_kg=np.full(two_N, float(inputs["mass_amu"]) * U),
        droplet_radii=np.full(two_N, float(inputs["droplet_radius_A"])),
        r0=np.zeros(N),
        E_kin_eV=np.zeros((two_N, 1)),
        E_pot_eV=np.zeros((two_N, 1)),
        E_initial_eV=np.zeros(N),
        E_dissip_eV=np.zeros((two_N, 1)),
        L_droplet_eV_ps=np.zeros((two_N, 1)),
        schema_version=_NEUTRAL_SCHEMA_VERSION,
    )


def main() -> None:
    inputs = json.loads(INPUTS_JSON.read_text())
    N = int(inputs["num_molecules"])
    num_steps = int(inputs["num_steps"])
    dt_ion = float(inputs["dt_ion_ps"])

    # Build cfg with all stochastic ion physics disabled.
    cfg = single_pulse_N2000(num_molecules=N)
    flags = inputs["cfg_flags"]
    cfg = replace(
        cfg,
        dt_ion=dt_ion,
        binding_energy_I_ion_eV=float(flags["binding_energy_I_ion_eV"]),
        potential_steepness=float(flags["potential_steepness"]),
        E_coulomb_scale=float(flags["E_coulomb_scale"]),
        geometric_scattering_crosssection_Iplus=float(
            flags["geometric_scattering_crosssection_Iplus"]
        ),
        sigma_dependent_on_v=bool(flags["sigma_dependent_on_v"]),
        mass_attach_probability=float(flags["mass_attach_probability"]),
        ion_scatter_angle_std_deg=float(flags["ion_scatter_angle_std_deg"]),
        additional_droplet_charges=int(flags["additional_droplet_charges"]),
        single_charge_ionization_allowed=bool(flags["single_charge_ionization_allowed"]),
        highly_charged_iodine=bool(flags["highly_charged_iodine"]),
        effusive_dynamics=bool(flags["effusive_dynamics"]),
    )

    neutral_ckpt = _build_neutral_checkpoint(inputs)

    # Step 0: build the initial ion state (corrected E_kin/E_pot per the
    # documented Python fixes vs legacy MATLAB t=0 bugs).
    ckpt = build_initial_ion_state(cfg, neutral_ckpt, num_steps_ion=1)
    state = ion_state_from_checkpoint_column(ckpt, 0)

    # Stochasticity is OFF via cfg, but the step API still requires an
    # rng. Seed deterministically.
    rng = np.random.default_rng(0)

    droplet_radii = ckpt.droplet_radii_angstrom
    charge = np.ones(2 * N, dtype=float)

    rows: list[tuple] = []

    def _record(step_idx: int, st) -> None:
        for atom in range(2 * N):
            rows.append((
                step_idx,
                float(st.time_ps),
                atom,
                float(st.x[atom]), float(st.y[atom]), float(st.z[atom]),
                float(st.vx[atom]), float(st.vy[atom]), float(st.vz[atom]),
                float(st.mass_kg[atom]),
                float(st.E_kin_eV[atom]),
                float(st.E_pot_eV[atom]),
                float(st.E_dissip_eV[atom]),
            ))

    _record(0, state)

    prev_dist: np.ndarray | None = None
    for step in range(1, num_steps + 1):
        new_state = ion_propagation_step(
            state,
            cfg=cfg,
            droplet_radii=droplet_radii,
            charge=charge,
            prev_distance_angstrom=prev_dist,
            rng=rng,
        )
        prev_dist = np.sqrt(
            (new_state.x - state.x) ** 2
            + (new_state.y - state.y) ** 2
            + (new_state.z - state.z) ** 2
        )
        state = new_state
        _record(step, state)

    header = ("step,t_ps,atom,"
              "x_A,y_A,z_A,"
              "vx_Aps,vy_Aps,vz_Aps,"
              "mass_kg,E_kin_eV,E_pot_eV,E_dissip_eV")
    with open(OUT_CSV, "w") as f:
        f.write(header + "\n")
        for r in rows:
            step_idx, t_ps, atom, x, y, z, vx, vy, vz, m, ek, ep, ed = r
            f.write(
                f"{step_idx},{t_ps:.16e},{atom},"
                f"{x:.16e},{y:.16e},{z:.16e},"
                f"{vx:.16e},{vy:.16e},{vz:.16e},"
                f"{m:.16e},{ek:.16e},{ep:.16e},{ed:.16e}\n"
            )

    print(f"Wrote {OUT_CSV} with {len(rows)} rows ({num_steps + 1} steps x {2*N} atoms).")
    # Quick sanity: total energy at step 0 vs final step.
    rows_arr = np.array([r[10:13] for r in rows])  # E_kin, E_pot, E_dissip
    n_atoms = 2 * N
    E_total_per_step = (rows_arr[:, 0] + rows_arr[:, 1] + rows_arr[:, 2]).reshape(-1, n_atoms).sum(axis=1)
    print(f"  E_total[0]   = {E_total_per_step[0]:.6f} eV")
    print(f"  E_total[end] = {E_total_per_step[-1]:.6f} eV")
    print(f"  drift        = {(E_total_per_step[-1] - E_total_per_step[0]):+.3e} eV")


if __name__ == "__main__":
    main()
