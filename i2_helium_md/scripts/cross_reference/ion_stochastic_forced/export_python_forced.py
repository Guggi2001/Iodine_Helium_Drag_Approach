"""Export Python ion forced-event trajectory for cross-reference.

CLAUDE.md validation target 5 (driver-level stochastic plumbing).
The collision kernel itself is covered by
``scripts/collision_comparison_test/``; this script verifies the
driver bookkeeping around it.

Uses the implemented modular functions: ``build_initial_ion_state``
+ ``ion_propagation_step``. Stochastic event probabilities are forced
to >= 1 entirely via cfg flags (no code paths bypassed):

    geometric_scattering_crosssection_Iplus = 1e6 A^2  (forces collision)
    sigma_dependent_on_v = True, sigma_ion_exponent = -2
    mass_attach_probability = 1.0                     (forces attachment)
    ion_scatter_angle_std_deg = 0.0                   (no extra smearing)

Run:
    python scripts/cross_reference/ion_stochastic_forced/export_python_forced.py

Output: python_forced.csv with one row per (step, atom). Columns:
    step, t_ps, atom,
    x_A, y_A, z_A, vx_Aps, vy_Aps, vz_Aps,
    mass_kg, E_kin_eV, E_pot_eV, E_dissip_eV, E_mass_attach_defect_eV,
    number_of_collisions, b_collision, b_attach, sigma_used_A2, depth_A
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from i2_helium_md.physics.collisions import velocity_dependent_cross_section
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
OUT_CSV = SCRIPT_DIR / "python_forced.csv"


def _build_neutral_checkpoint(inputs: dict) -> NeutralCheckpoint:
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
    two_N = 2 * N
    num_steps = int(inputs["num_steps"])
    dt_ion = float(inputs["dt_ion_ps"])

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
        sigma_ion_exponent=float(flags["sigma_ion_exponent"]),
        mass_attach_probability=float(flags["mass_attach_probability"]),
        scatter_mass_ion_amu=float(flags["scatter_mass_ion_amu"]),
        ion_scatter_angle_std_deg=float(flags["ion_scatter_angle_std_deg"]),
        v_limit_m_per_s=float(flags["v_limit_m_per_s"]),
        additional_droplet_charges=int(flags["additional_droplet_charges"]),
        single_charge_ionization_allowed=bool(flags["single_charge_ionization_allowed"]),
        highly_charged_iodine=bool(flags["highly_charged_iodine"]),
        effusive_dynamics=bool(flags["effusive_dynamics"]),
        hard_sphere_collision_mode=int(flags["hard_sphere_collision_mode"]),
    )

    neutral_ckpt = _build_neutral_checkpoint(inputs)
    ckpt = build_initial_ion_state(cfg, neutral_ckpt, num_steps_ion=1)
    state = ion_state_from_checkpoint_column(ckpt, 0)

    # Deterministic seeded RNG. With p_collision >= 1 and p_attach = 1
    # the actual draws don't influence event flags, but the RNG stream
    # still drives the impact-parameter, reference-direction, and
    # azimuth draws inside apply_collision -> the post-collision
    # velocities depend on the seed (and won't match MATLAB).
    rng = np.random.default_rng(0)
    droplet_radii = ckpt.droplet_radii_angstrom
    charge = np.ones(two_N, dtype=float)

    # sigma at t=0 (recorded for cross-language analytic check on the
    # v-dependent cross-section code path).
    v0_speed = np.sqrt(state.vx ** 2 + state.vy ** 2 + state.vz ** 2)
    sigma_t0 = velocity_dependent_cross_section(
        v0_speed,
        sigma_0_angstrom_sq=cfg.geometric_scattering_crosssection_Iplus,
        exponent=cfg.sigma_ion_exponent,
    )

    rows: list[tuple] = []

    def _record(step_idx: int, st, sigma_used: np.ndarray, b_collision: np.ndarray, b_attach: np.ndarray) -> None:
        depth = np.sqrt(st.x ** 2 + st.y ** 2 + st.z ** 2) - droplet_radii
        for atom in range(two_N):
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
                float(st.E_mass_attach_defect_eV[atom]),
                int(st.number_of_collisions[atom]),
                int(b_collision[atom]),
                int(b_attach[atom]),
                float(sigma_used[atom]),
                float(depth[atom]),
            ))

    # Step 0: no collision, no attachment, sigma at t=0 from initial v.
    b_zero = np.zeros(two_N, dtype=int)
    _record(0, state, sigma_t0, b_zero, b_zero)

    prev_dist: np.ndarray | None = None
    prev_n_coll = state.number_of_collisions.copy()
    prev_mass = state.mass_kg.copy()

    for step in range(1, num_steps + 1):
        # sigma_used at the START of this step (before any update),
        # using the pre-step velocity. This is what the driver feeds
        # into sample_collision_events.
        v_pre = np.sqrt(state.vx ** 2 + state.vy ** 2 + state.vz ** 2)
        if cfg.sigma_dependent_on_v:
            sigma_used = velocity_dependent_cross_section(
                v_pre,
                sigma_0_angstrom_sq=cfg.geometric_scattering_crosssection_Iplus,
                exponent=cfg.sigma_ion_exponent,
            )
        else:
            sigma_used = np.full(two_N, cfg.geometric_scattering_crosssection_Iplus)

        new_state = ion_propagation_step(
            state,
            cfg=cfg,
            droplet_radii=droplet_radii,
            charge=charge,
            prev_distance_angstrom=prev_dist,
            rng=rng,
        )

        # Recover per-step b_collision, b_attach by differencing.
        n_coll_delta = new_state.number_of_collisions - prev_n_coll
        b_collision_step = (n_coll_delta > 0).astype(int)
        # mass increases by exactly 4u per attachment per step here
        mass_delta = new_state.mass_kg - prev_mass
        b_attach_step = (mass_delta > 0.5 * 4.0 * U).astype(int)

        prev_dist = np.sqrt(
            (new_state.x - state.x) ** 2
            + (new_state.y - state.y) ** 2
            + (new_state.z - state.z) ** 2
        )
        prev_n_coll = new_state.number_of_collisions.copy()
        prev_mass = new_state.mass_kg.copy()
        state = new_state

        _record(step, state, sigma_used, b_collision_step, b_attach_step)

    header = (
        "step,t_ps,atom,"
        "x_A,y_A,z_A,vx_Aps,vy_Aps,vz_Aps,"
        "mass_kg,E_kin_eV,E_pot_eV,E_dissip_eV,E_mass_attach_defect_eV,"
        "number_of_collisions,b_collision,b_attach,sigma_used_A2,depth_A"
    )
    with open(OUT_CSV, "w") as f:
        f.write(header + "\n")
        for r in rows:
            (step_idx, t_ps, atom,
             x, y, z, vx, vy, vz,
             m, ek, ep, ed, edm,
             nc, bc, ba, sg, dp) = r
            f.write(
                f"{step_idx},{t_ps:.16e},{atom},"
                f"{x:.16e},{y:.16e},{z:.16e},"
                f"{vx:.16e},{vy:.16e},{vz:.16e},"
                f"{m:.16e},{ek:.16e},{ep:.16e},{ed:.16e},{edm:.16e},"
                f"{nc},{bc},{ba},{sg:.16e},{dp:.16e}\n"
            )

    print(f"Wrote {OUT_CSV} with {len(rows)} rows ({num_steps + 1} steps x {two_N} atoms).")
    # Per-side conservation invariant (must hold modulo Verlet drift):
    rows_arr = np.array([r[10:14] for r in rows])  # E_kin, E_pot, E_dissip, E_defect
    inv = (rows_arr[:, 0] + rows_arr[:, 1] + rows_arr[:, 2] + rows_arr[:, 3]).reshape(num_steps + 1, two_N).sum(axis=1)
    print(f"  E_kin+E_pot+E_dissip+E_defect at step 0   : {inv[0]:.6f} eV")
    print(f"  same at final step                        : {inv[-1]:.6f} eV")
    print(f"  drift                                     : {inv[-1] - inv[0]:+.3e} eV "
          f"({(inv[-1] - inv[0]) / inv[0] * 100:+.4f}%)")


if __name__ == "__main__":
    main()
