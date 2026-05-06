"""Export the Python ion t=0 state for the first ion-driver cross-reference.

Reads ``inputs.json`` (the shared single source of truth for both
sides), builds a minimal :class:`NeutralCheckpoint` of length T=1 with
the hand-crafted end-state, calls :func:`build_initial_ion_state`, and
writes ``python_t0.csv``.

No propagation step is taken. All stochastic ion physics is disabled
via the cfg flags listed in ``inputs.json``.

Run from anywhere:

    python scripts/cross_reference/ion_t0_state/export_python_t0.py

See ``README.md`` next to this file for the comparison protocol.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

# Make i2_helium_md importable when this script is run directly.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from i2_helium_md.physics.constants import U
from i2_helium_md.presets import single_pulse_N2000
from i2_helium_md.simulation.checkpoint import NeutralCheckpoint, _NEUTRAL_SCHEMA_VERSION
from i2_helium_md.simulation.ion_initial_state import build_initial_ion_state


SCRIPT_DIR = Path(__file__).parent
INPUTS_JSON = SCRIPT_DIR / "inputs.json"
OUT_CSV = SCRIPT_DIR / "python_t0.csv"


def _build_neutral_checkpoint_from_inputs(inputs: dict) -> NeutralCheckpoint:
    """Pack the hand-crafted t=0 state into a single-column NeutralCheckpoint.

    Only the fields actually read by ``build_initial_ion_state`` need
    physically meaningful values: positions, velocities, masses,
    droplet radii, num_molecules. The other dataclass fields
    (``r0``, ``E_initial_eV``, ``time_ps``, ``E_kin_eV``, ``E_pot_eV``,
    ``E_dissip_eV``, ``L_droplet_eV_ps``) are required by the schema but
    not consulted by the ion init builder; we set them to zeros.
    """
    N = int(inputs["num_molecules"])
    two_N = 2 * N
    T = 1   # single end-state column

    # Standard 2N layout: atom 1 of molecule i at row i; atom 2 at row N+i.
    a0 = inputs["atom_0"]
    a1 = inputs["atom_1"]

    pos_x = np.array([[a0["x_A"]], [a1["x_A"]]])
    pos_y = np.array([[a0["y_A"]], [a1["y_A"]]])
    pos_z = np.array([[a0["z_A"]], [a1["z_A"]]])
    vel_x = np.array([[a0["vx_Aps"]], [a1["vx_Aps"]]])
    vel_y = np.array([[a0["vy_Aps"]], [a1["vy_Aps"]]])
    vel_z = np.array([[a0["vz_Aps"]], [a1["vz_Aps"]]])

    # Mass: store in u in inputs, convert to kg here using Python's CODATA U.
    mass_kg = np.full(two_N, float(inputs["mass_amu"]) * U)
    droplet_radii = np.full(two_N, float(inputs["droplet_radius_A"]))

    return NeutralCheckpoint(
        num_molecules=N,
        time_ps=np.array([float(inputs["time_ps_at_t0"])]),
        positions_x=pos_x,
        positions_y=pos_y,
        positions_z=pos_z,
        velocities_x=vel_x,
        velocities_y=vel_y,
        velocities_z=vel_z,
        mass_kg=mass_kg,
        droplet_radii=droplet_radii,
        r0=np.zeros(N),
        E_kin_eV=np.zeros((two_N, T)),
        E_pot_eV=np.zeros((two_N, T)),
        E_initial_eV=np.zeros(N),
        E_dissip_eV=np.zeros((two_N, T)),
        L_droplet_eV_ps=np.zeros((two_N, T)),
        schema_version=_NEUTRAL_SCHEMA_VERSION,
    )


def main() -> None:
    inputs = json.loads(INPUTS_JSON.read_text())

    cfg = single_pulse_N2000(num_molecules=int(inputs["num_molecules"]))
    cfg_flags = inputs["cfg_flags"]
    cfg = replace(
        cfg,
        binding_energy_I_ion_eV=float(cfg_flags["binding_energy_I_ion_eV"]),
        potential_steepness=float(cfg_flags["potential_steepness"]),
        E_coulomb_scale=float(cfg_flags["E_coulomb_scale"]),
        single_charge_ionization_allowed=bool(cfg_flags["single_charge_ionization_allowed"]),
        additional_droplet_charges=int(cfg_flags["additional_droplet_charges"]),
        highly_charged_iodine=bool(cfg_flags["highly_charged_iodine"]),
        effusive_dynamics=bool(cfg_flags["effusive_dynamics"]),
        mass_attach_probability=float(cfg_flags["mass_attach_probability"]),
        sigma_dependent_on_v=bool(cfg_flags["sigma_dependent_on_v"]),
    )

    neutral_ckpt = _build_neutral_checkpoint_from_inputs(inputs)

    ion_ckpt = build_initial_ion_state(
        cfg,
        neutral_ckpt,
        num_steps_ion=int(inputs["num_steps_ion"]),
    )

    # Pull the t=0 column.
    rows = [
        ("x_A",              ion_ckpt.positions_x[0, 0],         ion_ckpt.positions_x[1, 0]),
        ("y_A",              ion_ckpt.positions_y[0, 0],         ion_ckpt.positions_y[1, 0]),
        ("z_A",              ion_ckpt.positions_z[0, 0],         ion_ckpt.positions_z[1, 0]),
        ("vx_Aps",           ion_ckpt.velocities_x[0, 0],        ion_ckpt.velocities_x[1, 0]),
        ("vy_Aps",           ion_ckpt.velocities_y[0, 0],        ion_ckpt.velocities_y[1, 0]),
        ("vz_Aps",           ion_ckpt.velocities_z[0, 0],        ion_ckpt.velocities_z[1, 0]),
        ("mass_kg",          ion_ckpt.mass_kg[0],                ion_ckpt.mass_kg[1]),
        ("droplet_radius_A", ion_ckpt.droplet_radii_angstrom[0], ion_ckpt.droplet_radii_angstrom[1]),
        ("time_ps",          ion_ckpt.time_ps[0],                ion_ckpt.time_ps[0]),
        ("E_kin_eV",         ion_ckpt.E_kin_eV[0, 0],            ion_ckpt.E_kin_eV[1, 0]),
        ("E_pot_eV",         ion_ckpt.E_pot_eV[0, 0],            ion_ckpt.E_pot_eV[1, 0]),
        ("E_dissip_eV",      ion_ckpt.E_dissip_eV[0, 0],         ion_ckpt.E_dissip_eV[1, 0]),
    ]

    with open(OUT_CSV, "w") as f:
        f.write("quantity,atom_0,atom_1\n")
        for name, v0, v1 in rows:
            f.write(f"{name},{v0:.16e},{v1:.16e}\n")

    print(f"Wrote {OUT_CSV} with {len(rows)} quantities.")
    for name, v0, v1 in rows:
        print(f"  {name:<18s}  atom_0 = {v0:+.6e}   atom_1 = {v1:+.6e}")


if __name__ == "__main__":
    main()
