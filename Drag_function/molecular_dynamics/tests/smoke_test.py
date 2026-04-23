"""Quick smoke test without pytest -- verifies Steps 1-3 wire up."""

import math
import sys
from pathlib import Path

# make package importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from iodine_he_sim import SimConfig, single_pulse_N2000
from iodine_he_sim.physics import EV, K_B, MASS_I_AMU, U


def close(a, b, rtol=1e-12):
    return math.isclose(a, b, rel_tol=rtol, abs_tol=1e-30)


failures = []
def check(name, cond):
    if cond:
        print(f"  OK   {name}")
    else:
        print(f"  FAIL {name}")
        failures.append(name)


print("=== Constants ===")
check("U == 1.66053907e-27", close(U, 1.66053907e-27))
check("EV == 1.602e-19", close(EV, 1.602e-19))
check("K_B == 1.380649e-23", close(K_B, 1.380649e-23))
check("MASS_I_AMU == 127", MASS_I_AMU == 127.0)

print("\n=== SimConfig derived quantities ===")
cfg = SimConfig()
expected_bind = 318.43 * K_B / EV
check(
    f"binding_energy_I_atom_eV == 318.43*k_B/eV ({expected_bind:.6e})",
    close(cfg.binding_energy_I_atom_eV, expected_bind),
)

expected_mol = 573.3 * K_B / EV * 1000.0
check(
    f"binding_energy_molecule_meV == 573.3*k_B/eV*1000 ({expected_mol:.6f})",
    close(cfg.binding_energy_molecule_meV, expected_mol),
)

check(
    "num_timesteps_neutral == 20000 for t_max=200, dt=0.01",
    cfg.num_timesteps_neutral == 20000,
)

check(
    "v_limit: 40 m/s -> 0.4 A/ps",
    close(cfg.v_limit_angstrom_per_ps, 0.4),
)

expected_emin = MASS_I_AMU * U * (40.0 ** 2) / 2.0 / EV
check(
    f"E_min_eV matches MATLAB formula ({expected_emin:.6e} eV)",
    close(cfg.E_min_eV, expected_emin),
)

print("\n=== single_pulse_N2000 preset ===")
p = single_pulse_N2000()
check("R0_GS_angstrom == 9.0", p.R0_GS_angstrom == 9.0)
check("num_molecules == 2000", p.num_molecules == 2000)
check("single_droplet_size == 2000", p.single_droplet_size == 2000)
check("geometric_scattering_crosssection_I == 30", p.geometric_scattering_crosssection_I == 30.0)
check("geometric_scattering_crosssection_Iplus == 2500", p.geometric_scattering_crosssection_Iplus == 2500.0)
check("binding_energy_I_ion_eV == 0.3", p.binding_energy_I_ion_eV == 0.3)
check("hard_sphere_collision_mode == 3", p.hard_sphere_collision_mode == 3)
check("T_particles_K == 0.4", p.T_particles_K == 0.4)
check("Xdip_active is True", p.Xdip_active is True)
check("lambda_pump_nm == 630", p.lambda_pump_nm == 630.0)
check("E_diss_eV == 1.556", p.E_diss_eV == 1.556)
check("v_limit_m_per_s == 40", p.v_limit_m_per_s == 40.0)
check("sigma_ion_exponent == -2", p.sigma_ion_exponent == -2.0)

print("\n=== overrides ===")
p2 = single_pulse_N2000(num_molecules=500, seed=42)
check("override num_molecules=500", p2.num_molecules == 500)
check("override seed=42", p2.seed == 42)
check("non-overridden stays at preset (R0=9)", p2.R0_GS_angstrom == 9.0)

print("\n=== validation ===")
try:
    single_pulse_N2000(num_molecules=0).validate()
    check("validate() rejects num_molecules=0", False)
except ValueError:
    check("validate() rejects num_molecules=0", True)

try:
    cfg_bad = single_pulse_N2000()
    cfg_bad.hard_sphere_collision_mode = 5
    cfg_bad.validate()
    check("validate() rejects collision mode 5", False)
except ValueError:
    check("validate() rejects collision mode 5", True)


print()
if failures:
    print(f"{len(failures)} FAILURES")
    for f in failures:
        print("  -", f)
    sys.exit(1)
else:
    print("ALL CHECKS PASSED")
