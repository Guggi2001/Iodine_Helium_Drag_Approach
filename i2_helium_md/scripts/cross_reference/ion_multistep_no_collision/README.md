# Ion driver cross-reference — deterministic multi-step (no collisions)

CLAUDE.md ion-stage validation targets **3** (several deterministic
ion steps with collisions disabled) and **4** (energy bookkeeping in
deterministic mode).

## What is tested

A side-by-side comparison of the Python ion driver vs. the legacy
MATLAB ion driver on a small **fully deterministic** multi-step
trajectory:

- 1 molecule (2 iodine ions) in the standard 2N layout.
- 100 leapfrog (velocity-Verlet) steps at `dt_ion = 0.01 ps` → 1 ps
  total, starting from the same hand-crafted state on both sides.
- Force model = ion-droplet potential + partner Coulomb only.
- All stochastic ion physics disabled at the **configuration** level
  (no implementation hacks):
  - `geometric_scattering_crosssection_Iplus = 0`
  - `sigma_dependent_on_v = false`
  - `mass_attach_probability = 0`
  - `ion_scatter_angle_std_deg = 0`
  - `additional_droplet_charges = 0`

The scenario is small enough to inspect by eye: at t=0 the atoms sit
inside a tiny droplet (R = 4 Å) with non-zero z and asymmetric per-
atom velocities, and Coulomb repulsion accelerates them apart over
the 1 ps run.

## What is *not* tested here

- Stochastic events (collisions, attachment, v-dependent σ, scatter
  angles): covered by the sibling
  `../ion_stochastic_forced/` task (driver plumbing) and by
  `scripts/collision_comparison_test/` (kinematics).
- t=0 bookkeeping bugs (E_kin / E_pot at the ion stage start):
  covered by `../ion_t0_state/`. Both this script and the MATLAB
  reference here use the **corrected** t=0 formulas, so t=0 is a
  normal CROSS_MATCH row and not a bug-tagged row.

## Files

| File | Role |
|---|---|
| `inputs.json` | Hand-crafted neutral end-state shared by both sides. Single source of truth for positions, velocities, mass, droplet radius, and disabled-physics flags. |
| `export_python_multistep.py` | Calls `build_initial_ion_state` + loops `ion_propagation_step`. Writes `python_multistep.csv` (one row per (step, atom)). |
| `export_matlab_multistep.m` | Inline reimplementation of `frog_step_ion.m` + the per-step bookkeeping block of `vmi_sim_3d_ion_propa.m` (lines 761, 765). Uses MATLAB legacy-rounded constants (`eV=1.602e-19`, `u=1.66053907e-27`, force-conversion literals `1.602e-9`, `9648.53322`). Writes `matlab_multistep.csv`. |
| `compare_multistep.py` | Loads both CSVs, prints per-quantity max abs / scaled-relative diffs plus per-side `E_total` drift. Non-zero exit only on disagreement beyond expected tolerance. |
| `python_multistep.csv` / `matlab_multistep.csv` | Generated outputs (committed for reference). |

## How to run

```bash
# Python
python scripts/cross_reference/ion_multistep_no_collision/export_python_multistep.py

# MATLAB
matlab -batch "cd('scripts/cross_reference/ion_multistep_no_collision'); export_matlab_multistep"

# Compare
python scripts/cross_reference/ion_multistep_no_collision/compare_multistep.py
```

## Quantities compared

Per (step, atom): `x, y, z, vx, vy, vz, mass_kg, E_kin_eV, E_pot_eV,
E_dissip_eV`, plus per-step system total energy `E_total_eV`.

## Pass / fail criteria

| Quantity | Category | Tolerance |
|---|---|---|
| `t_ps`, `E_dissip_eV` | **MATCH** | exact (must be 0) |
| `x, y, z, vx, vy, vz, mass_kg, E_kin_eV, E_pot_eV, E_total_eV` | **CONSTANT_ROUNDING** | scaled-relative ≤ 5×10⁻⁴ (500 ppm) |

The ~110 ppm budget reflects the difference between Python's CODATA
2022 constants (`EV` exact, `U = 1.66053906892e-27`) and the MATLAB
legacy 4-significant-figure values (`eV = 1.602e-19`, `u = 1.66053907e-27`)
plus the literal `1.602e-9` / `9648.53322` factors inside the force
conversions. The 500 ppm gate gives a 5× headroom on top of that
budget for compounding over the 100-step trajectory.

## Last verified result (Python side)

```
E_total[0]   = 2.222288 eV
E_total[end] = 2.222231 eV
drift        = -5.712e-05 eV   (~26 ppm of total — Verlet symplectic noise)
```

End-to-end comparison validation (against a synthesized MATLAB
stand-in with legacy constants):

```
positions  scaled-rel ≤ 1.27e-7   (~127 ppm)
velocities scaled-rel ≤ 4.16e-6
E_kin      scaled-rel ≤ 1.18e-4   (~118 ppm, dominated by EV literal)
E_pot      scaled-rel ≤ 1.01e-7
E_total    scaled-rel ≤ 3.65e-6
mass_kg    scaled-rel ≤ 6.5e-10   (CODATA U vs rounded u)
E_dissip,  t_ps                       exactly 0
VERDICT: PASS
```
