# Ion driver cross-reference — stochastic forced events (driver plumbing)

CLAUDE.md ion-stage validation target **5** (collision/statistical
behavior at the **driver level**, after the deterministic targets 1–4
are stable).

## What is tested

The **driver bookkeeping around** the hard-sphere collision kernel
in the full ion driver — specifically:

1. `ion_propagation_step` calls the collision sampler with the right
   per-step inputs (previous-step distance, depth, E₀, σ).
2. Velocity-dependent cross section is integrated correctly into the
   loop (`sigma_dependent_on_v = True`, `sigma_ion_exponent = -2`).
3. Mass attachment fires **after** the elastic-scattering kinematics.
4. The post-attachment `E_kin = ½ m_new v²` is recomputed with the
   new mass; the kinetic-energy bookkeeping closes via the
   `E_mass_attach_defect_eV` term added in schema v4 (Step 11e).
5. Per-step storage records `mass_history_kg`, cumulative
   `number_of_collisions`, `E_dissip_eV`, and `E_mass_attach_defect_eV`
   into the right checkpoint columns.

## What is *not* tested here

- The collision kernel itself (impact-parameter sampling, COM↔lab
  transform, ΔE/E₀ distribution): already covered statistically by
  `scripts/collision_comparison_test/` with 1000 atoms × 200 steps.
- A full single-pulse stochastic trajectory: explicitly out of scope
  per CLAUDE.md.
- Mode-1 / Mode-2 collision branches, Gaussian angle smearing,
  charged-droplet model, single-charge ionization: out of project
  scope.

## Approach: forced events, no RNG injection

Trying to compare full trajectories cross-language would require
identical Python/MATLAB RNG streams, which they don't share. Instead
the inputs are tuned so that the **stochastic events** are
deterministic:

- `geometric_scattering_crosssection_Iplus = 1e6 Å²` with
  `sigma_dependent_on_v = True, sigma_ion_exponent = -2` →
  `p_scatter = prev_dist · σ · ρ_droplet ≈ 58` ≫ 1 every step,
  so any uniform draw `trial < p_scatter` ⇒ collision fires.
- `mass_attach_probability = 1.0` → any uniform draw < 1 ⇒
  attachment fires after every collision.

Event **counts** and **mass history** are therefore identical on both
sides regardless of which RNG samples come out. Post-collision
**velocities** and the `E_dissip` / `E_mass_attach_defect`
accumulations remain RNG-dependent and are explicitly **not**
gated cross-language; only their per-side invariants are.

## Files

| File | Role |
|---|---|
| `inputs.json` | Hand-crafted state. 1 molecule, 5 ion steps, atoms deep inside R=10 Å droplet (depth ≈ −8 Å), moderate `\|v\| = 3 Å/ps` so E_kin stays well above the Landau cutoff for the entire run. |
| `export_python_forced.py` | Calls `build_initial_ion_state` + loops `ion_propagation_step`. All stochastic flags driven entirely by cfg — no library edits. Per-step `b_collision` and `b_attach` recovered by differencing `number_of_collisions` and `mass_history_kg` between consecutive states; `sigma_used` recomputed via `velocity_dependent_cross_section`. |
| `export_matlab_forced.m` | Inline reimplementation of `vmi_sim_3d_ion_propa.m` lines ~300–770 (frog_step_ion + Mode-3 sampling + scattering kinematics + mass attachment + defect bookkeeping). Uses legacy-rounded constants. No dependency on legacy globals/`inputfile`. |
| `compare_forced.py` | Three groups (CROSS_MATCH / PER_SIDE_INVARIANTS / NOT_COMPARED) with separate verdicts. |
| `python_forced.csv` / `matlab_forced.csv` | Generated outputs (committed for reference). |

## How to run

```bash
# Python
python scripts/cross_reference/ion_stochastic_forced/export_python_forced.py

# MATLAB
matlab -batch "cd('scripts/cross_reference/ion_stochastic_forced'); export_matlab_forced"

# Compare
python scripts/cross_reference/ion_stochastic_forced/compare_forced.py
```

## Quantities compared and verdict groups

### CROSS_MATCH — must agree across MATLAB and Python

| Quantity | Tolerance |
|---|---|
| `t_ps` | exact (1e-12) |
| `number_of_collisions` (cumulative) | exact integer |
| `b_collision` per step (derived) | exact integer |
| `b_attach` per step (derived) | exact integer |
| `mass_kg` per (step, atom) | 1e-30 kg (CODATA U vs rounded u) |
| `sigma_used_A2` at t=0 | scaled-relative ≤ 5×10⁻⁴ |
| `depth_A` at t=0 | scaled-relative ≤ 5×10⁻⁴ |

`sigma_used` and `depth` are checked **only at t=0**, when both sides
read positions/velocities directly from `inputs.json`. Later steps
depend on RNG-driven post-collision velocities and would diverge.

### PER_SIDE_INVARIANTS — must hold independently on each side

- `mass_history_kg` non-decreasing per atom.
- `b_attach` only where `b_collision` is True (no attachment without
  a preceding collision).
- `E_dissip_eV` non-decreasing per atom.
- Conservation:
  `(E_kin + E_pot + E_dissip + E_mass_attach_defect)[end] − [0]`
  drifts by less than 1% of the initial total energy. The
  `E_mass_attach_defect_eV` term is what closes the invariant when
  helium attaches; it was added to the Python checkpoint in Step 11e
  precisely to make this gating possible.

### NOT_COMPARED — RNG-dependent, informational only

`vx, vy, vz, E_kin_eV, E_pot_eV, E_dissip_eV,
E_mass_attach_defect_eV` after the first collision step. These will
differ between sides because the two RNG streams produce different
impact parameters, reference directions, and azimuths. The script
prints their max abs / scaled-relative diffs for inspection but does
not gate the verdict on them.

## Last verified result

End-to-end comparison validation (against a synthesized MATLAB
stand-in with legacy constants and a deliberately different RNG
seed):

```
CROSS_MATCH:
  t_ps                          0           PASS  (exact)
  number_of_collisions          0           PASS  (forced events)
  b_collision                   0           PASS  (forced events)
  b_attach                      0           PASS  (forced events)
  mass_kg                       1.5e-34     PASS  (constants rounding)
  sigma_used_A2 (step 0)        0           PASS
  depth_A      (step 0)         0           PASS

PER_SIDE_INVARIANTS:
  python: total drift  +0.0145%             PASS  (< 1% threshold)
  matlab: total drift  +0.0134%             PASS
  both:   mass-monotonic, attach-implies-collide, E_dissip-monotonic

NOT_COMPARED (informational):
  vx_Aps, vy_Aps, vz_Aps        large       (RNG-dependent, ignored)
  E_kin/E_pot/E_dissip/E_defect large       (RNG-dependent, ignored)

VERDICT: PASS
```

## Relationship to Step 11e (`E_mass_attach_defect_eV`)

This script's per-side conservation invariant **requires** the v4
schema field `E_mass_attach_defect_eV`. Without it, the recomputed
`E_kin = ½ m_new v²` after each helium attachment would overstate
the true kinetic energy by `½ Δm v²` and the conservation invariant
would drift far beyond the 1% threshold even with no integrator
error. Designing this comparison is what surfaced the missing field
in the original Python port; see `migration_log.md` Step 11e for the
fix.
