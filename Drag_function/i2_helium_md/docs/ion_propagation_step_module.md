# The `ion_propagation_step.py` module

## What problem does this file solve?

This module implements the **per-step physics body** of the ion
propagation. One call advances the entire system by one `dt_ion`
including:

- Leapfrog integration over Coulomb partner force + ion-droplet potential
- Mode-3 hard-sphere collision sampling with **velocity-dependent**
  cross section
- Elastic-scattering kinematics for collisions
- Mass attachment (helium atoms occasionally stick to ions)
- Energy bookkeeping (E_kin, E_pot, cumulative E_dissip,
  cumulative collision count)

It is the analogue of `propagation_step.py` for the ion stage. Both
share the design principle that the step function is **pure** — no
mutation of inputs, no side effects, no I/O. The driver
(`simulation/ion.py` — Step 11d) handles orchestration.

## Position in the dependency chain

```
IonStepState (input)  ──────────────► IonStepState (output)
                                      ▲
   ion_propagation_step               │
   ├── make_ion_step (leapfrog)       │
   │   └── partner_interaction_ion    │ pure function
   │       ion-droplet potential      │
   ├── velocity_dependent_cross_section
   ├── sample_collision_events (Mode 3)
   ├── apply_collision (elastic kinematics)
   └── mass attachment (rng.uniform < p)
```

## Public API

```python
from i2_helium_md.simulation.ion_propagation_step import (
    IonStepState,
    ion_propagation_step,
)

new_state = ion_propagation_step(
    state,                                      # IonStepState
    cfg=cfg,                                    # SimConfig
    droplet_radii=droplet_radii_angstrom,       # (2N,)
    charge=charge,                              # (2N,) -- all 1.0 in our scope
    prev_distance_angstrom=prev_distance,       # (2N,) or None for first step
    rng=rng,                                    # np.random.Generator
)
```

`IonStepState` is a frozen dataclass holding the per-atom dynamic
quantities. The big difference from `NeutralStepState` is that
**`mass_kg` is part of the state** — helium attachment changes it
during the run. `droplet_radii` and `charge` stay constant and are
passed to the function as separate arguments.

## Step sequence

1. **Leapfrog**: build `make_ion_step` closure with current mass,
   integrate one dt → new positions, new velocities, per-pair
   Coulomb potential.
2. **Depth**: `r_new − droplet_radius` per atom.
3. **Cross section**: per-atom `σ = σ_0 · v^exponent` if
   `cfg.sigma_dependent_on_v`, else constant `σ_0`.
4. **Collision sampling**: Mode 3, using the **previous step's**
   distance traveled (passed by the driver). On the first step,
   `prev_distance=None` and no collisions can occur.
5. **Collision kinematics**: `apply_collision` updates velocities
   for colliding atoms and returns ΔE per atom.
6. **Mass attachment**: `mass_attach_trial < p_attach` AND
   `b_collision` → mass += 4 amu (one He atom).  Random number
   drawn for ALL atoms (matches MATLAB rng pattern).
7. **Energy diagnostics**: `E_kin` uses **NEW** mass (matches MATLAB
   line 761). `E_pot` = ion-droplet + half-pair Coulomb. `E_dissip`
   accumulates ΔE. Collision count accumulates.
8. **Return** new `IonStepState` with `time_ps += dt_ion`.

## Mass-as-state design (Option C)

In the design discussion, three options were considered for handling
the changing mass:

- **A** rebuild step closure every iteration (simple but slow)
- **B** make `make_ion_step` accept mass as a runtime argument (clean
  but requires changing leapfrog API)
- **C** carry mass in `IonStepState`, rebuild closure inside the pure
  step function

We chose **C**. The closure overhead is small (a single function
construction per step, no per-atom Python loops) and it keeps the
leapfrog API unchanged.

## Energy conservation

A regression test in `test_ion_propagation_step.py
::TestEnergyConservation::test_drift_small_without_attachment` runs
30 steps with **mass attachment disabled** and asserts that

```
|ΔE_total| / |E_total_init| < 0.5%
```

Without attachment we observe ~0.002% drift over 50 steps — at the
leapfrog symplectic-error limit (~ppm per step).

**With** mass attachment, drift is non-zero by design: when a
helium atom attaches, its kinetic energy isn't tracked
(approximation). The legacy MATLAB has an `E_mass_attach_defect`
diagnostic for this; we omit it for now since it's purely a
diagnostic and not used in the final output.

## Out-of-scope branches

`ion_propagation_step` calls `_check_scope(cfg)` and raises
`ValueError` if cfg requests:

- `hard_sphere_collision_mode != 3`
- `effusive_dynamics = True`
- `single_charge_ionization_allowed = True`
- `additional_droplet_charges > 0`

The legacy MATLAB also has a `relative_energy_loss_ion` flag (an
alternative energy-loss model). It is not in our `SimConfig`
because both production input scripts leave it disabled. If we add
the field later, mirror the existing scope-check pattern.

## Testing

- `test_ion_propagation_step.py` — 24 unit tests across 6 classes
  (TestApi, TestFirstStep, TestReproducibility, TestEnergyBookkeeping,
   TestMassAttachment, TestEnergyConservation, TestScopeChecks,
   TestVelocityDependentSigma).
- `smoke_test_ion_propagation_step.py` — 30 sandbox checks.
