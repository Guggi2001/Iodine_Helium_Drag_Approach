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
- Per-step temperature diagnostic
  ``[<T'/T>_actual, <T'/T>_from_mass_ratio, <theta_lab>_rad]``
  averaged over the colliding atoms (legacy MATLAB
  ``diagnostic_array``, ``vmi_sim_3d_ion_propa.m:683``)

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
    baoab_propagation_step,   # drag-branch sibling (Slice 4)
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
5. **Collision kinematics**: `apply_collision(..., return_diagnostics=True)`
   updates velocities for colliding atoms, returns ΔE per atom, and
   exposes the COM-frame and lab-frame (post-smearing) cosines, mass
   ratio, and pre/post-collision energies needed to build the legacy
   temperature diagnostic. The recipe
   ``temperature_diagnostic_from_collision`` uses the lab-frame cosine
   for the angle column (matches MATLAB
   ``vmi_sim_3d_ion_propa.m:561`` where ``theta = acos(COStheta(b))``
   is built from the lab cosine, not the COM cosine) and reduces to a
   3-element row written into the new ``IonStepState`` field
   ``temperature_diagnostic`` (and later into
   ``IonCheckpoint.temperature_diagnostic`` -- see schema v5).
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

**With** mass attachment, recomputing
``E_kin = ½ m_new v²`` after a 4-amu helium atom attaches at the
atom's current velocity overstates the true kinetic energy of the
ion+helium system by ``½ Δm v²``: the helium contribution to mass
is added but the corresponding kinetic-energy bookkeeping isn't
provided by the model. Mirroring the legacy MATLAB diagnostic at
``vmi_sim_3d_ion_propa.m:762``, the step function now accumulates a
correction term per atom

```
E_mass_attach_defect[t+1] = E_mass_attach_defect[t]
                          - ½ (m_new − m_old) · v_post² · 100²/eV
```

(``v_post`` in Å/ps, mass in kg, factor ``100²/eV`` for the unit
conversion). Adding this term to the per-side total gives the
conservation invariant

```
E_kin + E_pot + E_dissip + E_mass_attach_defect ≈ const
```

modulo Verlet symplectic drift. The diagnostic is exposed as
``IonStepState.E_mass_attach_defect_eV`` (shape ``(2N,)``) and
persisted in ``IonCheckpoint.E_mass_attach_defect_eV`` (shape
``(2N, T)``) in schema v4.

## Temperature-diagnostic capture

``IonStepState`` carries an optional field
``temperature_diagnostic: np.ndarray | None`` of shape ``(3,)``:

```
[<T'/T>_actual, <T'/T>_from_mass_ratio, <theta_lab>_rad]
```

It is set from the ``CollisionDiagnostics`` returned by
``apply_collision`` via ``temperature_diagnostic_from_collision``.
On steps where no atom collided, the value is an all-NaN ``(3,)``
array. On states reconstructed from a checkpoint column (e.g. via
``ion_state_from_checkpoint_column``), the field is ``None`` because
the per-step measurement is not stored per atom.

The driver (``simulation/ion.py``) writes the row into
``IonCheckpoint.temperature_diagnostic[stored_step_idx, :]`` only when
that step is stored, so the array shape matches ``(num_stored_steps, 3)``
-- the same downsampling as MATLAB's
``diagnostic_array(1:reduction_timesteps:end, :)`` in
``vmi_sim_3d_ion_propa.m:883``.

## Drag branch (Slice 4): `baoab_propagation_step` + `_check_drag_scope`

The drag-model port adds a **third sibling** per-step function alongside
`neutral_propagation_step` and `ion_propagation_step`. The driver
(`simulation/ion.py`) dispatches once per run on `cfg.drag_coefficients is not
None`: present → drag path; absent → the collision path above, **bit-identical
and uncalled-for-drag**. The collision step is *not* rewritten in place — it
stays single-purpose so the hard-sphere regression surface is preserved.

### Shared scaffolding (behavior-preserving §3.1 lifts)

Three physics-free helpers were extracted from the collision step and are now
shared by both paths (CLAUDE.md rule 1), each guarded by
`test_ion_propagation_step.py` staying green:

- `_depth(x, y, z, droplet_radii)` — per-atom `r − r_droplet`.
- `_E_kin_eV(mass_kg, vx, vy, vz)` — `½ m v²` in eV.
- `_E_pot_per_atom(depth, E_pot_per_pair, cfg)` — ion-droplet binding + half
  partner Coulomb.

### `baoab_propagation_step(state, *, step, cfg, droplet_radii)`

The Tier-0 drag analog of `ion_propagation_step`: symmetric, but with the
*middle* replaced — **no** collision sampling, **no** `apply_collision`, **no**
mass attachment. The conservative B/A kicks and the dissipative O-step both live
inside the pre-built BAOAB closure `step` (built by the driver via
`make_ion_baoab_step`); this function does thin per-step accounting only:

1. `(pos', vel', E_pot_per_pair, ΔE_dissip) = step(pos, vel, cfg.dt_ion)`.
2. `depth = _depth(...)`; `E_kin = _E_kin_eV(...)` (fixed mass);
   `E_pot = _E_pot_per_atom(...)`.
3. **Dissipation eV conversion** — `ΔE_dissip` arrives from the stepper in
   `amu·Å²/ps²` (Slice 2's pure-mechanical handoff). It is converted with the
   *baseline idiom* (amu→kg via `U`, Å/ps→m/s via `100`, J→eV via `EV`):

   ```
   dE_dissip_eV = ΔE_dissip * U * 100**2 / EV
   ```

   the same conversion path as the mass-attach defect, so `E_dissip_eV` stays in
   a consistent eV with `E_kin_eV`/`E_pot_eV`. It accumulates into the carried
   `E_dissip_eV`.
4. **Tier-0 fills:** `mass_kg` unchanged (fixed mass), `E_mass_attach_defect_eV`
   and `number_of_collisions` carried unchanged (both 0 under the drag branch),
   `temperature_diagnostic` set to an all-NaN `(3,)` sentinel (no collisions to
   diagnose).

The §2.9 energy invariant reduces to `E_kin + E_pot + E_dissip ≈ const`, and
closure is **tight** (the O-step books dissipation analytically — exact — so the
only drift is the conservative Verlet drift).

### `_check_drag_scope(cfg, initial_mass_kg)`

The drag-branch analog of `_check_scope`. `_check_scope` demands collision mode
3, which is irrelevant under drag; this instead asserts the **Tier-0 runnability
envelope** (distinct from `check_drag_config`, which validates config *internal
consistency*). It raises `NotImplementedError` when:

```
noise_form != "none"            (active Langevin noise is Tier 3)
mass_scenario != "fixed"        (mass dynamics is Tier 1)
drag_form != "linear_cubic"     (other forms NotImplemented in drag.py)
effusive_dynamics / single_charge_ionization_allowed / additional_droplet_charges > 0
```

**Realized-mass trip-wire (Slice 4 fix, §6.5).** The final check reads the
**realized** initial ion mass `initial_mass_kg` — the array the stepper will
actually integrate, passed by the driver as `ckpt.mass_kg` *downstream* of the
`build_initial_ion_state` `m_eff` override — and raises if

```
max| initial_mass_kg / U − cfg.m_eff_amu |  >  _MASS_COEFFICIENT_CONSISTENCY_TOL_AMU  (~8 amu)
```

It reuses the §6.5 mass-insensitivity band (no new tolerance constant) and reads
the *realized* array, **not** the config field — reading the config would merely
echo what the override set and guard nothing. A ~127-amu mass (76 amu below
`m_eff`) trips it loudly. This is the trip-wire half of the defense-in-depth
pairing with the `build_initial_ion_state` override.

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
