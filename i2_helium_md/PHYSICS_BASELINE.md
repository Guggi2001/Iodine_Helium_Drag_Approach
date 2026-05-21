# Physics Baseline — `i2_helium_md`

**Date frozen:** 2026-05-20
**Repository state:** main @ `49e48af` ("19.05 finished post processing finished documentation")
**Purpose.** This document is the frozen reference for the **current** physics
implementation of the I₂ / I₂⁺-in-He-droplet molecular-dynamics code, captured
immediately before the planned replacement of the hard-sphere collision model
in `i2_helium_md/physics/collisions.py` by a TDDFT-calibrated drag-force model
for I⁺ in a helium bubble.

The drag-model port is the only scoped exception to the "do not change
collision physics" rule in `CLAUDE.md` ("Forbidden Without Explicit User
Approval" → "Active exception (scoped, time-limited)"). Everything else in this
document — neutral propagation, ion propagation outside the collision step,
checkpoint schema, constants, RNG draw order, presets, post-processing surface
— is **out of scope** for the upcoming work.

All line references are relative to the verified file revisions read on
2026-05-20 and quote the Python source — not the MATLAB legacy.

---

## 1. System overview

Single-pulse ionization of I₂ embedded in a He nanodroplet, followed by ion
dynamics on the dication (and optional singly-ionized) potential surfaces.
Neutral propagation thermalizes/relaxes the cluster; at t = t_neutral_end the
ion stage takes the final neutral state, switches charges and potentials, and
propagates the Coulomb-driven dissociation for 20 ps while particles collide
elastically with bubble helium atoms (the model being replaced).

Layered architecture:

```
            physics/constants.py
                    │
            physics/potentials.py            (Morse, droplet, Coulomb consts)
                    │
            physics/interactions.py          (pair force assembly, 2N layout)
                    │
            physics/leapfrog.py              (velocity-Verlet, model-agnostic)
                    │
   physics/collisions.py    ──┐              (the swap point for drag model)
                              │
            simulation/propagation_step.py   (neutral per-step orchestration)
            simulation/ion_propagation_step.py (ion per-step orchestration)
            simulation/{neutral,ion}.py      (drivers)
            simulation/checkpoint.py         (NeutralCheckpoint v2, IonCheckpoint v5)
                    │
            postprocess/*                    (diagnostics, comparison, plots)
```

The 2N layout convention pervades every array: for `N` molecules,
indices `0..N-1` are atom-1 of each pair and `N..2N-1` are the partner atom-2.
This is shared with MATLAB.

Internal units: **Å** for positions, **ps** for time, **Å/ps** for velocity
(1 Å/ps = 100 m/s), **eV** for energy (SI used internally only where forces
are computed), **kg** for mass.

---

## 2. Constants and unit conversions

File: `i2_helium_md/physics/constants.py`. Single source of truth for every
fundamental constant the simulation uses.

| Symbol / name                              | Value                          | Line | Notes |
|--------------------------------------------|--------------------------------|------|-------|
| `E_CHARGE`                                 | −1.602176634e−19 C             | 33   | exact (2019 SI); sign matches MATLAB |
| `EPSILON_0`                                | 8.8541878188e−12 F/m           | 34   | CODATA 2022 |
| `EV`                                       | 1.602176634e−19 J              | 36   | exact |
| `U`                                        | 1.66053906892e−27 kg           | 35   | atomic mass unit, CODATA 2022 |
| `K_B`                                      | 1.380649e−23 J/K               | 37   | exact |
| `HC`                                       | 1239.841984 eV·nm              | 38   | not used in dynamics path |
| `EV_PER_WAVENUMBER`                        | 1/8065.543937 eV / cm⁻¹        | 41   | spectroscopic helper |
| `BULK_DENSITY_HELIUM`                      | 0.0219 atoms/Å³                | 48   | liquid He bulk |
| `DENSITY_DROPLET`                          | 0.8 × 0.0219 ≈ 0.01752 atoms/Å³| 49   | effective density used by collision sampler |
| `MASS_I_AMU`                               | 127.0 u                        | 92   | iodine mass |
| `EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2`      | 1.602176634e−23                | 131  | single conversion shared by `leapfrog.py` and `interactions.py` |

Derived helpers in the same file:

- `droplet_radius_bulk_angstrom(N) = (3 N / (4 π · 0.0219))^(1/3)` (line 52)
- `coulomb_energy(r_Å)` returns Joules; `coulomb_velocity(r_Å, m_kg)` returns
  `sqrt(E_C / m)` (lines 98, 103).

### Things that are NOT in `constants.py`

- **He scatterer mass.** Set as `SimConfig.scatter_mass_neutral_amu = 4.0`
  (config.py:106) and `scatter_mass_ion_amu = 4.0` (config.py:107). The drag
  port should keep mass as a configurable, not hard-code it.
- **Coulomb prefactor.** Hardcoded literal `14.39964548 eV·Å` inside
  `ion_interaction_potential` at `interactions.py:114`. Equals `e²/(4πε₀)`
  expressed in eV·Å. Documented but not exported. Multiplied by
  `cfg.E_coulomb_scale` at the call site.

### Intentional drift from legacy MATLAB

The Python constants are full CODATA 2022 / SI 2019 values. The legacy MATLAB
`physical_constants.m` used 4-significant-figure approximations
(`eV = 1.602e-19`, etc.). This was a ~100 ppm error in the MATLAB values; the
Python correction introduces a ~10 ppm trajectory drift relative to MATLAB
that is regression-guarded but **not** to be reverted (see
`CLAUDE.md` → "Project Quality Principles" rule 10).

---

## 3. Potentials

File: `i2_helium_md/physics/potentials.py`.

### 3.1 Morse — neutral I₂ ground state ("X-state")

`morse_X(r, cfg)` — `potentials.py:163`.

```
U_X(r) = D_e · (1 − exp(−a · (r − R_e)))²    [eV, r in Å]
       − 0.9 · exp(−(r − 9)² / (2 · 0.3²))   if cfg.Xdip_active  (line 184)
```

Spectroscopic parameters in `I2_X_STATE` (line 148):

| D_e        | ω_e         | ω_e·x_e    | R_e       | a (1/Å) derived from `_morse_a()` (line 41) |
|------------|-------------|------------|-----------|---------------------------------------------|
| 1.556 eV   | 214.5 cm⁻¹  | 0.65 cm⁻¹  | 2.666 Å   | reduced mass μ = 127/2 × U                  |

Source: *J. Chem. Phys.* **107**, 9046 (1997), `doi:10.1063/1.475194`.

The `Xdip` Gaussian correction at r ≈ 9 Å is an empirical bump that mirrors
the avoided-crossing physics in the experimental fit; enabled by default
(`cfg.Xdip_active = True`, config.py:92).

### 3.2 Morse — I₂⁺ electronic states

Four ionic-state Morse curves in `I2PLUS_STATES` (line 151):

| state idx | D_e (eV) | ω_e (cm⁻¹) | ω_e·x_e (cm⁻¹) | R_e (Å) | IP_rel (eV) |
|-----------|----------|------------|----------------|---------|-------------|
| 0         | 2.70     | 240.0      | 0.69           | 2.61    | 0.00        |
| 1         | 2.03     | 230.0      | 0.29           | 2.61    | 0.63        |
| 2         | 1.26     | 141.0      | 0.32           | 2.95    | 1.68        |
| 3         | 0.56     | 117.0      | 0.38           | 2.95    | 2.44        |

`I2PLUS_IP_REL` at line 159, `I2PLUS_IP_0 = 9.36 eV` at line 160.

`morse_I2plus_state_select(r, state_ids)` at line 209 evaluates per-molecule
the chosen state's curve, shifted by `−V_morse(2.666; state) + IP_0 + IP_rel`
to give absolute ionic potential energy. Only active when
`cfg.single_charge_ionization_allowed = True` (default `False`; out of scope
for the drag port — `ion_propagation_step._check_scope` refuses to run if
this flag is on, see `ion_propagation_step.py:319`).

### 3.3 Droplet confining potential

`droplet_potential(r, steepness, binding_energy, offset=0)` at `potentials.py:73`.

```
V_drop(r) = ½ · (erf((r − offset)/steepness) + 1) · binding_energy   [eV]
```

The radial argument `r` is **depth relative to the droplet surface**:
`r_atom − droplet_radius`. Inside the droplet (`r ≪ 0`) V → 0; outside
(`r ≫ 0`) V → `binding_energy`.

The analytical derivative `droplet_force(r, …)` at `potentials.py:106`:

```
dV/dr = binding_energy / (steepness · √π) · exp(−((r−offset)/steepness)²)
```

This **replaces the MATLAB finite-difference form** for the droplet
contribution and is the only analytical force in the codebase — the I-I and
ion-ion partner forces still use the MATLAB-style finite difference for
byte-level compatibility (see §4).

Defaults from `SimConfig`:

| Field                            | Default     | Line | Used for                                |
|----------------------------------|-------------|------|-----------------------------------------|
| `potential_steepness`            | 14.2 Å      | 86   | atom-level erf width                    |
| `potential_steepness_molecule`   | 14.3324 Å   | 87   | molecule-level (post-processing only)   |
| `binding_energy_I_atom_K`        | 318.43 K    | 88   | neutral I atom; → eV via line 149       |
| `binding_energy_molecule_K`      | 573.3 K     | 89   | I₂ molecule; → meV via line 154         |
| `binding_energy_I_ion_eV`        | 0.3 eV      | 117  | I⁺ ion; consumed directly               |

The K→eV conversion uses `K_B / EV` (config.py:150).

### 3.5 Binding energies — physical meaning and drag-port implications

The three binding-energy fields above are physically distinct quantities,
not three names for the same thing. They enter the dynamics at different
sites and have different magnitudes:

| Field                                | Value               | Where it enters                                                                       |
|--------------------------------------|---------------------|---------------------------------------------------------------------------------------|
| `binding_energy_I_atom_K`            | 318.43 K ≈ 27.4 meV | `_droplet_acceleration` with `use_ion_binding=False` (`leapfrog.py:178`); neutral I |
| `binding_energy_I_ion_eV`            | 0.30 eV             | `_droplet_acceleration` with `use_ion_binding=True`  (`leapfrog.py:177`); I⁺        |
| `binding_energy_molecule_K`          | 573.3 K ≈ 49.4 meV  | **post-processing only** (`boltzmann_overlay.py`); not consumed by the dynamics loop  |

The erf-based droplet potential `V_drop(r) = ½ · (erf((r − offset)/steepness)
+ 1) · E_b` is **asymptotic to 0 deep inside the droplet** (r ≪ 0) and
**asymptotic to E_b outside** (r ≫ 0). E_b is the energy a particle pays
to leave. The ion binding (300 meV) is roughly **10× the neutral atom
binding (27 meV)** — the I⁺ sits in a much steeper trap because of
electrostriction / dielectric attraction by the surrounding helium.

Selection logic in `leapfrog.py:_droplet_acceleration` (lines 174-185):

```python
binding = (
    cfg.binding_energy_I_ion_eV if use_ion_binding
    else cfg.binding_energy_I_atom_eV
)
dU_dr = droplet_force_fn(depth, steepness=cfg.potential_steepness,
                         binding_energy=binding)
```

`use_ion_binding=False` for `_neutral_accel_fn` (`leapfrog.py:230-232`);
`True` for `_ion_accel_fn` (`:264-266`). So the moment the simulation
transitions from neutral to ion stage, the droplet well **deepens by ~10×
under identical positions/velocities** — this is the dominant geometric
change for the ion-stage initial condition.

Sanity guard. `cfg.validate()` (config.py:174-182) compares
`E_min_eV` (Landau cutoff) against `binding_energy_I_atom_eV` and emits
`"all neutrals will escape!"` `RuntimeWarning` if `E_min > binding`. With
the defaults (`v_limit_m_per_s = 40`, `MASS_I_AMU = 127`),
`E_min ≈ 1.06 × 10⁻⁷ eV` — far below the 27 meV neutral binding, so the
guard is dormant in production.

### 3.5.1 Why the drag port has to revisit binding

The 0.30 eV ion binding in `single_pulse_N2000` was **co-fit with the
hard-sphere cross-section and the velocity-exponent `−2`** to match the
HeDFT trajectories. Evidence that it is a fit parameter, not a fixed
physical constant:

- `single_pulse_N2000_18Angst` lowers `binding_energy_I_ion_eV` from
  0.30 eV to **0.05 eV** (presets.py:89) for the 18 Å droplet — a 6× drop
  for the same I⁺ species, justified by the different geometry.
- `single_pulse_droplet_distribution` keeps 0.30 eV but applies
  `E_coulomb_scale = 0.8` (presets.py:112) — i.e. the binding and the
  Coulomb screening are tuned together.

The drag-force model dissipates energy continuously rather than at
discrete scattering events; the effective escape barrier seen by an ion
on its way out couples to the **drag profile**, not just to V_drop. So
both `binding_energy_I_ion_eV` and the new drag parameters must be
**re-calibrated jointly** against the same HeDFT references. Treat the
current 0.30 / 0.05 eV values as starting points, not as ground truth.

Neutral atom binding (`binding_energy_I_atom_K = 318.43 K`) and droplet
steepness are likely **safe to keep** for the drag port — they govern
the neutral stage, which is unaffected by the collision-model swap.

### 3.4 Coulomb

`ion_interaction_potential(dr, q1, q2, cfg, *, state_ids=None)` at
`interactions.py:73`. For the canonical I⁺/I⁺ case (q₁ = q₂ = 1):

```
U_Coulomb(r) = cfg.E_coulomb_scale · 14.39964548 / r        [eV, r in Å]
```

`cfg.E_coulomb_scale = 1.0` by default (config.py:122); the
`single_pulse_droplet_distribution` preset overrides this to `0.8`
(presets.py:112) to model partial dielectric screening by the droplet.

If `cfg.single_charge_ionization_allowed = True` (default `False`), an
additional Morse-I₂⁺ term is added to pairs with `q₁ + q₂ == 1`
(interactions.py:122-126). This branch is **out of scope** and rejected by
the ion-step scope check.

---

## 4. Forces and 2N pair geometry

File: `i2_helium_md/physics/interactions.py`.

The integrator (`leapfrog.py`) consumes an `acc_fn(pos) → (acc, E_pot)`
callable, so the bridge from scalar U(r) to per-atom 3D acceleration lives
here. The bridge is identical for the neutral X-state Morse and the
ion Coulomb (+ optional Morse) potentials.

Per-step assembly (`partner_interaction_neutral` at `:244`,
`partner_interaction_ion` at `:284`):

1. `_pair_geometry` (`:151`) — compute pair distance `dr` (shape N) and unit
   vector `dr_unit` pointing from atom-2 to atom-1.
2. `_force_from_potential_fd` (`:175`) — finite-difference force
   `F = (U(r) − U(r+h)) / h`, with `h = 1e-4 Å`. **Kept for byte-level
   MATLAB compatibility**; the only analytical force in the codebase is the
   droplet contribution (§3.3).
3. `_acceleration_from_force` (`:204`) — broadcast scalar force to the 2N
   layout (atom-1 along `+dr_unit`, atom-2 along `−dr_unit`, Newton's 3rd
   law), divide by mass, multiply by `EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2`
   to land in Å/ps².

Public force entry points:

| Function                       | File:Line                | Used by                  |
|--------------------------------|--------------------------|--------------------------|
| `atom_interaction_potential`   | interactions.py:47       | neutral leapfrog branch  |
| `ion_interaction_potential`    | interactions.py:73       | ion leapfrog branch      |
| `partner_interaction_neutral`  | interactions.py:244      | `_neutral_accel_fn`      |
| `partner_interaction_ion`      | interactions.py:284      | `_ion_accel_fn`          |

The per-atom potential-energy convention for ions is `E_pot_per_pair / 2`
duplicated to both atoms (interactions.py:338), matching MATLAB.

---

## 5. Integrator

File: `i2_helium_md/physics/leapfrog.py`.

`velocity_verlet_step(pos, vel, acc_fn, dt)` at `leapfrog.py:74`. Symplectic,
second-order, kick-drift-kick:

```
(ax0, ay0, az0), _      = acc_fn(x0)                              # line 113
x1 = x0 + dt·v0 + 0.5·a0·dt²                                       # line 116-118
(ax1, ay1, az1), E_pot1 = acc_fn(x1)                              # line 121
v1 = v0 + 0.5·(a0 + a1)·dt                                         # line 124-126
return (x1, v1, E_pot1)
```

The function is **agnostic to the force model**: it consumes any callable
matching `(x, y, z) → ((ax, ay, az), E_pot[N])`. The drag-model port should
not need to touch this file.

The two force assemblers:

- `_neutral_accel_fn` at `:217` — droplet (atom binding) + optional I-I Morse.
- `_ion_accel_fn` at `:248` — droplet (ion binding) + ion-ion (Coulomb + opt.
  Morse).

Convenience factories that bind `cfg, mass, droplet_radii, charge,
state_ids` once and return a `step(pos, vel, dt) → (pos', vel', E_pot)`
closure:

- `make_neutral_step(cfg, mass, droplet_radii)` — `:287`
- `make_ion_step(cfg, mass, droplet_radii, charge, state_ids=None)` — `:318`

Note: the ion step closure is **rebuilt inside the ion per-step function**
on every iteration because mass changes via attachment
(`ion_propagation_step.py:188`).

---

## 6. Hard-sphere collision model — **THE BASELINE BEING REPLACED**

File: `i2_helium_md/physics/collisions.py`.

This module is **pure physics**: no `SimConfig` dependency, no I/O. The
neutral and ion drivers pull config values and pass them in as keyword
arguments. This is the swap point for the drag port.

### 6.1 Cross-section

`velocity_dependent_cross_section(v, *, sigma_0, exponent)` at
`collisions.py:149`:

```
σ(v) = σ_0 · v^exponent          [Å², v in Å/ps]
```

Defaults for ions (config.py):

- `geometric_scattering_crosssection_Iplus = 2500.0 Å²` (line 105)
- `sigma_ion_exponent = −2.0` (line 109)
- `sigma_dependent_on_v = True` (line 108) — if `False`, the constant
  `geometric_scattering_crosssection_Iplus` is used (ion-step line 218).

Neutrals always use a fixed cross-section
`geometric_scattering_crosssection_I = 30.0 Å²` (line 104).

At `v = 0` with `exponent < 0`, σ → +∞; the downstream
`trial < p_scatter` test always evaluates True, so an infinitely slow ion
always collides — this is the intended Landau-cutoff-paired behaviour.

### 6.2 Event sampling (Mode 3)

`sample_collision_events(…)` at `collisions.py:204`:

```
p_scatter = distance_travelled_Å · σ · ρ_droplet     (with ρ_droplet = 0.8 × 0.0219)
collide ⇔   uniform(0,1) < p_scatter
        ∧   depth < 0                                (inside droplet)
        ∧   E0_eV ≥ E_min_eV                         (Landau cutoff)
```

`distance_travelled` is the **previous step's** path length per atom (the
driver tracks this between calls — see ion-step §7 and neutral-step §8).

The Landau cutoff comes from `cfg.E_min_eV` (config.py:159):

```
E_min = ½ · (127 u) · (v_limit_m_per_s)² / EV     with v_limit = 40 m/s by default
```

Modes 1 (constant probability per step) and 2 (mean-free-path tracking) from
the MATLAB code are **not ported** (`collisions.py` module docstring; ion-step
`_check_scope` at line 310 rejects mode ≠ 3, neutral-step rejects same at
`propagation_step.py:138`).

### 6.3 Velocity update

`apply_collision(…)` at `collisions.py:305`. Per colliding atom:

1. **Impact parameter** by inverse-CDF: `b/R = √u`, `u ∼ Uniform(0,1)`
   (line 433).
2. **COM-frame scattering cosine** (hard sphere):
   `cos θ_COM = 2 (b/R)² − 1` (line 437).
3. **Mass ratio** `ρ = m_atom_amu / m_scatterer_amu` (line 445).
4. **Post-collision energy** (line 447):
   ```
   E₁ = E₀ · (1 + 2 ρ cos θ_COM + ρ²) / (1 + ρ)²
   ```
5. **Lab-frame angle** (line 458):
   ```
   cos θ_lab = (cos θ_COM + ρ) / √(1 + 2 ρ cos θ_COM + ρ²)
   ```
6. **Optional Gaussian smearing** on `θ_lab` only for colliders, std =
   `cfg.{neutral,ion}_scatter_angle_std_deg · π/180` (default 0 — disabled,
   config.py:111-112; lines 463-472).
7. **Reassemble velocity** with a random orthonormal basis perpendicular to
   the incoming velocity (lines 474-491) and a non-uniform azimuth sample
   `COSBETA = 2·u − 1` (line 494). Non-uniform-phi is a documented latent
   MATLAB quirk — mirrored exactly; an isotropy unit test in
   `tests/physics/test_collisions.py` verifies the resulting 3-D
   distribution is isotropic in the perpendicular plane.

Non-colliders keep their incoming velocity exactly (lines 517-520),
guarding against roundoff.

Returns `(vx_new, vy_new, vz_new, ΔE_eV)` and, with
`return_diagnostics=True`, a `CollisionDiagnostics` namedtuple
(`collisions.py:53`) used by the ion driver to record the legacy
temperature-diagnostic accumulator
`[⟨T'/T⟩_actual, ⟨T'/T⟩_from_mass_ratio, ⟨θ_lab⟩_rad]`
(`temperature_diagnostic_from_collision` at line 90).

### 6.4 What the drag model must reproduce or supersede

Validation surface (see §10). At minimum:

- HeDFT trajectory comparison at R₀ = 9 Å and 18 Å must remain inside the
  current tolerances of `compare_distance` / `compare_velocity_magnitude` in
  `i2_helium_md/postprocess/compare_trajectories.py`.
- The final-velocity histogram from the experimental-condition droplet
  preset must keep its agreement with the VMI reference under
  `i2_helium_md/postprocess/velocity_distribution.py`.
- Energy bookkeeping (§9) must remain balanced.

---

## 7. Ion propagation — per-step sequence

File: `i2_helium_md/simulation/ion_propagation_step.py`, function
`ion_propagation_step` at line 141.

Scope check (`_check_scope`, line 308) refuses to run if any of these are on:
`hard_sphere_collision_mode ≠ 3`, `effusive_dynamics`,
`single_charge_ionization_allowed`, `additional_droplet_charges > 0`.

Step body:

| # | Action                                                         | Line range |
|---|----------------------------------------------------------------|------------|
| 1 | Build per-step ion closure (mass may have changed); leapfrog   | 188-193    |
| 2 | Per-atom depth `r1 − droplet_radii`                            | 200-201    |
| 3 | E₀ at post-leapfrog velocity, eV                               | 206-208    |
| 4 | Per-atom σ via `velocity_dependent_cross_section` or constant  | 211-218    |
| 5 | `sample_collision_events` using **previous** step distance     | 220-232    |
| 6 | `apply_collision` with `return_diagnostics=True`               | 237-247    |
| 7 | Mass attachment: `Bernoulli(p_attach) ∧ collided` → +4 amu     | 249-256    |
| 8 | E_kin at new mass and velocity                                 | 261-262    |
| 9 | E_pot = droplet (ion binding) + ½·partner Coulomb per atom     | 264-273    |
| 10| Cumulative: E_dissip, n_collisions, E_mass_attach_defect       | 275-289    |

The defect term (`E_mass_attach_defect_eV`) compensates the spurious
kinetic-energy increase when 4 amu attaches at non-zero velocity:

```
ΔE_defect = − ½ · (m_new − m_old) · v_post²       (line 288)
```

so that `E_kin + E_pot + E_dissip + E_mass_attach_defect` is conserved
modulo Verlet drift. Mirrors `vmi_sim_3d_ion_propa.m:762`.

The per-step temperature diagnostic (3,) is recorded to
`IonStepState.temperature_diagnostic` (line 301) and persisted to the
checkpoint as the `(T, 3)` `temperature_diagnostic` array.

### 7.1 Mass attachment — extended detail

Mass attachment is the channel by which the I⁺ ion **gains weight over
time** by sticking individual He atoms (4 amu each) as it moves through
the droplet. It is **physically distinct** from the elastic collision that
sets the angular and energy update — although in the current model the
two share an RNG-coupled gate.

#### 7.1.1 The Bernoulli gate

`ion_propagation_step.py:249-256`:

```python
mass_attach_trial = rng.uniform(0.0, 1.0, size=n_atoms)
b_mass_attach    = (mass_attach_trial < cfg.mass_attach_probability) & b_collision
new_mass_kg      = state.mass_kg + b_mass_attach * 4.0 * U
```

Rules:

1. Drawn for **every** atom each step (not just colliders) to keep the
   RNG stream deterministic — mirrors `vmi_sim_3d_ion_propa.m:727`.
2. Attachment requires **both** `trial < cfg.mass_attach_probability`
   **and** `b_collision == True`. Non-colliding atoms cannot gain mass
   on a given step.
3. Mass increment is fixed at `4.0 × U` kg per attachment (one He).
   Multiple attachments compound over the run; `mass_history_kg[:, t_id]`
   records the trajectory.
4. The He mass is **hard-coded** at this site (4.0), unlike the
   collision-physics path where `cfg.scatter_mass_ion_amu` is used. This
   is an inconsistency to be aware of for the drag port.

Default probability `cfg.mass_attach_probability = 0.09` per collision
in `single_pulse_N2000`; **0.005** in `single_pulse_N2000_18Angst`
(presets.py:90). The drop with droplet size is empirical — fewer sticky
collisions at larger droplet radius.

#### 7.1.2 Why mass changes mean re-building the leapfrog closure

Because the per-atom mass is dynamic, the ion leapfrog closure
**must be rebuilt every step**. The full sequence at
`ion_propagation_step.py:184-193`:

```python
step_fn = make_ion_step(cfg, state.mass_kg, droplet_radii, charge)
(x1, y1, z1), (vx1, vy1, vz1), E_pot_coulomb_per_pair = step_fn(...)
```

`make_ion_step` (`leapfrog.py:318`) closes over `mass` via `_StepContext`;
that closure is single-mass-snapshot. Re-creating it every iteration is
the chosen "Option C" design (state carries mass; closure rebuilt inside
the pure step function). For the drag port: if drag introduces additional
state (e.g. a memory kernel, a He shell thickness), it should ride along
in `IonStepState` and re-close per step too — do **not** stash drag
state in module-level globals.

#### 7.1.3 The mass-attach kinetic-energy defect

When mass jumps from `m_old` to `m_new = m_old + Δm` at post-collision
velocity `v_post`, the recomputed kinetic energy at line 262

```
E_kin_new = ½ m_new v_post²
```

is spuriously larger than `½ m_old v_post²` by `½ Δm v_post²` —
attaching cold He at the ion's velocity injects fake kinetic energy
because the He atom is treated as if it were already comoving. The
defect term **subtracts that fictitious energy back out** so the global
invariant holds. From `ion_propagation_step.py:286-289`:

```python
mass_diff_kg            = new_mass_kg - state.mass_kg
dE_defect_eV            = -0.5 * mass_diff_kg * (v_post_sq * 100.0 ** 2) / EV
E_mass_attach_defect_new = state.E_mass_attach_defect_eV + dE_defect_eV
```

Conservation invariant (modulo Verlet drift):

```
E_kin + E_pot + E_dissip + E_mass_attach_defect ≈ const
```

`E_mass_attach_defect_eV` is **per atom, cumulative, in eV**, always ≤ 0.
It is the field that allows the energy-balance plots to close in the
presence of mass attachment.

#### 7.1.4 What the checkpoint exposes

`IonCheckpoint` (`simulation/checkpoint.py:140-210`, schema v5):

- `mass_kg (2N,)` — **initial** mass (at t=0 of ion stage), shape 2N.
- `mass_history_kg (2N, T)` — mass per atom per stored step. Strictly
  non-decreasing along T for each atom.
- `mass_final_kg (2N,)` — mass at end of run (= `mass_history_kg[:, −1]`).
- `E_mass_attach_defect_eV (2N, T)` — cumulative defect.
- `number_of_collisions (2N, T) int` — total collisions including
  sticky ones; not separated by collision kind.

There is **no per-step "sticky collision happened" boolean** in the
checkpoint. Reconstruction requires diffing `mass_history_kg` along T.

#### 7.1.5 Coupling that the drag port has to decide on

The current model bundles three behaviours into one Bernoulli step:

- **Energy & angle update** (always when `b_collision`),
- **Sticky vs elastic** (with probability `mass_attach_probability` for
  colliders),
- **He pickup as mass increase** (4 amu added to the ion).

Under a continuous drag model there is no discrete `b_collision` event;
the bundle has to be unpacked. The design choices:

1. **Drop attachment entirely.** Treat the I⁺ as carrying its He shell
   implicitly inside the drag profile. Then `mass_history_kg` becomes
   trivially constant and `E_mass_attach_defect_eV` stays at zero. The
   checkpoint can keep these fields (set to constant initial values) to
   avoid a schema bump.
2. **Keep attachment as an independent Poisson process** with rate
   `λ_attach [1/ps]` driven only by the drag step (e.g. proportional to
   ion speed through the bubble). Schema unchanged; the Bernoulli gate
   moves to a new code path that no longer reads `b_collision`.
3. **Solvation shell as a state variable.** Track shell radius or He
   count separately from the integer-amu mass jumps. Requires a
   `IonCheckpoint` v6 with a new array.

Whichever route is taken, the **energy-balance invariant must remain
closed**. Any new dissipation channel needs either a sibling accumulator
or to fold cleanly into `E_dissip_eV`. See §11.3 → §14.

---

## 8. Neutral propagation — per-step sequence

File: `i2_helium_md/simulation/propagation_step.py`, function
`neutral_propagation_step` at line 99.

| # | Action                                                         | Line range |
|---|----------------------------------------------------------------|------------|
| 1 | `make_neutral_step` closure (mass constant); leapfrog          | 147-152    |
| 2 | Per-atom depth                                                 | 155-156    |
| 3 | E₀ at post-leapfrog velocity                                   | 159-160    |
| 4 | `sample_collision_events` with `σ = cfg.geom_scatter_xsec_I`   | 162-174    |
| 5 | `apply_collision` (no diagnostics return)                      | 177-185    |
| 6 | E_kin                                                          | 188-189    |
| 7 | E_pot = droplet (atom binding) + ½·partner Morse per atom      | 191-197    |
| 8 | Cumulative: E_dissip, `L_droplet` (path length inside droplet) | 199-207    |

No mass attachment in the neutral stage (out of scope per module docstring).

---

## 9. Energy bookkeeping and invariants

Per-atom (2N, T) arrays stored in checkpoints:

| Array                       | Neutral | Ion | Meaning                                              |
|-----------------------------|:-------:|:---:|------------------------------------------------------|
| `E_kin_eV`                  | ✓       | ✓   | ½ m v², per-atom                                     |
| `E_pot_eV`                  | ✓       | ✓   | droplet + ½ partner (Morse for neutral / Coulomb for ion) |
| `E_dissip_eV`               | ✓       | ✓   | cumulative collision energy loss                     |
| `L_droplet_eV_ps`           | ✓       | —   | cumulative path length inside droplet (Å, despite the field name) |
| `E_mass_attach_defect_eV`   | —       | ✓   | accumulated negative compensator for mass-attach E_kin jumps |
| `number_of_collisions`      | —       | ✓   | cumulative count per atom                            |

Invariant: `E_kin + E_pot + E_dissip + E_mass_attach_defect` is conserved
to Verlet drift in the absence of dissipation channels. The drag port adds a
new dissipation channel; it must either reuse `E_dissip_eV` or introduce a
sibling field with its own checkpoint slot (§11.3).

Energy-balance verification scripts (post-processing):

- `scripts/post_processing/plot_neutral_energy_balance.py`
- `scripts/post_processing/plot_ion_energy_balance.py`
- `scripts/post_processing/plot_ion_temperature_diagnostic.py`

### Known intentional differences from MATLAB at t=0

(From `CLAUDE.md` → "Known MATLAB Bugs Not To Reproduce"; do **not** match
MATLAB on these — Python is correct.)

- Neutral-stage `E_pot` at t=0 includes the partner Morse contribution.
- Ion-stage `E_kin` at t=0 includes `vz` with the correct expression.
- Ion-stage `E_pot` at t=0 includes the `z` coordinate and the partner
  Coulomb term.

Regression tests guard these in `tests/physics/test_energy_balance.py` and
`tests/simulation/test_ion_propagation_step.py`.

---

## 10. `SimConfig` — drag-relevant fields

File: `i2_helium_md/config.py`. All tunable parameters are dataclass fields;
no globals. Defaults reproduce
`inputfiles_dft_comparison/single_pulse_N2000.m`.

| Field                                            | Default            | Line | Relevance to drag port |
|--------------------------------------------------|--------------------|------|------------------------|
| `dt_ion`                                         | 0.01 ps            | 55   | drag-force integration timescale |
| `ion_simulation_time`                            | 20.0 ps            | 54   | total ion-stage time |
| `T_particles_K`                                  | 0.4 K              | 73   | He bubble temperature proxy |
| `potential_steepness`                            | 14.2 Å             | 86   | droplet surface width (unchanged by port) |
| `binding_energy_I_ion_eV`                        | 0.3 eV             | 117  | ion-droplet well depth |
| `E_coulomb_scale`                                | 1.0                | 122  | dielectric screening knob |
| `v_limit_m_per_s`                                | 40.0 m/s           | 97   | Landau cutoff feed |
| `hard_sphere_collision_mode`                     | 3                  | 102  | will likely become irrelevant; keep field until removal approved |
| `geometric_scattering_crosssection_Iplus`        | 2500.0 Å²          | 105  | hard-sphere baseline; retained but unused under drag |
| `sigma_dependent_on_v`                           | True               | 108  | same as above |
| `sigma_ion_exponent`                             | −2.0               | 109  | same as above |
| `scatter_mass_ion_amu`                           | 4.0 u              | 107  | He mass; drag model must keep this configurable |
| `mass_attach_probability`                        | 0.09               | 118  | sticky-collision rate; out of scope if drag is collisionless |
| `ion_scatter_angle_std_deg`                      | 0.0                | 112  | angular smearing; drag-model-irrelevant |
| `binding_energy_I_atom_K`                        | 318.43 K           | 88   | unchanged |
| `Xdip_active`                                    | True               | 92   | unchanged |
| `seed`                                           | None               | 39   | RNG seed; preserve draw order |

Derived (read-only properties) at config.py:138-167:
`num_timesteps_neutral`, `v_limit_angstrom_per_ps`,
`binding_energy_I_atom_eV`, `binding_energy_molecule_meV`, `E_min_eV`.

---

## 11. Presets

File: `i2_helium_md/presets.py`. Start from a preset and override with
keyword arguments.

| Preset                                  | R₀ (Å) | N      | σ_I⁺ (Å²) | binding_I_ion (eV) | mass_attach | Other                       |
|-----------------------------------------|--------|--------|-----------|--------------------|-------------|-----------------------------|
| `single_pulse_N2000`           (line 18) | 9.0    | 2000   | 2500      | 0.30               | 0.09        | HeDFT canonical comparison  |
| `single_pulse_N2000_18Angst`   (line 70) | 18.0   | 2000   | 1600      | 0.05               | 0.005       | Fixed 18 Å droplet HeDFT    |
| `single_pulse_droplet_distribution` (line 95) | 2.666 | 8000 | 2500 | 0.30 | 0.09 | Variable droplet sizes, E_coulomb_scale = 0.8, `single_initial_position = False` |

`single_pulse_droplet_distribution` sets `R0_GS_angstrom = 2.666` (the
ground-state I₂ equilibrium distance), not 9 Å — distinct from the HeDFT
preset.

---

## 12. Checkpoint schemas

File: `i2_helium_md/simulation/checkpoint.py`. Both checkpoints are
`@dataclass` and serialized to `.npz` (`allow_pickle=False` on load).

### NeutralCheckpoint — schema version 2 (checkpoint.py:85, version at line 78)

Per-atom shape `(2N,)` static or `(2N, T)` per-step. Fields:
`time_ps (T,)`, `positions_{x,y,z}`, `velocities_{x,y,z}`, `mass_kg (2N,)`,
`droplet_radii (2N,)`, `r0 (N,)`, `E_kin_eV`, `E_pot_eV`,
`E_initial_eV (N,)`, `E_dissip_eV`, `L_droplet_eV_ps`.

### IonCheckpoint — schema version 5 (checkpoint.py:140, version at line 79)

History from v2 → v3 → v4 → v5 documented in `checkpoint.py:54-77`:

- v3: added `droplet_radii_angstrom`, `mass_history_kg`, `E_dissip_eV`.
- v4: added `E_mass_attach_defect_eV`.
- v5: added per-step `temperature_diagnostic (T, 3)`.

Fields (in addition to the neutral-equivalent set):
`positions_final_{x,y,z}`, `velocities_final_{x,y,z}`, `mass_final_kg`,
`mass_history_kg`, `b_ion_outside (N,) bool`, `relative_loss_per_ps`,
`number_of_collisions (int)`, `temperature_diagnostic (T, 3)`.

The loader rejects mismatched `schema_version`
(`_load_checkpoint`, lines 329-333). Bump the version whenever a field is
removed or its meaning changes; additions are backward-incompatible under
the current loader (it enforces `expected_fields ⊆ npz.files` at line 337).

**Drag port note.** If the drag model needs a new per-atom or per-step
diagnostic, bump `_ION_SCHEMA_VERSION` to 6 and add the field to
`IonCheckpoint`. Do **not** silently overload existing fields.

---

## 13. Validation surface (what the drag model must keep green)

### Trajectory comparison
- `i2_helium_md/postprocess/compare_trajectories.py` — `compare_distance`,
  `compare_velocity_magnitude`. References:
  `data/reference/9A_All_Data.csv`, `data/reference/18A_All_Data.csv`.
  Header convention `Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance`
  (CLAUDE.md → "Data Contracts").

### Final-velocity histograms
- `i2_helium_md/postprocess/velocity_distribution.py`. References:
  `data/reference/vmi_summary/vmi_iplus_he.csv`,
  `data/reference/vmi_summary/vmi_iplus_gas.csv`.
- Bin convention (CLAUDE.md → "Known Plotting Conventions"):
  internal bins at 0.04 Å/ps (= 4 m/s), 15-bin moving mean,
  display range up to 2800 m/s.

### Energy balance
- `scripts/post_processing/plot_neutral_energy_balance.py`,
  `plot_ion_energy_balance.py`,
  `plot_ion_temperature_diagnostic.py`.

### Run summary
- `scripts/post_processing/plot_run_summary.py` — driven by USER SETTINGS
  block at top of file (CLAUDE.md → "Post-Processing Workflow").

### Pytest highlights
- `tests/physics/test_potentials.py`
- `tests/physics/test_interactions.py`
- `tests/physics/test_collisions.py`
- `tests/physics/test_energy_balance.py`
- `tests/simulation/test_propagation_step.py`
- `tests/simulation/test_ion_propagation_step.py`
- `tests/postprocess/test_compare_neutral_to_hedft.py`

Canonical full-suite command (CLAUDE.md → "Testing"):

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q
```

---

## 14. Implications for the drag-model port

(Advisory — design decisions still belong to the user; this only catalogues
the unchanged surface around the swap point.)

- The integrator (`leapfrog.py`) is model-agnostic. A drag force can be
  added by extending `_ion_accel_fn` (`leapfrog.py:248`) with a third
  acceleration term, or by composing a new acceleration callable and
  feeding it to `velocity_verlet_step` directly.
- The call sites in the ion driver are:
  - cross-section calculation: `ion_propagation_step.py:211-218` (drop),
  - event sampling: `:220-232` (drop),
  - velocity update: `:237-247` (drop or replace with drag-tick),
  - mass attachment: `:249-256` (independent — decide whether sticky
    collisions still apply in the drag picture).
- If the drag step replaces the stochastic collision, the `E_dissip_eV`
  accumulator can keep being used as the catch-all dissipation channel; the
  per-step `temperature_diagnostic` becomes either NaN (no collision events
  per step) or repurposed.
- `SimConfig` should **grow** new fields for the drag model. Keep the
  hard-sphere fields in place until removal is explicitly approved
  (CLAUDE.md → "Forbidden Without Explicit User Approval").
- Checkpoint schema must bump if new diagnostics are added (§12). Do not
  silently overload existing fields.
- The neutral driver, neutral checkpoint schema, presets except for new
  drag fields, RNG draw order, physical-constants table, and
  post-processing surface remain **off-limits**.

---

## 15. Sources of truth

Python sources (all paths relative to repo root):

- `i2_helium_md/CLAUDE.md`
- `i2_helium_md/physics/constants.py`
- `i2_helium_md/physics/potentials.py`
- `i2_helium_md/physics/interactions.py`
- `i2_helium_md/physics/collisions.py`
- `i2_helium_md/physics/leapfrog.py`
- `i2_helium_md/simulation/propagation_step.py`
- `i2_helium_md/simulation/ion_propagation_step.py`
- `i2_helium_md/simulation/checkpoint.py`
- `i2_helium_md/config.py`
- `i2_helium_md/presets.py`

Documentation cross-references (the existing per-module docs):

- `docs/physics/constants_module.md`
- `docs/physics/potentials_visualization.png`
- `docs/physics/How the potentials work.md`
- `docs/physics/physics_background.md`
- `docs/physics/interactions_module.md`
- `docs/physics/collisions_module.md`
- `docs/physics/leapfrog_module.md`
- `docs/sampling/*.md`
- `docs/post_process/post_processing_strategy.md`

Reference data inventory:

- `data/reference/9A_All_Data.csv`, `18A_All_Data.csv` — HeDFT trajectories.
- `data/reference/vmi_summary/vmi_iplus_he.csv`, `vmi_iplus_gas.csv` —
  experimental VMI references.
- `data/reference/paper_v2/`, `paper_v3/`, `paper_v4/`, `paper_cov/` —
  frozen paper-figure exports.
