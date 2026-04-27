# Migration Log

Chronological record of decisions, deviations, and open questions in the
MATLAB → Python port. Keep this current whenever a design choice is made.

---

## Scope decisions

- **In scope:** single-pulse simulation, neutral + ion propagation, He-DFT
  comparison at R₀ = 9 Å.
- **Out of scope:** pump-probe, effusive dynamics, experimental VMI
  comparison, 18 Å HeDFT comparison (no data in repo), Abel inversion, all
  image-processing utilities.

The effusive-dynamics path existed in the original code as a validation
control (I₂ in vacuum without He droplet) to benchmark the MD against
gas-phase measurements. Since the goal here is droplet physics vs HeDFT,
we drop it. It could be re-added later as a single config flag that bypasses
the droplet potential and collisions.

---

## Architectural decisions

### 1. Replace MATLAB globals with a `SimConfig` dataclass

~36 MATLAB `global` variables were spread across `run_simulation.m`,
`physical_constants.m`, and the various `inputfiles_*/*.m` scripts. We
consolidate them into a single strongly-typed dataclass (`config.py`).
Every function that needs a parameter takes `cfg: SimConfig` explicitly.

**Rationale:**
- Single source of truth for all parameters
- IDE autocomplete + type checking catch typos
- Two simulations with different parameters can coexist in memory (impossible
  with globals)
- Derived quantities (`E_min_eV`, `num_timesteps_neutral`, etc.) stay
  consistent with their inputs automatically via `@property`

### 2. Presets as functions, not scripts

MATLAB input files (`single_pulse_N2000.m`) were scripts that mutated
globals. We replace each with a **function** (`single_pulse_N2000()`) that
returns a fully-configured `SimConfig`. Users override fields via keyword
arguments:

    cfg = single_pulse_N2000(num_molecules=500, seed=42)

**Rationale:**
- Explicit output (no hidden global state)
- Trivial parameter scans via loops
- Presets can build on other presets

### 3. Unit suffixes in field names

`T_particles_K`, `binding_energy_I_atom_eV`, `v_limit_m_per_s`, etc.
The MATLAB code relied on comments to indicate units, which were often
inconsistent or missing.

**Rationale:** makes unit errors visible at the call site, not buried in
definitions.

### 4. Analytical derivatives where MATLAB used finite differences

`droplet_potential.m` used `(U(r+h) - U(r))/h` with `h=1e-6` for the force.
We provide `droplet_force()` with the exact analytical derivative.

**Rationale:** no numerical noise at small h, and no ambiguity about
forward/backward/central differences. Unit-tested against finite difference.

### 5. Data files moved to `data/reference/`

All hardcoded Windows paths (`T:\github synchronized\...`) are removed.
Reference data files are copied once into `data/reference/` and referenced
via `cfg.data_dir`.

**Files to extract from legacy repo:**

| Legacy path | → | New path |
|---|---|---|
| `HeDFT_MD_comparison_neutral/custom_start_interpolating_functions.mat` | → | `data/reference/hedft_custom_start.mat` |
| `single_pulse_simulation/HeDFT_comparison/9Angström/data_vabs2.csv` | → | `data/reference/hedft_9A_velocity.csv` |
| `single_pulse_simulation/HeDFT_comparison/9Angström/R1-R2.csv` | → | `data/reference/hedft_9A_distance.csv` |

---

## Physics decisions

### Missing `g()` function (resolved)

`get_morse_potential_X.m` calls `g(9, 0.3, r)` but the function isn't
defined in the legacy repo. The user located its definition via MATLAB's
`type g` command:

    function res = g(mu, sig, x)
        res = exp(-(x - mu).^2 / (2*sig^2));
    end

Unit-height Gaussian with **argument order `(mu, sig, x)`**. Ported as
`_gaussian(mu, sig, x)` as a private helper inside `physics/potentials.py`,
since it's only used there.

Initial port mixed the argument order — user caught this before any wrong
physics was committed.

### Morse I₂⁺ state-select ported even though it's usually inactive

For the standard I⁺–I⁺ interaction, the code uses pure Coulomb `14.4/r`, not
the Morse I₂⁺ curves. The Morse state-select code is only triggered when
`single_charge_ionization_allowed = True` and exactly one atom of the pair
is ionized.

We port and test it anyway because:
- It's needed for the I–I⁺ mixed-charge mode
- It's small, self-contained, and already tested
- Removing it would be premature optimization

See `docs/physics_background.md` §2 for the detailed discussion.

### 1-based vs 0-based state indices

MATLAB uses `state_select_id ∈ {1, 2, 3, 4}`. Python uses
`state_ids ∈ {0, 1, 2, 3}`. Documented in the docstring of
`morse_I2plus_state_select`.

---

## Unresolved items

### External toolbox `T:\github synchronized\VMI_matlab`

Referenced by `setup_VMI_path_office_flir.m`. After inspection this is
purely lab-data analysis (Abel inversion, camera calibration) and is NOT
needed for the simulation itself. Confirmed out of scope; no port required.

### 18 Å HeDFT reference data

Referenced via `T:\Cloud\MATLAB iodine\MartiPiNotes\...` — not present in
the repo. We support only 9 Å comparison. Adding 18 Å later is a matter
of dropping the files into `data/reference/` and extending the loader.

---

## Testing strategy

- Every physics module has a `tests/test_<module>.py` covering:
  - Direct formula ports (compare to MATLAB formula on a grid)
  - Edge cases (r = 0, r = ∞, empty arrays)
  - Shape and error handling
  - Regression values (hand-computed baselines)
- `tests/smoke_test_*.py` exists for the sandbox environment where pytest
  isn't available. Can be deleted once the project is only run on a machine
  that has pytest (i.e. the user's own PC).

### Step 6 — leapfrog integrator

The MATLAB `frog_step_neutral.m` / `frog_step_ion.m` files are monolithic
functions that interleave the velocity-Verlet algorithm with force
computation (droplet + partner interactions + optional scattering +
optional droplet charges). We split this into three layers:

- `velocity_verlet_step(pos, vel, acc_fn, dt)` — pure algorithm, takes a
  callable that computes acceleration. Knows no physics.
- `_neutral_accel_fn` and `_ion_accel_fn` — assemble total acceleration
  from droplet + partner contributions.
- `make_neutral_step(cfg, mass, droplet_radii)` and `make_ion_step(...)` —
  factory functions that bind the context once and return a single-arg
  step function for the simulation loop.

**Algorithm note:** the MATLAB name "frog_step" is misleading. The
algorithm is velocity-Verlet in kick-drift-kick form, not classical
leapfrog. Velocity-Verlet stores position and velocity at the same
timestep, which is easier to checkpoint. Mathematically equivalent to
leapfrog for conservative forces.

**Features skipped:**
- `he_direction_scattering` (velocity randomization) — was zero in input
  files, never activated.
- `additional_droplet_charges` (charged droplet) — pump-probe scenario,
  out of scope.

**Energy conservation verified:** a vibrating I₂ molecule at 1 fs timestep
conserves energy to within 0.18% over 200 steps — well within the 1%
test tolerance. Drift is dominated by finite-difference force truncation,
not by the integrator itself.

### Test-tolerance lessons learned

- A discretization mistake in `test_vectorized` of `_gaussian` looked like a
  physics failure: `np.linspace(7, 11, 100)` gives step 0.0404 Å, so the
  nearest grid point to `mu=9.0` is at 8.98, and a σ=0.3 Gaussian evaluates
  there to 0.99774 — not 1.0 within 1e-3 tolerance.

  **Resolution:** use a loose-but-meaningful test (`0.99 < max ≤ 1.0` and
  location check). This catches amplitude bugs and center bugs without
  depending on grid alignment.

- Testing `partner_interaction_neutral` at R = R_e expected force ≈ 0, but
  the finite-difference form `F = (U(R) - U(R+h))/h` with h=1e-4 is only
  accurate to first order. Near a minimum where U'(R) = 0, it returns
  F ≈ -h/2 · U''(R_e) — small but not zero. For I₂ X state this is
  ~0.04 Å/ps², not 1e-2.

  **Resolution:** tolerance should be informed by the FD truncation error.
  Compare against the "real" forces at nearby points rather than against
  zero. Loosened tolerance to 1 Å/ps² which is still 100× tighter than the
  nearby true forces of 60–100 Å/ps².

  **Lesson:** wherever we use FD, expect truncation artifacts
  proportional to h and the second derivative. A future improvement is to
  replace FD with analytical derivatives where the potentials allow.

---

## Migration progress

See `README.md` for the live progress table.
