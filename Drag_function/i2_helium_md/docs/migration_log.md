# Migration Log

Chronological record of decisions, deviations, and open questions in the
MATLAB → Python port. Keep this current whenever a design choice is made.

---

## Project quality principles

These guide every porting decision and code review:

1. **No duplicate implementations of the same physics.** Whenever the same
   conversion, formula, or constant appears in two places, we consolidate it
   into a single source of truth in `constants.py` (or the appropriate
   shared module). Parallel implementations *will* drift over time —
   even if the numbers agree today.

2. **No dead code.** Unused imports, commented-out blocks, and "in case
   we need it later" branches all get removed. Git history preserves
   the past; the current source should reflect only the current intent.

3. **Units and conventions are encoded in names.** `mass_kg`, not just
   `mass`. `T_particles_K`, not just `T_particles`. `R0_GS_angstrom`, not
   `R0_GS`. This makes unit mistakes visible at every call site.

4. **Validation early, errors loud.** Wrong shape → `ValueError`. Wrong
   collision mode → `ValueError`. Wrong type → caught by the type system
   or asserted. We never let "silent broadcasting" or implicit conversions
   produce subtly wrong physics.

5. **Every public function has a docstring with units.** Inputs, outputs,
   shapes, units, edge cases. If reading the function requires reading a
   different file's comments, it's not documented enough.

6. **Modules are organised by concern, not by file size.** `physics/` is
   pure science. `sampling/` is randomness. `simulation/` is orchestration.
   Each layer depends only on the layers below it.

7. **Tests document intended behaviour.** A test that passes by accident is
   worse than no test. We choose tolerances based on the truncation error
   of the underlying numerical method, and we comment when a test is
   guarding against a specific MATLAB-vs-Python deviation.

8. **Audits are mandatory, not optional.** Periodically — and especially
   after any refactor — we re-read modules with fresh eyes looking for
   dead imports, drift between parallel implementations, and
   inconsistencies between code and docs. Two such audits already caught
   real issues (the `U` import in `leapfrog.py`; the misleading
   "all forces analytical" claim in the migration log).

Departures from these principles are documented as "open items" with a
plan to resolve.

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

### 4. Analytical derivatives — for the droplet force only (so far)

The MATLAB code uses finite differences `(U(r+h) - U(r))/h` for **all** forces:

- Droplet solvation force: `h=1e-6` in the neutral/ion propagation scripts
- Partner I-I Morse force: `h=1e-4` in `add_partner_interaction.m`
- Partner I⁺-I⁺ Coulomb force: `h=1e-4` in `add_partner_interaction_ion.m`

In Python we currently have **a mixed approach**:

| Force | Form | Where |
|---|---|---|
| Droplet solvation | ✅ analytical: `(E / (s·√π)) · exp(-((r-c)/s)²)` | `potentials.droplet_force` |
| Partner Morse  | ❌ finite-difference (h=1e-4) | `interactions._force_from_potential_fd` |
| Partner Coulomb | ❌ finite-difference (h=1e-4) | `interactions._force_from_potential_fd` |

We started with the droplet derivative because it has a clean closed form
and was easy to verify against FD. The Morse and Coulomb partner derivatives
were left as finite differences for two reasons:

1. To preserve byte-compatibility with MATLAB during the early port stages.
2. The Coulomb-plus-conditional-Morse mix in `ion_interaction_potential` is
   awkward to differentiate symbolically when single-charge ionization mode
   is enabled (per-molecule state-select).

**Next refactor on the queue:** make the partner Morse force analytical.
The derivative is straightforward:

    d/dr [ D_e * (1 - exp(-a(r-R_e)))² ]
        = 2·a·D_e·(1 - exp(-a(r-R_e)))·exp(-a(r-R_e))

And for the Xdip Gaussian dip:

    d/dr [ -A · exp(-(r-mu)²/(2σ²)) ]
        = A · (r-mu) / σ² · exp(-(r-mu)²/(2σ²))

Doing this would eliminate the FD truncation noise we currently see at
potential minima (~0.04 Å/ps² residual force at R_e).

The pure-Coulomb derivative is also trivial (`-q1·q2·14.4/r²`), so we should
do that at the same time.

The Morse-I2+ state-select case is the only awkward one, and is rarely
active in practice (only when `single_charge_ionization_allowed=True`).

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

### Unit-conversion refactor (post-Step 7 cleanup)

User audit caught two issues in `leapfrog.py`:

1. `U` (atomic mass unit) was imported but never used.
2. The same physical conversion ("force in eV/Å on a mass in kg → acceleration
   in Å/ps²") was implemented twice with different numerical forms:

   - **interactions.py**: `F / (mass / U) * 9648.533`
     (i.e. mass in u, factor 9648.5 with units (Å/ps²)·u/(eV/Å))
   - **leapfrog.py**: `F * 1.602e-9 / mass * 1e-14`
     (mass in kg, two-stage conversion via Newtons and m/s²)

   Both expressions give the same numerical result to 1 part in 10⁴ (the
   difference comes from rounding in 1.602e-9 vs the more precise constant).
   Functionally equivalent, but a maintainability hazard.

**Fix:** introduced a single shared constant
`EV_PER_ANGSTROM_PER_KG_TO_A_PER_PS2 = EV * 1e-4 = 1.602e-23` in
`physics/constants.py`. Both `interactions.py` and `leapfrog.py` now import
and use it. The kg-based form is more natural since `cfg`-derived masses are
already in kg.

**Verification:** all 80+ smoke checks still pass. The energy-conservation
test still shows 0.18% drift over 200 steps. Numerical results from the
neutral and ion step tests changed at the 4th significant digit (e.g.
21.308 → 21.306 Å/ps), exactly the size of the previous duplicate-constant
discrepancy of 1 part in 10⁴. The new version is internally consistent
because both force paths use the same constant.

**Lesson:** when porting cross-cutting physics constants, consolidate them
in `constants.py` from the start. Pulling them into "wherever they're used"
seems convenient but creates duplicates that drift over time.

### Step 7 — droplet-size sampling

Ports `generate_droplet_sizes.m` (the full pickup simulation) and
`get_dropletsize.m` (the empirical mean-size correlation). The MATLAB code
had two near-duplicate function files (`generate_droplet_sizes.m` and
`generate_droplet_sizes_simpler.m`) that differed only in whether the
pickup-cell simulation ran. We unified them into a single
`sample_droplet_sizes(cfg, mode=...)` function with `mode="raw"` or
`mode="post_pickup"`.

**Key decisions:**
- All `bayes_hist` plotting code stripped — visualisation belongs in
  postprocessing, not samplers.
- Post-pickup mode raises `ValueError` if too few one-pickup droplets are
  found, instead of MATLAB's confusing index error. Caller can either
  increase oversampling or fall back to `mode="raw"`.
- RNG is a per-call `np.random.Generator` derived from `cfg.seed`,
  replacing MATLAB's global `lognrnd` state.
- `_evaporation_per_pickup` returns negative `dN` (so `samples + dN`
  updates droplet sizes naturally), matching MATLAB's sign convention but
  documented explicitly.

**Note on use:** the default `single_pulse_N2000` preset sets
`use_single_droplet_size = True` and pins every molecule to N=2000, so
the sampler isn't actually invoked for our target scope. We port and test
it anyway because (a) the config exposes the realistic mode and (b) it's
a clean entry point into the sampling layer.

### Step 8 audit: thesis-text vs thesis-figure E_solv inconsistency

When verifying our pickup simulation against the supervisor's thesis
figure 3.2 (post-pickup droplet size distributions), an audit revealed
**an inconsistency within the thesis itself**:

- **Thesis text** explicitly states `E_solv = 30 meV` for iodine.
- **Thesis figure 3.2** can only be reproduced with `E_solv = 14 meV`.
  At 30 meV, the reduced-σ correction is so weak that reduced and
  normal distributions overlap; at 14 meV they separate as the figure
  shows.
- **Legacy MATLAB code** uses `E_solv = 14 meV` (matches the figure).

The first instinct was "trust the text, it must be the right value."
That gave a default of 30 meV, but a side-by-side comparison (which the
user pushed for) showed the resulting plot looks nothing like the
thesis figure. **The figure is the ground truth, not the text** -- the
figure was actually produced with the same code we have, while the
text was either revised later, copied from a different paper, or
referred to a different physical context.

**Resolution:**
- `E_solv_meV` is a kwarg on `sample_droplet_sizes()` and helpers.
- **Default = 14 meV** (matches both the figure and the legacy MATLAB).
  Pass `E_solv_meV=30.0` explicitly if you want the thesis-text value.
- Two regression tests lock in the figure-matching behaviour:
  - peak position at T=18 K should be ~2500-3500 atoms
  - reduced-σ distribution should shift right of normal-σ by >30%
- Inconsistency is documented prominently in the source comments.

**Quantitative verification (E_solv = 14, p = 25 mbar, d = 5 µm):**

| Feature | Thesis fig | Our reproduction |
|---|---|---|
| T=18 K normal peak | ~2500 | ~2500 ✓ |
| T=15 K normal peak | ~5000 | ~5000 ✓ |
| Reduced σ shifts right | yes, by ~5000 | yes, by ~5000-7000 ✓ |
| Reduced σ has sharp left cutoff | yes | yes ✓ |
| Means in panel (a) | normal/reduced split | reproduced ✓ |

**Lessons:**
1. Legacy code can disagree with the published paper that documents it.
2. **Within** a paper, the text and figures can disagree. Visual
   regression against figures beats following text quotations.
3. When the user pushes back on a "fix", verify with a side-by-side
   comparison rather than defending the change.

### Step 8 — radial position sampling

Ports `generate_radial_samples_3d.m`. Each molecule is placed at a radial
distance from its droplet center drawn from the Boltzmann density
`p(r) ∝ r² · exp(-U_drop(r-R) / k_B T)`.

**Key decisions:**

- **Cleaner energy unit handling.** MATLAB went meV → J → eV → meV in three
  separate places via the `[0:0.001:E_max]/1000*eV` idiom. We use eV
  throughout, with one explicit `K_B / EV` conversion to get `k_B` in eV/K.
- **Dropped the unused `p(E)` normalization.** MATLAB normalized the
  Boltzmann density `p(E)` once over an energy axis, then renormalized
  `p_radius(r)` over the radial axis. The first normalization had no effect
  on the rejection sampling and was confusing; we removed it.
- **Smarter rejection sampling.** Instead of MATLAB's hardcoded 1000-
  proposals-per-batch loop, we estimate the acceptance rate from a
  grid-sampled density and propose the right batch size in one go (with a
  5% safety factor and a top-up loop in case we under-shoot).
- **`np.unique` for heterogeneous droplet ensembles.** Each unique droplet
  radius gets its own rejection batch; results are scattered back to
  per-molecule slots via `return_inverse`. Faster than re-running the
  sampler for every duplicate radius.
- **Unphysical configuration is loud, not silent.** Density identically
  zero (e.g. T = 0 with infinite well) raises `RuntimeError` with a
  helpful message instead of dividing by zero somewhere downstream.

**Visual verification:** `radial_positions_visualization.png` shows
sampled histograms at T = 0.4, 4, 50 K tracking the analytical
density to within Monte Carlo noise.

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
