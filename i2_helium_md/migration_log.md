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

9. **Literal transliteration before clean refactor for any reference reproduction.**
   When porting code that produces a known reference output (a thesis
   figure, an experimental dataset, a published table), write the most
   literal possible line-by-line transliteration FIRST and verify against
   the target. Only then refactor for clarity. Going straight to a
   "clean" version skips the verification step that would catch subtle
   convention mismatches (density factors, normalisation choices,
   off-by-one indexing) immediately. The transliteration can be deleted
   once the clean version produces an identical numerical result. This
   principle was added after the thesis figure 3.2 reproduction took
   four wrong attempts because we kept refactoring before verifying.

10. **Don't keep bad features just to stay byte-identical to legacy.**
    Where the legacy MATLAB code uses 4-significant-figure constants
    or other approximations that we know to be inaccurate, the Python
    port uses the correct modern values (CODATA 2022 / SI 2019) and
    documents the choice prominently. "Reproduces the MATLAB output
    bit-for-bit" is not a goal in itself — it's a useful regression
    target *during a port*, but where the legacy was wrong, the new
    code should be right. Differences of order ~100 ppm in energy
    calculations are explicitly tolerated and noted. This principle
    was added after I downgraded a precise CODATA constant to match
    the legacy 4-sig-fig value, and the user correctly objected.

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

MATLAB input files (`single_pulse_N2000.m`,
`single_pulse_droplet_distribution.m`) were scripts that mutated globals. We
replace each with a **function** (`single_pulse_N2000()`,
`single_pulse_droplet_distribution()`) that returns a fully-configured
`SimConfig`. Users override fields via keyword arguments:

    cfg = single_pulse_N2000(num_molecules=500, seed=42)
    cfg = single_pulse_droplet_distribution(num_molecules=500, seed=42)

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

### Step 10 — neutral propagation (in progress)

The MATLAB `vmi_sim_3d_neutral_propa_HeDFT_mimic.m` is ~1000 lines of
mixed concerns: angle sampling, initial-state assembly, leapfrog
propagation, hard-sphere collisions, energy bookkeeping, debug
plotting. We deliberately **skip the literal-transliteration step**
for this file (per user direction) -- the MATLAB itself is the
reference. Instead we build the leaves first, in isolation, and
glue them together via a driver only after each leaf is fully
verified.

**Plan:**

1. ✅ `sampling/orientations.py` -- pure angle/bond-length sampling.
2. ⏳ `physics/collisions.py` -- hard-sphere collision physics
   (Braun thesis D.2). Mode 3 only for now.
3. ⏳ `simulation/neutral.py` -- driver that assembles initial state
   and runs the propagation loop.

**Step 10a: `sampling/orientations.py`**

Direct-but-clean port of the orientation-sampling block (MATLAB lines
~225-275). The module exports `sample_orientations(...)` returning a
frozen `MolecularOrientations` dataclass with five named arrays.

Key design choices:

- **Returns a dataclass, not five loose arrays.** Makes the contract
  explicit and guarantees consistent length.
- **Position angles always uniform on the sphere.** Only the molecular
  axis orientation is mode-dependent (anisotropic vs isotropic).
- **Conversion to atomic xyz is OUT of scope.** That depends on the 2N
  array layout convention owned by `simulation/neutral.py`.
- **Anisotropic sampler uses bounded batched rejection.** Cos² rejection
  has ~33% acceptance, so we propose ~3.5N candidates per batch and
  accumulate; cap at 20 batches as a safety guard.
- **Statistical regression tests:** `<(cos α sin δ)²> ≈ 0.6` in
  anisotropic mode and `≈ 1/3` in isotropic mode. These distinguish a
  correct cos²-weighted distribution from a bug-introduced uniform one.

**Step 10b: `physics/collisions.py`**

Mode-3-only hard-sphere collision physics, based on Andreas Braun's
PhD thesis section D.2. Two pure functions, no SimConfig dependency:

- `sample_collision_events(...)` -- given per-particle distance,
  depth, and energy, decide who collides via Poisson thinning
  with `p = d * sigma * rho_droplet`.
- `apply_collision(...)` -- sample impact parameter, compute new
  energy and lab-frame angle, build new 3D velocity vector with a
  random orthonormal basis perpendicular to the incoming velocity.

Key design choices:

- **No SimConfig coupling.** Driver pulls cfg fields and passes them
  in. Lets us test physics in isolation and reuse for both neutral
  and ion stages.
- **Modes 1 and 2 deliberately not ported.** They're legacy code
  paths; we'll add them if needed.
- **Azimuthal smearing convention preserved exactly from MATLAB.**
  The `COSBETA = uniform(-1, 1)` sampling in MATLAB is NOT a uniform
  azimuth, but the random reference direction (re-drawn per particle
  per step) makes the resulting scattered velocities isotropic in
  the perpendicular plane. We mirror this exactly and verify
  isotropy with a regression test on `<cos 2φ>` (catches 2-fold
  bias). Empirical `<cos 2φ> ≈ -0.002` with 200k samples.
  The convention is documented prominently in the module docstring
  with instructions for switching to uniform-phi if desired later.

Statistical regression tests lock in:

- Empirical collision rate matches `d × σ × ρ`.
- Mean fractional energy loss for I/He is ~5–7% per collision.
- Max fractional energy loss for equal masses approaches 1.0
  (full transfer).
- `<cos θ_lab>` matches the analytic 2-body transformation integral.
- Perpendicular-plane azimuth has zero `<cos φ>`, `<sin φ>`,
  and crucially `<cos 2φ>`.
- Speed self-consistency: `|v_new|² = 2E₁/m` to numerical precision.

**Step 10c-i: `simulation/initial_state.py`**

Builds the t=0 physical state and pre-allocates trajectory arrays for
the rest of the run. Single public function
`build_initial_state(cfg, num_steps, rng)` returning a fully-allocated
`NeutralCheckpoint` with column 0 populated.

Internally:

1. Sample droplet sizes (`sample_droplet_sizes` or `single_droplet_size`).
2. Convert N -> radius via `droplet_radius_bulk_angstrom`.
3. Sample radial positions of molecule centres.
4. Sample axis orientations + bond lengths.
5. Compute v0 from laser parameters (with `partner_interaction` branch
   for the energy budget).
6. Assemble per-atom xyz and velocities using the 2N layout.
7. Allocate full `(2N, num_steps)` arrays, fill column 0.
8. Return `NeutralCheckpoint`.

Key design choices:

- **Single function, not split into "sample then assemble".** The
  layout convention (2N) and the conversion from angles to xyz is
  the responsibility of this module. Splitting would have created
  redundant boundaries.
- **No DFT pre-fill.** The optional `custom_DFT_start` block from
  MATLAB is left out; if needed it'll be added as a separate function
  invoked by the driver.
- **Allocates full trajectory arrays up front.** Memory cost is
  computed by the driver; the auto-stride decision happens upstream
  before this function is called.
- **Variable name caveat documented.** The MATLAB `E_initial`
  variable is in joules despite its name. Our Python preserves the
  physics but uses explicit unit-suffixed names (`E_initial_J`).
  The per-molecule eV-valued field on the checkpoint is named
  `E_initial_eV` and stores `hc/lambda` in eV.
- **MATLAB's `sin(δ+π)`, `cos(δ+π)` simplifications.** We write
  `-sin(δ)` and `-cos(δ)` explicitly. The MATLAB form is preserved
  in a comment.

**Tests:** 17 pytest tests + 32 sandbox checks covering shapes, bond
length, atom-2 anti-alignment, single-pulse v=0, reproducibility,
photon energy, single-droplet-size mode, validation errors, save/load
round trip.

**Step 10c-iii: `simulation/neutral.py` driver (refactor + driver)**

Final piece of the neutral propagation. Two changes in this step:

1. **Refactored `propagation_step.py` to a pure function.** The
   original 10c-ii API took a `NeutralCheckpoint` and `t_id`,
   mutating column `t_id+1` in place. That coupling made auto-stride
   impossible without weird tricks. Refactored to:

   * `NeutralStepState` -- frozen dataclass holding the minimum-
     sufficient state (positions, velocities, cumulative diagnostics,
     time).
   * `neutral_propagation_step(state, *, cfg, mass_kg, droplet_radii,
     prev_distance_angstrom, rng) -> NeutralStepState` -- pure
     function. Does not mutate inputs. Doesn't know about checkpoints.
   * Helpers `state_from_checkpoint_column` and
     `write_state_to_checkpoint_column` for the driver to
     bootstrap/serialize.

   This decouples the physics step from the storage stride, which
   was the whole point of asking for a pure function. The same 11
   physics tests from the in-place version were rewritten for the
   pure API and all pass.

2. **Added `simulation/neutral.py` with `run_neutral_propagation`.**
   The driver:
   * Decides internal step count (single_pulse -> 2; otherwise
     `cfg.num_timesteps_neutral`).
   * Decides storage stride. If full-resolution checkpoint would
     exceed `max_bytes` (default 300 MB), only every K-th internal
     step is stored. Internal stepping always at `cfg.dt_neutral`
     so collision rate / leapfrog stability are preserved.
   * Calls `build_initial_state` (sized for stored steps).
   * Raises `NotImplementedError` if `cfg.custom_DFT_start` is True
     (the TD-HeDFT pre-fill is not yet ported).
   * Runs the inner loop using `neutral_propagation_step`, copying
     every K-th state into the checkpoint.
   * Always stores the final state in the last reachable column,
     even if `(num_internal - 1) % stride != 0` -- ensures
     trajectories end at the actual end-time.
   * Optionally saves via `RunDirectory`.

**Tests:** 15 pytest + 23 sandbox checks. End-to-end smoke includes
single-pulse, long run, strided run (verified storage size stays
under tight budget), DFT-stub, RunDirectory round trip, and energy
bookkeeping (~3% drift over 50 steps in collision-active run, <10%
tolerance).

**Memory budget verified empirically:**

| Run | Internal | Stride | Stored | Size |
|---|---|---|---|---|
| Single-pulse N=2000 | 2 | 1 | 2 | 1 MB |
| HeDFT 9 Å N=500 | 2000 | 1 | 2000 | 160 MB |
| HeDFT 9 Å N=2000 | 2000 | 3 | 667 | 214 MB |
| Long N=2000 | 20000 | 22 | 910 | 291 MB |

Step 10 of the migration plan is now complete.

**Bug fix follow-up (E_pot at t=0 missing partner Morse):**

The user reported a test failure: `test_total_energy_drift_small`
expected drift < 1% but got 1149%. Investigation revealed that
`build_initial_state` was computing `E_pot[t=0]` as the droplet
solvation potential ONLY, while `propagation_step.py` (and MATLAB
from t=1 onwards) computed it as droplet + half partner Morse. The
discontinuity between t=0 and t=1 was the full Morse pair energy,
~3 eV per molecule for R0=9 A initial conditions.

The legacy MATLAB code has the same bug: line 476 of
`vmi_sim_3d_neutral_propa_HeDFT_mimic.m` has
`E_pot(:,1) = droplet_potential_atom(...)` -- no partner term --
while line 885 has the full
`E_pot(:,t_id+1) = droplet_potential_atom(...) + [E_pot_partner;E_pot_partner]/2`.

Per project principle #10 ("don't preserve legacy approximations"),
we treat this as a bug and include the partner term at t=0.

After the fix, energy conservation in the test scenario is exact
(drift = 0 to machine precision, since single_pulse atoms don't
move). The long-run energy bookkeeping with collisions is unchanged
(2.83% drift over 50 steps), confirming the fix doesn't disturb the
collision physics.

Added two regression guards:

* `TestEPotIncludesPartner.test_E_pot_t0_includes_partner_morse`:
  E_pot[t=0] should be CONTINUOUS with E_pot[t=1] from a single-pulse
  step (jump < 1e-3 eV).
* `TestEPotIncludesPartner.test_E_pot_t0_at_R0_9_is_a_few_eV`:
  for R0=9 A and N=5, E_pot[t=0] > 1 eV (dominated by the Morse term).

These would have caught the bug immediately if I'd written them when
building `build_initial_state`. Lesson: when one module computes a
quantity using a particular formula and another module is supposed
to be consistent with it, **write a regression test that explicitly
checks continuity at the boundary**. Don't trust that "they probably
match because they look similar."

**The (now removed) verbose Step 10c-ii section was a duplicate of
content captured more accurately under Step 10c-iii's "refactor"
bullet, which describes the API as it actually exists (pure function,
not in-place).**

**Step 10c-prep: per-atom checkpoint shapes + droplet-radius utility**

Two prep changes before building the `simulation/neutral.py` driver:

1. **`NeutralCheckpoint` and `IonCheckpoint` energy/L_droplet diagnostics
   moved from `(N, num_steps)` to `(2N, num_steps)`** to match the
   legacy MATLAB convention. Schema version bumped from 1 to 2 (old
   files will fail to load with a clear error). The 2N layout lets us
   inspect per-atom energy bookkeeping for debugging; per-molecule
   values are recovered by summing/averaging atom 1 + atom 2 indices.
   Only `r0` and `E_initial_eV` remain per-molecule (they describe the
   molecule as a whole, not individual atoms).

2. **New utility `droplet_radius_bulk_angstrom(N)` in
   `physics/constants.py`.** The legacy MATLAB neutral propagation
   uses `R = 2.22 * N^(1/3)` -- a 3-significant-figure rounding of
   the formula derived from the bulk helium density. Per our "don't
   preserve legacy approximations" principle, we use the precise
   value (~2.2173 instead of 2.22). The 1200 ppm difference at
   N=2000 is documented in the docstring as expected.

   The function also documents the legacy MATLAB inconsistency: the
   neutral-propagation code uses bulk density (2.22 prefactor), but
   the pickup-cell sampler in `generate_droplet_sizes.m` uses 0.8x
   bulk (2.39 prefactor). Both are mirrored faithfully in the
   respective Python modules.

**Audit follow-up (validator shape coverage):** when the user asked
whether the schema-v2 change affected the `RunDirectory` wrapper, an
audit found that `RunDirectory` itself was shape-agnostic (it just
delegates to `save_neutral_checkpoint`/`load_neutral_checkpoint`).
**However**, the validator function `_validate_against_cfg` only
checked shapes for `mass_kg` and `droplet_radii`; it did NOT check
the new (2N, T) shapes for the energy arrays. A regression where
some future code wrote `E_kin_eV` with the old (N, T) shape would
have round-tripped silently rather than being caught.

Strengthened `_validate_against_cfg` to check:
- Static (2N,) arrays: `mass_kg`, `droplet_radii`, `mass_final_kg`,
  `positions_final_*`, `velocities_final_*`.
- Static (N,) arrays: `r0`, `E_initial_eV`, `b_ion_outside`.
- Trajectory (2N, T) arrays: all positions, velocities, energy
  diagnostics, plus consistency that they all share the same `T`.
- `time_ps` must agree with the trajectory `T`.

Added two pytest regression guards in `test_checkpoint.py` and
extended the `smoke_test_run_directory` to round-trip the actual
output of `build_initial_state` with cfg-validated load. This locks
in the contract: if the schema or the build function ever drifts,
either the file will fail to load or the regression test will fail.

**Audit follow-up (twice):** an early version of `collisions.py`
defined local constants `AMU_KG` and `EV_J` (with the more precise
CODATA value `1.602176634e-19`). I initially "fixed" this by
downgrading the local values to match the legacy MATLAB
4-significant-figure constant `EV = 1.602e-19`, on the rationale that
"matching MATLAB byte-for-byte is good." **The user pushed back: that
reasoning sacrifices physical accuracy for a false reproducibility
goal.** Legacy MATLAB constants are wrong by ~100 ppm; throwing away
~10 ppm of accuracy in our own code to "match" them is keeping a bad
feature.

**Correct fix:** updated `physics/constants.py` to use modern
CODATA 2022 / SI 2019 values throughout:

- `EV = 1.602176634e-19` (exact by 2019 SI redefinition)
- `E_CHARGE = -1.602176634e-19` (exact)
- `U = 1.66053906892e-27` (CODATA 2022)
- `EPSILON_0 = 8.8541878188e-12` (CODATA 2022)
- `K_B = 1.380649e-23` (already exact, no change)
- `HC = 1239.841984` (was 1240, ~130 ppm wrong)
- `EV_PER_WAVENUMBER = 1/8065.543937` (CODATA 2022)

This means our results will not be byte-identical to the legacy
MATLAB output — but they will be **physically more accurate**. The
discrepancy is at most 100 ppm for any individual energy
calculation, which is below the precision of any quantity we
actually compare against experiment.

If a future user needs to reproduce a legacy MATLAB run bit-for-bit,
they can override constants by reassignment before any other module
loads (documented in `docs/constants_module.md`).

**Lesson:** "byte-identical to the legacy code" is a useful
regression target during a port, but it is **not a goal in itself**.
Where the legacy code was sloppy, the new code should be better, and
documented as such.

### Step 9 — checkpoint I/O

Replaces MATLAB's bare `save('neutral_propagation_checkpoint', 'var1', 'var2', ...)`
calls with two clean dataclasses (`NeutralCheckpoint`, `IonCheckpoint`) plus
matching `save_*` / `load_*` helpers, persisting to `.npz` files. A second
module `run_directory.py` adds a `RunDirectory` convenience layer on top
that gives each simulation run a self-describing folder containing
`cfg.json`, `neutral.npz`, `ion.npz`.

**Key design choices:**

- **Two clean dataclasses** -- explicit field list, units in field names
  (`mass_kg`, `time_ps`, `E_kin_eV`), shape documented in the docstring.

- **Schema versioning.** Every checkpoint carries a `schema_version: int`.
  Loader refuses to load incompatible versions. Adding fields is
  backward-compatible; removing or renaming bumps the version.

- **No constants saved.** MATLAB save dumped `eV`, `u`, `mass`, etc. -- all
  recoverable from `constants.py`. Saving them creates drift risk.

- **No config flags saved.** MATLAB saved `effusive_dynamics`, mode flags,
  binding energies, ... -- these belong to `cfg`. The loader takes an
  optional `cfg` argument to validate shape consistency.

- **`.npz` instead of `.mat`.** No `scipy.io` dependency. Native to NumPy.
  `savez_compressed` gives 30-60% size reduction on smooth physical
  trajectories.

- **`allow_pickle=False` on load** -- refuses to load files containing
  pickled Python objects (security hygiene).

- **Path discipline.** Caller specifies the output path; the function
  creates parent directories as needed and adds `.npz` extension if
  missing. No `cd()` calls or hardcoded paths.

**Validation behaviour:**
- Missing file -> `FileNotFoundError`
- Wrong/missing `schema_version` -> `ValueError`
- Missing required fields -> `ValueError` with field list
- `cfg.num_molecules` mismatch -> `ValueError`

This is the first I/O module, and it sets the pattern for any future
disk-backed state in the codebase: explicit dataclass + schema version
+ shape validation.

**The `RunDirectory` convenience layer (`simulation/run_directory.py`):**

A second module adds a thin convenience class on top of the raw save/load
functions. The user picks one path -- the run directory -- and never types
filenames again:

```python
run = RunDirectory("data/runs/test01")
run.save_cfg(cfg)
run.save_neutral(neutral_ckpt)
# later, possibly different process:
run = RunDirectory("data/runs/test01")
cfg = run.load_cfg()
neutral = run.load_neutral()    # auto-validates against cfg.json
```

Why this layer:

- **Eliminates filename invention.** Every script uses `neutral.npz`,
  `ion.npz` -- always. No risk of typos like `Neutral.npz` vs `neutral.npz`.
- **Self-describing runs.** The cfg that produced a run lives next to the
  data as `cfg.json`. The legacy MATLAB code had no equivalent and that was
  a known pain point ("which inputfile produced this checkpoint?").
- **Two-script handoff is just a string.** Neutral-stage and ion-stage
  scripts agree on a directory path, not a set of filenames.
- **Forward-compatible.** Adding `pump.npz`, `probe.npz`, or
  `figures/` later doesn't change existing call sites.
- **No cwd dependency.** The legacy code used `cd()` to switch directories
  and relied on relative paths; we always use explicit absolute or
  relative paths.

Validation behaviours added:
- `cfg.json` with unknown fields -> `ValueError` (catches version skew).
- `cfg.json` is **never** silently overwritten by `save_neutral(cfg=...)`
  if it already exists -- this prevents a re-run with a different cfg from
  silently corrupting the run's record of what produced it.
- Loading a checkpoint with no explicit cfg auto-loads `cfg.json` for
  shape validation.

### Step 8 audit: thesis-figure 3.2 reproduction -- the full story

The thesis figure 3.2 reproduction took multiple rounds and several wrong
attempts. The final understanding:

1. The figure was produced by a separate didactic script in Treiber's
   thesis-writing workflow (file:
   ``conditional_droplet_size_distribution_simplified.m``), NOT by the
   production simulation code in this repo.

2. That didactic script uses a closed-form Poisson-convolution formula

       p_conditional(N) = P(k=1 | N) * p_lognormal(N) / Z

   where ``P(k | N)`` is the Poisson PMF with rate ``a * n_gas * sigma(N)``.
   No evaporation, no Monte-Carlo, no destroyed droplets -- a much simpler
   physical model than our ``_simulate_pickup``.

3. Two **subtle but critical** details in that script that I initially
   missed and only caught after the user explicitly disagreed with my
   first reproduction attempt:

   * **Density convention.** The script's ``R(N) = (3N/(4*pi*n_he))^(1/3)``
     uses the **bulk** helium density ``n_he = 2.18e28``, NOT the
     ``0.8 * n_he`` droplet density used elsewhere in the production
     codebase. The 1.077x scaling on R has a 1.16x effect on the cross
     section and dramatically shifts the threshold cutoff.

   * **Normalisation convention.** Both the normal-sigma and reduced-sigma
     conditional distributions are normalised by the **same**
     ``p_k_normalization`` (the normal-sigma value). The reduced-sigma
     distribution therefore does NOT integrate to 1; its integral is the
     fraction of one-pickup events that come from droplets above the
     kinetic-energy threshold. This is why in the thesis figure the
     dashed peaks at T=18 K are tiny (~0.3e-4) and at T=12 K are
     **taller than the solid peaks** (~1.45e-4 vs ~0.85e-4).

   I had originally re-normalised both branches to integrate to 1
   independently, which produced a plot that looked qualitatively similar
   but was quantitatively wrong on every peak height.

**Resolution:**

- Added :func:`conditional_size_distributions_analytical` to
  ``droplet_sizes.py``: literal port of the MATLAB script, returning
  BOTH p_normal and p_reduced as a tuple, with Treiber's exact density
  and normalisation conventions.
- Added :func:`plot_thesis_figure_3_2` to the diagnostics module:
  reproduces the thesis figure visually using the public function.
  Side-by-side comparison now matches pixel-for-pixel on every peak
  position, height, and width.
- Production sampler ``_simulate_pickup`` and its default
  ``E_solv = 14 meV`` are **unchanged**. The Monte-Carlo simulator is what
  feeds the actual initial conditions for our MD runs.
- Regression tests lock in the **specific** thesis-figure signatures, not
  just qualitative shape:
  - normal-sigma integrates to 1
  - reduced-sigma integrates to <1 (Treiber convention)
  - T=12: reduced peak HIGHER than normal peak
  - T=18: reduced peak much LOWER than normal peak (< 40%)
  - T=18: normal peak position in [1500, 4000]

**Lessons learned (this audit alone, in order):**

1. *First attempt:* read the thesis text, defaulted ``E_solv=30`` -- the
   user pushed back with a side-by-side comparison showing the dashed
   lines weren't separated.

2. *Second attempt:* reverted to ``E_solv=14``, declared "match" based on
   eyeballing solid lines only, ignored the dashed lines.

3. *Third attempt:* the user uploaded the actual Treiber script. I claimed
   to port it but missed the density and normalisation details. Generated
   another "match" that the user correctly identified as wrong.

4. *Fourth attempt:* after the user said "this doesn't match at all", I
   wrote a literal line-by-line MATLAB transliteration as a sanity check.
   The transliteration matched the thesis. Comparing it to my "clean"
   port revealed the two missed details. **Final fix gives pixel-perfect
   reproduction.**

**The lesson:** before claiming to reproduce a reference figure, write
the most literal possible transliteration first and verify against the
target. Only then refactor for clarity. Going straight to a "clean"
version skipped the verification step that would have caught the bugs
immediately. This is now reflected in the project quality principles
above.

### Step 8 audit follow-up: diagnostic plot count vs density

A second issue with the diagnostic plot was caught by the user comparing
it side-by-side with the MATLAB ``histogram()`` output: my Panel 2 (the
"raw vs one-pickup overlay") used ``density=True``, which renormalises
each histogram independently to integrate to 1. That made the one-pickup
histogram look like a peak shifted to N~16500, because that's where the
*conditional density* peaks.

MATLAB's ``histogram()`` defaults to **raw counts**, which is what the
reference image (``generate_droplet_sizes.m`` lines 226-229) shows. Raw
counts make the one-pickup histogram appear visibly smaller (because
most droplets don't end up with exactly one pickup) and at a peak only
slightly shifted right of the initial distribution.

**Fix:** Panel 2 now uses ``density=False`` and labels the y-axis
"count", matching MATLAB exactly. Side-by-side comparison with the
MATLAB plot confirms identical shape and ratio.

**Lesson:** when reproducing a reference plot, match the y-axis
convention (counts vs density vs probability mass) -- it's not just
visual styling, it changes which feature is at the peak. Always
side-by-side compare before declaring "looks right".

### Step 8 audit follow-up: diagnostics module

The legacy MATLAB `generate_droplet_sizes.m` has a `debug_pickup_plot=true`
flag that interleaves histograms and per-round `fprintf` reports into the
production sampler. To preserve that capability without polluting the
clean Python sampler, we added a new file
`i2_helium_md/sampling/droplet_sizes_diagnostics.py` with one public
function `diagnose_pickup(cfg, ...)` that:

- Re-runs `_simulate_pickup` with `return_diagnostics=True` (a new flag
  added to the production function -- when False, behaviour is byte-
  identical to before).
- Returns a `(diagnostics_dict, matplotlib_figure)` tuple.
- Prints a per-round physics report to stdout (suppressible).

The figure is a 4-panel layout including the **raw vs one-pickup overlay**
that MATLAB lines 226-229 produced, plus per-round yield bars and physics
diagnostics not in the MATLAB.

**Design choices:**
- Diagnostic mode skips the "too few singles" `ValueError` so users can
  diagnose *why* yield is low, instead of being blocked by the failure.
- `matplotlib.pyplot` is imported inside the figure-building function,
  not at module load. This keeps the module importable in headless
  environments.
- Production sampler's count check is preserved (still raises in
  `mode='post_pickup'`); only the diagnostics path is permissive.

**Use case verified:** the thesis figure 3.2 reproduction
(`thesis_comparison_corrected.png`) was generated using `diagnose_pickup`
with `oversample_factor=20`, allowing both normal and reduced-σ
distributions at all three temperatures (T=12, 15, 18 K) without
hitting the production count-check failure.

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

## Step 11a — Velocity-dependent cross section helper

Both production input scripts (`single_pulse_N2000.m` and
`single_pulse_droplet_distribution.m`) set `sigma_dependent_on_v = true`,
so the ion stage needs `sigma = sigma_0 * v ** sigma_ion_exponent` with
`sigma_ion_exponent = -2`.

Rather than entangling this in the collision sampler itself, I added a
small standalone helper `velocity_dependent_cross_section()` to
`physics/collisions.py`. The driver computes per-particle σ from
current speeds, then passes the array to the existing
`sample_collision_events()` (which already accepted a per-particle
σ array — that part of the API was forward-compatible from Step 10).

**Edge case at `v = 0`:** with negative exponent this yields `+inf`, not
NaN. After consideration we deliberately do NOT clamp, because:

1. `inf` propagates cleanly through `p_scatter = dist · σ · ρ_droplet`
   to give `p = +inf`.
2. The downstream check `trial < p_scatter` evaluates True for any
   finite `trial`, i.e. the ion always collides — physically correct
   for an "infinitely-slow" particle.
3. The Landau cutoff (`E0 < E_min`) provides the actual physical
   low-velocity floor; clamping σ would conflate the two mechanisms.

**Verification (sandbox, 10000-trial collision frequencies):**

| v [Å/ps] | predicted | observed |
|----------|-----------|----------|
| 0.0      | 100.00 %  | 100.00 % |
| 0.5      |  87.50 %  |  87.19 % |
| 5.0      |   8.75 %  |   8.62 % |
| 50.0     |   0.88 %  |   0.93 % |

Helper is covered by 7 unit tests in `test_collisions.py
::TestVelocityDependentCrossSection` and 4 sandbox checks in
`smoke_test_collisions.py`.

---

## Step 11b — Ion initial-state builder + IonCheckpoint v3 schema

### Schema bump v2 → v3

Added three fields to `IonCheckpoint` so postprocess and energy-
conservation diagnostics have everything they need:

- **`droplet_radii_angstrom: (2N,)`** — per atom, mirrors NeutralCheckpoint.
- **`mass_history_kg: (2N, T)`** — mass over time, since helium
  attachment changes mass during the run.
- **`E_dissip_eV: (2N, T)`** — cumulative energy dissipated per atom
  (matches NeutralCheckpoint convention).

Old v2 ion checkpoints cannot be loaded; the loader raises
`ValueError` with a clear message. Schema-version checks in
`tests/test_checkpoint.py` and `smoke_test_checkpoint.py` updated to
verify v3 fields round-trip correctly.

### `build_initial_ion_state` function

Takes a `NeutralCheckpoint` and produces an `IonCheckpoint` with
column 0 populated (positions, velocities, masses, energies) and
columns 1..T-1 zero-allocated for the driver to fill.

Reads from the neutral checkpoint:
- positions and velocities at column `start_id` (default -1 = last)
- `mass_kg` and `droplet_radii`

Computes at t=0:
- `E_kin = ½ m v²` (per atom, eV)
- `E_pot = ion-droplet + half-pair Coulomb` (per atom, eV)

### Bugs found in legacy MATLAB and fixed in Python

While porting, two t=0 bookkeeping bugs were found in
`vmi_sim_3d_ion_propa.m`:

**Line 289** -- E_kin formula:
```matlab
E_kin_ion(:,1) = mass_i.*(vx² + vy²)²/2/eV;
```
- Missing `vz`
- `(vx² + vy²)²` squares v² → produces `m * v⁴ / 2 / eV`

**Line 291** -- E_pot formula:
```matlab
E_pot_ion(:,1) = droplet_potential(sqrt(x² + y²) - R);
```
- Missing `z²` in the radial coordinate
- Missing the partner Coulomb term entirely (subsequent steps DO
  include it via `frog_step_ion`)

Both are "silent" because `E_pot_ion(:,1)` is only ever read for
diagnostic plotting at the end of the run -- the same pattern as
the t=0 E_pot bug we found in the neutral stage in Step 10. Python
port fixes both per principle #10. Regression tests in
`test_ion_initial_state.py::TestInitialEnergies` would catch either
bug if it returned (they verify v² and inclusion of vz, of z, and of
the Coulomb term).

### Out-of-scope features

`build_initial_ion_state` raises `NotImplementedError` at build time
if cfg requests any of:
- `effusive_dynamics`
- `single_charge_ionization_allowed`
- `additional_droplet_charges > 0`
- `highly_charged_iodine`

Both production input scripts (`single_pulse_N2000.m` and
`single_pulse_droplet_distribution.m`) leave all of these at their
default-disabled values. The early failure prevents silent wrong
physics if someone ever flips one of these on.

### Test coverage

- `test_ion_initial_state.py` -- 22 unit tests across 5 classes:
  TestApi (6), TestInheritance (5), TestInitialColumnsAreEmpty (3),
  TestInitialEnergies (4), TestScopeChecks (4), TestValidation (2)
- `smoke_test_ion_initial_state.py` -- 41 sandbox checks
- All 15 smoke test suites pass after the schema bump.

---

## Step 11c — Ion propagation step (pure function)

Mirrors the design pattern of Step 10c-ii (`propagation_step.py`):

- `IonStepState` frozen dataclass holds per-atom dynamic quantities
- `ion_propagation_step(state, *, cfg, droplet_radii, charge,
  prev_distance, rng) → IonStepState` is pure (no mutation, no I/O)

### Key design decision: mass carried in state (Option C)

Unlike neutral, ion mass changes per atom over time as helium attaches.
The cleanest approach is to make `mass_kg` part of `IonStepState` and
rebuild the `make_ion_step` closure inside the pure step function each
iteration. This keeps the leapfrog API unchanged and adds only one
function-construction per step (no per-atom Python loops).

### Step sequence

1. Leapfrog (rebuilds closure with current mass)
2. Depth into droplet
3. Per-atom cross section: v-dependent (`σ_0 · v^exponent`) if
   `cfg.sigma_dependent_on_v` (production default), else constant
4. Mode-3 collision sampling using previous step's distance
5. Apply elastic collisions (`apply_collision`)
6. Mass attachment: `rng.uniform < p_attach` AND `b_collision`
   → mass += 4 amu
7. Energy bookkeeping with **new** (post-attachment) mass
8. Return new state with `time_ps += dt_ion`

### Energy conservation result

Without mass attachment, `E_kin + E_pot + E_dissip` drift over 50
steps is **0.0022%** — at the leapfrog symplectic-error limit
(~ppm/step). With attachment enabled (production setting), recomputing
`E_kin = ½ m_new v²` after a 4-amu helium atom attaches at the atom's
current velocity overstates the true post-attachment kinetic energy by
`½ Δm v²`. The legacy MATLAB tracks this as a per-atom diagnostic
``E_mass_attach_defect`` (`vmi_sim_3d_ion_propa.m:762`).

The initial Step 11c port omitted that diagnostic. The omission was
re-evaluated during Step 11 cross-reference work and reversed in
**Step 11e** (see below): the field is now part of `IonStepState` and
`IonCheckpoint` (schema v4), and the conserved invariant is
`E_kin + E_pot + E_dissip + E_mass_attach_defect ≈ const` modulo
Verlet drift on each side.

### RNG draw pattern matches MATLAB

The mass-attachment trial is drawn for ALL atoms (not just colliders)
on every step, matching MATLAB line 727. This keeps the rng stream
deterministic in case anyone uses it for reproducibility checks
between Python and MATLAB later.

### Test coverage

- `test_ion_propagation_step.py` -- 24 unit tests across 8 classes
  (TestApi, TestFirstStep, TestReproducibility, TestEnergyBookkeeping,
   TestMassAttachment, TestEnergyConservation, TestScopeChecks,
   TestVelocityDependentSigma)
- `smoke_test_ion_propagation_step.py` -- 30 sandbox checks (all pass,
  including v-dependent σ behavior: slow atom 50/50, fast atom 1/50)
- All 16 smoke test suites pass

---

## Migration progress

See `README.md` for the live progress table.

---

## Step 11d — Full ion propagation driver

The full ion-stage driver has now been implemented in `simulation/ion.py`.

This completes the first functional Python ion-propagation path by connecting
the pieces introduced in Steps 11a–11c:

- velocity-dependent cross-section helper,
- `build_initial_ion_state`,
- pure `ion_propagation_step`,
- `IonCheckpoint`,
- `RunDirectory`.

The driver follows the same orchestration pattern as `run_neutral_propagation`
where appropriate:

1. build the initial ion checkpoint from the final or selected neutral state,
2. decide internal ion timestep count,
3. decide storage stride when needed,
4. propagate internally at the configured ion timestep,
5. write stored states into the checkpoint,
6. preserve the final state even when it does not align exactly with the
   storage stride,
7. optionally save through `RunDirectory`.

The implementation is intended to preserve existing Python corrections to
known MATLAB bookkeeping bugs rather than reintroducing them for byte-identical
legacy output.

Known corrections that remain relevant:

- ion-stage `E_kin` at `t=0` uses the full velocity vector,
- ion-stage `E_pot` at `t=0` uses the full radial coordinate,
- ion-stage `E_pot` at `t=0` includes the partner Coulomb contribution,
- Python uses the updated physical constants defined in `constants.py`.

The project is now in the ion-stage MATLAB/Python cross-reference phase.
The next goal is to validate the completed Python ion driver against small,
deterministic MATLAB reference cases before moving on to the public
single-pulse run script.

Recommended validation order:

1. ion `t=0` state copied from a tiny neutral checkpoint,
2. deterministic one-step ion propagation with collisions disabled,
3. deterministic multi-step ion propagation with collisions disabled,
4. energy-bookkeeping sanity in deterministic mode,
5. collision/statistical checks only after deterministic tests are stable.

Do not start with a full stochastic trajectory comparison, because collision
sampling, mass attachment, velocity-dependent cross sections, and RNG-stream
differences make such a comparison too entangled for a first reference test.

---

## Step 11e — `E_mass_attach_defect_eV` diagnostic ported (schema v3 → v4)

While planning the stochastic forced-event cross-reference (validation
target 5), the per-side conservation invariant
`E_kin + E_pot + E_dissip ≈ const` was found to be incorrect once
helium mass attachment is enabled: when 4 amu attaches to an atom
moving at `v`, recomputing `E_kin = ½ m_new v²` overstates the true
post-attachment kinetic energy by `½ Δm v²`. The legacy MATLAB tracks
the negative of this overstatement as a per-atom diagnostic
``E_mass_attach_defect`` at `vmi_sim_3d_ion_propa.m:762`, which makes
`E_kin + E_pot + E_dissip + E_mass_attach_defect` conserved up to
Verlet drift. The original Step 11c port omitted this term and
documented the omission as a "known approximation"; per project
principle #10 we reversed that omission so the conservation invariant
holds exactly on each side, which is a prerequisite for using the
invariant as a cross-language sanity check during the stochastic
comparison.

Changes:

- `IonCheckpoint` schema bumped v3 → v4. New field
  `E_mass_attach_defect_eV: (2N, T)` (per-atom, cumulative, eV).
- `IonStepState` gains the matching scalar-per-atom `(2N,)` field.
- `ion_propagation_step` accumulates the per-step increment

  ```
  dE_defect = -½ (m_new − m_old) · v_post² · 100²/eV
  ```

  using the post-collision, post-attachment velocity (verbatim from
  the legacy MATLAB).
- `build_initial_ion_state` allocates the array and seeds it with
  zeros at t=0.
- `_NUM_2N_T_ARRAYS_ION` bumped 12 → 13 in the storage-stride budget.
- Documentation updated in `docs/checkpoint_module.md`,
  `docs/ion_propagation_step_module.md`,
  `docs/ion_initial_state_module.md`, and `docs/ion_module.md`.

Verification:

- Full pytest suite (325 tests) passes after the change.
- The two existing deterministic cross-reference scripts
  (`scripts/cross_reference/ion_t0_state/` and
  `scripts/cross_reference/ion_multistep_no_collision/`) re-run to
  byte-identical output: their scenarios run with
  `mass_attach_probability = 0`, so the new field stays zero
  throughout. This confirms the diagnostic is purely additive when
  attachment is disabled.
- Forced-event smoke test (1 molecule, σ inflated, p_attach = 1)
  shows the defect accumulating negatively as expected with each
  attachment, partially compensating the spurious E_kin gain.

Old `IonCheckpoint` `.npz` files on disk become unreadable; the loader
reports `schema_version=3, this code expects 4. Re-run the simulation`,
which is the documented behaviour for schema bumps. No checkpoints are
present in `data/` — only in test scratch directories.

---

## Step 11f — Ion-driver MATLAB/Python cross-reference complete

The ion-stage MATLAB/Python cross-reference phase described in
`CLAUDE.md` is now complete.

Completed validation targets:

1. Ion `t=0` state copied from a tiny neutral checkpoint.
2. One deterministic ion step with collisions disabled.
3. Several deterministic ion steps with collisions disabled.
4. Energy bookkeeping in deterministic mode.
5. Collision/statistical behavior after the deterministic comparisons
   were stable.

Reference artifacts live under:

- `scripts/cross_reference/ion_t0_state/`
- `scripts/cross_reference/ion_multistep_no_collision/`
- `scripts/cross_reference/ion_stochastic_forced/`

The deterministic checks cover targets 1 through 4. The forced stochastic
case covers target 5 with collision and mass-attachment behavior isolated
enough to compare bookkeeping and event effects without starting from a
full production stochastic trajectory.

---

## Step 12 — Public single-pulse run script

Implemented `scripts/run_single_pulse.py` as the first user-facing entry point
for the validated neutral + ion pipeline.

The script:

- builds the canonical config from `single_pulse_N2000`,
- accepts small-run overrides (`--num-molecules`, `--seed`, and
  `--ion-simulation-time`),
- creates a `RunDirectory`,
- runs `run_neutral_propagation`,
- runs `run_ion_propagation` from the neutral checkpoint,
- writes `cfg.json`, `neutral.npz`, and `ion.npz` through the existing
  run-directory conventions,
- refuses to overwrite existing run artifacts unless
  `OVERWRITE_EXISTING_RUN = True`,
- prints concise progress and output paths.

The saved `cfg.json` keeps the preset's ion-stage
`sigma_dependent_on_v=True`. This is acceptable because the Python neutral
driver uses the neutral cross-section field directly, while the ion driver is
the stage that reads `cfg.sigma_dependent_on_v`. No stage-specific physics
config mutation is needed for the current Python APIs.

Usage example:

```bash
python scripts/run_single_pulse.py --run-dir results/single_pulse_test --num-molecules 10 --seed 123
```

Tests:

- Added `tests/test_run_single_pulse.py`.
- `pytest tests/test_run_single_pulse.py -q` passes: 3 tests.
- `pytest tests/test_run_directory.py tests/test_neutral.py tests/test_ion.py tests/test_run_single_pulse.py -q` passes: 63 tests.

Step 13 followed this by adding the first HeDFT loading and trajectory
comparison path in `postprocess/`.

---

## Step 13 — Post-processing comparison path

The first Python post-processing path is now implemented. The new code imports
the useful numerical and plotting pieces from the legacy HeDFT comparison
workflow without expanding into Abel inversion or full experimental VMI
analysis.

Implemented package modules:

- `i2_helium_md/postprocess/hedft_loader.py`
  - loads normalized 8-column HeDFT/TDDFT reference CSVs with header
    `Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance`,
  - infers the droplet radius from filenames such as `9A_All_Data.csv` and
    `18A_All_Data.csv`,
  - supersedes the split legacy 9 A import helpers for the current Python
    comparison contract.
- `i2_helium_md/postprocess/compare_trajectories.py`
  - computes mean MD I-I distance vs. time from an `IonCheckpoint`,
  - interpolates MD data onto the HeDFT time grid over the overlap interval,
  - reports RMSE and mean MD/HeDFT ratio,
  - applies the same overlap/interpolation contract to I1 and I2 velocity
    magnitudes.
- `i2_helium_md/postprocess/velocity_distribution.py`
  - loads exported VMI reference CSVs (`vmi_iplus_he.csv`,
    `vmi_iplus_gas.csv`),
  - computes final-velocity histograms from mass-selected, outside-filtered
    ion checkpoint atoms for simple simulation overlays.
- `i2_helium_md/postprocess/__init__.py`
  - exposes the post-processing dataclasses and helper functions as the public
    package API.

Implemented scripts:

- `scripts/plot_hedft_comparison.py`
  - loads an existing `RunDirectory` ion checkpoint,
  - loads the normalized HeDFT reference,
  - recreates the distance-trajectory figure,
  - recreates the velocity-vs-time figure,
  - optionally adds the VMI/reference velocity-distribution tile.
- `scripts/post_processing_comparison/compare.py`
  - remains as imported VMI-reference verification context for
    `vmi_iplus_he.csv` and `vmi_iplus_gas.csv`.

Reference data currently expected under `data/reference/`:

- `9A_All_Data.csv`
- `18A_All_Data.csv`
- `vmi_iplus_gas.csv`
- `vmi_iplus_he.csv`

Documentation added for the new modules and script:

- `docs/hedft_loader_module.md`
- `docs/compare_trajectories_module.md`
- `docs/velocity_distribution_module.md`
- `docs/plot_hedft_comparison_script.md`

Tests:

- `tests/test_hedft_loader.py`
- `tests/test_compare_trajectories.py`
- `tests/test_velocity_distribution.py`
- `tests/test_plot_hedft_comparison_smoke.py`

Verification in this documentation-update pass:

```bash
python -m pytest tests/test_compare_trajectories.py tests/test_hedft_loader.py tests/test_velocity_distribution.py tests/test_plot_hedft_comparison_smoke.py -q
```

Result: 36 tests passed.

Remaining work is no longer the first implementation of Step 13, but rather
post-processing validation and interpretation: run the numerical comparison
API on the current production run, record the RMSE/ratio values that matter,
and decide whether any outputs should become stable reference diagnostics.
Keep further work focused; do not expand into Abel inversion, pump-probe, or
full experimental VMI analysis without an explicit request.

---

## Run script generalized across input presets

`scripts/run_single_pulse.py` now supports a named `INPUT_PRESET` setting
instead of being hardwired to `single_pulse_N2000()`.

Supported presets:

- `single_pulse_N2000`
- `single_pulse_droplet_distribution`

`RUN_SIZE` keeps its previous meaning:

- `smoke` and `custom` apply `NUM_MOLECULES`, `SEED`, and `ION_TIME_PS` as
  overrides on top of the selected preset,
- `production` uses the selected preset exactly, except for explicit
  `PRODUCTION_*` overrides.

This keeps the editable-script workflow intact while letting the same public
pipeline run either migrated MATLAB input file.

Verification:

```bash
python -m pytest tests/test_run_single_pulse.py -q
```

Result: 8 tests passed.

---

## Legacy MATLAB live-debug and paper post-processing reproduction

Reproduces three classes of legacy MATLAB output that the Python port had
not yet covered:

1. **Neutral energy-balance debug figure**
   (`vmi_sim_3d_neutral_propa_HeDFT_mimic.m:965`) -- sum-over-atoms
   `E_kin`, `E_pot`, `E_dissip`, total `E_system` over the neutral time
   axis. All inputs were already in the existing `NeutralCheckpoint`,
   no physics changes required.
2. **Ion energy-balance debug figure**
   (`vmi_sim_3d_ion_propa.m:898`) -- per-molecule `E_kin`, `E_pot`,
   `E_dissip`, `E_mass_attach_defect`, `E_system`. All inputs already
   in `IonCheckpoint`, no physics changes required.
3. **Ion temperature-diagnostic figure**
   (`vmi_sim_3d_ion_propa.m:683` / `:883`) -- per-step
   `[<T'/T>_actual, <T'/T>_from_mass_ratio, <theta_lab>]` averaged over
   the colliding atoms. Required adding a new ion-checkpoint field.
4. **Minimal `post_process_single_pulse_paper_v3.m` reproduction**
   -- radial velocity comparison (already covered by
   `plot_experimental_comparison.py`), simulated azimuthal phi
   histogram, final ion mass spectrum, combined PDF export named
   `compare_simulation_and_measurement.pdf` to match legacy filename.

### IonCheckpoint schema v4 -> v5

Added `temperature_diagnostic: (num_steps, 3)` to `IonCheckpoint` -- the
only ion-checkpoint array whose leading dimension is `num_steps` rather
than `2N`. NaN where no atom collided in a stored step. Validation lives
outside the `(2N, num_steps)` loop in `_validate_against_cfg`.

`_load_checkpoint` already enforces the field set, so any pre-v5 file
fails at load time with the standard "schema_version=4, this code
expects 5; rerun the simulation" error. Existing `data/runs/*/ion.npz`
files must be regenerated.

### Capture path

`physics/collisions.py` now exposes a `CollisionDiagnostics` namedtuple
plus the recipe helper `temperature_diagnostic_from_collision`.
`apply_collision` accepts a keyword-only `return_diagnostics: bool =
False`; the default keeps the 4-tuple shape so all existing neutral-side
callers and tests remain unaffected. `IonStepState` carries an optional
`temperature_diagnostic` field set by `ion_propagation_step` from the
returned `CollisionDiagnostics`. `simulation/ion.py` writes the row into
`IonCheckpoint.temperature_diagnostic[stored_step_idx, :]` for every
stored step, mirroring the MATLAB
`diagnostic_array(1:reduction_timesteps:end, :)` downsampling.

### New post-processing modules

- `i2_helium_md/postprocess/energy_balance.py` -- pure recipe helpers:
  `neutral_energy_totals`, `ion_energy_totals`, `phi_histogram`,
  `mass_spectrum`. `mass_spectrum` uses bin edges at half-integer amu
  so 127.0 lands in the bin centred on 127 rather than on a bin edge.
- `i2_helium_md/postprocess/_smoothing.py` -- shared
  MATLAB-style `movmean` and trace normaliser; both
  `plot_experimental_comparison.py` and the new paper-figure script
  import from here instead of duplicating the helpers.

### New scripts (manual invocation)

Under `scripts/post_processing/`:

- `plot_neutral_energy_balance.py`
- `plot_ion_energy_balance.py`
- `plot_ion_temperature_diagnostic.py`
- `plot_paper_figure.py`

Each loads from a `RunDirectory` and writes under
`<run>/figures/`. None of them auto-run from
`scripts/run_single_pulse.py` -- the run script remains untouched.

### Out of scope (deferred)

The polar-VMI panels of `post_process_single_pulse_paper_v3.m` (cos^2
angular anisotropy fit, beta(v) function, 3-D `surf` of polar VMI image)
require a 2-D polar VMI image not present in `data/reference/`.
CLAUDE.md flags "full experimental VMI interpretation" as out of default
scope, so these panels are intentionally skipped pending an explicit
request and additional reference data.

### Verification

```bash
python -m pytest tests/test_energy_balance.py tests/test_collisions.py \
    tests/test_ion_propagation_step.py tests/test_ion.py \
    tests/test_ion_initial_state.py -q
```

Result: 106 tests passed. The schema bump was also smoke-tested with a
small end-to-end run (`N=50`, 2 ps ion stage): 197/200 stored steps
recorded a non-NaN diagnostic; `<T'/T>_from_mass_ratio` clusters at
~0.943, exactly the analytical asymptote `(1 + rho^2) / (1 + rho)^2`
for `rho = 127/4` (iodine on helium).
