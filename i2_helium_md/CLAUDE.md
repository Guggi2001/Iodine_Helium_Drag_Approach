# CLAUDE.md

Project guidance for Claude-style agents before editing or reviewing this
repository. This file holds **durable rules**, not live status — status,
decisions, and history live in the dedicated docs linked below.

## Project Snapshot

Python port of a legacy MATLAB molecular-dynamics codebase for iodine /
iodine-ion dynamics in helium nanodroplets. Main package: `i2_helium_md`.
Legacy MATLAB reference: `legacy_matlab_repository/`.

The MATLAB→Python transfer is complete: neutral + ion propagation, the
single-pulse presets, run-directory artifacts (`cfg.json`, `neutral.npz`,
`ion.npz`), cross-reference validation, the full in-scope post-processing
surface, and the frozen experimental reference exports under `data/reference/`.

The current phase is **drag-model physics**: replace the hard-sphere collision
model (`physics/collisions.py`) with a TDDFT-calibrated drag-force model for I⁺
in a helium bubble. The post-processing surface is the stable comparison layer
for the new physics. This phase is the scoped exception to the "do not change
collision physics" forbidden-list rule.

**Live status and full history are not duplicated here — consult the docs
(actual paths shown; root unless a directory is given):**

Repo root:
- `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` — live mass-model detail (locked mechanism, SQ1–SQ3 integrator/mass-jump split, R/A/OQ registers)
- `DRAG_PORT_DESIGN_DECISIONS.md` — frozen architecture choices
- `CALIBRATION_MAP.md` — cross-doc parameter index (Sourced/Bounded/Free/Derived classes, Tier anchors, identifiability)

`docs/drag_port/Tier0/`:
- `drag_migration_log_tier0.md` — Slices 1–4 + Tier-0 + extraction decision / withdrawal history
- `TIER0_FINDINGS.md` — Tier-0 verdict
- `METHOD_B_trajectory_matching_extraction.md` — extraction method + held-out validation
- `SLICE{1..4}_GOALS_*.md` — per-slice goal docs

`docs/drag_port/Tier1/`:
- `drag_migration_log_tier1a.md` — Tier-1a decision + delivery history (delivered)
- `TIER1A_IMPLEMENTATION_PLAN.md` — Tier-1a build plan (anchored kinematic mass-dynamics validation; delivered)

`docs/drag_port/Tier2/` (current phase — the active goal documents):
- `drag_migration_log_tier2.md` — Tier-2 decision + delivery log
- `TIER2_IMPLEMENTATION_PLAN.md` — program overview (11 slices / 6 phases A–F)
- per-phase detail plans: `TIER2_PHASE_{A,B,C,D,E,F}_IMPLEMENTATION_PLAN.md`
- `RESEARCH_QUESTIONS.md` — consolidated open physics questions (RQ1–RQ6:
  E_int(0) provenance, per-shed ε, suppressed fate/bare peak, ladder bottom,
  µs-flight channels, s_eff cross-check) — entry document of the
  literature-research / cross-validation phase (active since 2026-07-09)
- `TIER2_SENSITIVITY_ATLAS_PLAN.md` — **current goal document** (designed
  2026-07-23): the sensitivity-atlas study program
- `TIER2_PARAMETER_INFLUENCE.md` — compact per-knob influence reference
  (atlas D0; living doc — its GAP markers are the atlas targets)

`docs/matlab_port/`:
- `PHYSICS_BASELINE.md` — MD baseline

Compact state (verify against the log before relying on it): production law is
`shared_pure_cubic` (`γ = g·b·v²`), both presets wired to the shared bundle; drag
form is settled pure-cubic in-band (the production arbitration adds the
`capped_cubic` tail, cap **v_c ≈ 7.25–7.5** / `p_tail = −1` — the §4w
joint-landing basin); **Tier 0 is complete** (18 Å clean pass, 9 Å
non-radial flag). **Tier 1a is delivered** — anchored kinematic mass-dynamics
validation: the He shell schedule `n(t)` is read from the 9 Å TDDFT loss curve and a
controlled `fixed` vs `anchored_discrete` A/B was built (variable-mass integrator
SQ1–SQ3, v6 checkpoint, four-term ledger; continuous-velocity shed is the physical
path, cold-shed retained as a diagnostic bound). **Tier 2 is in its
OPTIMIZATION stage (build program A–F complete, F5 reconciled 2026-07-22; the
tier stays ACTIVE)** — the *generative* `biphasic_energy_gated` mass mechanism
(Poisson pickup + energy-gated RRK evaporation + `E_int` reservoir + Newton
cooling, **5-term** invariant, `IonCheckpoint` **v7** with `E_int`) is built
and arbitrated against the experimental I⁺Heₙ size distribution and fragment
mean-KE at production 2.70 eV. **Standing production point: `finc1v725`**
(`capped_cubic` v_c 7.25 / τ 3.2 / E₀ 0.27 / c1 rq4graded ladder / Landau-on
0.58), statistical reference = the pooled **N = 5000 battery** (5 × N = 1000
fresh seeds, findings §4cc): histogram-level landing seed-robust (n₁_solv
0.243, midHot 1.014, W₁ 0.571 ± 0.04). Figure surface: `plot_run_summary.py`
(MD window) + `plot_detection_summary.py` (detected ensemble, 14 sections) +
the pooled figures container. **Current goal (2026-07-24): the sensitivity
atlas** (`TIER2_SENSITIVITY_ATLAS_PLAN.md` + the D0 reference
`TIER2_PARAMETER_INFLUENCE.md`; results in
`TIER2_SENSITIVITY_ATLAS_FINDINGS.md`) — understanding-driven OAT influence
mapping around finc1v725. **Executed so far:** stage 2a (D4 Step 1
Method-B form table, artifact reuse, oracle bit-exact; Padé excluded by
arithmetic, subtractive gate fired/fit pending, sub-2.5 Å/ps TDDFT-blind)
and the §6.6 quadratic counterfactual **twin + MD spot-check → outcome
(b), MD-measured: the Tier-2 landing is form-blind** (an lq system with
its own artifacts lands all three N = 500 cells incl. the off-needle
control; the twin's needle was a frozen-chord artifact; form authority
rests solely on the Tier-0 traces where lq stays held-out-rejected;
the landing is never evidence for any form). New RQ11 signal
(NB-RQ11-12): the deep-KE lever sits in the **5–9 Å/ps mid-band**, not
sub-2.5 — with the lq E_bind-0.048 confound caveat pending. **§6.7
item 1 (lq sanity battery) EXECUTED 2026-07-24** (5 × N = 1000, qcc,
seeds 20260722–26 paired-by-seed vs `bigc1v725s{1..5}`; generator
`gen_tier2atlas_lqbattery.py`; oracles bit-exact): pre-registration
**3 MET / 2 MISSED** — histogram form-blindness CONFIRMED at battery
scale (W₁, midHot), but the §6.6 KE "lands better" is **REFUTED**
(paired lq χ²_med worse on 5/5 seeds, Δ+13.3±10.0; the §6.6 halving was
a seed-20260721/N = 500 artifact) and lq **over-suppresses** (supp
0.213 vs 0.187, Δ+0.026 ~6σ; NB-RQ11-13). Mid-band deep-KE warming
persists in direction (~300 counts) but buys no better KE landing.
**§6.7 item 2 (E_bind pair-separation scan) EXECUTED 2026-07-24**
(initial N = 500 single-seed + firm-up N = 1000 × 3 paired seeds;
generators `gen_tier2atlas_ebindscan.py`/`ebindseeds.py`; only the
climbed well overridden under the unvalidated-binding hatch, lq drag
stamp unchanged; oracle bit-exact): **the over-suppression is DISCHARGED
to the FORM, not the well** (matched-well 0.1168 Δsupp +0.034 ~8σ; well
adds only +0.010), and item-1's histogram match is a **(form, well)
co-compensation** — at a matched well lq is colder (midHot Δ−0.070 ~35σ)
and smaller (n̄ Δ−0.57); **W₁ alone stays form-blind even at matched
well**. Single-seed W₁ reads were seed-noise; χ²_med inconclusive even at
N = 1000; deep well 0.154 over-retains (handover-guard trips, 2/3 cells).
Net: forms are physically distinguishable; "form not limiting" holds only
for coarse W₁, not KE/size/fate. **Next (plan §6.7 item 3, designed NOT
triggered):** leave D → Axis A geometry grid → D2b → E₀/τ curves. RQ3/RQ5
reads and
the margin-3 Å pinned convention (I88) stay open in-tier. Tier-3 noise
stays next (second-moment under-dispersion; stubbed behind its enum,
NOT retired). New drag-program code stays behind the
`[PROCEED TO IMPLEMENTATION]` trigger.

## Current Scope

In scope:

- single-pulse neutral and ion dynamics,
- 9 A and normalized 18 A HeDFT comparison inputs (in `data/reference/`),
- VMI reference loading and final-velocity histogram helpers,
- consolidated post-processing diagnostics from finished run directories,
- authentic reproduction of legacy post-processing figures where reference data
  and run outputs are available,
- focused MATLAB/Python reference validation,
- drag-model physics calibrated against TDDFT for I⁺ in the helium bubble
  (replacing the hard-sphere collision model).

Out of scope unless the user explicitly asks: pump-probe support; effusive /
gas-phase dynamics; Abel inversion; full experimental VMI image interpretation;
broad experimental VMI analysis beyond current overlays; new physics branches;
broad refactors; live-debug 3D animations and visualization-only MATLAB utilities.

## Project Quality Principles

Use these for every porting decision, code review, and cleanup:

1. No duplicate implementations of the same physics. Shared conversions,
   formulas, and constants belong in `constants.py` or the appropriate shared
   module.
2. No dead code. Remove unused imports, commented-out blocks, and speculative
   branches. *Scoped exception:* the drag-config fields declared-but-not-yet-read
   (the rule-2 carries are recorded as prose entries in the **active phase's** drag
   migration log — currently `docs/drag_port/Tier2/drag_migration_log_tier2.md`;
   each carry is retired by the slice/tier that activates it. Convention: a field
   whose only reader is a config-load *guard* is still a carry — "guard-live" is
   not "physics-live").
3. Encode units and conventions in names: `mass_kg`, `time_ps`, `T_particles_K`,
   `R0_GS_angstrom`, etc.
4. Validate early and fail loudly. Wrong shape, unsupported collision mode,
   invalid type, or non-overlapping time axes should raise clear errors.
5. Public functions need docstrings with units, shapes, edge cases, inputs,
   and outputs.
6. Organize modules by concern: `physics/` is science, `sampling/` is
   randomness, `simulation/` is orchestration, `postprocess/` is analysis of
   finished runs.
7. Tests document intended behavior, units, tolerances, and known MATLAB/Python
   deviations.
8. Audit after refactors for dead imports, duplicate physics, and doc/code drift.
9. For any known reference output, literal transliteration comes before clean
   refactor. Reproduce the MATLAB behavior first, then refactor once verified.
10. Do not preserve bad legacy behavior for byte identity. Python uses corrected
    modern constants and fixes known MATLAB bookkeeping bugs unless the user
    explicitly asks for legacy behavior.

## Architecture Rules

Use `SimConfig` instead of globals (consolidated in `i2_helium_md/config.py`);
functions that need parameters receive `cfg: SimConfig` explicitly.

Use preset functions instead of scripts:

```python
from i2_helium_md import (
    single_pulse_N2000,
    single_pulse_N2000_18Angst,
    single_pulse_droplet_distribution,
)
```

Use `RunDirectory` for simulation artifacts. A run directory is self-describing:
`cfg.json`, `neutral.npz`, `ion.npz`, optional `figures/`.

Checkpoint I/O rules: checkpoints are explicit dataclasses; every checkpoint has
`schema_version`; incompatible versions fail at load time; no constants are saved
in checkpoints (config belongs in `cfg.json`); load with `allow_pickle=False`;
shape validation uses `cfg` when available.

Avoid changing checkpoint schema, physical constants, neutral propagation, ion
propagation, collision physics, or `scripts/run_single_pulse.py` behavior unless
a focused test reveals a real bug or the user explicitly asks.

## Drag-Model Port

This section governs the active drag-model phase. It sits above the
post-processing rules, which now describe the stable **comparison layer**.

### Working method

- **Physics-Definition vs. Software-Implementation boundary is strict.** Default
  to mathematical formulation (LaTeX) and logic trees. Do **not** write
  implementation code until the user gives the explicit trigger
  `[PROCEED TO IMPLEMENTATION]`. Discussion and specification come first.
- **Every model choice stays interchangeable behind a `SimConfig` enum** so the
  empirical cross-check against the TDDFT / VMI references can select among
  hypotheses. No drag form, noise model, mass scenario, or gate is hard-wired;
  each lives behind its own enum surface (`DRAG_PORT_DESIGN_DECISIONS.md`).
- **Before proposing a governing equation, perform and display a strict
  dimensional analysis.** Define every parameter's units (forces in `amu·Å/ps²`,
  velocity in `Å/ps`, etc.) and reject any formulation whose units do not balance.
- When brainstorming a model transition (discrete collisions → continuous drag),
  list the physical trade-offs explicitly: lost degrees of freedom, violated
  conservation laws, thermodynamic changes (e.g. loss of the Langevin
  fluctuation channel).

### Friction convention (unified — do not reintroduce the second one)

- **`γ(v)` is a force coefficient**, units **amu/ps**, defined
  `γ(v) = |F_drag(v)| / v`. The friction force is `γ(v)·v` — **no leading `m`**.
- The friction **rate** is `γ/m` [1/ps] and appears **only** inside the BAOAB
  damping exponent `e^(−γ·dt/m)`.
- The FDT noise amplitude `√(2·γ·k_B·T_eff)` uses `γ` [amu/ps] directly.
- Consequence: the drag-physics module is **mass-agnostic**. It never takes `m`.
  Mass enters only at the integrator's O-step, as one explicit division by `m(t)`.
  Do not pass mass into the drag module.

### Status, slices, and the field exception

All four implementation slices are delivered and reviewed (Slice 1
`physics/drag.py`, Slice 2 `physics/baoab.py`, Slice 3 `SimConfig` drag surface +
config-load guard + coefficient loader, Slice 4 ion-driver rewiring incl. the
`SLICE4_FIX` mass-consistency fix). Per-slice delivery records, the drag-config
declared-but-unread field exception (rule 2), the full Tier-0 diagnosis history
(including withdrawn readings), and the Method-B / shared-form / form-discrimination
records all live in `drag_migration_log_tier0.md` (Tier-1a history in
`drag_migration_log_tier1a.md`). **Consult the logs before assuming any status.**

### Validation hierarchy (sequential, not simultaneous)

Parameters are entangled; validate in tier order, fixing each tier's winner
before introducing the next unknown:

- **Tier 0 — COMPLETE.** Drag form, deterministic, fixed mass, in-window TDDFT
  traces. Locked: `shared_pure_cubic` (18 Å clean pass, 9 Å non-radial flag).
- **Tier 1a — DELIVERED.** Anchored kinematic mass-dynamics validation.
  The shell schedule `n(t)` is *anchored* to the 9 Å TDDFT loss curve
  (~21→19→14 He), not generated; the run is a controlled `fixed` vs
  `anchored_discrete` A/B that tests whether drag + variable mass reproduce `R(t)`,
  `|v(t)|` and whether the four-term §2.9 ledger closes. **OQ-independent** —
  picture/ladder/κ/ν/s are bypassed by the anchor. Built the variable-mass
  integrator upgrade (SQ1–SQ3: drag O-step under `m(t)`; post-jump `m⁺`), with the
  continuous-velocity shed as the physical path and cold-shed retained as a
  diagnostic bound. The predictive shell-timing variant ("1b") is **rejected**
  (TDDFT is not ground truth — experiment arbitrates at Tier 2). Plan + slices:
  `TIER1A_IMPLEMENTATION_PLAN.md`.
- **Tier 2 — ACTIVE, optimization stage (build A–F complete, F5 reconciled
  2026-07-22). Terminal I⁺(He)ₙ size distribution** vs. experimental detector
  data — the only observable that separates the mass scenarios *and* arbitrates
  the two genuinely-free knobs (ladder shape + electronic picture). The biphasic
  *generative* mechanism (Poisson pickup + energy-gated RRK evaporation +
  `E_int` reservoir + Newton cooling, **5-term** invariant; schema **v7**,
  `E_int` field) is built and arbitrated: standing point finc1v725, pooled
  N = 5000 battery reference, histogram-level landing seed-robust. In-tier
  optimization targets: RQ11 (deep-bin KE slope), RQ3, RQ5, the margin pin.
  Plan + slices: `TIER2_IMPLEMENTATION_PLAN.md`; F5 reconciliation entry in
  `drag_migration_log_tier2.md`.
- **Tier 3 — next (open, not started). Ensemble second moments** (noise) vs.
  VMI references — the under-dispersed cov/VMI panels are its measured
  motivation; stays stubbed behind its enum until its own plan + trigger.

A `mass_scenario`↔`drag_coefficients` consistency guard is enforced at config-load
(§6.5): constant-mass coefficients are self-consistent only with
`mass_scenario=fixed`. Both the Tier-1a `anchored_discrete` run and the Tier-2
`biphasic` production run trip this guard structurally and run under
`allow_inconsistent_mass_pairing=True`, defended by the §6.6 mid-window argument
(`m≈19 He = m_eff` mid-window; the `n=21`/`n=14` ends sit in the §6.7 free-zone).
Histogram comparisons (Tier 2/3) default to the Wasserstein metric.
Current tier status: `drag_migration_log_tier2.md` (active) /
`drag_migration_log_tier1a.md` / `TIER0_FINDINGS.md`.

## Post-Processing Comparison Layer

Prefer current Python APIs over ad hoc scripts:

1. Load a finished run with `RunDirectory`.
2. Load HeDFT references with `load_hedft_trajectory`.
3. Compute numerical diagnostics with `compare_distance` and
   `compare_velocity_magnitude` before changing trajectory plots.
4. Use `velocity_distribution.py` for VMI reference loading and mass-selected
   final-velocity histograms.
5. Use the focused post-processing helpers instead of rolling new histograms
   (`energy_balance.py`, `polar_velocity.py`, `velocity_2d.py`,
   `pair_correlation.py`, `time_resolved.py`, `boltzmann_overlay.py`).
6. Use `scripts/post_processing/plot_run_summary.py` for every in-scope
   diagnostic from a finished run.
7. Keep plot changes local to `scripts/post_processing/` unless a package API
   change is actually needed.
8. Add or update focused pytest coverage when behavior changes.

## General Request Handling

For investigation / audit / inspect / compare / explain requests: do not edit
files; read relevant docs and tests; report files inspected, conclusions, and
uncertainties; suggest the smallest safe next edit.

For implementation or fix requests: run `git status` first; preserve unrelated
user changes; make the smallest coherent change; do not touch unrelated files;
add or update tests when behavior changes; run relevant tests; show changed files
and remaining risks.

During the drag-model phase, "implementation" requests still require the
`[PROCEED TO IMPLEMENTATION]` trigger before any code is written.

## Post-Processing Workflow

Keep HeDFT and experimental workflows separate (the legacy combined MATLAB figure
mixed run settings):

- HeDFT trajectory comparison uses a HeDFT-comparison run; experimental VMI
  distribution comparison uses the experimental-condition run.
- Do not plot the HeDFT velocity-vs-time panel from an experimental run, or
  experimental VMI distributions from a HeDFT run.

Script entry points live in `scripts/post_processing/`. Each is configured
through its own `# USER SETTINGS` block (run dir, reference paths); read the
script header for invocation rather than duplicating it here.

The post-processing port is complete — changes are bug-fix-only. On any
disagreement with the legacy MATLAB figure, follow the **legacy-first rule**:
(1) inspect the exact MATLAB recipe, (2) literal-port the normalization / binning
/ smoothing / filtering / fit behavior, (3) verify numerically or visually,
(4) only then refactor. Do not read the full legacy plotting stack unless a
specific discrepancy requires it, and do not add out-of-scope analysis (Abel
inversion, pump-probe, effusive comparison, full VMI image interpretation).

## Data Contracts

Normalized reference data lives in `data/reference/`. Expected files:
`9A_All_Data.csv`, `18A_All_Data.csv`, `vmi_iplus_he.csv`, `vmi_iplus_gas.csv`.
HeDFT trajectory CSV header:
`Time_ps,V1_mag,V2_mag,V1_z,V2_z,V1_x,V2_x,R_distance`. (The earlier split 9 A
files `data_vabs2.csv` / `R1-R2.csv` are provenance, not the loader contract.)
Do not use hardcoded absolute MATLAB paths; verify whether an existing
`data/reference/*.csv` is VMI or HeDFT data before repurposing it.

When a legacy MATLAB script processes experimental data, that MATLAB path is the
source of truth: adapt it to export a small, inspectable reference (CSV preferred;
`.mat`/`.npz` matrix + JSON sidecar for 2-D VMI images), save under
`data/reference/`, and keep the exporter under `data/reference/scripts/`
(precedent: `export_vmi_reference_data.m`). Python loads and compares against the
export rather than reimplementing opaque extraction. Document full provenance:
source script/function, measurement IDs or input files, processing steps,
calibration/scaling factors, output columns + units, toolbox requirements, random
seed / molecule count / timestep count / enabled physics, and any known MATLAB
bug vs. intentional Python correction. For 2-D VMI image export specifics (axis /
grid conventions for `pcolormesh`, sidecar fields), follow the precedent exporter
and existing sidecars — do not invent a new format.

Reference data must be small, inspectable, reproducible (small CSV / JSON / NPZ /
text). Do not commit large checkpoints, large MAT files, generated figures, temp
debug output, or full run directories.

## Known Plotting Conventions

The exact numeric recipes — HeDFT velocity-trace cap (~30); experimental velocity
histogram bins / 15-bin moving mean / display range; the polar
`f(φ) = a + b·cos²(φ−φ0)` fit and `β = 2b/(2a+b)` recovery; angular pair-covariance
`θ = arctan2(vx,vy)+π` with zeroed diagonal; and the Boltzmann overlay
normalization — are encoded in the post-processing helpers and scripts that own
them. Those modules are the source of truth. If a panel disagrees with the legacy
figure, apply the legacy-first rule above before editing.

## Scientific-Code Caution

Clean code is not automatically correct physics. Before changing a formula, unit
conversion, force sign, random sampler, normalization convention, indexing
convention, or draw order: (1) locate the corresponding MATLAB source or previous
Python test, (2) explain the convention, (3) add a regression test or focused
numerical check, (4) then edit.

Do not start validation with a full stochastic trajectory comparison — too many
effects are entangled. Validate in this order:

1. direct formula comparison,
2. shape and unit checks,
3. one-step deterministic comparison,
4. multi-step deterministic comparison,
5. energy bookkeeping comparison,
6. stochastic statistical comparison,
7. full driver smoke test.

For deterministic tests, disable stochastic features (collisions off, mass
attachment off, fixed initial state, fixed molecule count, few timesteps, fixed
droplet radius, fixed seed). For the drag-model port the deterministic
configuration is the Tier-0 setup: `mass_scenario=fixed` at `m_eff`, noise zero,
`linear_cubic` drag. The seven-step order maps onto the tier hierarchy:
formula/shape/one-step/multi-step are Tier 0–1; energy bookkeeping closes the
§2.9 invariant; statistical comparison is Tier 2–3. For stochastic tests, prefer
distribution and moment checks over exact trajectory matching unless RNG identity
is guaranteed.

## Known MATLAB Bugs Not To Reproduce

Do not force Python to match these known MATLAB bookkeeping bugs:

- neutral-stage `E_pot` at `t=0` omitted partner Morse contribution,
- ion-stage `E_kin` at `t=0` used an incorrect velocity expression,
- ion-stage `E_kin` at `t=0` omitted `vz`,
- ion-stage `E_pot` at `t=0` omitted the `z` coordinate,
- ion-stage `E_pot` at `t=0` omitted the partner Coulomb term.

Tests should state whether Python is expected to match MATLAB, match within known
constant/unit differences, statistically match, or intentionally differ because a
MATLAB bug was corrected.

## Testing

Default: `pytest`. Python may not be on PATH; the absolute interpreter that works:

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q
```

After editing code: run the narrowest relevant test first, then broader tests if
multiple modules are affected; report exactly which tests ran; report failures
honestly; do not claim correctness without tests or a focused numerical check.

For post-processing changes, cover whichever apply: HeDFT CSV loader, comparison
with tiny synthetic checkpoints, missing-file / non-overlapping-time-axis
validation, VMI reference loader, mass-filtered final-velocity histogram,
non-interactive plotting smoke. Do not generate figures or production-sized
checkpoints in tests.

Justify numerical tolerances (tight for analytical ports; FD-step-consistent for
finite-difference forces; sample-size-based for Monte Carlo; constant-update
allowance for MATLAB/Python comparisons; a comment for any loose tolerance). Do
not update expected values blindly — first determine whether the code is wrong,
the test encoded a legacy bug, the tolerance is unreasonable, the model changed, a
constant difference is expected, the stochastic test is under-sampled, or the
reference came from the wrong MATLAB path.

## Forbidden Without Explicit User Approval

- deleting reference data,
- changing physical constants,
- changing checkpoint schema,
- changing random-number draw order,
- changing default simulation scope,
- changing neutral or ion propagation physics without a focused bug,
- changing collision physics without a focused bug,
- broad refactors,
- optimizing performance by changing numerical behavior,
- implementing out-of-scope MATLAB paths.

**Active scoped exception:** the user has approved replacing the hard-sphere
collision model with a TDDFT-calibrated drag model for I⁺ in the helium bubble.
This relaxes the "changing collision physics" rule only for that work. The
neutral driver, checkpoint schema, RNG draw order, default simulation scope, and
physical-constants table remain off-limits.

## Reporting

For post-processing work, report: (1) MATLAB files inspected, (2) Python files
changed, (3) data files or run directories used, (4) numerical or plotting
behavior changed, (5) tests run and results, (6) remaining risks or deferred
behavior. After a focused post-processing change, stop and report — do not
automatically continue into Abel inversion, full VMI interpretation, new physics
branches, or broad cleanup.
