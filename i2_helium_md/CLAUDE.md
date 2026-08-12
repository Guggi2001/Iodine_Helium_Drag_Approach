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
  (atlas D0; living doc — its GAP markers are the atlas targets; §17 is
  the physical-sensibility ledger: physics / convention / effective /
  scaffolding / missing, with what would retire each)
- `TIER2_MASS_SCENARIOS.md` — **open investigation** (2026-08-11): what
  mass flies through the Coulomb explosion, and why that is the KE₁
  deficit; M1–M6 open questions, M1 gates

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
sub-2.5 — with the lq E_bind-0.048 confound caveat pending. **§6.7 items 1–2 EXECUTED 2026-07-24** (lq sanity battery + E_bind
pair-separation scan; full record in `drag_migration_log_tier2.md` +
findings §6.7): cubic and quadratic are physically distinguishable — the
histogram (W₁) is form-blind but the KE/size/fate observables are not, and
the item-1 lq over-suppression is the FORM (the E_bind confound is
discharged). Nothing moves finc1v725; Tier-0 still rejects lq. **D2b
§4.1 provenance audit EXECUTED 2026-07-26** (read-only; as-built
geometry + which sampler ran in which production — D0 §15, plan §3.1
grid spec amended; parent-document geometry anchor recorded). **GEOMETRY
CORRECTION adopted in principle 2026-07-26** — the droplet geometry is
known wrong vs the experiment's source conditions; Axis A re-scoped as
the correction (atlas plan §3, staged G0–G4 in §3.5); D0 §17
physical-sensibility ledger added. **G0 FROZEN + stage 2b BUILT + Axis A G1
EXECUTED + the §3.5b retained-class arm DELIVERED 2026-07-26/27** — the
geometry grid is **complete at 11/11 cells** (the last five recovered
detection-only, zero MD; six reproduce bit-for-bit). The marginal-exclusion
bracket is **not tight**, so the R ≥ 49 Å deepKE/n̄ carry a stated range
(D0 §14.2, findings §G1.4). **D2b §4.3 grid re-weighting EXECUTED
2026-07-27** (zero MD): pre-registered oracle INADMISSIBLE (support hole
below R = 26.6 Å), post-hoc in-support oracle 7/7 (`nearest` adopted),
corrected-ensemble forecast on record (findings "D2b §4.3", D0 §15.7 —
the landing breaks at ensemble level). **Birth-depth lever decomposed
(D0 §14.3, 2026-07-27):** geometry-locked dressing channel (§4p: the T5
lever tripled supp — load-bearing, saturates at realistic depths) +
parameter-accessible transit channel. **G2 TAKEN 2026-07-27 (user):
the corrected geometry (legacy+raw ⟨N⟩ 12794 / Boltzmann 313.2 K /
E_bind 0.1168) is ADOPTED as the target; retained-policy sub-decision
deferred to G3; `finc1v725` stands until G4.** **G3 Step 1 EXECUTED 2026-07-27**
(zero MD, triggered): twin landmark re-issue at the corrected geometry —
S6 oracle bit-exact, 11-cell twin↔MD authority box measured (supp/n₁_solv
near-quantitative, n̄ residence-scale hot, trap a twin floor, center-pin
trap channel twin-invisible), corrected-ensemble twin row confirms the
D2b forecast cross-instrument (findings "G3 Step 1"; stage `g3landmarks`).
**G3 Step 2 twin scan EXECUTED 2026-07-27** (stage `g3scan`, zero MD,
6480 cells): oracles bit-exact, failure criterion NOT fired — **24
cells gate, all Tier-0-legitimate; fully-unflagged basin (v_c 5.5–6.0,
τ 4.8–6.4, E₀ 0.32–0.37) at the standing well; the standing chord
v_c 7.25 gates nowhere at the corrected geometry; landing
cascade-carried, not selection-carried** (findings "G3 Step 2", D0
§14.5). **G3 Step 3 MD ring EXECUTED 2026-07-27/28** (plan §3.5d; 14 cells ×
N = 500, seed 20260728): **GR-P1..P6 all land — the corrected-geometry
basin is MD-CONFIRMED** (gated: a037 v5.5/τ4.8/E₀0.37 n̄ 4.097; b031
v6.0/τ6.4/E₀0.31; d030 v5.5/τ6.4/E₀0.30; e154 deep-well; c50 shows
the basin nearly reaches v_c 5.0); **finc1v725 MD-measured broken at
realistic droplets** (f725: trap 0.437 incl. 0.159 marginal, n₁ 0);
instrument MD-calibrated (n̄ bias 0.24–2.66 small-end, n₁ ≤ 0.036);
marginal class ≈ 0 at the basin (retained-policy evidence). Adjudicated
(user): `exclude_all_coupled` interim, final call at G4; Axis B folded
into the ring. Findings "G3 Step 3" + D0 §14.5; scorer
`tier2atlas_g3ring_table.py`. **G4 Step 1 EXECUTED 2026-07-28** (plan
§3.5e; 15 × N = 1000): successor candidate `h405` (v_c 5.5 / τ 4.4 /
E₀ 0.405), the KE tension broken, and a measured W₁ floor ≈ 0.67 whose
cause — n₁ and n̄ not simultaneously matchable — localizes the residual
to the mechanism, not the drag surface. **2026-07-29 (low-n KE axis;
full records in the log + findings):** KE₁ re-anchored to the 1.00 eV
peak (user); §3.5g retro-scan — in-surface (v_c, τ, E₀) freedom
EXHAUSTED on KE₁, shallow-birth/mixture lever REJECTED unphysical
(user); §3.5h p_tail ring — placement confirmed but the kill fired
(velocity-only tail softening un-damps the cascade; single-knob p_tail
CLOSED; capped_cubic p_tail guard now [−4, 0]); band-limited tail
rejected as overfitting (user). The **s(n) drag state coupling was
BUILT and probe-EXECUTED 2026-07-29** (S1–S4 delivered, `off` default
bit-identical, full suite green): **all SC predictions REFUTED —
the coupling is measured GATE-CLIPPED** (no ion reaches n ≤ 8 inside
the droplet; ions exit fully dressed at n̄ 19), the axis is STOPPED,
and OQ-F cooling was REJECTED (user). **The §3.5j discussion was HELD
2026-07-29** (records: findings §3.5i.3 + §3.5j, plan §3.5j): the three
verification reads + two placement reads executed zero-MD — the
n = 20–21 spike is falsified ~30× at the reference AND identified as
the `suppressed` fate class (scoring scaffolding); the gas KER is
multi-channel (2.70 eV = point-Coulomb idealization); the experimental
n = 1 KED holds 48.7 % above the (A) cap 1.15 eV; the h405
counterfactual measures the (A) bracket ceiling KE₁ ≈ 0.708. Cold shed
and the CE-instantaneous-stripping variant DECLINED (bulk-refill
argument: net stripping only at the outbound crossing). **Measured
verdict: (A) exit stripping is necessary-not-sufficient; (B)
source-KER spread is required for the KE₁ position; the (C)
combination is the design target.** The **§3.5k budget-slope probe
EXECUTED 2026-07-29** (3 × N = 1000 at the h405 pins, CRN; plan §3.5k,
findings §3.5k, D0 §19): **BP-KILL NOT fired — S_k = 0.389 eV/eV**
(the source-side lever is alive; the pre-registered honest-residual
consequence is not forced); BP-P3 sign-inverted — suppression is
measured **E_int-driven** (∂supp/∂E₀ ≈ +1.7/eV; full proportional (B)
wiring parks fast ions in suppressed-bare with KE₁ at baseline → the
(A)/(B) complementarity is measured from both sides); needle persists
(width requires the mixture). **Domain input recorded (user paper,
RQ7/RQ8 NBs):** gas phase reproduced at 0.8·E_C with Q = 2/Q = 3 from
2.666 Å → channels 2.16 / 4.32 eV per I⁺ (fast peak = I⁺–I²⁺); the
standing 2.70 budget is ~25 % high vs this calibration; Q = 3 /
single-ionization exist as refused scope flags; the I²⁺ discriminator
is ANSWERED (covariance Fig. 6.2: Q3 share ~7 %/~20 % at
1.47/2.94×10¹⁴ W/cm²; retention-twin anchor recorded — per-ion
retention stochastic at fixed kinematics; RQ8 NB). **The §3.5l
pre-design state is CLOSED 2026-07-29**: reads (i)+(iii) executed
zero-MD (bare KED fast-fed + weakly two-lump → slow-bare is a new
kill axis; gas Q3 share 10/17/23 % at 160/300/600 mW; a
calibration-frame rider posted and withdrawn same-session — the Abel
gas export confirms 0.8·E_C in the repo frame), and the 0.8
provenance is CLOSED (Hatherly 1994, J. Phys. B 27 2993 — gas-phase
finite-pulse CE physics, channel-independent fraction, NOT droplet
screening; key numbers extracted to the docs, paper not repo-kept;
RQ7 NB). **Current goal: the (C)
build — `TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md` is the fresh-session
entry point** (design DRAFTED + ADJUDICATED 2026-07-29: channel
mixture {single, Q2 2.16, Q3 4.32} × shared f 0.8 via per-pair scale
emulation + Q3 partner mask, per-channel f_int, depth-graded exit
strip P₀(v)·G(j) retention-twin-anchored, ledger strip term;
OQ-A..K ALL closed — **OQ-K closed 2026-07-29 by user confirmation**
of the 300/600 mW ↔ 1.47/2.94×10¹⁴ mapping (w_Q3 ≈ 0.20
un-provisional); the bare-I⁺ m/q-127-row covariance decomposition is
**postponed (user) until the method proves successful**; checkpoint
v8 granted-scoped to three per-ion fields; weights never fitted on
the n = 1 KED — the circularity guard; dof accounting §6:
8 anchor-frozen / 3 P1-pinned / 1 probe-scanned). **P1 + P2 EXECUTED
2026-07-29 zero-MD** (committed instruments
`tier2atlas_ce_strip_prior.py` / `tier2atlas_ce_gas_widths.py`, all
oracles green incl. exact §3.5j crossing-band + variant reproduction:
strip prior box PINNED a = 2 [1.5, 2.5] / j₀ [1.5, 2] / w_j [0.5, 1],
a = 1 disfavored on W₁-parking + slow-bare; measured channel widths
σ_Q2 0.31 / σ_Q3 0.55 eV supersede the Hatherly priors, E_single
≈ 0.53 — findings "(C) pre-steps P1 + P2"). **P3 EXECUTED
2026-07-29 zero-MD** (`tier2atlas_ce_p3_forecast.py`; weights
(0.30, 0.50, 0.20) entering registration; C-full forecast KE₁ mean
0.92–1.04 / SD 0.51–0.55 / above-1.15 43–50 % vs ref 48.7 % / Q3
carries 53–61 % of the n = 1 bin; two registered tensions —
two-lump mode vs single-mode reference 0.891, and n₁ 0.18–0.23
short of the 0.31 reference; findings "(C) pre-step P3"). **P1–P3
COMPLETE and the PROBE REGISTRATION FROZEN 2026-07-29 (CP-1..8
bands user-approved — authoritative table: design doc §8).
PC-1..5 outcome pre-commitments ADJUDICATED (design §11) and the
probe BUILT + EXECUTED 2026-07-29 behind `[PROCEED TO
IMPLEMENTATION]`** (sampler + exit strip + checkpoint **v8**
(three OQ-E per-ion fields, silent v7→v8 shim) + pair-scale seam +
partner-mask scoring, off-mode byte-identical, full suite green;
4 cells cfull/aonly/bonly/cq3hi at the h405 pins, seed 20260729):
**REGISTRATION FAILED — CP-1..4 FAIL, kills CP-5/6/7 FIRED on
C-full (CP-6 at both f_int,Q3 ends). Measured: (B) alone places
KE₁ 0.920/needle broken but its E_single 0.53 channel traps 0.465;
(A) as specified over-tolls (ε × knock counts p90 20 → 0.5–0.7 eV
vs the 0.2 eV anchor) and FEEDS the suppressed gate (A-only supp
0.55 — Σ(n) collapses under untouched E_int). PC-3 TRIGGERED:
(A)-v1 stops, the current f_int wiring closes, the mixture
survives. Records: findings "(C) probe EXECUTED", D0 §20, design
§8/§9/§11, `atlas_ce_probe.csv`.** **(C) ADJUDICATED CLOSED
2026-07-30 (user): PC-3 ratified, ALL candidate follow-ups
DECLINED, the (C) line SHELVED as a measured boundary — the
covariance-anchored mixture cannot populate n = 1 at the reference
level (Q3 feeder ceiling ≈ 0.11 of the scored ensemble;
suppressed-conversion bound n₁ ≈ 0.22 < 0.31); the positive result
stands (KE₁ position is source physics, not drag/mass-mechanism);
nothing adopted, the scalar budget did NOT retire. Records:
findings "(C) adjudication CLOSED", D0 §20, log 2026-07-30.**
**Current arm (2026-07-30/31): the free-form linear counterfactual —
`TIER2_FREEFORM_LINEAR_TWIN_SWEEP_PLAN.md` is the fresh-session entry
point.** Step 0 (twin KE₁ ranking LICENSED ρ 0.998) + both twin arms
(35-cell basin incl. sub-plateau; ram term rejected, pure linear
selected in-family) + the §6.1 CRN MD ring all EXECUTED: **the ring
LANDS — 6/7 cells gate both-policy-clean, K-KE not fired, SUCCESS at
five cells (best MD KE₁ 0.903 at a = 27.5 vs h405p 0.637; `pure_linear`
form now in `physics/drag.py`; trap axis dead; cost = mid-band
overheating); nothing adopted. The §6.5 h405-clone battery (2026-08-12,
3 × N = 500) then measured that clone a **KE-equivalent only** — KE₁
matches h405 to −0.006 eV and midHot returns to 0.891 (so the mid-band
overheating is an *a*-coordinate, not a property of constant γ), but
W₁ lands 1.079 vs 0.765 and n₁ sits on the gate floor (`gate-marginal`);
twin W₁ is measured non-transferable outside the ring's box. The §6.6
(a, τ) joint ring (2026-08-12, 6 × N = 500) then closed the arm's central
question: **JR-P3 empty — KE₁ and the histogram are NOT simultaneously
improvable here.** The τ 4.8 a-curve is monotone in everything at once,
the family's W₁ minimum (0.785 at a ≈ 35) sits above h405's 0.767, τ 6.4
costs +0.23…+0.41 W₁ at every chord by *evaporating* the shells, and the
one cell beating h405 on KE₁ at in-band midHot (`j1`) pays W₁ 1.145 —
i.e. **KE₁ is a flight-mass observable and the size distribution pins the
flight mass**, so the residual is a mass question (`TIER2_MASS_SCENARIOS.md`),
not a drag one. NEXT USER GATE: routing the clone (accept-marginal /
re-site / recalibrate) + the §6.3 adoption discussion and escalation
battery.** **Atlas §6.5 Step 2 (E_bind) EXECUTED
2026-08-10** (twin, zero MD): `dKE₁/dE_bind = −0.54` — the exit toll is
half the naive ledger value (shared-erf well/density overlap refunds
~45 %), form- and depth-invariant, trap lever form-split, and the axis is
**closed as a KE lever**; **Step 3 MD-CONFIRMED** (1 × N = 500 CRN,
−0.5308, 2 % from the twin; MD-P5 withdrew the Step-2 twin-grading
caveat) — D0 §9.2/§9.3 + findings §6.5 Steps 2–3. **The "45 % refund" is
SOLVED 2026-08-11 as a mass-frame partition — the model is correct, the
ion pays 98.1 % of the well; RQ12's density-width leg is measured NULL
(+0.002 eV) and atlas §6.8 T1–T4 are closed — D0 §9.6.** **Current arm:
`TIER2_MASS_SCENARIOS.md`** — the flight mass, not the drag, carries the
KE₁ deficit (M1 gates). The
G4 adjudications (successor point, retained policy, ledger re-issue)
stay open; D2b A/B remainder, RQ3/RQ5 reads and the margin-3 Å pinned
convention (I88) stay open in-tier. Tier-3 noise
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
