# Tier 2 — Sensitivity Atlas (parameter-influence study program)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger, given per execution stage. Module
> descriptions are *interface contracts*; grids are *specifications*.
>
> **Status: DESIGNED 2026-07-23** (brainstorming session, user-approved
> section by section). No stage has been triggered.
>
> **Parent:** `TIER2_IMPLEMENTATION_PLAN.md` (Tier-2 optimization stage,
> post-F5-reconciliation). **Entry docs:**
> `TIER2_STAIRCASE_PROBE_FINDINGS.md` (§4a–§4ee, I1–I100 — the accumulated
> parameter knowledge this program systematizes);
> `drag_migration_log_tier2.md` (as-built state, standing point);
> `RESEARCH_QUESTIONS.md` (RQ11 NB register — the D4 form axis doubles as
> RQ11's form-discrimination step); `CALIBRATION_MAP.md` (calibration
> classes reused throughout).

---

## 0. Purpose and stance

The Tier-2 build program was *elimination-driven*: land the experimental
observables, then refute candidate owners of the one remaining deficit
(RQ11). This program is *understanding-driven*: build a systematic map of
what each model knob actually does to the observables — a **sensitivity
atlas**, with the accumulated but scattered knowledge of the findings doc
distilled first (D0) and its gaps filled by controlled studies (D2–D4).

**Stance (fixed for the whole program):**

- **Sensitivity atlas — reported, not adjudicated.** No parameter is
  re-tuned by this program. `finc1v725` stays the standing production
  point; nothing here discharges F5 or moves production. If an atlas
  result motivates moving the standing point, that is a *separate,
  pre-registered* decision outside this plan.
- **One-knob-at-a-time (OAT)** around the standing point everywhere,
  except *inside* Axis A, where the droplet-size × position factorial is
  itself the object of study.
- Any drag-law modification remains an interchangeable enum behind its
  own `[PROCEED TO IMPLEMENTATION]` trigger (`DRAG_PORT_DESIGN_DECISIONS.md`
  discipline).
- Dual purpose is declared, not hidden: D4 Step 1 (Method-B form re-fit)
  is simultaneously the RQ11 lever-hierarchy step 1 recorded in the RQ11
  NB register ("TDDFT is the only authority above ~1 Å/ps"). Its RQ11
  reading is reported into `RESEARCH_QUESTIONS.md`; its atlas reading
  lands here. One study, two write-ups.

**Center point (standing, from the log — verify against
`drag_migration_log_tier2.md` before execution):** `capped_cubic`
v_c 7.25 Å/ps / p_tail = −1, τ = 3.2 ps, E₀ = 0.27 eV, `rq4graded`
ladder, Landau v_limit 0.58 Å/ps, locked Tier-0 b, jointly-extracted
E_bind = `binding_energy_I_ion_eV` = 0.1168 eV, biphasic mass mechanism,
production 2.70 eV, 30 ps ion window, N = 1000 per battery member.

## 1. Common measurement protocol (governs every axis)

### 1.1 Observable vector

Every atlas cell reports the **committed scorer's** vector, unchanged —
no new metrics without pre-registration:

- fate split: suppressed weight, trapped fraction, bare weight,
- n̄ (detected mean size), n₁_solv (n = 1 solvated share),
- midHot (n = 2–3 KE ratio), W₁_solv (solvated-histogram Wasserstein),
- n₁ KE ratio (median-anchored, I77 convention),
- deep-bin KE ratios sim/ref over n = 10–17 (the RQ11 axis),
- χ²_med (median-anchored KE chi-square).

Per-axis figures are additive (Axis A: full detected histogram + the
n = 1 velocity distribution per cell), produced by the existing figure
surface (`plot_detection_summary.py` sections), never by new ad hoc
histogram code.

### 1.2 Significance yardstick

Influence is reported in **seed-SD units**, from the pooled N = 5000
battery (findings §4cc): W₁_solv SD ≈ 0.04, deep-bin KE read
0.0603 ± 0.0033, midHot seed-robust band [0.84, 1.14]. Convention:

- |Δ| < 1 seed-SD → "not resolved at cell N",
- 1–2 seed-SD → "suggestive" (eligible for N = 1000 confirmation),
- > 2 seed-SD → "resolved influence" (atlas headline material).

A cell run at N = 500 has ≈ √2 larger scatter than a battery member;
the doc's tables state cell N next to every number.

### 1.3 Cost ladder

1. **Twin counterfactual first** wherever the knob is twin-expressible
   (drag-form S(v) via the §4ee γ-tail seam, ladder scans, ε forward
   models). Oracle discipline: every twin session must reproduce the
   committed `h2b_s6_final` v7.25 row bit-for-bit before any new cell.
2. **Small-N MD** (N = 500, one fixed seed, full pipeline
   neutral → ion → relaxation → detection) where the twin cannot express
   the knob (droplet geometry, E_bind, mechanism-feedback legs of
   E₀/τ). Twin authority limits (frozen chord, channel (d), I69/I73)
   are restated in every twin-derived table.
3. **N = 1000 confirmation** only for effects > ~2 seed-SD, and MD
   spot-checks for drag forms only after D4 Steps 1–2 filter the list.

### 1.4 Discipline

- Scorer reuse mandatory (the F3/§4cc committed scorer; oracle rows
  re-verified at the start of every scoring session).
- All new code — generator scripts, report scripts, any drag-form enum —
  behind `[PROCEED TO IMPLEMENTATION]`, per stage.
- Run-directory namespace: atlas runs use a distinct tag (working name
  `_tier2atlas_`), following the probe-tag precedent (no literal
  `_tier2_` substring; lock test at implementation time; F3
  explicit-exclusion fallback).
- Findings are recorded **parameter→influence-style from day one** in
  the atlas findings doc (§2.3) — chronology goes to the migration log
  only. This is the anti-bloat rule: the 5000-line probe findings doc is
  the counterexample this program is structured to avoid.

## 2. Deliverable documents

### 2.1 `TIER2_PARAMETER_INFLUENCE.md` — compact reference (D0, written first)

Pure distillation of `TIER2_STAIRCASE_PROBE_FINDINGS.md` — **no new
analysis, zero MD**. The findings doc stays untouched as archive. One
section per knob, table-first:

| element | content |
|---|---|
| role | one sentence: where the knob enters the mechanism |
| class | Sourced / Bounded / Free / Derived (CALIBRATION_MAP convention) |
| influence | per observable: direction + magnitude, with the measured numbers (e.g. "ladder tail ±5 % → ∓3.5 pts suppressed weight, ~+5 % deep-KE, I93") |
| couplings | structural couplings (Σ-coupling I93; drag↔E_bind joint extraction §6.5.1; v_c↔τ↔E₀ arbitration basin §4w) |
| status | locked / standing / open, with § and I-number pointers |

Knob list (initial; the doc may split or merge): drag form + b, v_c /
p_tail, τ, E₀, κ, picture, f_int, f_ret, λ₀, ladder D₀(n) shape + tail +
Σ, Landau v_limit, per-shed ε, E_bind, droplet geometry (size ×
position), sampling laws (D2b chapter). Knobs with no measured influence
yet (droplet geometry, E_bind, sampling laws) carry an explicit **GAP**
marker — the GAP list is the atlas's target list.

Living doc: each atlas stage merges its results back in, replacing its
GAP marker.

### 2.2 This plan

`TIER2_SENSITIVITY_ATLAS_PLAN.md` — the design (this document).

### 2.3 `TIER2_SENSITIVITY_ATLAS_FINDINGS.md` — results

Created at first execution. Structured by axis, parameter→influence
tables + the per-axis figures, seed-SD-unit reporting, twin-authority
boxes where applicable. Delivery/decision records go to
`drag_migration_log_tier2.md` as usual.

## 3. Axis A — droplet size × initial position (D2; MD, most novel)

Everything so far sits at the production-sampled geometry; size and
position have never been controlled independently. The size–position
pair sets the He path length and hence the total drag/cooling exposure —
plausibly the largest unmapped influence on the histogram and on the
n = 1 velocity distribution.

### 3.1 Grid

**3 × 3 factorial, N = 500/cell, one fixed seed, full pipeline:**

- **R (droplet radius):** the production post-pickup size distribution's
  **q10, mean, q90**, converted N → R by the existing 0.8·ρ_bulk
  convention. The three concrete Å values are computed from the
  production sampler at execution (one scripted draw, recorded with the
  sampler's seed in the findings doc) — pinned by procedure here so the
  grid tracks the actual production distribution, not a guess.
- **r/R (fractional initial radial offset of the I₂ center):**
  {0 (center), 0.5, 0.8 (sub-surface)}, overriding the Boltzmann
  positional sampler with fixed offsets for grid cells (direction of the
  offset randomized per molecule as in production; only the radius is
  pinned).
- All other knobs at the standing point; droplet-size *sampling* is
  disabled per cell (fixed R), which is exactly the controlled variable.

### 3.2 Machinery

Sampling machinery exists (`sampling/droplet_sizes.py`,
`sampling/radial_positions.py`, `single_pulse_droplet_distribution`
preset). New pieces (behind the trigger): a generator mirroring
`gen_tier2_runs.py` (atlas namespace, overwrite/resume guards, budget
guard) and a pure-scorer report script (F3 idiom: knobs from `cfg.json`,
text table + CSV, figures behind `SHOW_FIGURE` only). Whether fixed-R /
fixed-offset cells need small additive config surface (e.g. a
fixed-offset override knob) is settled at implementation planning; any
such field follows the rule-2 carry convention if staged.

### 3.3 Questions the axis answers

- Does exposure alone (bigger droplet, deeper start) reorganize the
  detected histogram, and along which observables first?
- Is the n = 1 velocity distribution's shape a **geometry fingerprint**
  (distinguishable signatures of R vs r/R), or degenerate along an
  exposure-equivalent diagonal?
- Does geometry move the deep-bin cold tail (RQ11-relevant; reported
  under atlas stance, not adjudicated)?
- How much of the production ensemble's spread is geometry-inherited vs
  mechanism-stochastic?

### 3.4 Optional follow-up (not in the grid)

A binned read of the existing pooled N = 5000 battery by its sampled
(R, r) draws — free but statistically entangled; only meaningful as a
cross-check after the grid fixes the interpretation.

## 4. Axis A addendum — D2b: provenance audit + salvageability of the legacy sampling laws

The production sampling laws are theory-laden legacy ports. This
sub-axis audits them and measures what they change.

### 4.1 Provenance audit (zero MD, doc work)

Enumerate every theoretical ingredient with source + calibration class,
as the "sampling laws" chapter of `TIER2_PARAMETER_INFLUENCE.md`:

- **Size law** (`droplet_sizes.py`): nozzle correlation
  ⟨N⟩ = k₁·p^0.97·T^−3.88·d² (Lackner, empirical); log-normal width
  δ = 0.625 (Kornilov 2009); the pickup + evaporation chain (mean free
  path 4 Å, ε = 0.04 relative loss per collision, **E_solv = 14 meV
  with the documented thesis-text-vs-figure discrepancy against
  30 meV**); droplet density 0.8·ρ_bulk; the analytic D4 prior family
  (truncated Kornilov ln-normal, pickup-weighted ×N^(2/3),
  `cfg.droplet_size_prior`).
- **Position law** (`radial_positions.py`): thermal-equilibrium
  assumption p(r) ∝ r²·exp(−U_drop(r−R)/k_B T) at `T_particles_K` in
  the mean-field molecule–droplet well (`binding_energy_molecule_meV`,
  steepness) — an equilibrium ansatz, not a measured dopant location.

### 4.2 Distribution-level A/B (zero MD)

Compare the *sampled distributions themselves* under variants: raw vs
post_pickup vs analytic prior; E_solv 14 vs 30 meV; δ perturbations;
position law vs simple alternatives (center-pinned, uniform-in-volume).
Pure sampling plots; no propagation.

### 4.3 Observable-level sensitivity — via grid re-weighting (zero new MD), then ≤ 2 confirmations

The §3 grid is the transfer function {(R, r/R) → observable vector}.
Any candidate sampling law is evaluated by **re-weighting the grid
cells with its (R, r) density** — predicting the ensemble effect of
each theory variant without new MD (the coarse 3 × 3 support is a
stated limitation; interpolation convention recorded in the findings
doc). Only variants that matter after re-weighting get a full-sampler
MD confirmation (N = 500 each; expected: E_solv 14→30, and legacy
post_pickup vs analytic prior — ≤ 2 runs).

### 4.4 Salvageability verdict (defined, not vibes)

The legacy sampling laws are **salvageable** if the landed observables
stay inside the §1.2 seed-SD yardstick under all defensible variants
(class Bounded moves within their bounds; the E_solv discrepancy tested
at both published values). A variant that breaks the landing is flagged
as a real model risk in the atlas findings and in
`TIER2_PARAMETER_INFLUENCE.md` — reported, never silently retuned.

## 5. Axis B — E₀ and τ (D3; MD + partial twin)

The arbitration *sampled* these knobs (τ scans §4-series; the (v_c, τ,
E₀) joint-landing basin §4w); it never produced clean OAT influence
curves around the standing point.

- **Grids (OAT, N = 500/cell, fixed seed):**
  E₀ ∈ {0.17, 0.22, **0.27**, 0.32, 0.37} eV;
  τ ∈ {1.6, 2.4, **3.2**, 4.8, 6.4} ps (log-spaced; standing values
  bolded — center cells reuse existing battery data, so 4 + 4 new MD
  cells).
- **Twin-expressibility marking:** τ enters the twin through the cooling
  exposure K, E₀ through the gate — each knob's twin-expressible leg is
  run twin-first (cheap curve at dense resolution), the
  mechanism-feedback remainder (pickup/evaporation coupling, E_int
  reservoir dynamics) is what the MD cells measure. Twin-vs-MD
  disagreement is itself an atlas result (it localizes where mechanism
  feedback matters).
- **Deliverable:** per-knob influence curves (each observable vs knob
  value) with seed-SD bands; couplings noted against the §4w basin
  (E₀ and τ are not independent in arbitration space — the OAT curves
  quantify the local slopes that basin only implied).

## 6. Axis C — drag-law forms + E_bind (D4)

### 6.1 Form list and dimensional analysis

Friction convention unchanged (`γ(v) = |F|/v` [amu/ps], force `γ·v`, no
leading m; mass-agnostic module). ρ̂ = normalized local He density
(dimensionless). Candidate in-band forms:

| form | F(v) | γ(v) [amu/ps] | coefficient units | status tag / Tier-0 prior |
|---|---|---|---|---|
| pure cubic (baseline) | ρ̂·b·v³ | ρ̂·b·v² | b [amu·ps/Å²] | locked Tier-0 (§9); re-fit must reproduce locked b = 2.5154 (oracle) |
| power law (free n) | ρ̂·C·vⁿ | ρ̂·C·vⁿ⁻¹ | C [amu·Å^(1−n)·ps^(n−2)] | Tier-0 discriminator: shared n̂ = 2.927 ≈ 3, equivalent to cubic (§10) |
| quadratic | ρ̂·b₂·v² | ρ̂·b₂·v | b₂ [amu/Å] | Tier-0 **rejected-by-objective** (the lq fits collapse to this corner, a → 0; Δobj +0.034) |
| linear + quadratic | ρ̂·(b₁v + b₂v²) | ρ̂·(b₁ + b₂v) | b₁ [amu/ps], b₂ [amu/Å] | Tier-0 **rejected**: shared Δobj +0.034 (outside ±0.005 band), held-out 9 Å prediction FAILS (0.699 vs ≤ 0.45) |
| linear + cubic | ρ̂·(b₁v + b₃v³) | ρ̂·(b₁ + b₃v²) | b₁ [amu/ps], b₃ [amu·ps/Å²] | Tier-0 enum exists; shared fit drives a → 0 (= pure cubic) |

The **shifted cubic** F = ρ̂·b′·(v−v_L)³ (Free/de-anchored per I100) is
**dropped from Step 1 by user decision (2026-07-23)** — it is the only
form with no Tier-0 artifact, so including it would be the step's only
new fitting work. It stays **parked, re-addable** if the trace-tail
inspection shows genuine threshold structure; the machinery needs no
change.

**Parked candidates (added 2026-07-23; entry gate = the zero-cost γ(v)
arithmetic, fit only if that says they matter).** Two forms express
the S-shape the standing evidence sketches, as minimal extensions of
existing machinery:

| form | F(v) | γ(v) [amu/ps] | parameters | why promising |
|---|---|---|---|---|
| saturating cubic (Padé) | ρ̂·b·v³ / (1 + (v/v_s)³) | ρ̂·b·v² / (1 + (v/v_s)³) | b [amu·ps/Å²], v_s [Å/ps] | reduces to cubic in-band; high-v asymptote = constant force ρ̂·b·v_s³ — exactly the arbitrated p_tail = −1 behavior, smooth, with **one Free knob replacing (v_c, p_tail)**; the 9 Å band top (4.95) may let Method-B constrain v_s → tail Derived |
| subtractive (gated) cubic | ρ̂·b·v·(v² − v_f²) for v > v_f, else 0 | ρ̂·(b·v² − a)₊ | a [amu/ps], b [amu·ps/Å²]; v_f = √(a/b) [Å/ps] | the **existing `linear_cubic` family with a < 0** (sign-bound relaxation + zero-clamp, no new enum); at v_f = 1.5 it *is* the §4ee twin landing window (×0.44 at v = 2, ×0.91 at v = 5, → 1 high-v); suppression ×0.77–0.91 persists inside the 9 Å band, so the traces can genuinely constrain v_f |

Their γ(v) columns are computed first with **implied parameters from
existing results** (v_s ≙ the arbitrated 7.25; v_f from the twin
window 1.5–2.0) at locked b — zero cost, no fit. The quadratic family
is closed (Tier-0 rejection + wrong-direction low-v limb); memory /
Basset forces (unconstrainable by current traces) and density-exponent
variants (the gate axis, separately Derived) are explicitly not
candidates.

**Tier-3 bridge note (recorded here, not a Step-1 form):**
**discrete-emission drag** — Poisson momentum-loss events (roton-pair
sawtooth per NB-RQ11-9; or ps-scale vortex-ring shedding per
NB-RQ11-7, inside the MD window) with the *same mean* as the standing
law. Invisible to Method-B (traces constrain the mean) but adds
fluctuations — and the program's width-type residuals (RQ3 σ 6–10×
too narrow, under-dispersed cov/VMI panels, deep-bin scatter) are
second-moment observables. This is the physics-motivated starting
candidate for the Tier-3 noise-channel design, superseding generic
FDT noise as the default ansatz.

All units balance against v [Å/ps]; each form composes with the
production capped tail (v_c, p_tail) unchanged unless a form makes the
cap ill-posed, which the re-fit report must state.

**Anti-circularity note (user-raised, 2026-07-23).** The downstream
system (v_c, τ, E₀, ladder, E_bind pairing) was arbitrated
*conditional on cubic*, and those knobs trade against the drag through
K (I47/I72) — so the Tier-2 landing is **weak evidence about the
form** (a quadratic-based system might also have roughly landed after
re-arbitration) and must never be cited as confirming cubic. The
form's authority comes solely from the Tier-0 trace instruments,
which are downstream-independent: the free-exponent fit picked
n̂ ≈ 3 unforced, the shared objective rejects lq against TDDFT itself,
and the held-out 18 Å→9 Å prediction (the designed anti-circularity
test) fails lq outright. This privilege is **band-limited** —
outside 2.54–4.95 Å/ps the form was fit to experiment (the Free cap),
which is the standing critique the collaborator ask addresses.

**Registered priors and the framing finding.** The Tier-0 verdicts
(`TIER0_FINDINGS.md` §9/§10 tables) stand as **pre-registered priors**
for Step 1: a form rejected there stays rejected unless the re-fit
overturns it under the *same* objective with a stated reason (new data,
extended band, or a demonstrated objective defect) — never by switching
metrics until it passes. Motivating synthesis for the whole axis:
I72 (three observables prefer three v_c) + I40 (high-v over-drag) +
RQ11 (low-v over-drag) jointly describe an **S-shaped deviation from
cubic at both ends** — steeper-than-cubic (or gated) below ~2 Å/ps,
cubic in the TDDFT-pinned mid-band, sub-cubic/saturating above ~7 —
which no single polynomial or power law expresses. The capped tail is
the crude top-end half of that shape; the §4ee sub-form window is the
bottom-end half. The principled cap-remover is therefore **extending
TDDFT authority to where the cap lives** (production-kinematics traces
— the collaborator ask), not a different polynomial: with the fitted
Tier-0 coefficients, linear+quadratic softens the top end but carries
~2.5× the cubic's drag at v = 2 Å/ps (γ ≈ 26 vs ≈ 10 amu/ps) — the
wrong direction for both open deficits at low v — and even pure
quadratic still grows ×2.1 between v_c 7.25 and the production peak
10.5 Å/ps, where the arbitration wanted saturation (p_tail = −1).

### 6.2 Step 1 — Method-B form table (zero MD; dual-purpose; mostly reuse, not rerun)

Compile the **full Tier-0 three-mode protocol per form** against the
existing 9/18 Å TDDFT traces (`METHOD_B_trajectory_matching_extraction.md`).
**Existing Tier-0 fit artifacts are reused verbatim, not re-run**
(`data/reference/drag/shared/trajectory_matching/` + the per-case
`<variant>/fit_parameters.json` files already cover pure cubic, power
law, both lq variants, and linear+cubic in all three modes). With the
shifted cubic dropped (§6.1), **no new fits run at all** — only **one
oracle re-run** (shared pure-cubic must reproduce locked b = 2.5154)
to verify the machinery. Re-running everything is required only if the
objective or fit window changes — an explicit, logged decision, never
a default. The three modes, and what each answers:

1. **Shared joint fit** (both traces, one coefficient set + E_bind) —
   the physics claim (one law + ρ̂ scaling), and the *gating* mode: the
   Tier-0 objective (equal-weight mean of the per-case in-window |v₂|
   RMSEs + escape penalty) with the §9.4 bands and the ±0.005
   equivalence band, unchanged.
2. **Held-out** (fit 18 Å alone, predict the untouched 9 Å) — the
   trap-catcher for trajectory-matched-but-wrong laws (the Finding-3
   lesson: in-window match ≠ correct production behavior).
3. **Per-case single-curve** — **diagnostic only, never
   preset-wired**: measures identifiability (Tier-0 found it
   case-asymmetric — 9 Å alone pins n̂ = 2.651 ± 0.026, 18 Å alone
   rails to n = 4) and transferability (a large per-case-vs-shared
   coefficient gap = the shared-form ansatz failing for that form —
   itself a discrimination signal). Caveat carried: the 9 Å curve is
   transverse-contaminated (Tier-0 Finding 2), so per-case 9 Å
   coefficients absorb unrepresentable drift.

Plus the **low-v re-inspection of the trace tails** (do they show
coast/threshold behavior?), and one zero-cost diagnostic column: for
every re-fitted form, tabulate its implied γ(v) at **v = 2 / 7.25 /
10.5 Å/ps** against the arbitrated capped-cubic — so "would this form
still need a cap, and which way does its low-v limb bend relative to
the RQ11 deficit" is a computed column, not an argument.

Outputs: per-form × per-mode coefficients + residuals, the γ(v)
comparison table, the tail-inspection read. **This is simultaneously
RQ11 lever-hierarchy step 1** (form discrimination — TDDFT as the only
in-band authority, I100); the RQ11 reading goes to
`RESEARCH_QUESTIONS.md`. Stated limitation: the calibrated bands only
partially overlap the RQ11 deficit window (v ≈ 2–7 Å/ps), so low-v
discrimination power may be limited — which the trace-tail
re-inspection and the standing collaborator ask (TDDFT tails / extended
traces / production-kinematics traces) are there to supplement.

### 6.3 Step 2 — twin sweep per re-fitted form (zero MD)

Each Step-1-calibrated form runs through the §4ee S(v)-seam twin
machinery at the standing cell → observable-vector influence map.
Oracle discipline per §1.3; twin authority limits restated (the twin
under-expresses late low-v dissipation — channel (d) — so deep-bin
magnitudes are direction-only here).

### 6.4 Step 3 — MD spot-checks (gated)

Only for forms still interesting after Steps 1–2. Requires the drag-form
enum build (new `SimConfig` enum members + `physics/drag.py` branches),
behind its own `[PROCEED TO IMPLEMENTATION]`; N = 500 each.

### 6.5 E_bind influence study (zero-MD first, MD confirmation second)

**Learning goal (user-stressed, 2026-07-23):** does the
jointly-extracted E_bind actually influence the landed histogram and
the KE shape — i.e. how much of the landing is conditional on the
0.1168 eV member of the Tier-0 coupled pair?

**Scan values — the Tier-0 extracted spread, not arbitrary
multipliers.** Every Tier-0 form fit co-extracted its own E_bind; the
on-disk spread is the physically-motivated grid:
**{0.048 (lq shared), 0.071 (18 Å per-case cubic), 0.113 (power law),
0.1168 (standing), 0.154 (9 Å per-case cubic)}** eV. (The earlier
×0.5–×1.5 bracket is superseded; it spans the same range without the
provenance.)

**Step order (cost ladder §1.3 applied):**

1. **Scoring-level swap on the existing pooled battery (zero MD,
   §4dd forward-model idiom).** E_bind is the droplet exit well the
   ion climbs during propagation, so a literal re-run-free swap is
   impossible — but the first-order counterfactual is computable on
   stored states: escapers shift final KE by −δ (δ = E_bind_new −
   0.1168); fragments whose exit-side radial energy falls between the
   two barriers flip ejected↔trapped and are re-classed. Outputs the
   full observable vector per candidate E_bind. Stated limitation:
   weight-level — no feedback into cascade timing or near-surface
   drag exposure (second order).
2. **Twin counterfactual** at the standing cell (the S6 twin models
   the barrier — the §4ee un-trapping arm — and re-integrates chords,
   capturing the deceleration-under-drag coupling the scoring swap
   linearizes away).
3. **MD confirmation, N = 500** — only for values whose zero-MD effect
   exceeds ~1 seed-SD (expected ≤ 2–3 cells, not the full grid). Runs
   deliberately under `allow_unvalidated_binding_pairing` (the §6.5.1
   guard exists because drag ↔ E_bind is a jointly-calibrated pair;
   the flag's use is the documented, intentional exception, mirroring
   the §6.6 mass-pairing precedent).

Focus observables: trapped fraction, suppressed class, exit
deceleration → deep-bin KE (E_bind is the direct knob on the
trapped/exit boundary — the §4ee un-trapping echo makes this the
cleanest probe of the candidate-(iv)-adjacent margin class). Caveat
stated in every table: these are sensitivity reads of a broken joint
calibration, not candidate points; E_bind is *not* connected to the
He ladder D₀(n)/Σ (different bookkeeping surfaces: droplet mean-field
exit well vs per-atom shell energetics).

### 6.6 Parked test — quadratic counterfactual arbitration (twin, zero MD)

**Purpose (user-adjudicated 2026-07-23):** measure, rather than argue,
the anti-circularity question — could the *quadratic* form have landed
the Tier-2 observables if the downstream knobs had been re-arbitrated
around it? The Tier-0 trace verdict rejected lq with its own best
coefficients (shared Δobj +0.034; held-out 9 Å prediction 0.699 FAIL),
but no quadratic-based *downstream* arbitration has ever run — the
compensation capacity of (cap, τ, E₀) through K (I47/I72) makes that
an unknown, not a refuted claim.

**Design:** the lq form with its **own Tier-0 artifacts** — shared
coefficients (a ≈ 0, c ≈ 12.80) *and* its co-extracted
E_bind = 0.048 eV — run through the reconstructed S6 twin (§4dd/§4ee
machinery, oracle bit-exact first) with a **(cap, τ, E₀)
re-arbitration grid** around the standing cell; standard observable
vector, seed-SD yardstick. High-v treatment is swept, not inherited
(lq needs its own cap decision; its γ still grows linearly above the
band).

**Interpretation, pre-registered:** (a) if no lq cell approaches the
cubic landing under full downstream freedom → "quadratic fails
downstream too" becomes a *measured* statement and the anti-circularity
note closes; (b) if some lq cell lands comparably → the Tier-2
observables cannot discriminate forms, and the form choice rests
solely and explicitly on the Tier-0 trace instruments (held-out above
all). Both outcomes are informative; neither moves the standing point
(atlas stance). Twin authority limits (I69/I73, channel (d)) carried
as always — a near-landing in the twin would need an MD spot-check
before any strong claim.

**Conditional follow-up (only if (b) or a near-(b) fires):** Method-B
form discrimination re-run under the Tier-1a anchored variable-mass
m(t) instead of fixed m_eff — removes the largest Tier-0 architecture
conditionality (the 9 Å case sheds ~7 He in-window; estimated
exponent bias from the ±7 % mass drift is ~0.2 in n — real, but an
order of magnitude short of moving 3 → 2; the re-run turns that
estimate into a measurement).

**Status: EXECUTED 2026-07-24 (user-triggered after the Step-1 table) —
outcome (b) measured** (initial coarse-grid (a) read overturned the
same day by the user-caught fine (τ, E₀) rescan — the pass window sat
between the coarse τ points). Oracle bit-exact; 11 v_c chords; fine
scan 4,959 cells → **six four-way passes** at v_c 8.8–9.0 / τ 3.3–3.5 /
E₀ 0.27 (best cell beats the cubic base on W₁ and midHot centering).
Per pre-registration: **the Tier-2 observables cannot discriminate the
forms**; the form choice rests solely on the Tier-0 trace instruments
(held-out above all), and the anti-circularity note is upheld and
sharpened — the Tier-2 landing is measured to be form-blind, never
evidence for cubic. The lq landing is a needle (τ window ~0.3 ps,
single E₀ grid point) vs cubic's §4w basin — twin-geometry statement,
MD-unconfirmed. Pre-registered consequences now open (user decisions):
MD spot-check of an lq-passing cell (enum build + trigger) and the
variable-mass Method-B re-run. **MD spot-check EXECUTED 2026-07-24**
(3 × N = 500, `capped_linear_quadratic` enum + atlas-namespace
generator): all three cells land — including the off-needle control —
so (b) is **MD-confirmed** and the twin needle was a frozen-chord
artifact (the MD basin is wide). The variable-mass Method-B re-run
remains the open follow-up. Full record:
`TIER2_SENSITIVITY_ATLAS_FINDINGS.md` §6.6 section.

## 6.7 — §6.6 follow-on program (designed 2026-07-24, user-adjudicated; **items 1–2 EXECUTED 2026-07-24**, item 3 pending)

> **Item-1 (lq sanity battery): EXECUTED 2026-07-24.** 5 × N = 1000
> paired-by-seed; pre-registration **3 MET / 2 MISSED**. Histogram
> form-blindness confirmed at battery scale (W₁, midHot); the §6.6 KE
> "lands better" **refuted** (paired lq χ²_med worse on 5/5 seeds); lq
> over-suppresses.
>
> **Item-2 (E_bind pair-separation scan): EXECUTED 2026-07-24.** Initial
> N = 500 single-seed + firm-up N = 1000 × 3 paired seeds. **Resolves the
> confound: the over-suppression is the FORM, not the well** (matched-well
> Δsupp +0.034, ~8σ); item-1's histogram match is a **(form, well)
> co-compensation** — at a matched well lq is colder (midHot Δ−0.070) and
> smaller (n̄ Δ−0.57); **W₁ alone stays form-blind even at matched well**.
> The single-seed W₁ reads were seed-noise; χ²_med inconclusive even at
> N = 1000; deep well 0.154 over-retains (handover-guard trips). Full
> record: `TIER2_SENSITIVITY_ATLAS_FINDINGS.md` §6.7 item-2 RESULTS +
> `drag_migration_log_tier2.md`.
>
> **Item 3 (leave D → Axis A geometry grid → D2b → Axis B) stays behind
> its trigger.**

Agreed next steps after the §6.6 MD spot-check (all three lq cells
landed; findings §6.6 MD subsection). Each item stays behind
`[PROCEED TO IMPLEMENTATION]`.

1. **lq sanity battery — pooled N = 5000.** Base cell **qcc**
   (lq v_c 8.8 / τ 3.4 / E₀ 0.27; chosen on MD merit: W₁ 0.496,
   midHot 1.001). 5 × N = 1000 with **the cubic battery's exact seeds
   20260722–20260726** (paired-by-seed design: every observable
   becomes a per-seed paired difference vs the corresponding
   `bigc1v725s*` member — seed noise cancels; the form-blindness
   statement gets a paired-SD error bar). Free bonus: pooled deep bins
   (~300 fragments at n = 10) are the N = 1000-scale escalation of the
   NB-RQ11-12 mid-band read. Light BN-style pre-registration frozen
   before launch: pooled W₁ ∈ ~[0.46, 0.57]; midHot in the seed-robust
   band; deep slope persists; per-member χ²_med below the paired cubic
   member; supp within cell scatter of 0.187. Namespace
   `tier2atlas_conf270_qccbigs{1..5}`.
2. **E_bind scan on qcc — the pair-separation experiment (run first or
   alongside).** Two N = 500 runs, qcc config, seed 20260721, only
   E_bind overridden: **0.1168** (the standing cubic-pair value) and
   **0.154** (top of the Tier-0 extracted spread) — §6.5's
   physically-motivated grid. These deliberately break the §6.5.1
   joint pairing → run under the unvalidated-binding escape hatch (the
   §6.5 pre-registered exception), stamped per run. Headline question:
   **does the lq KE advantage (χ²_med 62–83 vs cubic 125.7) survive at
   E_bind 0.1168?** Collapse → the "lands better" finding re-attributes
   to the shallow exit well; survival → the mid-band-γ attribution
   firms up. Replaces the D0 §1/§2 confound caveats with a measurement.
   Focus reads per §6.5: trapped, suppressed, deep-bin KE, χ²_med.
3. **Then leave D** (parked: subtractive re-fit decision, variable-mass
   Method-B re-run, collaborator ask) and proceed per the §7 order:
   **Axis A geometry grid (2b)** → D2b audit → Axis B E₀/τ curves
   (whose center cells then reuse both batteries).

Cost: 5 × 1000 + 2 × 500 = 6,000 fragments ≈ 1.2× the cubic battery —
the atlas's largest single spend; hence the pre-registration and
seed-pairing discipline. Atlas stance unchanged: nothing here moves
`finc1v725`; adoption of any form/E_bind value stays outside this
program.

## 7. Execution order and budget

| stage | content | MD cost | gate |
|---|---|---|---|
| 1 | D0 compact reference doc | zero | none (doc work) |
| 2a | D4 Step 1 Method-B form table (artifact reuse + oracle; no new fits) + trace-tail re-inspection | zero | trigger (analysis scripts) |
| 2b | Axis A grid (§3) — runs parallel to 2a | 9 × 500 | trigger (gen + report scripts) |
| 3 | D2b audit + A/B + re-weighting (+ ≤ 2 confirmations) | ≤ 2 × 500 | trigger (sampling-variant runs) |
| 4 | Axis B E₀/τ curves | 8 × 500 | trigger |
| 5 | D4 Step 2 twin sweeps | zero | none (scratchpad twin, §4ee precedent) |
| 6 | D4 Step 3 spot-checks + E_bind zero-MD swaps/twin, then conditional MD confirms | ~0–7 × 500 | trigger (incl. enum build; swaps are scratchpad) |
| 7 | synthesis: merge all results into `TIER2_PARAMETER_INFLUENCE.md`, close GAP markers | zero | none |

Total new MD ≈ 19–26 cells × N = 500 ≈ 9.5–13k fragments ≈ 1.9–2.6×
the pooled battery — the explicit spend of the program (E_bind MD
cells now conditional on the zero-MD swap/twin reads). N = 1000
confirmations (only for > 2 seed-SD effects) are extra and
case-by-case.

## 8. Risks and boundaries

- **Twin authority:** frozen-chord limits (I69/I73, channel (d)) —
  every twin-derived number carries its authority box; closure-magnitude
  questions belong to MD.
- **Grid decoupling vs production coupling:** production samples (R, r)
  jointly; the factorial deliberately decouples them. The D2b
  re-weighting is the bridge back; its 3 × 3 support is coarse and the
  interpolation convention must be stated, not implied.
- **Guard flags:** E_bind cells run under
  `allow_unvalidated_binding_pairing`; biphasic cells continue under
  `allow_inconsistent_mass_pairing` (§6.6). Both usages are deliberate
  and logged per run.
- **No adoption path inside this program:** a drag form, sampling law,
  or E_bind value that "looks better" against experiment is an atlas
  finding. Adoption requires Method-B/TDDFT evidence where available
  (forms) and a separate pre-registered decision (everything else).
- **Anti-bloat:** findings recorded parameter→influence-style only
  (§1.4); if a session produces chronology, it goes to the migration
  log, not the findings doc.
- **Scorer drift:** oracle rows re-verified at the start of every twin
  and scoring session (bit-exactness precedent §4dd/§4ee).

## 9. Cross-links

- `TIER2_STAIRCASE_PROBE_FINDINGS.md` — source material for D0
  (§4a–§4ee, I1–I100); stays archival.
- `RESEARCH_QUESTIONS.md` — RQ11 NB register (D4 Step 1's second
  write-up target); RQ3/RQ5 unaffected by this program.
- `CALIBRATION_MAP.md` — calibration classes; the D0 doc and D2b audit
  feed newly-classed entries back to it.
- `METHOD_B_trajectory_matching_extraction.md` — the D4 Step 1
  machinery.
- `DRAG_PORT_DESIGN_DECISIONS.md` — enum-interchangeability rule for any
  Step-3 build.
- `drag_migration_log_tier2.md` — delivery/decision records per stage.
