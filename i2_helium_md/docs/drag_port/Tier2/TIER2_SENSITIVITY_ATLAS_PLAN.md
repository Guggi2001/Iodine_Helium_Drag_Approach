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

| form | F(v) | γ(v) [amu/ps] | coefficient units | status tag |
|---|---|---|---|---|
| pure cubic (baseline) | ρ̂·b·v³ | ρ̂·b·v² | b [amu·ps/Å²] | locked Tier-0; re-fit must reproduce locked b (oracle) |
| quadratic | ρ̂·b₂·v² | ρ̂·b₂·v | b₂ [amu/Å] | new |
| linear + quadratic | ρ̂·(b₁v + b₂v²) | ρ̂·(b₁ + b₂v) | b₁ [amu/ps], b₂ [amu/Å] | new |
| linear + cubic | ρ̂·(b₁v + b₃v³) | ρ̂·(b₁ + b₃v²) | b₁ [amu/ps], b₃ [amu·ps/Å²] | Tier-0 enum exists |
| shifted cubic | ρ̂·b′·(v−v_L)³ | ρ̂·b′·(v−v_L)³/v | b′ [amu·ps/Å²] | **Free / de-anchored (I100)** — phenomenological, kept for completeness |

All units balance against v [Å/ps]; each form composes with the
production capped tail (v_c, p_tail) unchanged unless a form makes the
cap ill-posed, which the re-fit report must state.

### 6.2 Step 1 — Method-B re-fit per form (zero MD; dual-purpose)

Re-extract each form's coefficients from the **existing 9/18 Å TDDFT
traces** with the Tier-0 Method-B machinery
(`METHOD_B_trajectory_matching_extraction.md`), over the calibrated
bands (18 Å: 2.54–3.02 Å/ps; 9 Å: 2.83–4.95 Å/ps), plus the low-v
re-inspection of the trace tails (do they show coast/threshold
behavior?). Outputs: per-form coefficients + fit residuals per band.
**This is simultaneously RQ11 lever-hierarchy step 1** (form
discrimination — TDDFT as the only in-band authority, I100); the RQ11
reading goes to `RESEARCH_QUESTIONS.md`. Stated limitation: the
calibrated bands only partially overlap the RQ11 deficit window
(v ≈ 2–7 Å/ps), so low-v discrimination power may be limited — which
the trace-tail re-inspection and the standing collaborator ask (TDDFT
tails / extended traces) are there to supplement.

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

### 6.5 E_bind scan (MD)

`binding_energy_I_ion_eV` OAT around the jointly-extracted 0.1168 eV:
**{0.0584 (×0.5), 0.0934 (×0.8), 0.1168, 0.1402 (×1.2), 0.1752 (×1.5)}**
— 4 new MD cells, N = 500 each. Runs deliberately under
`allow_unvalidated_binding_pairing` (the §6.5.1 guard exists because
drag ↔ E_bind is a jointly-calibrated pair; the flag's use is the
documented, intentional exception, mirroring the §6.6 mass-pairing
precedent). Focus observables: trapped fraction, suppressed class, exit
deceleration → deep-bin KE (E_bind is the direct knob on the
trapped/exit boundary — the §4ee un-trapping echo makes this the
cleanest probe of the candidate-(iv)-adjacent margin class). Caveat
stated in every table: these are sensitivity reads of a broken joint
calibration, not candidate points; E_bind is *not* connected to the
He ladder D₀(n)/Σ (different bookkeeping surfaces: droplet mean-field
exit well vs per-atom shell energetics).

## 7. Execution order and budget

| stage | content | MD cost | gate |
|---|---|---|---|
| 1 | D0 compact reference doc | zero | none (doc work) |
| 2a | D4 Step 1 Method-B re-fit + trace-tail re-inspection | zero | trigger (analysis scripts) |
| 2b | Axis A grid (§3) — runs parallel to 2a | 9 × 500 | trigger (gen + report scripts) |
| 3 | D2b audit + A/B + re-weighting (+ ≤ 2 confirmations) | ≤ 2 × 500 | trigger (sampling-variant runs) |
| 4 | Axis B E₀/τ curves | 8 × 500 | trigger |
| 5 | D4 Step 2 twin sweeps | zero | none (scratchpad twin, §4ee precedent) |
| 6 | D4 Step 3 spot-checks + E_bind scan | ~4–8 × 500 | trigger (incl. enum build) |
| 7 | synthesis: merge all results into `TIER2_PARAMETER_INFLUENCE.md`, close GAP markers | zero | none |

Total new MD ≈ 23–27 cells × N = 500 ≈ 11.5–13.5k fragments ≈ 2.3–2.7×
the pooled battery — the explicit spend of the program. N = 1000
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
