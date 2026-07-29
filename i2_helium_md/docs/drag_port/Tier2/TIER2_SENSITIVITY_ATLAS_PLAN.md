# Tier 2 — Sensitivity Atlas (parameter-influence study program)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger, given per execution stage. Module
> descriptions are *interface contracts*; grids are *specifications*.
>
> **Status: DESIGNED 2026-07-23** (brainstorming session, user-approved
> section by section). **Executed since:** stage 1 (D0), stage 2a
> (D4 Step 1), §6.6 counterfactual (twin + MD spot-check), §6.7 items
> 1–2, the **D2b §4.1 provenance audit**, the **Axis A pre-read**, the
> **G0 spec freeze**, **stage 2b (built)**, **Axis A / G1** and the **§3.5b
> retained-class arm + marginal bracket (2026-07-27)** — the geometry grid is
> **complete at 11/11 cells**, the last five recovered detection-only at zero
> MD — and the **D2b §4.3 grid re-weighting (2026-07-27, zero MD)**: the
> pre-registered oracle is INADMISSIBLE (53.8 % of the standing density sits
> below the grid support), the post-hoc in-support oracle passes 7/7 (method
> valid where support covers; `nearest` adopted), and the corrected-ensemble
> forecast is on record (findings "D2b §4.3"; D0 §15.7 — the landing breaks
> at ensemble level, trap 0.31–0.42, W₁ ≈ 9.0–9.8, deepKE 1.80–1.90 as Arm
> A–B ranges). **Next:** G2 (the pre-registered adoption decision) can now be
> taken on the direct G1 rows plus the ensemble forecast, with the §3.5b
> bracket caveat attached to the R ≥ 49 Å rows. Then the remaining D2b A/B →
> Axis B E₀/τ curves (sequencing question open: G3 re-locates (v_c, τ, E₀),
> so Axis B around the pre-G3 point may be deferred past G3/G4).
> **G2 TAKEN 2026-07-27 (user): the corrected geometry is ADOPTED as the
> target** (retained-policy sub-decision deferred to G3; `finc1v725`
> stands until G4). The **G3 twin-first scan** is designed along
> the Route A (cascade) / Route B (selection) split of §3.5 — mechanism
> basis in D0 §14.3 (birth-depth lever = geometry-locked dressing +
> parameter-accessible transit; the §4p receipts). **G3 Step 1 EXECUTED
> 2026-07-27 (zero MD, triggered):** the twin landmark re-issue at the
> corrected geometry — S6 oracle bit-exact, the 11 G1 cells re-run
> through the twin (the long-chord twin↔MD authority box: supp/n₁_solv
> near-quantitative, n̄ residence-scale hot, trap a twin floor with the
> center-pin trap channel twin-invisible, KE direction-only), and the
> corrected-ensemble twin row confirming the D2b §4.3 forecast
> cross-instrument (findings "G3 Step 1"; stage `g3landmarks` in the
> committed twin script). **G3 Step 2 DESIGNED + user-approved
> 2026-07-27 (§3.5c): the nested Route A/B twin scan** — chord surface
> (v_c × E_bind, incl. two outside-Tier-0-authority diagnostic v_c
> arms) × free surface (τ × E₀), zero-integration pre-scan first,
> pre-registered gate n₁_solv ∈ [0.19, 0.30] ∧ n̄ ∈ [4.4, 7.1] and
> failure criterion. **This is the next step**, awaiting its
> `[PROCEED TO IMPLEMENTATION]`.
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

**Purpose restated (user, 2026-07-26).** The goal is **model
understanding and physical soundness in general** — not RQ11. RQ11 is
one symptom among several. Two consequences, binding on the rest of the
program: (i) a study is worth running if it tells us what the model is
doing or whether an element is physically justified, *independently* of
whether it moves the deep-bin deficit; (ii) the program carries a second
deliverable beside the influence map — the **physical-sensibility
ledger** (`TIER2_PARAMETER_INFLUENCE.md` §17), which classes every
element as physics-constrained / convention / effective / scaffolding /
missing, records what would retire it, and is maintained as results
land. The atlas still adopts nothing (see the stance below), but it now
explicitly produces the **dossier** any later adoption decision would
read.

**Stance (fixed for the whole program):**

- **Sensitivity atlas — reported, not adjudicated.** No parameter is
  re-tuned by this program. `finc1v725` stays the standing production
  point; nothing here discharges F5 or moves production. If an atlas
  result motivates moving the standing point, that is a *separate,
  pre-registered* decision outside this plan.
- **Amended 2026-07-26 — geometry is the one declared exception.** The
  droplet geometry is *known* wrong from the experiment's own source
  conditions (D0 §15), so **Axis A is re-scoped as the geometry
  correction** (§3, staged G0–G4 in §3.5; user-adopted in principle).
  It is read for both purposes at once: what the corrected geometry
  does, and what the variables do to the model. Adoption still waits for
  the pre-registered G2 decision; every other knob remains under the
  stance above.
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
- midHot (**geometric** mean of sim/ref-mean KE over **n = 2–8** — corrected
  2026-07-26: this doc previously said "n = 2–3", which is not the recorded
  convention; both the band and the geometric aggregate are now committed
  code, `postprocess/tier2_confirmation.midhot_ratio`), W₁_solv
  (solvated-histogram Wasserstein),
- n₁ KE ratio (median-anchored, I77 convention),
- deep-bin KE ratio sim/ref over n = 10–17 (the RQ11 axis; **arithmetic** mean
  of the per-bin ratios — `deep_bin_ke_ratio`),
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

## 3. Axis A — droplet geometry: influence map **and** the correction program (G0–G4)

Everything so far sits at the production-sampled geometry; size and
position have never been controlled independently. The size–position
pair sets the He path length and hence the total drag/cooling exposure —
the largest unmapped influence on the histogram and on the n = 1
velocity distribution.

**Re-scoped 2026-07-26 (user).** The provenance audit (§4.1) and the
parent-document anchor established that the geometry is not merely
unmapped, it is **known wrong**: the standing point runs at ⟨N⟩ = 2000
(R̄ 26.6 Å, mean solvation depth 9.0 Å) against the experiment's own
source conditions of ⟨N⟩ ≈ 12794 (R̄ ≈ 49.7 Å, depths 29–40 Å). Axis A
therefore carries **two jobs at once**, both declared:

1. **influence** — what droplet size and birth position do to the model
   (the atlas job), and
2. **correction** — the scoping stage of fixing the geometry, staged
   G0–G4 in §3.5.

The stance amendment in §0 applies to this axis only: geometry is
corrected under G0–G4; every other knob stays reported-not-adjudicated.

### 3.0 Why the correction simplifies the model

Five rows of the D0 §17 ledger classed **S** (scaffolding) retire
together at the corrected geometry, with no new knob introduced:

| element | why it exists today | fate after correction |
|---|---|---|
| ⟨N⟩ = 2000 pin | Tier-0 TDDFT droplet → S2-D4 freeze → twin D4 family; its only source-consistency check was evaluated at 23 K instead of the preset's 14 K | replaced by the source correlation |
| `uniform_volume` birth law | T7 twin parity, because the parent's thermal law center-pins at R ≈ 28 Å | replaced by the parent's own law |
| margin 3 Å | free convention, load-bearing (I88: 3.4σ on n₁_solv at 6 Å) | disappears with the birth law |
| T5 `density_tied` | dressing for under-dressed shallow births | structurally inert (ρ̂ saturates) |
| T6 `sigma_proportional` | onset scaled by Σ(n₀)/Σ(n*) | structurally inert |

### 3.1 Grid — **option (B), adopted 2026-07-26 (user decision)**

Axis A is the **standing-vs-anchored-geometry test**, not a local
sensitivity ring: the R axis spans from the inherited geometry to the
parent document's own droplet ensemble, and the position axis is run as
competing *laws* rather than fixed offsets.

**3 × 3 factorial + 2 large-R cells = 11 cells, N = 500/cell, one fixed
seed shared by every cell, full pipeline:**

- **R (droplet radius), fixed per cell:**

  | cell | R [Å] | N [He] | provenance |
  |---|---|---|---|
  | R1 | **26.6** | 1727 | the standing point's realized mean (⟨N⟩ = 2000 prior) |
  | R2 | **34.0** | 3605 | the parent document's *smallest* droplet — externally anchored (its quoted mean solvation depth 29 Å) |
  | R3 | **49.4** | 11059 | the parent ensemble's mean radius (raw ln-normal at ⟨N⟩ = 12794; realized ⟨R⟩ = 49.7, median 48.6 Å) |

  Conversion is the bulk form 2.2173·N^(1/3) throughout (D0 §15.3).
  **4th/5th cells FUNDED 2026-07-26 (user): R = 68.3 Å (N = 29227, the
  parent's largest) on L2/L3** — they are the only cells that close the
  D2b re-weighting extrapolation (the parent's q95 otherwise sits outside
  the grid support, §4.3) and the trapped-explosion question (§3.6) lives
  at large R.

- **Birth law, per cell** (replaces the r/R offsets):

  | cell | law | configuration |
  |---|---|---|
  | L1 | **center-pin** | `single_initial_position=True` (r₀ ≡ 0) — the Tier-0 convention and the small-R limit of the parent law |
  | L2 | **Boltzmann, parent-consistent** | `birth_position_law="boltzmann"` with `binding_energy_molecule_K=313.2` (β₂ = 26.99 meV, the DFT fit the parent used) |
  | L3 | **`uniform_volume`, m = 3 Å** | the standing production law |

  L2 at R2 reproduces the parent's quoted 29 Å depth (29.33 Å) — the
  cell that makes the grid externally checkable. The as-built well
  573.3 K = 49.4 meV is a **robustness check, not a cell** (+1.3 Å in
  depth, D0 §15.4); **G0 decision: 313.2 K is a per-run override, the
  config default stays 573.3 K** (legacy port fidelity; the arm is unused
  in production, so a default change buys nothing).

  **L1/L2 near-degeneracy, kept deliberately (2026-07-26).** Their mean
  depths and chords nearly coincide at every R (26.6/26.6, 34.0/33.8,
  49.4/47.9 Å) and both saturate the dressing, so the pair differs
  essentially only in birth-depth **spread** — which is precisely §3.3's
  "how much spread is geometry-inherited" contrast. User decision: keep
  all three L1 cells for a balanced factorial rather than reallocating
  L1/R2.

- **Pre-registered exposure table** (mean isotropic chord from birth
  point to surface — the quantity drag integrates; computed at design
  time, so the MD reads can be checked against it):

  | | L1 center | L2 Boltzmann | L3 uniform m3 |
  |---|---|---|---|
  | **R1 = 26.6 Å** | 26.6 | 26.6 (depth 24.8) | **21.4** (depth 8.9) ← standing geometry |
  | **R2 = 34.0 Å** | 34.0 | 33.8 (depth 29.3) | 26.9 (depth 10.7) |
  | **R3 = 49.4 Å** | 49.4 | 47.9 (depth 35.1) | 38.3 (depth 14.5) |

- **Pre-registered structural expectation** (from the §3.4 pre-read,
  findings §C): the birth-law column is the strong axis — the escaping
  fragment's path *is* its birth depth — while the R row acts mainly
  through trapped/suppressed selection, with a residual n̄ slope at deep
  births.
- **Pre-registered numeric prediction for the anchored cells** (frozen
  before launch): at L2/R1 the dressing saturates (ρ̂ = 0.994,
  n₀ = 20.9 of 21) ⇒ supp ≈ 0, n̄ > 6.3, n₁_solv ≲ 0.1, W₁ ≳ 1.6 — i.e.
  **the landing breaks at the anchored geometry, through the
  dressing/suppression channel, not through drag exposure**. A miss that
  matches this prediction is a confirmation, not a surprise.
- **Scoring rule (mandatory, from pre-read §D).** Pooled W₁ 0.571 is
  better than *every* geometry bin (best 0.639). The landing is a
  mixture property, so fixed-geometry cells are compared **to each other
  and to a re-weighted mixture reconstruction**, never directly against
  the committed acceptance. Scoring nine fixed cells against the
  acceptance would manufacture nine spurious failures.
- **Pre-registered observable resolution (added 2026-07-26).** §6.7
  item 2's finding 4 measured that at N = 500 with a single seed **W₁
  swings ≈ 0.1 on seed alone** — the refuted single-seed reads in that
  scan were W₁ reads. Therefore: every cell uses **the same seed**
  (common random numbers, so cross-cell differences are partly paired);
  per-cell W₁ is **reported but not read as a discriminator** below
  |Δ| ≈ 0.15; the axis's discriminating observables are **supp, trap, n̄,
  n₁_solv, midHot and deep-bin KE**, which resolved at 8–35σ in that same
  scan. This is consistent with §3.1's mixture rule — per-cell W₁ was
  never the target.
- All other knobs at the standing point. Droplet-size *sampling* is
  disabled per cell (fixed R) — that is the controlled variable.
- **E_bind held at 0.1168 eV in every cell** (G0 decision, §3.6): its
  R-dependence is bounded analytically at ≤ 0.009 eV over the whole span
  (D0 §9.1), with an MD bracket cell **armed as a rider** on the G1 read.

### 3.1b Grid-centering decision

**(B), 2026-07-26 (user):** span the anchor rather than ring the
standing point; the superseded option-(A) grid and the full option
comparison are in the migration log.

### 3.1c Dressing / onset / cooling control cells — rejected

**Rejected 2026-07-26 (user):** each requires the *less* physical arm of
its pair, and the three legs all read one shared surface
(`rho_he_ratio(depth, steepness=drag_gate_steepness(cfg))` —
`ion_initial_state.py:204`, `detection_stage.py:454`, `ion.py:211`), so
they are one physical fact and not independently variable; the
correction (§3.5) makes them saturate instead.

### 3.2 Machinery — **BUILT 2026-07-26** (stage 2b, triggered)

Sampling machinery already existed (`sampling/droplet_sizes.py`,
`sampling/radial_positions.py`, `single_pulse_droplet_distribution`
preset). Delivered under the trigger:

- **`scripts/gen_tier2atlas_geometry.py`** — the 11-cell generator (atlas
  namespace, skip-complete/overwrite guards, `--dry-run`, process pool),
  seed **20260727** shared by every cell, each cfg diffed field-by-field
  against the standing battery member.
- **`scripts/post_processing/tier2atlas_geometry_table.py`** — pure-scorer
  report: the pooled-battery drift oracle first, then the grid table
  (geometry columns + committed observable vector), CSV behind a setting.
- **`droplet_size_sampler_mode`** config field + plumb (see below).
- **`postprocess/tier2_confirmation.ke_band_ratio` / `midhot_ratio` /
  `deep_bin_ke_ratio`** — midHot and deep-bin KE were scratchpad-only
  conventions until now; committing them (with the pooled oracle reproducing
  midHot 1.0110 and deep-KE 0.6313) removes the drift risk the atlas tables
  would otherwise carry.
- **`scripts/tier2_common.cfg_diff_vs_reference`** — the previously
  duplicated "only these knobs may differ from the reference run" guard,
  now one implementation shared by all five atlas generators, with the
  forward-compatibility rule a new `SimConfig` field requires.

**Every G1 cell is expressible with existing fields**
(`single_droplet_size` + `use_single_droplet_size`,
`droplet_size_prior="legacy"`, `birth_position_law`,
`single_initial_position`, `initial_position_margin_angstrom`,
`binding_energy_molecule_K`); the fixed-offset override knob that option
(A) would have required is **not** built. Note
`binding_energy_molecule_K` is currently read only by the Boltzmann arm —
L2 is its first non-default use in the drag branch.

**Three guards are exercised deliberately** and must be verified at
generator-test time, not discovered at runtime (config line numbers as
read 2026-07-26):

1. T8 analytic-vs-fixed-size: the analytic prior arms require
   `use_single_droplet_size=False` (`config.py:1142`) ⇒ every fixed-R
   cell must run `droplet_size_prior="legacy"`.
2. `uniform_volume` requires `single_initial_position=False`
   (`config.py:738`) — L3.
3. **`birth_position_law="boltzmann"` requires
   `initial_position_margin_angstrom == 0.0`** (`config.py:731`) — so L1
   and L2 cells cannot carry the 3 Å margin. Newly noted.

**One new config field (G0-1 decision, 2026-07-26 — amends the earlier
"no new config surface" claim).** The corrected *ensemble* arm needs
`mode="raw"`, which is unreachable from config today:
`simulation/initial_state.py:83` hardcodes `mode="post_pickup"` for the
`legacy` arm. Adopted: a `droplet_size_sampler_mode ∈ {raw,
post_pickup}` selector defaulting to `post_pickup`, so no existing run
changes (no default-scope change) — built with the stage-2b generator
under its trigger, though G1's fixed-R cells do not use it. Rationale in
D0 §15.6: under `legacy` the mean is not a knob, the nozzle correlation
supplies ⟨N⟩ = 12794 from the preset's own 40 mbar / 14 K, so the
corrected ensemble needs **no new number** — and the analytic prior
family retires with the pin rather than needing its truncation window
widened.

### 3.3 Questions the axis answers

- **Does the landing survive at the anchored geometry?** (the option-(B)
  headline: R3/L2 is the parent document's own droplet ensemble and
  birth law, at ≈ 2.2× the standing point's He chord. Reported, never
  re-tuned.)
- Does exposure alone (bigger droplet, deeper start) reorganize the
  detected histogram, and along which observables first?
- Is the n = 1 velocity distribution's shape a **geometry fingerprint**
  (distinguishable signatures of R vs birth law), or degenerate along an
  exposure-equivalent diagonal?
- Is the margin-3 Å pin (I88) still a live lever once the birth law is
  the parent's own — i.e. does L2 make the margin convention moot?
- Does geometry move the deep-bin cold tail (RQ11-relevant; reported
  under atlas stance, not adjudicated)?
- How much of the production ensemble's spread is geometry-inherited vs
  mechanism-stochastic?

### 3.4 Battery pre-read — EXECUTED 2026-07-26 (zero MD)

Formerly "optional follow-up", run first because it is free: the pooled
N = 5000 battery binned by its own sampled (R, depth). Committed scorer
reused, oracle reproduced. Results in findings "Axis A pre-read" §A–§G;
the three that changed this plan: R is a *selection* knob and birth
depth the physics knob (§3.1 pre-registration rewritten); the histogram
landing is a **mixture** property (§3.1 scoring rule); and the deep-bin
KE deficit closes ×2.5 with birth depth.

### 3.5 The correction, staged (G0–G4)

**G0 — freeze the corrected-geometry spec (doc, zero MD). FROZEN
2026-07-26.** All four decisions taken:

1. **Size prior → `legacy` + `raw`.** The parent's quoted quantiles match
   **raw** (q05 34.9 / q95 68.7 Å vs its 34 / 68.3); post-pickup gives
   37.7 / 74.1. The analytic arm is now *rejected*, not merely blocked: at
   ⟨N⟩ = 12794 with δ = 0.625 roughly a quarter of the ln-normal mass sits
   above its [250, 16000] He window, so it would need the window widened
   **and** the twin's D4 family re-pinned. Under `legacy` the mean is not
   a knob at all — the nozzle correlation supplies ⟨N⟩ = 12794 from the
   preset's own 40 mbar / 14 K, so the corrected ensemble needs **no new
   number**, and the analytic prior family retires together with the ⟨N⟩
   pin (a sixth **S** row; D0 §15.6 + ledger §17). Cost: the one
   `droplet_size_sampler_mode` selector of §3.2, since `mode="raw"` is
   unreachable from config today (`initial_state.py:83`).
2. **Birth law + well depth → `boltzmann` at 313.2 K as a per-run
   override**; the config default stays 573.3 K (legacy MATLAB value —
   changing it is a default-scope change with no benefit while the arm is
   unused in production). The ~1.3 Å systematic stays a provenance-defect
   row in the ledger; 573.3 K is a robustness check, not a cell.
3. **E_bind → defer, bounded analytically.** All cells run at 0.1168 eV.
   The R-dependence is bounded at **≤ 0.009 eV (7.6 %) across the whole
   grid span including R = 68.3 Å**, by splitting the well into an
   R-independent local (snowball / electrostriction) term and a Born
   far-field term `c/R`, `c = (q²/8πε₀)(1 − 1/ε) = 0.390 eV·Å` at
   liquid-He ε = 1.0572 — the local term cancels in the *difference*,
   which is the only quantity claimed. Against the measured lever
   (trap ≈ 0.85 per eV) that is Δtrap ≈ 0.006, versus the 0.0003 → 0.158
   R-selection effect the grid measures, and the sign is one-sided
   (deeper at large R). Free by-product: the same formula gives 0.096 eV
   at a = 3.5 Å against the extracted 0.1168 eV — the first independent
   number on E_bind, ~20 %. **Rider armed, not funded:** an MD bracket
   cell at R3 fires only if G1 shows the trapped channel near a boundary
   (trapped > ~0.4, or a detection-handover-guard trip like the one
   0.154 eV produced in §6.7 item 2). Derivation: D0 §9.1.
4. Source conditions p = 40 mbar / T = 14 K stay (they are the preset's
   own values and reproduce the parent ensemble).

**G1 — scoping study = the §3.1 grid (MD).** Read twice per cell: the
corrected-geometry observables, and the influence map. This is the only
MD this axis spends before a decision.

**G2 — pre-registered adoption decision (geometry only). TAKEN 2026-07-27
(user): ADOPTED.** Registered in advance: the landing is *expected* to
break at the anchored geometry (§3.1 prediction), and **a broken landing
is not a reason to keep the wrong geometry** — it is the reason G3
exists. Decision record: the corrected geometry (the G0 spec — `legacy`
+ `raw` sizes at the preset's own 40 mbar / 14 K ⇒ ⟨N⟩ = 12794,
Boltzmann 313.2 K births, E_bind 0.1168 eV) is the **target
configuration**; G3 is authorized as the next stage (behind its own
trigger); `finc1v725` stands until G4 delivers its successor. Dossier
read: the external anchor (D0 §15), the confirmed G1 prediction, the
§4.3 ensemble forecast, the §14.3 mechanism decomposition.
**Sub-decision deferred (user):** the detection retained-ion policy is
NOT fixed into the corrected-geometry definition at G2 — it stays a
per-run choice (`exclude_all_coupled` for atlas/G3 scoping runs, with
the bound/marginal decomposition and the item-6 caveat) and is decided
at **G3**, when the re-arbitrated drag point shows how large the coupled
class really is.

**G3 — re-arbitration at the corrected geometry.** (v_c, τ, E₀) were
fitted at the wrong exposure and must be re-located, not carried over.
Cost ladder as always: the twin's K landmarks are re-issued at the
corrected geometry (they are all center-pinned ⟨N⟩ = 2000 numbers
today), the basin located there, then a small MD confirmation ring —
twin free, MD ~10–20 cells × 500. **Step 1 (landmark re-issue + G1
deep-cell oracle check) EXECUTED 2026-07-27** — stage `g3landmarks`,
findings "G3 Step 1": the re-issued center-pin landmark table is
law-tagged (the recorded K 0.74460 is pure-cubic; the capped tail is
what makes the anchored radii traversable), the twin↔MD authority box
is re-measured on all 11 G1 cells, and the corrected-ensemble twin row
confirms the D2b §4.3 forecast cross-instrument. Future corrected-
geometry twin sessions oracle against `h2b_g3_corrected_row.csv`.

**G3 design sharpened (2026-07-27 discussion; mechanism basis D0 §14.3).**
The working hypothesis G3 tests, held by the user: *a parameter point
exists at the realistic geometry that reproduces the landed observables.*
The birth-depth lever decomposes into a geometry-locked channel (partial
dressing — saturated at realistic depths, no knob reaches it) and a
parameter-accessible channel (transit/cascade — exactly the (v_c, τ, E₀,
E_bind) surface), so the twin scan is read along two routes:

- **Route A (cascade):** does any (v_c, τ, E₀) corner shed ~10 more
  He/ion and cool deep arrivals? Energetically plausible (≈ 0.05–0.15 eV
  extra E_int vs E₀ = 0.27); the deepKE sign flip (0.63 shallow →
  1.8–1.9 deep, crossing at depth ≈ 11–15 Å) means the arbitration
  passes through the right value. The small measured in-tier lever
  slopes were all short-transit numbers and do not bound this.
- **Route B (selection):** does a trap-boundary corner (E_bind, v_c —
  the E_bind trap lever is measured at 0.85/eV) make predominantly the
  shallow-born minority detectable, with the marginal class's
  at-or-below-support KE (§3.5b item 9) as the physical detector-side
  argument? This is the emergent version of what the ⟨N⟩ = 2000 prior
  was hand-doing, and the physically-attractive outcome.

Acceptance framing for the scan: n₁_solv / n̄ / W₁_solv and the KE curve
are the hard targets; **supp is a softer, RQ3-coupled target** (the
experimental n = 0 bin is a two-channel mixture; the solvated histogram
renormalizes to n ≥ 1), so supp → 0 at saturated dressing is not by
itself disqualifying. Failure of both routes is itself the atlas
result: it localizes the corrected-geometry failure to the mechanism
(dressing law, pickup, ladder) and reads the standing point as an
effective model of the detected subset.

**G4 — re-baseline.** New standing point, fresh 5 × N = 1000 battery,
figure surface re-generated, D0 §17 re-issued with the retired **S**
rows struck.

**Protected — not re-fit by the correction:** Tier-0 drag form and b
(calibrated per unit density behind the ρ̂ gate; the interior is bulk at
any R, so the extraction transfers unchanged), Tier-1a, the committed
scorer, checkpoint schema, RNG draw order, the constants table.

**Discipline:** nothing adopts before G2; `finc1v725` stands until then;
each stage behind its own `[PROCEED TO IMPLEMENTATION]`.

### 3.5b The handover-staging decision — SPEC FROZEN + AMENDED + **DELIVERED 2026-07-27** (arm + bracket; grid complete at 11/11, zero MD)

> **Outcome:** arm built, all 11 cells scored, six reproduced bit-for-bit
> (item-10 oracle). The item-8 bracket came out **NOT TIGHT** at r3l1 / r4l2 /
> r4l3, so item 6's caveat **stays** and the bracket columns stay; the
> `deepKE`/`n̄` values at R ≥ 49 Å are conditional on the exclusion. The
> retained class is confirmed slow and high-n, and measured to be **not** a
> null direction — it is the RQ11 direction. Results: findings §G1.4, D0
> §14.1–§14.2.

G1 is complete for R ≤ 34 Å and **blocked** above it: all five anchored-radius
cells trip the detection P1–P3 handover guard (mechanism and numbers in D0
§14.2 / findings §G1.3). This sub-section records the adjudication and the
frozen spec so the build can proceed without re-deriving either.

**User's interpretation (2026-07-27), which the numbers are to be read under.**
Big droplets combined with an **over-strong drag law** produce a large retained
population; those ions' trajectories are *not interesting* and should not be
computed; this is a **long-standing design error, not a discovery**. The only
quantity that matters about them is **the percentage that gets stuck**.

**Frozen spec.**

1. **New arm on the existing enum surface:**
   `detection_droplet_retained_policy = "exclude_all_coupled"`. Every ion still
   helium-coupled at handover is excluded from the free-flight loop and
   counted; none is integrated further; no error. Default stays `refuse`;
   the atlas geometry cells opt in, so **no existing run changes**.
2. **The retained class stays decomposed** — two different claims must not be
   merged:
   - **bound** — provably cannot escape (exact conservative test: total energy
     vs the effective potential incl. angular momentum). Physics; the existing
     V0-2 class, `state_reason = "droplet_retained"`.
   - **marginal** — energetically able to escape but still inside helium at
     handover. A *modelling* exclusion, conditional on the drag law and the
     window length; its own `state_reason` value.
   The scorer counts both as retained through **one shared constant** (not two
   code paths), and the report prints `trap_bound` / `trap_marginal`
   separately. Rationale: collapsing them makes the retained fraction
   unreadable later — no one could tell physics from convention.
3. **No checkpoint schema change** — a new *value* in the existing
   `state_reason` string field; no new fields, no dtype change. It does extend
   the reason vocabulary, so `postprocess/tier2_confirmation` must be updated
   in the same change or it would silently mis-handle the new value.
4. **Zero new MD.** The five blocked cells failed at the *detection* step and
   their `relaxation.npz` (a v7 ion checkpoint) is on disk, so recovery is a
   detection-only re-run — minutes, not hours.
5. **Cannot retroactively alter the six scored cells.** At R ≤ 34 Å every
   violator was provably bound — that is why those cells passed — so the
   marginal class is empty there and the R1/R2 numbers are untouched.
6. **Caveat travelling with every table:** the marginal fraction is
   conditional on a cubic law extrapolated ~3× beyond its 9/18 Å calibration
   band. It is a diagnostic of the model at this geometry, **not** a
   prediction of nature. Item 8's drop rule is what can retire this caveat.

**Amendment 2026-07-27 — items 7–10.**

7. **The arm buys no compute** (correction to the recorded interpretation).
   "Should not be computed" is the right *intent*, but the arm cannot deliver
   it: it touches only the **detection** step, and the retained ions' cost was
   already paid in the neutral + 30 ps ion MD, which ran before the guard
   fired. Two candidate savings were checked and both fail:
   - the relaxation-stage early exit is **not** blocked by the retained class —
     "frozen" there means the *evaporation cascade* froze (`n == 0` or
     `E_int < D₀(n)`), not spatial decoupling. (**Corrected at build time
     2026-07-27:** an earlier draft of this item cited the config *default*
     `cooling_spatial_gate = "none"`; these cells actually run
     `"density_scaled"`. The conclusion is unchanged and in fact firmer — under
     `density_scaled` it is the *ejected* fragments whose cooling stops as
     ρ̂ → 0, so they are the ones that can fail to freeze, while retained ions
     sit at ρ̂ ≈ 1 and freeze **fastest**.)
   - an MD-stage early abort on provably-bound ions is not free either: the
     class is not identifiable until late in the window, and the stage is
     vectorized, so a real saving would need array compaction.

   Consequence for sequencing: **do not re-run any cell to test the arm.** The
   five blocked cells hold `relaxation.npz` (v7), so the arm is a
   **detection-only re-run — minutes, zero MD** (§3.5b item 4). A fresh N = 500
   would cost ~44 min/cell and measure nothing the stored states do not
   already contain.

8. **One-off bracket diagnostic — the exclusion's own bias check.** Scoring the
   arm proves nothing on its own: the excluded ions are absent from every
   observable *by construction*, so the scored vector cannot report whether the
   exclusion biased it. The check is nearly free because the arm already
   computes the split — `_conservatively_bound` (`detection_stage.py:640`)
   evaluates each ion's `E_tot` against `V_eff(r') = U(r'−R) + L²/2mr'²`, and
   because the E2 relaxation is **zero-gamma** (`relaxation_stage.py:43`) that
   same `E_tot` **is** a marginal ion's exact conservative asymptotic KE. No
   integration, no new physics.

   - **Arm A (headline):** marginals excluded — the committed observable vector.
   - **Arm B (diagnostic):** marginals injected into the *scored* ensemble at
     `(n = n_handover, KE = E_tot_handover)`. Bound ions never enter Arm B.
     Scoring-side only — it writes no checkpoint and fabricates no detection
     rows. `(n, KE)` is the complete input to every observable in §1.1, so the
     injection is exact rather than approximate.

   **Pre-registered before the run** (otherwise "tight" is post-hoc):
   - *Direction:* Arm B adds low-KE, mid-to-high-n fragments ⇒ deep-bin KE
     **down**, n̄ slightly up, midHot ≈ flat, supp unchanged.
   - *Excluded from the test:* the **fate split**. Arm B reclassifies
     marginal → detected by construction, so that Δ is definitional, not bias.
     The test runs on the detection-conditional shape observables: n̄,
     n₁_solv, midHot, deep-bin KE, W₁_solv, χ²_med.
   - *Tight ⇔* bracket width < 1 seed-SD on each of those, using the pooled
     N = 5000 SDs scaled by √(N_pool / N_cell_detected), **and**
     |ΔW₁| < 0.15 (§3.1's single-seed W₁ resolution floor).
   - *Drop rule:* **tight** ⇒ record in the findings doc, delete the
     diagnostic, and retire item 6's travelling caveat. **Not tight** ⇒ the
     diagnostic stays, the offending observable is flagged, and the follow-up
     is the kinetic forward model (conservative orbit from the stored state
     with live Poisson pickup at λ₀ρ(depth(t)) and RRK evaporation) — which
     also tests whether pickup **mass-loads marginals into the bound class**,
     i.e. whether the conservative split *understates* trapping.
   - *As built (2026-07-27):* the diagnostic is its own script,
     `scripts/post_processing/tier2atlas_retained_bracket.py`, so the drop
     rule deletes a file rather than a column — the grid report never carries
     bracket columns. The `trap_bound` / `trap_marg` columns it *does* carry
     are item 2's permanent decomposition and are **not** covered by the drop
     rule. Outcome: **NOT TIGHT**, so the script stays.

   Expected split of the outcome, stated so it is not read as a surprise: the
   claim "the retained class does not shift the result" and the claim "it is
   just slow, high-n ions" pull against each other — slow and high-n places
   them in exactly the deep-n / low-KE bins where RQ11 lives and where r3/r4
   currently read deepKE 1.38–1.49 (too hot). A near-zero Δ on W₁/n̄/supp with
   a resolved Δ on deep-bin KE is therefore the *likeliest* result, and the
   bracket must be read per observable, never aggregated.

9. **Third output, free from the same array — the experimental-acceptance
   check.** Compare the marginals' `KE_asym` distribution against the
   experimental reference's KE support. If it sits below that support, these
   fragments are undetected **in the experiment** too, and exclusion stops
   being a modelling convention: it becomes the correct model of the
   measurement. This is the one result that could retire item 6's caveat on
   physical rather than statistical grounds.

10. **Arm oracle.** Re-running detection on the six already-scored R ≤ 34 Å
    cells under the new arm must reproduce them **bit-for-bit** — every
    violator there was provably bound, so the marginal class is empty at those
    radii (item 5). This is the arm's own regression test, in the scorer-oracle
    idiom of §1.4.

**What the spec deliberately does not resolve.** Whether the drag law
over-dissipates over 25–50 Å paths. That is the standing collaborator ask
(extended-range TDDFT), not something this axis can settle; the arm makes the
geometry grid readable while leaving the question open and flagged. The
distinguishing test, if it is ever wanted cheaply: a geometric retained
fraction is insensitive to (v_c, b), an over-dissipation artifact is not.

### 3.5c G3 Step 2 — the Route A/B twin scan (DESIGNED 2026-07-27, user-approved; **EXECUTED 2026-07-27** — 24/6480 cells gate, basin Tier-0-legitimate; findings "G3 Step 2")

Design frozen from the Step-1 numbers (findings "G3 Step 1"); scan spec
adjudicated by the user 2026-07-27.

**Structural basis (measured, Step 1).** In the twin, trap and the KE
curve depend only on the chord integration — the **(v_c, E_bind)
surface**; τ (exact K rescale) and E₀ (fate map) are **free
post-processing**. So Route B = the chord surface, Route A = the free
surface at fixed chords, and one nested factorial covers both at the
cost of the chord families alone. Twin caveat pre-registered: in MD, τ
feeds back on trap through the mass mechanism — the clean decomposition
is itself a twin approximation (authority box, D0 §14.4).

**Fixed:** the committed corrected master draw (seed 20260727, m =
20000), rq4graded ladder, p = 1, p_tail = −1, detection = the
non-trapped chord read.

**Grids.**

- **Chord surface:** v_c ∈ {5.0, 5.5, 6.0, 6.5, 7.25, 8.0, 9.0, 10.0}
  (legitimate range; the floor respects the TDDFT band top 4.95 Å/ps)
  **plus v_c ∈ {3.5, 4.25} labeled "outside Tier-0 authority —
  diagnostic only"** (they override the calibrated band; if only these
  land, the finding points at the collaborator ask, not at a parameter
  point). E_bind ∈ {0.048, 0.1168, 0.154} eV (the Tier-0 extracted
  spread; the §6.5/§6.7 documented joint-pairing exception, stamped per
  table). ≈ 30 integrations ≈ 40–55 min.
- **Free surface:** τ ∈ {2.4, 3.2, 4.8, 6.4, 9.6, 12.8} ps (log-ish;
  any landing that *requires* τ ≫ the sourced 6.55 is flagged against
  its calibration class, never adopted silently); E₀ ∈ [0.17, 0.52]
  step 0.01. ≈ 10⁴ scored cells total.

**Block order.**

0. **Oracles:** `h2b_g3_corrected_row.csv` re-derived bit-exact at the
   standing cell (the Step-1 landmark convention) + the S6 machinery
   oracle.
1. **Zero-integration pre-scan:** from the *stored* standing chord, the
   required-exposure-reduction map — for every (τ, E₀), the detected
   fraction below the strip ceiling, and hence how much X-reduction the
   chord axes must deliver. Analytic Route-A kill criterion before any
   new integration: if X → 0.25·X cannot produce n₁_solv ≈ 0.24 at any
   (τ, E₀), Route A is dead at this surface.
2. Chord families, then the full nested scoring.

**Per-cell outputs:** the committed observable vector + det_yield +
**the detected-subset R quantiles vs the source** (the Route-B
selection diagnostic — Step 1 showed the corrected birth law has no
shallow-birth minority; the selectable axis is droplet size).

**Pre-registered acceptance (twin level, authority box applied):**
hard gate = n₁_solv ∈ [0.19, 0.30] AND n̄ ∈ [4.4, 7.1] (target 4.07 +
the residence-conditional bias bracket [+0.3, +3]); W₁ reported, not
gating (twin-bias-loaded); KE direction-only readout (deepKE response
to Route A is not hand-predictable); supp soft (RQ3 mixture); trap
reported as a **floor** with the composition beside it.

**Pre-registered failure criterion:** no cell inside the gate anywhere
on the grid *including the diagnostic arms* ⇒ Route A and Route B fail
at the (v_c, τ, E₀, E_bind) surface — itself the atlas result: the
corrected-geometry failure localizes to the mechanism (pickup
re-filling, per-shed ε, ladder shape — the knobs the twin does not
carry), and the standing point reads as an effective model of the
detected subset.

**After the scan:** basin (if any) → MD confirmation ring ~10–20 ×
N = 500 behind its own trigger; the retained-policy sub-decision and
the Axis-B deferral question are decided on the scan's evidence (G2
record).

**Execution record (2026-07-27, stage `g3scan`, zero MD):** oracles
bit-exact; Route-A kill NOT fired (141/216 free cells reach the n₁ band
at f ≥ 0.25) but pure rescaling never lands the joint gate — the
landing needs chord reshaping. **24/6480 cells gate, all inside the
Tier-0 v_c range (diagnostic arms: 0 — the failure criterion did not
fire).** Fully-unflagged basin at the standing well: (v_c 5.5, τ 4.8,
E₀ 0.36–0.37) + (v_c 6.0, τ 6.4, E₀ 0.32–0.33); the standing chord
v_c 7.25 gates at 0/216. Landing is cascade-carried (det_yield
0.87–0.99, no R-selection). Full record: findings "G3 Step 2", D0
§14.5; the MD confirmation ring is §3.5d.

### 3.5d G3 Step 3 — the MD confirmation ring (DESIGN FROZEN + user-approved 2026-07-27; **EXECUTED 2026-07-27/28** — GR-P1..P6 all land, 4 cells gate in MD, basin CONFIRMED; findings "G3 Step 3")

The twin basin (§3.5c) confirmed or refuted in real MD, with the ring
placed by the **measured twin→MD transfer rules**, not centered naively
on the twin-gated cells:

1. n₁_solv/supp transfer near-quantitatively (§14.4 box ±0.02);
2. n̄ is twin-hot by the residence-scaled bias **+0.2…+3.2 He** and
   falls ≈ 0.3 He per 0.01 eV of E₀ ⇒ each sub-basin gets a downward
   **E₀ ladder** covering the whole bias bracket;
3. twin gates can be too narrow (§6.6 frozen-chord lesson) ⇒
   twin-*fail* controls at the basin edges;
4. trap is a twin floor and the deep-birth trap channel is
   twin-invisible ⇒ trap predicted one-sided, bound/marginal
   decomposition recorded per cell.

**Fixed:** corrected geometry (G0/G2: `legacy`+`raw` sizes at the
preset's own conditions, Boltzmann births at 313.2 K per-run override,
margin 0), standing pins otherwise (rq4graded tabulated ladder,
capped_cubic p_tail −1, co-moving shed, density_tied shell, p = 1
partition, Landau-on 0.58, budget 2.70), **N = 500 per cell, one fresh
shared seed 20260728** (CRN-paired), retained policy
`exclude_all_coupled` (interim, decomposition recorded). Namespace
`tier2atlas_conf270_g3r*`; generator `gen_tier2atlas_g3ring.py`; scorer
`tier2atlas_g3ring_table.py` (pooled-battery drift oracle first, §1.4).

**The 14 cells with their frozen twin rows** (from the committed
`h2b_g3scan_predictions.csv`; the GR-P1 oracle re-reads these
string-exact before any MD):

| cell | v_c | well | τ | E₀ | twin trap | supp | n̄ | n₁_solv | W₁ | twin gate |
|---|---|---|---|---|---|---|---|---|---|---|
| a031 | 5.5 | eb1168 | 4.8 | 0.31 | 0.0535 | 0.0008 | 6.717 | 0.0371 | 1.8778 | 0 |
| a033 | 5.5 | eb1168 | 4.8 | 0.33 | 0.0535 | 0.0124 | 5.930 | 0.1051 | 1.1757 | 0 |
| a035 | 5.5 | eb1168 | 4.8 | 0.35 | 0.0535 | 0.0559 | 5.211 | 0.1700 | 0.7258 | 0 |
| a037 | 5.5 | eb1168 | 4.8 | 0.37 | 0.0535 | 0.1339 | 4.564 | 0.2097 | 0.5026 | 1 |
| b029 | 6.0 | eb1168 | 6.4 | 0.29 | 0.1294 | 0.0033 | 6.046 | 0.0857 | 1.4151 | 0 |
| b031 | 6.0 | eb1168 | 6.4 | 0.31 | 0.1294 | 0.0392 | 5.224 | 0.1775 | 0.8761 | 0 |
| b033 | 6.0 | eb1168 | 6.4 | 0.33 | 0.1294 | 0.1318 | 4.487 | 0.2255 | 0.6786 | 1 |
| c50  | 5.0 | eb1168 | 4.8 | 0.34 | 0.0111 | 0.0479 | 4.586 | 0.1858 | 0.7477 | 0 |
| c65  | 6.5 | eb1168 | 6.4 | 0.33 | 0.2132 | 0.0572 | 6.153 | 0.1298 | 1.8528 | 0 |
| d030 | 5.5 | eb1168 | 6.4 | 0.30 | 0.0535 | 0.0284 | 4.374 | 0.2029 | 0.7529 | 0 |
| d036 | 6.0 | eb1168 | 4.8 | 0.36 | 0.1294 | 0.0511 | 6.391 | 0.1365 | 1.8973 | 0 |
| e0482 | 5.5 | eb0482 | 4.8 | 0.36 | 0.0111 | 0.0913 | 5.356 | 0.1870 | 1.0057 | 0 |
| e154 | 5.5 | eb154 | 4.8 | 0.36 | 0.0761 | 0.0919 | 4.651 | 0.1959 | 0.6407 | 1 |
| f725 | 7.25 | eb1168 | 3.2 | 0.27 | 0.3080 | 0.0000 | 16.994 | 0.0000 | 12.1053 | 0 |

(arm A = a03x E₀ ladder through sub-basin 1; arm B = b0xx through
sub-basin 2; arm C = c50/c65 twin-fail edge controls; arm D =
d030/d036 τ-crosses, one per sub-basin; arm E = well bracket at the
a-arm center, `eb0482`/`eb154` under the documented pairing hatch;
arm F = f725 the standing point = the broken-landing continuity anchor. a037 is also arm A's twin-gated
top; twin n̄ spans 6.7→4.6 (A) and 6.0→4.5 (B), so the MD 4.07 crossing
sits inside each arm for any bias in the +0.2…+3.2 bracket.)

**MD-level acceptance (bias-free, pre-registered):** a cell *lands* iff
MD n₁_solv ∈ [0.19, 0.30] ∧ MD n̄_det ∈ [3.77, 4.37] (4.07 ± 0.3).
W₁ reported, not gating (single-seed N = 500 scatter ≈ 0.1–0.15); KE
(midHot/deepKE/χ²_med) reported, direction-only vs the twin; supp soft
(RQ3 mixture).

**Pre-registered predictions:**

- **GR-P1 (oracle):** the 14 twin rows above reproduce string-exact
  from the committed scan CSV; the scorer reproduces the pooled-battery
  recorded row within 0.002. Fails ⇒ session invalid.
- **GR-P2 (basin transfer — the headline test):** ≥ 1 cell in
  arms A ∪ B lands the MD acceptance. This requires the
  near-quantitative n₁ transfer AND a small-end n̄ bias (≈ +0.5…+1.1)
  to co-occur at the same E₀ — genuinely falsifiable; a large-end bias
  (+3) pushes the landing below every twin-n₁-viable E₀ and refutes
  the basin at MD level (itself the atlas result).
- **GR-P3 (bias direction):** MD n̄ < twin n̄ at every paired cell,
  with Δn̄ ∈ [0.2, 3.2] at ≥ 70 % of cells.
- **GR-P4 (n₁ transfer):** |MD n₁_solv − twin| ≤ 0.05 at every
  non-control cell with ≥ 100 scored solvated ions.
- **GR-P5 (controls fail as predicted):** c50 misses low (n̄ under the
  MD floor), c65 misses on n₁ (< 0.19); f725 reproduces the broken
  landing (n₁_solv < 0.05, n̄ > 10, trap ≥ 0.17 = twin 0.308 minus the
  max floor bias).
- **GR-P6 (policy evidence):** trap ordering eb0482 < eb1168 < e154 at
  the arm-E triple; bound/marginal decomposition recorded per cell for
  the G4 retained-policy call.

**Adjudications taken with this design (user, 2026-07-27):**
retained-policy — `exclude_all_coupled` stays the *interim* convention
for the ring with the decomposition recorded; the **final policy call
moves to G4** on the ring's measured class sizes (the scan's
(v_c, well)-sensitivity already points at the over-dissipation-artifact
branch). Axis B — **declared folded** into this ring (arms A/B/D are
its MD E₀/τ content at the corrected geometry) + the committed §3.5c
free-surface maps; any standalone Axis-B study is deferred past G4 and
re-opens only if the ring's MD curves disagree with the twin maps
beyond the §14.4 box.

**Execution record (2026-07-27/28):** all 14 cells complete (detached
relaunch after two external kills; oracles passed on every launch).
**GR-P1..P6 all land — the basin is MD-CONFIRMED**: gated cells
**a037 (n̄ 4.097) / b031 / d030 / e154**; the MD basin spans both τ
values at v5.5 and the deep well; c50 shows it nearly reaches v_c 5.0
(fails on n₁ by 0.005, not on the predicted n̄ clause); f725 measures
finc1v725 broken at the corrected geometry (trap 0.437 = 0.278 bound
+ 0.159 marginal, n₁ 0). n̄ bias 0.24–2.66 He small-end (14/14
twin-hot), n₁ transfer ≤ 0.036, well-trap lever ≈ 0.8/eV. Marginal
class ≈ 0 at the basin vs 0.159 at the standing chord (retained-policy
evidence, MD-grade). KE: no gated cell holds midHot and deepKE at once
(RQ11 successor question at v_c 5.5–6.0). Full record: findings
"G3 Step 3", D0 §14.5; table CSV committed. **G4 is live** (successor
point + retained policy + ledger re-issue; user adjudications).

### 3.5e G4 Step 1 — the fine ridge sweep (DESIGNED + user-approved 2026-07-28; **EXECUTED 2026-07-28** — Blocks 0/1/3 + a user-requested ladder arm; successor candidate **h405**, W₁ floor ≈ 0.67, deficit localized to the size-distribution *shape*; findings "G4 Step 1")

The G3 ring confirmed the basin but left the successor point
underdetermined: **no gated cell holds both KE axes**, and the two clean
sub-basins were *adjacent corners of a grid that was never scanned
between them*. This stage refines the cheap axes before any structural
conclusion is drawn (the standing anti-refutation rule: never declare a
surface exhausted while a post-hoc-rescorable axis is coarser than the
incumbent's known basin width).

**The resolution gap, stated.** Step 2 scanned v_c at step 0.5
({5.0, 5.5, 6.0, 6.5, …}) and τ at ratio ≈ 1.33 ({3.2, 4.8, 6.4, 9.6}).
The fully-unflagged basin is (v_c 5.5, τ 4.8) and (v_c 6.0, τ 6.4) — a
**diagonal pair with no scanned cell between them**; c50 missing the MD
gate by 0.005 on n₁ says the same from the low-v_c side. E₀ (step 0.01)
was already fine and is not the gap.

**Objective (user, 2026-07-28):** joint score with the **KE tension as
the thing to break** — the n₁/n̄ hard gate stays, surviving cells are
ranked on W₁ + midHot + deepKE.

#### Block 0 — the twin's ranking authority (zero cost, **gating**)

§14.4 licenses the twin for n₁ (quantitative), n̄ (biased), trap (floor)
— **not** for W₁ ("bias-loaded") or KE ("direction-only"), which are
exactly the ranking observables. The G3 ring bought the fix for free:
**14 paired (twin, MD) cells** never read as a transfer measurement.

- Per observable (W₁, midHot, deepKE, χ²_med): Spearman ρ between twin
  and MD across the 14 cells, **with CI**.
- **Permission gate (pre-registered): ρ ≥ 0.7 ⇒ the twin may *rank* on
  that observable in Block 1; ρ < 0.7 ⇒ gate-only/reported and MD must
  rank that axis.** Used as permission only — never as a correction to
  a twin W₁/KE value.
- **n̄ bias model** Δn̄ = f(residence / trap) fitted from the same 14
  cells, replacing the crude [+0.3, +3] bracket that forced the Step-2
  gate to n̄ ∈ [4.4, 7.1]. This tightens the twin gate toward the MD
  band itself and sharpens the scan more than the grid refinement does.
- Verify the **deep-bin KE seed-SD** from the pooled battery: §1.2
  records "0.0603 ± 0.0033", which does not reconcile with the pooled
  deepKE 0.631. The score's deepKE normalizer is provisional until this
  is resolved.
- Output: `h2b_g4_transfer.csv`; D0 §14.4 extended.
- Risk: n = 14 over a narrow span ⇒ weak ρ estimates. If Block 0
  licenses nothing, Block 1 degrades to a gate-only scan and the MD
  budget shifts toward the two-stage pattern.

#### Block 1 — the fine ridge scan (twin, zero MD)

Hypothesis under test: **the two clean sub-basins are opposite corners
of one connected diagonal ridge in (v_c, τ), and the KE tension is
minimized between them** — where nothing was ever evaluated.

| axis | Step 2 | fine | cost |
|---|---|---|---|
| v_c | {5.0, 5.5, 6.0, 6.5} step 0.5 | {5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5} step 0.25 | 9 new chord integrations (3 v_c × 3 E_bind) ≈ 15 min |
| τ | {3.2, 4.8, 6.4, 9.6} ratio 1.33 | {4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8} step 0.4 | free (post-processing) |
| E₀ | [0.17, 0.52] step 0.01 | [0.26, 0.42] step **0.005** | free |
| E_bind | {0.048, 0.1168, 0.154} | **unchanged** | — (no new value ⇒ no new pairing exception) |

≈ 5540 cells, same order as the 6480 of Step 2. τ = 6.8 carries the
existing > 6.55 calibration-class flag. Everything else stays at the
§3.5d pins (corrected geometry, rq4graded, p = 1, p_tail −1, co-moving
shed, Landau 0.58, budget 2.70).

Gate: n₁_solv ∈ [0.19, 0.30] ∧ **bias-corrected** n̄ ∈ [3.77, 4.37].
Rank: the joint score below, on the Block-0-licensed observables only.

**The old grid is an exact sub-lattice of the new one** — which is what
makes G4-P1 a strong oracle.

#### Block 2 — p_tail, conditional (twin, zero MD)

Fires **only** if no Block-1 cell holds midHot ∈ [0.85, 1.15] ∧ deepKE
≥ 0.6 simultaneously. Then p_tail ∈ {−0.5, −1, −1.5, −2} on the top two
chords (6–8 integrations).

Rationale: the D0 verdict "p_tail is not a second lever" was measured at
**v_c 7.25**. The cap has moved to 5.5–6.0 — *into* the 5–9 Å/ps
mid-band that NB-RQ11-12 identified as the deep-KE lever — so the tail
exponent now shapes γ exactly where the deep-KE sensitivity lives and
the old verdict does not automatically transfer.

Either outcome is a result: it breaks the tension ⇒ the deep-KE lever is
the tail *exponent*, not v_c (a direct NB-RQ11-12 answer); it does not
⇒ the tension is **structural across the whole drag surface**, which
localizes it to the mechanism (ladder / pickup / ε) — a far stronger
statement than "no better cell was found".

#### Block 3 — MD finalists (6 × N = 1000, one fresh shared seed, CRN-paired)

MD's role has changed from breadth ("does the basin exist", N = 500 was
right) to **discrimination** among finalists on W₁/KE differences of
0.1–0.3 against N = 500 single-seed scatter of ±0.1–0.15. Budget is
therefore spent on precision, not coverage (user decision 2026-07-28).

- 3 new ridge cells (or 2 + one p_tail cell if Block 2 fired),
- **a037 and b031 replicated at N = 1000** — mandatory: without them the
  new N = 1000 numbers are not comparable to the ring's N = 500 numbers,
  and they measure the N = 500 scatter the current ranking sits inside,
- 1 pre-registered off-ridge fail control.

Acceptance: the same bias-free MD gate; ranking by the joint score.

#### Block 4 — G4 proper

5 × N = 1000 pooled battery (§4cc pattern) at the winner, then the user
adjudications: successor point, retained-policy final call
(`exclude_all_coupled` interim; ring evidence = marginal class ≈ 0 at
the basin vs 0.159 at the standing chord), D0 §17 ledger re-issue with
the retired **S** rows struck, figure surface regenerated.

#### The joint score (pre-registered)

Hard gate `n₁ ∈ [0.19, 0.30] ∧ n̄ ∈ [3.77, 4.37]`, then the **Pareto
front** over (W₁, |ln midHot|, |ln deepKE|), with one scalar for
reporting and tie-break:

```
S = W₁/0.571  +  |ln midHot|/ln(1.15)  +  |ln deepKE|/ln(1.15)
```

Equal weights, declared **convention, not truth**; the deepKE
normalizer is provisional pending the Block-0 seed-SD check. χ²_med is
reported, not scored (scale varies 44–1322 across the ring and it is
largely redundant with midHot + deepKE). Evaluated on the existing
numbers:

| point | W₁ term | midHot term | deepKE term | **S** |
|---|---|---|---|---|
| standing `finc1v725` (old geometry, pooled N = 5000) | 1.00 | 0.08 | 3.29 | **4.37** |
| a037 (N = 500) | 1.45 | 0.44 | 6.04 | **7.94** |
| e154 | 1.73 | 0.29 | 6.04 | **8.06** |
| b031 | 1.86 | 4.96 | 2.35 | **9.17** |
| d030 | 2.05 | 1.51 | 9.36 | **12.9** |

Reading: the incumbent still wins, **and its own score is dominated by
the deepKE deficit** — RQ11 was the largest single defect before the
geometry correction; the correction worsened it at a037 and improved it
at b031. The sweep therefore has one explicit target: **beat 4.37.**

#### Pre-registered predictions

- **G4-P1 (oracle):** the fine scan reproduces all 24 Step-2 gated cells
  **bit-exact** on the shared sub-lattice, and the transfer table
  re-reads the 14 committed ring rows string-exact. Fails ⇒ session
  invalid, no number read.
- **G4-P2 (ridge connectivity):** the gated set at the standing well
  forms a *connected* path from (5.5, 4.8) to (6.0, 6.4). Falsifiable —
  they may be two islands.
- **G4-P3 (tension):** ≥ 1 ridge cell holds midHot ∈ [0.85, 1.15] **and**
  deepKE ≥ 0.6 at once. Its negation is the Block-2 trigger.
- **G4-P4 (score):** ≥ 1 ridge cell reaches S < 4.37 — the corrected
  geometry gets at least as close to experiment as the wrong geometry
  did. Legitimate comparison: W₁ and both KE ratios are distances to the
  *same* experimental reference.
- **G4-P5 (MD):** the Block-3 winner's MD S beats a037's re-measured
  N = 1000 S.

#### Risks and boundaries

1. **Thin transfer sample** (14 cells, narrow span) — Block 0 may
  license nothing; fallback above.
2. **Twin τ-feedback:** in MD, τ feeds back on trap through the mass
  mechanism, so the free-surface decomposition is itself a twin
  approximation (§3.5c caveat). A 0.4-ps τ step may be finer than the
  twin's fidelity in τ; the d030/d036 τ-crosses partially measure this,
  and if it fails τ reverts to a 0.8-ps step.
3. **Refinement may buy only noise-level gains** — the N = 1000
  finalists are the guard.
4. Nothing adopts before the G4 adjudication; **`finc1v725` stands**
  throughout. Tier-0 form authority, b, the committed scorer, checkpoint
  schema, RNG draw order and the constants table are untouched.

#### Execution record (2026-07-28)

Blocks 0, 1 and 3 all executed the same day; Block 2 did **not** fire.
Total spend: zero-MD for Blocks 0–1, **15 cells × N = 1000** for Block 3
(6 finalists + a 9-cell ladder arm the user requested after seeing the
first six). All oracles passed on every launch. Full record: findings
"G4 Step 1" (three sections) and D0 §14.4/§14.5 + the §2/§3/§4 knob
entries.

- **Block 0** — W₁ ρ +0.824 and midHot ρ +0.996 **licensed**, deepKE
  ρ +0.327 **not** ⇒ the deep-KE axis is not twin-scannable at any
  resolution. n̄ bias model resid SD 0.28 He. Two plan §1.2 provenance
  corrections recorded (the "W₁ SD 0.04" is the SEM of the pooled mean;
  the deep-bin "0.0603 ± 0.0033" is not the deepKE ratio).
- **Block 1** — G4-P1 oracle bit-exact on 408 shared cells; 574/5544
  gate (172 clean); **G4-P2 CONFIRMED**, the ridge is one connected
  diagonal band; G4-P4 not met on the licensed axes once the transfer is
  applied as a regression; **G4-P3 proved NOT EVALUABLE at twin level**
  — a pre-registration defect found by G4's own Block 0, since its
  deciding observable is the unlicensed one.
- **Block 3** — G4F-P3: the W₁ deficit is **not** sample size (a037
  N = 500 → 1000 moved −0.034 against a per-seed SD of 0.095).
  **G4F-P4: the mid-vs-deep KE tension is BROKEN** (deepKE 0.438 →
  0.614 across f2/f1/f3 at near-identical n̄/n₁; f3 holds midHot 1.023
  *and* deepKE 0.614) — so **Block 2 (p_tail) was never triggered** and
  the G3-ring "no cell holds both axes" reading was a coarse-E₀
  artifact. G4F-P2 refuted by 0.001–0.003 on n₁ (design gap: n₁ was
  gated uncorrected; its ±0.02 transfer scatter decides the band edge).
- **Ladder arm** — E₀ slope measured (per +0.005 eV: W₁ −0.02…−0.04,
  n₁ +0.004…+0.007, n̄ −0.11…−0.13, midHot −0.02…−0.03, monotone in all
  three τ arms, gate window ~0.01 eV). **W₁ floor ≈ 0.67**, measured
  twice independently (arm asymptote 0.669; Block-0 regression intercept
  0.671). **Root cause: n₁ and n̄ are not simultaneously matchable** —
  n₁ = 0.243 needs E₀ ≈ 0.435 where n̄ ≈ 3.3, i.e. 0.77 He below
  experiment. The residual is a *shape* mismatch in the detected size
  distribution, localizing it to the mechanism (ladder / pickup / ε),
  not the drag surface. τ ordering is the twin's reverse
  (4.4 < 4.8 < 5.2 on W₁). Deep-KE ceiling ≈ 0.656 at v_c 5.25, paid in
  midHot.
- **Successor candidate: h405** (v_c 5.5, τ 4.4, E₀ 0.405) — gated,
  best gated S 1.559, better than a037 on every axis, supp 0.183 vs the
  standing 0.187, marginal class 0.002.

**Still open (the G4 adjudications, user):** successor point, retained
policy, ledger re-issue, and whether the W₁ floor is recorded as the
corrected geometry's honest residual or chased into a new mechanism
axis. Verification before re-baselining: a pooled 5 × N = 1000 battery
at the chosen point.

### 3.5f G4 Step 2 — h405 pooled battery + W₁ residual anatomy (DESIGNED + user-approved 2026-07-28; **EXECUTED 2026-07-28** — GV-P1/P2 CONFIRMED, GV-P3 REFUTED (pooled S 1.734 > 1.683); anatomy vs the frozen fingerprints: **no knob matches → honest-residual branch**, F1-partial + trapped-tail caveats recorded; low-n KE axis found mid-stage and pre-registered hard; findings "G4 Step 2" ×3)

G4 Step 1 left one residual axis — the W₁ floor ≈ 0.67, measured twice
independently and root-caused: **n₁ and n̄ are not simultaneously
matchable on the (v_c, τ, E₀) surface**, a *shape* mismatch in the
detected size distribution that localizes to the mechanism (ladder
shape, pickup, per-shed ε), not the drag surface. This stage does two
things in parallel: **verify h405** (the designed pooled battery), and
**diagnose the floor** before deciding whether it is chased into a new
mechanism axis or recorded as the corrected geometry's honest residual
(user decision 2026-07-28: *diagnose first, then decide* — depth
option B, decomposition + pre-registered fingerprints, no MD probes).

The stage's product is a **decision input**, not an adoption: either a
named mechanism knob with a matched fingerprint (→ that knob is the
next designed axis), or an evidence-backed "record 0.67 as the honest
residual" verdict. Atlas stance holds throughout: nothing adopts
mid-stage, `finc1v725` stands until the G4 adjudications fire.

#### Block V — the pooled battery at h405 (MD: 5 × N = 1000, fresh seeds)

§4cc pattern at h405 (v_c 5.5, τ 4.4, E₀ 0.405; all other knobs at the
§3.5d pins, corrected geometry). Oracles before any number is read:
the pooled-battery scorer-drift oracle, and the committed h405 Block-3
row rescored **string-exact**.

Pre-registered:

- **GV-P1 (gate):** h405 gates on the pooled read
  (n₁ ∈ [0.19, 0.30] ∧ n̄ ∈ [3.77, 4.37]).
- **GV-P2 (floor consistency):** pooled W₁ ∈ [0.64, 0.78]; per-seed SD
  comparable to the standing battery's 0.095.
- **GV-P3 (score):** pooled S beats a037's re-measured 1.683 and lands
  near the Block-3 single-seed 1.559.

GV-P1/P2 failure ⇒ h405 is not seed-robust ⇒ back to the ridge with
the battery evidence; no adjudication fires. All three pass ⇒ the three
near-mechanical adjudications become fireable (user calls, recorded at
G4 proper): **successor = h405**; **retained policy finalized as
`exclude_all_coupled`** (ridge marginal class 0.001–0.003 ≈ empty);
**ledger re-issue** (D0 §17, retired S rows struck) + figure-surface
regeneration at the successor.

#### Block D — W₁ residual anatomy (zero MD, scorer-only)

W₁ on a 1-D histogram is exactly ∫|F_mod − F_exp|, so it decomposes
into **signed per-bin CDF-gap contributions**. Compute the
decomposition (per-bin sign + magnitude + cumulative share) for:

1. **pooled h405** (Block V output — the object of interest),
2. **pooled `finc1v725`** (standing N = 5000, wrong geometry): what the
   0.579 landing looked like bin-wise — the shape the corrected
   geometry lost,
3. **the τ 4.4 E₀ arm** (h405 → h410 → h415, committed CRN-paired
   N = 1000 rows): the finite-difference **bin signature of the E₀
   lever** — the measured template of what a knob signature looks
   like, and the direct read of which bins are **E₀-inaccessible**
   (those carry the floor).

#### Block F — pre-registered mechanism fingerprints (doc work; FROZEN before any Block-D number is read)

For each candidate knob — **ladder shape (c1 / rq4graded grading; RQ4),
pickup (Poisson rate), per-shed ε (RQ2)** — derive from the mechanism
(analytically / from the cascade structure, not by eyeball) its
expected signature: sign on n₁, sign on n̄, and the bin-pattern of the
CDF-gap change. The measured deficit constrains the derivation target
in advance: at matched n₁ (0.243 at E₀ ≈ 0.435), n̄ runs ≈ 0.77 He
light — the model is missing solvated/high-n weight *relative to*
n = 1, so each fingerprint must state whether the knob can add
solvated-shoulder mass **without paying it back out of n₁**.

**Match rule (pre-registered):** a knob is a live candidate iff its
signature moves, with the correct sign, the bins carrying ≥ 70 % of the
residual W₁, without a wrong-signed prediction on n₁ or n̄.

**Decision table (pre-registered):**

- exactly one knob matches → it becomes the next designed axis (own
  plan section, twin-authority caveat carried: twin W₁ cannot
  discriminate ≲ 0.05, so that axis budgets MD from the start);
- multiple match → cheapest-first discrimination (analytic/twin
  fingerprint sharpening) before any MD;
- none match → the W₁ floor is recorded as the corrected geometry's
  **honest residual** (D0 §17 + findings), the G4 adjudications
  complete, and the tier proceeds to the open items (RQ3/RQ5 reads,
  margin-3 Å pin I88, D2b A/B remainder, then Tier-3).

#### Machinery and cost

Behind `[PROCEED TO IMPLEMENTATION]`: a battery generator for h405
(adapting the §4cc / `gen_tier2atlas_g4finals.py` pattern), one new
decomposition scorer under `scripts/post_processing/` (scorer-only;
reads committed run dirs + the experimental reference; emits a
committed CSV), focused tests. Block F is doc work in the findings/plan
before the scorer's output is read. MD spend: 5 × N = 1000 (~45–90 min
wall); everything else zero-MD.

#### The low-n KE axis — PRE-REGISTERED as a hard axis (user adjudication 2026-07-28)

Found mid-stage (user observation, verified on the committed artifacts
before the Block-V battery was read): the corrected-geometry ridge
**halves the n = 1 mean KE** (h405 0.637 / f3 0.668 eV vs the incumbent's
1.034 against ref mean/median/mode 1.302/1.128/0.891) and takes ~20 % off
n = 2 (0.551/0.580 vs ref 0.706), while n ≥ 3 lands (h405 n = 3: 0.471 vs
0.496). The scored surface was blind to it: n = 1 KE enters only the
median-anchored, *reported-not-scored* χ²_med (which did register:
242 → 350/426), and midHot's n = 2–8 geo-mean dilutes the bad n = 2
behind six good bins. Mechanism: the corrected geometry's ~9 → ~35 Å
birth depth makes the fastest, most-stripped fragments pay the full
transit-drag toll — the incumbent's n = 1 KE match was partly a
shallow-birth artifact (the same wrong-geometry compensation as its
W₁ 0.571). Note: p_tail (Block 2) was dismissed on the *deep*-KE axis;
the n = 1 KE axis is where the high-v tail exponent actually acts, so
that dismissal does not carry over.

**Pre-registration (governs every future ridge/successor ranking):**

- Observables: **KE₁** = sim mean detected KE in bin n = 1 vs the
  reference **median** 1.128 eV (the I77 n = 1 anchor convention);
  **KE₂** = sim mean in bin n = 2 vs the reference **mean** 0.706 eV.
- Score terms, same form as the §3.5e S convention:
  `|ln(KE₁/1.128)|/ln(1.15) + |ln(KE₂/0.706)|/ln(1.15)` — **added to S**
  in any future ranking.
- **Hard gate:** band half-widths set at 2× the per-seed SD measured by
  the Block-V battery (the deepKE-normalizer precedent: axis frozen now,
  norm frozen from the first measured scatter, recorded in the findings
  before any ranking uses it).
- The Block-V scorer reports KE₁/KE₂ (mean/median/mode, per member +
  pooled) from this stage on.

#### Risks and boundaries

1. **Fingerprint degeneracy** — the three knobs may produce
   sign-identical bin patterns at this resolution; the decision table's
   "multiple match" branch absorbs this without extra MD.
2. **Battery non-robustness** — GV-P2's [0.64, 0.78] band is wide
   enough that a failure is informative, not noise.
3. Nothing adopts mid-stage; Tier-0 form authority, b, the committed
   scorer, checkpoint schema, RNG draw order and the constants table
   are untouched.

#### Execution record (2026-07-28)

Built and executed same-day behind the user's trigger; build plan
`TIER2_G4_STEP2_IMPLEMENTATION_PLAN.md` (7 tasks, 8 commits). Every
oracle passed on every launch/read. Mid-stage events: (a) first battery
launch lost 4/5 members to disk-full — recovered by the user-approved
trajectory strip (194 scored dirs, committed scorers re-verified
bit-exact) and a same-seed relaunch; (b) the **low-n KE finding**
(user-observed, verified on committed artifacts, pre-registered above
before the battery was read).

- **Block V:** GV-P1 CONFIRMED (pooled n₁ 0.2094 / n̄ 3.877; 5/5 members
  gate), GV-P2 CONFIRMED (pooled W₁ 0.7653 — the floor's third
  independent measurement; per-seed SD 0.0337, 2.8× tighter than the
  standing 0.0954), **GV-P3 REFUTED** (pooled S 1.734 > a037's 1.683;
  the Block-3 1.559 was a favorable draw) ⇒ no adjudication fires
  automatically. Low-n KE measured at 5× stats: KE₁ ratio 0.568 ± 0.003,
  KE₂ 0.778 ± 0.002 — physics, not noise.
- **Block D:** residual = **both-ends deficit** (n = 1 −0.101; tail
  n ≥ 14 −0.055) vs core n = 2–12 excess +0.156 — under-dispersion.
  E₀-lever measured n = 1–4-only ⇒ bins 7–20 (≈ 45 % of residual)
  E₀-inaccessible. **Provenance correction:** "experimental 0.243/4.07"
  in the Step-1 narrative are the incumbent's sim values; the true
  solvated reference is n₁ 0.3103 / n̄ 4.889 (gates are
  incumbent-anchored conventions; W₁/floor unaffected; the
  not-simultaneously-matchable conclusion strengthens).
- **Match rule:** F2/F3 refuted as pre-registered; F1 partial (low side
  addressable, tail Σ-locked) ⇒ **no match → honest-residual branch**,
  with the trapped-tail shadow (missing tail mass ≈ the excluded
  trap_bound 0.0768, measured slow-and-high-n) recorded as the
  quantitative stake on the open retained-policy adjudication.

**Open after this stage (all user):** successor point (GV-P3 refuted —
h405 unverified at battery level), retained policy (now carrying the
trapped-tail stake), ledger re-issue, and the W₁-floor + low-n-KE
residual pair (recorded honest; candidate levers noted: retained
policy, p_tail on the *low-n KE* axis, p occupancy exponent).

### 3.5g The low-n KE retro-scan (DESIGNED 2026-07-29; zero MD — the in-surface freedom must be measured before any new axis; **EXECUTED 2026-07-29** — O1–O3 passed; R5: in-gate bound +0.19 vs required +0.31 eV → **in-surface freedom EXHAUSTED, p_tail motivated**; R4: the n = 1 channel is **bifurcated** — deep-born cells produce zero n = 1, every standing-chord n = 1 is shallow-born at 1.04–1.09 eV = the reference peak scale → the birth-law/mixture lever is the second live candidate; findings "§3.5g")

**Motivation (user direction, 2026-07-29).** Before p_tail (or any new
knob) is opened on the low-n KE axis, the freedom the *existing*
(v_c, τ, E₀) surface has on KE₁/KE₂ must be measured at the resolution
already on disk — the refine-before-refute rule (the §6.6 lesson, now a
standing working principle). Every committed corrected-geometry run
retains its `detection.npz` (the trajectory strip removed only
`ion.npz`/`relaxation.npz`), so KE₁/KE₂ are retro-scorable over the full
committed surface with zero MD. KE₁ is a high-SNR observable (Block-V
per-seed SD ≈ 0.003 eV at N = 1000, expected ≈ 0.02–0.03 eV per N = 500
cell from the n = 1 occupancy) — single-seed cells rank cleanly, unlike
W₁.

**Convention amendment (user adjudication 2026-07-29, supersedes the
§3.5f anchor for scoring; the old anchor stays flagged-not-overwritten
per the §1.2 precedent).** The experimental n = 1 KED shows a wide
*upper* tail read as a second process (not the drag-cascade channel this
model owns); the n = 1 target is therefore the **peak, ≈ 1.00 eV** —
between the reference mode 0.891 and median 1.128, deliberately not
chasing the tail-inflated mean 1.302. From this section on: **KE₁ is
scored against 1.00 eV** (score term `|ln(KE₁/1.00)|/ln 1.15`); ratios
vs the old median anchor stay reported for continuity. **KE₂ keeps the
reference mean 0.706 eV** (no second-process signature adjudicated on
n = 2). Current distance at h405 (pooled): KE₁ 0.641 → deficit
−0.36 eV against 1.00.

**Scope (52 rows, all committed runs, no new MD):**

1. the **h405 battery** (5 × N = 1000 + pooled) — the anchor;
2. the **G4 finals + ladder arm** (15 × N = 1000, CRN-paired seed
   20260729) — the τ arms 4.4/4.8/5.2 at v_c 5.5, the E₀ ladders, v525
   (v_c 5.25) and b031 (v_c 6.0);
3. the **G3 ring** (14 × N = 500, seed 20260728) — wider (v_c, τ)
   coverage incl. c50/c65 edge cells, the e-well bracket, and **f725**
   (standing chord at the corrected ensemble);
4. the **geometry grid** (11 × N = 500, fixed geometry, standing chord)
   plus the **incumbent battery** (5 × N = 1000 + pooled, standing
   geometry) — KE₁/KE₂ vs (R, depth) at *fixed* drag parameters: the
   depth-provenance read (does any shallow/small-droplet cell produce
   reference-scale-fast n = 1?).

Excluded: the lq runs (form authority is Tier-0's; the E_bind confound
family) and the pre-correction N = 500 `finc*` probes (superseded by the
incumbent battery at 10× the statistics).

**Oracles (§1.4 — all before any new number is read):**

- **O1 scorer-drift:** the standing pooled battery reproduces its
  recorded observable row (the standard `check_oracle`).
- **O2 committed-KE:** the five h405 battery members + pooled, rescored
  through this scorer's own code path, reproduce the committed
  `atlas_g4step2_battery.csv` KE₁/KE₂ columns (mean/med/mode) to
  1e-9 relative — the KE extraction is bit-compatible with Block V.
- **O3 incumbent-KE:** the incumbent pooled battery lands on the
  recorded KE₁ 1.034 / KE₂ 0.754 (3-decimal comparison; the Block-V
  scorer's recorded constants).

**Pre-registered readings (computed only after O1–O3 pass):**

- **R1 (τ):** ∂KE₁/∂τ and ∂KE₂/∂τ at v_c 5.5 — across the three
  CRN-paired arms, both raw and at matched n̄ (the gated chord
  h405/h375/h345). This is the axis the user named; it has never been
  read.
- **R2 (v_c):** KE₁/KE₂ at v_c 5.0/5.25/5.5/6.0/6.5/7.25 from
  c50/v525/(5.5-cells)/b031/c65/f725 — the cap sets where the
  uncalibrated tail begins, so this is the strongest in-surface
  candidate.
- **R3 (E₀):** within-arm slopes per +0.005 eV (CRN-paired ladders);
  the h405↔f3 pair already suggests ≈ −0.03 eV KE₁ per +0.01 eV E₀
  (composition effect — to be confirmed on the full arms).
- **R4 (geometry at fixed chord):** KE₁/KE₂ vs (R, depth_mean) over the
  11 grid cells + the incumbent battery + f725 — decomposes the deficit
  into transit-toll vs mixture-weights and settles the n = 1 birth-depth
  provenance question (2026-07-29 discussion).
- **R5 (reachability verdict, the section's deliverable):** linearize
  the measured slopes about h405 and ask whether **any** in-gate move
  (n₁ ∈ [0.19, 0.30], n̄ ∈ [3.77, 4.37], midHot ∈ [0.85, 1.15] soft)
  reaches KE₁ ≥ 0.95 (1.00 − 1σ of the linearization). **If yes** —
  the identified chord becomes the next MD arm (own registration,
  before any p_tail talk). **If no** — the in-surface freedom is
  recorded exhausted on the low-n KE axis, and the p_tail axis is
  motivated *with this record as its evidence* (the same structure by
  which the E₀ arm established the W₁ floor).

**Boundaries.** Pure scorer, runs nothing, mutates nothing; nothing
adopts from this section (successor/retained/ledger stay the open G4
adjudications); Tier-0 form authority, the committed scorers, checkpoint
schema, RNG draw order and the constants table untouched. N = 500 KE₁
carries ≈ 0.02–0.03 eV single-seed scatter — ring/grid rows are read at
that resolution, never over it.

### 3.5h The p_tail ring — the low-n KE axis opened on the drag tail (DESIGNED + user-approved 2026-07-29; **EXECUTED 2026-07-29** — placement CONFIRMED in-bracket (pt15 0.917 / pt20 1.151), but **PT-P3 KILL FIRED** (midHot 1.9–4.4 at every KE₁ ≥ 0.95 cell; deepKE to 12.5×; PT-P4 sign-inverted — residence/relaxation channel dominates; PT-P5 refuted) → **single-knob p_tail CLOSED**; measured shape suggestion: band-limited tail v_c2 ≈ 9.5; findings "§3.5h")

**Adjudications recorded (user, 2026-07-29, post-§3.5g discussion):**

1. **The shallow-birth/mixture lever is REJECTED on physics.** The central
   parent-Boltzmann birth law *is* the physics (solvated I₂ sits at the
   droplet center); the §3.5g R4 shallow cells were an artifact of the old
   wrong geometry, not a channel reality offers. The R4 bifurcation is
   re-read diagnostically: if real births are central and the real n = 1
   peak is ≈ 1.0 eV, then real deep-born n = 1 ions ARE fast — **the
   current tail over-drags**, and R4 becomes evidence *for* weak high-v
   drag, not for a mixture lever.
2. **The p_tail axis is approved** as the next designed axis (this
   section), with the §3.5g record + the composition ceiling as its
   motivation.

**The composition ceiling (exploratory read on the committed h405 battery,
2026-07-29; to be re-scored under this section's oracle chain):** the
pooled n = 1 KE distribution is a needle — max **0.734 eV**, p99 0.724,
**zero of 1556 ions above 1.0 eV** (n = 2: max 0.643). No re-selection /
composition knob (λ₀, s_eff, f_ret, ladder, retained policy, occupancy p)
can raise a mean to 1.0 over a population capped at 0.73 — only levers
that *add kinetic energy to an ion* can, and with KER experimental and
the geometry physically fixed, the force law above the TDDFT band is the
unique remaining energy-side lever.

**Guard widening (the §3.5h scoped code change):** the config-load guard
pinned `capped_cubic` `p_tail` to the Step-1c set {0, −1} with the
explicit message "any other exponent needs a fresh adjudication" — this
is that adjudication. Widened to **−4 ≤ p_tail ≤ 0** (dissipative
softening only; positive exponents stay refused). The drag module itself
was already general (`(v/v_c)**p_tail`); no physics code changed. Tests
updated as adjudication-driven (documented in-file), plus an explicit
p_tail = −2 closed-form exponent-law regression the generic identities
cannot catch. `capped_linear_quadratic` keeps {0, −1}.

**Placement (1-D anchored transit integral, zero MD — placement
authority only, no ranking):** anchored to reproduce the measured h405
toll (start 14.33 Å/ps = 1.35 eV; exit 9.87 = 0.641 eV at p_tail = −1
⇒ L_eff = 26.1 Å at ρ̂ = 1 under the extraction mass 202.95; the
calibration absorbs the mass — p* is bracket-invariant over m ∈
[127, 250]):

| p_tail | −1.0 | −1.5 | −1.75 | −2.0 | −2.5 | −3.0 |
|---|---|---|---|---|---|---|
| predicted KE₁ [eV] | 0.641 (anchor) | 0.889 | 0.984 | 1.062 | 1.173 | 1.241 |

**KE₁ = 1.00 eV at p_tail ≈ −1.80.** Known omissions: the heating
back-reaction (weaker tail → less drag work → less E_int → shallower
cascade → the *population* ending at n = 1 changes) and the ensemble
spread — MD measures the truth.

**Cells (3 × N = 1000 + 1 conditional; seed 20260729 = the finals seed,
so every comparison against the committed h405 row is CRN-paired):**

| cell | p_tail | role |
|---|---|---|
| pt15 | −1.5 | below-target bracket |
| pt20 | −2.0 | at/above-target bracket (placement 1.06) |
| pt30 | −3.0 | far probe: overshoot + midHot-damage measurement |

Everything else = the committed `g4fh405` pins verbatim (corrected
geometry, v_c 5.5 / τ 4.4 / E₀ 0.405, exclude_all_coupled). The
conditional 4th cell (E₀ recenter along the measured §3.5g slopes)
fires only per PT-P4 below.

**Oracles (§1.4, before any MD and by --dry-run):** each cell's cfg
diffs against the committed `g4fh405` `cfg.json` in **exactly**
`{"drag_coefficients"}`, and inside the bundle only `p_tail` differs
(b, v_c, stamps bit-identical); the scorer re-runs the standing
scorer-drift oracle and reproduces the committed finals h405 row to 4
decimals before any new number is read.

**Pre-registered predictions:**

- **PT-P1 (placement band):** KE₁ rises monotonically with |p_tail|;
  per-cell bands = placement ± 0.10 eV: pt15 ∈ [0.79, 0.99], pt20 ∈
  [0.96, 1.16], pt30 ∈ [1.14, 1.34]. Outside-band ⇒ the 1-D placement
  model is wrong in a way worth diagnosing before any further tail work.
- **PT-P2:** KE₂ rises toward 0.706 with the same ordering; KE₃ rises
  mildly (≤ +0.10 at pt20).
- **PT-P3 (the kill criterion, v525-informed):** if **every** cell with
  KE₁ ≥ 0.95 shows midHot > 1.15, the tail exponent alone fails the
  axis exactly as the v_c cap move did, and the axis moves to a joint
  (p_tail × v_c/τ/E₀) re-tune or the honest-residual branch — no
  further single-knob tail cells.
- **PT-P4 (back-reaction):** n₁ falls and n̄ rises with |p_tail| (less
  drag work → less E_int → shallower cascades). If a cell with KE₁ ∈
  [0.95, 1.16] lands n₁ < 0.19, the conditional recenter cell fires:
  E₀ raised along the measured §3.5g arm slopes (+0.004…+0.007 n₁ per
  +0.005 eV), one cell, same seed, registered by this clause.
- **PT-P5 (orthogonality):** deep bins are untouched — deep-n ions live
  near/below v_c, where the law is byte-identical; deepKE within ±0.05
  of the h405 battery's 0.50 (per-seed SD 0.02). A violation means the
  tail reaches deeper into the cascade than the exit-speed picture
  says, and the §3.5g slope arithmetic must be redone with p_tail live.
- **Tier-0 legitimacy by construction:** the in-band branch (v ≤ 5.5)
  is byte-identical to the Tier-0-locked pure cubic; the TDDFT traces
  end at 5.58 / 3.36 Å/ps, so no trace constrains any of these cells'
  tail. No Tier-0 re-run is owed.

**Scoring:** the full standing surface + the §3.5g KE terms (KE₁ vs the
1.00 anchor, KE₂ vs 0.706). Success shape: a cell holding KE₁ ∈
[0.95, 1.16], midHot ≤ 1.15, gate (n₁/n̄) intact — such a cell becomes
the candidate for a pooled verification battery (the GV precedent)
before any successor talk. **Nothing adopts from this section;
`finc1v725` stands; the G4 adjudications stay open.**

**MD spend:** 3 × N = 1000 (~45–60 min at concurrency 3) + at most one
conditional cell. Disk checked: ≈ 0.8 GB/cell against 17 GB free.

### 3.5i The drag state coupling s(n) — design ADJUDICATED + probe REGISTERED (2026-07-29; own document)

The low-n KE axis continues off the γ(v) surface: after the §3.5h kill
(velocity-only softening un-damps the cascade) and the user's rejection
of the band-limited tail as overfitting, the adjudicated route is the
**drag–shell state coupling** `γ(v,d,n) = g(d)·s(n)·γ_form(v)` — the
drag learns the ion has stripped. Full design, dimensional analysis,
nine closed adjudications (OQ-A..I) and the registered 3-cell probe
(sa22/sa30/sa44: ρ_shell across its Bounded range at the h405 pins,
seed 20260729; needle-break signature SC-P2; midHot kill SC-P4) live in
**`TIER2_DRAG_STATE_COUPLING_DESIGN.md`** — the axis's own document.
Build slices S1–S4 wait for `[PROCEED TO IMPLEMENTATION]` (fresh
session). Nothing adopts; `finc1v725` stands; G4 adjudications open.

### 3.6 Open questions this axis must answer or explicitly defer

- **Is E_bind R-dependent? — ANSWERED at G0 (2026-07-26): yes, but
  negligibly.** The well splits into an R-independent local
  (snowball/electrostriction) term and a Born far-field term; only the
  latter carries R, and helium's ε = 1.0572 makes it small — **≤ 0.009 eV
  (7.6 %) over the whole span**, i.e. ≈ ⅛ of the smallest step §6.7
  item 2 measured as a live lever (trap +0.058 per the 0.0482 → 0.1168
  step ⇒ ≈ 0.85/eV, hence Δtrap ≈ 0.006 here). All cells run at
  0.1168 eV; an MD bracket cell is **armed as a rider** on the trapped
  read below. Full derivation D0 §9.1; G0 item 3.
- **Does the trapped channel explode at large R?** Trapped rises
  0.000 → 0.184 across the battery's R quintiles (R̄ 19.6 → 35.0 Å); at
  R ≈ 50–68 Å it may dominate, in which case the *detected* ensemble is a
  small-R-biased subset of the source ensemble. A G1 read, not a separate
  study — and now also the **trigger for the E_bind rider** (trapped
  > ~0.4 or a handover-guard trip). Secondary consequence to watch: a
  large trapped fraction shrinks the *detected* count a cell delivers, so
  the §3.1 resolution pre-registration must be re-checked against the
  realized detected N at R3 and R = 68.3 Å before those rows are read.
- **Do T5/T6 stay in the codebase once inert? — yes.** Keep the enums
  (delivered, tested, needed if the geometry question reopens) and mark
  them inert-at-production in the ledger — guard-live arms with recorded
  provenance, not dead code under rule 2. Reinforced by the control
  rejection: no physically realizable control can decompose them anyway
  (one shared `rho_he_ratio` surface, §3.1c).

## 4. Axis A addendum — D2b: provenance audit + salvageability of the legacy sampling laws

The production sampling laws are theory-laden legacy ports. This
sub-axis audits them and measures what they change.

### 4.1 Provenance audit — **EXECUTED 2026-07-26**

Delivered as D0 §15 (the authoritative chapter: which sampler ran in
which production, the bulk N → R convention, the parent-document anchor,
and the inherited ⟨N⟩ = 2000 pin) — it re-ordered this sub-axis, putting
the ⟨N⟩ pin ahead of δ and E_solv.

### 4.2 Distribution-level A/B (zero MD)

Compare the *sampled distributions themselves* under variants, in this
priority order (re-ordered 2026-07-26 by the §4.1 audit):

1. **the ⟨N⟩ pin: 2000 (standing) vs 12794 (source correlation) vs
   16.3k (pickup-weighted)** — the dominant, previously unlisted axis;
2. raw vs post_pickup vs analytic prior at a *matched* ⟨N⟩ (so the
   shape difference is not confounded with the pin) — **partly
   pre-answered 2026-07-26**: the parent document's quoted R range
   (34–68.3 Å) matches the **raw** ln-normal's q05/q95 (34.9 / 68.7 Å),
   not post-pickup's (37.7 / 74.1 Å), so the parent ensemble is raw;
3. position law vs alternatives (center-pinned, uniform-in-volume at
   margins 0/3/6 Å, Boltzmann) — noting that Boltzmann's own output is
   R-dependent (§4.1 item 3), so it must be compared at each R;
4. δ perturbations (0.40 / 0.625 / 0.80);
5. E_solv 14 vs 30 meV — **demoted**: it only enters the legacy pickup
   chain, which production does not run.

Pure sampling plots; no propagation. Every N → R conversion uses the
bulk form (§4.1 item 2).

### 4.3 Observable-level sensitivity — via grid re-weighting (zero new MD), then ≤ 2 confirmations — **EXECUTED 2026-07-27** (oracle INADMISSIBLE on the full standing density — support hole below R = 26.6 Å; post-hoc in-support oracle 7/7, `nearest` adopted; corrected-ensemble forecast recorded in findings "D2b §4.3" / D0 §15.7; the below-support R ≈ 20 Å × L3 cell is the designed first confirmation candidate)

The §3 grid is the transfer function {(R, r/R) → observable vector}.
Any candidate sampling law is evaluated by **re-weighting the grid
cells with its (R, r) density** — predicting the ensemble effect of
each theory variant without new MD (the coarse 3 × 3 support is a
stated limitation; interpolation convention recorded in the findings
doc). Only variants that matter after re-weighting get a full-sampler
MD confirmation (N = 500 each; **expected after the §4.1 audit: the
⟨N⟩ 12794 source-condition pin first, then legacy post_pickup vs
analytic prior — ≤ 2 runs**; E_solv is no longer a candidate).

**Support note.** The grid spans R = 26.6 → 49.4 Å, so the
source-condition pin is inside the support and re-weightable; the
parent's q95 tail (to 68.3 Å) extrapolated under the 9-cell design and is
**closed by the two funded R = 68.3 Å cells** (§3.1).

**Mandatory re-weighting oracle (added 2026-07-26).** Before any candidate
law is re-weighted, the method must reproduce a *known* ensemble:
re-weight the grid with the standing point's own sampled (R, depth)
density and check it recovers the pooled N = 5000 battery row (n₁_solv
0.243, W₁ 0.571, supp 0.187, n̄ 4.07, trapped 0.067) within the stated
interpolation error. If the reconstruction of a known ensemble fails, no
reconstruction of an unknown one is admissible — same discipline as the
scorer oracle (§1.4), zero MD. The interpolation convention (cell
weights, nearest-neighbour vs bilinear in (R, depth)) is fixed by this
oracle and then recorded in the findings doc.

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
| 2b | Axis A grid (§3) = geometry-correction stage **G1**; 3 × 3 + 2 **funded** R = 68.3 Å cells (controls rejected, §3.1c) | 11 × 500 | trigger (gen + report scripts + the `droplet_size_sampler_mode` selector) |
| 3 | D2b audit (**§4.1 DONE 2026-07-26**) + A/B + re-weighting (**§4.3 DONE 2026-07-27, zero MD**) (+ ≤ 2 confirmations; first candidate = the below-support R ≈ 20 Å × L3 cell) | ≤ 2 × 500 | trigger (sampling-variant runs) |
| 4 | Axis B E₀/τ curves | 8 × 500 | trigger |
| 5 | D4 Step 2 twin sweeps | zero | none (scratchpad twin, §4ee precedent) |
| 6 | D4 Step 3 spot-checks + E_bind zero-MD swaps/twin, then conditional MD confirms | ~0–7 × 500 | trigger (incl. enum build; swaps are scratchpad) |
| 7 | synthesis: merge all results into `TIER2_PARAMETER_INFLUENCE.md`, close GAP markers | zero | none |
| G3 | geometry re-arbitration: Step 1 landmarks + Step 2 twin scan (zero MD) + Step 3 MD ring (**DONE** 2026-07-27/28) | 14 × 500 | trigger (per step) |
| G4-1 | fine ridge sweep (§3.5e): Blocks 0/1 zero-MD + Block 3 finalists **and** the ladder arm (**DONE** 2026-07-28) | 15 × 1000 | trigger (per block) |

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
