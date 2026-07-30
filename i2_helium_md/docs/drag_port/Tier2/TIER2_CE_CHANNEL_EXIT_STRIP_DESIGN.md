# Tier-2 (C) Design — CE Channel Mixture (B) + Depth-Graded Exit Stripping (A)

**Status: DESIGN ADJUDICATED (2026-07-29, user: "I follow your
recommendations") — physics definition only. No code until
`[PROCEED TO IMPLEMENTATION]`. OQ-A..K ALL CLOSED: OQ-K was closed
2026-07-29 by user confirmation of the power↔intensity mapping
(300 mW ↔ 1.47×10¹⁴ / 600 mW ↔ 2.94×10¹⁴ W/cm² — the 600 mW KED
reference pairs with the ~20 % Q3 covariance share). P1 + P2
EXECUTED 2026-07-29 (zero MD; findings "(C) pre-steps P1 + P2"):
the strip prior box is PINNED (a = 2 [1.5, 2.5], j₀ ∈ [1.5, 2],
w_j ∈ [0.5, 1]; a = 1 disfavored) and the width priors are
MEASURED (σ_Q2 0.31 / σ_Q3 0.55 eV upper bounds; E_single ≈ 0.53) —
the §3.5 table below carries the updated values. The m/q-127-row
ask is postponed (user) until the method proves successful. **P3
EXECUTED 2026-07-29** (weight vector entering registration
(0.30, 0.50, 0.20); C-full placement forecast + two registered
tensions — mode two-lump, n₁ 0.18–0.23; findings "(C) pre-step
P3"). **PROBE REGISTRATION FROZEN 2026-07-29 (user: "I approve the
bands")** — the CP-1..8 magnitudes and frozen inputs are in §8; the
probe cells in §9. **Outcome pre-commitments PC-1..5
ADJUDICATED 2026-07-29 (§11); `[PROCEED TO IMPLEMENTATION]` ISSUED;
the probe was BUILT and EXECUTED 2026-07-29 (4 × N = 1000, seed
20260729; findings "(C) probe EXECUTED", `atlas_ce_probe.csv`) —
REGISTRATION FAILED: CP-1..4 FAIL, kills CP-5/6/7 FIRED on C-full,
CP-6 at both coupling-arm ends. Measured decomposition: (B) alone
places KE₁ = 0.920 with the needle broken (its slow single channel
traps 0.465); (A) as specified over-tolls (ε-dominated 0.5–0.7 eV
at knock counts p90 = 20 — the ε box is inconsistent with its own
0.2 eV anchor) and FEEDS the suppressed gate (supp 0.55 in A-only:
Σ(n) collapses under an untouched E_int). The §11 PC-3 consequences
are TRIGGERED — (A)-v1 stops, the current f_int wiring closes, the
mixture survives as measured; next-step adjudication is with the
user.** Nothing here is
adopted; `finc1v725` stands; h405 candidacy and the G4
adjudications stay open.**

**ADJUDICATED CLOSED 2026-07-30 (user): the PC-3 consequences are
ratified, ALL candidate follow-ups are DECLINED (ε → ~0 gated strip /
E_int co-strip / trap-class pricing), and the (C) line is SHELVED as
a measured boundary — the covariance-anchored mixture cannot populate
n = 1 at the reference level (Q3 feeder ceiling ≈ 0.11 of the scored
ensemble; suppressed-conversion bound n₁ ≈ 0.22 < 0.31). Re-opening
requires new external evidence, not a re-parameterization. Records:
findings "(C) adjudication CLOSED", D0 §20, log 2026-07-30. The
positive result stands: KE₁ position is source physics (B-only
0.920, needle broken) — not drag or mass-mechanism physics.**

Precedent and template: `TIER2_DRAG_STATE_COUPLING_DESIGN.md` (the
s(n) axis — designed, probed, refuted, stopped). This document designs
the **(C) combination** adjudicated in findings §3.5j: (B) supplies
the source-KER spread (multi-channel Coulomb explosion), (A) converts
fast + dressed exiters into small-at-existing-KE at the droplet
boundary. Entry point context: plan §3.5l.

---

## 1. Motivation — the evidence chain (pointers, not re-derivation)

- **(A) alone is necessary-not-sufficient** — counterfactual bracket
  ceiling KE₁ ≈ 0.71–0.77 eV vs the 1.00 eV anchor; 48.7 % of the
  experimental n = 1 KED sits above the (A) cap (findings §3.5j).
- **(B) alone cannot place KE₁** — under full proportional E₀ wiring
  the fast ions exit through the suppressed-bare door with KE₁ at
  baseline (bud411, findings §3.5k); suppression is **E_int-driven**
  (∂supp/∂E₀ ≈ +1.7 /eV).
- **The source-side lever is measured alive:** S_k = 0.389 eV/eV
  (§3.5k, BP-KILL not fired; D0 §19).
- **The needle requires the mixture:** KE₁ SD 0.034–0.047 at any
  single budget (BP-P4) vs experimental σ 0.697 eV.
- **Channels are paper-anchored:** Hatherly et al, J. Phys. B: At.
  Mol. Opt. Phys. 27 (1994) 2993–3003 (consulted 2026-07-29; not
  kept in-repo — key numbers extracted here) — gas-phase finite-pulse CE,
  peak KER a *channel-independent* fraction of point-charge E_C
  (0.75/0.65·E_C at 90/200 fs), channels carry intrinsic width
  (total-KER FWHM 2.6 eV for I⁺+I⁺, 5.8 eV for I⁺+I²⁺ at 200 fs).
  The group's own gas spectrum is reproduced at 0.8·E_C with Q = 2 /
  Q = 3 from 2.666 Å (RQ7/RQ8 NBs); the Abel export confirms peaks
  2.26/4.11 in the repo's own calibration.
- **The Q3 weight is measured, not fitted:** I²⁺Heₙ ↔ I⁺He event
  covariance ~7 % (1.47×10¹⁴) / ~20 % (2.94×10¹⁴ W/cm²) (thesis
  Fig. 6.2; RQ8 NB).
- **Per-ion retention is stochastic at fixed kinematics:** the
  retention-twin anchor — a third of n = 1 fragments at low intensity
  have a fully-stripped same-KER twin (RQ8 NB follow-up round). This
  is the (A) closure's target distribution.
- **The suppressed fate class is scaffolding:** all 334 suppressed
  fragments in the h405 battery carry mechanical shell n = 21 and
  cross the surface at 10.66–11.69 Å/ps — the model already re-labels
  its fast dressed exiters at scoring level; (A) replaces that rule
  with physics (findings §3.5j structural identification; D0 §17).
- **New kill axis from read (i):** the experimental bare KED holds
  ≤ 2 % below 1.0 eV — slow-bare overpopulation is scoreable
  (findings "§3.5l reads (i)+(iii)").

## 2. Physical picture

At the pump pulse each I₂ is multiply ionized into one of a small set
of CE channels. Which channel fires is a per-molecule stochastic event
tied to the local intensity, not to the droplet: {single ionization
(I₂⁺ dissociation, slow), Q = 2 (I⁺+I⁺), Q = 3 (I⁺+I²⁺, fast)}. Every
channel releases a *fraction* f ≈ 0.8 of its point-charge Coulomb
energy from R₀ = 2.666 Å (finite-pulse sequential ionization,
Hatherly) with intrinsic per-channel width. Inside the droplet the ion
is always dressed (bulk-refill argument, §3.5j — impulse stripping
self-heals where there is bath), so the drag/mass machinery is
untouched in-flight. The one place net stripping exists is the
**outbound surface crossing**: the shell He are knocked out
stochastically, outer shell eagerly, innermost 1–3 shadowed/protected
(pure binding-threshold → all-bare; pure Poisson-count → mid-tail
parking; both measured wrong, §3.5j). A fast Q3-derived ion therefore
exits under-dressed or bare *at its existing fast KE* — the n = 1 KED
becomes the drag-tolled image of the fast KER channels, and the
experimental bare peak (mean 3.706 eV) is the under-dressed Q3
transit.

## 3. Governing forms and strict dimensional analysis

### 3.1 (B) — the channel mixture

Per molecule m (one CE pair), draw independently of geometry:

    c_m ~ Categorical(w_single, w_Q2, w_Q3),      Σ w_c = 1
    E_m ~ TruncNormal(E_c, σ_c;  E_m > E_min)     [eV per I⁺ fragment]

with channel means (point-charge E_C per I⁺ from R₀, times the shared
fraction f):

    E_ref  = 14.39964548 / (2 · R₀_Å) = 2.7006 eV      (R₀ = 2.666 Å)
    E_Q2   = f · 2 · E_ref · (1/2) · 2 = f · 2 · E_ref / 2 · 2  — see below
    written explicitly:
      Q2 (q₁q₂ = 1): total f·E_C = f·5.4012 eV → per I⁺  f·2.7006  (= 2.16 at f = 0.8)
      Q3 (q₁q₂ = 2): total f·2·E_C            → per I⁺  f·5.4012  (= 4.32 at f = 0.8)
      single:        E_single  Bounded prior [0.3, 0.8] eV (Abel gas slow
                     structure 0.27/0.49/0.70; not Coulombic — no f)

f is **one shared Bounded parameter** across channels (Hatherly:
fraction channel-independent at fixed pulse), prior 0.80, range
[0.65, 0.90]. Width priors from Hatherly 200 fs, per-I⁺ σ:
σ_Q2 ≈ 0.55 eV, σ_Q3 ≈ 1.23 eV (FWHM/2.355 of the per-ion half of
2.6/5.8 eV); to be sharpened by the zero-MD Abel width read (§9 P2).

**Wiring (the emulation route, recommended — OQ-A):** the pair's
Coulomb drive is scaled per molecule,

    U_m(r) = s_m · q₁q₂ · 14.39964548 / r,   s_m = E_m / E_ref   [dimensionless]

with both fragments kept at unit charge. For point charges with equal
masses this is **kinematically exact** for the I⁺ fragment: the force
∝ q₁q₂/r² enters only through the product, so scale s at q₁q₂ = 1
reproduces the trajectory of the true charge pair. What it does NOT
reproduce is the partner's identity: in a Q3 event the partner is an
I²⁺ (detected on a different m/q row), while the emulation produces a
second fast "I⁺". **Partner masking** is therefore part of the design:
in Q3-labeled molecules one fragment is tagged `q3_partner` and
excluded from every I⁺-scored observable (histogram, KED, fates).
Without the mask the ensemble gains a spurious fast-I⁺ population of
weight ≈ w_Q3/2 — a designed-in scoring error. (§10 OQ-I.)

The standing scalar budget `coulomb_available_eV = 2.7006` (scale 1.0)
**retires with the mixture**: there is no single budget under (C);
the per-ion stamp becomes E_m. The known-high-25 % RQ7 fact is thereby
discharged structurally, not by a uniform re-anchor (bud226k measured
what a uniform re-anchor does).

### 3.2 Per-channel E₀ coupling

The S2 onset generalizes from the scalar `f_int` to a per-channel
coupling (T6 p-law preserved verbatim):

    E_int,i(0) = f_int,c(m) · E_m · (Σ(n₀,i)/Σ(n*))^p     [eV]

Measured calibration pair (findings §3.5k): full proportional coupling
on Q3 (f_int 0.15 → E₀ 0.6165) parks the fast ions in suppressed-bare
with KE₁ at baseline; pinned E₀ (f_int 0.0985) keeps them alive but
overheats the cascade at *uniform* Q3 weight. Under the mixture the
overheating dilutes with w_Q3; the Q3 coupling is the design's main
*shape* lever on the fast branch:

    f_int,Q2 = 0.15 (standing; h405 pins unchanged on the dominant channel)
    f_int,Q3 ∈ [0.0985, 0.15]   — Bounded by the measured bracket
    f_int,single = 0.15 default (weakly identified; slow channel)

All f_int,c dimensionless.

### 3.3 (A) — the depth-graded exit strip

At an **outbound surface crossing** with exit speed v_x [Å/ps] and
mechanical shell n_x, each rung j = 1…n_x (j = 1 innermost) is
independently knocked out with probability

    P_knock(j) = P₀(v_x) · G(j)
    P₀(v)      = min(1, (v / v_strip)^a)                 [dimensionless]
    G(j)       = 1 / (1 + exp(−(j − j₀) / w_j))          [dimensionless]

- `v_strip` **Sourced** at ≈ 9.9 Å/ps (the equal-mass max-transfer
  full-strip threshold vs the production `rq4graded` rungs — the CF-3
  convention, findings §3.5j); `a` Bounded (steepness of the velocity
  gate), prior a ≈ 2 (ram ∝ v²).
- `j₀` (protection depth, rungs) Bounded prior [1, 3]; `w_j` (grading
  width, rungs) Bounded prior ≈ 1. G encodes "outer eager, inner
  shadowed" — the element §3.5j measured as *required* (sharp
  threshold → all-bare at 0.484 eV; Poisson-count → W₁ 1.5–1.6
  mid-tail parking).
- Independent Bernoulli draws at fixed (v_x, n_x) implement the
  **retention-twin stochasticity**: same-kinematics ions leave with
  different survivor counts. Anchor: at the n = 1-producing exit
  speeds (~10.2 Å/ps), P(all stripped) must be same-order as
  P(one survivor) (a third of n = 1 fragments have a bare twin at
  1.47×10¹⁴ — RQ8 NB). This anchors (j₀, w_j, a) *without touching
  the n = 1 KED* (§7).
- **n = 0 outcomes allowed** (OQ-L adjudicated, §3.5j) — feeds RQ3.

**Energetics and ledger (the §3.5j-round OQ-K — NOT this document's
§10 OQ-K condition mapping — adjudicated):** stripping is re-labeled
surface-crossing drag work with an explicit term — no new energy
source or sink. Per knocked rung j the ion loses

    ΔE_j = D₀(j) + ε_carry        [eV]

(D₀(j) the rung binding from the resolved ladder; ε_carry Bounded
[0, ~0.05] eV per He, the kinetic carry-off), applied to the outbound
velocity at the crossing (direction preserved); the mass drops by
m_He per knock through the existing variable-mass machinery (post-jump
m⁺ convention, SQ1–SQ3). The 5-term invariant is unchanged in total:
E_strip = Σ ΔE_j is booked as a labeled sub-term of the dissipation
channel at the crossing. Total expected scale ≈ 0.2 eV ≈ the measured
post-exit slowdown (§3.5j). The fate of E_int on a full strip is
§10 OQ-H.

### 3.4 What the suppression rule becomes

The `suppressed` scoring class (self-unbound at detection, scored as
bare-equivalent with mechanical n = 21) is the scaffolding this design
retires (D0 §17): its population — fast dressed exiters — now passes
through the strip at the crossing and lands in physical n = 0…small-n
bins at existing KE. The scoring rule itself is NOT removed in v1 (no
behavior change off-mode; and residual self-unbound ions must still be
handled) — the *prediction* is that its occupancy → ~0 under (C)-on
(§8 CP-1).

### 3.5 Dimensional table (every new parameter)

| parameter | units | class | prior / range (P1/P2-updated 2026-07-29) | anchor |
|---|---|---|---|---|
| w_single, w_Q2, w_Q3 | — (Σ=1) | Bounded | w_Q3 ≈ 0.20 @ 600 mW (OQ-K confirmed); w_single open | I²⁺ covariance, power series, gas slow band (§7) |
| f | — | Bounded | 0.80 [0.65, 0.90] | own gas fit + Hatherly + Abel peaks |
| E_single | eV | Bounded | 0.53 (σ ≈ 0.42) in [0.3, 0.8] | P2 measured (Abel gas slow component) |
| σ_Q2, σ_Q3 | eV | Bounded | **0.31 / 0.55** (measured upper bounds; Hatherly 0.55/1.23 superseded) | P2 measured (§9); window-stable ± 0.02 |
| f_int,Q3 | — | Bounded | [0.0985, 0.15] | bud411k/bud411 bracket |
| f_int,single | — | Free (weak) | 0.15 | — |
| v_strip | Å/ps | Sourced | 9.9 | equal-mass max transfer vs rq4graded |
| a | — | Bounded | **2 [1.5, 2.5]** (a = 1 disfavored: W₁ parking + slow-bare) | P1 measured (§9); ram scaling |
| j₀ | rungs | Bounded | **[1.5, 2]** (was [1, 3]; ends killed by twin ratio) | P1 measured; retention-twin |
| w_j | rungs | Bounded | **[0.5, 1]** | P1 measured; retention-twin |
| ε_carry | eV/He | Bounded | [0, 0.05] | 0.2 eV total-scale constraint |

Honest count: ~11 new degrees of freedom, the largest single freedom
addition in Tier 2. Mitigation: 7 carry external anchors independent
of the arbitration observables; the kill set (§8) is pre-registered;
the circularity guard (§7) is absolute.

**Dof accounting (added 2026-07-29, user challenge "11 seems
large").** The raw count is the honest price of replacing *two*
idealizations at once — a delta-function source (1 knob: the scalar
budget) becomes a physical 3-channel source (8 of the new dof:
2 weights, f, 2 widths, E_single, 2 couplings), and a scoring rule
(the suppressed relabeling — itself a hidden dof-equivalent hard-coded
at scoring level) becomes a physical boundary process (4 dof: a, j₀,
w_j, ε_carry). The *effective* freedom at probe time is far smaller
than the raw count:

| tier | parameters | count | status at probe launch |
|---|---|---|---|
| frozen at external anchors | f, w_single, w_Q2, w_Q3, E_single, σ_Q2, σ_Q3, f_int,Q2, f_int,single, ε_carry | 8-ish | fixed values, anchor-derived (gas fit, covariance, P2 read, 0.2 eV bound); behave as measured inputs with error bars, not fit knobs |
| pinned by the zero-MD P1 step | a, j₀, w_j | 3 | calibrated detection-only on committed checkpoints against the retention-twin + bare/solvated split, then frozen |
| scanned by the probe | f_int,Q3 | 1 | the coupling arm (bracket ends) |

So the probe launches with **one** scanned knob, scores ~12
pre-registered observables (n̄, n₁, supp, trap, W₁, midHot, deepKE,
KE₁, KE₁ SD, KE₂, bare weight, slow-bare), and carries three kills a
wrong survival shape cannot tune around (W₁ parking, slow-bare,
midHot). dof ≪ observables with bands frozen in advance — the
CALIBRATION_MAP identifiability stance is preserved: Bounded-with-
anchor is not Free.

## 4. Tier-0 / Tier-1a legitimacy

- Both mechanisms live behind their own enums, default **off**;
  off-mode is **byte-identical** (zero new RNG draws when off — the
  streams are untouched, not merely unused; the s(n) S1–S4 precedent).
- The drag module is untouched: no new coupling, still mass-agnostic,
  γ(v) convention unchanged. Tier-0 (fixed mass, deterministic,
  channels off) is unaffected by construction.
- The strip requires shell bookkeeping → v1 is **biphasic-only**
  (config-load guard: `exit_strip_mode != "off"` requires
  `mass_scenario = "biphasic"`; same for `ce_channel_mode`).
- The §6.5/§6.6 mass-pairing stance is unchanged (mid-window m_eff
  argument does not involve the source or the crossing).

## 5. Architecture placement (minimal violence)

- **sampling/** owns the channel draw (new focused module; dedicated,
  separately-seeded named RNG stream appended after all existing
  streams — the forbidden-list "random-number draw order" is honored:
  existing draws are bit-identical, channels-off consumes nothing).
- **ion_initial_state** consumes the per-molecule (c_m, E_m): pair
  Coulomb scale s_m, per-ion onset E_int(0) (§3.2). `E_coulomb_scale`
  (global) remains and multiplies s_m (oracle: global × per-ion).
- **the ion driver / relaxation stage** applies the strip at the
  outbound gate crossing it already detects (the drag spatial gate);
  mass jump through the existing SQ machinery; ledger sub-term booked
  at the crossing.
- **detection scoring**: partner mask (§3.1); suppressed rule kept,
  occupancy predicted → 0.
- **Checkpoint**: per-ion channel label, sampled E_m, strip count are
  new per-ion state → **schema v8, a forbidden-list change requiring
  explicit user approval** (§10 OQ-E). No constants in the checkpoint;
  weights/f/etc. live in cfg.json as always.

Enum surface sketch (names final at implementation; rule-2 carries
recorded at build): `ce_channel_mode {"off","sampled"}`,
`ce_channel_weights`, `ce_fraction_f`, `ce_channel_sigma_eV`,
`ce_single_ker_eV`, `ce_q3_partner_mask (default True)`,
`internal_energy_partition_fraction_per_channel`,
`exit_strip_mode {"off","depth_graded"}`, `exit_strip_v_ref`,
`exit_strip_exponent`, `exit_strip_protect_j0`,
`exit_strip_width_rungs`, `exit_strip_carry_eV`.

## 6. Trade-offs and lost physics (working-method disclosure)

- **The Q3 partner is not propagated as I²⁺**: kinematics of the
  scored I⁺ exact (point-Coulomb argument §3.1), but the partner's own
  dressing/detection physics (the I²⁺Heₙ covariance row) is *not*
  modeled — the model cannot predict the I²⁺ observables it is
  anchored on. Accepted for v1; un-refusing `highly_charged_iodine`
  is the upgrade path (OQ-A).
- **No Q4 (I⁺+I³⁺)**: the bare upper-lump reach above 4.32 eV (read
  (i): mode 4.758) may contain it (per-I⁺ ≈ 6.5 eV at f = 0.8, at the
  detector edge). Refused for v1, recorded as the bare-lump rider
  (OQ-G) — bare is a soft anchor only.
- **Channel draw ⊥ droplet geometry**: no intensity–position
  correlation inside the ensemble (focal-volume physics not modeled).
- **Strip is instantaneous at the crossing**: no finite surface
  thickness; multiple crossings (trapped orbits) are OQ-J.
- **Single channel emulated as a low-KER budget**, not I₂⁺ PECs — the
  refused `single_charge_ionization_allowed` branch stays refused in
  v1.
- The neutral driver, checkpoint RNG order, constants table: untouched.

## 7. Anchoring and the circularity guard

Channel weights and strip parameters are **never fitted on the n = 1
KED** (the 2026-07-29 meta-discussion guard — the n = 1 KED is the
arbitration observable). Anchors, in priority order:

1. **w_Q3**: the I²⁺ covariance shares (~7 %/~20 %) at the mapped
   condition (600 mW ↔ 2.94×10¹⁴ provisional — §10 OQ-K), cross-checked
   by the read-(iii) power ordering (10.2/16.8/23.1 %).
2. **w_single vs w_Q2**: the gas slow band + the low-intensity Q2
   partner structure (RQ8 NB); the m/q-127-row covariance ask would
   close it cleanly (open rider).
3. **f, σ_c**: the own-experiment gas fit (0.8) + Abel peak/width
   reads; Hatherly brackets.
4. **(A) parameters (a, j₀, w_j)**: the retention-twin distribution
   (same-KER dressing splits) + the bare-vs-solvated *weight* split —
   explicitly not the n = 1 KED shape.
5. **Bare KED**: soft anchor (lump positions; the ≤ 2 % slow-bare
   bound is a kill, not a fit target).

## 8. Sign-level predictions — REGISTRATION FROZEN (user approval 2026-07-29)

**Frozen bands (the P3 proposal, user: "I approve the bands"). These
are the registered magnitudes; the sign-level statements below stand
as their rationale. Scored on the C-full cell unless stated; kills
are kills wherever they fire.**

| prediction | frozen band | basis |
|---|---|---|
| CP-1 | supp ≤ 0.05 | conversion by construction (from 0.183) |
| CP-2 | KE₁ SD ∈ [0.30, 0.70] | P3 forecast 0.51–0.55 |
| CP-3 | KE₁ **mean** ∈ [0.85, 1.20] AND above-1.15 share ∈ [0.35, 0.55] | P3 forecast 0.92–1.04 / 0.43–0.50; mean-level per Tension 1 (the two-lump mode is a named discriminator, not a gate) |
| CP-4 | n₁_solv ∈ [0.18, 0.26] | P3 forecast band + one seed-SD; NOT the 0.31 reference (Tension 2) |
| CP-5 (kill) | ensemble slow-bare ≤ 0.01 (bare-row share below 1 eV ≤ 5 %) | experimental 1.9 % |
| CP-6 (kill) | midHot ∈ [0.85, 1.15] | h405 neighborhood (§4cc seed-robust band) |
| CP-7 (kill) | W₁_solv ≤ 0.80 | baseline 0.709 + one seed-SD |
| CP-8 (soft) | bare KED two-lump qualitative | read (i); Q4 rider acknowledged |

**Frozen inputs:** weights (w_single, w_Q2, w_Q3) = (0.30, 0.50,
0.20) (OQ-C closed; OQ-K confirmed); f = 0.80; σ_Q2/σ_Q3 =
0.31/0.55 eV, E_single = 0.53 (P2); strip box a = 2 [1.5, 2.5],
j₀ ∈ [1.5, 2], w_j ∈ [0.5, 1] (P1); f_int,Q3 = the scanned arm
(0.0985 / 0.15); v_strip 9.9 Sourced; ε_carry [0, 0.05].
Cells: C-full / A-only / B-only / coupling arm, N = 1000, CRN,
`drag_state_coupling = "off"`. Circularity guard absolute: nothing
above was fitted on the n = 1 KED.

Original sign-level statements (rationale record):

- **CP-1**: suppressed-class occupancy → ~0 under (C)-on (the
  scaffolding retires; replaced by physical n = 0…2 at fast KE).
- **CP-2**: the KE₁ needle breaks — SD > 0.08 (vs 0.034–0.047 at any
  single budget) and the n = 1 KED upper-half weight rises toward the
  48.7 % reference via the selection-amplified Q3 share.
- **CP-3**: KE₁ mode moves into a band around the 1.00 eV anchor
  (band frozen at registration; the zero-MD mixture arithmetic puts
  the blend at ~0.85–1.2).
- **CP-4**: n₁_solv rises toward 0.31 (bracket precedent 0.353;
  mixture gate-survivors alone ~0.17 — the delta is (A)'s attribution).
- **CP-5 (kill)**: slow-bare stays ≤ few % below 1.0 eV (read-(i)
  axis; the band frozen at registration).
- **CP-6 (kill)**: midHot/deepKE stay in the h405 neighborhood — the
  w_Q3-diluted, f_int-coupled fast branch must NOT re-fire the §3.5h
  cascade (midHot band frozen at registration).
- **CP-7 (kill)**: no mid-tail parking — W₁ must not degrade above
  the h405 baseline (the §3.5j Poisson-mode failure signature).
- **CP-8**: bare KED two-lump structure qualitatively reproduced
  (soft; Q4 rider acknowledged).

## 9. THE PROBE — EXECUTED 2026-07-29 (built behind `[PROCEED TO IMPLEMENTATION]`; cells cfull/aonly/bonly/cq3hi via `gen_tier2atlas_ce_probe.py`, scored by `tier2atlas_ce_probe_table.py` — REGISTRATION FAILED, see the status header + findings "(C) probe EXECUTED")

Zero-MD pre-steps (before any launch):

- **P1 — strip-form prior calibration, detection-only — EXECUTED
  2026-07-29** (committed instrument `tier2atlas_ce_strip_prior.py`,
  artifact `atlas_ce_strip_p1.csv`; findings "(C) pre-steps P1 + P2"):
  the §3.5j counterfactual re-run with the actual P₀·G form, all five
  oracles passed (incl. exact reproduction of the §3.5j crossing band
  and the sharp/bracket variant rows). Prior box pinned — see the
  §3.5 table; the retention-twin ratio is a-independent at the twin
  speed and selects (j₀, w_j) alone.
- **P2 — Abel gas width read — EXECUTED 2026-07-29** (committed
  instrument `tier2atlas_ce_gas_widths.py`; peak oracle 2.26/4.11
  passed): σ_Q2 0.31 / σ_Q3 0.55 eV (measured upper bounds — narrower
  than the Hatherly priors), E_single ≈ 0.53 (σ 0.42). The n = 1 KED
  width must therefore come mostly from channel separation + strip
  selection, not intrinsic channel width.
- **P3 — weight-composition forecast table — EXECUTED 2026-07-29**
  (committed instrument `tier2atlas_ce_p3_forecast.py`, artifact
  `atlas_ce_p3_forecast.csv`; findings "(C) pre-step P3"): weight
  vector entering registration **(0.30, 0.50, 0.20)**; C-full
  forecast KE₁ mean 0.92–1.04 / SD 0.51–0.55 / above-1.15 43–50 %
  (ref 48.7 %) / Q3 carries 53–61 % of the n = 1 bin; **two
  registered tensions**: the blend is two-lump (mode parks at the
  Q2 toll lump 0.46 vs single-mode reference 0.891 — CP-3 frozen at
  blend-MEAN level, bimodality a named probe discriminator) and n₁
  forecasts 0.18–0.23, short of the 0.31 reference (CP-4 frozen at
  the forecast band). Registration-band proposal in the findings —
  the freeze is the next user gate.

**Standing exclusion (user-confirmed 2026-07-29): every (C) cell runs
`drag_state_coupling = "off"`** — the §3.5i s(n) axis is STOPPED
(all SC predictions refuted; coupling measured GATE-CLIPPED: no ion
reaches n ≤ 8 inside the droplet). (C) does not resurrect it: the
in-droplet gate-clipping measurement is unchanged by exit stripping,
which acts at the outbound crossing where the spatial drag gate ends.
"off" is the config default and bit-identical (S1–S4); the cfg-diff
oracle vs the committed h405 reference catches any stray activation.
The module stays stubbed behind its enum per the architecture rule
(interchangeable enums; Tier-3 precedent) — excluded from runs, not
deleted.

MD cells (N = 1000, CRN on all pre-existing streams; the channel
stream is new and documented as such — pairing with the standing
batteries holds for everything except the channel draws):

- **C-full**: h405 pins + mixture (nominal weights) + strip (prior
  box center).
- **A-only**: single channel at the retiring 2.7006 + strip on —
  attributes (A) in isolation against the §3.5j bracket.
- **B-only**: mixture on, strip off — attributes (B); expected to
  reproduce the §3.5k complementarity (fast branch parks suppressed).
- **Coupling arm**: C-full with f_int,Q3 at the bracket ends
  (0.0985 vs 0.15) — the midHot/supp trade measured in-mixture.

Basin note: the (C) probe is a *mechanism* probe at the h405 pins; a
drag-surface re-tune (v_c/τ/E₀) is expected downstream if the
mechanism lands (bud226k measured the uniform-re-anchor breakage; the
mixture changes the game but is not presumed to preserve the pins).

## 10. Open questions — ADJUDICATED (user, 2026-07-29: "I follow your recommendations")

- **OQ-A (route) — CLOSED: emulation v1.** Per-pair scale emulation
  (§3.1 exactness argument); `highly_charged_iodine` /
  `single_charge_ionization_allowed` stay refused (no scope change);
  un-refusing recorded as the upgrade path if the partner observables
  are ever promoted to targets.
- **OQ-B (channel energies) — CLOSED: adopted.** f = 0.80 → Q2 2.16 /
  Q3 4.32 as the in-mixture channel means; the scalar 2.7006 budget
  retires with the mixture; the basin re-tune is downstream scope
  (the probe runs at the h405 pins and reports what breaks).
- **OQ-C (weights) — CLOSED as procedure; VALUES freeze at
  registration.** Anchor set approved (I²⁺ covariance → w_Q3 ≈ 0.20
  at the reference condition; gas slow band + Q2 partner structure
  for the w_single/w_Q2 split; provisional split (w_single, w_Q2,
  w_Q3) ≈ (0.30, 0.50, 0.20) with wide bounds). Final numbers frozen
  at probe registration after P2/P3; the m/q-127-row ask remains the
  clean closer.
- **OQ-D (strip form) — CLOSED: family approved.** Logistic G(j) +
  power-law P₀(v) with the prior box (a ≈ 2, j₀ ∈ [1, 3], w_j ≈ 1);
  P1 may move values inside the box, not the family.
- **OQ-E (checkpoint v8) — GRANTED, scoped.** Schema v8 approved
  for exactly three per-ion fields: channel label, sampled E_m
  [eV], strip count. Load-time incompatibility fails loudly per the
  checkpoint rules; v7 files stay readable via the existing
  versioned-load path. (Recorded as an explicit forbidden-list
  exception, granted under the blanket adjudication — flagged so the
  user can revoke; scope creep beyond these three fields re-opens
  the OQ.)
- **OQ-F (Q3 coupling) — CLOSED: probe arm.** f_int,Q3 measured at
  the bracket ends (0.0985 / 0.15) as the coupling arm; no single
  value pre-committed.
- **OQ-G (Q4) — CLOSED: stays refused for v1.** The bare upper-lump
  rider stands (read (i)); bare remains a soft anchor.
- **OQ-H (E_int on full strip) — CLOSED: discard-with-ledger-label.**
  The complex ceases to exist; the residual E_int is booked as a
  labeled marginal (A8-collective precedent), never silently dropped.
- **OQ-I (partner mask) — CLOSED: True, non-optional in v1.**
- **OQ-J (multiple crossings) — CLOSED: every outbound crossing.**
  Bulk-refill re-dresses on re-entry, so each outbound pass strips
  again; the probe reports the trapped-class sensitivity.
- **OQ-K (condition mapping) — CLOSED (user confirmation,
  2026-07-29).** "I confirm the peak intensities are 300 mW and
  600 mW equivalently" — the mapping 300 mW ↔ 1.47×10¹⁴ /
  600 mW ↔ 2.94×10¹⁴ W/cm² stands as a confirmed fact. The 600 mW
  `ihe_ked` reference pairs with the ~20 % Q3 covariance share;
  **w_Q3 ≈ 0.20 is no longer provisional** and the OQ-C weight
  freeze is un-gated (values still freeze at registration, after
  P3). Rider recorded: the m/q-127-row covariance decomposition
  (the clean w_single/w_Q2 closer) is **postponed by the user until
  the (C) method proves successful**.

## 11. Outcome pre-commitments — ADJUDICATED (user, 2026-07-29: "I agree with your recommendations")

Recorded **before** the probe build so outcome interpretation cannot be
post-hoc rationalized:

- **PC-1 (Tension 1 interpretation).** If C-full lands the CP-3 mean
  but the n = 1 KED stays clearly bimodal against the single-mode
  0.891 reference, the pre-named candidate deficiencies are
  (a) the designed-out intensity–position correlation (channel draw ⊥
  geometry; focal-volume averaging would smear channel weights
  continuously and fill the valley) and (b) missing channel structure
  (the single-channel idealization; the Q4 rider). (a) is the
  physically favored valley-filler. No other explanation may be
  introduced without a new registered read.
- **PC-2 (sequencing on a pass).** Probe passes CP-1..8 →
  (i) mechanism provisionally adopted at the h405 pins,
  (ii) drag-surface re-tune (v_c/τ/E₀) with the mixture ON (G3-ring
  machinery reused), (iii) only then the G4 adjudications close
  (successor point, retained-policy finalization from the
  `exclude_all_coupled` interim, ledger re-issue) against the
  re-tuned point. G4 does NOT close before the re-tune.
- **PC-3 (failure semantics).** CP-5 firing at a = 2 indicts the
  P₀(v) *family* (a = 1 already disfavored in-box); CP-7 firing
  indicts the G(j) family (P1 may only move values inside the box,
  OQ-D) — either family-level kill stops the (A) strip design and
  returns to the §3.5j alternatives. CP-6 firing at **both** coupling-
  arm ends closes the current (C) f_int wiring but not the mixture
  itself; at one end only, the other end is the surviving arm.
- **PC-4 (m/q-127 trigger).** "Method proves successful" = a CP-band
  pass on the probe (not post-re-tune adoption): the postponed
  bare-I⁺ m/q-127-row covariance ask fires at probe-pass, so the
  re-tune runs with properly anchored w_single/w_Q2 instead of
  re-tuning twice.
- **PC-5 (CRN scope note, build-level).** CRN pairing with the
  standing batteries holds on all pre-existing streams, but C-full /
  B-only break *per-molecule* pairing at the source (every E_m
  differs). Seed-level comparisons against the h405 battery are clean
  only for A-only (single channel at the retiring 2.7006; strip
  Bernoullis aside). The probe scorer must state this explicitly.
- **Tension-2 note.** An MD n₁ landing *above* the CP-4 forecast band
  is diagnostic (kernel saw first crossings only; OQ-J re-crossings
  can convert trap-class ions), not a miss; the residual toward the
  0.31 reference is expected to wait for the PC-2 re-tune.

## 12. Cross-links

- Plan: `TIER2_SENSITIVITY_ATLAS_PLAN.md` §3.5j–§3.5l.
- Findings: "§3.5j", "§3.5k", "§3.5l reads (i)+(iii)".
- Influence: `TIER2_PARAMETER_INFLUENCE.md` §19 (budget), §17
  (suppression scaffolding ledger).
- RQs: RQ7 (budget provenance, closed 2026-07-29), RQ8 (channel
  assignment + covariance anchors), RQ3 (bare/n = 0 fate, fed by
  OQ-L-allowed outcomes), RQ11 (deep-bin KE, untouched here).
- Paper: Hatherly, Stankiewicz, Codling, Frasinski & Cross,
  J. Phys. B: At. Mol. Opt. Phys. 27 (1994) 2993–3003
  (doi:10.1088/0953-4075/27/14/032; consulted 2026-07-29, not kept
  in-repo — the design-relevant numbers are quoted in §1/§3.5).
