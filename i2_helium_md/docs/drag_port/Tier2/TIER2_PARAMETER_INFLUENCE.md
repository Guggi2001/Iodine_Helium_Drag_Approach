P# Tier 2 — Parameter Influence Reference (compact)

> **What this is.** The D0 deliverable of `TIER2_SENSITIVITY_ATLAS_PLAN.md`:
> a **pure distillation** of `TIER2_STAIRCASE_PROBE_FINDINGS.md` (§4a–§4ee,
> I1–I100) into per-knob influence statements. **No new analysis** — every
> number here is a pointer-backed quote of a measured result. The findings
> doc stays the archive/provenance; this doc is the navigation layer.
>
> **How to read an entry.** *Role* = where the knob enters the mechanism.
> *Class* = CALIBRATION_MAP calibration class (Sourced / Bounded / Free /
> Derived). *Influence* = measured effect on the observables, with
> magnitudes and insight pointers (I-numbers → the findings-doc register;
> § → findings-doc sections). *Couplings* = structural entanglements with
> other knobs. *Status* = locked / standing / open / **GAP** (no measured
> influence yet — atlas target).
>
> **Living doc.** Each sensitivity-atlas stage merges its results in and
> closes its GAP marker. Last built: 2026-07-23 (distillation of I1–I100).
>
> **Standing point (context for every "standing" value):** finc1v725 —
> `capped_cubic` v_c 7.25 / p_tail −1, τ 3.2 ps, E₀ 0.27 eV, `rq4graded`
> ladder, Landau v_L 0.58 Å/ps, locked Tier-0 b, E_bind 0.1168 eV,
> margin 3 Å, N = 5000 pooled battery reference (§4cc).

---

## 0. Summary table

| knob | class | headline influence | status |
|---|---|---|---|
| drag form + b (in-band) | Derived (Method-B/TDDFT) | owns KE scale + fragment composition; sole surviving RQ11 owner-candidate | locked in-band; form question OPEN (RQ11) |
| v_c, p_tail (capped tail) | Free | v_c owns the whole mid-bin KE curve; joint-landing basin v_c ∈ [7.25, 7.5]; p_tail is not a second lever | standing 7.25 / −1 |
| τ (cooling clock) | Bounded | descent clock; joint closure only as a (v_c, τ) pair; race coordinate with f_int under the gate | standing 3.2 ps |
| E₀ = E_int(0) | Bounded [0.2, 0.5] eV | fate-cliff position (bare↔shell split); re-landed at solvation scale 0.22–0.27 under full geometry | standing 0.27 eV |
| κ (ladder steepness) | Free | near-dead staircase lever (inverted + normalization-capped) | standing (rq4graded supersedes) |
| picture (electronic) | Free (selection) | magnitude-degenerate ≤ 4 %; discrimination belongs to the size distribution | standing `statistical_mixture` family |
| f_int | Derived coordinate | superseded by absolute E₀; gated: only the race margin Δ× is physical | retired as free knob |
| f_ret | Bounded (prior small) | pinned 0.1 throughout; never swept | **GAP** |
| λ₀ pickup (+ p cap) | Sourced + Bounded (p conditional-Free) | re-filling ≈ 7 ions/100; p = 1 de-suppresses to the twin and lifts n₁ KE toward experiment | standing p = 1 |
| s_eff (RRK dof) | Bounded | cascade *rate*, arm-conditional (8 ungated / 30 gated); detector-compressed except the n = 1 bin | standing 8-family |
| ladder D₀(n), Σ | rungs Sourced; shape RQ4-graded | owns n₁/small-n shape; rq4graded > floor1 ≫ 2σ; deep tail Σ-coupled ×20 too stiff to move deep KE | standing rq4graded |
| Landau v_L | Sourced 0.58 Å/ps | bit-flat on scored surface across 0.30–0.58; acts only on the retained class | standing 0.58, quiet |
| per-shed ε | Bounded small (NB-RQ23-1) | refuted as RQ11 owner; ε ≈ 1–2 meV/shed sub-dominant compatible | ε = 0 standing |
| shed momentum convention | convention (two-valued) | cold-shed injects ×1.611 KE over a full strip; histogram convention-blind; bare-bin KE reads the fragmentation convention | co-moving basis for twin parity; RQ3-coupled |
| birth margin | pinned convention (3 Å) | **the sensitive robustness lever**: 6 Å moves n₁_solv −0.171 (3.4σ) | standing 3 Å; in-tier open item (I88) |
| E_bind (ion–droplet well) | Derived (joint Method-B) | gates the trapped class (5–11 %); never swept | **GAP** |
| droplet geometry (R × r) | sampled, never controlled | fate is droplet-size + birth-position ordered; prior decoupled from KE axis | **GAP** (Axis A) |
| sampling laws (size + position) | theory-laden legacy ports | unaudited (E_solv 14 vs 30 meV discrepancy known) | **GAP** (D2b) |

---

## 1. Drag law — in-band form + coefficient b

| | |
|---|---|
| role | continuous friction γ(v) = ρ̂·b·v² (force ρ̂·b·v³) inside the bubble; owns all in-window dissipation |
| class | Derived — Method-B extraction from 9/18 Å TDDFT traces (bands: 18 Å 2.54–3.02, 9 Å 2.83–4.95 Å/ps); locked Tier 0 |

**Influence (measured):**

- Owns the KE *scale* and the fragment *composition*: n₁/bare occupants
  are near-edge births losing 0.8–1.2 eV to drag vs 2.6–2.8 eV for deep
  bins (I60). Geometry owns the KE *shape*; the envelope [current-law,
  ballistic] brackets experiment everywhere (I43).
- Over-dissipates at production speeds: detected model speeds undershoot
  the VMI I⁺He peak ~2–2.5× (I40); measured K₂.₇₀ = 0.746 with 41.7 % of
  exposure beyond the calibrated band (I37).
- The RQ11 deep-bin cold tail (sim/ref 0.81 → 0.38 over n = 10–17, a
  slope not a scale, I92) accrues in the **in-band pure-cubic phase at
  v ≈ 2–7 Å/ps, t ≈ 2–14 ps, in-bubble** (I90); with candidates
  (i)/(ii)/(iii) refuted, the in-band law at low v is the surviving
  owner-candidate (I96) — ~1.5 % of the deep survivors' dissipation.
- Low-v form space mapped in the twin (§4ee): hard/erf floors are
  perfectly speed-selective and act by un-trapping the marginal retained
  class (trap 0.042 → 0); the sub form at v_f ≈ 1.5–2.0 is the twin's
  landing window (deep → 1.06…0.83, midHot ≤ 1.14); the parameter-free
  shifted-cubic at locked b over-corrects globally (all bins ×1.2–1.4,
  midHot 1.40 — though W₁ improves 0.678 → 0.502) (I98/I99).
- Literature: above ~1 Å/ps in a droplet the drag *shape* is
  unconstrained — TDDFT is the only authority; the (v−v_L)³ law is
  near-threshold/negative-ion/mean-of-sawtooth and does not transfer
  (I97, I100, NB-RQ11-1..9).
- Atlas D4 Step 1 (2026-07-23, artifact reuse, zero MD): the three-mode
  form table stands as compiled — lq rejected (shared Δobj +0.0339,
  held-out 9 Å FAIL) with its low-v limb ×2.54 the production drag at
  v = 2 (wrong direction for RQ11); power law cubic-equivalent
  (n̂ = 2.927). Tail re-inspection: **no drag-dominated trace samples
  below v ≈ 2.5 Å/ps** — the sub-2.5 shape is TDDFT-blind pending the
  collaborator ask. The subtractive gated form is trace-constrainable
  in-band (×0.65–0.38 at 2.54 for v_f 1.5–2.0 vs the whole 18 Å window
  at 2.54–3.02); its re-fit is a pending decision. Full tables:
  `TIER2_SENSITIVITY_ATLAS_FINDINGS.md` (NB-RQ11-10).

**Couplings:** jointly extracted with E_bind (§9 below); KE↔histogram
anti-correlate through the exposure integral K (I47) — any drag change
re-opens the (v_c, τ, E₀) arbitration (I99).

**Status:** locked in-band; the *form* question is RQ11's open axis
(lever hierarchy: TDDFT re-inspection / Method-B form re-fit ≫ candidate
(iv) ≫ Free phenomenological suppression). After atlas D4 Step 1: priors
unchanged; Padé parked (see §2), subtractive fit-decision pending,
shifted cubic's re-add condition not met. After the §6.6 quadratic
counterfactual (2026-07-24, corrected same day by the user-caught fine
rescan): **outcome (b) — the Tier-2 observable vector is measured to be
form-blind.** A re-arbitrated lq system (v_c 8.8–9.0, τ 3.3–3.5,
E₀ 0.27) lands the histogram + KE observables comparably to production
cubic in the twin (six four-way passes; best cell beats base W₁ and
midHot centering; n₁KE ~10 % colder; needle-width basin vs cubic's §4w
basin). Consequence: the form choice rests **solely on the Tier-0
trace instruments** (lq stays rejected there: held-out 0.699 FAIL,
n̂ = 2.927), and the Tier-2 landing must never be cited as evidence
for cubic — the anti-circularity caveat is upheld as a measured fact.
Pre-registered follow-ups open: MD spot-check of an lq-passing cell;
variable-mass Method-B re-run.

## 2. Capped tail — v_c and p_tail

| | |
|---|---|
| role | production arbitration's above-band extension: cap at v_c with tail exponent p_tail |
| class | Free (the program's genuinely fitted knobs, arbitrated by experiment) |

**Influence (measured):**

- **v_c owns the entire mid-bin KE curve**: v6.5 → v7.5 at p−1 moves
  midHot 1.74 → 0.95 with n2–n8 within ~10 % (I78). The KE χ² argmin was
  right-censored at the bracket top until the taper was understood (I71).
- **Joint landing is a basin, not a knife-edge**: v_c ∈ [7.25, 7.5] at
  (τ3.2, E₀0.27, rq4graded); W₁_solv bowls at 7.25 (0.424), n₁_solv
  crosses 0.31 at ≈ 7.25, midHot crosses 1.0 at ≈ 7.2; closes above
  ~7.75 (whole curve ×0.68–0.88 at v8.0) (I80/I81).
- **p_tail is the wrong second lever**: a more-negative exponent re-heats
  the whole KE curve *and* re-strips — KE and histogram move together,
  not orthogonally (I78). p = 1 tails excluded outright (I51).
- Three observables prefer three v_c (KE χ² ≥ 8.5, n₁ anchor ≈ 6.2–6.3,
  W₁ ≈ 7.25) — a drag-form-*shape* statement, not a calibration error
  (I72).
- Hard cutoff (tail removed) excluded: 99.7 % suppression → bare at
  v_c ≤ 10 (I48).

**Couplings:** with τ (joint closure exists only as a pair, I49); with
the ladder — orthogonal (ladder is KE-neutral: midHot flat across the
ladder family at every v_c, I79).

**Status:** standing v_c 7.25 / p_tail −1; every sensitivity-ring
perturbation at the blessed point degrades W₁ (+0.18..+0.33) — the
standing cell sits at the basin optimum (I88). The candidate smooth
cap-replacement (Padé saturating cubic, one knob v_s) is **excluded by
zero-cost arithmetic** (atlas D4 Step 1): in-band cubicity forces
v_s ≳ ~15 Å/ps, which forfeits the p_tail = −1 saturation at the
production peak — the principled cap-remover remains extending TDDFT
authority to production kinematics (the collaborator ask), not a
different closed form.

## 3. τ — Newton cooling clock

| | |
|---|---|
| role | E_int Newton-cooling time constant; sets the evaporative-descent clock and the in-bubble leak |
| class | Bounded (GAH25 band; re-classed from the 6.55 pin, §I.8) |

**Influence (measured):**

- Alone, cannot fix the KE↔histogram anti-correlation; **(v_c, τ)
  jointly close** what neither does alone (18-cell basin at τ ≈ 4.1;
  refined to τ 3.2 at the c1 winner, τ 3.8 at the c4 basin) (I49, I82).
- Under the gate, τ and f_int are **not independent** — only the race
  margin Δ× = t_eject − t× is physical; the gated arm is binary in the
  t×↔ejection race (cliff, re-opening at τ 16.5 × f_int 0.25) (I7, I16,
  I18).
- s_eff↔τ separate through first-shed timing (I4).

**Couplings:** v_c (joint basin); f_int/E₀ (race coordinate); the
cooling exposure K is the single integral through which KE and stripping
trade off (I47).

**Status:** standing 3.2 ps.

## 4. E₀ = E_int(0) — internal-energy budget

| | |
|---|---|
| role | absolute internal energy at onset [eV]; feeds the RRK gate via the Σ(21) crossing |
| class | Bounded — RQ1 band [0.2, 0.5] eV, floor ≈ 0.22 (row 16; replaced f_int as the physical variable, I25) |

**Influence (measured):**

- Sets the **fate cliff**: closed-form bare ⇔ K < K* = ln(E₀/Σ(21));
  the (E₀, K) two-parameter fate map reproduces measured columns to
  0.84 He (I34). At a congruent single-point ensemble the suppressed
  fraction is a delta step — coexistence requires ensemble heterogeneity
  (I26).
- The experimental histogram inverts into a smooth unimodal p(E₀)
  coinciding with the RQ1 band untuned (W₁ = 0.086 at the probe basis)
  (I29–I31).
- **Scale is geometry-dependent**: pinned-droplet reads wanted
  E₀ ≈ 0.38–0.41 (I39); full position geometry re-lands the optimum at
  the solvation scale 0.22–0.27 eV (I45). Standing 0.27.
- Budget-corner regime change: c1's suppressed channel closes entirely
  at E₀ = 0.23 (supp ≈ 0), inverting orderings (I82).
- The scenario *budget* (0.80 vs 2.70 eV) is bookkeeping-only in the
  delivered model — the gated map is budget-invariant in absolute E₀
  (I28).

**Couplings:** τ (race), droplet-K axis (the cliff position moves with
exposure: E*(K) = Σ(21)·e^K, I39), sampling width δ (bare-fraction
trade-off, I35/I36).

**Status:** standing 0.27 eV.

## 5. κ — ladder steepness (Form U)

| | |
|---|---|
| role | sigmoid steepness of the D₀(n > 1) ladder |
| class | Free (one of the two original Free knobs) |

**Influence:** near-dead: inverted *and* normalization-capped for the
21→14 staircase (Form-U floors D₀(21) at ≈ 0.53·D₀(1)) (I2); Σ(21) is
nearly κ-independent (~11 % over the pure-κ range, CALIBRATION_MAP row
21). Superseded in practice by the rq4graded ladder family (§8).

**Status:** standing within rq4graded; not an active lever.

## 6. picture — electronic-state selection

| | |
|---|---|
| role | selects the D₀ curve family (statistical_mixture / x2_only / cooling_relaxed) |
| class | Free (selection) |

**Influence:** magnitude-degenerate — ≤ 4 % spread from x-invariance;
its only staircase signal is gate-open timing via Σ(21), degenerate with
f_int (I3); landing band picture-robust (Wave 3). Discrimination belongs
to the terminal size distribution, i.e. Tier-2 arbitration, where the
mixture family stands.

**Status:** standing `statistical_mixture` family (c1 bundle).

## 7. f_int, f_ret, λ₀, p — budget fraction, retention, pickup

- **f_int** (Derived coordinate = E₀/E_avail): parametrization broken in
  both directions (no mechanical partition; literal coupling ~1 % would
  be inert) → replaced by absolute E₀ (I25). Under the gate it was a
  joint timing + budget knob (I15); only the race margin Δ× is physical
  — independent τ × f_int grids sample the race incidentally (I18).
  Retired as a free knob.
- **f_ret** (Bounded, prior small): pinned 0.1 through the whole probe +
  campaign chain; **never swept — GAP** (atlas candidate if ever
  motivated; no measured influence on record).
- **λ₀ pickup** (Sourced + Bounded, ~0.7–1.1/ps): live Langmuir
  re-filling is real and measurable — MD suppression sits ≈ 7 ions/100
  below the frozen twin at leg B, ≈ 1–3/100 at p = 1 (I63, I65).
- **p occupancy-cap exponent** (conditional-Free, standing p = 1):
  p = 0 → 1 de-suppresses ≈ 3× to the twin's prediction (best twin↔MD
  W₁ 0.31–0.47 of the chain) and lifts solvated n₁ KE to ≈ 1.1 eV —
  closest to the experimental n₁ on the solvated branch; touches only
  the onset, never chord dynamics (trapped p-invariant) (I65/I66).

## 8. Ladder — D₀(n) rungs, tail, Σ

| | |
|---|---|
| role | per-shell binding ladder: RRK gate, suppression criterion (E_ej > Σ(n₀)), descent bookkeeping, e_bind E_pot fold |
| class | D₀(1) Sourced (±3 cm⁻¹, IHe05); shape RQ4-graded (2.2 : 1.5 : 1.3 diagnostic); floor Sourced (bulk µ_He) |

**Influence (measured):**

- **Owns the small-n histogram shape**: n = 1 is a one-rung-wide E_ej
  window — the taper controls n₁; equal-deep-rung families revert n₁/n₂
  to the flat chord-K ratio (I42). The solvated histogram independently
  demands the RQ4 target ratios through a purely geometric forward model
  (I42).
- **rq4graded > floor1, decisively**: N = 500 verdict on every histogram
  read (W₁ 0.524/0.657 vs 0.843, ≫ 2σ; c4's KE curve also fails
  structurally) (I86). Cheap-bottom ladders strip *past* n = 1 into bare
  and *lower* n₁_solv — rq4graded's 20 meV D₀(1) barrier is what piles
  weight at n = 1 (I79).
- The ladder is **KE-neutral** (midHot flat across the family at every
  v_c) — orthogonal to the drag axis (I79).
- **The deep tail cannot move deep-bin KE**: deep-KE and histogram are
  rigidly Σ-coupled, ~×20 too stiff — ±5 % deep rungs ↔ ∓3.5 pts
  suppressed weight; ×1.05 already bleeds supp 0.202 → 0.167, ×1.2
  wrecks midHot (1.46); Σ-preserving transfers blow W₁ to 1.1–2.0.
  Corollary: the seed-robust histogram landing is itself evidence the
  deep rungs are approximately right (I93, §4dd).

**Couplings:** Σ(21) sets the gate threshold (κ-insensitive); the
suppression criterion and KE re-mapping read the *same* tail (the I93
stiffness); picture selects the curve family (§6).

**Status:** standing rq4graded; RQ4 external arbitration remains the
blocking authority for the taper's *physics* (I42).

## 9. E_bind — ion–droplet surface well (0.1168 eV)

| | |
|---|---|
| role | depth of the mean-field droplet exit barrier for the ion (`binding_energy_I_ion_eV`); decides eject vs trapped |
| class | Derived — jointly extracted with the drag coefficients (Tier-0 Method-B); §6.5.1 guard enforces the exact pairing |

**Influence (measured, indirect only):**

- Defines the **trapped droplet-retained class**: 5–11 % of fragments
  (inward partners of off-center births, chord K up to ≈ 18) dissipated
  below the 0.117 eV barrier, never ejecting; includes centrifugal
  resonances; fragile at meV scale (1–3 ions/config flip at a
  convention change) (I44, I56, I61).
- The §4ee low-v twin sweep moved this boundary *from the drag side*:
  un-trapping the marginal class (trap 0.042 → 0) was the floors'
  dominant effect — an in-model echo of RQ11 candidate (iv).
- **Never swept directly — GAP** (Axis C §6.5 of the atlas plan; runs
  require `allow_unvalidated_binding_pairing`).

**Couplings:** drag law (jointly calibrated pair — "correct drag traps
the ions", TIER0_FINDINGS); *not* connected to the ladder D₀(n)/Σ
(different bookkeeping surfaces: droplet exit well vs per-atom shell
energetics).

## 10. Landau v_L — dissipation threshold

| | |
|---|---|
| role | γ = 0 gate below v_L (0.58 Å/ps); also a latent shared reader (sets the neutral-stage E_min threshold, I84) |
| class | Sourced — bulk roton Landau value, droplet-persistent down to ~1000 atoms (NB-RQ11-1..3) |

**Influence (measured):** **bit-flat on every scored observable** across
v_L ∈ {0.30, 0.40, 0.58} (Δn₁_solv = ΔW₁ = 0.0000) while demonstrably
acting on the retained class (retained KE 0.277 → 0.005–0.023 eV,
monotone) — the scored surface is structurally insulated, not the arm
inert (I84, re-confirmed at the blessed point I88). The Landau floor
(≈ 0.003 eV) sits an order of magnitude below the scored deep bins —
part of why RQ11 candidate (i) fails (I89).

**Status:** standing 0.58, quiet on the scored surface; literature
deep-read confirms the value but shows the *gate shape* above v_L is
unconstrained (I100).

## 11. Per-shed ε — evaporative recoil

| | |
|---|---|
| role | translational energy carried per shed He (currently ε = 0; OQ-F origin, I23) |
| class | Bounded small (NB-RQ23-1: ~0.5–0.6 meV prescription) |

**Influence (measured):** required ε_req(n) to close RQ11 is
sign-changing (−3…−4 meV/shed at n = 2–3, +3 → +9 over n = 10–17); flat
ε ≈ 6–8 meV closes the deep slope (χ²_med 242 → 66) only by breaking
midHot (1.25–1.33 vs [0.84, 1.14]); physical drag-retention is
*inverted* (deep-bin sheds happen in-bubble, 43 % retention vs 95 % at
n₁) — retention-weighted forward model floors at χ²_med ≈ 198.
**Refuted as RQ11 owner; ε ≈ 1–2 meV sub-dominant stays compatible**
(I94, §4dd).

## 12. Shed momentum convention (cold-shed vs co-moving)

| | |
|---|---|
| role | what the evaporated He carries away: at-rest-in-lab (cold-shed, momentum-conserving) vs co-moving (ε ≈ 0) |
| class | convention — RQ2/RQ3-coupled (OQ-J) |

**Influence (measured):** cold-shed injects +0.91 eV over a full
post-exit strip (×1.611 = m₂₁/m₁ in KE) — it *was* the twin−MD n₁ KE
divergence, closed entirely by the co-moving counterfactual (I59, I61).
The histogram is convention-blind (W₁(A″, A′) = 0.14–0.20); the
suppressed-class *fragmentation* read is convention-decided: pooled
conventions give intact 1.96 / co-moving 1.36 / cold-shed 2.82 eV vs
ref bare mean 3.706 — direction momentum-conserving, all 6–10× too
narrow in σ, and the bare mean is kinematically un-sourceable from the
production channel (cap ≈ 3.37 eV pre-drag) → n = 0 is a two-channel
mixture; the suppressed class belongs to the low-KE shoulder (I60, I64,
I95).

## 13. Birth margin — 3 Å pinned convention

| | |
|---|---|
| role | minimum birth depth below the droplet surface for the sampled I₂ center |
| class | pinned convention (I88) |

**Influence (measured):** **the sensitive robustness lever** at the
blessed point: margin 3 → 4.67 / 6.0 Å moves n₁_solv −0.094 / −0.171
(≈ 1.9σ / 3.4σ), collapses the ratio, halves suppression — the
histogram landing depends on this pin (I88). In-tier open item.

**Couplings:** position sampling (§15) — the margin truncates the
Boltzmann radial law near the surface; the position axis is a two-sided
race lever (I55).

## 14. Droplet geometry — size R × birth position r **(GAP — Axis A)**

Never controlled independently; everything known is inference from
sampled ensembles:

- The droplet-radius axis has a **quantified demand**: ≥ 9.6 % of the
  ensemble below the 0.22 eV solvation floor requires a larger-K
  (larger-droplet) tail (I30); width attribution: the ~21-bin span must
  be carried by the E₀ and/or droplet-K axes — per-point reads are
  ≤ 1 bin wide (I33).
- The **position axis alone opens the two-sided race** (first MD weight
  at n = 0/1 with no dressing; suppressed/bare ordering transfers
  exactly from the twin) (I55); center-pinning was half of the T3 park
  (I52).
- The **droplet prior and the KE axis are decoupled**: flipping the
  kornilov prior moves the histogram (n̄ ≈ −1 He) while solvated n₁ KE
  moves ≤ 0.12 eV (I68); prior axis asymmetric at the blessed point
  (δ 0.40 / pickup ≈ 1.5σ; δ 0.80 noise-level) (I88).
- **Fate is birth-dressing + droplet-size ordered** (deep survivors:
  larger droplets 30.1 Å, fuller dressing 18.4; retained: 34.0 Å)
  (I90); frozen-geometry twin bias (channel (d)) scales with in-droplet
  residence, 0 → −0.8 He across the v_c bracket (I73), (τ, E₀)-blind
  (I74).

**Atlas target:** the §3 factorial grid (R × r/R) + D2b re-weighting.

## 15. Sampling laws — size distribution + radial positions **(GAP — D2b)**

Theory-laden legacy ports, unaudited:

- **Size law** (`sampling/droplet_sizes.py`): nozzle correlation
  ⟨N⟩ = k₁·p^0.97·T^−3.88·d² (Lackner); log-normal δ = 0.625 (Kornilov);
  pickup + evaporation chain (mean free path 4 Å, ε = 0.04/collision,
  **E_solv = 14 meV vs thesis-text 30 meV — documented discrepancy**);
  0.8·ρ_bulk; analytic D4 prior family (`droplet_size_prior`).
- **Position law** (`sampling/radial_positions.py`): thermal-equilibrium
  ansatz p(r) ∝ r²·exp(−U_drop(r−R)/k_B T) at T_particles_K
  (`binding_energy_molecule_meV`, steepness).

**Atlas target:** D2b provenance audit + distribution A/B + grid
re-weighting; salvageability = landed observables stay inside the
seed-SD yardstick under all defensible variants.

## 16. Quiet / structural surfaces (for completeness)

- **s_eff (RRK dof, Bounded):** a *rate*, not an amount (I13);
  arm-conditional landing (ungated ≈ 8, gated ≈ 30 — the ungated prior
  does not carry, I10); the detector compresses it (arrival state, not
  arrival n) *except* the n = 1 bin, which selects the kinetics band
  (I20, I32). Standing: the 8-family via the production bundle.
- **E2 relaxation stage:** KE-frozen at handover on every solvated bin
  (E_dissip gain 0.0000; Landau-gated E2 drag touches only the retained
  class); detection-stage sheds KE-negligible — the post-window chain is
  causally disconnected from the KE curve (I89); E2 checkpoint stride is
  scored-read-neutral (I67).
- **Detection time:** Sourced 8.53 µs; ungated reads are
  detector-converged (n_detect ≡ n_relaxed); the gated long-time
  structure is the Klots evaporative-ensemble regime (I19, I23).
- **Scenario budget (0.80/2.70 eV):** bookkeeping-only in the delivered
  model; the gated map is budget-invariant in absolute E₀ (I28).
- **Guard flags:** biphasic runs under `allow_inconsistent_mass_pairing`
  (§6.6 mid-window defense); E_bind sweeps will use
  `allow_unvalidated_binding_pairing` (§6.5.1). Deliberate, logged.

---

## Cross-references

- Archive/provenance: `TIER2_STAIRCASE_PROBE_FINDINGS.md` (§4a–§4ee,
  I1–I100, boundaries §6, open questions §7, run inventory §8).
- Program that fills the GAPs: `TIER2_SENSITIVITY_ATLAS_PLAN.md`.
- Calibration classes + identifiability: `CALIBRATION_MAP.md`.
- Open physics questions: `RESEARCH_QUESTIONS.md` (RQ1–RQ11 + NB
  registers).
- Decision history: `drag_migration_log_tier2.md`.
