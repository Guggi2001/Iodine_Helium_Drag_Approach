1C
# Tier 2 — Staircase / Capability Probe Findings

> **Status:** consolidated findings record, written 2026-07-07 after the fourth
> probe wave (cooling-spatial-gate total-strip A/B) executed; last updated
> **2026-07-21 after the n=1-deweight KE re-score + E-2/ladder twin scans +
> the (v7.5) MD joint-landing + the v_c sensitivity ring (§4u–§4w, I77–I81 —
> the n₁ 1.302 eV anchor is a mixture-mean artifact but de-weighting it leaves
> the high-drag pull with the mid-bins n2–n8; the E-2 (v_c, p_tail) corner has
> no joint cell because p_tail re-couples the axes, while raising v_c alone
> lands the whole mid-bin KE; the ladder is KE-orthogonal but saturated at
> rq4graded, so the v7.5 histogram recovery is a twin→MD strip bias — and the
> (v7.5, τ3.2, E₀0.27) MD spot-check (§4v, I80) CONFIRMED the joint landing
> exists (n₁_solv 0.291 / n₁/n₂ 2.30 / W₁ 0.496 / midHot 0.90 / χ²_med 14,
> the first MD cell to land both axes; §4t's "no joint cell in capped_cubic"
> was a sparse-v_c-grid artifact), and the v_c ring (§4w, I81) mapped it as a
> *basin* v_c ∈ [7.25, 7.5] with the optimum refining to ≈ 7.25 on the
> histogram — the c1 N=500 finalist center), then the §I.11.5 endgame legs
> S4 + S5 (§4x–§4y, I82–I84): the c4/floor1 joint region EXISTS in MD at
> (v7.0–7.25, τ3.8, E₀0.25) — n₁_solv 0.312/0.303 dead-on, W₁ ≈ 0.63 vs c1's
> 0.42–0.50, so the N=500 ladder contest is real with the c4 finalist at
> (v7.0, τ3.8, E₀0.25) — and the Landau v_L gate DISSOLVES (scored
> observables bit-flat across v_L ∈ {0.30, 0.40, 0.58} while the arm acts on
> the retained class; N=500 runs Landau-on at 0.58, external re-pinning
> non-blocking)**; **then S6 (§4z, I85–I87): the N=500
> finalists land NO cell under the frozen joint bars — N=500 unmasks a
> systematic deep-bin KE miss (χ²_med 69–202) while the robust
> histogram metrics carry and the ladder contest resolves for
> c1/rq4graded; riders winner-gated, not fired**;
> re-pilot Stage 2 is §4t (I74–I76 — the (τ, E₀) grid at two carried v_c;
> the twin is quantitative at v65; each axis lands at different cells and no
> joint cell exists in capped_cubic); re-pilot Stage 1 is §4s (I71–I73 — the
> v_c bracket in real MD: χ² argmin right-censored twin+1-step; the
> speed-selective tension survives; channel (d) scales with drag exposure);
> T9 leg D is §4r
> (I68–I70 — the kornilov droplet prior; DP-P3/P4 confirmed, DP-P1/P2 split);
> T9 leg C is §4q (I65–I67 — the T6 `sigma_proportional` p-law
> de-suppresses to the twin; W₁ 0.31–0.47, the chain's best; n₁ KE moves toward
> experiment); T9 leg B is §4p (I62–I64 — the T5 dressing transfers
> quantitatively; pickup re-filling measured; the bare-bin KE gap is RQ3
> bookkeeping); T9 leg A″ is §4o (I61); the I.11.2 item-1
> n₁-composition re-read is §4n (I58–I60, OQ-J); T9 leg A′ is §4m
> (I55–I57); the Slice-T3 pilot
> first look is §4l (I52–I54); Wave 11 (the RQ7 production-kinematics
> probe — K₂.₇₀ measured) is §4i; Wave 10's width decomposition + K
> forward model is §4h; Wave 9's E₀-mixture inversion is §4g;
> Wave 8's cliff anatomy is §4f. This document
> collects **every insight and numerical result** of the pre-F5 probe program;
> the H.2b analytic feasibility pass (zero-MD, 2026-07-11) is §4j (I41–I45)
> in one place; the decision/delivery history stays in
> `drag_migration_log_tier2.md` (entries of 2026-07-05/06/07) and the design
> contracts stay in `TIER2_STAIRCASE_PROBE_PLAN.md` (+ Addendum A) and the
> cooling-gate plan referenced there.
>
> **Stance (unchanged throughout):** existence/capability probes, **reported,
> not auto-adjudicated**. TDDFT is not ground truth (the rejected "1b" stance
> holds; the 9 Å non-radial flag stands); experiment arbitrates at Tier 2 via
> the terminal I⁺Heₙ size distribution. **Nothing here discharges the F5 gate**
> — that is resolved only by the N=500 Stage-1/2 campaign at 0.80 eV.

---

## 0. What the probe program is

Before the F5 production switch (2.70 eV) and the N=500 campaign spend, four
small-N probe waves asked, in sequence:

| Wave | Runs | Question | Outcome |
|---|---|---|---|
| 1. Full-grid capability probe (2026-07-05) | 45 | Does any in-band (κ, picture, τ) land the anchored 9 Å staircase 21→19→14? How flexible is the terminal shell? | **Miss everywhere** — freeze at n ≈ 20; RRK-dof mechanism-level OQ fired |
| 2. `s_eff` mini-probe (2026-07-06) | 12 | Is the freeze the classical RRK dof convention `s = 3n−3`? | **Yes** — constant s_eff ≈ 8 lands magnitude *and* timing |
| 3. Picture cross-check (2026-07-06) | 10 | Does the 45-run picture-near-degeneracy prediction hold under reduced s? | **Confirmed** — ≤ 4 % magnitude spread; s_eff landing band picture-robust |
| 4. Cooling-gate total-strip A/B (2026-07-06/07) | 24 | Can the mechanism express total stripping (the experimental 43 % bare-I⁺ peak)? What does density-gating K2 cooling change? | **Near-total strip reachable** (n̄ ≈ 2.7) with the gate on at mid-band τ; gated arm is **all-or-nothing in τ**; bare n = 0 not reached (floor n = 2) |
| 5. Gated landing re-location (2026-07-07) | 5 | Under the production gate, where does the s_eff landing move (the ungated [8,12] prior is gate-conditional)? | **Landing shifts to s_eff ≈ 30** (gate ~4× the in-window shedding); **but the gated in-window staircase and the terminal read decouple** — s_eff sets the *rate*, not the endpoint, so the terminal is flight-time-dependent |
| 6. f_int probe (2026-07-07) | 10 | The one knob every wave pinned: what does f_int control under the gate — timing only, or the race and the budget? Can it re-align staircase timing / reach bare? | **f_int is a joint timing + effective-budget knob gated** (the in-bubble leak over [t×, t_eject]); cliff located (0.50 < f_int* < 0.65 at τ=6.55); the dead τ=16.5 arm re-opens at 0.25; timing re-aligns at ≈ 0.40–0.42; **bare unreached at any f_int** |
| 7. Detected-read re-score (2026-07-08) | 0 (re-read of all 102) | What is the terminal read at the Sourced detector time t_detect = 8.53 µs, per (gate, s_eff, τ, f_int)? Does the gated landing point arrive at the floor or mid-shell? | **Ungated ≡ relaxed (detector-converged)**; the gated f_int = 0.50 arm is **s_eff-compressed at the detector** (n̄ 2.7–4.1 over s_eff 1→30; the landing point arrives 100 % live at n̄ = 4.09); early-open leak points survive as **converged mid-shell weight**; dead arms are literal n = 21 detector weight; **bare unreached** (min n = 2) |
| 8. Suppressed-fraction cliff anatomy (2026-07-09) | 6 (+ zero-MD step 1) | OQ-B weight-level: does a high-enough E_int(0) split the ensemble ≈ 43.5 % suppressed (bare-candidate) / 56.5 % opened? How wide is the Wave-6 cliff? | **The cliff is a delta at E\* = 0.4612 eV** (ensemble kinematically congruent → no split at any E_int(0): 0 % or 100 % only); the **sub-rung sliver below E\* expresses n = 1** (all ions, first time); the **budget knob is bookkeeping-only** (K₂₇₀ ≡ K₀₈₀; absolute-E_int(0) map budget-invariant; OQ-H fired) |
| 9. E₀-mixture inversion (2026-07-09) | 8 (+ pure-fit) | Invert the experimental histogram into an implied p(E₀): does a non-negative mixture of single-point detected distributions reproduce it, and what E₀ provenance does the fit imply? | **Clean inversion** — simplex fit of 13 `s_eff = 8` columns reproduces the histogram to **W₁ = 0.086 bins** (top 6 bins < 0.006); implied p(E₀) **smooth, unimodal, concentrated at/above E\* ≈ 0.46 eV**, in the RQ1 band without tuning; deep tail forces **≥ 9.6 % below the 0.22 eV floor** (droplet-R demand); **s_eff detector-identifiable via n = 1** (W9-P4); no multi-modal (branching) signature |
| 10. Steps 1+1.5: width decomposition + K forward model (2026-07-10) | 0 (re-read of all 41 gated dirs) | Which axis carries the histogram width? Can a narrow E₀ ≈ 0.28 eV + droplet-K spread source the 43.5 % bare peak — at which K-scale? | **Fixed-(E₀, droplet) columns are ≤ 1 bin wide** (W10-P1) — width is E₀/K-axis only; the (E₀, K) closed form is MD-validated (0.84 He); **probe K-scale excluded** for bare at narrow E₀ (needs δ ≈ 12, ~19× Kornilov — W10-P3); **production K-bracket lands bare 42.6 % / W₁ = 0.99 untuned** at (0.28 eV, δ = 0.8) — W10-P5 split verdict confirmed, **RQ7 co-requisite** |
| 11. RQ7 production-kinematics probe (2026-07-11) | 2 (first off-9 Å MD) | Measure K₂.₇₀ instead of bracketing it; does the closed form transfer? | **K₂.₇₀ = 0.746 measured** (ballistic bracket 0.49 wrong — drag eats the channel speed); 42 % beyond-band validity caveat (over-drag reading, VMI-corroborated ~2×); closed form transfers (0.24 He); narrow-E₀ landing re-calibrates to **E₀ ≈ 0.38–0.41 eV** pinned-droplet |
| H.2b analytic feasibility pass (2026-07-11) | 0 (zero-MD forward model) | Are the solvated targets (n₁ 31 %/n₂ 14.2 % + the mean-KE curve) reachable inside the all-bounded lever set (L1 position × L2 dressing × L3 floor ladder × priors × E₀ band)? | **Outcome (c) at the frozen bar — no cell passes T1–T4 at either exposure bracket**; but the miss is *structured*: bins 2–8 and the KE envelope/shape land; the residuals are exactly **n₁** (ladder-taper-controlled — the slid-X2 family *worsens* it; the RQ4-graded 2.2:1.5:1.3 diagnostic nearly lands, W₁ 0.575) and the **KE scale** (speed-selective ×1.3–1.5, the v_c direction); new **trapped droplet-retained class** 6–11 % |

**Common pins:** 9 Å case, `coulomb_available_eV = 0.80`, N = 50 single fixed
seed, 30 ps ion window, Tier-0 dt/drag bundle, `f_int = 0.5` (bridge pin —
timing-aligned: t× ≈ 5 ps on the anchored first shed), `f_ret = 0.1`,
`λ₀ = 0.9/ps`, E2 relaxation enabled. Wave 4 sets the relaxation cap to
1000 ps (a gated-off self-unbound complex never freezes — the 8.53·10⁶ ps
flight cap would not terminate).

**Standing wiring guarantees (all 87 probe dirs):** the (κ=1, mixture,
τ=6.55, s_eff=per-n, gate=none) point reproduces the Phase-D bridge exactly
(Δn̄ = 0.67, n_end = 20.33) — including after the cooling-gate build (the
`none`-arm byte-identity promise verified end-to-end); the 5-term
`ion_ledger_closure` residual is uniform ≈ 2.23·10⁻⁵ eV across every run.

**Anchored targets (Tier-1a comparator `build_shell_schedule`, t★ = 5 ps):**
7 sheds in-window → n_ion_end = 14, first shed at ≈ 5 ps.

---

## 1. Wave 1 — 45-run capability probe: the mechanism freezes at n ≈ 20

Grid: κ ∈ {0.5, 1, 2, 4, 8} × picture ∈ {`statistical_mixture`, `x2_only`,
`cooling_relaxed`} × τ ∈ {2.6, 6.55, 16.5} ps. **Outcome (b) of the probe
plan §6 — no in-band point lands:**

- Max reach **Δn̄ = 1.7 sheds** (κ=0.5, τ=16.5, `cooling_relaxed`) vs the
  anchored 7.
- Relaxed terminal n confined to **[19.3, 20.95]** across the whole grid —
  the pinned mechanism expresses essentially one shell.

### Structural findings (from the scoreboard + ladder/evaporation code)

1. **The κ lever is inverted** (corrects bridge findings §2 lever 1). Δn̄
   *decreases* monotonically with κ everywhere. The Form-U cliff centre is
   pinned at n\*+½ = 21.5 and the stripping range 21→14 lies *below* it:
   σ(21) = sigmoid(−κ/2) ≤ ½ always, so large κ drives the in-shell rungs
   *up* toward D₀(1). The κ→0 direction is also normalisation-capped —
   **Form-U structurally floors D₀(21) at ≈ 0.53·D₀(1)** (shallow minimum
   near κ ≈ 0.25). No κ makes the near-edge rungs cheap.
2. **The picture knob is near-degenerate for the staircase.** D₀(n) and Σ(n)
   both scale with (D₀(1) − D_floor), so the RRK argument
   x = D₀(21)/Σ(21) ≈ 0.032 is nearly picture-invariant — the three picture
   blocks of the scoreboard are almost identical.
3. **τ is the strongest but capped lever.** Δn̄ grows ≈ linearly over the band
   (0.2 → 1.7), but the gate-open time t× = τ·ln(f_int·E/Σ(21)) grows with τ;
   the ≈ 80 ps that 7 sheds would extrapolate to puts gate-open past the
   30 ps window. Dead end in-band and out.

### Kinetic diagnosis

k = ν·(1−x)^(s−1) with x ≈ 0.032 and s−1 = 59 gives k ≈ 0.15·ν ≈ 0.36/ps —
the freeze is **entirely the RRK exponent acting on a small x**. The swept
levers move x at most ~2× (linear); s acts exponentially. Energetics are
self-sustaining by construction (per shed, E_int and Σ(n) drop by the same
D₀(n): gate margin shed-invariant; 7 sheds cost 0.06 of the 0.188 eV
budget). This fired the **pre-registered RRK-dof mechanism-level OQ**
(bridge findings §2 lever 3): `s = 3n−3` treats the floppy quantum He₂₁
shell as 60 fully-coupled classical oscillators.

**Tier-2-blocking consequence:** the crossing construction caps the cascade
budget at Σ(21) at *any* `coulomb_available_eV` (budget moves only t×), and
the experimental abundance holds 43 % bare I⁺ with < 1 % weight at n ≥ 19 —
a mechanism frozen at n ≈ 20 cannot approach the arbitration observable at
0.80 **or** 2.70 eV. The probe saved the N=500 campaign from a dead band.

---

## 2. Wave 2 — `s_eff` mini-probe: the dof convention was the freeze

Grid: `s_eff ∈ {per-n, 30, 20, 12, 8, 5}` × τ ∈ {6.55, 16.5} ps at the
bridge pins (κ=1, mixture), via the already-plumbed `cfg.evap_rrk_dof`
override. **Outcome (a) of Addendum A.4 — a constant s_eff lands the
staircase in magnitude AND timing.** The τ = 6.55 ps arm:

| s_eff | Δn̄ | n_ion_end_mean | t_first_shed_ps | n_traj_MAD |
|---|---|---|---|---|
| per-n (60 at n=21) | 0.67 | 20.33 | 6.26 | 3.99 |
| 30 | 2.08 | 18.92 | 5.97 | 3.02 |
| 20 | 3.51 | 17.49 | 5.70 | 2.23 |
| 12 | 5.61 | 15.39 | 5.51 | 1.25 |
| **8** | **7.41** | **13.59** | **5.42** | **1.00** (grid argmin) |
| 5 | 9.35 | 11.65 | 5.37 | 2.41 |

The τ = 16.5 ps arm lands *magnitude only*: s_eff=20 reaches 5.93 sheds →
n_end 15.07 but with t_first ≈ 13.2 ps (timing killed by t× ∝ τ).
`frac_ions_shed = 1` at every set-s_eff point.

### Insights

1. **The kinetic diagnosis is quantitatively confirmed.** Dropping s from 60
   to 8 moves Δn̄ from 0.67 to 7.41 with *nothing else moving*.
2. **s↔τ separability works:** first-shed timing stays gate-open-governed
   (≈ 5.4–6.0 ps at τ=6.55; ≈ 12.9–13.5 ps at τ=16.5, tracking t×),
   magnitude is s-governed. **Landing region: s_eff ∈ ≈ [8, 12] × mid-band
   τ** — the anchored staircase alone nearly pins both, the best
   identifiability result of the program so far.
3. **The flexibility question is answered.** Relaxed terminal n spans
   [7.3, 20.9] across the sweep (monotone in s and τ; spreads 1.3–1.5 He).
   The mechanism was never inexpressive; the classical dof convention froze
   it. *(The report headline's per-picture ranges are a design artifact —
   only `statistical_mixture` received s_eff overrides in this wave.)*
4. **Constant s suffices at mid-band τ** (no in-window overshoot past 14 at
   s_eff=8): the n-dependent scaled convention `s = α·(3n−3)` (α ≈ 0.13) is
   *not demanded* — kept as a documented alternative arm, not a build.

**Physics interpretation.** The landing s_eff ≈ 8–12 says the effective heat
bath is ~an order of magnitude smaller than the 60 classical modes — the
*expected* direction for a cold, floppy, quantum He shell (most modes
quantum-frozen or uncoupled on the sub-ps shed timescale; cf.
evaporative-ensemble treatments of quantum clusters). The RRK form, energy
gate, ladder, and Σ(21) crossing budget all survive — the "gate/ladder
structure suspect" arm did **not** fire.

---

## 3. Wave 3 — picture cross-check: the landing band is picture-robust

10 runs: picture ∈ {`x2_only`, `cooling_relaxed`} × s_eff ∈ {30, 20, 12, 8,
5} at the landing arm (κ=1, τ=6.55). **Magnitude near-picture-invariance
CONFIRMED** — Δn̄ per picture (mixture / cooling_relaxed / x2_only):

| s_eff | Δn̄ (mix / cool / x2) | n_ion_end_mean (mix / cool / x2) |
|---|---|---|
| 30 | 2.08 / 2.03 / 2.02 | 18.92 / 18.97 / 18.98 |
| 20 | 3.51 / 3.34 / 3.22 | 17.49 / 17.66 / 17.78 |
| 12 | 5.61 / 5.59 / 5.33 | 15.39 / 15.41 / 15.67 |
| 8 | 7.41 / 7.46 / 7.13 | 13.59 / 13.54 / 13.87 |
| 5 | 9.35 / 9.11 / 9.10 | 11.65 / 11.89 / 11.90 |

- Shed-magnitude spread ≤ ~0.3 sheds (**≤ 4 %**) at every s_eff — the
  x = D₀(21)/Σ(21) picture-invariance prediction holds under reduced s.
- The one picture-sensitive read is **timing**: t_first at s_eff=8 is
  5.42 / 4.22 / 3.07 ps (mix / cool / x2; anchor 5 ps) — Σ(21) differs per
  picture, shifting t×. This is weak discrimination, not a picture
  selection: t× is **degenerate with the Bounded f_int** (a modest f_int
  shift re-aligns any picture's gate-open). Picture stays a near-flat co-fit
  dimension, entangled with f_int only through timing; its discrimination
  remains the size distribution's job.
- Trajectory-MAD argmin over all scored runs: (mixture, κ=1, τ=6.55,
  s_eff=8) at 1.00 He; both other pictures reach MAD ≈ 1.36 at s_eff=8.

**Consequence (user decisions, 2026-07-06):** s_eff promoted **Derived →
Bounded** (band ≈ [5, 20], landing prior [8, 12]; classical `s = 3n−3`
demoted to the classical-limit arm; `cfg.evap_rrk_dof` *is* the selection
surface — no new enum). F2 Stage 1 re-scoped from the κ×picture co-fit to an
**s_eff×τ co-fit** (18 runs, 9 Å / 0.80 eV; κ pinned with a sensitivity
spot-check; picture + f_int ride as a reported timing-degenerate pair).
Implemented 2026-07-06.

---

## 4. Wave 4 — cooling spatial gate + total-strip A/B (the last runs)

### 4.1 Motivation and arm

The experimental abundance reference peaks at **43 % bare I⁺**; the question
was whether the locked mechanism can express total (or near-total) stripping
at all. A structural trace found **no fixed-asymptote wall** (`E_inf(N)→0`
was built for OQ6 reachability); the limiter is **K2 Newton cooling**, which
drains E_int toward 0 and quenches shedding once E_int < D₀(n) — terminal n
floors at n ~ few, reached within a few τ. The user then flagged the real
asymmetry: **drag (γ ∝ ρ_He) and pickup (λ ∝ ρ_He) both gate off outside the
bubble via the shared erf-complement surface, but K2 cooling was applied
ungated everywhere** — including after ejection into vacuum, where there is
no droplet bath (the GAH25 τ was fit for a near-droplet *growing* Na⁺
shell). Hence the new interchangeable arm (`cfg.cooling_spatial_gate`,
default `none` byte-identical):

- `none`: Ė_int|K2 = −E_int/τ (status quo).
- `density_scaled`: Ė_int|K2 = −(ρ_He/ρ_bulk)·E_int/τ — cooling → 0 outside
  the bubble, sharing the drag/pickup gate boundary (no new steepness knob).

Grid (24 runs): gate ∈ {`none`, `density_scaled`} × s_eff ∈ {1, 2, 3, 5}
(s_eff = 1 is the **max-kinetics bound**: the RRK factor is identically 1) ×
τ ∈ {6.55, 16.5, 30} ps, at κ=1 / mixture / probe pins; relaxation cap
1000 ps.

### 4.2 Results (scored 2026-07-07; mixture, N = 50)

**τ = 6.55 ps arm — the gate turns shell-retention into near-total strip:**

| s_eff | gate | Δn̄ | n_ion_end | t_first [ps] | n_relax_mean | spread | n_relax_min | frac_frozen |
|---|---|---|---|---|---|---|---|---|
| 1 | none | 12.46 | 8.54 | 5.31 | 8.54 | 1.35 | 5 | 1.00 |
| 1 | density_scaled | 18.29 | **2.71** | 6.68 | 2.71 | 0.50 | **2** | 1.00 |
| 2 | none | 11.66 | 9.34 | 5.33 | 9.34 | 1.34 | 6 | 1.00 |
| 2 | density_scaled | 18.17 | 2.83 | 6.68 | 2.73 | 0.49 | 2 | 1.00 |
| 3 | none | 10.89 | 10.11 | 5.36 | 10.11 | 1.38 | 7 | 1.00 |
| 3 | density_scaled | 18.01 | 2.99 | 6.69 | 2.84 | 0.42 | 2 | 0.98 |
| 5 | none | 9.35 | 11.65 | 5.37 | 11.65 | 1.48 | 8 | 1.00 |
| 5 | density_scaled | 17.52 | 3.48 | 6.69 | 3.02 | 0.35 | 2 | 0.89 |

**τ = 16.5 and 30 ps arms — the gated arm sheds NOTHING:**

| τ [ps] | s_eff | gate | Δn̄ | n_ion_end | t_first [ps] | n_relax_mean | n_relax_min | frac_frozen |
|---|---|---|---|---|---|---|---|---|
| 16.5 | 1 | none | 16.31 | 4.69 | 12.88 | 4.69 | 2 | 1.00 |
| 16.5 | 2 | none | 15.74 | 5.26 | 12.88 | 5.26 | 3 | 1.00 |
| 16.5 | 3 | none | 15.02 | 5.98 | 12.90 | 5.96 | 3 | 1.00 |
| 16.5 | 5 | none | 13.62 | 7.38 | 12.90 | 7.34 | 4 | 1.00 |
| 16.5 | 1–5 | density_scaled | **0.00** | **21.00** | — | 21.00 | 21 | **0.00** |
| 30 | 1 | none | 16.20 | 4.80 | 23.16 | 3.13 | 2 | 1.00 |
| 30 | 2 | none | 15.02 | 5.98 | 23.17 | 3.51 | 2 | 1.00 |
| 30 | 3 | none | 13.79 | 7.21 | 23.18 | 4.01 | 3 | 1.00 |
| 30 | 5 | none | 11.89 | 9.11 | 23.21 | 5.16 | 4 | 1.00 |
| 30 | 1–5 | density_scaled | **0.00** | **21.00** | — | 21.00 | 21 | **0.00** |

Ledger residual stays uniform 2.23·10⁻⁵ eV across all 24 runs (5-term
closure holds on both arms — E_dissip books the actual gated drain).

### 4.3 Mechanism trace of the τ cliff (verified in the trajectories)

The gated arm is governed by a **race between the gate-open crossing and
droplet ejection**. With f_int·E = 0.40 eV and Σ(21) ≈ 0.188 eV (mixture),
the ungated crossing time is t× = τ·ln(0.40/0.188) ≈ 0.757·τ → 4.96 / 12.5 /
22.7 ps at τ = 6.55 / 16.5 / 30 — matching the measured ungated first sheds
(5.31 / 12.88 / 23.16 ps). The ion exits the ~27.9 Å droplet at **~6 ps**
(r̄ = 25.7 Å at 5 ps, 35.7 Å at 8 ps; drag, pickup, and now cooling all gate
off there):

- **τ = 6.55 (crossing ≈ ejection):** cooling completes the crossing just as
  the ion leaves (first shed 6.68 ps, ~1.3 ps later than ungated — cooling
  slows during exit). After ejection the quench is gone: the self-bound
  complex sheds freely down to n̄ ≈ 2.7 by ~20 ps (E_int held ≈ 5 meV,
  below D₀ → frozen). In-window near-total strip.
- **τ ≥ 16.5 (crossing after ejection):** cooling shuts off at ejection with
  E_int held at **0.280 eV > Σ(21) ≈ 0.188 eV** — the self-bound gate
  **never opens**. Zero sheds; n = 21 flat through the 30 ps window *and*
  the full 1000 ps relaxation; `frac_frozen = 0.0` (not converged — held
  above the gate forever, terminated only by the cap). This is exactly the
  freeze-termination loss the code review flagged: under `density_scaled`,
  "E_int monotone non-increasing" does **not** imply termination; the finite
  cap (plus the config step-budget guard) is what bounds it.

### 4.4 Findings

1. **Near-total stripping is expressible — with the gate, at mid-band τ.**
   `density_scaled` at τ = 6.55 strips the whole ensemble to n̄ ≈ 2.7–3.5
   (relaxed min n = 2, tight spreads 0.35–0.50, `frac_ions_shed` = 1). The
   post-ejection K2 quench was indeed the shell-retention limiter, exactly
   as the structural trace predicted.
2. **The ungated arm confirms the shell-retaining floor.** Even at the
   max-kinetics bound (s_eff = 1, RRK factor ≡ 1), `none` at mid-band τ
   floors at n̄ ≈ 8.5; longer τ digs deeper (n̄ ≈ 4.7 at 16.5 ps; relaxed
   n̄ ≈ 3.1 at τ = 30 with the cascade continuing into the relaxation
   stage) but at the price of first sheds at 13–23 ps — an order of
   magnitude off the anchored 5 ps. Ungated deep stripping and staircase
   timing are irreconcilable.
3. **The gated arm is all-or-nothing in τ** (a structural cliff, not a
   tuning artifact): deep strip if t× ≲ t_eject (~6 ps), frozen at n = 21
   forever if t× > t_eject. There is no intermediate. Two consequences:
   - the gated arm's viability is controlled by the **t×↔ejection race**,
     i.e. by (τ, f_int, Σ(21)/picture) jointly — the f_int degeneracy now
     has a *hard* consequence, not just a timing shift;
   - a residual out-of-bubble cooling floor (`rho_min`, the deferred OQ)
     would soften the cliff into a slow late opening — the A/B brackets the
     two extremes (`none` = floor 1, `density_scaled` = floor 0).
4. **Bare I⁺ (n = 0) is NOT reached anywhere in the grid.** The deepest
   shell is **n = 2** (gated τ = 6.55, all s_eff; ungated only at s_eff ≤ 2
   with τ ≥ 16.5). The floor is energetic: the Σ(21)-crossing construction
   hands the cascade exactly Σ(21) at gate-open, so any drain after the
   crossing (residual in-bubble cooling between sheds, the f_ret retention
   loss) leaves it short of the full 21-rung cost; E_int lands at a few meV
   < D₀(n) around n ≈ 2–3 and freezes (frac_frozen = 1.0 at s_eff ≤ 2).
   **Expressing the experimental 43 % bare-I⁺ peak therefore remains open**
   at 0.80 eV / f_int = 0.5 — a finding for the F5/production discussion,
   not a retune target.
5. **`frac_frozen` earned its place as a first-class column.** It cleanly
   separates converged terminal reads (1.0), near-threshold RRK slow-down at
   the cap (0.98 / 0.89 at gated s_eff = 3 / 5 — larger s is slower near
   threshold), and the never-opened gate (0.0) — without it the τ ≥ 16.5
   gated rows would silently read as converged n = 21.
6. **The total-strip corner and the staircase corner are disjoint.** The
   staircase argmin (ungated, s_eff ≈ 8, τ = 6.55, MAD 1.00) sheds 7; every
   deep-strip point destroys staircase fidelity (gated τ = 6.55 MAD
   5.97–7.53; ungated long-τ first sheds 13–23 ps). No single (s_eff, τ,
   gate) point does both.
7. **Every probe point produces a narrow terminal distribution** (spreads
   0.35–1.5 He), while the experimental abundance spans bare through n ≳ 20.
   Whatever reproduces the broad experimental shape must come from ensemble
   heterogeneity and/or the regime axis, not from a single knob point's
   intrinsic spread — an identifiability note for F4/F5, reported not
   adjudicated.

---

## 4b. Wave 5 — gated landing re-location: the landing moves to s_eff ≈ 30, and the terminal decouples

Grid (5 runs): s_eff ∈ {8, 12, 16, 20, 30} at τ = 6.55 ps,
`cooling_spatial_gate = density_scaled`, probe pins (κ=1, mixture,
f_int = 0.5, 9 Å / 0.80 eV / N=50, 1000 ps relaxation cap). Motivation: the
Addendum-A staircase landing s_eff ∈ [8, 12] was established **ungated**;
Wave 4 showed the gate removes the post-ejection quench and roughly doubles
in-window shedding at s_eff ≤ 5, so the ungated prior could not be assumed to
carry to a gated production arm. Wiring: the gated s_eff=5 point re-scored
unchanged from Wave 4 (Δn̄ = 17.52); ledger residual uniform 2.23·10⁻⁵ eV.

**Gated τ = 6.55 ps sweep (`density_scaled`), full bracket incl. Wave-4 s ≤ 5:**

| s_eff | Δn̄ | n_ion_end | t_first [ps] | MAD | n_relax_mean | n_relax_min | frac_frozen |
|---|---|---|---|---|---|---|---|
| 1 | 18.29 | 2.71 | 6.68 | 7.53 | 2.71 | 2 | 1.00 |
| 2 | 18.17 | 2.83 | 6.68 | 7.20 | 2.73 | 2 | 1.00 |
| 3 | 18.01 | 2.99 | 6.69 | 6.81 | 2.84 | 2 | 0.98 |
| 5 | 17.52 | 3.48 | 6.69 | 5.97 | 3.02 | 2 | 0.89 |
| 8 | 16.16 | 4.84 | 6.71 | 4.70 | 3.11 | 3 | 0.88 |
| 12 | 14.30 | 6.70 | 6.80 | 3.25 | 3.86 | 3 | 0.17 |
| 16 | 12.46 | 8.54 | 6.90 | 2.04 | 4.20 | 3 | 0.01 |
| 20 | 10.62 | 10.38 | 6.99 | 1.13 | 4.79 | 4 | 0.00 |
| **30** | **7.20** | **13.80** | **7.28** | **0.96** | 6.21 | 5 | **0.00** |

**Ungated τ = 6.55 ps sweep (`none`) for contrast** (the Wave-2 landing at
s_eff ≈ 8; here `n_relax` = `n_ion_end` exactly — frac_frozen = 1.0 throughout):

| s_eff | Δn̄ | n_ion_end | t_first [ps] | MAD | n_relax_mean | frac_frozen |
|---|---|---|---|---|---|---|
| 1 | 12.46 | 8.54 | 5.31 | 4.76 | 8.54 | 1.00 |
| **8** | **7.41** | **13.59** | **5.42** | **1.00** | 13.59 | 1.00 |
| 30 | 2.08 | 18.92 | 5.97 | 3.02 | 18.92 | 1.00 |
| per-n | 0.67 | 20.33 | 6.26 | 3.99 | 20.33 | 1.00 |

### Findings

1. **The gated staircase landing is at s_eff ≈ 30, not 8** — prediction
   confirmed. The gate removes the post-ejection quench, so to hold the
   in-window cascade to the anchored 7 sheds (n_ion_end 14) you need slower
   near-threshold kinetics: gated s_eff=30 gives Δn̄ 7.20, n_ion_end 13.80,
   MAD 0.96 (grid best). The gate multiplies the effective in-window shedding
   by ≈ 3.75× (ungated s_eff=8 ↔ gated s_eff=30 both land). **Running the F2
   campaign on the ungated [8,12] prior under a gated production arm would
   have missed badly** (gated s_eff=8 → n_ion_end 4.84, MAD 4.70) — exactly
   the mistake Wave 5 existed to prevent.
2. **The gated in-window staircase and the terminal read DECOUPLE — the
   central Wave-5 result.** Ungated, the cascade quenches by sim-end so
   n_ion_end = n_relaxed exactly (frac_frozen = 1.0 at every s_eff): one read.
   Gated, at the staircase-landing s_eff=30, frac_frozen = 0.00 — the
   relaxation ran the full 1000 ps cap without a single ion freezing, n
   drifting 13.80 (30 ps) → 6.27 (1000 ps), E_int 0.104 → 0.034 eV, slope
   still −0.0012 n/ps and decelerating logarithmically. **s_eff sets the
   *rate* of the post-ejection cascade, not its endpoint.** The endpoint is
   the same energetic floor n ≈ 2 for all s_eff (reached in-window at s_eff=1,
   frozen at 30 ps; reached only after ≫ 1000 ps at s_eff=30). Confirmed by
   the relaxation lengths: gated s_eff=1 relaxation stored **2 timesteps**
   (froze immediately); gated s_eff=30 stored the **full 50001** (never froze).
3. **Consequence — the gated terminal size distribution is flight-time
   dependent at the staircase-landing s_eff.** The observable that arbitrates
   Tier 2 is the *terminal* distribution, but under the gate at s_eff ≈ 30 the
   terminal is not converged within a tractable window: the physical ~8.5 µs
   flight is 8500× the 1000 ps cap, and the cascade is still live. Either the
   flight time must be resolved (intractable at this dt), or the terminal is
   cap-artifactual, or the gated s_eff must be pushed even higher so the
   cascade freezes shell-retaining in-window (but that breaks the staircase).
   This is a genuine model prediction, not a numerical nuisance: a gated
   complex with a positive self-bound margin G = E_int − Σ(n) at ejection
   keeps evaporating in vacuum until G closes, and the closing rate is what
   s_eff controls.
4. **Timing lands late by ~2.3 ps under the gate** (Addendum B.4 outcome (b)).
   Every gated point sheds first at 6.7–7.3 ps vs the anchored t★ = 5 ps: the
   gated crossing lags because cooling slows as the complex exits the bubble,
   extending t×. A modest f_int increase re-aligns it (t× = τ·ln(f_int·E/Σ));
   this is the campaign's f_int dimension doing its job, not a miss.
5. **The Σ(21)-crossing energetic floor is unchanged by the gate** (n ≈ 2,
   not bare) — the gate changes *how fast and how completely* the cascade
   reaches the floor within a window, not the floor itself. Wave-4 finding 4
   (bare I⁺ unexpressed) stands.

### What Wave 5 means for the arm choice

The clean result of Waves 4+5 together: **the two arms are terminal-regime
opposites.** Ungated freezes fast → a *shell-retaining* terminal at the
s_eff-selected n (well-defined, converged). Gated does not freeze at the
staircase-landing s_eff → a *deep-stripping* terminal heading to the n ≈ 2
floor over the flight (rate-limited, flight-time-dependent). Neither pure arm,
at a single knob point, produces the **broad bimodal** experimental
distribution (43 % bare *and* substantial shell-retaining weight to n ≳ 20).
That reinforces Wave-4 finding 7: the experimental shape must come from
**ensemble heterogeneity across the t×↔ejection race** (some trajectories
gate-open before ejection and strip deep; some cross after ejection and retain
their shell) and/or a genuine regime distribution — not from any one
(gate, s_eff, τ, f_int) point. The gate is what *makes both sides reachable*;
the spread across them is the campaign's real target.

---

## 4c. Physical interpretation — what s_eff really is (the cooling race)

Synthesis of the post-Wave-5 discussion (2026-07-07). This resolves the
recurring "won't high s just freeze at high n again?" question and pins down
the physical meaning of each swept knob. Grounded in `evaporation.py`: the
shed band is `D_0(n) < E_int < Σ(n)`, each shed drains `dE_int = −D_0(n)`
(K1) and drops `n → n−1`, so `Σ(n)` falls by the same `D_0(n)`.

**s_eff is a rate knob, not an amount knob.** It enters only the RRK rate
`k = ν·(1 − D_0(n)/E_int)^(s−1)` — the *exponent of the evaporation rate*. It
appears in no energy balance, so it has **no direct effect on terminal n**.
Every apparent effect of s on "how many He come off" is indirect and runs
through a race.

**What drains E_int, and what closes the gate:**
- **Evaporation** drains `−D_0(n)` per shed, but `Σ(n)` drops by the same
  amount, so the self-bound margin `G = E_int − Σ(n)` is **shed-invariant**.
  Evaporation alone never closes its own gate — it is self-sustaining and, left
  alone, strips to the energetic floor where `E_int < D_0(n)`.
- **Cooling (K2)** drains `−E_int/τ` continuously and is `Σ`-blind, so it
  drives `G` negative. **Cooling is the only thing that closes the gate at
  high n.**

**Terminal n is therefore an evaporation-vs-cooling race** (rate `k` vs `1/τ`;
what matters is roughly `k·τ`). High s → slow evaporation → cooling zeroes
`E_int` while n is still high → energetic freeze at high n. Low s → fast
evaporation → many sheds before cooling wins → freeze at low n. This is why s
*looked* like an amount knob ungated, and why s and τ are entangled ungated
(longer τ strips deeper at fixed s; higher s strips shallower at fixed τ — both
in the data, §4b / Wave 2).

**The two freezes are physically different (answers the recurring question):**
- **Ungated high-s, high-n = an *energetic* freeze.** Cooling closed the gate;
  `E_int < D_0(n)`; it is done for all time. Robust.
- **Gated high-s, high-n = a *kinetic pause*.** After ejection there is no
  cooling competitor, `G` reverts to shed-invariant, the gate stays
  energetically open — `E_int` is still above `D_0(n)`, the cascade is merely
  slow. It has *not* frozen; it is still descending to the same low floor, just
  not within the window. Fragile — purely a function of wait time.

Data proof of the distinction: gated s=1 stored **2** relaxation timesteps
(froze at n≈2 instantly); gated s=30 stored the **full 50001** (never froze,
n 13.8 → 6.27, still shedding). Same energetic destination (~n≈2, set by the
brief in-bubble/crossing cooling, roughly s-independent); wildly different
arrival time. **Gated, s does not set the terminal n — it sets the terminal n
*at a given detection time*.**

**Per-knob roles (unified):**

| Knob | Physical role | Sets terminal n? |
|---|---|---|
| **s_eff** | evaporation rate (RRK exponent); pure kinetics | only via the cooling race → **yes ungated, no gated** (only the descent rate) |
| **τ** | cooling rate `1/τ`; also pre-crossing timing `t× = τ·ln(f_int·E/Σ)` | ungated yes (it *is* the competitor); gated only via `t×` |
| **f_int** | `E_int` at gate-open vs `Σ` → timing (`t×`) only; never the budget | no |
| **cooling gate** | whether the cooling competitor persists after ejection | it is the *switch* that flips s from an endpoint-lever (ungated) to a rate-only-lever (gated) |

**One sentence:** cooling converts internal energy into a *closed gate*;
evaporation alone is self-sustaining and always strips to the floor; s only
changes evaporation's *speed*, so s affects the endpoint only for as long as
cooling is present to race against. The gate removes that race after ejection
— which is precisely why it changes what s *means*. For understanding the
simulation, the **cooling gate is the more fundamental knob than s**: it
decides whether s is an outcome parameter or a clock parameter.

### On the staircase itself — a capability showcase, not a meaningful anchor

**The literal reproduction of the 21→19→14 staircase carries no physical
meaning as a target.** It only demonstrates *capability*: that by choosing
parameters within this 30 ps timeframe the mechanism *can* express that
in-window shedding profile. It says nothing about the physically relevant
outcome, because — as §4b makes explicit — we do not know what happens after
the window: under the gate the cascade is still descending (flight-time
dependent), and even the "landing" s_eff ≈ 30 is a snapshot of a live
process, not a converged terminal. The staircase is therefore a
**showcase/existence check** (the mechanism is expressive enough), never a
calibration anchor. The only anchor that carries meaning is the experimental
**terminal** size distribution, read at an explicit, physically-motivated
detection time (the production-run stance, B.1/§4b caveat b). TDDFT is not
ground truth here in a second, sharper sense than the standing 9 Å
non-radial flag: it constrains the *in-window* trace, and the in-window trace
is not the observable.

---

## 4d. Wave 6 — the f_int probe: under the gate, f_int is a budget knob, not a timing nicety

Grid (10 new runs, executed 2026-07-07 under the trigger; Addendum C.1;
scratchpad-driver route, zero repo-code change): `density_scaled`, τ = 6.55:
s_eff ∈ {8, 30} × f_int ∈ {0.24, 0.30, 0.35, 0.65}; τ = 16.5: s_eff ∈ {8, 30}
× f_int = 0.25; probe pins otherwise. The f_int = 0.50 companions (Waves 4/5)
are the wiring oracle — all three re-scored exactly (s_eff = 5/8/30);
ledger residual uniform ≈ 2.23·10⁻⁵ eV over all **102** dirs.

**Gated τ = 6.55 ps, f_int sweep (`density_scaled`; f_int = 0.50 rows = Wave 4/5):**

| s_eff | f_int | Δn̄ | n_ion_end | t_first [ps] | MAD | n_relax_mean | n_relax_min | frac_frozen |
|---|---|---|---|---|---|---|---|---|
| 8 | 0.24 | 9.08 | 11.92 | 0.68 | 2.93 | 10.96 | 8 | 0.51 |
| 8 | 0.30 | 10.72 | 10.28 | 2.20 | 3.56 | 9.13 | 7 | 0.56 |
| 8 | 0.35 | 12.12 | 8.88 | 3.22 | 4.03 | 7.66 | 7 | 0.64 |
| 8 | 0.50 | 16.16 | 4.84 | 6.71 | 4.70 | 3.11 | 3 | 0.88 |
| 8 | 0.65 | **0.00** | **21.00** | — | 4.50 | 21.00 | 21 | 0.00 |
| 30 | 0.24 | 2.87 | 18.13 | 1.84 | 3.00 | 14.88 | 13 | 0.00 |
| 30 | 0.30 | 3.75 | 17.25 | 3.22 | 2.50 | 12.85 | 11 | 0.00 |
| 30 | 0.35 | 4.63 | 16.37 | 3.88 | 2.10 | 11.20 | 10 | 0.00 |
| 30 | 0.50 | 7.20 | 13.80 | 7.28 | **0.96** | 6.21 | 5 | 0.00 |
| 30 | 0.65 | **0.00** | **21.00** | — | 4.50 | 21.00 | 21 | 0.00 |

**Gated τ = 16.5 ps at f_int = 0.25 (the Wave-4 dead arm):**

| s_eff | Δn̄ | n_ion_end | t_first [ps] | MAD | n_relax_mean | n_relax_min | frac_frozen |
|---|---|---|---|---|---|---|---|
| 8 | 14.17 | 6.83 | 1.61 | 5.91 | 5.68 | 4 | 0.48 |
| 30 | 6.35 | **14.65** | 2.27 | 1.46 | 8.55 | 7 | 0.00 |

(Gated `n_relax` rows with `frac_frozen < 1` are cap-truncated snapshots —
the §6 I11 boundary applies throughout.)

### Prediction verdicts (pre-registered in Addendum C.1)

- **P1 (timing re-alignment) — qualitatively CONFIRMED, quantitatively
  refined.** The gated first shed is strongly f_int-tunable and monotone
  (0.7 → 7.3 ps over the sampled band — it spans the anchored t★ = 5 ps).
  But f_int = 0.35 lands 3.88 ps, not the predicted ≈ 5: the anchored
  timing sits at **f_int ≈ 0.40–0.42**, not 0.35. The miss localizes an
  assumption: the ~2 ps gated crossing *lag is not f_int-independent* — it
  grows with t× (lag ≈ 1.3 ps at f_int = 0.35 vs ≈ 2.3 ps at 0.50), because
  a crossing closer to ejection happens under already-weakening cooling.
- **P2 (cliff, closed side) — CONFIRMED exactly.** f_int = 0.65 at τ = 6.55
  sheds nothing at either s_eff (n = 21 flat, frac_frozen = 0: held above
  the gate forever). The in-band cliff sits between f_int 0.50 and 0.65 at
  τ = 6.55, consistent with the derived edge ≈ 0.59 minus the exit-lag
  correction.
- **P3 (re-opening the dead arm) — CONFIRMED.** f_int = 0.25 re-opens
  τ = 16.5 (zero sheds at 0.50 in Wave 4): s_eff = 8 deep-strips
  (n_end 6.83), and s_eff = 30 lands **near-staircase magnitude**
  (n_end 14.65, MAD 1.46) with early timing — the race demonstrated from
  the f_int side, and a second staircase-magnitude region exposed at long τ.
- **P4 (floor depth vs the in-bubble leak) — CONFIRMED, with the sign made
  vivid.** The relaxed floor is monotone in f_int: *earlier* gate-open
  (lower f_int) means a **shallower** strip (s_eff = 8: n̄ 10.96 at
  f_int = 0.24 vs 3.11 at 0.50; min-n 8 vs 3), because the cascade runs
  in-bubble under active cooling and the leak over [t×, t_eject] eats the
  Σ(21) budget. The deepest strip sits where the crossing hugs ejection —
  the sampled optimum is f_int = 0.50 (t× ≈ t_eject), and the trend says
  the true optimum hugs the cliff edge from below. **Bare I⁺ is not reached
  at any f_int**: even the minimum-leak corner floors at n ≈ 2–3, and every
  earlier opening makes it *worse*. The OQ-B gap is now quantified from the
  f_int side: within the Σ(21)-crossing construction, no (s_eff, τ, f_int,
  gate) reaches n = 0 at 0.80 eV.

### Structural finding — f_int's role changes under the gate

The standing "f_int is timing-only (the budget is Σ(21) by the crossing
definition)" note is **ungated-only**. Under `density_scaled`, f_int sets
*how much of the crossing budget survives to ejection*: the budget at
crossing is always Σ(21), but the in-bubble cooling drain over [t×, t_eject]
scales with how early the gate opens. Gated, f_int is therefore a **joint
timing + effective-budget knob** — it cannot re-align staircase timing
without also moving magnitude (s_eff = 30 / f_int = 0.35: timing improves to
3.9 ps but n_end rises to 16.4). Landing both at τ = 6.55 needs a *joint*
(s_eff, f_int) co-fit (rough interpolation: s_eff ≈ 25, f_int ≈ 0.42) —
exactly the campaign's Stage-1 shape, now with a quantified prior.

---

## 4e. Wave 7 — the detected read: s_eff compresses away, the race margin owns the detector

Executed 2026-07-08 (under the `[PROCEED TO IMPLEMENTATION]` trigger;
scratchpad driver through the delivered Slice-DS `run_detection_stage` —
zero repo-code change, **zero new MD runs**): every one of the **102** probe
dirs re-read at the Sourced detector time **t_detect = 8.53 µs**
(CALIBRATION_MAP row 24), seeding from its existing `relaxation.npz`;
`detection.npz` written per dir; full re-score with the delivered report
(`n_detect_*` + state-reason columns populated on all 102 rows).

**Execution note (recorded, not a code change).** The 63 pre-Wave-4 dirs are
stamped with the legacy relaxation cap 8.53·10⁶ ps (chosen before Slice DS
existed), so the config-load *nominal*-window bound rejects
t_detect ≤ 30 + 8.53·10⁶ ps even though their E2 runs terminated at
all-frozen after tens of ps. The transient per-dir cfg view therefore carried
the **realized** relaxation duration (read from the artifact) alongside the
detection fields, so the nominal and realized bounds state the same physical
fact; on-disk `cfg.json` was never modified, and the stage's realized-t_h
re-check and the P1–P3 handover guard ran on every dir — **zero violations**
(as predicted: all ions sit at erfc-underflowed ρ_He at handover).

**Wiring:** ungated no-op identity exact on all **75** ungated dirs
(`n_detect ≡ n_relaxed` elementwise, zero events drawn, all ions `frozen`) —
W7-P1 confirmed; the bridge point re-scores unchanged (Δn̄ 0.67 /
n_ion_end 20.33 / n_relaxed 20.33 / **n_detect 20.33**); ledger residual
uniform ≈ 2.23·10⁻⁵ eV across all 102 rows.

### Gated τ = 6.55 / f_int = 0.50 — the s_eff sweep at the detector

(`n_relax` = the 1000 ps cap read of §4b, for contrast; `frozen`/`texh` =
detector state-reason fractions, `texh` = `time_exhausted` = cascade live at
the detector.)

| s_eff | n_relax (cap) | n̄_detect | detect min–max | frozen | texh |
|---|---|---|---|---|---|
| 1 | 2.71 | 2.71 | 2–4 | 1.00 | 0.00 |
| 2 | 2.73 | 2.73 | 2–4 | 1.00 | 0.00 |
| 3 | 2.84 | 2.82 | 2–4 | 1.00 | 0.00 |
| 5 | 3.02 | 2.92 | 2–4 | 0.99 | 0.01 |
| 8 | 3.11 | 3.05 | 2–4 | 0.94 | 0.06 |
| 12 | 3.86 | 3.08 | 3–4 | 0.95 | 0.05 |
| 16 | 4.20 | 3.18 | 3–4 | 0.86 | 0.14 |
| 20 | 4.79 | 3.79 | 3–4 | 0.24 | 0.76 |
| **30** | 6.21 | **4.09** | 4–5 | **0.00** | **1.00** |

### Gated f_int sweep at the detector (the Wave-6 grid re-read)

| τ [ps] | s_eff | f_int | n_relax (cap) | n̄_detect | detect min–max | frozen | supp | texh |
|---|---|---|---|---|---|---|---|---|
| 6.55 | 8 | 0.24 | 10.96 | 10.60 | 8–13 | 0.87 | 0 | 0.13 |
| 6.55 | 8 | 0.30 | 9.13 | 8.79 | 7–11 | 0.90 | 0 | 0.10 |
| 6.55 | 8 | 0.35 | 7.66 | 7.37 | 6–9 | 0.93 | 0 | 0.07 |
| 6.55 | 8 | 0.50 | 3.11 | 3.05 | 2–4 | 0.94 | 0 | 0.06 |
| 6.55 | 8 | 0.65 | 21.00 | 21.00 | 21 | 0 | 1.00 | 0 |
| 6.55 | 30 | 0.24 | 14.88 | 13.02 | 12–14 | 0.00 | 0 | 1.00 |
| 6.55 | 30 | 0.30 | 12.85 | 10.89 | 10–12 | 0.00 | 0 | 1.00 |
| 6.55 | 30 | 0.35 | 11.20 | 9.21 | 8–11 | 0.00 | 0 | 1.00 |
| 6.55 | 30 | 0.50 | 6.21 | 4.09 | 4–5 | 0.00 | 0 | 1.00 |
| 6.55 | 30 | 0.65 | 21.00 | 21.00 | 21 | 0 | 1.00 | 0 |
| 16.5 | 8 | 0.25 | 5.68 | 5.27 | 4–7 | 0.89 | 0 | 0.11 |
| 16.5 | 30 | 0.25 | 8.55 | 6.64 | 6–8 | 0.00 | 0 | 1.00 |

The ten gated dead arms (τ ∈ {16.5, 30} × s_eff ∈ {1, 2, 3, 5} at
f_int = 0.50; τ = 6.55 × s_eff ∈ {8, 30} at f_int = 0.65) all arrive **100 %
`suppressed` at exactly n = 21** — W7-P3 confirmed; the delivered gate
semantics ride to the detector verbatim (OQ-B's fragmentation question is now
literal weight in the observable, still unresolved by construction).

### Detection-time sensitivity band (log-spaced, from stored event times)

n̄(t) at t = 10³ / 10⁴ / 10⁵ / 10⁶ / 8.53·10⁶ ps — a pure re-read of
`detection.npz` (design §3.3); linear fractions of t_detect would have shown
three near-identical numbers, the descent is logarithmic:

| τ [ps] | s_eff | f_int | 10³ | 10⁴ | 10⁵ | 10⁶ | 8.53·10⁶ |
|---|---|---|---|---|---|---|---|
| 6.55 | 30 | 0.50 | 6.21 | 5.31 | 4.94 | 4.43 | 4.09 |
| 6.55 | 20 | 0.50 | 4.79 | 4.17 | 3.99 | 3.97 | 3.79 |
| 6.55 | 16 | 0.50 | 4.20 | 4.00 | 3.77 | 3.44 | 3.18 |
| 6.55 | 30 | 0.35 | 11.20 | 10.36 | 9.83 | 9.48 | 9.21 |
| 6.55 | 30 | 0.24 | 14.88 | 14.03 | 13.66 | 13.31 | 13.02 |
| 6.55 | 8 | 0.24 | 10.96 | 10.79 | 10.70 | 10.64 | 10.60 |
| 16.5 | 30 | 0.25 | 8.55 | 7.77 | 7.20 | 6.92 | 6.64 |

### Findings

1. **W7-P2 refuted — in the informative direction.** The gated staircase
   landing (s_eff = 30, τ = 6.55, f_int = 0.50) does **not** reach the n ≈ 2
   energetic floor by the detector: it arrives **100 % `time_exhausted` at
   n̄ = 4.09** (91 % at n = 4, 9 % at n = 5), still descending at
   ≈ 0.4–0.9 He per *time decade* (band row 1). At s_eff = 30 the floor is
   asymptotic — effectively unreachable at any laboratory flight time. I11's
   "flight-time dependent" now has a number, and the Wave-5 open end is
   resolved: the gated landing point contributes **deep-strip-adjacent weight
   (n ≈ 4), not mid-shell weight**.
2. **The detector read compresses s_eff — the headline.** At f_int = 0.50 the
   whole s_eff ∈ [1, 30] range lands n̄_detect ∈ [2.7, 4.1], while the same
   points span 2.7–13.8 in-window and 2.7–6.2 at the 1000 ps cap. What s_eff
   decides at the detector is the **arrival state** (frozen by s ≤ 8 vs live
   at s = 30; the frozen fraction falls 1.00 → 0.00 over s_eff 1 → 30), not
   the arrival n. This is I13 measured at the physical read: s is a clock,
   and 8.53 µs is long enough that the clock mostly runs out.
3. **Identifiability consequence (sharpens I18).** The two reads are
   near-orthogonal: the in-window staircase is s_eff-sensitive and
   race-margin-degenerate; the detector distribution is race-margin- and
   leak-sensitive (f_int at fixed τ moves n̄_detect 3.05 → 10.60 at s_eff = 8)
   and nearly s_eff-blind. A campaign co-fit gets s_eff from the staircase
   prior and the race margin Δ× from the terminal distribution — the two
   targets constrain different knobs, which is the good case.
4. **W7-P4 confirmed — the gated mid-shell channel is real and converged.**
   The early-gate-open leak points arrive mostly *energetically frozen* at
   mid-shell (s_eff = 8: n̄_detect 10.60 / 8.79 / 7.37 at f_int = 0.24 / 0.30
   / 0.35 with frozen fractions 0.87–0.93). Combined with the deep-strip
   corner (n ≈ 3–4) and the suppressed class (n = 21), **the gated arm alone
   spans n_detect ≈ 3–13 continuously along the race margin, plus the n = 21
   class** — the ensemble-heterogeneity target (I12/I18) is expressible
   within a single arm at the detector.
5. **Bare I⁺ is unreached at the detector** (global min n_detect = 2, and
   only via the fast-kinetics s_eff ≤ 5 corner). OQ-B stands at the
   physically meaningful read; at s_eff ≳ 20 the floor is *additionally*
   kinetically out of reach (finding 1), so flight time cannot resolve it
   either.

### Physical interpretation (post-Wave-7 discussion, 2026-07-08) — why nothing freezes, why s_eff compresses, and where the two cracks are

Synthesis of the post-Wave-7 discussion (user: "I expected all to freeze").
The result is simultaneously *textbook* — in a way that partially validates
the mechanism's long-time structure — and *constructed* — in a way that
localizes the remaining qualitative gaps in two bookkeeping conventions,
not in any swept knob.

**1. The post-ejection cascade is an exactly closed system.** Under the
gate, after ejection, each shed removes `D₀(n)` from `E_int` *and* the same
`D₀(n)` from the remaining ladder cost `Σ(n)`, so the self-bound margin
`G = E_int − Σ(n)` [eV] is exactly shed-invariant (§4c) — evaporation can
never close its own gate. The crossing construction hands the cascade
`E_int = Σ(21)` at gate-open, so `|G|` equals precisely the in-bubble
cooling leak between crossing and ejection, and the energetic floor sits at
the rung where the remaining ladder below costs less than `|G|` (n ≈ 2 at
f_int = 0.50, higher at earlier openings — the P4/I17 leak map). Almost
nothing reaches that floor, because of:

**2. The kinetic wall, and why it double-log-compresses.** Descending, the
RRK argument `x = D₀(n)/E_int` steepens like a hyperbola: `E_int ≈ Σ(n)`
shrinks roughly as n·D₀(1) (the low rungs all sit near D₀(1) below the
Form-U cliff) while `D₀(n)` grows toward D₀(1), so `x(n) ≈ 1/n` near the
bottom. The next shed lands within the remaining flight iff `k·t ≳ 1`,
i.e. `(s−1)·|ln(1−x)| ≲ ln(ν·t_detect) ≈ 17`; that defines a kinetic wall
`x_c(s) = 1 − e^(−17/(s−1))` (≈ 0.44 at s = 30, ≈ 0.91 at s = 8). Both
compressions follow at once: in **s**, x_c moves enormously but the steep
`x(n) ~ 1/n` maps it to Δn ≈ 1–2 (and at s ≤ 8 the energetic floor
intervenes first — hence frozen at 2.7–3.1 vs kinetically parked at 4.1);
in **t**, flight time enters only through ln(ν·t), so a decade buys a few
percent of x_c — the measured 0.4–0.9 He per decade. This is the **Klots
evaporative-ensemble regime**: an ensemble observed at time t sits at the
rung where `k(E, n)·t ≈ 1` regardless of where it started, drifting
logarithmically. A µs cluster cascade landing there is what cluster physics
says should happen — in this narrow sense Wave 7 *validates* the
mechanism's long-time structure. (W7-P2 failed by linearly extrapolating a
logarithmic law: four more decades bought ~2 He, not the descent to the
floor.)

**3. Crack 1 — the perpetual cascade is a construction (→ OQ-F).** The
reason evaporation is *exactly* self-sustaining is that a shed drains only
the binding energy: `dE_int = −D₀(n)` and nothing else (verified in the
delivered `internal_energy_budget.dE_int_shed_eV`; the cold-shed kick is
mechanical bookkeeping, not an E_int drain). A real statistical evaporation
also carries away a translational release per shed, of order the
equipartition share `ε ~ (E_int − D₀(n))/s ≈ 5–20 meV` at
E_int ~ 0.1–0.19 eV and s ~ 8–30 — over a ~17-shed cascade that integrates
to ~0.1–0.3 eV, **the same order as the entire Σ(21) ≈ 0.188 eV budget**.
Not a correction: first-order. With ε included, G strictly decreases per
shed → cascades genuinely self-terminate, terminals converge (the
`time_exhausted` class largely vanishes), freezes land at higher n, the
deep-strip reach shrinks — and OQ-B gets *worse* (stripping becomes
harder). Same epistemic class as the RRK-dof convention Wave 2 resolved: a
statistical-mechanics convention, not a calibration knob, silently
controlling qualitative behavior.

**4. Crack 2 — the suppressed class is physically a fragmentation channel
(→ OQ-B, sharpened into a hypothesis).** A `suppressed` ion rides 8.5 µs
carrying `E_int > Σ(n)` — *more internal energy than the total binding of
its entire shell*. The model parks it at its handover n by deliberate OQ-B
non-resolution; physically, a net-self-unbound complex does not arrive
intact. If the suppressed class fragments, the detector map inverts in
meaning: **bare I⁺ becomes the crossed-after-ejection side of the Δ× race**
(not the tail of the cascade side, where every kinetic lever has failed to
produce it), the experimental **bimodality falls out of the cliff structure
for free** (one ensemble distributed in Δ× splits into
suppressed→fragmented→bare and cascade→broad shell weight n ≈ 3–13), and
the B.1(2) production prediction inverts: at 2.70 eV / f_int = 0.5 the
whole ensemble is suppressed — previously "sheds nothing, dead arm", it
would read "entirely bare-candidate", on the same axis as the 43 % bare
peak.

**One sentence:** the cascade side of the model behaves like a textbook
evaporative ensemble whose detector read is owned by the race margin and
the ladder geometry (s_eff demoted to a clock), and the two features that
look wrong or missing — the eternal cascade and the bare peak — both point
at the same two energy-bookkeeping conventions (no per-shed KE release;
suppressed = inert), not at the swept knobs: the probe program has
exhausted the calibration space and now presses on mechanism conventions,
which is exactly what it was built to localize.

### The top two experimental bins (post-Wave-7 discussion continued, 2026-07-08) — can the never-opened side also explain I⁺He?

The experimental reference's two largest bins are bare **43.5 %** and
I⁺He (n = 1) **17.5 %**, then a smooth decreasing tail (8.0 / 5.2 / 4.0 /
3.3 % for n = 2–5). The delivered ladder (mixture, κ = 1) has a **flat
bottom**: D₀(n) = 9.22 meV for every low rung, Σ(n) ≈ n·9.22 meV.

**Structural fact — the cascade side cannot make n = 1 either.** The
cascade-side terminal rung is a direct readout of the in-bubble leak in
units of one rung (freeze at the smallest n with Σ(n−1) < |G|):

- **bare requires |G| = 0 exactly** (measure zero — *why* every kinetic
  lever failed: any positive leak leaves the cascade one increment short,
  by G-invariance);
- **n = 1 requires |G| < 9.22 meV** (a leak under one rung), *and* that
  sliver is kinetically strangled — ending at n = 1 means shedding at
  n = 2 with E_int barely above D₀, i.e. x → 1 where the RRK bracket
  collapses. Doubly suppressed; hence Wave 7's global minimum n = 2 with
  n = 1 absent in 10 200 ion detections.

The experiment's two largest bins are therefore **both structurally
outside the cascade's reach** — they must come from the never-opened side
(or from ladder physics the flat bottom lacks, below).

**The never-opened side under the three fragmentation specs:**

- **(a) Inert (delivered semantics):** contributes only n = 21. No.
- **(b) Gateless RRK boil-off, K1 = D₀-only kept:** G > 0 is
  shed-invariant, so no energetic floor is ever hit, and at n = 1 the
  delivered rate is the *direct* channel k = ν (no RRK barrier) — the
  complex boils to **exactly bare, every time** (at G₀ ≈ 0.09 eV the
  rates stay ~0.2 ps⁻¹ down the ladder; done in tens of ps). The
  cliff-edge fringe (G₀ ≲ one rung) parks *kinetically* at n = 2 —
  never at n = 1, which has no barrier to hide behind. Explains the bare
  peak; gives **zero** at n = 1.
- **(c) Boil-off with per-shed release ε (OQ-F):** G decreases by ~ε per
  shed; the excess burns out in ~G₀/ε sheds and the complex freezes
  energetically a few rungs later — **including at n = 1**, now generic
  rather than fine-tuned. The ensemble's G₀ distribution is one-sided
  starting at 0 (cliff-edge trajectories), so the never-opened class
  splits into **large-G₀ bulk → bare** plus a **small-G₀ fringe → a
  decreasing small-n tail with n = 1 populated**. Quantitative anchors:
  G₀ ≈ 0.092 eV at 0.80 eV (Wave 4's held E_int = 0.280 eV vs
  Σ(21) = 0.188 eV); G₀ ≈ 0.35 eV estimated at 2.70 eV — against a
  full-ladder ε cost 21·ε ≈ 0.1–0.4 eV, so **bare is marginal at the
  validation budget and robust at production**: a bare peak that
  strengthens with budget is a genuine, testable prediction of this
  picture.

**Answer: the never-opened side explains n = 1 only jointly with OQ-F.**
Neither crack alone produces both top bins — spec (b) gives bare without
n = 1; ε without fragmentation gives neither. OQ-B's endpoint
distribution *is set by* OQ-F's ε: they are one mechanism discussion, not
two.

**The honest alternative for n = 1 (→ OQ-G):** the elevated n = 1 / n = 2
ratio (2.2×, then smoothly 1.5×, 1.3×) could be thermodynamic rather than
race statistics — if the *real* first rung is much deeper than the outer
ones (ion-induced-dipole binding of the last He on I⁺, plausibly tens of
meV vs the flat 9.22 meV Form-U bottom), n = 1 is a natural
"last survivor" of any cascade, a magic-number-like stability the
flat-bottom ladder cannot express. Identifiability consequence: **the
small-n abundance tail is the first observable found that reads the
ladder's *bottom*** (every in-window quantity probes only the top rungs
near n = 21, and the ladder shape is one of the two genuinely-free knobs).
n = 1 is arguably the single most mechanism-discriminating bin in the
distribution: it separates "fragmentation fringe with ε" (tail shape set
by the G₀ distribution; budget-dependent) from "deep first rung" (tail
shape set by D₀(1)/D₀(2); budget-robust) — and the two are combinable.

### The f_int parametrization (post-Wave-7 discussion continued, 2026-07-09) — OQ2 fires

User objection: "if f_int = 0.5, half of the Coulomb explosion goes into
internal energy, then the explosion is only half as violent — this can't be
right; actual values should be ~1 % or lower." Verified against the
delivered code, the objection holds in a **sharper** form:

- **The model implements no partition at all.** `E_int(0) = f_int·E_avail`
  is deposited once at t = 0 (S2 onset, `ion_initial_state.py` →
  `e_int_onset_eV`) and **nothing is subtracted from the mechanics** — the
  Coulomb dynamics and initial velocities are f_int-blind. At f_int = 0.5
  the model books 0.40 eV (validation) / 1.35 eV (production, ≈ 7× the
  total shell binding) of shell heat with no mechanical source. A *true*
  partition at 0.5 would slow the fragments by √2 and would have broken the
  Tier-0/1a drag validation and any VMI comparison — the delivered velocity
  physics is consistent only because the "partition" isn't one. (Row 14's
  soft upper ~0.2 "advisory" already flagged the unease.)
- **But the literal ~1 % kills the mechanism.** f_int = 0.01 →
  E_int(0) = 8 meV at 0.80 eV — barely one rung (D₀(21) = 5.97 meV), far
  below Σ(21) = 188 meV: never suppressed, ~1 shed, frozen at n ≈ 20
  (27 meV → 3–4 sheds at 2.70 eV). No crossing, no race, no deep strip, no
  bare. **There is no physically comfortable literal value — the
  parametrization is broken, not the number.**
- **The row-14 floor identity says the natural variable is absolute.** The
  scenario-keyed floors multiply out to the *same energy*:
  0.235 × 0.80 ≈ 0.188 eV and 0.065 × 2.70 ≈ 0.176 eV — both ≈ Σ(21). The
  "scenario-keyed f_int floor" is the single statement
  `E_int(0) ≥ Σ(21)` refracted through a fraction-of-budget
  parametrization; in absolute `E_int(0)` [eV] the scenario-keying
  dissolves, and the band the whole probe program found effective is
  E_int(0) ∈ [~0.19, ~0.5] eV at both budgets.
- **Physical sources that do not touch the KER** (and hence evade the
  objection): (i) solvation **reorganization on vertical ionization** — the
  shell equilibrated to neutral I₂ suddenly sits displaced in the ionic
  potential; deposit up to ~Σ(21) ≈ 0.19 eV, budget-independent; (ii)
  **electronic / spin–orbit relaxation of nascent I⁺** (~0.7–0.9 eV
  fine-structure scale, ¹D ~1.7 eV) partially degrading into the shell —
  budget-independent, big enough, and it ties the hitherto-free **picture
  knob to the E_int(0) provenance**; (iii) retained drag heating (bounded
  by the drag share of E_dissip — tens of meV, a supporting term); (iv) the
  literal KER coupling, ~1 % = 8–27 meV (the user's estimate — real but
  minor).
- **Working hypothesis (pending literature validation — NOT adjudicated):**
  `E_int(0)` ≈ 0.2–0.5 eV **absolute**, essentially budget-independent,
  sourced from reorganization + electronic relaxation. Consequences if
  confirmed: the B.1(2) production prediction inverts (t× stops scaling
  with budget; what changes at 2.70 eV is the earlier t_eject), the
  bare-peak budget-dependence (§4e above) re-routes through earlier
  ejection (less cooling time → larger G₀ → *more* bare at production —
  same sign, different mechanism), the Δ× race coordinate survives
  untouched, and f_int reclassifies from a dimensionless Bounded fraction
  to an absolute Bounded energy [eV].
- **This is CALIBRATION_MAP OQ2 ("KE_shed / partition") firing.** It
  couples to OQ-F (what a shed drains), OQ-B (what suppression means), and
  OQ-G (the ladder bottom) — all four live in the same energy-bookkeeping
  layer at the ionization/shed/self-unbound boundaries. → The consolidated
  register for the literature-research phase is
  `RESEARCH_QUESTIONS.md` (user decision, 2026-07-09).

---

## 4f. Wave 8 — the suppressed-fraction cliff is a delta; the sub-rung sliver expresses n = 1; the budget knob is bookkeeping-only

Executed 2026-07-09 (under the `[PROCEED TO IMPLEMENTATION]` trigger;
Addendum D of the probe plan; scratchpad drivers through the delivered
generator/detection/report pipelines — zero repo-code change). Step 1 is a
**zero-MD analysis** of the on-disk f_int = 0.65 dirs; step 2 added **6 new
MD runs** (not the planned 7 — the step-1 delta collapsed the three
"ECDF-crossing" points into a two-point bracket), all with `detection.npz`;
full re-score over all **108** dirs (ledger residual uniform ≈ 2.23·10⁻⁵ eV;
bridge + gated fi0.50/fi0.65 wiring oracles exact).

### Step 1 — the per-ion critical value is one number, not a distribution

Pre-gate-open the ride is shed-free and pickup-free (n ≡ 21 at Langmuir
saturation), so the gated K2 drain is source-free and linear:
`E*_i = Σ(21)·E₀/E_int,i(∞)` read directly from the stored histories.
Result: **E*_i = 0.461218816 eV with ensemble width 4.2·10⁻¹³ eV** — a
delta function, not an ECDF. `f_int* = 0.5765` at 0.80 eV. Root cause,
verified in the artifacts: **the probe ensemble is 100 kinematically
congruent replicas** — every molecule sits at its droplet center
(|r|₀ = 4.5 Å), onset speeds agree to 10⁻¹⁹ Å/ps, one droplet radius
(27.936 Å), so the radial trajectories (2·10⁻¹¹ Å agreement at 5 ps) and
hence the per-ion cooling-exposure integrals are identical:
K_tot = ln(E₀/E_end) = **0.898297** for every ion. Analysis validation:

- from-positions replication of the exposure (per-step
  `exp(−dt·ρ_erfc(r−R)/τ)` through the shared gate surface, depth = r − R):
  max |E_pred − E_stored| = 1.9·10⁻¹⁵ eV;
- pre-open E_int linearity across f_int (f0.50 vs f0.65 dirs): residual
  9.4·10⁻¹⁶ eV over 63,300 ion-steps; pre-open positions **bit-identical**;
- opening-time law `K(t_open) = ln(E₀/Σ(21))`: predicts the observed
  f0.50 opening at 6.33 ps exactly.

### Step 2 — six bracket runs, every pre-registered number hit

All gated `density_scaled`, τ = 6.55, probe pins; E₀ = f_int·E_avail [eV]:

| run | E₀ [eV] | predicted | measured |
|---|---|---|---|
| b080 f_int 0.57, s30 | 0.4560 | opens 10.83 ps; 0 % supp; \|G\| ≤ 2.12 meV | 100/100 open **at 10.83 ps**; 0 % supp; G ∈ [−2.12, −1.79] meV; n_relax 3.87; **n_detect 2.00** (100 % texh) |
| b080 f_int 0.57, s8 | 0.4560 | same opening (s-blind); **energetic floor n = 1** | opens 10.83 ps identically; all 100 frozen at **n = 1** by 633 ps; **n_detect 1.00** (100 % frozen) |
| b080 f_int 0.58, s30 | 0.4640 | never opens (+1.13 meV) | 0/100 open; 100 % suppressed at n = 21; K_tot identical |
| b270 anchor, s30 | 0.5200 | suppressed, shed-free (exposure source) | 100 % suppressed; **K_tot = 0.898297 ≡ the 0.80 eV value** |
| b270 f_int 0.17, s30 | 0.4590 | opens 11.96 ps; \|G\| ≤ 0.90 meV | opens **at 11.96 ps**; G ∈ [−0.90, −0.73] meV; n_relax 3.78; n_detect 1.98 (2 ions **frozen at n = 1**, 98 texh) |
| b270 f_int 0.18, s30 | 0.4860 | never opens (+10.1 meV) | 0/100 open; 100 % suppressed |

### Prediction verdicts

- **W8-P1 (MD lands on the step-1 prediction) — CONFIRMED exactly:**
  fractions, opening times (10.83 / 11.96 ps), and |G| windows all match;
  the pre-open separability assumption is airtight (bit-identical pre-open
  trajectories across f_int).
- **W8-P2 (finite cliff width) — REFUTED.** Width 4·10⁻¹³ eV: machine
  noise, not physics. There is no transition region.
- **W8-P3 (in-band 43.5 % crossing) — REFUTED.** The suppressed fraction is
  a step 0 → 1 at E* = 0.4612 eV: **no E_int(0) splits the ensemble**; the
  only available fractions are 0 % and 100 %. (E* itself is in-band —
  inside the sourced 0.2–0.5 eV window.)
- **W8-P4 (opened side ≠ the broad tail) — CONFIRMED, sharpened.** Near the
  cliff the opened side is a near-delta at n = 1–2 (not n ≈ 3–5 as
  predicted — the sub-rung leak was not anticipated), even farther from the
  broad experimental tail than predicted.
- **W8-P5 (bare grows with budget via earlier t_eject) — REFUTED IN-MODEL,
  at the mechanism level.** `coulomb_available_eV` has exactly **one**
  physics reader: the S2 onset deposit `E_int(0) = f_int·E_avail`
  (`ion_initial_state.py`). The Coulomb mechanics, trajectories, ejection,
  and exposure are **budget-blind** — K₂₇₀ ≡ K₀₈₀ to 10⁻¹³, and the MD
  bracket confirms the same absolute step (E* between 0.459 and 0.486 eV at
  2.70 eV). In absolute E_int(0) coordinates the delivered model is exactly
  budget-invariant; "production" currently changes the E_int bookkeeping
  and nothing else.

### Findings

1. **OQ-B's "high enough f_int" question is answered: yes, but all-or-none.**
   Above E* = 0.4612 eV the *entire* ensemble arrives `suppressed`
   (100 % bare-candidate under the RQ3 fragmentation hypothesis) — never
   43.5 %. A single-knob-point ensemble cannot express the experimental
   bare/shell coexistence, not approximately, not at any value: the Δ×
   race margin is ensemble-uniform by construction. I12/I18's "the shape
   must come from ensemble heterogeneity" upgrades from an inference to a
   structural theorem about the current preset: **heterogeneity must be
   injected** (droplet-radius distribution, thermal onset spread, an
   E_int(0) distribution) before the bare fraction can take any value
   between 0 and 1.
2. **The sub-rung sliver expresses n = 1 — the previously unreachable
   second-largest experimental bin.** For E₀ ∈ (E* − D₀(1), E*) the leak
   |G| is under one rung, so the opened cascade's energetic floor is
   **n = 1** (not the n ≈ 2–3 of every earlier wave, whose leaks were
   multi-rung). MD: all 100 ions frozen at n = 1 (s_eff 8; reached by
   633 ps, detector-stable), and even at s_eff 30 the deeper-sub-rung 2.70
   point put 2/100 at n = 1 within flight. The I24 "n = 1 is doubly
   suppressed" statement is now refined: the kinetic strangulation applies
   only as |G| → one full rung; a *shallow* sub-rung leak is kinetically
   open. Consequence for the tail: the flat-bottom ladder **can** populate
   n = 1 without OQ-G's deep first rung — but only from a window hugging
   the cliff from below that is one rung wide in *leak* units, i.e.
   **D₀(1)·e^(K_tot) = 22.6 meV wide in E_int(0)** (E₀ ∈ 0.4386–0.4612 eV;
   ~7.5 % of the sourced band). Feeding the 17.5 % experimental n = 1
   weight from this window constrains any heterogeneity model; OQ-G depth
   remains the competing (combinable) route. *(Correction 2026-07-09: an
   earlier revision stated the window as 9.22 meV in E₀ — that width is
   correct in |G| units only; the e^K = 2.455 map factor was dropped.)*
3. **The budget knob is bookkeeping-only — the production scenario's
   kinematics are unmodeled (new OQ-H).** The 2.70 eV channel physically
   comes with different Coulomb-explosion mechanics (higher KER = faster
   fragments, earlier ejection, shorter exposure); the delivered model
   implements none of that. Two standing statements void in-model: the
   B.1(2) production prediction in *both* its versions (original t×
   scaling and the RQ1-inverted earlier-t_eject route), and the
   bare-vs-budget discriminator separating RQ3 from RQ4 — currently
   untestable in-model on the kinematic route. What survives: in absolute
   E_int(0) the whole Wave-4..8 gated map transfers to production
   verbatim (one number, E*, for every budget).
4. **The delta cliff is preset-conditional, not mechanism-intrinsic.** The
   congruence traces to the deterministic single-droplet center-placed
   onset of the probe preset. The heterogeneity axes exist in the codebase
   (droplet-size distribution preset; thermal sampling) — the cliff
   *anatomy* under a heterogeneous ensemble (each ion carrying its own
   E*_i via R_i-dependent exposure) is exactly what a future
   heterogeneity wave / the campaign would map. The step-1 method
   (per-ion E*_i from stored histories) applies unchanged there.

### The heterogeneity program (post-Wave-8 discussion, 2026-07-09) — the histogram as a margin-distribution readout, and the two injection axes

Wave 8's negative result (I26: no single-point split) turned constructive
in discussion: because the per-ion terminal bin is a deterministic
function of the leak |G| = Σ(21) − E₀·e^(−K), **the experimental
histogram is a linear readout of the ensemble's margin distribution**
through fixed windows — at the 9 Å pins, 22.6 meV of E₀ per rung
(D₀(1)·e^K, flat-bottom regime), bare above E* = 0.4612 eV, n = 1 at
0.4386–0.4612 eV, one rung per window below. The Wave-6/7 points already
validate the mid-map (E₀ = 0.40 → n̄_detect 3.05; 0.28 → 7.37; 0.192 →
10.6, the deviation being the high-n rung shrink). The implied p(E₀) from
the experimental bins is smooth and unimodal: ≈ 43.5 % above E*, density
≈ 0.77 %/meV at the cliff declining 2.2× / 1.5× / 1.3× per window — mode
≈ 0.45 eV, scale ~50 meV — **inside the RQ1-sourced 0.2–0.5 eV band
without tuning**. That convergence (the sourced provenance range and the
histogram-required range agree) is the strongest circumstantial support
the biphasic mechanism has received. Wave 9 (Addendum E; approved
2026-07-09) inverts this properly, with the fitted p(E₀) treated as a
**physics claim** on E_int(0) provenance (user decision (b)): a
multi-modal fit would read as electronic fine-structure branching (RQ1).

**The two injection axes, and why the droplet-size distribution is the
physically-grounded future one.** Margin heterogeneity can enter through
E₀ (the reservoir; Wave 9's axis — cheap, zero new physics via run
mixtures) or through the **exposure K** — and K is set by the droplet
radius (larger droplet → longer in-bubble transit → larger K → larger
per-ion E*_i and deeper leak at fixed E₀). The experimental droplet
ensemble is a **log-normal size distribution** — this axis is not a
modeling option but a known physical fact the current single-radius
preset suppresses; the codebase already carries
`single_pulse_droplet_distribution` on the neutral side. The two axes
own different parts of the histogram:

- the **near-cliff structure** (bare / n = 1 / n = 2 ratios) is
  E₀-resolution-limited — the E₀ axis owns it;
- the **deep-shell tail** (n ≳ 11) requires E₀ below the 0.22 eV
  solvation floor under pure-E₀ heterogeneity (the W9-P2
  pre-registration) — the droplet-R tail supplies exactly this weight
  with in-band E₀ (bigger K, deeper leak);
- the **near-absence of n ≥ 19** (< 1 % experimentally) *constrains* the
  large-R tail: a large-enough droplet fully quenches the cascade inside
  (the ungated-like shell-retaining corner), so the experimental gap at
  high n bounds how much large-R weight the ensemble can carry.

A future droplet-distribution wave/slice is therefore the expected
follow-up once Wave 9 quantifies the below-floor demand — it is a real
build (wiring the distribution preset into the Tier-2 ion pipeline;
onset-shell size vs R is an open design question), unlike Wave 9's pure
post-processing. Joint (E₀, R) fitting is the natural endpoint; the
step-1 per-ion E*_i extraction method applies unchanged to a
heterogeneous ensemble.

**s_eff identifiability — possible, not fixed (user decision (c)).** The
n = 1 window populates its bin only if the cascade reaches the floor
(s_eff = 8: 100 %; s_eff = 30: kinetic wall parks it at n = 2 — the
measured fi0.57 pair). If the mixture picture holds, the small-n bins
would make s_eff detector-identifiable for the first time (Wave 7 found
the detector nearly s-blind). This is **recorded as a possibility only**
— it awaits the W9-P4 read and adjudication; no campaign re-scope is
implied by it.

### On the run count and the b270 staircase columns

The addendum's 7-run matrix executed as 6 (the delta has one crossing; a
2-point bracket verifies it — the third transition point had nothing left
to measure). The b270 rows carry staircase-metric columns scored against
the 9 Å / 0.80 eV anchored comparator; they are **N/A physically** (the
generator-level budget guard exists precisely to prevent this pairing in
campaign use; the scratchpad route bypassed it deliberately, for exposure
measurement, not staircase comparison) — read only their suppression /
terminal / detected columns.

---

## 4g. Wave 9 — the E₀-mixture inversion: the experimental histogram is a smooth in-band p(E₀); the top bins reproduce to W₁ = 0.086 with in-band + above-cliff mass alone

Executed 2026-07-09 (under the `[PROCEED TO IMPLEMENTATION]` trigger;
Addendum E of the probe plan; scratchpad drivers through the delivered
generator/detection pipelines + a pure post-processing fit — zero repo-code
change). **8 new MD runs** (7 primary `s_eff = 8` columns + 1 `s_eff = 30`
companion), each with `detection.npz` at the Sourced t_detect = 8.53 µs →
**116 probe dirs**. Wiring oracles exact: ledger residual uniform
2.231·10⁻⁵ eV across all 8 new dirs; the on-disk fi0.50 s8 / s30 and fi0.65
points re-score unchanged (n_detect 3.05 / 4.09 / 21.00 — the Wave-7 values).

**User decisions carried into execution (post-design, 2026-07-09):** (1) the
fit is a **simplex-constrained** least squares (non-negative *and* Σw = 1) —
p(E₀) is a probability density, not raw NNLS; (2) the suppressed class
collapses to a single **bare** column under RQ3 spec (b), implemented per-ion
(`state_reason == "suppressed"` → n = 0); (3) the deep-tail shortfall is
**reported, not forced** — the basis is not extended to cover bins it cannot
reach; (4) W9-P4 is a **localized column-swap demonstration**, no campaign
re-scope.

### The basis lands on the window map (E₀ = f_int · 0.80 eV)

The 7 new `s_eff = 8` columns arrive exactly where the flat-bottom window map
(§4f) places them, all **frozen/converged** (the floors are *reached* at
s = 8, Wave-8 I27):

| f_int | E₀ [eV] | n̄_detect | arrival | window |
|---|---|---|---|---|
| 0.53 | 0.424 | 2.04 | 100 % frozen | n = 2 |
| 0.48 | 0.384 | 3.38 | 88 % frozen | n = 3–4 |
| 0.45 | 0.360 | 4.31 | 91 % frozen | n = 4–6 |
| 0.42 | 0.336 | 5.22 | 96 % frozen | n ≈ 5–6 |
| 0.20 | 0.160 | 12.26 | 91 % frozen | n ≈ 12 |
| 0.15 | 0.120 | 14.48 | 87 % frozen | n ≈ 14 |
| 0.10 | 0.080 | 16.62 | 89 % frozen | n ≈ 16–17 |

The `s_eff = 30` companion at f_int = 0.53 arrives **n̄_detect = 3.02, 100 %
`time_exhausted`** (vs the s = 8 sibling's n = 2, frozen) — the kinetic wall
measured a third time (with the on-disk fi0.57 s8/s30 pair n = 1 vs 2), the
W9-P4 lever.

### Primary fit (13 `s_eff = 8` columns) — W9-P1 CONFIRMED

The experimental histogram is reproduced to **L2 residual 0.016,
Wasserstein-1 = 0.086 bins**. The fit is near-diagonal (each near-cliff
column owns one bin; bare is fed **solely** by the suppressed→bare column,
weight 0.436 ≈ experimental 0.435):

| bin | exp | fit | resid |
|---|---|---|---|
| bare | 0.435 | 0.436 | +0.001 |
| n = 1 | 0.175 | 0.176 | +0.001 |
| n = 2 | 0.080 | 0.081 | +0.001 |
| n = 3 | 0.052 | 0.054 | +0.002 |
| n = 4 | 0.040 | 0.040 | −0.000 |
| n = 5 | 0.033 | 0.038 | +0.005 |
| n = 6 | 0.028 | 0.016 | −0.012 |
| n = 7–17 | (tail) | (tail) | \|resid\| ≤ 0.003 |
| n = 18–20 | 0.010 tot | 0.001 | −0.009 tot (basis floor) |

The top six bins match to < 0.006. Near-cliff **density on the physical
22.6 meV rung window** (D₀(1)·e^K, the convention W9-P1's prediction lives
in): n = 1 → **0.780 %/meV** (predicted 0.77), n = 2 → 0.375 %/meV,
ratio **2.08×** (predicted 2.2×). **The top bins are reproducible with in-band
+ above-cliff E₀ mass alone** — W9-P1 confirmed quantitatively.

### Implied p(E₀), W9-P2, and the provenance claim (W9-P3)

The fitted weights, read as a density over E₀, give:

- **mass above E\* (0.4612 eV, the bare/suppressed class): 0.436** — a point
  mass (all E₀ > E* collapse to bare under spec (b)), matching the
  experimental bare bin;
- **mass in the sourced band [0.22, E\*]: 0.467** — the resolved shell bins;
- **mass below the 0.22 eV solvation floor: 0.096** (**W9-P2 CONFIRMED**):
  9.6 % of the ensemble must sit below the RQ1 floor (columns E₀ = 0.192 /
  0.16 / 0.12 / 0.08 eV) to feed the deep tail — *and* n = 18–20 (~1 %
  experimentally) lies **beyond the basis floor** (needs E₀ < 0.08 eV), an
  additional below-floor demand. This quantifies, in advance, how much weight
  the droplet-radius axis (§4f) must supply with in-band E₀ (larger K → deeper
  leak at fixed E₀).

**W9-P3 (the physics claim) — no multi-modal signature; mode at the cliff,
not an interior bump.** The implied p(E₀) is **unimodal** (no genuine second
mode), but the mass does **not** sit in a symmetric ~50 meV bump at 0.45 eV as
pre-registered: it **rises monotonically toward the cliff** and concentrates
in the ≥ E* class (0.436), with a declining in-band tail through ≈ 0.28 eV and
the below-floor tail (~0.10). So the provenance read is **E_int(0) concentrated
at/above ≈ 0.46 eV with an in-band declining tail** — inside/around the RQ1
[0.2, 0.5] band but pressed to its **upper edge and above E\***, consistent
with reorganization + electronic-relaxation sourcing near the top of the band
(RQ1/OQ2). The interior density carries roughness (a weight spike at
E₀ = 0.384 / n = 4, zeros at the collinear 0.36 / 0.40 columns) that is a
**basis-collinearity artifact** (fi0.45 / 0.48 / 0.50 all peak at n ≈ 3–4),
**not** physical multi-modality — the fit does **not** demand electronic
fine-structure branching. (A branching signature would have to survive column
de-collinearization; this one does not.)

### W9-P4 (s_eff identifiability) — CONFIRMED, recorded as possibility only

Re-fitting with the two near-cliff columns swapped to `s_eff = 30`
(n_detect 1 → 2 and 2 → 3, the kinetic wall) leaves the **n = 1 bin unfed**
(no s = 30 column reaches n = 1): fit n = 1 = 0.000, resid **−0.175**; L2
jumps 0.016 → 0.182 and Wasserstein 0.086 → **1.398 (16×)**. The detector's
small-n bins **select the kinetics band** — s = 8 (floors reached) fits, s = 30
(kinetic wall) cannot, so **s_eff is detector-identifiable via the n = 1 bin**
for the first time (Wave 7 found the detector near s-blind at the deep-strip
corner). Per user decision (c): **recorded as a possibility, no campaign
re-scope**.

### What Wave 9 establishes

The Wave-8 negative result (no single knob point splits the ensemble; the
cliff is a delta, I26) is now turned fully constructive: because ions never
interact, a weighted mixture of single-point runs **is** an ensemble with
p(E₀) on the grid, and the experimental histogram **inverts cleanly** into a
smooth, unimodal, in-band-plus-above-cliff p(E₀) (top bins to W₁ = 0.086)
**without tuning** — the sourced provenance range and the histogram-required
range coincide (the strongest circumstantial support the biphasic mechanism
has received, now quantitative). The two standing caveats survive intact and
are quantified: the deep tail demands ≥ 9.6 % below-floor mass (the
droplet-R axis's job, W9-P2), and the whole inversion is **conditional on
RQ3 spec (b)** (suppressed → cleanly bare; an ε-type channel would add
suppressed-side small-n weight and re-shape the near-cliff fit — OQ-B/OQ-F).
Wave 9 builds **no p(E₀) sampling surface** and does **not** discharge the F5
gate; the mixture lives entirely in post-processing.

### Physical-justifiability of the landed parameters + interpretation risks (post-execution discussion, 2026-07-09)

A post-execution discussion assessed whether the Wave-9 operating point sits
in a physically defensible range. Recorded as caveats on the *interpretation*,
not new results.

**s_eff = 8 — order-of-magnitude defensible, but a rate knob and a
methodological pick.** Physically it says the He₂₁ shell acts *kinetically*
like ~8 active oscillators, not the classical `3n−3 = 60`; the ~10× reduction
is the expected direction for cold quantum He (quantum mode-freezing — only
soft delocalized intershell/surface modes participate on the sub-ps shed
timescale; cf. Klots evaporative-ensemble treatments of weakly-bound/quantum
clusters; an effective heat capacity ~10–20 % of Dulong–Petit `3Nk_B` is
plausible for He — consistent, not precise). Three risks on reading "8" as a
physical constant: (i) s_eff is fundamentally a **rate** knob (I13) — under the
gate it does not set the terminal n, and the detector read compresses it (I20),
so the exact value does less work than it looks; (ii) **s = 8 was chosen for
the Wave-9 basis because at s = 8 the columns reach their *energetic* floor**
(E₀ → n deterministic, the window map) while s = 30 stalls at the kinetic wall
— a methodological choice, not a fit; the only physical signal on the value is
W9-P4 (the small-n bins prefer the fast-kinetics side); (iii) constant s is a
diagnostic simplification — the physical form scales with size
(`α·(3n−3)`, OQ-C), so "8" is best read as the effective bath near the *top* of
the ladder (n ≈ 21, where the in-window sheds happen).

**E_int(0) — band defensible; the upper-edge concentration is a hard physical
demand.** The fit concentrates **43.6 % of the mass above E\* = Σ(21)·e^K =
0.461 eV**, i.e. E_int(0) ≳ 2.4× the *total* shell binding for the bare class.
Solvation reorganization on vertical ionization caps at ~Σ(21) ≈ 0.19 eV
(budget-independent, the Na⁺-benchmark route), so the > E* mass **requires**
an additional ≳ 0.27 eV deposit from I⁺ **electronic / spin–orbit relaxation**
(³P fine-structure scale ~0.94 eV) degrading into the shell. Wave 9 therefore
effectively claims **the bare peak is electronically sourced**, tying the
hitherto-free picture knob to E_int(0) provenance (RQ1) — falsifiable. Risk:
the ≥ 9.6 % below-0.22-eV-floor mass (I30) is **not** E₀-justifiable — it must
come from the **droplet-radius axis** (larger K), not the reservoir; pure-E₀
heterogeneity cannot supply the deep tail.

**τ = 6.55 ps — plausible, factor-~2.5 uncertain** (geometric-mean pin; GAH25
Na⁺ *growing-shell* provenance, a different scenario than an ejected shedding
complex). **κ = 1 / ladder bottom — a config default, not a data-justified
value** (κ is near-dead for the staircase — inverted + normalization-capped;
the flat-bottom ladder `D₀ ≈ 9.22 meV` is physically suspect at the bottom, the
real last-He-on-I⁺ rung plausibly much deeper — OQ-G; one of the two
genuinely-free knobs, read by the n = 1 bin). **gate = density_scaled — the
most physically-grounded choice** (bath-gated cooling; shared ρ_He boundary
with drag/pickup).

**Interpretation risks carried from execution.** The whole inversion is
conditional on **RQ3 spec (b)** (suppressed → cleanly bare; an ε channel /
OQ-F would add suppressed-side small-n weight and reshape the near-cliff fit);
N = 50 single seed at 1 % granularity; the interior p(E₀) carries
**basis-collinearity roughness** (the fitted *weights* are robust, the
*density* readout is window-convention-sensitive); no p(E₀) sampling surface is
built and nothing here discharges the F5 gate.

---

## 4h. Wave 10 (Steps 1 + 1.5) — the width lives on the E₀/K axes; narrow E₀ fails at the probe K-scale and lands the bare peak at the production K-scale bracket

Executed 2026-07-10 (under the `[PROCEED TO IMPLEMENTATION]` trigger;
scratchpad route `wave10_step1_step15.py` — **zero new MD, zero repo-code
change, no run artifacts touched**): the two zero-MD steps of Addendum F *as
amended* (K-reframing F.2b; Step 1.5 F.3b; corrected W10-P3; K-scale leg
W10-P5). Step 2 (the droplet slice) and F.5 (the n-graded ladder arm) remain
un-built and gated. Inputs: all 41 gated (`_cgds`) probe dirs re-read
(detection.npz per column; stored trajectories of the fi = 0.65 suppressed
dir for exposure integrals), the F.2b closed form on the delivered ladder
(Σ(21) = 0.18783720 eV), and `integrated_i_he_abundance.csv`. Scoreboards:
`wave10_step1_columns.csv`, `wave10_step15_K_of_R.csv`,
`wave10_step15_scan.csv` (scratchpad).

### Step 1 — width decomposition at fixed (E₀, droplet) — W10-P1 CONFIRMED

- **Wiring oracle.** The per-ion exposure integral
  `K_i = (1/τ)·∫ρ̂(|r_i(t)|−R₀)dt` over the stored ion + relaxation
  trajectories reproduces the Wave-8 anchor: mean K = 0.897541 vs 0.898297
  (0.08 %, trapezoid-vs-per-step discretization); per-ion spread
  (max−min) = 9.0·10⁻¹³. The implied per-ion `E*_i = Σ(21)·e^{K_i}` has
  spread **4.15·10⁻¹³ eV — the Wave-8 delta (I26) re-derived from
  trajectories** rather than from the bisection.
- **Opened side.** Across all 41 gated columns the detected read at fixed
  (E₀, droplet) is **≤ 1 bin wide**: n_detect MAD 0.04–0.88 He (every
  column < 1), min–max spans 1–5 bins, modal-bin fraction 0.36–1.00.
  Near the cliff the columns are essentially delta (fi ≥ 0.48: MAD ≤ 0.47;
  the fi = 0.57 column is n = 1 for all 100 ions).
- **Verdict (W10-P1).** Poisson pickup + trajectory heterogeneity carry
  ≈ ≤ 1 bin of the ~21-bin experimental span. At fixed droplet the entire
  Wave-9 tail width came from the E₀ mixture axis; **under a narrow E₀ the
  droplet-K axis is the only remaining width carrier** — the Step-2 target
  is now quantified as the full span.
- **Bonus re-confirmation (I28).** The b270 columns line up in absolute E₀:
  fi = 0.17 (E₀ = 0.459 eV < E\*) opens and lands n̄ = 1.98; fi = 0.18/0.19
  (0.486/0.513 eV > E\*) arrive 100 % suppressed — budget-invariance holds
  in the detected read.

### Step 1.5a — the closed form is MD-validated; K(R) is steep

The F.2b fate map (bare ⇔ K < K\* = ln(E₀/Σ(21)); else no-shed leak +
exact ε = 0 ladder descent from n = 21) scored against the 13 gated
`s_eff = 8` columns at the pinned K: **mean |diff| = 0.84 He over the opened
columns, +0.05/+0.04/±0.00 at fi = 0.50/0.53/0.57** — i.e. ≈ exact in the
cliff-adjacent region that owns the bare and n = 1 bins; worst −1.7 He
mid-shell (the no-shed leak approximation slightly over-books the in-bubble
drain when the gate opens early). The pinned-trajectory exposure curve
K(R = R₀·(N/2000)^{1/3}) has slope d ln K/d ln R ≈ 1.23:

| N [He] | R [Å] | K | E\* = Σ(21)·e^K [eV] |
|---|---|---|---|
| 250 | 14.0 | 0.385 | 0.276 |
| 500 | 17.5 | 0.498 | 0.309 |
| 1000 | 22.1 | 0.664 | 0.365 |
| **2000** | **28.0** | **0.899** | **0.462** |
| 4000 | 35.4 | 1.229 | 0.642 |
| 8000 | 44.2 | 1.657 | 0.985 |
| 16000 | 55.9 | 2.231 | 1.748 |

### Step 1.5b — the fate-map scan — W10-P3 CONFIRMED (probe K), W10-P5 CONFIRMED (the split verdict)

Sharp E₀ ∈ [0.24, 0.32] eV × stated log-normal droplet priors
(Kornilov δ = 0.625 about ⟨N⟩ = 2000; δ = 0.40 / 0.80 sensitivity;
pickup-weighted ∝ N^{2/3} variant) × two K-scales (probe 0.898;
production bracket ×0.545 ≈ 0.49 from fragment-speed scaling √(2.70/0.80)):

- **Probe K-scale: the droplet axis is excluded for the bare peak at narrow
  E₀.** At E₀ = 0.28 eV every stated prior yields ≤ 2 % bare (max anywhere:
  12.7 % at the E₀ = 0.32/δ = 0.80 corner); analytically, 43.5 % bare needs
  σ_lnK ≈ 4.95 → **log-normal δ ≈ 12, ~19× the Kornilov width** — not a
  tuning miss but an order-of-magnitude exclusion. W₁ ≥ 2.7 bins everywhere
  at this scale: at narrow E₀ the pinned kinematics put *all* mass
  mid-shell, so the tail *shape* fails too (the W10-P2 analytic pre-verdict
  is negative at the probe scale).
- **Production K-scale bracket: the bare peak lands inside the RQ1 band.**
  bare crosses 43.5 % within the scanned E₀ band for every prior; best
  sampled point (E₀ = 0.28 eV, δ = 0.80): **bare 42.6 %, W₁ = 0.993 bins,
  untuned** (δ = 0.625 gives 33.1 % at 0.28 eV / 55.9 % at 0.30 eV — the
  cliff sits mid-ensemble, K\*(0.28) ≈ 0.40 vs K₀ ≈ 0.49). The small-n
  envelope is monotone as observed; the mid-tail (n = 2–4) overshoots
  (11.6/10.1/8.0 % vs 8.0/5.2/4.0 %); the n = 1 bin under-fills
  (13.6 % vs 17.5 % — model n1/n2 step 1.18× vs experimental 2.18×).
- **W10-P4 pre-read.** The residual n = 1 deficit factor is **1.29× ≤ the
  ~1.4× F.5 electronic-taper cap** — coverable by the capped X₂-first
  graded picture without fitting the ladder to the bins (existence-level;
  the n = 2–4 overshoot is not addressed by the taper and stays open).
- **W10-P5 CONFIRMED — the split verdict, exactly as pre-registered.**
  Narrow E₀ fails at the probe K-scale and succeeds at the production
  K-scale bracket. Consequence: the resolution is **"narrow E₀ + RQ7
  kinematics"** — the amended F.8 outcome-(a)-via-RQ7 exit, *not* the RQ1
  super-solvation route; **the RQ7 kinematics arm is now co-requisite** for
  any Step-2 adjudication, and the droplet slice's job shrinks to the
  width/shape of the K-distribution, not the bare-peak location.

### Boundaries of the Wave-10 zero-MD read

1. The production K-scale is a **ballistic speed-scaling bracket**
   (×√(0.80/2.70) applied multiplicatively to K(R)), not model output —
   RQ7 stays unmodeled; no 2.70 eV MD was run.
2. K(R) is the **pinned-trajectory approximation** (stored r(t) under the
   R₀ drag environment re-integrated against shifted surfaces); the slope
   bias is direction-mixed and small vs the 19× probe-scale margin.
3. **n_eject = 21 for all R** (the delivered model dresses ions fully at
   t = 0); physically small droplets dress less → smaller Σ(n_eject) →
   suppression *easier* — a bias in *favour* of the narrow-E₀ reading.
4. The closed form carries the −1.7 He mid-shell no-shed bias and the
   ≤ 1.4 He detected-read s-dependence (I20) — both ≈ 0 at the cliff.
5. All of this is **analytic pre-verdict**: Step 2 (the MD droplet slice)
   arbitrates; conditional on RQ3 sequential-shed (suppressed → bare) and
   ε ≈ 0. Reported, not auto-adjudicated; nothing here discharges the F5
   gate.

---

## 4i. Wave 11 — K₂.₇₀ measured at 0.746: the ballistic bracket was wrong, the pre-registered E₀ band fails, and the narrow-E₀ landing re-calibrates to E₀ ≈ 0.38–0.41 eV (still in-band) with a 42 %-beyond-band validity caveat

Executed 2026-07-11 (under the `[PROCEED TO IMPLEMENTATION]` trigger given
same day; Addendum G with the opened companion promoted to default).
Scratchpad route `wave11_gen_runs.py` / `wave11_measure_K.py` /
`wave11_step3_scan.py`; **two new MD dirs** (the first off 9 Å geometry:
`..._tier2probe_b270_..._cgds_R2.67`, suppressed primary E₀ = 0.52 eV and
opened companion E₀ = 0.28 eV at `R0_GS_angstrom = 2.666`), zero repo-code
change. Scoreboards `wave11_K_scoreboard.csv`, `wave11_step3_scan.csv`
(scratchpad).

### Step 0 — opening checks (all pass, one correct guard refusal)

- Config-load **accepts** `R0_GS_angstrom = 2.666` (no geometry bound;
  the standing §6.5 pairing warning is the only trip). Both tags free.
- The pure-cubic loader **correctly refuses b = 0** (`linear_cubic with
  a == 0 requires b > 0`) — the drag-off ballistic oracle therefore ran at
  b = 10⁻¹² (γ·t/m ~ 10⁻¹¹, drag-free to machine precision).
- **Ballistic dt oracle:** at the ~11× steeper Coulomb onset, the fixed
  dt = 0.01 ps integrator reproduces the analytic
  v(r) = √(2·[14.4/2.666 − 14.4/r_sep]/2 − E_bind)/m_eff) — including the
  bundle's effective binding well 0.116758 eV, which the fragment climbs
  out of on exit — to **rel. error 8.8·10⁻⁵**. (First oracle attempt
  omitted E_bind and "failed" at 2.2 % ≡ exactly the well depth — a
  correct-physics reminder, not an integrator problem.)

### Step 1 — the measurement

- **W11-P1 (wiring oracle) — CONFIRMED exactly.** The re-implemented
  exposure script reproduces the Wave-10 trajectory-route control at
  **K = 0.897541** (0.8975410665…, per-ion spread 9.0·10⁻¹³ — the I26
  congruence delta re-derived).
- **K₂.₇₀ = 0.74603** (suppressed primary; per-ion spread 7.4·10⁻¹³;
  exposure ∫ρ̂dt = 4.886 ps; t_exit = 4.62 ps vs 5.63 ps at 9 Å). The
  independent E_int cross-read agrees: ln(0.52/0.24642) = 0.7468 (0.1 %,
  trapezoid-vs-per-step discretization). Suppressed ride confirmed
  (n ≡ 21 through ion end, relaxation, and detection).
- **W11-P2 — REFUTED as registered.** K₂.₇₀ = 0.746 sits **above** the
  pre-registered [0.5, 0.7]: the scale factor is S_K = 0.831, nowhere
  near the ballistic ×0.545. Drag eats almost all of the extra channel
  speed: the fragment peaks at 10.51 Å/ps mid-flight and exits at
  3.82 Å/ps — slower transit than ballistic scaling assumed, hence the
  larger exposure. E\*(K₂.₇₀) = Σ(21)·e^0.746 = **0.396 eV**.
- **W11-P3 (congruence + m(t) feedback) — CONFIRMED.** The suppressed
  production ensemble is kinematically congruent (10⁻¹³); the opened
  companion (cascade lightens the ride 21→~7 He in-window) measures
  K = 0.7299 ± 0.0111 — the m(t)-feedback on K is **−2.2 %**, and it is
  the first non-degenerate per-ion K spread ever measured in the probe
  program (0.0222 — channel-RNG-driven, still ≪ 1 rung in E\* terms).

### Step 2 — validity decompositions (the load-bearing caveat) + VMI

- **W11-P4 — SPLIT.** Overlap share (pair separation < 2×4.67 Å):
  **8.0 %** — inside the pre-registered ≲ 10–15 %. Beyond-band share
  (exposure accrued at speed above the 9 Å control's in-window max
  5.23 Å/ps): **41.7 %** (42.7 % on the companion) — **far above** the
  pre-registration. The production fragment spends the first half of its
  exposure at up to **2.0×** the speed the pure-cubic law was calibrated
  on. Both identified biases push K up (cubic force ∝ v³ is steeper than
  form drag; the booked inter-ion density is an over-count), so
  **K₂.₇₀ = 0.746 is an upper-bound-flavored measurement**; the hard
  lower bracket remains the ballistic 0.49.
- **Fragment speeds vs VMI (the independent discriminator).** Detected
  asymptotic speeds: suppressed→bare class 4.11 Å/ps, companion (n ≈ 6)
  5.31 Å/ps — vs the experimental `vmi_iplus_he.csv` peak at
  **10.1 Å/ps (1010 m/s)**. Ballistic references: 20.3 (bare), 15.7
  (dressed n = 21). The experiment sits *between* the model and
  ballistic: the locked drag **over-dissipates production fragments by
  ~2–2.5× in speed** — the first VMI-side kinematics constraint, and its
  direction is *consistent with the over-drag reading* of the
  beyond-band caveat (true K₂.₇₀ < 0.746). Mass-convention caveat noted
  (model speeds are complex speeds; RQ3 sequential shed at ε ≈ 0
  preserves speed, so the comparison is convention-clean to first order).

### Step 3 — closed form validated at production kinematics; the fate-map scan at measured K

- **The F.2b closed form holds at the new geometry:** companion MD
  n_detect = **5.76** (min 5 / max 8; relaxed 6.16) vs the closed-form
  ladder descent of E₀·e^(−K) = 0.1328 eV → **n = 6**. 0.24 He agreement
  — I34's two-parameter (E₀, K) reduction transfers to production
  kinematics unchanged.
- **The K(R) machinery wiring-oracles cleanly:** the probe-scale curve
  re-derived here matches the Wave-10 §4h table to ≤ 0.007 in K
  (slope 1.226); the production curve is steeper *and lower*:
  K(N) = 0.302/0.399/0.540/**0.746**/1.041/1.449/1.977 for
  N = 250/500/1000/2000/4000/8000/16000 (slope d lnK/d lnR ≈ 1.303).
- **W11-P5 — REFUTED as registered.** In the pre-registered
  E₀ ∈ [0.24, 0.32] no stated prior reaches the bare peak at the
  measured K (best anywhere in-band: 18.4 % bare, W₁ = 2.33 at
  E₀ = 0.32/δ = 0.80). The Wave-10 "landing at 0.28 eV" was an artifact
  of the wrong ballistic bracket.
- **But the construction lands *better than Wave 10* one step up-band
  (still inside RQ1's sourced [0.2, 0.5] eV):** bare crosses 43.5 %
  within E₀ ≈ 0.37–0.43 for **every** stated prior; bare-pinned rows:
  δ = 0.625 → (E₀ = 0.38, bare 42.7 %, n₁ 7.7 %, W₁ 0.785);
  δ = 0.80 → (0.38, 44.3 %, 6.0 %, 0.694); pickup-weighted →
  (0.42, 44.1 %, 7.7 %, 0.723); δ = 0.40 → (0.385, ~47 %, **12.0 %**,
  1.19 — best n₁, worst mid-tail: n₂₋₄ 36 % vs 17.2 %). Global best
  untuned: **W₁ = 0.496** at (E₀ = 0.41, δ = 0.80) — better than
  Wave 10's 0.993.
- **Robustness of the landing to the K uncertainty:** the bare-crossing
  scalar E₀ tracks E\*(K₀) = Σ(21)·e^(K₀); over the *entire* bracketed
  range K₂.₇₀ ∈ [0.49 (ballistic), 0.746 (measured upper-bound)] that is
  E₀ ≈ 0.30 → 0.41 eV — **inside the RQ1 band everywhere**. The
  beyond-band caveat therefore moves the calibration point, not the
  qualitative verdict.
- **The K-scale is near-degenerate in W₁:** even the probe-scale curve
  lands at (E₀ = 0.46, W₁ = 0.570). The histogram alone barely
  discriminates the kinematics — which is precisely why measuring
  K₂.₇₀ (this wave) rather than fitting it was the right move: the
  kinematics question is now settled by measurement, and the scalar E₀
  sweep inherits a *measured* anchor E\* = 0.396 eV.
- **The n = 1 bin worsens:** at bare-pinned production-K points the n₁
  deficit is 1.46× (δ = 0.40) / 2.3× (δ = 0.625) / 2.9× (δ = 0.80) —
  beyond the ~1.4× F.5 electronic-taper cap except at the narrow prior,
  which in turn overshoots n₂₋₄ by 2×. Within this construction the
  n = 1 bin cannot be repaired by (E₀, prior) alone → the sharpened
  RQ4 ladder-bottom question (and/or the tens-of-meV E₀ smear leg)
  inherits it.

### Boundaries of the Wave-11 read

1. **The beyond-band share (42 %) is the load-bearing caveat**: K₂.₇₀ =
   0.746 is measured *inside the locked model*, whose pure-cubic law is
   exercised at up to 2× its calibrated speed band here. Both identified
   biases point up → treat 0.746 as the upper edge and 0.49 as the hard
   lower bracket. A drag law validated at production speeds (TDDFT at
   2.666 Å kinematics, or a velocity-capped form) is the only way to
   shrink this — recorded as a candidate research item, **not** silently
   extrapolated over.
2. The VMI speed comparison is existence/direction-level (complex-speed
   vs detected-fragment-speed convention; single knob point; no
   ensemble). Its direction (model too slow by ~2×) independently favors
   the over-drag reading.
3. n_eject = 21 for all geometries stands (§4h boundary 3); at 2.666 Å
   the two ions start inside one shared first shell (separation < 4.67 Å)
   — the full-dressing convention is maximally strained here, though its
   exposure-side cost is bounded by the 8 % overlap share.
4. The scan inherits every F.3b approximation (pinned-trajectory K(R),
   no-shed leak, ε = 0, RQ3 sequential shed, detected ≈ energetic floor
   at s_eff = 8) — MD-validated at the two measured points (0.84 He at
   9 Å, 0.24 He at production) but analytic in between; N = 50, single
   seed, single droplet radius per dir.
5. All of this is **reported, not auto-adjudicated**: the pre-registered
   W11-P2/P5 are honestly refuted; whether the re-calibrated landing
   (E₀ ≈ 0.38–0.41) supersedes the Wave-10 narrow-E₀ endorsement, and
   whether the drag-recalibration item fires, is the user's call.
   Nothing here discharges the F5 gate.

---

## 4j. H.2b analytic feasibility pass — the bounded lever set cannot land the solvated targets; the miss localizes to the RQ4 taper + the v_c scale; a trapped droplet-retained class appears

Executed 2026-07-11 (under the `[PROCEED TO IMPLEMENTATION]` trigger given
same day; Addendum H §H.2b as frozen, decisions D1–D7). Scratchpad route
`h2b_feasibility.py` — **zero MD, zero repo-code change, no run artifacts
touched**; a self-contained 1D two-body chord forward model importing the
repo's own ladder (`dissociation_ladder`), density gate (`helium_density`),
bundle (b = 2.5153509, E_bind = 0.1167578 eV), erf steepness 14.2 Å, and
R(N) = 2.2173·N^(1/3). 20 000 molecules → 40 000 fragment chords per
exposure bracket (importance-reweighted over all priors × margins; seed
20260711; dt = 0.01 ps, t_end = 150 ps); fate map = the unified F.2b closed
form `E_ej = E₀ᵢ·e^(−Kᵢ)`, suppressed iff `E_ej > Σ(n_eject)`, else exact
ε = 0 descent (detected ≈ energetic floor, the validated s_eff = 8
convention). Scoreboards: `h2b_scan.csv` (7 308 cells),
`h2b_lever_interaction_map.csv`, `h2b_feeder_map.csv` (scratchpad).

### Wiring oracles — the 1D re-implementation reproduces the MD landmarks

| oracle | model | reference |
|---|---|---|
| Σ(21) mixture κ=1 | 0.18783720 eV | 0.18783720 (exact) |
| K production center-pin | 0.74460 | 0.74603 (Wave 11; 0.19 %) |
| K 9 Å center-pin | 0.89767 | 0.89754 (Wave 10/11; 0.014 %) |
| t_exit / v_peak / v_detect (prod) | 4.60 ps / 10.53 / 4.13 Å/ps | 4.62 / 10.51 / 4.11 |
| ballistic v_inf vs analytic | rel. 4.3·10⁻⁴ | (Wave-11 oracle form incl. E_bind) |
| E\*(K) / companion descent | 0.3955 eV / n = 6 | 0.396 / n = 6 (MD 5.76) |

dt-halving moves K by 3·10⁻⁴ — integration-converged.

### Lever-interaction map (the D5 deliverable) — L2 is margin-floored, and the D2 pre-derivation holds

At the code's actual 14.2 Å erf width, `n_eject = round(21·ρ̂)` at the
shallowest allowed birth is **13 / 14 / 15** at margin 3 / 4.67 / 6 Å
(N = 2000 reference: ρ̂ = 0.62 / 0.68 / 0.72). The weighted fraction with
`n_eject ≤ 10` is **0.0 at every (prior, margin)** — there is **no true
surface class** inside the firm band: born-bare is impossible and
Σ(n_eject) ≥ 0.12 eV always. Consequently the **gate-open-at-birth region
is empty in-bounds** for both E₀-law brackets (p = 1 requires
E₀ < Σ(21) = 0.188 eV — below the RQ1 band; p = 0 requires
E₀ < Σ(n_eject) ≤ 0.188 eV — same) — the D2 pre-derivation is confirmed
with the margin floor: **the W12b-P1 "fast n = 1–3 surface feeder" does not
exist under any bounded E₀ law.** The n = 1–3 weight is fed by the
near-cliff chord-K band instead (feeder map: n = 1 comes from births at
depth −6…−10 Å, n_eject 15–17, K ≈ 0.17–0.29 — just above the cliff).
Shell-averaged-ρ̂ sensitivity: with shell radius 4.67 ≪ steepness 14.2, the
averaging correction is ≤ (4.67²/6)·|ρ̂″| ≲ 0.01 in ρ̂ → < 0.25 atoms —
sub-rounding, dispositioned analytically.

### The scan — no landing at the frozen bar, at either bracket

7 308 cells (E₀ ∈ [0.20, 0.50] × p ∈ {0, 1} × 4 ladders + the RQ4
diagnostic × 4 priors × 3 margins × 2 exposure brackets): **zero cells pass
T1∧T2∧T3.** The ballistic bracket is *worse* everywhere (best W₁ = 1.94 —
under-cooled: bare-heavy, flat solvated), so the verdict is **not
law-conditional** (the D3 protection). Best cells per ladder (current law):

| ladder | best cell | W₁ | n₁ | n₁/n₂ | bare | trapped |
|---|---|---|---|---|---|---|
| flat | d080, m3, p1, E₀ 0.21 | 0.989 | 0.160 | 1.14 | 0.05 | 0.080 |
| floor1 (knob-free X₂) | d080, m3, p1, E₀ 0.22 | 0.813 | 0.223 | **1.83** | 0.11 | 0.080 |
| slid2 | d080, m3, p1, E₀ 0.22 | 0.725 | 0.214 | 1.25 | 0.06 | 0.080 |
| slid3 | d080, m3, p1, E₀ 0.23 | 0.722 | 0.214 | 1.29 | 0.12 | 0.080 |
| **rq4graded (diagnostic)** | d080, m6, p1, E₀ 0.25 | **0.575** | **0.267** | **1.76** | 0.13 | 0.065 |

Targets: W₁ ≤ 0.5, n₁ ∈ [0.26, 0.36], ratio ∈ [1.75, 2.6]. Structural
reads:

1. **p = 1 wins everywhere; p = 0 dies by over-suppression** — exactly the
   D2 pre-derivation (constant E₀ over shrunken Σ(n_eject) makes the
   off-center population suppressed → bare-heavy). Under p = 1 the descent
   target scales with Σ(n_eject), so the margin-floored dressing *does*
   feed small-n — just not n = 1 specifically.
2. **The n₁ bin is ladder-taper-controlled, and sliding the X₂ transition
   up makes it WORSE** (slid2/slid3 ratio 1.25/1.29 vs floor1's 1.83):
   n = 1 is a one-rung-wide window in E_ej, so equal deep rungs widen the
   n = 2/3 windows equally and the *ratio* reverts to the (flat) chord-K
   density ratio ≈ 1.1–1.3. What steepens n₁/n₂ is a genuine **taper**
   (rung 1 > rung 2 > rung 3). The knob-free X₂ floor (1.44×) reaches
   1.83; the experimental 2.18 needs ≈ 2× — beyond the cap, as
   pre-registered.
3. **The RQ4-graded diagnostic (2.2 : 1.5 : 1.3 — outside the bounded set,
   reported as a prediction FOR the external calculation) nearly lands:**
   T2 + T3 pass, bins 2–8 match to ≲ 0.01–0.02 absolute
   (pickup/m3/E₀ 0.23: 0.139/0.142, 0.105/0.092, 0.074/0.071, 0.066/0.058,
   0.056/0.049, 0.055/0.039, 0.054/0.036; bare 3.8 %), W₁ = 0.575–0.582 vs
   the 0.5 bar. Through a purely geometric forward model, **the solvated
   histogram now independently demands the RQ4 target ratios.**
4. **T4 (KE): the envelope holds at every top cell** — the experimental
   mean-KE curve lies inside [current-law, ballistic] for all n = 1–12 —
   but the single-λ scale fails (worst bin ×1.83–2.36 vs the ×1.5
   tolerance). The needed lift is **speed-selective in the v_c direction**:
   ×2.2 at n = 1 falling to ×1.7 at n = 12 (speed ×1.5 → ×1.3), and the
   λ-blend is a conservative proxy (a real velocity-capped law lifts fast
   fragments and leaves the slow n ≥ 13 tail at its in-band current-law
   values, where the blend artificially fails the < 0.1 eV band). One
   marginal band violation is real: the small-weight near-cliff fast class
   at n ≈ 19 arrives at 0.122 eV.
5. **E₀ re-lands at the solvation scale.** The optima sit at
   E₀ ≈ 0.22–0.27 eV ≈ E_solv — not the pinned-droplet 0.38–0.41 (Wave 11)
   — because the position axis supplies the low-K mass that the E₀ scan
   previously had to buy with a higher cliff. With bare un-targeted, the
   narrow-E₀ story becomes *more* physical under the full geometry.
6. **A trapped, droplet-retained ion class appears (new).** 6.5–11.3 %
   (weighted; prior/margin-dependent) of fragments — inward-going partners
   of off-center births with chord exposures up to K ≈ 18 — are dissipated
   by the current law below the 0.117 eV solvation barrier and **never
   leave the droplet** (parked near the far surface where the residual
   Coulomb push ≈ the well gradient). They are excluded from every
   detected read (they would appear as large cluster masses, not I⁺Heₙ
   bins). The class is t_end-conditional (150 ps read; the park is
   quasi-static) and ballistic-bracket-absent — a current-law prediction,
   MD-arbitrable by the W12 leg.

### Verdict (frozen outcome shapes)

**Outcome (c) at the strict bar — reachability inside the all-bounded
lever set is refuted, at both exposure brackets.** But the miss is not
diffuse; it decomposes into exactly two quantified, already-named items:

- **the RQ4 ladder taper** (the histogram side: the bounded floor variant
  reaches ratio 1.83 / W₁ 0.81; the RQ4-graded ratios close bins 2–8 and
  reach W₁ 0.575) — the external many-body calculation is confirmed as the
  **blocking arbiter**: if it returns ≈ 2 : 1.4 : 1.2 or steeper, the
  histogram closes inside physics; if it returns a plateau, the F.5 escape
  clause (non-zero ε / missing mechanism) fires with the geometric
  alternatives now exhausted;
- **the W13 v_c speed scale** (the KE side: envelope ✓, needed lift
  ×1.3–1.5 speed-selective — W13-P1's premise pre-confirmed with numbers).

The slid-X₂ family is **rejected as the repair** (worsens the ratio); the
W12b dressing arm's pre-registered fast-feeder fingerprint (W12b-P1) is
**predicted absent** under both bounded E₀ laws; the W12/W12b MD legs'
role shrinks from exploration to verifying the geometric components the
pass says are load-bearing (the chord-K density near the cliff, the
trapped class, the KE shape). Reported, not auto-adjudicated — sequencing
is the user's call; nothing here discharges the F5 gate.

### Boundaries of the H.2b read

1. Pinned-mass straight-line chords (no re-pickup, no bending; the
   measured m(t)-feedback scale on K is −2.2 %, Wave 11); the no-shed leak
   is baked into the unified `E_ej = E₀ᵢe^(−K)` form (−1.7 He mid-shell
   bias at early openings, ≈ 0 at the cliff).
2. Detected ≈ energetic floor (s_eff = 8 convention; ≤ 1.4 He
   s-dependence, I20); RQ3 spec-(b) suppressed → bare; ε = 0.
3. The trapped classification is a 150 ps operational read of a
   quasi-static park; leak-out over µs flight is not modeled.
4. The λ-blend T4 operationalization is a conservative stand-in for a
   velocity-capped law (per-fragment v_c physics is W13's).
5. Analytic/MC read (no N = 50 quantization; MC noise ≪ the 0.075+ W₁
   margins); single master seed; mean-⟨N⟩ = 2000 log-normal convention
   stated (D4).
6. The KE targets are the D7 user-supplied means, provenance-pending (the
   H.2 export prerequisite stands).

---

## 4k. Addendum I Steps 0–1 — the corrected reference re-draws the T4 read: envelope holds everywhere, the required lift is near-uniform ×1.6–2.2 (n ≥ 2) with an n=1-only excess

Executed 2026-07-15 under `[PROCEED TO IMPLEMENTATION]` (Addendum I as
frozen, decisions I-D1–I-D6). Zero MD, zero repo-code change. The recovered
H.2b forward model (`h2b_feasibility.py` + cached 40 000-fragment master
ensemble, seed 20260711) revalidated first — **all six wiring oracles
reproduce exactly** (Σ(21) = 0.18783720 eV; K_prod = 0.74460; K_9Å = 0.89767;
t_exit/v_peak/v_detect = 4.60/10.53/4.13; ballistic rel 4.3·10⁻⁴;
E\*(K) = 0.3955 / n = 6). The Addendum-I tail integrator reproduces the
cached current law at v_c = ∞ to machine precision (max |ΔK| = 4.6·10⁻¹⁴
over a 2 000-molecule subset, all four tails). Driver:
`addendum_i_steps01.py` (scratchpad); scoreboards
`addendum_i_step0_envelope.csv`, `addendum_i_sweep.csv`.

### Step 0 — D7 retirement (the T4 targets were the legacy moment convention)

The D7 table is **retired**. Corrected reference
(`data/reference/ihe_ked/IHe_KED_reference.csv`) vs the D7 values the H.2b
T4 verdict was scored against:

| n | D7 [eV] | corrected [eV] | ratio | σ_pt [eV] |
|---|---|---|---|---|
| 0 | 2.900 | 3.7057 | ×1.278 | 0.038 |
| 1 | 0.974 | 1.3017 | ×1.336 | 0.071 |
| 2 | 0.545 | 0.7058 | ×1.295 | 0.027 |
| 3 | 0.390 | 0.4964 | ×1.273 | 0.020 |
| 6 | 0.248 | 0.2670 | ×1.077 | 0.003 |
| 9 | 0.154 | 0.1635 | ×1.061 | 0.004 |
| 12 | 0.105 | 0.1070 | ×1.019 | 0.004 |

The correction is **not** the coherent +25–32 % everywhere: it is ×1.27–1.34
for n ≤ 3 but shrinks to ×1.02–1.08 by n ≥ 10 (the D7 tail values were
already near the corrected measure). Consequence: the D7-based
"speed-selective ×2.2 → ×1.7" lift reading (I43) had the right sign but the
wrong shape — see below. Error treatment from here on (I-D4): per-point
stat ⊕ sys; the 4 % calib and 6 % condition bands profiled as one coherent
factor f ∈ [0.902, 1.102].

### Step 0 — envelope verdict: (c′) does NOT fire

Self-consistent brackets (binning *and* speeds from the same law) at the
five §4j best-per-ladder cells: **every bin n = 1…17 lies inside the
[current-law, ballistic] envelope at every cell, already at f = 1** — the
ballistic curve sits at 1.71–2.15 eV across bins (≥ ×1.3 above the n = 1
target), the current-law curve at 0.33–0.50 eV (n = 1) falling to
0.03–0.06 eV (n ≥ 13). A drag-form solution is not excluded anywhere;
Addendum I proceeds to Step 1. **I-P1 CONFIRMED** (n = 1 inside the
ballistic bracket without needing the coherent bands).

Corrected required-lift curve (exp / current-law KE, f = 1):

- flat cell: ×2.60 at n = 1, ×1.57–1.94 for n = 2…17 (near-uniform);
- floor1 / slid2 / slid3: ×2.8–3.2 at n = 1, ×1.7–2.2 for n ≥ 2;
- rq4graded: ×4.0 at n = 1, ×2.0–2.6 for n ≥ 2.

Current law vs the I-D5 bar: **0/12 bins pass at every cell** (worst
×2.34–×3.61 after profiling f) — the corrected reference makes the
current-law undershoot *worse* than the stale §4j read (×1.83–2.36).
The n ≥ 13 band (< 0.1 eV) holds for the current law at all cells.

**Supersessions (D7-conditional numbers, corrected here):** §4j structural
read 4 ("T4: envelope holds … needed lift ×2.2 at n = 1 falling to ×1.7 at
n = 12") and the §4j verdict's "W13 v_c speed scale ×1.3–1.5" premise are
superseded: the corrected lift is **×2.6–4.0 at n = 1** and **near-uniform
×1.6–2.2 for n = 2…17** (cell-dependent). The "speed-selective" shape
claim weakens — selectivity concentrates in the *single* n = 1 bin, the
rest of the curve needs an almost flat lift. §4j boundary 6 (KE targets
provenance-pending) is discharged. I43 carries a supersession note.

### Step-1 pre-registered predictions (registered 2026-07-15 before any sweep read-out; sweep launched, output unread)

- **I-P2 (ordering/monotonicity):** per-bin detected mean KE increases
  monotonically as v_c decreases for every tail, and at fixed v_c the tail
  ordering is KE(cut) ≥ KE(p = −1) ≥ KE(p = 0) ≥ KE(p = 1) ≥ current law
  (pointwise γ ordering; strict at the center pin, approximate per bin
  under binning reshuffles).
- **I-P3 (the n = 1 residual — the sharp discriminator):** no single
  (p, v_c) supplies the n = 1 excess (×2.6–4.0) *and* the flat n ≥ 2 lift
  (×1.6–2.2) simultaneously: the best cells land n = 2…12 and leave n = 1
  low, **unless** the weaker law's binning shift re-feeds n = 1 from
  faster near-cliff chords. If a tail passes all 12 bins including n = 1,
  W13-P1's "v_c owns the scale" reading strengthens; if n = 1 alone fails
  everywhere, the n = 1 class carries extra physics (RQ8-adjacent
  channel mixing or dressing) and the miss is *localized*, not a drag
  failure.
- **I-P4 (band):** any tail that lands n = 2…12 keeps every populated
  n ≥ 13 bin below the 0.1 eV band (the implied ~×1.7 lift takes the
  0.03–0.06 eV current-law tail to ≤ 0.09 eV).

### Step 1 — tail-family sweep: no landing (0 / 10 416); the scale is fixable, the slope is not; KE and histogram anti-correlate through K

28 tail integrations (p ∈ {1, 0, −1, cut} × v_c ∈ {5.3, 6, 7, 8.5, 10, 12,
15}) × the focused fate-cell subset (floor1/rq4graded × d080/pickup ×
margins {3, 4.67, 6} × E₀ 0.20–0.50, p_couple = 1) = **10 416 scored
cells** (`addendum_i_sweep.csv`).

1. **Zero cells pass the I-D5 bar.** Max per-cell pass count is 9 — and
   every ke_pass ≥ 8 cell is *degenerate*: it *empties* the n = 1 bin
   (model weight ≤ 0.01 % vs the experimental 31 %, dropping n = 1 out of
   scoring via the 0.2 % weight floor), pins E₀ at the 0.20 grid edge
   (87/89 cells) and the coherent factor at its 0.902 lower edge (88/89),
   and destroys the histogram (W₁ = 1.31–2.36). The weight floor's
   bin-evasion loophole is exposed and reported — those cells are not
   candidates.
2. **The global scale is deliverable; the slope is not.** The genuine
   compromise cell (p = 0, v_c = 6.0, floor1, pickup, m6, E₀ 0.21) lands
   **n = 1…6 inside ×1.25** (n = 1 at ×0.87 with 20 % bin weight — the
   binning shift does re-feed n = 1) with W₁ = 0.845, n₁ = 0.200,
   ratio 1.32, bare 10⁻³ — but n ≥ 7 undershoots progressively:
   ×0.81 (n = 7) → ×0.49 (n = 12) → ×0.26 (n = 17). The experimental
   ⟨E⟩(n) decays much more slowly than any model curve: deep bins select
   the most-exposed chords (K-selection) and their late, slow dynamics
   ride the **in-band** cubic law — which is TDDFT-locked and untouchable
   by every admissible tail. The miss localizes to the deep-bin slope,
   n ≥ 7.
3. **KE and histogram anti-correlate structurally.** Every KE-improving
   tail weakens the exposure (center-pin K 0.745 → K_q50 0.18–0.37),
   which starves the evaporative descent feeding small n: among ke_pass
   ≥ 8 cells W₁ ≥ 1.31 (vs the current-law best 0.575–0.813); among
   W₁ ≤ 0.9 cells ke_pass ≤ 7. Under the frozen fate map
   (E_ej = E₀·e^(−K), ε = 0) and the pinned cooling clock (τ = 6.55 ps),
   *both* observables ride the single integral K — the drag lever cannot
   serve them simultaneously.
4. **The hard cutoff is excluded as a physical law** (as designed — it
   was the diagnostic bracket): at v_c ≤ 10 it makes the ensemble
   near-ballistic (K_q50 ≈ 0.18) → 99.7 % suppressed → bare, W₁ ≈ 3.67;
   at v_c = 12 > peak speed it is numerically identical to the current
   law.
5. **E₀-edge tension (recorded):** the KE fit pulls E₀ to (and against)
   the 0.20 eV RQ1-band bottom, while the histogram-good cells sit at
   0.21–0.25 — a second, milder expression of the same K-coupling.

### Prediction verdicts

- **I-P1 CONFIRMED** (Step 0: n = 1 inside the ballistic bracket at f = 1).
- **I-P2 CONFIRMED** — center-pin v_inf strictly monotone in v_c per tail
  and strictly ordered cut ≥ p = −1 ≥ p = 0 ≥ p = 1 ≥ current law at every
  v_c (e.g. v_c = 6: 15.36 / 10.01 / 5.45 / 4.36 / 4.13 Å/ps); K_q50
  ordering identical. The drag lever is well-behaved — the failure is
  structural, not a search artifact.
- **I-P3 CONFIRMED, first branch, with a twist:** no (p, v_c) supplies
  the n = 1 excess and the flat n ≥ 2 lift simultaneously. The
  anticipated binning re-feed of n = 1 *does* occur (compromise cell:
  n = 1 passes at 20 % weight) — so the "n = 1 alone fails everywhere"
  alternative is **refuted**; the localized residual is the **deep-bin
  slope (n ≥ 7)**, not n = 1.
- **I-P4 vacuously satisfied** (no cell landed all of n = 2…12);
  band_n13_ok held at every near cell — consistent, untested. Noted: the
  < 0.1 eV band bar is generous relative to the actual n ≥ 13 reference
  values (0.066–0.096 eV) — the fit metric sees the deep-bin undershoot
  the band bar cannot.

### Verdict (Addendum I outcome shapes) — *superseded same day by the Step-1b verdict below (joint closure); retained per program convention*

**Outcome (b): scale lands, shape fails, for every admissible tail.**
Step 2's entry condition (a bar-passing cell) is **not met — no MD build
is recommended from this sweep.** The conditional shortlist, *if* the
K-decoupling question (OQ-I) is resolved in favour of a joint
recalibration: **p = 0, v_c ≈ 6 Å/ps** (linear-force tail; best
two-observable cell) with **p = −1, v_c ≈ 7–8.5** as the alternate.
Reported, not auto-adjudicated — sequencing is the user's call; nothing
here discharges the F5 gate.

### Step 1b — joint (v_c, τ) mini-sweep (§I.8, user-approved same day): JOINT CLOSURE at τ = 4.1 ps on the rq4graded ladder — outcome (a)

Executed 2026-07-15 immediately after Step 1 (design + predictions
J-P1–J-P3 registered in plan §I.8 pre-execution). Mechanics: τ enters the
dynamics nowhere, so K(τ) = K(6.55)·(6.55/τ) per chord — 13 cached tail
integrations × a free (τ, E₀, cell) rescan: **28 899 cells**
(`addendum_i_joint.csv`; τ ∈ {6.55, 5.2, 4.1, 3.3, 2.6, 2.0}).

1. **Joint closure exists and is contiguous: 18 cells pass the full
   I-D5 KE bar (12/12 within ×1.25) AND W₁ ≤ 0.85 — all at τ = 4.1 ps,
   all on rq4graded**, spanning (p = 0, v_c = 6–7) and (p = −1,
   v_c = 7–8.5), both priors, all three margins, E₀ = 0.22–0.26
   (interior — the Step-1 grid-edge pinning resolves). KE-bar passes
   (any W₁) span τ ∈ {3.3, 4.1, 5.2} peaked at 4.1 (2/74/40).
2. **One cell lands everything the program has ever scored:**
   (p = −1, v_c = 7, τ = 4.1, rq4graded, pickup, m3, E₀ = 0.25) passes
   **T1 ∧ T2 ∧ T3 ∧ the 12/12 KE bar**: W₁ = 0.272 (vs the 0.5 bar and
   the H.2b in-bounds best 0.575), n₁ = 0.280 ∈ [0.26, 0.36],
   n₁/n₂ = 1.861 ∈ [1.75, 2.6], worst KE ratio ×1.243, deep-bin worst
   ×0.804⁻¹, bare 0.149, trapped 0.054, f = 0.916. The first full
   two-observable landing in the Tier-2 probe program — where H.2b
   concluded in-bounds unreachability, opening one pinned knob (τ) plus
   the bounded tail closes both observables *on the RQ4-anticipated
   taper*.
3. **Both knobs are load-bearing.** The v_c = 15 current-law control
   never joint-closes at any τ (W₁ ≥ 1.53 among KE-decent cells) — τ
   alone cannot do it; Step 1 showed v_c alone cannot either. The
   KE–histogram anti-correlation (I47) is *resolved by the pair*, exactly
   the OQ-I arm-(a) mechanism: v_c owns the speed integral, τ re-owns
   the descent clock.
4. **Form discrimination emerges** (the I-D4 curve-fit upgrade paying
   off): p = −1 wins the histogram ratio (T2/T3 pass) but its far tail
   sags (n ≥ 13 KE ratios 0.59–0.78); p = 0 gives the flattest KE shape
   ever seen in the program (ratios 0.85–1.04 across *all* n = 1…17 at
   its standout cell, deep bins fully recovered) but tops out at ratio
   ≈ 1.66 on T2. The two shortlist forms are experimentally
   distinguishable — by the n ≥ 13 mean-KE tail and the n₁/n₂ ratio.
5. **The deep-bin slope was K-coupling after all:** at the matched clock
   the deep bins recover to ×0.80–1.04 — the Step-1 undershoot was
   binning selection (weak-K chords flooding deep bins), not
   in-band-locked late dynamics. OQ-I arm (c) does not need to open.

**Prediction verdicts:** **J-P1 CONFIRMED** (matched clock: the
W₁-minimizing τ per (p, v_c) tracks 6.55·K_q50(tail)/K_q50(current),
τ\* = 4.1 at (p = 0, v_c = 6) as computed; E₀ optima interior; W₁
recovers *below* the predicted 0.58–0.85 — to 0.20–0.43).
**J-P2 REFUTED in its conservative branch** — the deep bins do clear the
bar at the matched clock (worst ×0.80–0.91); informative: the slope
belonged to K-coupling, not to the locked in-band law.
**J-P3 CONFIRMED** (no joint closure without a tail; KE per fixed
binning τ-flat by construction).

### Verdict (superseding the Step-1 verdict above): Step-1b outcome (a)

**The W13 Step-2 MD build re-opens with a two-knob (v_c, τ) arm and a
defined target region** — (p = 0, v_c ≈ 6–7) and (p = −1, v_c ≈ 7–8.5)
at τ ≈ 4.1 ps, E₀ ≈ 0.23–0.26, rq4graded ladder. Three adjudications
belong to the user before any build: (i) **the τ re-classification
event fires as pre-registered** — τ = 4.1 ps is ×0.63 of the
GAH25-sourced 6.55 ps pin; the pin becomes a fitted knob and its
sourcing must be re-argued or reclassified (CALIBRATION_MAP propagation
at build time); (ii) **the closure is rq4graded-conditional** — the
ladder taper remains the RQ4 external calculation's to confirm; the
stakes sharpen: with (v_c, τ) open, the taper closes *both*
observables, and a returned plateau now breaks the KE side too;
(iii) **form choice** — implement both shortlist tails for MD
discrimination, or pick one (the n ≥ 13 tail + n₁/n₂ discriminate).
Reported, not auto-adjudicated; nothing discharges F5.

### Step 1c — closure-basin refinement (§I.9): the closure is a plateau; refined targets (v_c, τ) ≈ (7.5, 3.8) / (6.0–6.5, 3.8–4.0); the taper buys the full histogram bar, not the joint closure itself

Executed 2026-07-15 (design + K-P1–K-P3 registered in §I.9 pre-execution;
RQ9 opened and parked in parallel per user decision). 8 new tail
integrations (fine v_c grid) + free τ-rescale: **39 312 cells**
(`addendum_i_basin.csv`; τ = 3.0–5.6 step 0.2, E₀ = 0.20–0.32).

1. **The closure is a calibratable plateau, not a knife-edge
   (K-P1 CONFIRMED):** 225 joint closures. Per form: p = 0 basin
   v_c ∈ [6.0, 7.5] × τ ∈ [3.2, 4.8] × E₀ ∈ [0.21, 0.25] (96 cells);
   p = −1 basin v_c ∈ [7.0, 8.5] × τ ∈ [3.0, 4.8] × E₀ ∈ [0.21, 0.27]
   (113 cells) — τ width ≈ 1.6 ps, v_c width ≈ 1.5 Å/ps, both an order
   wider than predicted minimums.
2. **24 full-house cells** (joint KE ∧ T1 ∧ T2 ∧ T3 — the Step-1b
   single cell was the edge of a region, not a fluke), τ ∈ [3.4, 4.2],
   all rq4graded, mostly d080 prior / 3 Å margin. Best:
   (p = −1, v_c = 7.5, τ = 3.8, d080, m3, E₀ = 0.25) with **W₁ = 0.196**,
   n₁ = 0.289, ratio 1.861, 12/12 KE (worst ×1.225), bare 0.137.
   **Refined calibration targets: p = −1: (v_c, τ) ≈ (7.5, 3.8);
   p = 0: (6.0–6.5, 3.8–4.0)** — the τ center shifts 4.1 → ≈ 3.8–4.0
   (the Step-1b grid quantization). ×0.58–0.61 of the GAH25 pin (RQ9
   numbers update accordingly).
3. **K-P2 CONFIRMED:** p = 1 (Newton tail) — zero joint closures on the
   fine grid; the form family truncates to {p = 0, p = −1}.
4. **K-P3 REFUTED for the joint criterion — a finding:** floor1 (the
   knob-free, in-bounds X₂ ladder) holds **16 joint closures**
   (12/12 KE + W₁ down to 0.575) at slightly larger τ (4.0–4.8) and
   lower E₀ (0.21–0.23). But **zero full houses**: floor1's n₁ tops at
   0.23 (< T3's 0.26) and its ratio at 1.83 only marginally.
   Consequence — the RQ4-conditionality statement *refines*: **bounded
   physics + (v_c, τ) now reaches the entire KE curve plus a
   W₁ ≈ 0.58 histogram** (the level H.2b's best achieved with *no* KE
   match); the RQ4 taper is what buys the last stretch — the n₁
   weight, the n₁/n₂ ratio, and W₁ 0.58 → 0.20. A returned RQ4 plateau
   no longer breaks the KE side (Step-1b's sharpened stake is walked
   back to: it breaks T2/T3, as in H.2b).

### Boundaries of the Steps 0–1(b,c) read

1. All §4j 1D-model boundaries carry verbatim (pinned-mass straight
   chords, no-shed leak, detected ≈ energetic floor s_eff = 8, ε = 0,
   150 ps trapped read, mean-⟨N⟩ = 2000 priors).
2. The v_c grid is coarse (7 points; Step 1b adds none) and E₀ stops at
   the RQ1-band bottom 0.20 — the sweeps bracket, they do not optimize;
   conclusions rest on the bracketing pattern, not on any single best
   cell.
3. The KE score is mean-to-mean per the export convention; the n = 0…4
   P(E) curve shapes were not scored (validation-only per I-D4).
4. p_couple = 0 was not re-swept (H.2b: dies by over-suppression;
   weakened drag only worsens that direction).
5. The joint closure's deep-bin passes hug the ×1.25⁻¹ edge (0.80–0.91)
   — MD confirmation can break them; the T2 pass sits at 1.861 in a
   [1.75, 2.6] band (margin ≈ 6 %); the τ grid quantizes at 4.1
   (neighbors 3.3/5.2 hold KE-bar passes but not joint closures — the
   τ landing is sharper than the v_c one).
6. τ rescaling is exact within the model (τ appears only in the
   exposure bookkeeping); any physics in which τ feeds back on the
   dynamics (e.g. temperature-dependent drag) is outside this read.

---

## 4l. Slice-T3 pilot first look — the undressed MD parks at n ≈ 7–9; the closure miss localizes to the ensemble geometry, not the calibrated law

Read 2026-07-16 (scratchpad first-look over the four delivered
`…_tier2probe_conf270_c{1..4}` detection artifacts; delivery record:
log entry "Slice T3 DELIVERED"; the formal ihe_ked scoring pass has not
run — this is the structural read that re-sequences it).

**Result.** All four configs produce a **narrow mid-shell cluster with
zero weight below n = 5** (100 ions each, 90–95 % `frozen`, 0 %
suppressed — E₀ ≥ 0.23 eV opens the gate in-bubble at every birth):

| | n range | peak bins | bare | solvated n₁ | n̄_detect |
|---|---|---|---|---|---|
| C1 | 6–12 | n8 0.34 / n9 0.40 | 0.00 | **0.00** | 8.63 |
| C2 | 5–12 | n8 0.31 / n9 0.43 | 0.00 | **0.00** | 8.77 |
| C3 | 6–10 | n6 0.29 / n7 0.50 | 0.00 | **0.00** | 6.97 |
| C4 | 6–11 | n8 0.42 / n9 0.41 | 0.00 | **0.00** | 8.39 |

S2-P3 (the rq4graded/floor1 n₁ split) is **structurally silent** —
rungs 1–3 are never reached; the experimental targets (solvated
n₁ = 0.31, n₁/n₂ = 2.18) are out of reach at *any* (drag, τ, ladder)
value of the matrix.

**Diagnosis.** The Step-1b/1c closure cells were scored on the full
H.2b ensemble — L1 birth margin (3–6 Å), **L2 depth-dressing
n_eject(d) = round(21·ρ̂)**, the D2 E₀-coupling p = 1, and the D4
droplet priors. The repo MD carries only the position axis:
`ion_initial_state.py` starts **every** biphasic ion at the full
n₀ = 21 shell (`ANCHOR_N_START` — the Tier-1a validation convention),
so off-center births change exposure only, every ion must descend the
full Σ(21) ladder, and E_ej buys only ~12 rungs → the n ≈ 7–9 park.
The §I.10 margin-0 "recorded caveat" is hereby **measured to be
load-bearing**, and the H.3b park premise ("building the dressing would
confirm an inertness prediction") is overturned in the direction its
own revival clause reserved.

**What does transfer.** Cross-config ordering follows the F.2b
race/leak picture exactly (shorter τ → larger K → smaller E_ej → fewer
post-ejection sheds → higher terminal n: 8.63 (τ3.8) > 8.39 (τ4.4) >
6.97 (τ6.55)); the capped tails lift detected KE as designed (C1/C2/C4
×1.3–2.0 vs the experimental means in their populated bins, current-law
C3 ×0.85–1.15) — though the bin-level KE read is geometry-conditional
until the shell question is fixed. All four stages ran clean at
production kinematics (first in-repo), and the dt-halving spot check
passed (drifts ≤ 0.30× SEM).

Consequence: the **Step-2c geometry-closure slice plan** (plan §I.11,
slices T5–T9) — build the missing ensemble levers in MD and re-run the
confirmation under them. The §I.10 T4 scoring gate is absorbed there.

---

## 4m. T9 leg A′ — the position axis alone opens the two-sided race; twin-parity is quantitative on the class/histogram axis; the KE composition diverges through the real cascade

Executed 2026-07-16 (T9 oracle-chain leg A′, all four C-configs, under
the leg trigger; delivery record: log entry "T9 leg A′ EXECUTED").
Exactly **one physics lever** flipped vs the delivered T3 dirs:
`birth_position_law = "uniform_volume"` at the Step-1c full-house
margin 3 Å (Slice T7) — plus the two consequences the lever forces
(the V0-2 `exclude` retained policy; relaxation cap 8000 ps). Twin
predictions were **pre-registered before any MD** (AP-P1..P4 + CSVs;
log entry of the same day). Dirs: `…_tier2probe_conf270_apc{1..4}`.

### The A/B (N = 50 → 100 ions; twin m = 20000 chords)

Conventions: `droplet_retained` excluded from the read (trapped
column, fraction of all ions); MD `suppressed` counted at bin 0 (the
RQ3 spec-b bare-candidate class, exactly the twin's fate-map
bookkeeping); W₁ over the full 0–21 histograms [bins].

| config | trapped (twin/MD) | supp→bare (twin/MD) | n₁ (twin/MD) | n̄_det (twin/MD) | W₁ |
|---|---|---|---|---|---|
| c1 | 0.033 / 0.050 | 0.089 / **0.095** | 0.172 / **0.147** | 6.22 / 5.80 | 0.69 |
| c2 | 0.035 / 0.060 | 0.035 / **0.032** | 0.170 / **0.160** | 6.55 / 5.97 | 0.69 |
| c3 | 0.067 / 0.110 | 0.336 / **0.315** | 0.069 / **0.079** | 4.48 / 4.12 | 0.55 |
| c4 | 0.033 / 0.050 | 0.118 / **0.116** | 0.138 / **0.126** | 6.05 / 5.75 | 0.58 |

### Verdicts (pre-registered reads)

1. **AP-P1 CONFIRMED.** Every MD histogram broadens from the T3 ±2-bin
   park (I52: zero weight below n = 5) to the twin's near-full-range
   two-sided shape; n₁ is the top solvated bin on all rq4graded/floor1
   configs — the first MD weight ever at n = 1 and n = 0.
2. **AP-P2 CONFIRMED, sharpened.** The twin−MD offset is −0.3 to
   −0.6 He (twin slightly high) — smaller than leg A's +1.4 — and the
   *shape* agreement is W₁ = 0.55–0.69 bins over 22 bins.
3. **AP-P3 CONFIRMED.** The droplet-retained class is real in MD in
   every config (5/6/11/5 ions), ~1.6× the twin's 150 ps chord read,
   with c3 the largest on both sides.
4. **AP-P4 CONFIRMED.** The suppressed/bare ordering
   c3 ≫ c4 > c1 > c2 transfers **exactly** (twin 0.336/0.118/0.089/
   0.035 → MD 0.315/0.116/0.095/0.032) — the Δ×-race weight physics
   crosses the 1D→3D boundary quantitatively.

### The divergence — histogram parity is NOT composition parity

The one axis that does **not** transfer is per-bin mean detected KE at
small n: MD n₁ ≈ 2.42–2.48 eV (near-ballistic; per-fragment budget
2.70 eV) vs twin 0.97 (c1) — ×2.5; c3: MD 0.666 vs twin 0.261. In the
twin's ε = 0 fate map, a bin's population and its speed ride the same
chord exposure; in the MD, the real RRK cascade + per-shed kicks +
live E_int dynamics let *different trajectories* feed the same bin.
Same histogram, different occupants — visible only on the KE axis.
This is the S2c-P4 *listed* fate-map↔real-cascade channel firing, not
a model-structure surprise, but it is load-bearing: the (n, mean-KE)
curve is composition-sensitive, so KE reads off undressed A′ bins are
geometry-conditional, and leg B (T5 dressing) re-decides which
trajectories feed small n. (The experimental n₁ mean 1.302 eV sits
between the undressed-MD 2.48 and the twin 0.97 — noted, not scored.)

### Structural findings of the execution itself (all guard-caught)

1. **The trapped class collides with the detection decoupling
   contract** → the `detection_droplet_retained_policy` arm (V0-2
   convention; `refuse` default byte-inert).
2. **Marginal ions fly conservative scattering orbits in E2** (the
   relaxation "coulomb" mode is zero-gamma): apc3's "ion 9"
   (+1.14 meV radial margin) turned out to be a **centrifugal
   resonance** — above the radial threshold, below its effective-
   potential barrier, drifting ~5 Å per 3000 ps. No relaxation cap
   resolves such an ion; the `exclude` bound criterion was upgraded to
   the exact conservative-mechanics condition
   (E_tot < max_path [U + L²/2mr′²]).
3. **The E2 zero-gamma convention is an open physics question**: a
   real ion orbiting at ρ ≈ 0.1–0.5 feels drag and would be captured;
   the conservative convention makes the resonance long-lived instead.
   Recorded for user adjudication (a drag-live E2 arm) — required
   before T9's N = 500, where marginal ions are guaranteed.

### Boundaries

N = 50 single seed (class fractions quantize at 1 %); one detection
RNG realization per dir; the undressed KE read is geometry-conditional
(above); cap 8000 ps is leg-specific adequacy, not a pin; the twin
comparison is at p_couple = 0 / undressed n_eject = 21 by construction
— nothing here reads on the experimental targets (that is the
re-centered T9 re-pilot's job); nothing discharges F5.

---

## 4n. I.11.2 item 1 — the n₁-composition re-read: the KE excess is the cold-shed momentum convention, not the Coulomb share; the twin−MD KE gap closes under a co-moving counterfactual (OQ-J fired)

Executed 2026-07-17 (plan §I.11.2 item 1; zero new MD, scratch read over
the four on-disk apc dirs; delivery record: log entry "I.11.2 item 1
EXECUTED"). Method: per ion — birth radius, chord cosine, first-shed time,
droplet-exit time, mass-at-exit — plus the **exact per-fragment Coulomb
work** `W_i = ∫ F_coul·v_i dt` from the stored trajectories (pair-sum
renormalized against the identity `W_A + W_B = U_c(0) − U_c(end)`;
post-window residual split by the impulse rule). The per-ion detected-KE
ledger then closes to ≤ 0.06 eV on every live class:

`KE_det = E_sym + coul_x − drag + boost + resid`, with `E_sym = 2.70 eV`
(symmetric half of the 5.40 eV pair release), `coul_x = W_i − E_sym` (the
mass-asymmetric share), `drag = E_dissip_detected`,
`boost = −E_mass_transfer_detected` (the net KE injected by
momentum-conserving mass events), `resid ≈ well work + numerics`.
(Empirical ledger identity verified first:
`E_kin + E_pot + E_dissip + E_mass_transfer + E_int` is conserved on the
stored series — mean end-drift −0.0008 eV.)

### The decomposition (c1; c2/c4 within a few 10 meV throughout)

| class | N | KE_det | v_det | coul_x | drag | boost | KE_cf | r̄_birth | μ̄ | t_exit | n_exit | t_shed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| supp(bare) | 9 | 1.959 | 13.38 | +0.094 | 0.812 | 0.000 | 1.959 | 22.99 | +0.79 | 0.48 | 21 | — |
| n = 1 | 14 | 2.479 | 19.10 | +0.084 | 1.187 | **+0.908** | **0.955** | 20.33 | +0.62 | 0.83 | 21 | 2.00 |
| n = 2–4 | 23 | 1.434 | 14.09 | +0.056 | 1.799 | +0.492 | 0.628 | 16.07 | +0.45 | 1.44 | 20.0 | 1.32 |
| n = 5–9 | 25 | 0.357 | 6.53 | −0.034 | 2.466 | +0.137 | 0.215 | 16.44 | −0.20 | 2.76 | 17.6 | 1.44 |
| n ≥ 10 | 24 | 0.083 | 2.84 | −0.080 | 2.663 | +0.067 | 0.066 | 16.32 | −0.68 | 5.27 | 16.1 | 1.26 |
| retained | 5 | 0.015 | 1.11 | −0.107 | 2.755 | +0.028 | 0.015 | 22.32 | −0.93 | 9.48 | 16.6 | 2.15 |

(`KE_cf` = the co-moving-shed counterfactual, below. c3 — current law,
τ = 6.55 — same structure at lower speed: n₁ KE 0.666, boost +0.220,
KE_cf 0.256.)

### Verdict on the pre-registered hypothesis — refuted as dominant

The §I.11.2 argument ("v = 19.1 > 15.4 Å/ps requires the mass-asymmetric
Coulomb split") does **not** survive measurement:

1. **The measured Coulomb-share excess is +0.08 eV (≈ 3 %), not ~1 eV.**
   n₁ ions ride the acceleration essentially mass-symmetric: they exit the
   droplet at 0.7–0.9 ps with the **full n = 21 shell** (mean masses over
   the first 2 ps: 210 vs 209 amu self/partner) and shed only *after*
   ejection (t̄_shed ≈ 1.7–2.2 ps > t_exit). The hypothesized compounding
   ("early shed ⇒ light while the force acts ⇒ bigger share + less drag")
   does not operate — there is no light-early route in the data.
2. **The dominant term is the cold-shed momentum convention (+0.91 eV).**
   The delivered Tier-2 evaporation channel composes
   `mass_jump.cold_shed_velocity_components` (`physics/evaporation.py`):
   each shed leaves the He **at rest in the lab frame**, the complex keeps
   its full momentum, so `v → v·m/m′` and KE rises by `KE·Δm/m′` per shed
   (verified event-by-event in the stored series: v_ratio ≡ m_ratio at
   every shed; pickup is the exact mirror capture). A full post-exit strip
   21→1 multiplies KE by `m₂₁/m₁ = 1.611`. `E_mass_transfer` books exactly
   the injected energy (−0.908 eV for n₁), which is how the ledger closes.
3. **The co-moving counterfactual reproduces the twin to ≤ 2 % per config.**
   For post-exit sheds momentum conservation gives the exact counterfactual
   `KE_cf = ½·m_det·(v_det·m_det/m_exit)²` (co-moving sheds keep v, not p):
   n₁ KE_cf = **0.955 / 0.952 / 0.256 / 0.933 eV** (c1/c2/c3/c4) vs the
   twin's **0.97 / — / 0.261 / —**. The twin's fate map *is* the co-moving
   convention in disguise — **I57's "histogram parity without composition
   parity" is entirely the shed-frame convention**, not a residual-geometry
   mystery. The composition axis (near-edge outward births, drag deficit
   1.19 vs 2.66 eV) decides *which* ions occupy n₁; the convention decides
   their *speed*.

### Consequences

- **OQ-J fired** (§7): the shed-frame convention is a first-order
  observable-level choice at production kinematics. The A8 "cold-shed"
  grounding ([Nat23], Na⁺) is an **at-rest** result — for a complex at
  rest, "He at ≈ 0 KE" and "He co-moving" coincide; at 10–19 Å/ps they
  diverge by ×1.6 in n₁ KE. Tier-1a adjudicated continuous-velocity as the
  physical path for the anchored channel; the Tier-2 generative channel
  silently composed the cold operator (pickup-capture symmetry). Physical
  evaporation leaves He ~co-moving plus an isotropic thermal recoil —
  i.e. continuous shed + RQ2's ε; OQ-J and RQ2 are one coupled discussion.
- **The experimental n₁ mean (1.302 eV) sits between the two conventions**
  on the capped-tail configs (cold 2.42–2.48; co-moving 0.93–0.96) — the
  bracket is convention-spanned, not geometry-spanned. Leg B's dressed
  re-read (smaller shell at birth ⇒ smaller boost factor `m_birth/m_det`)
  shrinks the gap mechanically.
- **Calibration-transfer warning:** the §4k (v_c, τ) joint closure passed
  the KE bar in the *twin*, i.e. under co-moving bookkeeping. A cold-shed
  MD inherits a ×1.6 n₁ lift the twin never sees — the same (v_c, τ) does
  **not** land both conventions. Any leg-B/T5 KE pre-registration must
  state its shed-convention basis explicitly (the item-3 amendment,
  sharpened).
- **RQ3 discriminator (noted):** the suppressed/bare class sheds nothing,
  so its detected read is convention-free (v ≈ 13.4 Å/ps intact-complex) —
  but the *fragmentation* read is convention-decided: co-moving break-up
  keeps v (bare KE ≈ 1.18 eV at c1 speeds), momentum-conserving-bare gives
  `v·m₂₁/m_I` ≈ 22.3 Å/ps (≈ 3.3 eV). The experimental bare-bin mean KE is
  a direct fragmentation-convention discriminator.

### Boundaries

Zero new MD; N = 50 single seed per config; the counterfactual is exact
only for post-exit sheds (n₁/bare: exact; deep bins: first-order — their
boost is ≤ 0.15 eV anyway); `resid` absorbs well work and share-integral
numerics (≤ 0.06 eV live classes); nothing here adjudicates the
convention — that is OQ-J's user decision; nothing discharges F5.

---

## 4o. T9 leg A″ — the co-moving convention lands the twin's KE axis in real MD; twin parity is now two-axis; the histogram is convention-blind as claimed

Executed 2026-07-17 (leg trigger "[PROCEED TO IMPLEMENTATION] leg A″";
delivery + pre-registration record: log entries "T9 leg A″ TRIGGERED" /
"EXECUTED"). Exactly **one lever** flipped vs the A′ dirs:
`evaporation_shed_convention = "co_moving"` (the OQ-J adjudicated working
convention; the one-lever lock test asserts the field-diff is exactly
that). Dirs: `…_tier2probe_conf270_apcmc{1..4}`, all five artifacts,
N = 50 bridge seed, margin 3 Å, `exclude` policy, 8000 ps cap.

### The A/B/twin read (conventions as §4m; twin = the §4m pre-registration)

| config | supp→bare (A″/A′/twin) | n₁ (A″/A′/twin) | n̄_det (A″/A′/twin) | W₁(A″,A′) | W₁(A″,twin) |
|---|---|---|---|---|---|
| c1 | 0.096/0.095/0.089 | 0.149/0.147/0.172 | 5.89/5.80/6.22 | **0.14** | 0.68 |
| c2 | 0.032/0.032/0.035 | 0.151/0.160/0.170 | 6.14/5.97/6.55 | **0.20** | 0.71 |
| c3 | 0.326/0.315/0.336 | 0.081/0.079/0.069 | 4.08/4.12/4.48 | **0.17** | 0.51 |
| c4 | 0.117/0.116/0.118 | 0.117/0.126/0.138 | 5.84/5.75/6.05 | **0.15** | 0.58 |

| config | n₁ KE: A″ | §4n counterfactual | twin | A′ (cold) | n₁ v̄ [Å/ps] |
|---|---|---|---|---|---|
| c1 | **1.005** | 0.955 | 0.972 | 2.479 | 12.16 |
| c2 | **1.025** | 0.952 | 0.963 | 2.473 | 12.28 |
| c3 | **0.312** | 0.256 | 0.261 | 0.666 | 6.77 |
| c4 | **0.991** | 0.933 | 0.949 | 2.423 | 12.08 |

### Verdicts (pre-registered A″-P1..P4)

1. **A″-P1 CONFIRMED — the histogram is convention-blind.** W₁(A″, A′)
   = 0.14–0.20 bins (vs the 0.55–0.69 twin distance); the suppressed/bare
   ordering c3 ≫ c4 > c1 > c2 transfers exactly; class fractions move by
   ≤ 3 ions. The I59 claim survives real MD.
2. **A″-P2 CONFIRMED — I57 is closed in the real pipeline.** n₁ mean KE
   lands at 1.005/1.025/0.312/0.991 eV vs twin 0.972/0.963/0.261/0.949 —
   the cold ×2.5 divergence collapses to ×1.03–1.20, and the **whole
   solvated per-bin KE curve now matches the twin** (n = 2: 0.74/0.74;
   n = 3: 0.59/0.57; n = 5: 0.36/0.36 on c1; same quality on c2–c4).
   Twin↔MD parity is now **two-axis** (histogram + KE). A small uniform
   residual (+0.04–0.06 eV above the §4n counterfactual, all configs)
   remains — composition-level (in-bubble sheds ride outside the
   counterfactual's post-exit-momentum assumption), not convention-level.
3. **A″-P3 CONFIRMED in signature, strict sub-claim refined.** Every shed
   ion books `E_mass_transfer` > 0; the ion-stage 5-term closure drift is
   unchanged vs the cold arm (max transient ≈ 0.11–0.14 eV, mean end-drift
   −8·10⁻⁴ eV). The "suppressed/retained KE identical to A′" oracle holds
   only to **≤ 10–20 meV**, not bitwise: ion pairs are Coulomb-coupled, so
   a never-shedding ion still feels its opened partner's
   convention-shifted trajectory. Corollary: **1–3 marginal near-barrier
   ions per config flip into the retained class** (c3: 11 → 14) — the
   same near-barrier fragility the A′ execution record flagged, now
   demonstrated at meV-scale perturbations (sharpens the E2-dissipation
   adjudication's stakes for N = 500).
4. **A″-P4 CONFIRMED.** n₁ detected speeds 12.1–12.3 Å/ps (pre-registered
   ≈ 11.9; was 19.1 cold).

Class-composition note (expected, not scored): the bare-candidate bin's
KE stays the intact-complex read (1.93–1.98 eV vs the twin's fragmented
1.16–1.19) — suppressed ions never shed, so their read is
convention-free; that bin's KE is RQ3's file, not OQ-J's.

### Boundaries

Undressed geometry (n_eject = 21, p_couple = 0) throughout — nothing here
reads on the experimental targets; N = 50 single seed, one detection RNG
realization; the twin's KE basis and the MD's are now the same convention
by construction (the item-3 amendment's stated-basis requirement is
satisfied for leg B); nothing discharges F5.

---

## 4p. T9 leg B — the T5 dressing transfers quantitatively: the suppressed/bare class triples along the twin's prediction, the tail truncates at the dressed band, and every deviation is a pre-listed channel

Executed 2026-07-18 (leg trigger "Go ahead with leg B run and compare to
twin"; delivery + pre-registration record: log entries "T9 leg B
TRIGGERED" / "EXECUTED"). Exactly **one lever** flipped vs the certified
apcm (A″) baseline: `initial_shell_model = "density_tied"` (Slice T5; the
one-lever lock test asserts it). Dirs: `…_tier2probe_conf270_bc{1..4}`,
all five artifacts, N = 50 bridge seed, margin 3 Å, `exclude` policy,
8000 ps cap, co-moving shed convention (the A″ basis). Twin: `stage_legb`
re-score at the dressed configuration (m = 20000; dressed chord mass;
`b_undressed` anchor rows equal to the A′ twin **exactly** — the in-stage
wiring oracle).

### The A/B/twin read (conventions as §4m; KE co-moving on both sides)

| config | supp→bare (B / A″ / twin_b) | n̄_det (B / A″ / twin_b) | n₁ (B / A″ / twin_b) | W₁(B, twin_b) | W₁(B, A″) |
|---|---|---|---|---|---|
| c1 | **0.322** / 0.096 / 0.389 | 4.13 / 5.89 / 3.94 | 0.100 / 0.149 / 0.101 | **0.49** | 1.76 |
| c2 | **0.267** / 0.032 / 0.349 | 4.40 / 6.14 / 4.17 | 0.122 / 0.151 / 0.105 | **0.51** | 1.74 |
| c3 | **0.464** / 0.326 / 0.533 | 3.08 / 4.08 / 2.78 | 0.071 / 0.081 / 0.049 | **0.41** | 1.00 |
| c4 | **0.378** / 0.117 / 0.436 | 3.97 / 5.84 / 3.67 | 0.067 / 0.117 / 0.076 | **0.46** | 1.87 |

Seed verification: MD dressed n₀ q05 = 13.9 / mean = 17.18 per config
(twin n_eject q05 = 13.0 / mean = 16.71; min single-atom n₀ = 12 — the
per-atom ±R0/2 offset reaching just below the twin's molecule-center
margin floor). Trapped: B 10/10/16/10 per 100 vs twin
0.055/0.060/0.093/0.055 and A″ 6/7/14/6.

| config | n₁ KE: B | twin_b | A″ (undressed) | bare-class KE: B (intact complex) | B rescaled to bare mass | twin_b (bare) |
|---|---|---|---|---|---|---|
| c1 | **0.661** | 0.549 | 1.005 | 1.598 | **1.084** | 0.999 |
| c2 | **0.626** | 0.512 | 1.025 | 1.461 | **0.999** | 0.891 |
| c3 | **0.263** | 0.176 | 0.312 | 0.585 | **0.392** | 0.352 |
| c4 | **0.601** | 0.499 | 0.991 | 1.516 | **1.023** | 0.945 |

### Verdicts (pre-registered BP-P1..P4)

1. **BP-P1 CONFIRMED — with the pre-registered caveat firing in the
   stated direction.** The suppressed/bare class roughly triples
   (c1 0.096 → 0.322; c2 0.032 → 0.267; c4 0.117 → 0.378) and the
   ordering c3 > c4 > c1 > c2 transfers exactly. MD lands systematically
   **0.06–0.08 below the twin** on every config — the listed
   pickup-re-filling channel (under-dressed shells re-fill in-bubble, Σ
   grows before gate-open, fewer ions stay self-unbound), now measured
   at ≈ 7 ions/100.
2. **BP-P2 CONFIRMED.** n̄_det lands 4.13/4.40/3.08/3.97 (pre-reg ≈
   3.9/4.2/2.8/3.7; MD +0.2–0.3 above twin — the same re-filling
   direction); the deep tail truncates near the dressed band (top
   detected weight ends by n ≈ 12–14 vs A″ weight to n ≈ 18+);
   **W₁(B, twin_b) = 0.41–0.51 bins — the best twin↔MD histogram
   agreement of the whole oracle chain** (A′/A″: 0.51–0.71), while the
   lever's own move is 2–4× larger (W₁(B, A″) = 1.0–1.9).
3. **BP-P3 CONFIRMED on the solvated curve; the bare bin resolves into
   RQ3 bookkeeping.** n₁ mean KE drops from the A″
   1.005/1.025/0.312/0.991 to **0.661/0.626/0.263/0.601** — onto the
   dressed twin's 0.549/0.512/0.176/0.499 at ×1.20/×1.22/×1.50/×1.20
   (c3's ×1.50 is a 6-ion bin; the mid-bins n5–n12 straddle the twin
   within ~±30 % with no systematic sign). The apparent ×1.6 bare-bin
   excess (1.598 vs 0.999 at c1) is **not dynamics**: the MD reports the
   intact dressed complex (m̄ ≈ 187–191 amu) while the twin books the
   RQ3 spec-b bare fragment at m_I; rescaling the MD read by
   m_I/m_complex (the co-moving break-up value) collapses it to
   **×1.08–1.12** on every config — exactly the §4n "fragmentation read
   is convention-decided" note, pre-listed, now with numbers.
4. **BP-P4 CONFIRMED.** Trapped 10/10/16/10 per 100 = ×1.7–1.8 the
   twin's 150 ps chord read (the A′ precedent was ×1.6), c3 largest on
   both sides; the dressing raises the trapped class vs A″ (6/7/14/6)
   in the twin-predicted direction.

**S2c-P4 holds:** every observed deviation is one of the three listed
channels (pickup re-filling; per-atom vs molecule-center dressing —
measured at ≈ +0.5 He mean and a single n₀ = 12 atom; trapped-class
dynamics beyond the chord read). No unlisted divergence appeared.

### Boundaries

N = 50 single seed, one detection RNG realization per dir (class
fractions quantize at ~1 %); the dressed read is at the C-matrix knob
values — the Step-1c closure basin was located under the 1D ensemble
with droplet priors, so **nothing here scores the experimental
targets** (that is the re-centered T9 re-pilot's job after T6/T8);
the bare-bin KE remains RQ3/RQ8-gated (the rescale above is a
bookkeeping identity, not a fragmentation model); T6 (p-law) not
flipped — E_int(0) is constant across the dressed ensemble; nothing
discharges F5.

---

## 4q. T9 leg C — the T6 `sigma_proportional` p-law lands on the twin: p = 1 de-suppresses the over-suppressed side to the twin's prediction (W₁ 0.31–0.47, the chain's best), n̄ and the solvated KE curve transfer, and every deviation is a pre-listed channel

Executed 2026-07-19 (leg trigger "Go ahead with the leg-C
pre-registration" → "can we also make the actual C runs?"; delivery +
pre-registration record: log entries "T9 leg C TRIGGERED" / "EXECUTED").
Exactly **one lever** flipped vs the certified `bc` (leg-B) baseline:
`internal_energy_partition_law = "sigma_proportional"` (Slice T6, p = 1;
the one-lever lock test asserts it — cumulative on leg B's `density_tied`
+ co-moving). Dirs: `…_tier2probe_conf270_cc{1..4}`, all five artifacts,
N = 50 bridge seed, margin 3 Å, `exclude` policy, 8000 ps cap, co-moving
shed. Twin: `stage_legc` re-score at p = 1 (m = 20000; the `c_p0` anchor
rows equal the leg-B dressed twin **exactly** — the in-stage wiring oracle,
verified at m = 20000).

**Execution-provenance note (session-specific, physics-neutral).** The
`cc` `relaxation.npz` are ~2 MB (vs the `bc` dirs' ~380 MB): the E2
checkpoint byte-budget was shrunk for this session's execution window (the
full-trajectory compressed save otherwise exceeded the process lifetime).
This coarsens only the *intermediate* relaxation trajectory; detection
seeds from the **terminal column only** (`detection_stage.py`
`ion_state_from_checkpoint_column(seed_ckpt, -1)`) and the stage always
stores the true final state last (`relaxation_stage.py` 439–440), so every
scored read below (all terminal/detected) is **bit-exact** — the coarser
stride touches no CP-P1..P4 quantity. Per-stage resume via a scratchpad
driver over the delivered stage functions (zero repo-code change).

### The A/B/twin read (conventions as §4m; KE co-moving on both sides)

| config | supp→bare (C / twin_c / B) | n̄_det (C / twin_c / B) | n₁ (C / twin_c / B) | W₁(C, twin_c) | trapped (C / twin) |
|---|---|---|---|---|---|
| c1 | **0.101** / 0.111 / 0.322 | 4.87 / 4.98 / 4.13 | 0.191 / 0.212 / 0.100 | **0.44** | 0.110 / 0.055 |
| c2 | **0.022** / 0.051 / 0.267 | 5.27 / 5.31 / 4.40 | 0.202 / 0.217 / 0.122 | **0.47** | 0.110 / 0.060 |
| c3 | **0.301** / 0.343 / 0.464 | 3.76 / 3.77 / 3.08 | 0.108 / 0.084 / 0.071 | **0.31** | 0.170 / 0.093 |
| c4 | **0.124** / 0.143 / 0.378 | 4.87 / 4.88 / 3.97 | 0.146 / 0.172 / 0.067 | **0.45** | 0.110 / 0.055 |

| config | n₁ KE: C | twin_c | B (bc) | A″ undressed | bare-class KE: C (intact) |
|---|---|---|---|---|---|
| c1 | **1.077** | 0.998 | 0.661 | 1.005 | 1.990 |
| c2 | **1.104** | 0.964 | 0.626 | 1.025 | 1.953 |
| c3 | **0.317** | 0.253 | 0.263 | 0.312 | 0.694 |
| c4 | **1.089** | 0.973 | 0.601 | 0.991 | 1.939 |

The solvated per-bin KE curve tracks the twin (c1 n2 0.74/0.67, n3 0.53/0.47,
n4 0.37/0.36, n5 0.28/0.29; same quality c2–c4).

### Verdicts (pre-registered CP-P1..P4)

1. **CP-P1 CONFIRMED — the de-suppression headline.** p = 1 drops the
   suppressed/bare class to roughly a third of leg B (c1 0.322 → 0.101;
   c2 0.267 → 0.022; c3 0.464 → 0.301; c4 0.378 → 0.124) and lands **on the
   twin p = 1** (0.111/0.051/0.343/0.143), ordering c3 > c4 > c1 > c2 exact.
   The pre-registered pickup-re-filling channel fires in the **corrected**
   direction: MD lands **≈ 1–3 ions/100 below** the twin ("fewer stay
   self-unbound" — under-dressed shells re-fill, Σ grows before gate-open;
   the TRIGGERED entry's word "above" was a slip, the parenthetical
   mechanism was right, and leg B measured the same MD-below-twin sign at a
   larger ≈ 7/100 gap — the gap shrinks at p = 1 because there is less
   suppression to re-fill).
2. **CP-P2 CONFIRMED.** n̄_det rises to 4.87/5.27/3.76/4.87 (from `bc`
   4.13/4.40/3.08/3.97), matching the twin to ≤ 0.11 He; n₁ ≈ doubles vs
   `bc`; the deep tail truncates near the dressed band (top weight to
   n ≈ 12–14). **W₁(C, twin_c) = 0.31–0.47 bins — the best twin↔MD histogram
   agreement of the whole oracle chain** (A′ 0.55–0.69, A″ 0.51–0.68, B
   0.41–0.51), while the lever's own move (`bc` → `cc`) is far larger.
3. **CP-P3 CONFIRMED.** n₁ mean KE **rises** to 1.077/1.104/0.317/1.089 eV
   (from `bc` 0.661/0.626/0.263/0.601) — the de-suppressed n₁ occupants ride
   shallower descents. It lands ≈ the A″ undressed (1.005/1.025/0.312/0.991)
   and **closest to the committed experimental n₁ = 1.302 eV of any leg** on
   the solvated branch; MD sits ×1.08–1.15 above twin (the same small uniform
   composition residual as A″). Bare-class KE (1.99/1.95/0.69/1.94, intact
   dressed complex) stays RQ3/RQ8 bookkeeping — renormalised out.
4. **CP-P4 CONFIRMED.** Trapped 0.110/0.110/0.170/0.110 (≈ `bc`
   0.10/0.10/0.16/0.10 — p-invariant, as pre-registered: the p-law never
   touches the chord dynamics), c3 largest, ≈ 1.8× the twin's 150 ps chord
   read (the A′/leg-B precedent).

**S2c-P4 holds:** every deviation is one of the three listed channels
(pickup re-filling, now measured ≈ 1–3 ions/100; per-atom vs
molecule-center dressing; trapped-class dynamics). No unlisted divergence.

### Boundaries

N = 50 single seed, one detection RNG realization per dir (class fractions
quantize at ~1 %); the shrunk-checkpoint provenance note above; the dressed
read is at the C-matrix knob values — the Step-1c basin was located under
the 1D ensemble, so **nothing here scores the experimental targets** (the
re-centered T9 re-pilot's job after T8); the bare-bin KE stays RQ3/RQ8-gated;
the E2 Landau-gated drag arm (§I.11.2 item 2) is still unbuilt — required
before N = 500, not the N = 50 legs; nothing discharges F5.

## 4r. T9 leg D — the T8 kornilov droplet prior: the axis transfers directionally with the suppressed ordering exact and the KE axis prior-quiet as the twin claimed; but the frozen-geometry twin under-models the axis — W₁ regresses to 0.54–0.73 and the histogram softens ≈ 0.5 He below the twin (channel (d)'s first measured bite)

Executed 2026-07-20 (log entries "T9 leg D TRIGGERED" / "EXECUTED"). One
semantic lever vs the certified `cc` baseline:
`droplet_size_prior = "kornilov_lognormal"` at the D4 primary (⟨N⟩ = 2000,
δ = 0.625; the guard-paired `use_single_droplet_size=False` is the same
lever) — legs A″/B/C carry over cumulatively (co-moving shed, `density_tied`,
p = 1). N = 50 bridge seed, margin 3 Å, exclude policy, 8000 ps cap,
`zero_gamma` (T8-D4: the Landau arm stays the N = 500 gate). The four `dc`
dirs ran under the I67 shrunk-checkpoint accommodation (per-config parallel
resume drivers; ~2 MB `relaxation.npz`; ≈ 14 min wall-clock for all four).
The scorer re-implementation was oracle-locked first: it reproduces the
recorded §4q `cc` numbers exactly on every column.

### The A/B/twin read (conventions as §4m; KE co-moving on both sides)

Pre-registered twin: `h2b_leg_d_predictions.csv` / `h2b_leg_d_ke.csv`
(committed; the `d_delta` wiring oracle reproduces the leg-C twin verbatim
at m = 20000). Realized prior N_q05/q50/q95 = 596/1638/4636.

| config | supp (D / twin_d / C) | n̄_det (D / twin / C) | n₁ (D / twin) | W₁(D, twin) | n₁ KE (D / twin / C) | trapped (D / twin) |
|---|---|---|---|---|---|---|
| c1 | 0.087 / 0.127 / 0.101 | 3.89 / 4.42 / 4.87 | 0.228 / 0.236 | **0.71** | 1.037 / 0.992 / 1.077 | 0.080 / 0.051 |
| c2 | 0.077 / 0.058 / 0.022 | 4.19 / 4.77 / 5.27 | 0.154 / 0.246 | **0.73** | 0.982 / 0.958 / 1.104 | 0.090 / 0.056 |
| c3 | 0.393 / 0.383 / 0.301 | 2.82 / 3.35 / 3.76 | 0.067 / 0.087 | **0.54** | 0.317 / 0.245 / 0.317 | 0.110 / 0.081 |
| c4 | 0.098 / 0.163 / 0.124 | 3.93 / 4.36 / 4.87 | 0.196 / 0.191 | **0.68** | 1.038 / 0.967 / 1.089 | 0.080 / 0.051 |

### Verdicts (pre-registered DP-P1..P4)

1. **DP-P1 SPLIT — ordering + twin band CONFIRMED, the net-rise-vs-`cc`
   arm REFUTED on c1/c4 through the listed channel (a).** The suppressed
   ordering transfers exactly (0.393 > 0.098 > 0.087 > 0.077 = c3 > c4 >
   c1 > c2). c2/c3 rise as predicted (0.022 → 0.077, 0.301 → 0.393, both
   within ~2 ions/100 of the twin — statistical at ~90 fragments). But
   c1/c4 *fall* vs `cc` (0.101 → 0.087, 0.124 → 0.098), landing **below**
   the twin by 4.0/6.5 ions/100 — inside the pre-listed re-filling band
   (1–7 ions/100): where the small-droplet suppression push is modest,
   in-bubble pickup re-filling overcompensates it.
2. **DP-P2 SPLIT — direction CONFIRMED, magnitude ≈ 2× the prediction, the
   W₁-scale claim MISSED.** n̄_det falls on every config, but by 0.94–1.08 He
   (predicted 0.42–0.55) — MD sits a **uniform 0.43–0.58 He below the twin**.
   n₁ lands in the predicted band on c1/c4 (0.228/0.196 vs 0.236/0.191), low
   on c2 (0.154 vs 0.246) and c3 (0.067 vs 0.087). **W₁(D, twin) = 0.54–0.73
   — the chain-best trend breaks** (A′ 0.55–0.69 → A″ 0.51–0.68 → B
   0.41–0.51 → C 0.31–0.47 → D 0.54–0.73). Statistical context: the
   per-config n̄ deficit is ≈ 1 SE (σ_n ≈ 4.5 over ~90 detected fragments →
   SE(n̄) ≈ 0.47), but it is same-signed on all four configs — a systematic
   softening, attributed to the *listed* channel (d): the twin's chord is
   frozen at birth geometry while the MD's well depth, gate depth, and
   pickup exposure follow each ion's R_i through the live cascade. The
   droplet axis is exactly where that freeze bites first.
3. **DP-P3 CONFIRMED — the structural claim of the leg.** The KE axis is
   prior-quiet in MD exactly as the twin claimed: n₁ mean KE moves ≤ 0.12 eV
   under the prior flip (1.077 → 1.037 / 1.104 → 0.982 / 0.317 → 0.317 at
   the printed precision / 1.089 → 1.038), and MD tracks the twin at
   ×1.03–1.07 (c3 ×1.29 — a 0.07 eV absolute gap in the low-KE regime).
   The solvated KE curve stays owned by the drag law + onset, not the
   droplet axis; distance to the experimental n₁ = 1.302 eV stays ≈ leg C's.
4. **DP-P4 CONFIRMED (ratio refined).** Trapped drops under the prior
   exactly as the twin said (0.110/0.110/0.170/0.110 → 0.080/0.090/0.110/
   0.080), c3 largest, at ×1.4–1.6 the twin's 150 ps chord read (predicted
   1.7–1.9): small-droplet marginal ions clear the shallower well even a
   little more readily in real MD.

**S2c-P4 holds:** every deviation lands in a listed channel — (a)
re-filling on the suppression axis, (d) frozen-vs-dynamic geometry on the
histogram axis, N = 50 statistics elsewhere. No unlisted divergence.

### What leg D means for the endgame

The §I.11 oracle chain A′→D is **complete**: every physics arm (position,
shed convention, dressing, p-law, droplet prior) transfers directionally to
real MD. The droplet axis is the first where the 1D twin's *quantitative*
authority ends — the W₁ regression and the uniform softening turn the
post-leg-B stance ("the twin is a basin-locator, not a predictor") from a
caveat into a measurement. Consequence for the re-pilot: the
confirmation-matrix re-centering must be MD-driven (as already decided),
with the twin's droplet-axis output used for ordering/direction only; the
N = 500 winner run (Landau arm on) is the statistical arbiter of the
0.5 He softening.

### Boundaries

N = 50 single seed (~90 detected fragments per config; class fractions
quantize at ~1 %, SE(n̄) ≈ 0.47 — the softening is a four-config joint
read, not per-config significant); the `dc` dirs use the I67
shrunk-checkpoint accommodation; the knob values are the C matrix —
**nothing here scores the experimental targets** (the re-centered re-pilot's
job); the bare bin stays RQ8-gated; nothing discharges F5.

---

## 4s. Re-pilot Stage 1 — the v_c bracket in real MD: the KE axis pins v_c only one-sidedly (χ² argmin right-censored at the bracket top, twin+1-step on every arm), the speed-selective tension survives the real cascade, and channel (d) is measured to scale with drag exposure (Δn̄ 0 → −0.8 He across the bracket)

Executed 2026-07-20 (log entries "Re-pilot Stage-1 pre-registration
EXECUTED" / "Stage-1 MD sweep EXECUTED"). The §I.11.4.2 Stage-1 KE-pin
sweep: 15 cells = the frozen v_c brackets (c1/c4: 6.5–8.5, c2: 5.5–7.5,
step 0.5) at the full leg-D configuration, one knob (v_c) moving. 13 new
N = 50 dirs (`s1c1v65`-style conf tags) + the certified `dc2`/`dc4`
reused as bracket centers; `s1c1v75` re-run fresh as the
**pipeline-identity oracle** — it reproduces the certified `dc1` row
exactly on every column (full-stack determinism, not just scorer
identity). Scored by the delivered Stage-0 scorer (`min_count = 2`,
sim-SE widened, committed error model).

**Execution wrinkle (records the I70 pattern amendment):** the four
4-cell worker processes were killed by the session process window
mid-sweep (8 cells complete, 4 partial-through-ion, 1 unstarted); all 13
cells completed via **one-cell-per-process** per-stage resume drivers
(the leg-C I67 route: E2+detection resumed from the intact `ion.npz`;
partial artifacts verified loadable before resuming). One cell per
process is the pattern going forward.

### The A/B/twin read (conventions §4r; χ² = profiled committed error model)

| cell | v_c | trap | supp (MD/twin) | n̄ (MD/twin) | Δn̄ | n₁ | W₁(MD,twin) | n₁KE (MD/twin) | χ²_prof | W₁_solv | n₁/n₂ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| s1c1v65 | 6.5 | 0.010 | 0.091/0.140 | 3.48/3.53 | −0.04 | 0.303 | **0.209** | 1.240/1.229 | 102.9 | 1.055 | **2.14** |
| s1c1v70 | 7.0 | 0.050 | 0.095/0.134 | 3.60/3.96 | −0.36 | 0.295 | 0.438 | 1.113/1.114 | 95.8 | 0.912 | 2.15 |
| s1c1v75 | 7.5 | 0.080 | 0.087/0.127 | 3.89/4.42 | −0.53 | 0.228 | 0.709 | 1.037/0.992 | 54.4 | 0.771 | 1.31 |
| s1c1v80 | 8.0 | 0.100 | 0.089/0.119 | 4.22/4.88 | −0.65 | 0.178 | 0.862 | 0.927/0.863 | 43.1 | **0.720** | 0.94 |
| s1c1v85 | 8.5 | 0.110 | 0.079/0.110 | 4.48/5.30 | −0.82 | 0.157 | 1.038 | 0.830/0.731 | **27.6** | 0.811 | 0.70 |
| s1c2v55 | 5.5 | 0.050 | 0.074/0.066 | 3.75/4.11 | −0.36 | 0.189 | 0.491 | 1.134/1.129 | 49.9 | 1.096 | 0.90 |
| s1c2v60 | 6.0 | 0.050 | 0.074/0.062 | 4.09/4.44 | −0.35 | 0.179 | 0.483 | 1.042/1.042 | 93.0 | 0.789 | 0.89 |
| s1c2v65 | 6.5 | 0.090 | 0.077/0.058 | 4.19/4.77 | −0.58 | 0.154 | 0.730 | 0.982/0.958 | 24.1 | 0.801 | 0.70 |
| s1c2v70 | 7.0 | 0.100 | 0.067/0.054 | 4.31/5.09 | −0.78 | 0.156 | 0.916 | 0.917/0.877 | 43.7 | 0.839 | 0.67 |
| s1c2v75 | 7.5 | 0.110 | 0.056/0.049 | 4.63/5.40 | −0.77 | 0.146 | 0.912 | 0.885/0.800 | **23.6** | 0.878 | 0.68 |
| s1c4v65 | 6.5 | 0.010 | 0.101/0.178 | 3.59/3.53 | +0.06 | 0.253 | **0.258** | 1.250/1.210 | 159.7 | 1.032 | **2.08** |
| s1c4v70 | 7.0 | 0.050 | 0.095/0.172 | 3.75/3.92 | −0.17 | 0.221 | 0.457 | 1.161/1.093 | 171.1 | 0.949 | 1.40 |
| s1c4v75 | 7.5 | 0.080 | 0.098/0.163 | 3.93/4.36 | −0.42 | 0.196 | 0.679 | 1.038/0.967 | 115.9 | 0.877 | 1.20 |
| s1c4v80 | 8.0 | 0.100 | 0.100/0.153 | 4.27/4.78 | −0.52 | 0.167 | 0.753 | 0.913/0.835 | 49.2 | 0.879 | 0.88 |
| s1c4v85 | 8.5 | 0.110 | 0.101/0.140 | 4.57/5.20 | −0.62 | 0.135 | 0.839 | 0.797/0.702 | **22.9** | 1.030 | 0.71 |

(Experimental references: solvated n₁/n₂ = 0.310/0.142, ratio 2.18;
n₁ KE anchor 1.302 eV.)

### Verdicts (pre-registered S1-P1..P5)

1. **S1-P1 CONFIRMED twice over.** The twin-side oracle passed at
   generation (bracket centers ≡ leg-D `d` rows, m = 20000), and the MD
   side adds full pipeline identity: the fresh `s1c1v75` regeneration
   reproduces the certified `dc1` scored row exactly on every column.
2. **S1-P2 CONFIRMED with one edge break.** n₁-KE is monotone
   (v_c ↑ ⇒ KE ↓) on all three arms; the c1/c4 KE degeneracy holds
   (split 0.001–0.048 eV ≤ the registered 0.06). The ×1.10 n₁ band
   holds on 12/15 cells and breaks **only at the bracket-top v85 cells**
   (×1.13–1.14): the MD-above-twin bias *grows with v_c* — at the low
   edge it vanishes (×1.01). Deep-bin (n ≥ 5) band compliance is
   partial at N = 50 (2–6 of 6–9 bins) — thin-bin statistics, reported.
3. **S1-P3 CONFIRMED — with a right-censoring caveat.** The MD
   profiled-χ² argmin lands at the twin's cell **plus exactly one step
   toward higher v_c on every arm** (c1: v85, c2: v75, c4: v85 — all in
   the registered sets). Caveat: the c1/c4 (and c2) argmins sit at the
   bracket **edge with χ² still falling** — the whole-curve v_c\* is
   right-censored; the bracket does not close the optimum from above.
4. **S1-P4 CONFIRMED — the form discrimination is now two-sided.** The
   p = 0 arm cannot reach the 1.302 eV n₁ anchor in-bracket (MD cap
   1.134 eV at v55, on the twin's 1.129 — the registered claim, exact).
   The p = −1 arms reach 1.240/1.250 eV at v65 — *short* of the
   registered 1.27–1.35 window because the ×1.03–1.07 bias is absent at
   low v_c; reaching the anchor extrapolates to **v_c ≈ 6.2–6.3**,
   below the frozen bracket.
5. **S1-P5 — the registered non-uniformity check FIRED (the structure
   finding of the sweep).** Directions confirmed (n̄ ↑ and trapped ↑
   with v_c, monotone on every arm; suppressed falls end-to-end but is
   non-monotone within the bracket — 1–2-ion N = 50 quantization, the
   listed statistical channel). But **Δn̄(MD − twin) is NOT uniform: it
   runs monotonically from ≈ 0 at v65 to −0.62…−0.82 He at v85 on every
   arm.** The leg-D "uniform ≈ 0.5 He softening" was a single-v_c
   snapshot: **channel (d) scales with drag exposure** (more drag →
   longer in-droplet residence → larger frozen-vs-live geometry
   divergence). Twin histogram parity tracks it: W₁(MD, twin) 0.21–0.26
   at v65 (beating the leg-C chain-best) → 0.84–1.04 at v85.

### The Stage-1 physics headline — the speed-selective tension survives real MD

The registered S1-P4 question ("where does MD resolve the H.2b
speed-selective residual?") is answered: **it does not.** The two
experimental axes pull v_c in opposite directions within the
capped-cubic family:

- the **whole-curve KE χ²** improves monotonically toward high v_c
  (argmin right-censored at v85: more drag fits the deep bins);
- the **n₁ anchor** (1.302 eV) and the **solvated ratio** (n₁/n₂ = 2.18)
  land at/below the low edge — v65 gives ratio 2.14/2.08 (nearly exact)
  and n₁ = 0.303/0.253 (the H.2b n₁ ≈ 0.31 target, first time in MD),
  while the χ² there is 4–7× the arm minimum (mid-bins overshoot).

No single v_c lands both — a **drag-form-shape statement** (the taper
family, not the calibration) measured in real MD, echoing I45/H.2b
outcome (c). Note also the W₁_solv column: its per-arm minimum
(v75–v80) coincides with *neither* axis's preference — three
observables, three different v_c optima.

### Boundaries

N = 50 single seed (SE(n̄) ≈ 0.47, class fractions ~1 %; the Δn̄ trend
is a 15-cell joint read); the (τ, E₀) knobs sit at the C-matrix values —
the histogram-side numbers here are *not* Stage-2 reads (τ/E₀ are the
histogram's own levers and have not moved yet); ledger closure not
re-checked per cell (the delivered pipeline's standing residual);
`zero_gamma` E2 throughout (RP-D7); nothing here discharges F5.

---

## 4t. Re-pilot Stage 2 — the (τ, E₀) grid at two carried v_c: the twin is quantitatively exact at v65 across the whole budget plane, the KE and histogram axes each land the experimental target at different cells, and the joint landing does not exist inside the swept capped-cubic family

Executed 2026-07-20/21 (log entries "Stage-2 … pre-registration EXECUTED" /
"Stage-2 MD sweep EXECUTED"). Adjudication (b): 62 new N = 50 cells +
5 Stage-1/`dc2` center reuses = the frozen 67-cell grid (c1/c4 ×
v_c ∈ {6.5, 8.5} × τ ∈ {3.2, 3.8, 4.4, 5.0} × E₀ ∈ {0.21, 0.23, 0.25,
0.27} + the c2 spot diagonal). Execution: one-cell-per-process heal
drivers (run/resume/regenerate per artifact state — several session
kill-sweeps were absorbed at zero rework beyond the in-flight stage).
Scored by the Stage-0 scorer (conventions §4r; χ² = profiled committed
error model, min 2 ions/bin, sim-SE widened).

### Verdicts (pre-registered S2s-P1..P5)

1. **S2s-P1 CONFIRMED** (twin generation oracle; recorded at
   pre-registration).
2. **S2s-P2 CONFIRMED.** n̄_det monotone decreasing in E₀ at fixed τ and
   in τ at fixed E₀, and suppressed monotone increasing in E₀, on **all
   four** (arm, v_c) columns — no exceptions in 32 column-checks.
3. **S2s-P3 CONFIRMED emphatically — the twin is quantitative at v65.**
   Δn̄(MD − twin) at the v65 cells: **32/32 within ±0.25 He** (mean
   −0.036, range [−0.17, +0.08]) with *no* systematic growth along
   either budget axis: channel (d) has no budget-dependence at low
   drag. Combined with S2s-P4 (**31/32** v85 cells inside the carried
   [−0.9, −0.5] He band, mean −0.715, one marginal at −0.91): the
   frozen-geometry bias is a **pure drag-exposure function** across the
   entire (τ, E₀) plane — I73 extended from a line to the full plane.
4. **S2s-P5 SPLIT.** The KE non-factorization is **real in MD** (c1v65
   n₁ KE runs 1.73 → 1.20 eV across the grid — fate-map composition,
   exactly as the twin predicted); the per-cell twin band holds at v65
   (**34/35** within ×1.10) but at v85 only **21/32** sit within ×1.15
   (breaks ×1.14–1.16 — the drag-end bias, again). The Stage-1 KE pin
   demonstrably does not transfer: the §I.11.4.2 guard fired as
   designed.

### The joint-landing read (the leg's real question)

Each experimental axis is **individually landed** — at different cells:

| cell | τ | E₀ | n₁ KE [eV] | χ²_prof | W₁_solv | n₁_solv | n₁/n₂ |
|---|---|---|---|---|---|---|---|
| **s2c1v85t32e23** (KE-best) | 3.2 | 0.23 | **1.310** (anchor 1.302) | **13.9** | 1.520 | 0.101 | 1.00 |
| s2c1v85t32e25 | 3.2 | 0.25 | 1.009 | **13.0** | 1.292 | 0.098 | 0.53 |
| **s2c1v85t38e27** (W₁-best) | 3.8 | 0.27 | 0.518 | 62.9 | **0.594** | **0.320** | 3.43 |
| **s2c1v65t32e27** (joint compromise) | 3.2 | 0.27 | 1.197 | 162.4 | 0.736 | **0.326** | 2.80 |
| s2c4v65t32e27 (best c4) | 3.2 | 0.27 | 1.099 | 227.1 | 0.793 | 0.229 | 2.00 |

(References: n₁ KE 1.302 eV; solvated n₁ = 0.310, n₁/n₂ = 2.18. The
registered v65 × τ3.2 × E₀{0.25, 0.27} candidate region did exactly what
the twin said: the n₁ anchor is *crossed* inside it — 1.355 at e25,
1.197 at e27 — and e27 lands n₁_solv 0.326 ≈ the H.2b 0.31 target, the
first MD cells ever to reach it. The Stage-1 χ² right-censoring is
**resolved interior** by the budget axes: the KE optimum sits at
(v85, τ3.2, E₀ 0.23–0.25), not beyond the v_c bracket — the deferred
extension (a) is likely moot.)

**But no cell lands both.** KE-best cells miss the histogram by ×2–3
(n₁ 0.10 vs 0.31); histogram-best cells miss the KE axis (χ² ≥ 63; the
W₁-best cell's n₁ KE is 0.52 eV — 2.5× low). The I72 drag-form-shape
statement is now **cell-resolved across a 3-knob family**: within
capped_cubic, (v_c, τ, E₀) cannot land the KE curve and the solvated
histogram simultaneously anywhere in the swept volume. The miss
direction is stable: cells that strip enough to populate n₁ ≈ 0.31
arrive there too slow (or, at v85, with the whole curve dragged down);
cells that keep the KE curve keep too much shell.

**Ladder discrimination (I51/S2c-P3 direction, N = 50-soft):** rq4graded
(c1) beats floor1 (c4) on the histogram axis at every matched knob
(top-4 W₁_solv are all c1; c4's best is 0.793); formal verdict stays the
N = 500 finalists' job (RP-D6).

### NB (2026-07-21, post-§4t discussion) — the KE-miss anatomy at the landed cells: the two histogram-landing cells fail the KE axis in opposite ways, and the v65 miss is a mid-bin bulge, not an n₁ miss

Per-bin read of the three key cells (sim/ref mean-KE ratio per bin;
counts 2–28 ions/bin at N = 50):

- **(v65, τ3.2, E₀0.27) — histogram landed, low drag.** n₁ is ×0.92
  (1.197 vs 1.302 eV — within the correlated bands' reach) and the deep
  tail n ≥ 9 sits **on** the reference (×0.93–0.99). The whole χ² = 162
  is a **mid-bin bulge**: n₂–n₈ run ×1.40–1.90 (n₃: 0.855 vs 0.496 eV).
  The experimental curve falls much more steeply from n₁ to n₈ than the
  model's.
- **(v85, τ3.8, E₀0.27) — histogram landed, heavy drag.** A coherent
  whole-curve depression: ×0.40–0.83 on every bin, worst at n₁ (0.518).
- **(v85, τ3.2, E₀0.23) — KE landed.** ×0.84–1.23 on all 13 scored
  bins — the curve *shape* is expressible; this cell just retains shell.

**Physics reading.** At weak drag the model's fragment speed barely
differentiates by retention class — moderately-stripped fragments keep
nearly the n₁-class speed — while the experiment couples retention to
slowness strongly. The data asks for a **bin-selective** correction:
leave the fastest fragments (n₁, already right at v65) and the slowest
(deep tail, already right) untouched, and slow specifically the
*intermediate-speed* band that populates n₂–n₈.

**E-2 design consequence.** A more negative `p_tail` alone acts only
*above* v_c — i.e. on the n₁-class fragments that need no fixing at
v65. The promising joint move is **v_c raised toward 7.0–7.5** (pulling
the mid-speed band back under full cubic drag → damping n₂–n₈) **with
p ∈ {−2, −3}** (keeping the above-cap n₁ population fast) — sharper
speed selectivity from both sides of the threshold; the locked
pure-cubic `b` below the cap stays untouched throughout. Whether a
joint cell exists in that (v_c, p_tail, τ, E₀) volume is the twin-first
E-2 scan's question (I74 licenses the twin as the quantitative search
tool at low drag). Caveat: single mid-bin ratios are soft (5–11
ions/bin), but a same-signed ×1.4–1.9 excess across seven consecutive
bins is not noise.

### Boundaries

N = 50 single seed (n₁/n₂ ratio reads carry SE ≈ ±0.5–0.7 — n₁_solv is
the stable histogram read); the c2 spot leg confirms directions only
(3 cells); `zero_gamma` E2 (RP-D7); ledger closure not re-checked per
cell; nothing here discharges F5.

---

## 4u. The n=1-deweight KE re-score + the E-2 taper-corner twin scan: the n₁ anchor is a mixture-mean artifact (mode 0.89 / median 1.13 / mean 1.30, bg-structural), but de-weighting it does NOT move the high-drag pull — the mid-bins n2–n8 own it; and the E-2 (v_c, p_tail) corner has no joint cell because p_tail re-couples the two axes, while raising v_c alone lands the whole mid-bin KE curve at the cost of histogram over-retention (KE↔histogram anti-correlation through K, re-confirmed)

Executed 2026-07-21 (zero-MD, pure post-processing over the on-disk
Stage-1/Stage-2 scored cells + the committed leg-D twin primitives;
`tier2_confirmation_score` and `tier2_h2b_forward_model` imported
verbatim, nothing in the repo mutated). Motivated by the user's n₁
mixture-mean theory (the two-channel bare reading, RQ8, extended into the
n = 1 bin).

### Part 1 — the n=1-deweight re-score (the ihe_ked reference is self-flagging)

The committed `IHe_KED_reference.csv` fingerprints the n = 1 bin exactly
as the theory predicts: **mode 0.891, median 1.128, mean 1.302 eV**
(mean/mode = 1.46; σ = 0.697), and n = 1 is the **only** low-n fragment
whose `dominantError = bg-structural`, with `bgOffShift = 0.378 eV`
(29 % of the mean — COLUMNS.md: "the bg choice matters essentially only
at n = 1"). The mean is a two-population mixture mean; the solvated core
sits at the median/mode, and the MD n₁ (≈ 1.20 eV) has no Coulomb-fast
tail (§4n/I58: the MD n₁ KE is the co-moving shed convention, Coulomb
share +0.08 eV) — so scoring MD-mean (core) against ref-mean
(core + structural tail) is not like-for-like.

Re-score (oracle: baseline χ²_prof reproduces the recorded values exactly
— 102.94 vs 102.9 / 22.95 vs 22.9 / 162.42 vs 162.4 / 13.90 vs 13.9 /
62.85 vs 62.9). Four n = 1 treatments: **base** (n ≥ 1, mean), **drop**
(`n_min = 2`), **median** (ref n₁ → 1.128), **median+bg** (median with sys
err widened by the 0.378 bgOffShift), **mode** (→ 0.891):

| cell | base | drop | median | med+bg | n₁ share of χ² |
|---|---|---|---|---|---|
| s2c1v65t32e27 (histogram-landed) | 162.4 | 90.3 | 121.4 | 91.7 | **44 %** |
| s1c1v65 | 102.9 | 65.9 | 78.6 | 66.6 | 36 % |
| s2c1v85t32e23 (KE-best) | 13.9 | 12.6 | 13.4 | 12.6 | 9 % |

- **The n₁ anchor was inflated** — de-weighting removes 36–44 % of the
  low-drag KE χ² (the user's theory holds as a data fact). Median+bg ≈ drop
  to <2 %, so the correction is robust. The MD n₁ (1.20) matches the ref
  **median (×1.06)** far better than the **mean (×0.92)**; mode (0.89) is
  worse (×1.34) — the MD sits at the experimental *core*, as the mixture
  reading requires.
- **But the high-drag pull survives every treatment.** Stage-1 χ² argmin
  stays at the bracket top (c1/c4 → v85, c2 → v75) under base/drop/median/
  med+bg; Stage-2 KE-best stays at **(v85, τ3.2)** under all four. c1
  mean-χ² at v65 vs v85: 200 → 39 (base), 151 → 35 (drop) — v85 still wins
  ~4×. **The two-axis tension is driven by n2–n8, not n=1.**

Consequence: adopt median-anchored (or n=1-excluded) KE scoring going
forward — a like-for-like correctness fix (n=0 is already excluded under
I-D4) that removes a 40 %-of-χ² confound and *isolates* the residual to
the mid-band. It does **not** buy out of the taper decision.

### Part 2 — the E-2 taper-corner twin scan (p_tail is the wrong second lever)

Twin re-score at the histogram-landing budget (leg-D config: rq4graded/c1,
kornilov prior, density_tied dressing, p_onset = 1, co-moving KE; τ = 3.2,
E₀ = 0.27), sweeping the §4t-NB E-2 corner v_c ∈ {6.5, 7.0, 7.5} ×
p_tail ∈ {−1, −3}. m = 6000 (oracle: the leg-D c1 `d` config reproduces
the committed m = 20000 row — n̄ 4.428 vs 4.424, n₁ KE 0.9927 vs 0.9921).
`midHot` = geo-mean(sim/ref-mean) over n2–8 (want ≈ 1):

| v_c | p_tail | n̄ | n₁_solv | W₁ | midHot | n₁ KE |
|---|---|---|---|---|---|---|
| 6.5 | −1 | 3.69 | 0.263 | 0.451 | 1.744 | 1.165 |
| 7.0 | −1 | 4.14 | 0.253 | 0.481 | 1.240 | — |
| **7.5** | **−1** | 4.65 | **0.238** | 0.939 | **0.953** | **0.906** |
| 6.5 | −3 | 2.53 | **0.304** | 1.545 | **4.697** | 1.636 |
| 7.0 | −3 | 2.65 | 0.297 | 1.428 | 4.187 | — |
| 7.5 | −3 | 2.87 | 0.288 | 3.421 | 1.480 | 1.480 |

(ref: n₁_solv 0.310, n₁/n₂ 2.18; ref mean-KE n1..8 = 1.30/0.71/0.50/0.39/
0.32/0.27/0.23/0.19.) Per-bin at **v7.5/p−1**: n2 0.637 (ref 0.706), n3
0.469 (0.496), n4 0.366 (0.390), n5 0.307 (0.319) — the **entire mid-band
lands within ~10 %**, n₁ at the core (0.906 ≈ median 1.13).

Findings:

1. **p_tail is the wrong second lever — it re-couples the two axes.** More
   negative p (−3) re-heats the *whole* curve (every fragment passes the
   above-cap fast phase during ejection → less drag there lifts all bins:
   midHot 1.74 → 4.70 at v6.5) *and* re-strips (n̄ 3.69 → 2.53). KE and
   n₁_solv move **together**, not orthogonally. **No joint cell exists in
   the (v_c, p_tail) corner** — the E-2 move as hypothesised does not
   reconcile the axes; it relocates the miss.
2. **But raising v_c alone lands the mid-bin KE** (v6.5 → v7.5 at p−1:
   midHot 1.744 → 0.953, the full n2–n8 curve on the reference). **The
   mid-bin bulge is a real, fixable mid-band drag deficit** — I72's
   taper-shape statement, now shown constructively: the taper *magnitude*
   (cap position v_c) is the lever, not the tail *exponent* p.
3. **The residual at v7.5 is histogram over-retention** (n₁_solv 0.238 vs
   0.310), a *stripping-at-fixed-KE* deficit. KE and stripping are coupled
   through the cooling exposure K (more drag → more K → cooler KE **and**
   lower E_ej → less shedding → higher n): the **KE↔histogram
   anti-correlation through K** (Addendum-I Step 0–1) re-confirmed in the
   taper corner. The taper cannot break it because both effects are K.
4. **The decoupling lever is the ladder bottom (RQ4), not the taper.** The
   ladder moves n_det at ≈ fixed v_inf (KE-neutral to first order) — the one
   knob that can re-strip v7.5's cooled fragments back toward n₁_solv 0.31
   without re-heating. The two "genuinely free knobs" map cleanly onto the
   two axes: **v_c ↔ the KE curve, ladder ↔ n₁_solv** — the good
   identifiability case, *if* a cheap-enough ladder bottom can shed at
   v7.5's low E_ej. Next twin scan: (v_c ≈ 7.5 × ladder family).
5. **The twin's known bias favours MD here.** The twin *under-strips* at
   high v_c (I73/I74: MD − twin ≈ −0.4…−0.7 He), so MD n₁_solv at v7.5 sits
   **above** the twin's 0.238, toward 0.31 — the histogram gap is smaller in
   MD. Corroboration on disk: Stage-1 s1c1v75 (v_c 7.5, τ3.8/E₀0.25) already
   gave MD n₁_solv 0.250 with χ² halved to 54.4. The **un-scored
   (v7.5, τ3.2, E₀0.27) MD cell** — between the v65/v85 Stage-2 columns — is
   the twin-identified best joint candidate and the cheapest MD spot-check.

### NB (2026-07-21) — the (v_c × ladder) follow-up: orthogonality holds, but the ladder is already saturated at rq4graded; the histogram recovery is a twin→MD strip bias, not a ladder move

Finding-2's proposed decoupler (I78: ladder ↔ n₁_solv at v_c-set KE) was
scanned — same leg-D twin, p_tail = −1, τ 3.2 / E₀ 0.27, v_c ∈ {6.5, 7.0,
7.5, 8.0} × ladder ∈ {rq4graded, flat, floor1, slid2, slid3} (chords cache
per v_c; the ladder touches only the fate map). Bottom rungs D0(1..3)
[meV]: rq4graded 20.3/13.8/12.0, flat 9.2/9.2/9.2, floor1 13.3/9.2/9.2,
slid2 13.3/13.3/9.2, slid3 13.3/13.3/13.3.

At v7.5 (midHot ≈ 0.95 — the KE-landing v_c):

| ladder | n₁_solv | n₂_solv | ratio | bare | midHot |
|---|---|---|---|---|---|
| **rq4graded** | **0.238** | 0.127 | **1.87** | 0.198 | 0.953 |
| floor1 | 0.171 | 0.095 | 1.80 | 0.313 | 0.957 |
| slid2 | 0.164 | 0.138 | 1.18 | 0.284 | 0.956 |
| slid3 | 0.164 | 0.138 | 1.18 | 0.250 | 0.969 |
| flat | 0.125 | 0.102 | 1.23 | 0.341 | 0.967 |

1. **Orthogonality CONFIRMED.** `midHot` is flat across all five ladders at
   every v_c (v7.5: 0.95–0.97; v6.5: 1.74–1.77; v8.0: 0.82–0.83) — the
   ladder is genuinely KE-neutral. v_c ↔ KE and ladder ↔ histogram *are*
   separable axes (the good identifiability case).
2. **But the ladder is already saturated at rq4graded.** It gives the
   **highest** n₁_solv (0.238) and the best ratio (1.87 → target 2.18);
   every cheaper-bottom ladder strips *past* n = 1 into **bare** (bare
   0.20 → 0.34), *lowering* n₁_solv. rq4graded's expensive D0(1) = 20 meV
   barrier is load-bearing — it stops fragments at n = 1 instead of n = 0.
   I78's "cheap-bottom ladder re-strips to n₁" is **refuted**: cheap bottoms
   overshoot to bare. No ladder in the family reaches n₁_solv 0.31 at
   midHot ≈ 1.
3. **The histogram recovery is a twin→MD strip bias, not a ladder change.**
   The twin under-strips at high v_c, and the gap is *measured* at v6.5:
   twin n₁_solv 0.263 vs the on-disk MD cell s2c1v65t32e27 = **0.326**
   (+0.063 at matched n̄ — the MD's frozen-vs-live channel-(d) sharpens the
   n = 1 peak). Extrapolated to v7.5 (where channel (d) is larger, I73), MD
   n₁_solv should exceed the twin's 0.238 toward ≈ 0.30 — landing the
   histogram *and* the KE (midHot 0.95) at one cell. **rq4graded is already
   the correct ladder; the only ladder route to more n₁ is a still-more-
   bottom-heavy shape, which is RQ4-external territory (parked).**

**Consequence:** the decisive test is unchanged and now sharper — the MD
spot-check at **(v7.5, τ3.2, E₀0.27, rq4graded/c1)**. The twin says KE
lands, the ladder is right, and the twin→MD bias should carry n₁_solv from
0.238 to ≈ 0.31 there. If it does, that is the first joint (KE + histogram)
MD landing; if it doesn't, the two-axis miss is real inside the full
(v_c, p_tail, ladder) lever set and RQ4 (a more bottom-heavy ladder) or a
non-drag lever is required. Twin-only, direction-read; nothing discharges
F5.

### Boundaries

Zero-MD; reported, not auto-adjudicated. The twin KE axis is
**direction-only** (I74: KE is not twin-factorizable); `midHot` reads the
*sign and rough size* of the mid-band move, not a calibrated χ². Twin
m = 6000 (means stable; the leg-D oracle holds). n₁ re-anchoring is a
scoring-convention change proposed here, not yet adopted in the committed
scorer. Nothing here discharges F5.

---

## 4v. The (v7.5, τ3.2, E₀0.27) MD spot-check — the joint landing exists: the twin-identified cell lands the histogram AND the KE curve simultaneously, the first MD cell to do both; §4t's "no joint cell" was a sparse-v_c-grid artifact plus the n₁ anchor inflation

Executed 2026-07-21 under `[PROCEED TO IMPLEMENTATION]`. One MD cell,
`s2c1v75t32e27` — the c1 leg-D configuration (rq4graded, p_tail = −1,
kornilov prior, density_tied dressing, p_onset = 1, co-moving KE,
τ = 3.2, E₀ = 0.27, N = 50, production 2.70 eV, 8000 ps E2 cap) with only
v_c overridden 6.5 → 7.5 vs the on-disk v65 cell. Built by reusing the
committed `gen_tier2_md_confirmation` pins + `build_biphasic_cfg` (zero
repo-code change); **oracle: the saved cfg.json is byte-identical to
s2c1v65t32e27 except drag v_c** (asserted before compute). Full pipeline
neutral → ion → E2 → detection; 8/100 droplet-retained (excluded),
suppressed 0.141.

### The v_c progression at c1 / τ3.2 / E₀0.27 (MD, §4r scoring)

| cell | n̄_det | n₁_solv | n₁/n₂ | W₁_solv | midHot(n2–8) | n₁ KE | χ²_base | χ²_drop | χ²_med |
|---|---|---|---|---|---|---|---|---|---|
| v65 (histogram, no KE) | 3.63 | 0.326 | 2.80 | 0.736 | 1.669 | 1.197 | 162.4 | 90.3 | 121.4 |
| **v75 (NEW — both)** | **3.97** | **0.291** | **2.300** | **0.496** | **0.904** | **0.915** | **24.0** | **11.5** | **14.4** |
| v85 (KE, no histogram) | 4.74 | 0.167 | 0.765 | 0.951 | 0.792 | 0.663 | 34.1 | 23.4 | 27.2 |

(experimental targets: n₁_solv 0.310, n₁/n₂ 2.18; KE anchor mean 1.302 /
median 1.128 / mode 0.891.) Per-bin sim/ref-mean at v75: n2 ×1.01, n3
×0.95, n4 ×0.94, n5 ×0.84, n6 ×1.00, n7 ×0.90, n8 ×0.73, n9 ×0.87 — the
n2–n7 band lands within ~10 %.

### Findings

1. **The joint landing exists — at v7.5, the un-sampled v_c.** v75 lands
   the histogram (n₁_solv 0.291 ≈ 0.310; n₁/n₂ 2.30 ≈ 2.18; W₁_solv 0.496 —
   the best in the whole re-pilot at this budget) **and** the KE curve
   (midHot 0.904; n₁ 0.915 at the core; under the like-for-like
   median-anchored convention χ²_med = 14.4, essentially the grid-best
   13.9 — but *with* the histogram landed, not sacrificed). It is the first
   MD cell to satisfy both axes. §4t's "no (v_c, τ, E₀) cell lands both in
   capped_cubic" was a **sparse-grid artifact**: Stage-2 sampled only
   v_c ∈ {6.5, 8.5}; the joint optimum sits between them, exactly where the
   §4u twin + strip-bias pointed.
2. **The twin's quantitative prediction + the strip-bias correction were
   right.** The twin (I78/I79) predicted v75 n₁_solv 0.238, midHot 0.953,
   n₁ KE 0.906; MD delivered **0.291 / 0.904 / 0.915**. The +0.053
   histogram lift is the twin→MD under-strip bias (I79), the predicted
   direction and ≈ magnitude — the twin's ordering authority + the measured
   bias jointly located the cell before it was run.
3. **The mid-bin bulge is gone.** midHot 1.669 (v65) → 0.904 (v75): raising
   v_c pulls the mid-speed band under full cubic drag and cools n2–n8 onto
   the reference, exactly the §4u mechanism — and the histogram survives
   because the strip bias holds n₁_solv near target. The two axes are *not*
   irreconcilable in capped_cubic; they meet at v7.5.
4. **Half the raw KE χ² is still the n₁ mixture-mean artifact** (χ²_base
   24.0 → χ²_drop 11.5): n₁ contributes ~52 % under the raw-mean anchor,
   ~0 under the median anchor — I77 re-confirmed, and the reason the
   like-for-like read (χ²_med 14.4) is the honest one.

### Boundaries

N = 50 single seed — a **pilot** spot-check; the formal joint-landing
verdict is the N = 500 finalist read (RP-D6). The joint landing is under
the median-anchored n₁ KE convention (I77, proposed not yet wired). A few
deep/thin bins run cold (n8 ×0.73, n10 ×0.78, n11 ×0.63 — 2–6 ions/bin).
`zero_gamma` E2, 8000 ps cap; ledger not re-checked. Nothing here
discharges F5.

---

## 4w. The v_c sensitivity ring around 7.5 — the joint landing is a *basin* [7.25, 7.5], not a knife-edge; the optimum refines toward v_c ≈ 7.25 on the histogram (exactly as pre-registered); the N=50 χ²_med is thin-bin-noisy and must be read against the robust W₁/midHot metrics

Executed 2026-07-21 under `[PROCEED TO IMPLEMENTATION]`. The pre-registered
S3r ring: 4 new N=50 cells v_c ∈ {7.0, 7.25, 7.75, 8.0} (tags v700/v725/
v775/v800), one-cell-per-process, each cfg-oracle byte-matched to the v65
cell except v_c (all 4 passed). With the on-disk 6.5/7.5/8.5 anchors, the
7-point curve at c1/τ3.2/E₀0.27/leg-D (median-anchored n₁, I77):

| v_c | n₁_solv | n₁/n₂ | W₁_solv | midHot | n₁ KE | χ²_med | χ²_base |
|---|---|---|---|---|---|---|---|
| 6.50 | 0.326 | 2.80 | 0.736 | 1.669 | 1.197 | 121.4 | 162.4 |
| 7.00 | 0.341 | 3.11 | 0.586 | 1.162 | 1.056 | 37.5 | 52.5 |
| **7.25** | **0.309** | **2.27** | **0.424** | **0.987** | 0.987 | (66)* | (73)* |
| **7.50** | 0.291 | 2.30 | 0.496 | 0.904 | 0.915 | **14.4** | 24.0 |
| 7.75 | 0.256 | 2.00 | 0.611 | 0.947 | 0.853 | 20.7 | 30.6 |
| 8.00 | 0.256 | 2.00 | 0.649 | 0.878 | 0.804 | (139)* | (157)* |
| 8.50 | 0.167 | 0.77 | 0.951 | 0.792 | 0.663 | 27.2 | 34.1 |

(targets: n₁_solv 0.310, n₁/n₂ 2.18, midHot ~1. *= N=50 thin-bin χ² spikes,
off the smooth W₁/midHot trend — see finding 3.)

### Findings vs the pre-registered S3r predictions

1. **A joint BASIN exists, not a knife-edge (P1 confirmed on robust
   metrics; P3 confirmed).** On the noise-robust reads — W₁_solv, n₁_solv,
   n₁/n₂, midHot — both axes land smoothly across **v_c ∈ [7.25, 7.5]** (and
   marginally to 7.0/7.75): W₁_solv is a clean bowl minimizing at 7.25
   (0.424), midHot slides monotonically 1.67 → 0.79 crossing ~1.0 at ≈ 7.2,
   n₁_solv slides 0.33 → 0.17 crossing the 0.31 target at ≈ 7.25. v75 is at
   the center of the basin, not its edge.
2. **The joint optimum refines to v_c ≈ 7.25 on the histogram (P2 confirmed
   exactly).** n₁_solv = 0.309 at v7.25 — dead on the 0.31 target, with the
   best W₁_solv (0.424) and n₁/n₂ 2.27 ≈ the 2.18 target, and midHot 0.987
   (KE on-target). The histogram optimum sits just below 7.5, precisely the
   P2/P3 prediction. v7.5 remains the χ²-best (14.4) and is jointly strong;
   the two straddle the optimum.
3. **The N=50 χ²_med is thin-bin-noisy — read the robust metrics.** χ²_med
   is non-monotone with off-trend spikes at v7.25 (66) and v8.0 (139) that
   break the otherwise smooth curve. Per-bin inspection confirms these are
   **not** real KE misses: v7.25's resolved bins scatter tightly around 1.0
   (n2 ×1.07, n3 ×1.06, n5 ×1.00, midHot 0.987) — the χ² is inflated by
   2–4-ion deep bins (n6 ×1.33/4 ions, n7 ×0.60/3 ions) weighted by the
   reference's small errors. At N=50 the squared-deviation χ² is dominated
   by thin-bin scatter; W₁_solv and midHot (outlier-robust) are the reliable
   discriminants until N=500 stabilises the deep bins. **The strict
   pre-registered χ²_med ≤ 30 gate therefore under-counts the window (flags
   only v7.5) — a criterion artifact, not a knife-edge.**
4. **v8.0+ is genuinely over-dragged.** Distinct from the noise: v8.0's
   whole per-bin curve runs cold (×0.68–0.88) and n₁_solv falls to 0.256 —
   the histogram degrades and the KE overshoots to the cold side. The basin
   closes above ~7.75.

### Consequence — the c1 N=500 finalist

The joint landing is robust across **v_c ∈ [7.25, 7.5]** (outcome (a)+(c):
basin confirmed, optimum refined below 7.5). The c1 N=500 finalist should
center **v_c ≈ 7.25–7.5** (7.375 ± the ring as its sensitivity leg), read
on W₁_solv/midHot, with χ² re-adjudicating 7.25 vs 7.5 once N=500 stabilises
the deep bins. Adopting the median-anchored n₁ KE convention (I77) into the
committed scorer remains the prerequisite.

### Boundaries

N = 50 single seed — W₁/n₁_solv/midHot are the stable reads; χ²_med is
thin-bin-noisy (finding 3) and per-cell differences within [7.25, 7.5] are
below the N=50 resolution. median-anchored n₁ (I77); `zero_gamma` E2,
8000 ps cap. Nothing here discharges F5.

---

## 4x. S4 — c4/floor1 twin-first finalist location: a c4 joint region EXISTS at (v7.0–7.25, τ3.8, E₀0.25) — n₁_solv lands dead-on in MD (0.312/0.303) with midHot in-band, W₁ stays ≈ 0.63–0.65 (c1's basin: 0.42–0.50) so the N=500 ladder contest is real, not a strawman; the twin→MD strip lift is budget-dependent on the floor1 arm (+0.10 at E₀0.25, ≈ 0 at the high-suppression E₀0.27 corner)

Executed 2026-07-21 under `[PROCEED TO IMPLEMENTATION]` (S4 + S5 in
parallel; plan §I.11.5.1). **Twin scan** (m = 6000, stage_repilot2 draw
discipline verbatim, scratchpad driver importing the committed module —
zero repo-code change): v_c ∈ {6.5..8.5, 0.25-spaced} × τ ∈ {3.2, 3.8} ×
E₀ ∈ {0.23, 0.25, 0.27} on floor1/p−1, with the rq4graded twin riding the
same nine chords (the ladder enters only the fate map). Committed:
`h2b_s4_c4scan_{predictions,ke}.csv` (m = 6000, 108 rows) +
`h2b_s4_c4_final_{predictions,ke}.csv` (m = 20000 refreeze). Then
**4 MD confirm cells** (N = 50, τ3.8, zero_gamma E2, leg-D arms; each cfg
byte-oracle vs its on-disk s2c4v65t38e{25,27} sibling — v_c-only diffs,
all passed).

### Verdicts (pre-registered S4c-P1..P4)

1. **S4c-P1 CONFIRMED.** The m = 6000 anchor cell reproduces the
   committed m = 20000 dc4 row: n̄ 4.354 vs 4.359 (Δ −0.005), n₁ KE
   0.9688 vs 0.9668 (Δ +0.0020) — tighter than the §4u c1 precedent.
2. **S4c-P2 CONFIRMED (KE-axis ladder-blindness).** max |midHot(c4) −
   midHot(c1)| ≤ 0.045 over the whole grid, shrinking with v_c (0.045 at
   v6.5 → 0.010 at v8.5). The c4 midHot ≈ 1 crossing sits at the same
   v_c as c1's.
3. **S4c-P3 SPLIT.** The c4 < c1 n₁_solv ordering holds on the
   E₀ ∈ {0.25, 0.27} region but **inverts at the E₀ = 0.23 corner** —
   there c1's suppressed channel closes entirely (supp ≈ 0.000 on every
   rq4graded E₀0.23 cell) and its n₁_solv collapses to ≈ 0.16 with
   ratio ≈ 1.0; floor1 keeps a small suppressed population and sits
   above it. A budget-corner regime change, not sampling noise.
   The registered *expectation* (bias-corrected c4 peak below the
   window ⇒ outcome (b)) was **wrong**: candidates exist — at τ3.8, not
   the c1-like τ3.2.
4. **S4c-P4 REFINED.** The candidate region is contiguous but sits at
   **(v_c ≈ 6.75–8.0, τ3.8)** with the best cells at v7.0–7.25 —
   overlapping the registered v7.25–7.5 guess only at its low edge (the
   floor1 histogram needs the longer τ to populate n₁).

### The MD confirm quartet (§4r scoring, median-anchored n₁)

| cell | supp | ret | n̄_det | n₁_solv | n₁/n₂ | W₁_solv | midHot | n₁ KE | χ²_med |
|---|---|---|---|---|---|---|---|---|---|
| s2c4v65t38e27 (anchor) | 0.424 | 1 | 2.68 | 0.211 | 1.20 | 0.797 | 1.275 | 0.935 | 47.8 |
| s2c4v65t38e25 (anchor) | 0.212 | 1 | 3.41 | 0.282 | 2.44 | 0.811 | 1.504 | 1.153 | 147.9 |
| s4c4v675t38e27 | 0.427 | 4 | 2.58 | 0.236 | 1.62 | 0.878 | 1.056 | 0.845 | 167.6 |
| s4c4v700t38e27 | 0.432 | 5 | 2.79 | 0.185 | 1.11 | 0.873 | 0.894 | 0.776 | 36.8 |
| **s4c4v700t38e25** | 0.189 | 5 | 3.49 | **0.312** | 3.00 | 0.651 | **1.121** | 1.019 | 61.9 |
| **s4c4v725t38e25** | 0.191 | 6 | 3.64 | **0.303** | 3.29 | **0.630** | **0.975** | 0.945 | 77.1 |

(targets: n₁_solv 0.310, n₁/n₂ 2.18, midHot ~1, n₁ KE median-anchor
1.128.)

### Findings

1. **The c4 joint region exists in MD — outcome (a) at the MD level.**
   (v7.0, τ3.8, E₀0.25) lands n₁_solv 0.312 (target 0.310, dead-on) with
   midHot 1.121; (v7.25, τ3.8, E₀0.25) lands 0.303 / 0.975. Both pass
   the robust joint bars (n₁_solv ∈ [0.26, 0.36] ∧ midHot ∈ [0.80, 1.20]);
   both fail the strict χ²_med ≤ 30 gate at N = 50 — the §4w
   thin-bin-noise pattern (I81), read on W₁/midHot until N = 500.
   Under the frozen tie-break (lower χ²_med subject to the robust bars)
   the **S4-located c4 finalist is (v7.0, τ3.8, E₀0.25)** — χ² 61.9 vs
   77.1 and the exact n₁_solv — with v7.25/e25 the adjacent basin cell.
2. **The ladder contest at N = 500 is real, not a strawman — and c1
   still leads.** c4's best W₁_solv is 0.630–0.651 vs c1's basin
   0.424–0.496, and its ratio overshoots (3.0–3.3 vs 2.18 — n₂
   under-filled where c1's v7.25 gives 2.27). c4 now enters the RP-D6
   arbiter at its own basin-located best, so a c1 win at N = 500 is a
   defensible bounded-physics verdict (I51/S2c-P3), not a sparse-grid
   artifact.
3. **The twin→MD under-strip lift is NOT universal (channel-(d)-family
   structure).** Measured lifts on floor1: **+0.100 / +0.097** at the
   E₀0.25 cells (twin 0.212/0.206 → MD 0.312/0.303) — ≈ double the
   c1-measured +0.05..0.065 band — but **+0.019 / −0.022** at the
   E₀0.27 cells, where the suppressed channel is heavily populated
   (supp ≈ 0.43). The lift collapses when suppression is large. The
   §I.11.5.1 discipline (bias bands as *search guidance only*, never a
   correction) is vindicated: the band found the region, MD placed the
   cells.
4. **Drag exposure moves the retained class on c4 too:** retained 1/100
   at the v65 anchors → 4–6/100 at v6.75–7.25, the §4v/§4w c1 pattern.

### Boundaries

N = 50 single seed on the MD quartet (ratio reads carry SE ≈ ±0.5–0.7;
n₁_solv is the stable read); twin scan m = 6000 refrozen at m = 20000 for
the named cells only; the E₀0.23 corner regime change is a twin-level
read (no MD cell there); `zero_gamma` E2 (the S5 outcome licenses the
N = 500 Landau configuration separately); the finalist designation
follows the frozen outcome-(a) rule and is **reported** — the S6
pre-registration freezes it. Nothing here discharges F5.

---

## 4y. S5 — the Landau v_L bracket at the §4v cell: the gate DISSOLVES — every scored observable is bit-flat across v_L ∈ {0.30, 0.40, 0.58} while the arm demonstrably acts on the retained class (E_dissip +0.58–0.60 eV; retained KE 0.277 → 0.005–0.023 eV, monotone in v_L); pre-E2 identity holds by hash despite the E_min shared reader

Executed 2026-07-21 under `[PROCEED TO IMPLEMENTATION]` (plan §I.11.5.2).
Three N = 50 cells at the s2c1v75t32e27 pins with
`relaxation_dissipation = "landau_gated_drag"`, v_L ∈ {0.30, 0.40, 0.58}
Å/ps (`v_limit_m_per_s` 30/40/58), vs the on-disk zero_gamma baseline.
Scratchpad driver reusing the committed generator (zero repo-code
change).

**A latent shared reader was identified at build time and measured:**
`v_limit_m_per_s` also feeds the *neutral-stage* collision threshold
through the derived `cfg.E_min_eV` (propagation_step.py /
collisions.sample_collision_events), so a v_L change could in principle
perturb pre-E2 physics. At these kinematics it does not fire: **on all
three arms `neutral.npz` and `ion.npz` are sha256-identical to the
baseline** (the eV-scale fragment energies never cross the 1.05 ↔ 2.21
meV threshold window). S5L-P2 therefore holds **by measurement**, and
the S5L-P1 cfg oracles passed on every arm (l40 diff =
{relaxation_dissipation} only; l30/l58 + {v_limit_m_per_s}).

### The bracket read (§4r scoring, median-anchored n₁)

| arm | ret | supp | n̄_det | n₁_solv | n₁/n₂ | W₁_solv | midHot | χ²_med | E_dissip(ret) [eV] | E_kin(ret) [eV] |
|---|---|---|---|---|---|---|---|---|---|---|
| zero_gamma (base) | 8 | 0.141 | 3.97 | 0.291 | 2.30 | 0.496 | 0.909 | 14.4 | 21.179 | 0.277 |
| landau v_L 0.30 | 8 | 0.141 | 3.97 | 0.291 | 2.30 | 0.496 | 0.909 | 14.4 | 21.784 | 0.0053 |
| landau v_L 0.40 | 8 | 0.141 | 3.97 | 0.291 | 2.30 | 0.496 | 0.909 | 14.4 | 21.777 | 0.0093 |
| landau v_L 0.58 | 8 | 0.141 | 3.97 | 0.291 | 2.30 | 0.496 | 0.909 | 14.4 | 21.761 | 0.0229 |

### Verdicts (pre-registered S5L-P1..P4)

1. **S5L-P1 CONFIRMED** (cfg byte-oracles, above).
2. **S5L-P2 CONFIRMED by hash** (pre-E2 identity; the shared-reader
   hazard measured quiet).
3. **S5L-P3 CONFIRMED, stronger than registered.** Not merely "below
   N = 50 resolution": Δn₁_solv = 0.0000, ΔW₁ = 0.0000, Δretained = 0 on
   every arm; ΔmidHot = −0.0002 (the listed residual coupling — the erf
   spatial gate's tail just outside the surface — four orders below the
   bar). χ²_med 14.4 on all four arms, which simultaneously re-confirms
   the recorded §4v value through the **committed** median-anchor scorer.
4. **S5L-P4 CONFIRMED — quiet is a dissolved gate, not a silent arm.**
   The retained class books +0.605/+0.598/+0.582 eV of E2 drag
   dissipation and its residual KE collapses 0.277 → 0.0053/0.0093/
   0.0229 eV — monotone in v_L exactly as the physics demands (lower
   v_L = wider gate = colder endpoint).

### Findings

1. **Outcome (a): the `v_L` gate dissolves.** The scored observables are
   structurally insulated from the Landau arm (ejected class:
   g → 0 outside; retained class: excluded from every IHe_n read), and
   the bracket confirms the insulation is airtight at the detected
   surface across the full credible v_L range. Consequence (per the
   frozen outcome shape): **N = 500 runs Landau-on at v_L = 0.58** (the
   physics-motivated bulk ceiling), the legacy 0.40 becomes the
   winner-cell sensitivity spot (RP-D7 mitigation, now data-backed), and
   the external re-pinning reply is demoted to a recorded spot value —
   it no longer gates anything.
2. **The arm's physical output is the retained-class endpoint** (colder
   by ~0.25–0.27 eV, v_L-graded) — invisible to Tier-2 scoring but the
   physically-correct configuration for the marginal class the N = 500
   arbiter carries.

### Boundaries

One cell (v7.5), one seed — the bracket licenses the *gate*, not a v_L
calibration; nothing measures which v_L is physically right (that
remains the external question, now non-blocking). The retained-class
endpoint state is diagnostic-only under the exclude policy. The E_min
shared-reader quietness is measured **at these kinematics** — a future
low-energy channel would need the hash check repeated (recorded as an
execution convention, not a code change). Nothing here discharges F5.

---

## 4z. S6 — the N = 500 finalists: NO cell lands the frozen joint acceptance — N = 500 unmasks a systematic KE-curve miss that the N = 50 sim-SE-dominated sigmas had absorbed (χ²_med RISES everywhere, refuting S6f-P3's direction, with the mechanism identified); the robust histogram carry holds (P2: v7.25 fully in-band) and the RP-D6 ladder contest resolves decisively for c1/rq4graded (P4, ≫ 2σ); the winner-gated riders do not fire

Executed 2026-07-21 under the S6 `[PROCEED TO IMPLEMENTATION]`
(pre-registration plan §I.11.5.5 committed before any MD). Three N = 500
cells, fresh seed 20260721, Landau-on v_L = 0.58, production 2.70 eV,
8000 ps E2 cap, 8.53 µs detection; every cfg oracle passed (diff exactly
the amended `{num_molecules, seed, relaxation_dissipation,
v_limit_m_per_s}` set); scored by the committed median-anchor scorer
(mean-legacy alongside).

### The scored read (targets: n₁_solv 0.310, n₁/n₂ 2.18, midHot ≈ 1, χ²_med ≤ 30)

| cell | scored | trapped | supp | n̄_det | n₁_solv | n₁/n₂ | W₁_solv | midHot | χ²_med | χ²_mean | bars |
|---|---|---|---|---|---|---|---|---|---|---|---|
| finc1v725 | 943 | 0.057 | 0.161 | 4.003 | **0.2718** | 2.048 | 0.5238 | **0.9905** | 125.7 | 132.0 | n₁ ✓ midHot ✓ χ² ✗ |
| finc1v750 | 935 | 0.065 | 0.162 | 4.252 | 0.2551 | 1.923 | 0.6572 | 0.9095 | **68.6** | 73.2 | n₁ ✗(−0.005) midHot ✓ χ² ✗ |
| finc4v700 | 952 | 0.048 | 0.249 | 3.537 | 0.2322 | 1.953 | 0.8426 | 1.1018 | 202.5 | 219.1 | n₁ ✗ midHot ✓ χ² ✗ |

### Verdicts (pre-registered S6f-P1..P5)

1. **S6f-P1 CONFIRMED** (cfg oracles in-process + independent diff check;
   scorer lock pre-registered in §I.11.5.5).
2. **S6f-P2 SPLIT.** n₁_solv and midHot carry at every cell (v7.25
   0.2718 ∈ [0.206, 0.411]; every midHot within ± 0.15 of its N = 50
   value). W₁_solv exceeds its ± 0.12 band at v7.5 (0.657 vs
   [0.376, 0.616]) and c4 (0.843 vs [0.531, 0.771]) — and the worsening
   is same-direction at all three cells (0.424 → 0.524, 0.496 → 0.657,
   0.651 → 0.843): an N-resolution tail structure emerging, not pure
   seed luck (v7.25, fully in-band, is the seed-robustness read that
   holds).
3. **S6f-P3 REFUTED — with the mechanism identified.** χ²_med rose at
   every cell (66.3 → 125.7, 14.4 → 68.6, 61.9 → 202.5). The N = 50
   per-point sigma was **dominated by the sim-side SE** (0.03–0.08 eV);
   at N = 500 it collapses to the reference's own error (0.010–0.03 eV)
   and the same fractional miss costs ~5–10× more χ², while five more
   deep bins enter the fit (ke_npts 12 → 17). The "thin-bin spike"
   diagnosis (I81) had the mechanism backwards: N = 500 does not smooth
   χ² down — it **unmasks a systematic KE-curve miss** the pilot
   statistics could not see.
4. **S6f-P4 CONFIRMED, decisively — the ladder verdict.** c1 beats c4 on
   every histogram read at N = 500: W₁_solv 0.524/0.657 vs 0.843,
   χ²_med 125.7/68.6 vs 202.5, n₁_solv 0.272/0.255 vs 0.232
   (c1 − c4 ≈ +0.04 ≈ 2σ at SE ≈ 0.015 — the I76 direction resurfaces
   at N = 500 after the S4 N = 50 near-tie). The registered direction
   (W₁: c1 < c4) lands far beyond 2σ; the I51/S2c-P3 bounded-physics
   claim is formally read: **rq4graded > floor1**.
5. **S6f-P5: 2/3 in-band.** Δn̄(MD − twin): v7.25 −0.384 ∈
   [−0.63, −0.33] ✓; c4 −0.264 ∈ [−0.46, −0.16] ✓; v7.5 −0.395
   **outside** its own [−0.83, −0.53] band but **inside** the original
   §I.11.5.3 drag-exposure interpolation guess (−0.3..−0.6): the N = 50
   measurement (−0.68) was the outlier and N = 500 regresses to the
   interpolation. Tracked never corrected (I69/I83).

### The KE-miss anatomy (per-bin z², committed error model, median anchor)

- **c1 cells — coherent deep-bin undershoot.** Sim mean KE at n ≥ 12
  sits 30–60 % below the reference (v7.25: n = 16 sim 0.028 vs ref
  0.069 eV, z² 59.7; the n ≥ 12 tail carries ≈ 91 of 125.7; v7.5:
  ≈ 37 of 68.6), plus a mild n = 2–4 overshoot at v7.25 (z² 11/9/6).
  The n = 1 bin is clean under the median anchor (z² 0.1/0.0) — the
  I77 convention did its job.
- **c4 — whole-curve rotation.** Hot at n = 2–4 (z² 28/31) and cold
  from n ≥ 8 (z² 8–24 per bin): the floor1 KE curve fails
  structurally, corroborating the histogram-side rejection.
- **Physics read:** the model's deeply-solvated survivors arrive **too
  cold** at the detector — invisible at N = 50 (2–6 ions/bin), now the
  sharpest experimental constraint the program has produced. Whether
  this is an E2/exposure artifact, a ladder-tail artifact, or a real
  relaxation-channel gap is NOT adjudicated here.

### Formal outcome

Under the frozen rule (winner = lower χ²_med subject to all three
bars) **no finalist lands and no winner exists; the winner-gated riders
(Stage-3 sensitivity ring + the 0.40 Landau spot) do not fire.**
finc1v725 is the best cell on the robust bars (both n₁_solv and midHot
in-band; midHot 0.990 dead-on); finc1v750 holds the lowest χ²_med
(68.6) but sits 0.005 below the n₁_solv window. Reported, not
auto-adjudicated: the F5 reconciliation (discharge vs re-scope) and any
deep-bin KE follow-up are the next user adjudication.

### Boundaries

Single seed per cell (the P2 carry is the seed read); Landau-on per the
S5 license (zero_gamma untested at N = 500; S5 measured the scored
surface insulated at N = 50); the n ≥ 14 sim means still carry 2–6
ions/bin at N = 500; the suppressed/bare channel stays RQ8-gated.
Nothing here discharges F5.

---

## 4aa. RP-D4 riders at the blessed production point (D1–D2 adjudication) — the 0.40 Landau spot is bit-flat (RP-D7 mitigation complete); every sensitivity-ring perturbation DEGRADES W₁ (+0.18..+0.33: the pinned leg-D configuration sits at the basin optimum); the birth margin is the sensitive lever (6 Å moves n₁_solv −0.17, out of window); N = 50 χ² deltas carry the I85 caveat

Executed 2026-07-21 under `[PROCEED TO IMPLEMENTATION]` (the D1–D4
adjudication block, log entry "PRODUCTION POINT ADJUDICATED"). Six
N = 50 cells, **one lever each** off the on-disk `s2c1v725t32e27`
baseline (same seed 20260604 — paired A/B); all six cfg oracles passed
with exact single-key diffs (`r725l40`'s diff is
`{relaxation_dissipation}` alone — v_limit 40 is the config default).

### The ring read (§4r scoring, median anchor; deltas vs baseline)

| cell (lever) | supp | n̄_det | n₁_solv | n₁/n₂ | W₁_solv | midHot | χ²_med |
|---|---|---|---|---|---|---|---|
| baseline s2c1v725t32e27 | 0.138 | 3.904 | 0.3086 | 2.273 | 0.4240 | 1.0163 | 66.3 |
| r725d040 (δ 0.40) | 0.138 | 4.372 | 0.2346 | 1.583 | 0.5993 | 1.0886 | 32.6 |
| r725d080 (δ 0.80) | 0.137 | 3.695 | 0.3293 | 2.250 | 0.6081 | 1.0008 | 9.5 |
| r725pick (pickup prior) | 0.143 | 4.297 | 0.2308 | 1.636 | 0.6617 | 1.1179 | 26.8 |
| r725m467 (margin 4.67 Å) | 0.116 | 4.347 | 0.2143 | 1.385 | 0.6006 | 1.1038 | 31.8 |
| r725m600 (margin 6.0 Å) | 0.084 | 4.632 | 0.1379 | 0.571 | 0.7513 | 1.1563 | 40.3 |
| r725l40 (Landau v_L 0.40) | 0.138 | 3.904 | 0.3086 | 2.273 | 0.4240 | 1.0161 | 66.6 |

### Findings

1. **The 0.40 Landau spot is bit-flat — RP-D7 mitigation complete.**
   Δn₁_solv = Δratio = ΔW₁ = Δn̄ = 0.0000; ΔmidHot −0.0002 (the §4y
   erf-tail residual); Δχ²_med +0.36. Combined with S5's three-point
   bracket at v7.5 and the N = 500 Landau-on production run, the
   `v_L` axis is measured quiet at the winner cell family on the
   scored surface across its full credible range.
2. **Every perturbation degrades W₁_solv** (+0.175 δ0.40 / +0.184
   δ0.80 / +0.238 pickup / +0.177 m4.67 / +0.327 m6.0): the pinned
   leg-D configuration sits at/near the basin optimum of every
   robustness lever swept. Reported as bands, never re-fit (RP-D4).
3. **The birth margin is the sensitive lever.** n₁_solv −0.094 at
   4.67 Å and −0.171 at 6.0 Å (≈ 1.9σ / 3.4σ at the N = 50
   SE ≈ 0.05), with the ratio collapsing (2.27 → 0.57) and the
   suppressed channel halving (0.138 → 0.084). **Recorded caveat: the
   leg-D margin 3 Å is a pinned convention and the histogram landing
   depends on it at the multi-σ level.**
4. **The droplet-prior axis is asymmetric.** Narrower (δ 0.40) and the
   pickup-weighted family are material (n₁_solv −0.074 / −0.078,
   ≈ 1.5σ); wider (δ 0.80) is noise-level (+0.021) — the D4 primary
   (kornilov δ 0.625) is not knife-edge in the wide direction.
5. **N = 50 χ²_med deltas are direction-only** (the I85 sim-SE
   unmasking caveat): δ0.80's 9.5 is NOT evidence a ring cell beats
   the baseline at N = 500 statistics.

### Boundaries

N = 50, single seed, one lever at a time (no interaction terms);
bands never re-fit; the ring rode `zero_gamma` (clean one-lever cfg
diffs — licensed by finding 1 + §4y). Nothing here discharges F5.

---

## 4bb. RQ11 exposure diagnostic (read-only) — the cold tail is set inside the 30 ps MD window: candidate (i) E2 exposure is REFUTED (E_dissip gain 0.0000 on every solvated bin; KE frozen at handover), the deficit accrues in the in-band cubic phase ≈ 2–14 ps at v ≈ 2–7 Å/ps, and fate is birth-dressing/droplet-size-ordered

Executed 2026-07-22, read-only (two scratchpad scripts over the on-disk
`finc1v725` npz surface; no repo code, no MD — the RQ11 "candidate
resolution: onset-n anatomy + velocity-class decomposition, read-only"
pass). Convention check: raw `state_reason`-masked per-n means reproduce
the committed scorer's deep bins at printed precision (n = 16: 0.0285 vs
the recorded 0.028 eV; census frozen 742 / time_exhausted 49 /
suppressed 152 / droplet_retained 57 = the §4z scored 943). NB the
on-disk `n_detected` rides at handover n for the suppressed/retained
classes — the 0/NaN forcing lives in `DetectedEnsembleView`; a raw read
that ignores `state_reason` pollutes n ≥ 12 with both classes.

### The stage-localization read (relaxation.npz + detection.npz)

Per detected solvated bin (means): KE at handover ≈ KE at E2-end ≈ KE at
detection everywhere — n ≥ 12: 0.0485 → 0.0557 → 0.0556 eV (the small
early rise is shed recoil); KE settle times 30–70 ps into the 8000 ps
E2 window; **zero** fragments still KE-decaying at the cap; cumulative
`E_dissip` gain across E2 = **0.0000 eV on every solvated bin**. The E2
Landau-gated drag acts only on the in-droplet retained class (as §4y
measured); escaped ions coast. Detection-stage sheds are a low-n
phenomenon: all 140 events sit on fragments detected at n ≤ 15,
overwhelmingly n ≤ 10, and move KE negligibly (n = 1: 1.0340 → 1.0299).
The Landau floor is irrelevant to the scored deep bins (KE at v_L 0.58
≈ 0.0030–0.0036 eV for n = 12–20 vs their ≈ 0.055 eV) — only the
excluded retained class sits on it (0.0025 eV, n ≈ 16.7), which
*explains* the §4y/§4aa v_L-quietness of the scored surface. E2 trims
the mapping near-diagonally: deep survivors go handover n 14.4 →
detected 13.4 (≈ 1 shed in 8.5 µs); n = 1 arrives from handover n 2.0.

### The in-window read (ion.npz, 0–30 ps)

All fate groups launch identically (peak KE 1.69–1.80 eV, ≈ 13 Å/ps at
0.5 ps) and differentiate by dressing: birth n_shell 15.5 (future
n = 1) → 18.4 (future n ≥ 12), initial droplet radius 25.7 → 30.1 Å
(retained: 34.0 Å) — fate is exposure-ordered end to end (total
E_dissip 1.258 → 2.576 → 2.676 eV monotone across the groups). For the
60 future n ≥ 12 survivors the above-cap phase ends at t ≈ 2.1 ps
(spread 1.94–2.31; 79 % of their E_dissip accrues above v_c 7.25, vs
99 % for the n = 1 group), but they fall below the experimental
deep-bin KE level (0.069 eV, the n = 16 reference) only at mean
13.5 ps (range 6.2–30, 60/60 cross): **the deficit-critical cooling is
the in-band pure-cubic segment ≈ 2–14 ps at v ≈ 2–7 Å/ps**, followed
by a flat coast (mild recoil recovery 0.046 → 0.049 eV after 20 ps).
Every solvated fragment is outside its droplet by 30 ps (deep group at
r ≈ 80 Å ≈ 2.7 droplet radii); the retained class is 77 % inside.

### What this does to RQ11's candidate list

- **(i) E2 exposure/over-cooling — REFUTED.** The E2 cap is causally
  disconnected from the KE curve; a longer/shorter cap changes nothing
  for detected survivors. A large-N re-run needs no cfg change on this
  axis (diff stays `{num_molecules, seed}`).
- **(ii) ladder tail / freeze-out ordering — alive, reframed:** the
  mapping is near-diagonal post-window, so the lever is the *in-window*
  (KE, n) exit correlation, not late-time re-mapping.
- **(iii) missing recoil/relaxation channel — alive, sized:** deep
  survivors shed ≈ 5 He total (18.4 → 13.4); closing the 0.02–0.04 eV
  deficit needs only ~4–8 meV recoil per shed (RQ2's per-shed ε).
- **(iv) NEW — deep-bin population channel:** experiment's n ≥ 12 mass
  may draw on the excluded droplet-retained class via µs droplet
  evaporation (RQ5); in the MD frame those are colder still
  (0.0025 eV), so this explanation requires a beam-frame/detection
  argument before it can be scored.

### Boundaries

Single run, single seed, N = 500 (60 deep solvated fragments); raw-npz
means, not a `DetectedEnsembleView` re-derivation (spot-checked against
the committed scorer at n = 16 and the class census only); the
0.069 eV crossing threshold applies the n = 16 reference level
group-wide; scripts scratchpad-only per the S4/S6 driver precedent.
Nothing here discharges F5.

---

## 4cc. The N = 5000 battery (5 × N = 1000 fresh seeds) — BN-P2 fully in-band (the histogram carry is seed-robust), the W₁ scare was seed scatter (SD ± 0.095 measured; pooled 0.571 = the converged branch), the blessed cell's draw was mildly favorable (≈ 2σ on n₁/supp — winner's curse quantified, not dramatic), and the RQ11 deficit is a SLOPE in n (sim/ref 0.81 → 0.38 across n = 10–17), not a uniform scale

Executed 2026-07-22 under the BN `[PROCEED TO IMPLEMENTATION]`
(pre-registration frozen and committed before any MD; launch record in
the log). Five sequential N = 1000 cells, seeds 20260722–26, cfg diff
exactly `{num_molecules, seed}` off the on-disk finc1v725 cfg
(oracle-checked at every launch); ≈ 26–38 min/cell. Execution note:
the harness background-task layer killed the driver four times on a
per-cell cadence (machine healthy — no OOM, no crash, clean event
log); cells resumed loss-free from their fixed seeds, and the final
two cells ran in a harness-detached process. Scored by the committed
median-anchor scorer per cell; the pooled read via
`read_confirmation_detection` on the concatenated
(state_reason, n_detected, E_kin) arrays (rule 1).

### The scored read (per cell + pooled; §4r conventions, midHot per the pinned launch-record definition)

| cell | scored | trap | supp | n̄_det | n₁_solv | n₁/n₂ | W₁_solv | midHot | χ²_med | n≥12 KE |
|---|---|---|---|---|---|---|---|---|---|---|
| s1 | 1862 | 0.069 | 0.192 | 4.189 | 0.2292 | 1.675 | 0.7284 | 1.0355 | 154.8 | 0.0576 |
| s2 | 1873 | 0.064 | 0.185 | 4.187 | 0.2462 | 1.741 | 0.6062 | 1.0226 | 114.2 | 0.0618 |
| s3 | 1858 | 0.071 | 0.188 | 4.054 | 0.2485 | 1.682 | 0.5192 | 1.0212 | 129.2 | 0.0588 |
| s4 | 1869 | 0.066 | 0.192 | 3.917 | 0.2462 | 1.706 | 0.5595 | 1.0109 | 136.3 | 0.0657 |
| s5 | 1868 | 0.066 | 0.180 | 3.991 | 0.2462 | 1.604 | 0.4822 | 0.9805 | 123.0 | 0.0583 |
| **pooled** | **9330** | **0.067** | **0.187** | **4.068** | **0.2433** | **1.680** | **0.5713** | **1.0139** | **242.0** | **0.0603** |

Per-seed SDs (the first seed-scatter measurement on this surface):
n₁_solv ± 0.0079, supp ± 0.0048, trapped ± 0.0030, W₁ ± 0.0954,
midHot ± 0.0207, χ²_med ± 15.4, deep KE ± 0.0033.

### Verdicts (pre-registered BN-P1..P6)

1. **BN-P1 CONFIRMED at launch** (all five cfg oracles; scorer lock
   reproduced the §4z row column-for-column; midHot definition pinned
   — see the launch record).
2. **BN-P2 CONFIRMED — all four bars in-band pooled.** n₁_solv 0.2433
   ∈ [0.230, 0.313] and well inside the frozen window [0.206, 0.411];
   midHot 1.0139 ∈ [0.84, 1.14]; suppressed 0.187 ∈ [0.130, 0.193];
   trapped 0.067 ∈ [0.038, 0.076]. The histogram-level landing is
   seed-robust.
3. **BN-P3 CONFIRMED on the converged branch.** Pooled W₁ 0.5713 ∈
   [0.404, 0.644]. The measured seed SD (± 0.095) retro-diagnoses the
   §4z/s1 "drift" scare: 0.424 (N = 50), 0.524 (N = 500) and 0.728
   (s1) are all within ≈ 2σ of the seed mean 0.579 — the
   "N-resolution tail structure" reading of S6f-P2 was seed scatter,
   not an N-effect. Best estimate of the true model-vs-experiment
   W₁_solv: **≈ 0.57 ± 0.04**.
4. **BN-P4 CONFIRMED — and the shape read is the battery's physics
   headline.** Pooled n ≥ 12 KE 0.0603 ∈ [0.043, 0.069]; χ²_med > 30
   at every seed (114–155) and pooled (242.0 — the I85 σ-shrink
   continues at pooling; no rescue, exactly as registered). The
   deficit **shape is a monotone slope, not a scale**: pooled sim/ref
   mean-KE ratio falls 0.81 → 0.75 → 0.74 → 0.67 → 0.59 → 0.58 →
   0.55 → 0.38 across n = 10–17 (N/bin 301 → 17). Both surviving
   RQ11 candidates are slope-compatible (drag exposure grows with
   dressing; per-shed recoil heating shrinks with terminal n since
   deep survivors shed least) — but uniform-scale stories are dead.
5. **BN-P5 SPLIT (registered primary in, S6 carry marginally out).**
   Pooled Δn̄(MD − twin 4.387) = −0.319 ∈ [−0.72, −0.05]; the tighter
   S6 interpolation carry [−0.63, −0.33] is missed by 0.011 — the
   N = 500 value (−0.384) regresses mildly toward the twin.
6. **BN-P6 CONFIRMED in substance (letter: 23/25 bin-readings at
   0.0000, two at 0.0001).** The 0.1 meV residuals trace to ~10–16
   deep fragments per run (n_det 10–18, ≤ 5 meV each — late escapers
   catching a sliver of in-droplet Landau-gated E2 drag): ≲ 1 % of
   fragments, two orders below the 30–60 meV deficit, and
   energy-*removing* — I89's causal disconnection of the E2 stage
   from the KE curve stands.

### Winner's curse, quantified

Every fresh seed sits below the blessed cell on n₁_solv (0.229–0.249
vs 0.2718, ≈ +1.8σ favorable draw) and above it on suppression
(0.180–0.192 vs 0.161, ≈ −2.2σ) — finc1v725's histogram numbers were
mildly seed-favorable, as regression-to-the-mean predicts for a
selected winner. But the effect is ≈ 2σ, and **every pooled bar stays
in-band**: the standing-result structure (histogram-level landing +
characterized KE miss) survives out-of-sample. The pooled battery is
now the better statistical reference for the standing point.

### Boundaries

One cell family (v7.25 only — no cross-cell contest at N = 5000); the
pooled χ²_med (242.0) carries the same committed error model, so its
absolute scale rides the I85 σ-mechanics; the RQ3 bare-peak
distribution read (≈ 1745 pooled suppressed fragments) is available
in the run dirs but not scored here; the deep-bin slope read thins to
17 fragments at n = 17. Nothing here discharges F5.

---

- **I1 (Wave 1).** In-band (κ, picture, τ) cannot land the staircase: freeze
  at n ≈ 20, max 1.7 sheds. Kinetic, not energetic — the RRK exponent
  (s−1 = 59) on x ≈ 0.032.
- **I2 (Wave 1).** κ is inverted *and* normalisation-capped for 21→14
  (Form-U floors D₀(21) at ≈ 0.53·D₀(1)); κ is a near-dead staircase lever.
- **I3 (Waves 1+3).** Picture is magnitude-degenerate (≤ 4 % spread, from
  x-invariance); its only staircase signal is gate-open timing via Σ(21),
  degenerate with f_int. Discrimination belongs to the size distribution.
- **I4 (Wave 2).** The freeze was the classical dof convention: constant
  s_eff ≈ 8 lands magnitude and timing simultaneously (MAD 1.00 He);
  landing region s_eff ∈ [8, 12] × mid-band τ; s↔τ separate through
  first-shed timing. → s_eff promoted to **Bounded**; F2 Stage 1 re-scoped
  to the s_eff×τ co-fit.
- **I5 (Wave 2).** Effective RRK bath ~8–12 modes, not 60 — consistent with
  quantum mode-freezing / weak coupling in a cold He shell. Gate, ladder,
  RRK form, and Σ(21) crossing all survive.
- **I6 (Wave 4).** The shell-retention limiter is the *ungated* K2 quench
  after ejection, not any asymptote wall. Density-gating K2 (shared bubble
  boundary with drag/pickup) unlocks near-total stripping (n̄ ≈ 2.7) at
  mid-band τ.
- **I7 (Wave 4).** The gated arm is binary in the t×↔ejection race: τ (with
  f_int and Σ(21)) selects between deep strip and a permanently closed gate
  — a structural cliff the campaign must respect; `rho_min` (deferred OQ)
  is the softening lever.
- **I8 (Wave 4).** Bare I⁺ is not reachable in this grid (floor n = 2,
  energetic under the Σ(21)-crossing budget). The 43 % bare experimental
  peak is still unexpressed — open item for F5/production.
- **I9 (program-level).** The probes cost ≈ a third of one campaign stage
  and re-pointed the entire Stage-1 design twice (κ×picture → s_eff×τ; gate
  arm now bracketed). Capability probes before campaign spend paid off
  exactly as intended.
- **I10 (Wave 5).** The s_eff staircase landing is **arm-conditional**:
  ungated s_eff ≈ 8, gated s_eff ≈ 30 (the gate ~4× the in-window shedding
  by removing the post-ejection quench). A gated campaign must re-anchor
  s_eff — the ungated prior does not carry.
- **I11 (Wave 5).** Under the gate, the in-window staircase read (ion-end n)
  and the terminal read (relaxed n) **decouple**: s_eff sets the cascade
  *rate*, not the endpoint (energetic floor n ≈ 2 for all s_eff). At the
  staircase-landing s_eff ≈ 30 the terminal is not converged in 1000 ps
  (frac_frozen = 0, still shedding) → the gated terminal is **flight-time
  dependent**. Ungated has no such split (frac_frozen = 1 everywhere;
  ion-end = terminal).
- **I12 (Waves 4+5).** The two arms are **terminal-regime opposites**
  (ungated = shell-retaining/converged; gated = deep-stripping/rate-limited);
  neither single knob point yields the broad bimodal experimental
  distribution. The observable's shape must come from **ensemble
  heterogeneity across the t×↔ejection race**, which the gate makes
  two-sided — the campaign's real target, not a single point fit.
- **I13 (interpretation, §4c).** s_eff is fundamentally a *rate* (RRK
  exponent), not an amount; it controls terminal n only through the
  evaporation-vs-cooling race (`k·τ`). Remove the race (gate off after
  ejection) and s reverts to a pure clock — same low energetic floor for all
  s, different arrival time. Ungated "freeze high" is energetic (permanent);
  gated "high n at finite time" is kinetic (transient). Cooling, not
  evaporation, is what closes the self-bound gate (evaporation is
  `G`-invariant, self-sustaining).
- **I14 (stance, §4c).** Literal staircase reproduction is a **capability
  showcase, not a calibration anchor** — it shows the mechanism *can* express
  the in-window profile, but the in-window trace is not the observable
  (we don't know the post-window terminal). Only the experimental *terminal*
  size distribution at an explicit detection time carries meaning.
- **I15 (Wave 6).** Under the gate, **f_int is a joint timing +
  effective-budget knob** — "timing-only" is an ungated-only statement. The
  budget at crossing is always Σ(21), but the in-bubble cooling leak over
  [t×, t_eject] scales with how early the gate opens, so f_int moves
  magnitude and timing together; the gated crossing lag also *grows* with
  t× (≈1.3 ps at f_int = 0.35 → ≈2.3 ps at 0.50). Staircase timing
  re-aligns at f_int ≈ 0.40–0.42; landing timing *and* magnitude needs a
  joint (s_eff, f_int) co-fit (interpolated: s_eff ≈ 25, f_int ≈ 0.42).
- **I16 (Wave 6).** The t×↔ejection race, demonstrated from the f_int side:
  the cliff is located in-band (0.50 < f_int* < 0.65 at τ = 6.55; sheds
  nothing above), and the Wave-4 dead τ = 16.5 arm **re-opens** at
  f_int = 0.25 — exposing a second staircase-magnitude region
  (s_eff = 30: n_end 14.65, MAD 1.46, early timing). τ and f_int trade off
  along the race; neither is independently identified by the staircase.
- **I17 (Wave 6).** The floor-vs-leak map: the deepest strip sits where the
  crossing hugs ejection (sampled optimum f_int = 0.50: n̄ ≈ 3.1, min 3);
  *earlier* opening is **shallower** (n̄ ≈ 11 at f_int = 0.24 — the leak
  eats the budget), and later opening is the closed gate. **Bare I⁺ is
  unreached at any f_int** — the OQ-B gap is now quantified from the f_int
  side: within the Σ(21)-crossing construction, no (gate, s_eff, τ, f_int)
  reaches n = 0 at 0.80 eV.
- **I18 (conclusion, user-adjudicated 2026-07-07).** **τ and f_int are not
  independent knobs under the gate — only their race coordinate is
  physically meaningful.** Every gated observable organizes along the
  crossing-vs-ejection margin Δ× = t_eject − t× (with t× = τ·ln(f_int·E/Σ);
  Waves 4+6: the cliff, the re-opening, the leak, the timing). Consequence
  for any future campaign grid: sample the race margin Δ× deliberately
  (f_int derived per τ as f_int = (Σ/E)·e^{t×/τ}), with τ retained as its
  own dimension only for what it independently controls (the in-bubble
  leak/quench strength). Independent τ × f_int grids "make no sense" —
  they sample the race incidentally and unevenly.
- **I19 (Wave 7).** The detected read exists for every probe dir at the
  Sourced t_detect = 8.53 µs (zero new MD; zero P1–P3 guard violations).
  Ungated is **detector-converged**: `n_detect ≡ n_relaxed` exactly (the
  no-op identity, 75/75 dirs). The gated dead arms are **literal detector
  weight at n = 21** (100 % `suppressed`) — OQ-B's fragmentation question is
  now visible in the observable itself. The I14 stance has its observable:
  the terminal distribution at an explicit detection time is now computed,
  not extrapolated.
- **I20 (Wave 7).** **The detector read compresses s_eff.** At f_int = 0.50
  the whole s_eff ∈ [1, 30] range arrives at n̄_detect ∈ [2.7, 4.1] (vs
  2.7–13.8 in-window); s_eff decides the *arrival state* (frozen vs
  `time_exhausted`), the race margin + in-bubble leak decide the arrival *n*.
  The staircase (s-sensitive) and the detector distribution (s-blind,
  Δ×/f_int-sensitive) are near-orthogonal reads — the campaign co-fit gets
  each knob from its own target.
- **I21 (Wave 7).** The gated arm alone spans **n_detect ≈ 3–13 continuously
  along the race margin, plus the suppressed n = 21 class**; the mid-shell
  weight is *converged* (energetic leak freezes, 87–93 % frozen at
  s_eff = 8), not a cap snapshot. Ensemble heterogeneity across Δ× (I12/I18)
  is expressible within a single arm at the detector.
- **I22 (Wave 7).** Bare I⁺ is unreached at the detector (global min
  n_detect = 2); at s_eff ≳ 20 the n ≈ 2 floor is additionally *kinetically*
  asymptotic (≈ 0.4–0.9 He per time decade at s_eff = 30), so neither budget,
  f_int, kinetics, nor flight time expresses the 43 % bare peak within the
  Σ(21)-crossing construction — OQ-B in its sharpest form yet.
- **I23 (post-Wave-7 interpretation, §4e).** The detected read is the
  **Klots evaporative-ensemble regime** (the ensemble sits at k·t ≈ 1 and
  drifts logarithmically) — the mechanism's long-time structure is
  textbook-validated — but the two remaining qualitative gaps are localized
  in two *bookkeeping conventions*, not in any swept knob: (1) the eternal
  cascade traces to the shed draining only D₀(n) (no per-shed translational
  release ε ~ (E_int−D₀)/s ≈ 5–20 meV, first-order over a cascade → OQ-F);
  (2) the missing bare peak plausibly *is* the suppressed class
  (E_int > Σ(n) complexes physically fragment; bare = the
  crossed-after-ejection side of the Δ× race, bimodality = the cliff → the
  OQ-B fragmentation hypothesis). Both are mechanism-convention OQs of the
  same epistemic class the RRK-dof OQ was (resolved by Wave 2).
- **I24 (post-Wave-7 interpretation, §4e continued).** The experiment's two
  largest bins (bare 43.5 %, I⁺He 17.5 %) are **both structurally outside
  the cascade side's reach** (bare needs |G| = 0 exactly; n = 1 needs a
  sub-rung leak *and* sits behind the collapsing RRK bracket). The
  never-opened side explains them **only jointly with OQ-F**: gateless
  boil-off without ε yields exactly bare and zero n = 1 (G > 0 invariant;
  the n = 1 direct channel has no barrier); with ε, the one-sided G₀
  distribution gives bare (bulk) + a decreasing small-n tail (fringe),
  predicting a budget-dependent bare peak. Competing/combinable
  alternative for n = 1: a deep first rung (OQ-G). The small-n abundance
  tail is the first observable that reads the **ladder bottom**; n = 1 is
  the most mechanism-discriminating bin in the distribution.
- **I25 (post-Wave-7 interpretation, §4e continued — OQ2 fired).** The
  f_int parametrization is broken in both directions: the model implements
  **no partition** (E_int(0) = f_int·E_avail is booked on top of the
  untouched Coulomb mechanics — 0.40/1.35 eV with no mechanical source at
  f_int = 0.5), while the physically defensible literal coupling (~1 %,
  8–27 meV) leaves the mechanism inert. The row-14 scenario-keyed floors
  multiply out to the same absolute energy (≈ Σ(21) at both budgets) — the
  natural variable is an **absolute E_int(0) [eV]**, working hypothesis
  0.2–0.5 eV budget-independent from ionization reorganization +
  electronic/spin–orbit relaxation (which would tie the picture knob to
  E_int(0) provenance). Pending literature validation; with OQ-B/F/G this
  completes the set of energy-bookkeeping conventions the probe program
  localized → `RESEARCH_QUESTIONS.md`.

- **I26 (Wave 8).** The gated suppression cliff at a probe knob point is a
  **delta function**: E*_i = Σ(21)·e^(K_tot) = 0.461219 eV with ensemble
  width 4·10⁻¹³ eV, because the probe ensemble is **kinematically
  congruent** (deterministic center-placed onset, one droplet radius →
  identical per-ion cooling exposures K_tot = 0.898297). The suppressed
  fraction is a step 0 → 1: **no single-point ensemble expresses the
  43.5 % bare / 56.5 % shell coexistence at any E_int(0)** — the Δ× race
  heterogeneity the observable needs (I12/I18) is structurally absent from
  the preset and must be *injected* (droplet-radius distribution, thermal
  onset, E_int(0) spread). All per-ion outcome spread at a knob point is
  post-opening channel RNG only.
- **I27 (Wave 8).** The **sub-rung sliver** E₀ ∈ (E* − D₀(1), E*) has a
  leak |G| under one rung, so the opened cascade's energetic floor is
  **n = 1** — expressed in MD (100 % of ions at n_detect = 1 at s_eff = 8;
  2 % even at s_eff = 30 on the deeper-sub-rung point). The experiment's
  second-largest bin is reachable by the flat-bottom ladder without OQ-G
  depth, but only from a window one rung wide in leak units =
  **22.6 meV wide in E_int(0)** (D₀(1)·e^K; ~7.5 % of the sourced band) —
  a strong constraint on any heterogeneity model that wants to feed
  17.5 % of the ensemble through it; I24's "kinetically strangled" holds
  only as |G| → one full rung.
- **I28 (Wave 8).** **The scenario budget is bookkeeping-only** — its single
  physics reader is the S2 onset deposit; Coulomb mechanics, ejection, and
  exposure are budget-blind (K₂₇₀ ≡ K₀₈₀ exactly; MD-verified same absolute
  step at 2.70 eV). In absolute E_int(0) the delivered model is exactly
  budget-invariant: the entire gated map transfers to production verbatim,
  the B.1(2) production prediction is void in both versions, and the
  bare-vs-budget discriminator (RQ3 vs RQ4) has **no in-model mechanism**
  until production Coulomb kinematics are modeled (→ OQ-H).

- **I29 (Wave 9).** The experimental I⁺Heₙ histogram **inverts cleanly** into a
  smooth, unimodal implied p(E₀): a simplex-constrained mixture of the 13
  `s_eff = 8` detected columns reproduces it to **L2 = 0.016, Wasserstein-1 =
  0.086 bins**, the top six bins to < 0.006 (bare fed solely by the
  suppressed→bare column, weight 0.436 ≈ 0.435; near-cliff density on the
  22.6 meV rung window 0.780 %/meV at n = 1, matching the sourced-histogram
  prediction 0.77). Because ions never interact, a run-weighted mixture **is**
  a p(E₀) ensemble — the I26 delta cliff turned constructive. The
  histogram-required E₀ range coincides with the RQ1-sourced [0.2, 0.5] band
  without tuning (strongest quantitative circumstantial support for the
  biphasic mechanism).
- **I30 (Wave 9).** The deep-shell tail forces **≥ 9.6 % of the ensemble below
  the 0.22 eV solvation floor** (columns E₀ ≤ 0.192 eV), plus n = 18–20 (~1 %)
  beyond the s_eff = 8 basis floor (E₀ < 0.08 eV) — W9-P2 confirmed. This is
  the **droplet-radius axis's quantified demand**: the single-radius preset
  cannot supply it from in-band E₀, but a larger-K (larger-droplet) tail
  supplies deep-strip weight at in-band E₀. In-band [0.22, E*] mass 0.467;
  above-E* (bare class) 0.436.
- **I31 (Wave 9).** The implied p(E₀) shows **no multi-modal signature** — mass
  rises monotonically to the cliff and concentrates at/above E* ≈ 0.46 eV, an
  in-band declining tail below it; interior roughness is a basis-collinearity
  artifact (fi0.45/0.48/0.50 co-peak at n ≈ 3–4), not electronic-branching
  structure. Under the W9-P3 decision (fitted p(E₀) = a physics claim on
  E_int(0) provenance), the read is **E_int(0) concentrated at the upper edge
  of / above the RQ1 band**, consistent with reorganization + electronic
  relaxation (RQ1/OQ2) — **not** a demand for fine-structure channels.
- **I32 (Wave 9).** **s_eff becomes detector-identifiable via the n = 1 bin**
  (W9-P4): swapping the near-cliff columns to s_eff = 30 (kinetic wall,
  n_detect 1→2 / 2→3) leaves n = 1 unfed, L2 0.016 → 0.182, W₁ 0.086 → 1.398
  (16×). The small-n bins select the kinetics band (s = 8 floors reached vs
  s = 30 parked) — the first detector-side s_eff handle (Wave 7 found the
  deep-strip detector near s-blind). **Recorded as a possibility only; no
  campaign re-scope** (user decision (c)). Conditional on RQ3 spec (b)
  throughout.
- **I33 (Wave 10, Step 1).** **Width attribution is settled at fixed
  droplet:** every gated column's detected read at fixed (E₀, droplet) is
  ≤ 1 bin wide (n_detect MAD 0.04–0.88 He; per-ion E\*_i spread
  4.15·10⁻¹³ eV, the Wave-8 delta re-derived from trajectories). The
  ~21-bin experimental span must be carried by the E₀ and/or droplet-K
  axes — Poisson + trajectory heterogeneity contribute ≈ nothing.
- **I34 (Wave 10, Step 1.5).** **The gated detected read is two-parameter
  to first order:** the closed-form (E₀, K) fate map — bare ⇔
  K < K\* = ln(E₀/Σ(21)), else no-shed leak + exact ε = 0 descent —
  reproduces the 13 measured s_eff = 8 columns to 0.84 He mean (≈ exact at
  the cliff). The full histogram machinery can be forward-modelled on
  paper before any droplet slice is built.
- **I35 (Wave 10, Step 1.5).** **Probe-K exclusion:** at the pinned 9 Å
  kinematics (K₀ = 0.898), a narrow E₀ ≈ 0.28 eV cannot source the bare
  peak from droplet width — 43.5 % bare needs a log-normal δ ≈ 12,
  ~19× the Kornilov 0.625; every stated prior gives ≤ 2 % bare and
  W₁ ≥ 2.7 bins. The corrected-direction W10-P3 is confirmed at this
  scale by an order of magnitude, not a margin.
- **I36 (Wave 10, Step 1.5).** **Production-K landing:** at the
  speed-scaling bracket K₀ ≈ 0.49 the suppression cliff sits mid-ensemble
  and (E₀ = 0.28 eV, δ = 0.80) lands **bare 42.6 %, W₁ = 0.993 bins
  untuned**, with a monotone small-n envelope; residual n = 1 deficit
  1.29× ≤ the F.5 taper cap, mid-tail (n = 2–4) overshoot open. W10-P5's
  split verdict holds → **"narrow E₀ + RQ7 kinematics" is the live
  resolution and RQ7 is co-requisite for Step-2 adjudication.**
  Conditional on RQ3 sequential-shed and ε ≈ 0 throughout; analytic
  pre-verdict, Step 2 arbitrates.
- **I37 (Wave 11).** **K₂.₇₀ = 0.746 measured** (first MD off 9 Å
  geometry; suppressed shed-free ride at `R0_GS = 2.666 Å`, congruent to
  10⁻¹³, E_int cross-read agrees to 0.1 %). The ballistic ×0.545 bracket
  is **wrong**: drag eats almost all the extra channel speed (peak
  10.5 Å/ps decelerated to 3.8 Å/ps by exit), S_K = 0.831. **Validity:
  41.7 % of the exposure accrues beyond the calibrated speed band**
  (overlap share 8 %); both biases push K up → 0.746 is the upper edge,
  0.49 the hard lower bracket. W11-P2 refuted as registered.
- **I38 (Wave 11).** **The (E₀, K) closed form transfers to production
  kinematics unchanged:** opened companion MD n_detect = 5.76 (5–8) vs
  closed-form n = 6 (0.24 He); the m(t)-feedback of an opened cascade on
  K is −2.2 % — the fate map needs no geometry-specific correction, and
  the K(R) machinery wiring-oracles against the Wave-10 table to ≤ 0.007.
- **I39 (Wave 11).** **W11-P5 refuted as registered, landing
  re-calibrated in-band:** at measured K the pre-registered
  E₀ ∈ [0.24, 0.32] cannot source the bare peak (≤ 18 % everywhere), but
  every stated prior crosses 43.5 % at E₀ ≈ 0.37–0.43 — inside RQ1's
  [0.2, 0.5] — with untuned best W₁ = 0.496 (E₀ = 0.41, δ = 0.80),
  better than Wave 10's 0.993. The bare-crossing E₀ tracks
  E\*(K₀) = Σ(21)·e^(K₀) ≈ 0.30–0.41 eV over the whole K bracket → the
  K uncertainty moves the scalar calibration, not the verdict; the
  K-scale itself is near-degenerate in W₁ (probe curve lands at 0.570)
  — settled by measurement, not fittable from the histogram. The n = 1
  deficit worsens to 1.5–2.9× (δ-dependent trade-off against the n₂₋₄
  overshoot) — beyond the F.5 taper cap except at δ = 0.40 → the bin
  moves to RQ4 / the E₀-smear leg.
- **I40 (Wave 11).** **First VMI-side kinematics constraint:** detected
  model speeds at production kinematics (4.1–5.3 Å/ps) undershoot the
  experimental I⁺He peak (10.1 Å/ps) by ~2–2.5×, with ballistic at
  15.7–20.3 — the locked drag over-dissipates production fragments, in
  the same direction as the beyond-band over-drag reading (true
  K₂.₇₀ < 0.746). Candidate research item: a drag law validated at
  production speeds (TDDFT at 2.666 Å kinematics / velocity-capped
  form); existence-level, single knob point.
- **I41 (H.2b).** **The all-bounded lever set cannot land the solvated
  targets** (7 308-cell scan, zero T1∧T2∧T3 passes, both exposure
  brackets — not law-conditional). The miss is *localized*, not diffuse:
  bins 2–8 and the KE envelope/shape land; the residuals are exactly the
  n₁ bin (histogram side) and the KE scale (speed side). Outcome (c) at
  the strict bar, decomposing into the two already-named items (I42,
  I43).
- **I42 (H.2b).** **The n₁ bin is ladder-taper-controlled, and the
  slid-X₂ family is rejected as its repair:** n = 1 is a one-rung-wide
  E_ej window, so equal deep rungs (slid2/3, ratio 1.25/1.29) revert
  n₁/n₂ to the flat chord-K density ratio, while the knob-free X₂ floor
  reaches 1.83 (cap ≈ 1.44×) and the **RQ4-graded diagnostic
  (2.2 : 1.5 : 1.3) nearly lands** (T2+T3 pass, W₁ 0.575 vs the 0.5 bar,
  bins 2–8 to ≲ 0.02). Through a purely geometric forward model the
  solvated histogram **independently demands the RQ4 target ratios** —
  the external many-body calculation is the blocking arbiter; ε/missing-
  mechanism fires only if it returns a plateau.
- **I43 (H.2b).** **Geometry owns the KE shape; the scale is v_c's:**
  the experimental mean-KE curve sits inside the [current-law,
  ballistic] envelope at every top cell, and the required lift is
  speed-selective (×2.2 at n = 1 → ×1.7 at n = 12 in energy; ×1.5 → ×1.3
  in speed) — exactly the velocity-cap direction. W13-P1's premise
  (v_c owns scale, not histogram shape) is pre-confirmed with numbers.
  *[Superseded in shape by §4k (2026-07-15): the D7 targets were the
  legacy moment convention. Against the corrected reference the envelope
  conclusion stands (re-confirmed at f = 1), but the lift is ×2.6–4.0 at
  n = 1 and near-uniform ×1.6–2.2 for n = 2…17 — the selectivity lives in
  the n = 1 bin alone.]*
- **I44 (H.2b).** **A trapped, droplet-retained ion class exists under
  the current law** (new structural class): 6–11 % of fragments — inward
  partners of off-center births, chord exposures up to K ≈ 18 — are
  dissipated below the 0.117 eV solvation barrier and never eject
  (quasi-static park at the far surface). Absent in the ballistic
  bracket; invisible to the I⁺Heₙ bins; t_end-conditional; MD-arbitrable
  by the W12 leg. Every prior wave's exposure map (droplet-R axis, K(R))
  implicitly assumed ejection — position heterogeneity breaks that.
- **I45 (H.2b).** **E₀ re-lands at the solvation scale under full
  geometry** (optima at 0.22–0.27 eV ≈ E_solv vs the pinned-droplet
  0.38–0.41 of I39): the position axis supplies the low-K mass the E₀
  scan previously bought with a higher cliff, and bare is un-targeted
  under the two-channel reading. Also confirmed with the margin floor:
  no gate-open-at-birth class exists in-bounds (n_eject ≥ 13,
  Σ(n_eject) ≥ 0.12 eV) — the W12b-P1 fast surface feeder is **predicted
  absent** under both bounded E₀-law brackets p ∈ {0, 1}.
- **I46 (Addendum I, Step 0/1).** **The corrected reference converts the
  KE miss from a scale problem into a slope problem:** the required lift
  is near-uniform ×1.6–2.2 for n = 2…17 (D7's "speed-selective ×1.3–1.5"
  was the legacy-convention bias), the tail family delivers that scale at
  n ≤ 6, but the experimental ⟨E⟩(n) decays far more slowly than any
  model curve — the deep bins (n ≥ 7) undershoot progressively to
  ×0.26–0.49 by n = 12–17. Deep bins select the most-exposed chords whose
  late slow dynamics ride the TDDFT-locked **in-band** cubic law — no
  admissible high-v tail touches them.
- **I47 (Addendum I, Step 1).** **KE and histogram anti-correlate through
  the single exposure integral K** under the frozen fate map
  (E_ej = E₀·e^(−K), ε = 0) and the pinned cooling clock (τ = 6.55 ps):
  every KE-lifting tail cuts K_q50 from 0.75 to 0.18–0.37 and starves the
  evaporative descent (ke_pass ≥ 8 ⇒ W₁ ≥ 1.31; W₁ ≤ 0.9 ⇒ ke_pass ≤ 7).
  The drag lever alone cannot serve both observables — decoupling needs a
  second knob (τ, or a fate-map-level change; → OQ-I).
- **I48 (Addendum I, Step 1).** **The tail family is well-behaved and the
  hard cutoff is excluded:** center-pin v_inf strictly monotone in v_c
  and strictly ordered across tails (I-P2 confirmed) — the Step-1 failure
  is structural, not a search artifact; the cut bracket produces
  99.7 % suppression → bare (W₁ ≈ 3.7) at v_c ≤ 10 and degenerates to the
  current law above the peak speed.
- **I49 (Addendum I, Step 1b).** **Joint closure exists: (v_c, τ)
  together do what neither knob does alone.** 18 contiguous cells at
  τ = 4.1 ps × (p = 0, v_c 6–7 / p = −1, v_c 7–8.5) on the rq4graded
  ladder pass the full ×1.25 KE bar (12/12) with W₁ ≤ 0.85; the best
  cell (p = −1, v_c = 7, pickup, m3, E₀ = 0.25) is the **first full
  two-observable landing of the program** (T1 ∧ T2 ∧ T3 ∧ KE bar;
  W₁ = 0.272, n₁ = 0.280, ratio 1.861, worst KE ×1.24). The I47
  anti-correlation resolves by the pair: v_c owns the speed integral, τ
  re-owns the descent clock; the deep-bin slope was K-selection, not
  in-band physics. Conditional on the RQ4 taper (stakes sharpened: a
  plateau now breaks both observables) and on re-classifying the GAH25
  τ pin (4.1 = 0.63 × 6.55 — fires as pre-registered in §I.8).
- **I50 (Addendum I, Step 1b).** **The two surviving tails are
  experimentally distinguishable:** p = −1 (constant-force tail) wins
  the histogram ratio (T2/T3) but sags in the far KE tail (n ≥ 13
  ratios 0.59–0.78); p = 0 (linear-force tail) yields the flattest KE
  shape yet (×0.85–1.04 across all n = 1…17) but tops at ratio ≈ 1.66.
  Discriminators for the MD build / experiment: the n ≥ 13 mean-KE tail
  and n₁/n₂.
- **I51 (Addendum I, Step 1c).** **The closure is a plateau, and the
  taper's role sharpens to exactly the last stretch:** the joint basin
  spans ≈ 1.5 Å/ps × 1.6 ps per form (225 cells; 24 full houses; best
  W₁ = 0.196 at (p = −1, v_c = 7.5, τ = 3.8)); p = 1 is excluded
  outright. K-P3 refuted: the knob-free floor1 ladder joint-closes too
  (12/12 KE + W₁ 0.58) — **in-bounds physics + (v_c, τ) now reproduces
  the full KE curve and an H.2b-best-level histogram**; the RQ4 taper
  specifically buys n₁ = 0.26+, ratio ≥ 1.75, and W₁ 0.58 → 0.20. RQ4's
  arbitration is therefore about the histogram's small-n shape only —
  the KE side stands either way.
- **I52 (Slice T3, §4l).** **The undressed MD cannot express the
  closure's small-n structure at any knob value:** with every ion
  starting at the full n₀ = 21 shell, all four §I.10 configs park in a
  narrow n ≈ 7–9 cluster (zero weight below n = 5; bare and n₁
  identically zero). The 1D↔MD gap is the *ensemble geometry* (L2
  dressing + E₀ law + margin + prior), not the calibrated (v_c, τ,
  ladder) values.
- **I53 (Slice T3, §4l).** **The mechanism transfer itself is
  confirmed at production kinematics:** τ/E₀ ordering across C1–C4
  follows the F.2b race/leak closed form; the capped tails lift
  detected KE (×1.3–2.0 vs current-law ×0.85–1.15 in populated bins);
  90–95 % frozen arrivals (converged reads); dt-halving drifts
  ≤ 0.30× SEM.
- **I54 (Slice T3, §4l).** **`ANCHOR_N_START` = 21 is a
  validation-era convention doing production work:** the biphasic seed
  hard-codes the Tier-1a 9 Å full-shell start regardless of birth
  position — physically indefensible for surface births at 2.666 Å
  kinematics, and exactly the H.3b "n_eject(depth)" surface that was
  designed and parked.
- **I55 (leg A′, §4m).** **The position axis alone is a two-sided race
  lever, and the twin transfers quantitatively on the class axis:**
  undressed MD with only `uniform_volume` births (margin 3 Å)
  reproduces the pre-registered twin histograms at W₁ = 0.55–0.69 bins,
  with the suppressed/bare ordering across configs transferring exactly
  (0.315/0.116/0.095/0.032 vs twin 0.336/0.118/0.089/0.035). Corollary:
  the T3 park (I52) was center-pinning as much as missing dressing —
  the first MD weight at n = 0 and n = 1 needed no shell dressing at
  all.
- **I56 (leg A′, §4m).** **The droplet-retained class is real in MD
  (5–11 % per config) and includes centrifugal resonances** — marginal
  ions above the radial escape threshold but below their
  effective-potential barrier, effectively permanent under the
  conservative (zero-gamma) E2 relaxation. Handled by the V0-2
  `exclude` policy with the exact conservative-mechanics bound
  criterion; whether E2 should carry drag (which would physically
  capture such ions) is an **open question for user adjudication**
  before any N = 500 run.
- **I57 (leg A′, §4m).** **Histogram parity is not composition
  parity:** MD small-n bins are fed by near-ballistic fragments (n₁
  mean KE 2.42–2.48 eV vs the twin's 0.97; c3 0.67 vs 0.26) — the
  fate-map↔real-cascade channel (S2c-P4, listed) decouples a bin's
  weight from its occupants. The (n, mean-KE) curve is
  composition-sensitive, so KE conclusions require the full dressed
  geometry (leg B onward), not histogram agreement alone.
- **I58 (I.11.2 item 1, §4n).** **The pre-registered mass-asymmetric
  Coulomb-split hypothesis is refuted as I57's driver:** the exact
  per-fragment Coulomb-work integral gives a share excess of only
  +0.08 eV (≈ 3 %) for n₁ ions — they exit the droplet full-shell at
  0.7–0.9 ps, mass-symmetric during the acceleration, and shed only
  after ejection. There is no light-early route in the data.
- **I59 (I.11.2 item 1, §4n).** **I57's mechanism is the cold-shed
  momentum convention:** the delivered evaporation channel leaves each
  shed He at rest in the lab frame (`v → v·m/m′` per shed, verified
  event-by-event), injecting +0.91 eV over a full post-exit strip
  (×1.611 = m₂₁/m₁ in KE). The co-moving counterfactual
  `KE_cf = ½·m_det·(v_det·m_det/m_exit)²` reproduces the twin per
  config to ≤ 2 % (0.955/0.952/0.256/0.933 vs 0.97/—/0.261/—) — the
  twin *is* the co-moving convention, and the twin−MD KE divergence is
  entirely convention, zero residual mystery. → **OQ-J**; couples to
  RQ2 (ε) as one "what does the evaporated He carry away" discussion;
  the §4k (v_c, τ) KE closure is co-moving-based and does not transfer
  to a cold-shed MD at small n.
- **I60 (I.11.2 item 1, §4n).** **Composition map at production
  kinematics:** n₁ and bare occupants are near-edge outward births
  (r̄_birth 20–23 Å, μ̄ +0.6–0.9, drag loss 0.8–1.2 eV vs 2.6–2.8 eV
  for deep bins) — the drag deficit is the composition axis; the
  convention (I59) is the speed axis. The suppressed/bare detected
  read is convention-free (no sheds), but its RQ3 *fragmentation* read
  is convention-decided (co-moving ≈ 1.18 eV vs momentum-conserving
  bare ≈ 3.3 eV at c1 speeds) — the experimental bare-bin mean KE
  discriminates the fragmentation convention.
- **I61 (leg A″, §4o).** **The co-moving convention closes I57 in real
  MD and makes twin↔MD parity two-axis:** with only the shed convention
  flipped vs A′, the whole solvated per-bin KE curve lands on the twin
  (n₁ ratio ×2.5 → ×1.03–1.20; W₁(A″, A′) = 0.14–0.20 bins — the
  histogram is convention-blind as claimed; suppressed ordering exact).
  Two refinements: no ion is *strictly* convention-isolated (pair
  Coulomb coupling moves never-shed ions by ≤ 20 meV), and the marginal
  near-barrier class is fragile at meV scale (1–3 ions/config flip
  retained; c3 11 → 14) — raising the E2-dissipation adjudication's
  stakes for N = 500. The leg-B baseline is the A″ dirs, on the same
  convention basis as the twin by construction.
- **I62 (leg B, §4p).** **The T5 dressing transfers quantitatively
  across the 1D→3D boundary:** one lever (density_tied) triples the
  suppressed/bare class along the twin's prediction (ordering exact),
  truncates the deep tail at the dressed band, and lands
  W₁(MD, twin) = 0.41–0.51 bins — the chain's best histogram agreement
  — while moving the histogram itself by 2–4× that (W₁ vs A″ 1.0–1.9).
  The dressing, not the position axis alone, is what feeds the
  bare-candidate class at scale.
- **I63 (leg B, §4p).** **Pickup re-filling is real and measurable:**
  the MD suppressed fraction sits uniformly 0.06–0.08 below the twin
  (and n̄_det 0.2–0.3 above) — the live Langmuir channel re-dresses
  under-dressed shells in-bubble before gate-open, a channel the 1D
  twin lacks. First direct quantification of the dressing↔pickup
  interplay (the H.3b boundary), ≈ 7 ions/100 at these knob values.
- **I64 (leg B, §4p).** **The bare-bin KE gap is bookkeeping, not
  dynamics:** MD reports the suppressed class as an intact dressed
  complex (m̄ ≈ 190 amu); the twin books the RQ3 spec-b bare fragment.
  The co-moving break-up rescale m_I/m_complex collapses the ×1.6 gap
  to ×1.08–1.12 on every config. The experimental bare-bin mean KE
  therefore reads directly on the fragmentation convention (sharpens
  I60); the solvated KE curve needs no such caveat — it lands on the
  twin at ×1.2.
- **I65 (leg C, §4q).** **The T6 `sigma_proportional` p-law transfers
  quantitatively and de-suppresses to the twin.** Flipping p = 0 → 1 on the
  certified `bc` baseline drops the suppressed/bare class to ≈ a third
  (0.322/0.267/0.464/0.378 → 0.101/0.022/0.301/0.124), landing on the twin
  p = 1 (0.111/0.051/0.343/0.143) with ordering c3 > c4 > c1 > c2 exact and
  **W₁(C, twin) = 0.31–0.47 bins — the best twin↔MD agreement of the whole
  oracle chain**. The §4j "p = 0 over-suppresses" verdict is now demonstrated
  in real MD: the under-dressed-birth onset regularised to Σ(n₀)/Σ(n*) is the
  physical direction. Pickup re-filling puts MD ≈ 1–3 ions/100 **below** the
  twin (smaller than leg B's ≈ 7/100 — less suppression to re-fill).
- **I66 (leg C, §4q).** **p = 1 moves the solvated n₁ KE toward experiment.**
  The de-suppressed n₁ occupants ride shallower descents, so n₁ mean KE rises
  1.08/1.10/0.32/1.09 eV (from `bc` 0.66/0.63/0.26/0.60) — ≈ the A″ undressed
  and the **closest of any leg to the committed experimental n₁ = 1.302 eV**
  on the solvated branch; the whole solvated KE curve tracks the twin at
  ×1.08–1.15. Trapped is p-invariant (0.11/0.11/0.17/0.11 ≈ `bc`), confirming
  the p-law touches only the onset, never the chord dynamics.
- **I67 (leg C, §4q — method).** **The E2 relaxation checkpoint is
  scored-read-neutral under stride.** Detection seeds from the terminal column
  only and the stage always stores the true final state last, so shrinking the
  checkpoint byte budget (here ~380 MB → ~2 MB) leaves every terminal/detected
  read bit-exact while only coarsening the unused intermediate trajectory — a
  reusable accommodation when the full-trajectory compressed save exceeds an
  execution window (used for the `cc` dirs; the `bc` dirs keep the full
  trajectory).
- **I68 (leg D, §4r).** **The droplet and KE axes are decoupled.** Flipping
  the T8 kornilov prior on the certified `cc` baseline moves the histogram
  (n̄ down ≈ 1 He, suppression re-ordered within the re-filling band,
  ordering c3 > c4 > c1 > c2 exact) while the solvated n₁ KE moves ≤ 0.12 eV
  (c3 unchanged at the printed precision) — exactly the twin's
  prior-invariance claim (DP-P3). The solvated per-n KE curve is a pure
  drag-law/onset observable; the droplet prior buys histogram freedom
  without re-opening the KE calibration.
- **I69 (leg D, §4r).** **The chain-best W₁ trend breaks at the droplet
  axis, and the miss is the twin's frozen birth geometry.** W₁(D, twin) =
  0.54–0.73 (vs leg C's 0.31–0.47) from a uniform 0.43–0.58 He MD-below-twin
  softening (≈ 1 SE per config, same-signed on all four): the MD's well
  depth, gate depth, and pickup exposure follow R_i through the live
  cascade while the twin's chord is frozen at birth geometry — channel (d),
  listed since the trigger, now measured. On the droplet axis the twin's
  authority is ordering/direction, not magnitude: the basin-locator stance
  (post-leg-B) is now a measurement, and the re-pilot re-centering must be
  MD-driven.
- **I70 (leg D, §4r — method).** **The I67 accommodation composes with
  per-config parallelism.** Four per-config resume drivers (shrunk
  relaxation checkpoints) ran the whole leg in ≈ 14 min wall-clock / ≈ 1
  CPU-hour, with the scorer oracle-locked against the recorded §4q numbers
  before any leg-D data was read — the leg-execution pattern for the
  re-pilot matrix.
- **I71 (re-pilot Stage 1, §4s).** **The KE axis pins v_c only
  one-sidedly.** The MD profiled-χ² argmin lands at the twin's cell + one
  step toward higher v_c on every arm (S1-P3), but sits at the bracket
  edge with χ² still falling on all three arms — the whole-curve v_c\* is
  right-censored at v_c = 8.5 (p = −1) / 7.5 (p = 0). The twin's
  ordering authority held; its one-step bias direction was exactly right.
- **I72 (re-pilot Stage 1, §4s).** **The speed-selective tension survives
  the real cascade — a drag-form-shape finding.** Within capped_cubic no
  v_c lands the whole KE curve and the n₁ anchor/solvated ratio
  simultaneously: χ² pulls to v_c ≥ 8.5 while n₁ = 1.302 eV extrapolates
  to v_c ≈ 6.2–6.3 and the v65 cells land n₁/n₂ = 2.14/2.08 (ref 2.18)
  with n₁ ≈ 0.30 (the H.2b target, first time in MD). Three observables
  (KE χ², n₁ anchor, W₁_solv) prefer three different v_c. The taper
  *shape* — not the calibration — is what's short (I45 measured in MD).
- **I73 (re-pilot Stage 1, §4s).** **Channel (d) scales with drag
  exposure.** Δn̄(MD − twin) runs monotonically 0 → −0.6…−0.8 He across
  each bracket (v65 → v85), and W₁(MD, twin) runs 0.21 → 1.04: the leg-D
  "uniform 0.5 He softening" (I69) was a single-v_c snapshot. At the
  low-drag end the frozen-geometry twin is nearly exact (W₁ 0.21–0.26,
  the chain's best parity anywhere); the bias is a *function of
  in-droplet residence time*, not a constant offset — pre-registered
  check S1-P5's non-uniformity arm, fired as designed.
- **I74 (re-pilot Stage 2, §4t).** **The twin is a quantitative
  histogram predictor at low drag across the whole budget plane.** At
  v65, Δn̄(MD − twin) sits within ±0.17 He on all 32 (τ, E₀) cells with
  no budget trend; at v85 the softening is a uniform −0.72 He. Channel
  (d) is a pure drag-exposure function — (τ, E₀)-blind. Consequence:
  low-drag twin sweeps can be trusted near-quantitatively for histogram
  design work; the KE axis cannot be factorized the same way (the
  fate-map composition moves n₁ KE by ~0.5 eV across the grid at fixed
  v_c — S2s-P5, predicted by the twin and confirmed in MD).
- **I75 (re-pilot Stage 2, §4t).** **Each experimental axis is landed —
  at different cells — and the joint landing does not exist inside the
  swept capped-cubic family.** KE: (v85, τ3.2, E₀0.23) hits the n₁
  anchor exactly (1.310 vs 1.302 eV; profiled χ² 13.9, grid-best band)
  with a badly shell-retaining histogram (n₁ 0.10). Histogram: (v85,
  τ3.8, E₀0.27) and (v65, τ3.2, E₀0.27) land n₁_solv 0.320/0.326 ≈ the
  H.2b 0.31 target (first MD cells ever) with the KE axis missed (0.52
  eV / χ² 162). The registered candidate region behaved exactly as the
  twin said (the anchor is crossed inside it), the Stage-1 χ²
  right-censoring resolves *interior* on the budget axes (extension (a)
  moot) — and still no (v_c, τ, E₀) point lands both. I72 sharpened to
  a three-knob, cell-resolved statement about the taper *shape*.
- **I76 (re-pilot Stage 2, §4t).** **rq4graded > floor1 on the histogram
  axis at matched knobs** (top-4 W₁_solv all c1; c4 best 0.793 vs c1
  0.594) — the I51/S2c-P3 direction expressed across the full grid,
  still N = 50-soft; the formal ladder verdict remains the N = 500
  two-finalist read (RP-D6).
- **I77 (n=1-deweight re-score, §4u).** **The experimental n₁ mean-KE
  anchor (1.302 eV) is a two-population mixture mean, not the solvated
  core.** The ihe_ked reference self-flags it: mode 0.891 / median 1.128 /
  mean 1.302 (mean/mode 1.46), and n = 1 is the sole low-n fragment with
  `dominantError = bg-structural` (bgOffShift 0.378 eV, 29 %). The MD n₁
  (≈ 1.20) has no Coulomb-fast tail (§4n) and matches the ref **median**
  (×1.06), not the mean (×0.92) — scoring MD-core against ref-mixture-mean
  is not like-for-like. De-weighting n = 1 removes **36–44 %** of the
  low-drag KE χ² (oracle: baseline reproduces §4s/§4t exactly), but the
  χ² argmin stays at high drag under every treatment (drop / median /
  median+bg): **the high-drag pull is owned by the mid-bins n2–n8, not by
  n=1.** Adopt median-anchored (or n=1-excluded) KE scoring — a
  like-for-like correctness fix (n=0 already excluded, I-D4) that isolates
  the residual; it does not resolve the two-axis tension.
- **I78 (E-2 taper-corner twin scan, §4u).** **The E-2 (v_c, p_tail)
  corner has no joint cell — p_tail is the wrong second lever.** A
  more-negative tail exponent re-heats the whole KE curve (all fragments
  cross the above-cap fast phase during ejection) *and* re-strips, so KE
  and n₁_solv move together, not orthogonally. But **raising v_c alone
  lands the entire mid-bin KE curve** (v6.5 → v7.5, p−1: midHot 1.74 →
  0.95, n2–n8 within ~10 %, n₁ at the core) — the mid-bin bulge is a real,
  fixable mid-band drag deficit (I72 constructive). The v7.5 residual is
  histogram over-retention (n₁_solv 0.238 vs 0.310), a
  stripping-at-fixed-KE deficit: KE and stripping are coupled through the
  cooling exposure K (Addendum-I anti-correlation re-confirmed), which the
  taper cannot break. **The decoupling lever is the ladder bottom (RQ4),
  not the taper** — v_c ↔ KE curve, ladder ↔ n₁_solv is the orthogonal
  factorization to test next (twin v_c ≈ 7.5 × ladder family). The twin
  *under-strips* at high v_c (I73/I74), so MD n₁_solv at v7.5 sits above
  0.238 toward 0.31 — the **un-scored (v7.5, τ3.2, E₀0.27) MD cell** is the
  twin-identified best joint candidate (Stage-1 s1c1v75 already gave MD
  n₁_solv 0.250, χ² halved to 54).
- **I79 (v_c × ladder twin scan, §4u NB).** **The v_c↔KE / ladder↔n₁_solv
  factorization is orthogonal but the ladder is already saturated at
  rq4graded — the histogram recovery at v7.5 is a twin→MD strip bias, not a
  ladder move.** `midHot` is flat across the ladder family at every v_c
  (v7.5: 0.95–0.97) — the ladder is KE-neutral (orthogonality confirmed).
  But rq4graded gives the **highest** n₁_solv (0.238) and best ratio (1.87);
  cheaper-bottom ladders (flat/floor1/slid2/slid3) strip *past* n = 1 into
  bare (0.20 → 0.34), *lowering* n₁_solv — I78's cheap-bottom decoupler is
  **refuted** (rq4graded's 20 meV D0(1) barrier is what piles weight at
  n = 1). No twin cell reaches n₁_solv 0.31 at midHot ≈ 1. The gap closes
  via the twin's under-strip bias (measured at v6.5: twin 0.263 vs MD 0.326,
  +0.063 at matched n̄), larger at v7.5 → MD n₁_solv ≈ 0.30 expected. The
  decisive test is the **(v7.5, τ3.2, E₀0.27, rq4graded/c1) MD cell**; a
  more-bottom-heavy ladder is RQ4-external (parked).
- **I80 (v7.5 MD spot-check, §4v).** **The joint landing exists — the
  twin-identified (v7.5, τ3.2, E₀0.27) cell lands the histogram AND the KE
  curve simultaneously, the first MD cell to do both.** s2c1v75t32e27
  (c1 leg-D, only v_c 6.5 → 7.5 vs the on-disk v65 cell; cfg byte-oracle
  passed): **n₁_solv 0.291** (target 0.310), **n₁/n₂ 2.30** (2.18),
  **W₁_solv 0.496** (the re-pilot's best), **midHot 0.904** (n2–n7 within
  ~10 %), n₁ KE 0.915 (at the core), like-for-like **χ²_med 14.4** ≈ the
  grid-best 13.9 — but with the histogram *landed*, not sacrificed. §4t's
  "no joint cell in capped_cubic" was a **sparse-v_c-grid artifact** (only
  6.5/8.5 sampled) plus the n₁ anchor inflation (I77); the joint optimum
  sits at the un-sampled v7.5, exactly where the §4u twin + strip-bias
  (I78/I79) pointed — the twin predicted 0.238/0.953/0.906, MD delivered
  0.291/0.904/0.915 (the +0.053 lift is the predicted twin→MD strip bias).
  N = 50 pilot; the formal verdict is the N = 500 finalist (RP-D6).
- **I81 (v_c sensitivity ring, §4w).** **The joint landing is a basin
  v_c ∈ [7.25, 7.5], not a knife-edge; the optimum refines to v_c ≈ 7.25 on
  the histogram, exactly as pre-registered (S3r-P2/P3).** The 7-point curve
  (v_c 6.5→8.5 at c1/τ3.2/E₀0.27) shows W₁_solv a clean bowl minimizing at
  7.25 (0.424; n₁_solv 0.309 ≈ target 0.31, n₁/n₂ 2.27 ≈ 2.18, midHot 0.987),
  midHot sliding monotonically through ~1.0 at ≈ 7.2, and n₁_solv crossing
  0.31 at ≈ 7.25. The N=50 **χ²_med is thin-bin-noisy** (off-trend spikes at
  v7.25=66 and v8.0=139 that per-bin inspection confirms are 2–4-ion deep-bin
  scatter, not KE misses) — so the strict χ²≤30 gate under-counts the window;
  read on the outlier-robust W₁/midHot the basin is clear. The basin closes
  above ~7.75 (v8.0 genuinely over-dragged, whole curve ×0.68–0.88). c1 N=500
  finalist: center v_c ≈ 7.25–7.5, χ² re-adjudicates once N=500 stabilises
  the deep bins.
- **I82 (S4 c4 location, §4x).** **The c4/floor1 joint region exists in
  MD — at (v7.0–7.25, τ3.8, E₀0.25), not the c1-like τ3.2.** MD lands
  n₁_solv 0.312/0.303 (target 0.310) with midHot 1.121/0.975 — both
  robust joint bars pass — while W₁_solv stays 0.630–0.651 vs c1's
  0.424–0.496 and the ratio overshoots (3.0–3.3 vs 2.18, n₂
  under-filled). The S4-located c4 finalist is **(v7.0, τ3.8, E₀0.25)**
  (frozen tie-break); the RP-D6 ladder contest at N = 500 is thereby
  real — c4 enters at its own basin-located best and c1 still leads on
  the histogram shape. Registered-expectation miss recorded: S4c-P3
  predicted no c4 window (outcome (b)); candidates existed at the longer
  τ. Sub-finding: c1's suppressed channel closes entirely at the
  E₀ = 0.23 corner (supp ≈ 0), inverting the c4 < c1 n₁_solv ordering
  there — a budget-corner regime change.
- **I83 (S4, channel-(d) family).** **The twin→MD under-strip lift is
  budget- and arm-dependent, not universal:** +0.100/+0.097 on the
  floor1 E₀0.25 cells (≈ 2× the c1-measured +0.05..0.065 band) but
  +0.019/−0.022 at the high-suppression E₀0.27 corner (supp ≈ 0.43) —
  the lift collapses when the suppressed channel is heavily populated.
  Bias bands remain *search guidance only* (the I69 discipline); any
  quantitative twin claim near a high-suppression corner must be
  MD-confirmed cell-by-cell.
- **I84 (S5 Landau bracket, §4y).** **The `v_L` gate dissolves:** every
  scored observable is bit-flat across v_L ∈ {0.30, 0.40, 0.58}
  (Δn₁_solv = ΔW₁ = 0.0000, ΔmidHot = −0.0002, retained 8/100 unchanged,
  χ²_med 14.4 on all arms) while the arm demonstrably acts on the
  retained class (E_dissip +0.58–0.60 eV; retained KE 0.277 →
  0.005–0.023 eV, monotone in v_L) — quiet because the scored surface is
  structurally insulated, not because the arm is inert. N = 500 runs
  Landau-on at v_L = 0.58 with 0.40 as the winner spot; the external
  re-pinning is non-blocking. Execution finding: `v_limit_m_per_s` is a
  latent shared reader (it also sets the neutral-stage collision
  threshold `E_min_eV`) — measured quiet here by neutral/ion sha256
  identity, to be re-checked at any lower-energy channel.
- **I85 (S6 no-landing + the sigma unmasking, §4z).** **No N = 500 cell
  lands the frozen joint acceptance** — every cell fails χ²_med ≤ 30
  (125.7 / 68.6 / 202.5) and only finc1v725 holds both robust bars
  (n₁_solv 0.2718, midHot 0.9905). S6f-P3 is refuted with its mechanism:
  the N = 50 χ² was sim-SE-dominated (σ 0.03–0.08 eV vs the reference's
  0.010–0.03), so N = 500 does not resolve a spike — it unmasks a
  systematic KE-curve miss. The winner-gated riders did not fire; the
  F5 reconciliation is the next adjudication.
- **I86 (S6 ladder verdict, §4z).** **The RP-D6 contest resolves for
  c1/rq4graded at N = 500 on every histogram read:** W₁_solv 0.524/0.657
  vs c4's 0.843, χ²_med 125.7/68.6 vs 202.5, n₁_solv 0.272/0.255 vs
  0.232 (≈ 2σ; the I76 direction resurfaces after the S4 N = 50
  near-tie). The I51/S2c-P3 bounded-physics claim is formally read:
  rq4graded > floor1 — c4's KE curve also fails structurally
  (whole-curve rotation, hot shallow / cold deep).
- **I87 (S6 deep-bin KE miss, §4z).** **The model's deeply-solvated
  survivors arrive too cold:** sim mean KE at n ≥ 12 sits 30–60 % below
  the reference on the c1 cells (n ≥ 12 carries ≈ 72 % / 54 % of
  χ²_med at v7.25 / v7.5), and W₁_solv worsens same-direction at all
  three cells vs N = 50 — an N-resolution tail structure, the sharpest
  experimental constraint the program has produced. Provenance
  (E2/exposure vs ladder tail vs missing relaxation channel) is
  unadjudicated. Sub-finding: v7.5's N = 50 Δn̄(MD − twin) = −0.68 was
  the outlier — N = 500 regresses to the −0.3..−0.6 interpolation band.
- **I88 (RP-D4 riders at the blessed point, §4aa).** **The 0.40 Landau
  spot is bit-flat** (RP-D7 mitigation complete; the v_L axis is quiet
  on the scored surface across its credible range at the winner cell
  family) and **every sensitivity-ring perturbation degrades W₁_solv**
  (+0.18..+0.33) — the pinned leg-D configuration sits at the basin
  optimum of the swept robustness levers. **The birth margin is the
  sensitive lever** (n₁_solv −0.094 / −0.171 at 4.67 / 6.0 Å,
  ≈ 1.9σ / 3.4σ; ratio collapses; suppression halves): the leg-D
  margin 3 Å is a pinned convention the histogram landing depends on.
  Droplet-prior axis asymmetric (δ 0.40 and pickup material ≈ 1.5σ;
  δ 0.80 noise-level). N = 50 χ² deltas direction-only (I85 caveat).
  Bands recorded, never re-fit.
- **I89 (RQ11 stage localization, §4bb).** **The deep-bin cold tail is
  set inside the 30 ps MD window.** KE is frozen at handover on every
  solvated bin (settle 30–70 ps into E2; E2 `E_dissip` gain 0.0000 —
  the Landau-gated E2 drag touches only the in-droplet retained class);
  detection-stage sheds are n ≤ 15 / KE-negligible; the Landau floor
  (≈ 0.003 eV) is an order of magnitude below the scored deep bins
  (≈ 0.055 eV) and hosts only the excluded retained class — explaining
  the §4y/§4aa v_L-quietness. RQ11 candidate (i) E2 exposure is
  **refuted**; the E2 cap is causally disconnected from the KE curve,
  so a larger-N re-run keeps the cfg diff at `{num_molecules, seed}`.
- **I90 (RQ11 in-window anatomy, §4bb).** **The deficit accrues in the
  in-band cubic phase, and fate is birth-dressing-ordered.** All fate
  groups launch identically (≈ 1.7–1.8 eV peak at 0.5 ps); the future
  n ≥ 12 survivors leave the above-cap regime at ≈ 2.1 ps (79 % of
  their dissipation) yet cross below the experimental deep-bin KE
  (0.069 eV) only at ≈ 6–30 ps (mean 13.5): the deficit-critical
  segment is in-band pure-cubic drag at v ≈ 2–7 Å/ps. Birth dressing
  (n_shell 15.5 → 18.4) and droplet radius (25.7 → 30.1 Å; retained
  34.0 Å) order the fate classes; the post-window n-mapping is
  near-diagonal (deep: ≈ 1 E2 shed), so the remaining RQ11 levers are
  the in-window (KE, n) exit correlation (ii), per-shed recoil ε
  (~4–8 meV/shed suffices, iii), or a beam-frame-dependent retained
  population channel (iv, RQ5-coupled).
- **I91 (battery seed-robustness + winner's curse, §4cc).** **The
  finc1v725 histogram landing is seed-robust at pooled N = 5000**:
  BN-P2 all in-band (n₁_solv 0.2433, midHot 1.0139, supp 0.187,
  trapped 0.067) and W₁ lands the converged branch (pooled 0.5713;
  true value ≈ 0.57 ± 0.04). First seed-scatter measurement:
  W₁ SD ± 0.095 — the §4z "N-resolution tail structure" (S6f-P2) and
  the s1 = 0.728 scare were seed scatter, not an N-effect. The
  blessed cell's draw was ≈ 2σ favorable on n₁/supp (winner's curse,
  quantified, mild); the pooled battery supersedes the single N = 500
  run as the statistical reference for the standing point.
- **I92 (RQ11 shape, §4cc).** **The deep-bin KE deficit is a monotone
  slope in n, not a uniform scale**: pooled sim/ref falls 0.81 → 0.38
  across n = 10–17 (301 → 17 ions/bin). χ²_med > 30 at every seed
  (114–155) and pooled (242.0) — no statistical rescue, per I89.
  Slope-compatible survivors: in-window exposure growing with
  dressing (ii) and per-shed recoil shrinking with terminal n (iii);
  uniform-scale explanations are excluded. Δn̄(MD − twin) regresses
  to −0.319 (in the registered band; 0.011 outside the tighter S6
  carry). E2-quietness holds to ≤ 0.1 meV at bin level (a ≲ 1 %
  deep-fragment residual of ≤ 5 meV late in-droplet drag,
  sign-irrelevant) — I89 stands.

---

## 6. Boundaries — what these results do NOT establish

- **N = 50, single seed:** means are stable (bridge-validated); distribution
  tails are not. Distribution-shape claims wait for the N=500 campaign.
- **9 Å / 0.80 eV / f_int = 0.5 / f_ret = 0.1 / κ = 1 / mixture pins**
  everywhere except where explicitly swept. In particular the Wave-4 cliff
  position (τ ≈ 6.55 strips, τ ≥ 16.5 freezes) is **f_int-conditional**
  (t× = τ·ln(f_int·E/Σ)); the floor pin would shift it.
- **TDDFT is not ground truth** — the staircase is a prior/anchor (9 Å
  non-radial flag stands); the arbitration observable is the experimental
  size distribution.
- **The staircase carries no meaning as a target (I14, §4c).** Reproducing
  21→19→14 in-window is a *capability showcase* only — the in-window trace is
  not the observable, and under the gate the cascade is unfinished at the
  window edge. Never treat a staircase "landing" as a calibration success;
  only the experimental terminal distribution at an explicit detection time
  does that.
- **No fidelity verdicts:** everything above is reported; adjudication is
  the user's (reporting-gate stance). The F5 gate is resolved only by the
  0.80 eV N=500 campaign.
- **The cooling gate is probe-scoped:** wired into the staircase-probe path
  (+ a latent `_cgds` campaign tag suffix from the review fix); the campaign
  generator and F3/F4 reports do not sweep it yet — promotion into the
  campaign is a separate decision.
- The relaxed reads at the 1000 ps cap with `frac_frozen < 1` (gated
  s_eff = 3, 5 at τ = 6.55) are cap-truncated, not converged; the gated
  τ ≥ 16.5 rows (frac_frozen = 0) are *by construction* never-converging.
- **The gated terminal read is flight-time dependent above the freeze-out
  s_eff (I11).** At the gated staircase-landing s_eff ≈ 30, `n_relaxed`
  (6.27 at 1000 ps) is a *snapshot of a live cascade*, not a terminal; the
  true endpoint (~n 2 over the ~8.5 µs flight) is not computed. Gated
  terminal-n numbers at s_eff ≳ 12 are cap-truncated and must never be read
  as converged. Only the *ungated* terminal reads, and the gated reads at
  low s_eff where frac_frozen ≈ 1, are converged.
- **NB (post-Wave-7, 2026-07-08) — the two cap-truncation items above are
  historical.** The detected read at the Sourced t_detect now exists for
  every probe dir (§4e); `n_relaxed` is a convergence diagnostic only
  (design §3.5), and the anticipated "~n 2 endpoint over the flight" was
  itself wrong at s_eff = 30 (the descent is logarithmic; the detector
  catches it at n̄ = 4.09, `time_exhausted`). Detected reads inherit the
  N = 50 / single-seed boundary plus a single detection-stage RNG
  realization per dir — distribution-shape claims still wait for the N=500
  campaign. Over the µs flight the detected read is RRK-evaporation-only by
  construction (OQ-E, design §4).

- **NB (post-Wave-8, 2026-07-09).** The single-seed boundary is now sharper
  than "tails are unstable": at a probe knob point the ensemble is
  kinematically **congruent** (I26), so *no* amount of N or seeds changes
  the delta cliff under this preset — the missing spread is structural
  (deterministic congruent onsets), not statistical. Wave-8 detected reads
  additionally verify that the whole pre-open analysis layer (exposure,
  opening times, |G| leaks) is deterministic and exactly reproducible from
  stored artifacts.

## 7. Open questions raised (documented, not built)

- **OQ-A (rho_min):** a residual out-of-bubble cooling floor (dragged He
  cloud) — would soften the Wave-4 cliff; deferred at the gate design.
- **OQ-B (bare-peak channel):** what expresses the 43 % bare I⁺ — the
  2.70 eV budget (moves only t× under the crossing construction, so *not*
  via more cascade budget), a lower f_int (earlier gate-open → longer
  post-ejection cascade under the gate), kinetic-energy assist, or a
  mechanism extension? Belongs to the F5/production discussion.
  *Post-Wave-7 sharpening (I22):* unreached at the detector too (min
  n_detect = 2); flight time cannot resolve it (the floor is kinetically
  asymptotic at s_eff ≳ 20), and the `suppressed` n = 21 class puts the
  fragmentation question directly into the observable.
  *Post-Wave-7 hypothesis (§4e interpretation, item 4):* the bare peak may
  *be* the suppressed class — a crossed-after-ejection complex carries
  E_int > Σ(n) and physically fragments rather than riding intact; the
  experimental bimodality then falls out of the Δ× cliff structure (one
  ensemble, two race outcomes). Document-level; the fragmentation
  channel's spec (statistical boil-off vs prompt, with/without OQ-F's ε)
  belongs to the OQ-B resolution, behind the trigger.
  *n = 1 corollary (§4e continued, I24):* the fragmentation channel covers
  the 17.5 % I⁺He bin **only with OQ-F's ε** (gateless boil-off without ε
  goes exactly to bare — G > 0 invariant, no barrier at the n = 1 direct
  channel); OQ-B and OQ-F are therefore one coupled mechanism discussion.
  *Post-Wave-8 resolution of the weight question (§4f, I26/I27):* "high
  enough f_int" does put the ensemble on the bare-candidate side — but
  **all-or-none** (100 % suppressed above E* = 0.4612 eV, 0 % below; the
  cliff is a delta at a knob point). The 43.5 % *fraction* is therefore
  not a knob outcome at all: it requires injected ensemble heterogeneity
  across E*_i (droplet-radius / onset / E_int(0) distributions). Two
  cascade-side corrections: the sub-rung sliver below E* expresses
  **n = 1** (I24's strangulation holds only at a full-rung leak), and the
  bare-vs-budget discriminator is void until OQ-H is resolved.
- **OQ-E (µs-flight channels, design §4):** over the 8.53 µs continuation
  the only active channel is RRK evaporation — no radiative cooling, no
  electronic relaxation, no residual-gas collisions. The detected read
  (Wave 7) leans on this by construction; a domain-expert question for the
  production discussion.
- **OQ-F (per-shed kinetic-energy release, post-Wave-7, 2026-07-08):** the
  delivered K1 drains only the binding energy per shed
  (`dE_int = −D₀(n)`); a statistical evaporation additionally releases
  translational energy ε ~ (E_int − D₀(n))/s ≈ 5–20 meV per shed —
  first-order over a ~17-shed cascade (~0.1–0.3 eV vs the entire
  Σ(21) ≈ 0.188 eV budget). Consequences if included: G strictly decreases
  → cascades self-terminate (no eternal `time_exhausted` class), terminals
  converge at higher n, deep-strip reach shrinks, OQ-B worsens. A
  mechanism-convention OQ of the RRK-dof class (§4e interpretation,
  item 3); document-only, no implementation. **Coupled to OQ-B:** the
  fragmentation channel's endpoint distribution (bare vs the small-n tail,
  incl. the n = 1 bin) is *set by* ε — resolve together (I24).
- **OQ-G (ladder-bottom depth, post-Wave-7, 2026-07-08):** the Form-U
  bottom is flat (D₀(1..5) = 9.22 meV, mixture κ = 1), but the real
  I⁺–He first rung is plausibly much deeper (ion-induced dipole, first
  shell) — which would make n = 1 a thermodynamic "last survivor" and
  explain its elevated 2.2× step over n = 2 without race statistics. The
  small-n abundance tail is the only identified observable that reads the
  ladder bottom; discriminates from the OQ-B/OQ-F picture via budget
  dependence (§4e continued). Literature/domain-expert question
  (I⁺–Heₙ binding energies); document-only.
- **OQ2 (S2-onset provenance — CALIBRATION_MAP register, FIRED
  2026-07-09):** the f_int parametrization `E_int(0) = f_int·E_avail` is
  physically indefensible in both directions (no actual partition
  implemented; literal ~1 % coupling leaves the mechanism inert — §4e
  "OQ2 fires"). Working hypothesis: absolute, budget-independent
  `E_int(0)` ≈ 0.2–0.5 eV from ionization reorganization +
  electronic/spin–orbit relaxation. Pending literature validation, then a
  reclassification decision (fraction → absolute eV; row-14 scenario
  keying dissolves).

- **OQ-H (production Coulomb kinematics, post-Wave-8, 2026-07-09):** the
  scenario budget `coulomb_available_eV` is bookkeeping-only (single
  physics reader: the S2 onset deposit) — the Coulomb-explosion mechanics,
  fragment speeds, ejection time, and cooling exposure are budget-blind
  (I28: K₂₇₀ ≡ K₀₈₀ exactly). Physically the 2.70 eV channel means faster
  fragments / earlier ejection / shorter exposure; none of that is
  modeled. Until it is (initial-separation or charge-state kinematics, an
  interchangeable arm behind the trigger), every budget-dependent
  prediction — B.1(2) in both versions, the RQ3-vs-RQ4 bare-vs-budget
  discriminator — is void in-model, while the absolute-E_int(0) map
  transfers across budgets verbatim. Couples RQ1 (the E_avail provenance)
  and the F5/production discussion.

**The consolidated register for all of these (OQ-B/E/F/G/H + OQ2 + the
resolved-but-uncrosschecked RRK-dof band) is `RESEARCH_QUESTIONS.md`
(created 2026-07-09) — the entry document of the literature-research /
domain-expert cross-validation phase.**
- **OQ-J (shed-frame / momentum convention of the evaporation channel,
  fired 2026-07-17 by the I.11.2 item-1 re-read, §4n):** the delivered
  Tier-2 generative evaporation composes the **cold-shed** operator
  (`mass_jump.cold_shed_velocity_components` in `physics/evaporation.py`):
  He left at lab rest, complex keeps its momentum, KE × m/m′ per shed.
  A8's grounding ([Nat23], Na⁺) is an *at-rest* result where "cold" and
  "co-moving" coincide; at production fragment speeds the conventions
  diverge first-order (+0.9 eV / ×1.61 on the n₁ bin), and Tier-1a had
  adjudicated **continuous-velocity (co-moving)** as the physical path
  for the anchored channel. Physical evaporation = co-moving + isotropic
  thermal recoil ε ⇒ OQ-J and RQ2/OQ-F are one coupled resolution.
  Consequences until adjudicated: every MD small-n KE read is
  convention-inflated relative to the twin/closure basis (I59); the §4k
  (v_c, τ) calibration does not transfer to a cold-shed MD at small n;
  leg-B KE pre-registrations must state their convention basis. Candidate
  resolution: an interchangeable shed-convention enum (cold /
  continuous_velocity, per the existing mass_jump operators) behind its
  own trigger, plus the RQ2 ε discussion. User adjudication required
  before the T9 endgame; → `RESEARCH_QUESTIONS.md` RQ10.
  **ADJUDICATED (user, 2026-07-17): co-moving is the working convention
  for the twin-parity legs (T5 onward).** Grounds: evaporation is
  thermal in the complex rest frame (u_thermal ≈ 0.5–1 Å/ps ≪ v_lab —
  co-moving is the correct zeroth order; cold shed needs a directed
  ~75 meV/He backward launch with no energy source, defensible only as
  a bound); [Nat23] is at-rest, where the conventions coincide — no
  contradiction with A8; Tier-1a already adjudicated
  continuous-velocity as the physical path (the Tier-2 channel's cold
  composition was a silent divergence, not a decision); the §4k
  closure/calibration is co-moving-based and transfers. Shape: an
  interchangeable `evaporation_shed_convention ∈ {cold, co_moving}`
  enum, **`cold` byte-inert default, legs stamp `co_moving`** (the T7
  birth-law precedent); build behind its own trigger — **DELIVERED
  2026-07-17** (log entry "Shed-convention enum DELIVERED";
  CALIBRATION_MAP row 26; full suite 2385 passed). **RQ10 stays
  open** as the physical-resolution question: the true convention is
  co-moving + isotropic thermal recoil — the ε → 0 limit vs the
  maximal-kick bound, one coupled RQ2+RQ10 discussion.
- **OQ-C (n-dependent s):** the scaled convention `s = α·(3n−3)` stays a
  documented alternative arm; not demanded by any current data.
- **OQ-D (conditional triggers, standing):** R6 9-Å re-extraction and the
  p↔κ occupancy-cap split remain document-only (firing criteria in the
  Phase-F plan §F5).
- **OQ-I (K-decoupling / joint (v_c, τ) calibration — Addendum I Step 1,
  2026-07-15):** under the frozen fate map both the detected KE and the
  evaporative descent ride the single exposure integral K, so no drag
  tail can lift the KE curve without collapsing the histogram (I47). The
  candidate decouplings, in increasing invasiveness: (a) **joint
  (v_c, τ) recalibration** — τ = 6.55 ps is a probe pin (GAH25,
  Bounded); a faster cooling clock restores descent at weak drag while
  v_c owns the KE scale; (b) fate-map-level change (per-shed ε > 0 /
  RQ2, or a non-exponential E_ej(K) form) — the F.5 escape clause; (c)
  the deep-bin slope specifically may implicate the *detected-KE side*
  of the deep bins instead (late in-band dynamics, Tier-0-locked — any
  revisit there needs new TDDFT input, not a knob). Documented, not
  built; the user adjudicates which arm (if any) opens.
  **[Arm (a) ANSWERED in-model same day (Step 1b, §4k): joint closure
  at τ = 4.1 ps — I49. Arms (b)/(c) do not need to open for this
  purpose. What remains of OQ-I is the τ re-classification question
  (GAH25 sourcing of the 6.55 ps pin vs the fitted 4.1 ps) and the
  RQ4-conditionality of the closure — both user adjudications recorded
  in the §4k Step-1b verdict.]**

## 8. Run inventory and reproduction

**116 probe dirs** under `data/runs/` (`*_tier2probe_*` namespace, disjoint
from the F3 campaign glob): 45 (Wave 1) + 10 (Wave 2; the two per-n controls
reuse Wave-1 dirs) + 10 (Wave 3, scratchpad-driven through the delivered
pipeline — no repo-code change) + 22 new in Wave 4 (24-point grid; the two
ungated s_eff=5 points at τ ∈ {6.55, 16.5} share dirs with Wave 2) + 5 new
in Wave 5 (gated s_eff ∈ {8,12,16,20,30} at τ=6.55, `density_scaled`;
scratchpad-driven, no repo-code change) + **10 new in Wave 6** (gated f_int
sweep: τ=6.55 s_eff ∈ {8,30} × f_int ∈ {0.24, 0.30, 0.35, 0.65} + τ=16.5
s_eff ∈ {8,30} × f_int = 0.25; scratchpad-driven, no repo-code change; the
f_int dimension rides the existing `_fiX.XX` tag) + **6 new in Wave 8**
(the cliff bracket, all gated τ=6.55: b080 f_int ∈ {0.57, 0.58} at s_eff=30
+ f_int=0.57 at s_eff=8, and the first-ever **b270** dirs — the 2.70 eV
anchor f_int=0.192593 (tag `fi0.19`) + the pair f_int ∈ {0.17, 0.18}, all
s_eff=30; scratchpad-driven, no repo-code change; `detection.npz` written
at generation time via the Wave-7 route) + **8 new in Wave 9** (the
E₀-mixture basis columns, all gated τ=6.55 b080: `s_eff=8` at
f_int ∈ {0.53, 0.48, 0.45, 0.42, 0.20, 0.15, 0.10} + one `s_eff=30` companion
at f_int=0.53; scratchpad-driven, no repo-code change; `detection.npz`
written at generation time). All gated dirs are stamped
`relaxation_time_ps = 1000.0` and carry the `cooling_spatial_gate` field in
`cfg.json`.

**Wave 7 (2026-07-08) added `detection.npz` to every one of the 102 dirs**
(zero new MD runs): a scratchpad driver through the delivered
`run_detection_stage`, seeding each dir from its `relaxation.npz` under a
transient detection-enabled cfg view (`detection_time_ps = 8.53e6`;
plus, for the 63 pre-Wave-4 dirs stamped with the legacy 8.53·10⁶ ps
relaxation cap, the *realized* relaxation duration so the config-load
nominal-window bound and the stage's realized-t_h check state the same fact
— §4e execution note). On-disk `cfg.json` files are unmodified; the detected
read is regenerable deterministically (stage-private RNG stream keyed on
`cfg.seed`).

- Generator: `scripts/gen_tier2_staircase_probe.py` (active USER SETTINGS =
  the Wave-4 A/B grid; the 12-run mini-probe and 45-run full grid are
  preserved as restore comment blocks).
- Scorer: `scripts/post_processing/tier2_staircase_probe_report.py` (pure
  scorer; table + optional CSV/figures; the numbers in §4.2 come from a
  full re-score of all 87 dirs on 2026-07-07; the §4e detected numbers from
  the 2026-07-08 re-score of all 102 dirs with `detection.npz` present —
  ion/relaxed columns unchanged, the wiring identity; the §4f numbers from
  the 2026-07-09 re-score of all 108 dirs — bridge + gated fi0.50/fi0.65
  oracles exact, ledger residual uniform. NB: the b270 rows' staircase
  columns are scored against the 9 Å / 0.80 eV anchor and are physically
  N/A — §4f run-count note).
- Stale-artifact policy (bridge findings §3) applies: any later
  config-surface change invalidates the dirs; regeneration is the recovery
  path.

**Cross-links:** `TIER2_STAIRCASE_PROBE_PLAN.md` (+ Addendum A);
`TIER2_PHASE_D_BRIDGE_FINDINGS.md` (the standing miss + levers, lever-3
resolution NB); `CALIBRATION_MAP.md` (s_eff row 10 Bounded, rows 19/20
annotations, row 5a cooling gate); `MASS_DYNAMICS_LOCKED_energy_gated_
evaporation.md` (§4 dof NB, A11 resolution, §6 K2 gate, §6.11 regime);
`TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (F2 Stage-1 re-scope NB, F5 gate);
`drag_migration_log_tier2.md` (delivery + decision history).
