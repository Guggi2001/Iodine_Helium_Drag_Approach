1C
# Tier 2 — Staircase / Capability Probe Findings

> **Status:** consolidated findings record, written 2026-07-07 after the fourth
> probe wave (cooling-spatial-gate total-strip A/B) executed; last updated
> 2026-07-11 after Wave 11 (the RQ7 production-kinematics probe — K₂.₇₀
> measured, §4i; Wave 10's width decomposition + K forward model is §4h;
> Wave 9's E₀-mixture inversion is §4g;
> Wave 8's cliff anatomy is §4f). This document
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

## 5. Consolidated insight register

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
- **OQ-C (n-dependent s):** the scaled convention `s = α·(3n−3)` stays a
  documented alternative arm; not demanded by any current data.
- **OQ-D (conditional triggers, standing):** R6 9-Å re-extraction and the
  p↔κ occupancy-cap split remain document-only (firing criteria in the
  Phase-F plan §F5).

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
