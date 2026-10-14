# Calibration Map — I⁺/He Drag-Port + Mass-Dynamics

> **Entry point:** `DRAG_PORT_DESIGN_DECISIONS.md` (§0 document map). This is the **cross-doc parameter index**; mass-model detail is in `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`.

One table mapping every model parameter to its **classification**, its
**tier anchor**, and its **cross-check**. Spans both
`DRAG_PORT_DESIGN_DECISIONS.md` (drag, gate, noise) and
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (pickup, evaporation,
cooling, ladder, early window). Physics-Definition only; no code.

> **Change history.** The reasoning behind each row's class lives in its
> source: the MASS doc revision notes, `DRAG_PORT_DESIGN_DECISIONS.md`,
> `docs/drag_port/Tier2/RESEARCH_QUESTIONS.md`,
> `docs/drag_port/Tier2/THESIS_FINAL_GOALS.md` and
> `docs/drag_port/Tier2/drag_migration_log_tier2.md`. This file's dated
> update log (2026-06-17 … 2026-07-10) was removed on 2026-10-01 and is
> kept in git history. **2026-10-14:** re-keyed from h405 to the thesis
> production law B2; h405 moved to its own contrast column; superseded
> classes kept as "History" inside each row.

## Classification scheme

- **Locked** — fitted to TDDFT (Tier-0 Method B), then fixed; not
  re-fit downstream. TDDFT authority covers only the fitted 18 Å trace
  band, 2.54–3.02 Å/ps; the 9 Å trace (2.83–4.95 Å/ps) enters only as a
  one-sided exit bound, not a fit (G4.9 / G4.12). The former
  2.54–4.95 Å/ps was the union of both traces (pre-G4.9). (Harmonized 2026-10-01 with D0, which used to say
  "Derived" for the same quantities.) Under B2 the Tier-0 output is a
  **one-parameter ridge**, not a point: the ridge is TDDFT-fixed, the
  position on it is Free (rows 1–3; THESIS_FINAL_GOALS G4.12 / G11.1).
- **Sourced** — value taken from external literature / ab initio, carried with a
  prior (and its transfer caveat).
- **Derived** — computed from other quantities; not independently fit.
- **Bounded** — physical bounds + sweep band; reported as a range/regime, not a
  point fit.
- **Free** — genuinely determined by fitting to data (the only true knobs).
  In production (B2) this means **selected by the experiment at G5**
  (Tier-2 detected size distribution + KE observables, `post_pickup`
  sampler): the ridge position $c$ and $E_\text{int}(0)$ only
  (THESIS_FINAL_GOALS G5 Outcome, G11.1).
- **Inactive** (added 2026-10-01) — declared and selectable in config, but
  not read by production: either fixed into a precomputed table or switched
  off by the selected arm. Its provenance is the value baked in, not a fit.

## Validation tiers (sequential, ordered by separability; DESIGN §6.4)

- **Tier 0 (COMPLETE)** — drag form + effective binding, deterministic, fixed
  mass, scored *in-window* against the TDDFT traces (Method B
  trajectory-matching). Thesis result: 18 Å-only extraction at
  m = I⁺ + 21 He and the 3.5 Å gate → ridge $F=a\,v+c\,v^2$; 9 Å is a
  one-sided exit bound (THESIS_FINAL_GOALS G4.9, G4.12, G11.2, G11.4).
- **Tier 1 (RETIRED 2026-10-13; historical)** — mass scenario,
  deterministic, full post-transient trajectory + the time-resolved shell
  trajectory (~21→19→14 He). With the 18 Å-only extraction at the true
  shell mass nothing is left to validate; the mass question is a
  discussion section (G11.3 / G11.5; `drag_migration_log_tier1a.md`).
- **Tier 2 (COMPLETE for the thesis)** — terminal I⁺Heₙ **size
  distribution** (+ KE observables) vs experiment. The arbiter for the
  Free set. *(NB 2026-07-01, Tier-2 Phase-E scope: the delivered
  observable is the **integer-n size distribution** — per-fragment
  velocity histograms were **cut** (the `vmi_summary/*.csv` references are
  aggregate I⁺He / I⁺-gas distributions, not per-n). See
  `docs/drag_port/Tier2/TIER2_PHASE_E_IMPLEMENTATION_PLAN.md`.)*
- **Tier 3 (open, not started)** — ensemble second moments: final-velocity
  histogram **width** + (if applicable) VMI angular spread. Isolates the
  noise model.

## Flag legend

`R#` risk · `A#` assumption · `OQ#` open question (author-contact/provenance).
All in the MASS doc unless marked **[D]** = DESIGN doc.

---

## The map

> **Production column (re-keyed 2026-10-14).** "Production (B2)" is the
> thesis production law (user, 2026-10-09): the value the generator
> `scripts/gen_tier2_detfix_postpickup.py --cell lq --b2` runs, read off
> the committed run
> `9A_drag_srho3p5A18A_linear_quadratic_N1000_detfix_conf270_g5lqppc160e3297s1_th500`
> (`cfg.json`). Reference = pooled N = 6000
> `…N6000_detfix_conf270_g5lqppc160e3297pooled_th500`; observables:
> THESIS_FINAL_GOALS G5 Outcome table. "Contrast (h405)" is the
> identifiability contrast law (standing point 2026-09-30 → 10-08;
> `gen_tier2_detfix_battery.py`, run
> `9A_drag_shared_pure_cubic_N1000_detfix_conf270_h405s1_th500`); "= B2"
> means the same value. The effective high-tail variant **B3** differs
> from B2 only in rows 2, 3 and 16 (a 19.8, E_bind 0.0482 eV, E₀ 0.3389
> eV; not Tier-0-consistent, G10.5). Per-knob influence:
> `TIER2_PARAMETER_INFLUENCE.md` (B2); "D0 §x" refers to the frozen
> `TIER2_PARAMETER_INFLUENCE_h405_finc1v725.md`.

| # | Domain | Parameter (units) | Production (B2) | Contrast (h405) | Role | Class | Primary anchor (Tier) | Cross-check | Flag |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Drag | drag law + ridge coefficient $c$ (amu/Å); $F=g\,(a\,v+c\,v^2)$ | `linear_quadratic`, **c 1.6**, uncapped (bundle `srho3p5A18A_linear_quadratic`, mass stamp 19 He) | `capped_cubic`, $b$ 2.5154 amu·ps/Å² ($F=g\,b\,v^3$ below $v_c$; bundle `shared_pure_cubic`) | drag form and its one free shape coefficient | **Free on the TDDFT ridge (B2, 2026-10-09)** — Tier 0 (18 Å only, m = I⁺ + 21 He, 3.5 Å gate) pins a one-parameter ridge: every c 0–7.2 passes 18 Å, c ≤ 5 passes the 9 Å exit bound (G4.12, G11.4). The experiment selects c 1.6 under `post_pickup` (G5; the raw sampler selects c 2.0 = A5). History: cubic $b$ 2.5154 **Locked** (Method B, jointly with $a$ and $E_\text{bind}$, 14.2 Å gate) — now only h405's base; the G4 refit at 3.5 Å gave $b$ 2.2300 (G4.8, superseded by G4.9) | Tier 0 ridge (18 Å) → Tier 2 (G5) | 18 Å RMSE 0.118 at m = 21 (band 0.19), 9 Å exit margin +0.82 Å/ps (G11.4); h405 fails the 9 Å bound (−1.21 Å/ps, over-dissipates); G4.11 cross-case cubic preference contested by the radial 9 Å runs (G11.1); history: `power_law` $n̂{=}2.93$ recovered cubic | 9 Å transverse flag **[D]**; identifiability (G5): no unique drag law |
| 2 | Drag | $a$ (amu/ps) | **18.007** | 0 (pure cubic below $v_c$) | linear drag coefficient | **Derived (ridge-tied)** — $a=a(c)$ from the 18 Å Tier-0 fit at fixed $c$ (`g412_ridge/c1p600.json`, m = 19 ridge; at m = 21 $a(1.6)=17.757$, not adopted — B2 passes Tier 0 unchanged at m = 21, G11.4). History: $a=0$ **Locked** (Method B drove $a\to0$, pure cubic) | Tier 0 ridge (18 Å) | same $F$ as A5 at 3 Å/ps, 7 % lower at 9.7 Å/ps (G5) | §1.2 **[D]** |
| 3 | Drag | $E_\text{bind}$ (eV) | **0.11727** | 0.1168 (bundle well `eb1168`) | effective ion binding, co-fit with drag | **Derived (ridge-tied)** — co-fitted with $a$ at fixed $c$ on the 18 Å trace (G4.12), then fixed; TDDFT pins the well / drag split below ≈ 0.08 eV (G10.5). History: **Locked** 0.1168 (cubic co-fit); G4 3.5 Å co-fit 0.1277 (G4.8, superseded); **identifiability flag** (h405-era, see "Flagged identifiability item") | Tier 0 ridge (18 Å) | +0.0004 eV at m = 21 (G11.4); B3's shallow well 0.0482 eV fails 18 Å (RMSE 0.580, G10.5); integrated ladder $\sum_i D_0$ (§6.5.1) | OQ1 (electronic provenance) |
| 4 | Drag | $v_c$ (Å/ps) | — (uncapped) | **5.5** | high-$v$ drag cap — `capped_cubic` tail onset (in-band $v\le v_c$ is the **byte-identical** pure cubic) | **Inactive in B2** (the lq ridge law has no cap; the high-speed shape is what TDDFT leaves open, G11.1). h405: **Free, fitted to the Tier-2 observables at G4, v_c = 5.5** (reclassified 2026-10-01; corrected-geometry basin 5.5–6.0, G3/G4; finc1v725's 7.25 breaks there). History: **Bounded (2026-07-16, §I.10 Slice T1; was contingent R10-(b))** — band $5.3\lesssim v_c\lesssim15$ (H.4; the ≥5.3 floor protects 0.80 eV in-window byte-identity, measured max speed 5.23); Step-1c refined targets $p{=}{-}1$: 7.5 / $p{=}0$: 6.0–6.5 → **Derived** at the T4 winner | Tier 2 $(n,\langle E\rangle)$ mean-KE curve (I-D4 full-curve fit) | Slice-T1 oracles: $v_c{=}\infty$ / $v_c\ge v_\text{max}$ byte-identity vs `linear_cubic(a{=}0)`; delivered 0.80 eV bridge probe dir reproduced **byte-identically** at $v_c{=}5.3$ (both tails) | R10, §I.10 |
| 4b | Drag | $p_\text{tail}\in\{0,-1\}$ (dimensionless) | — (uncapped) | −1 (saturated constant-force tail) | `capped_cubic` tail exponent, $\gamma_\text{tail}=b\,v_c^2(v/v_c)^{p}$ ($p{=}0$ Stokes-like linear force, $p{=}{-}1$ saturated constant force) | **Inactive in B2.** h405: **Free choice (2 arms, 1 selection; config-load-restricted 2026-07-16, §I.10 Slice T1)** — $p{=}1$ (Newton tail) **excluded** (Step-1c K-P2: zero joint closures); hard cutoff **excluded as physics** (I48) | Tier 2 $(n,\langle E\rangle)$ + solvated histogram (C1/C2 form discrimination) | discrimination read (I50): $p{=}{-}1$ wins $n_1/n_2$ ratio, sags $n\ge13$ KE tail; $p{=}0$ flattest KE curve of the program | §I.10, I50 |
| 5 | Gate | $g(\text{depth})$ (dimensionless) | **`erf_independent`, $s_\rho$ = 3.5 Å** (`drag_gate_steepness`; shared surface with rows 5a, 8, 27; = the bundle's extraction gate) | `density_proportional`, erf width 14.2 Å | drag spatial gate, G4→G2 | **Sourced (borrowed), adopted with B2** — 10–90 % width 1.81·$s_\rho$ = 6.3 Å vs 5.7 Å DFT / 6–8 Å experiment for ⁴He droplets (Harms et al., PRB 58, 3341 (1998); THESIS_FINAL_GOALS G2); the solvation potential keeps 14.2 Å. History: **Derived** (erf-tied reuse of the 14.2 Å potential width; density-width scan NULL at h405, D0 §9.6.5) | external (droplet density profile) | Tier-0 ridge extracted at the same 3.5 Å gate (G4.9/G4.12); `h405g35` (h405 at 3.5 Å: W₁ 0.818, KE₁ 0.643, G5) | §5 **[D]** |
| 5a | Gate | `cooling_spatial_gate` (arm) | **`density_scaled`** (not the default `none`) | = B2 (on the 14.2 Å surface) | K2 cooling spatial gate: `none` (locked, ungated $-E_\text{int}/\tau$) / `density_scaled` ($\tau_\text{eff}{=}\tau/\rho_\text{ratio}$, off outside bubble) | **Free choice** (2 arms, 1 selection; default `none` **byte-identical**; added 2026-07-06) | Tier 2 size dist (total-strip A/B capability probe) | reuses **row-5** drag gate surface (`drag_gate_steepness`, one bubble boundary); ejection-regime total-vaporization lever (shell-retaining ↔ near-bare); $\rho_\text{min}$ floor **deferred OQ** | MASS §6 K2/§6.11 |
| 6 | Noise | $T_\text{eff}$, FDT amplitude | **off** (`noise_form="none"`; Tier 3 not started) | = B2 | multiplicative local-FDT bath kick (N2) | **Bounded/Derived** (tied to $\gamma(v)$; $\propto\sqrt{\gamma g k_BT_\text{eff}}$) | Tier 3 (ensemble width / VMI) | strict-FDT shown dynamically null (§1.3a) | §1 **[D]** |
| 7 | Pickup | $\lambda_\text{attach}(\rho_\text{He},n)$ (ps⁻¹) | $\lambda_0=0.9$ ps⁻¹ (live through all of the 500 ps Stage I) | = B2 | Poisson He capture rate (density-gated, occupancy-capped) | **Sourced + Bounded** — central ~0.7–1.1/ps from GAH25 (García-Alfonso, Barranco, Pi, Halberstadt, J. Chem. Phys. 163, 144309 (2025)) **Table II**, slope A of $n_1(t)=A\,t$ ($n_1\le5$; ⁴He-TDDFT, Ak⁺@He₂₀₀₀, sudden surface ionization): Na⁺ 1.3, K⁺ 1.0, Rb⁺ [1.1, 0.5] (two values from n₁ oscillations), Cs⁺ 0.7. Transfer to I⁺ via **Table I** ion–He pair parameters: Rb⁺ $D_e$ 204.1 K / $R_e$ 3.10 Å, Cs⁺ 168.9 K / 3.37 Å; I⁺ $X_2$ (IHe05) 143.9 cm⁻¹ ≈ 207 K / 3.25 Å → between Rb⁺ and Cs⁺, nearest Rb⁺. Band = Cs⁺ 0.7 … Rb⁺ upper 1.1; the Rb⁺ lower value 0.5 is **not** included (no recorded reason). **Not** 2.0: Na⁺ experiment (Albrechtsen et al., J. Chem. Phys. 162, 174309 (2025), quoted GAH25 p. 6: 2.04 ± 0.13 /ps at N = 9000 → 1.65 at N = 3500). Transfer caveat: surface-born ion solvating in a 2000-atom droplet vs a Coulomb-explosion I⁺ fragment; ±factor-2; ×Langmuir cap $(1-n/n^*)_+^{\,p}$ (A12, 2026-06-21) | Tier 2 size dist | 9/18 Å density contrast; cap inert for ejection; untested pickup-side n̄ lever under `post_pickup` (D0 §15.9) | R1 (Na⁺→I⁺), R7 ($v$-dep), A12 |
| 7p | Pickup | $p\ge0$ | 1 (`pickup_occupancy_cap="langmuir"`) | = B2 | occupancy-cap sharpness (Langmuir exponent, A12) | **Free (conditional)** — Phase-B default **$p=1$ fixed** (NB 2026-07-01: the $p\!\leftrightarrow\!\kappa$ tie is **inverse** — rigid shell = large $\kappa$ = smaller $p$; **NOT** $p=\kappa$; match cutoff slopes if tied); split only if size dist. demands | Tier 2 size dist | [Nat23] resting-ion saturation bounds $n^*,p$; $p\!\leftrightarrow\!\kappa$ coupling (inverse); Phase-B §3.2 | A12, OQ4/OQ8 |
| 8 | Pickup | $\rho_\text{He}(\text{depth})$ (Å⁻³) | analytic erf-complement $\hat\rho$ (`rho_he_ratio`), **width 3.5 Å** (row-5 surface) | width 14.2 Å (the solvation-*potential* width reused for the density; D0 §17, RQ12) | density profile gating capture | **Sourced** — profile shape baseline, width per row 5 (G2). History: the 14.2 Å reuse (convention, D0 §17 **C**) | baseline + G2 | 9/18 Å density contrast | — |
| 9 | Evap | $\nu$ (ps⁻¹) | 2.42 (`evap_rate_prefactor_per_ps`) | = B2 | RRK prefactor (I⁺–He stretch freq) | **Sourced (pinned 2.42, 2026-06-17)** — $\omega_e{=}80.6$ cm⁻¹ from IHe05 EPAPS $V''(R_e){=}748.1$ cm⁻¹/Å² | Tier 2 size dist | I2-notes cascade timing (OQ5); near-threshold spacing → effective $s$ | OQ5, A11 |
| 10 | Evap | $s_\text{eff}$ (dimensionless) | 8 (constant) | = B2 | RRK effective vibrational DOF | **Bounded (promoted 2026-07-06; was Derived)** — constant $s_\text{eff}$, band **≈[5, 20]** (staircase landing [8, 12] at mid-band τ; picture-robust ≤4 %); classical mode-count $s{=}3n{-}3$ (corr. 2026-06-21 from $3n-6$) demoted to the classical-limit arm of the dof-convention selection; $n{=}1$ direct dissociation $k=\nu$ and the $s\ge1$ guard unchanged | 9 Å anchored staircase (magnitude; s↔τ separated by first-shed timing; TDDFT prior) + Tier 2 size dist | classical band $[3n{-}3,3n]$ superseded (landing far below: quantum mode-freezing / blurred shell — the A11 down-drift); κ-joint resolved **weak** | A11 (resolution NB 2026-07-06) |
| 11 | Cooling | $\tau_\text{dissip}$ (ps) | **6.4** | 4.4 (fitted at G4; ×0.60 of GAH25's first-shell 7.3) | Newton-cooling time of $E_\text{solv.struct}$ | **Bounded; value selected by preliminary experimental calibration, frozen at G5** — external sources, **Na⁺ only**: GAH25 (J. Chem. Phys. 163, 144309 (2025)) **Table III**, Newton fit (Eq. 12) of $E_\text{solv.struct}(t)$, ⁴He-TDDFT Na⁺@He₂₀₀₀: first shell τ = 7.3 ± 0.2 ps (t₀ 6.53, E∞ −3424 K), second shell τ = 16.5 ± 0.6 ps (t₀ 5.0, E∞ −4144 K); Na⁺ experiment τ = 2.6 ± 0.4 ps at N = 3600 (Albrechtsen et al., J. Chem. Phys. 162, 174309 (2025), quoted GAH25 p. 10, which calls the gap "puzzling"). Band $[2.6,16.5]$ = experiment … GAH25 second shell. **6.55 ps is not a GAH25 value**: it is the band's geometric mean $\sqrt{2.6\cdot16.5}$, a config default (`drag_migration_log_tier2.md`, "τ default = geometric mid 6.55 ps"), not to be confused with GAH25's onset t₀ = 6.53. No Rb⁺/Cs⁺-like cooling time exists in GAH25; it cools a growing shell (ours shrinks). The value 6.4 is **not** externally prescribed: it was selected against the same Tier-2 observables at a different law (G4.13 τ arms 4.8–6.4 → G4.14 Step-1 τ lines, 74 MD runs, fitted KE₁ vs W₁ trade-off → amendment A5 cell c 2.0 / τ 6.4, KE-led end of the curve; c 2.0, `raw` sampler). 6.4 sits at ×0.88 of GAH25's first-shell 7.3 and ×0.98 of the 6.55 geometric mean — noted after A5, not the selection basis. B2 (c 1.6, `post_pickup`) inherits τ 6.4 unchanged; **never re-scanned under B2** (G5 scanned c and E₀ only). Not a G5 knob (G11.1: τ trades the histogram tail against KE₁ at every c, G4.14). History: h405 **Free**, fitted at G4 τ = 4.4 ps (reclassified 2026-10-01); before that **Bounded** with **RQ9 reclassification pending (2026-07-16)** — the §I.8/§I.9 joint $(v_c,\tau)$ closure at $\tau\approx3.0$–$4.8$ ps | external (GAH25 Table III / exp) | dn̄/dτ ≈ −1.1 /ps under `post_pickup`, restoring n̄ costs midHot (D0 §15.9); $t_\times$ vs GAH25 $t_0$ | R8 (Na⁺-only), RQ9 |
| 12 | Cooling | $E_\infty(N)=-\|S(N)\|$ (eV) | $\|S\|=0.308$ (`solv_struct_asymptote_eV`) | = B2 | K2 asymptote, **occupancy-resolved** | **Sourced (split, 2026-06-17)** — full-shell $\|S_{\mathrm{I^+}}\|{=}0.308$ eV [I2-notes]; GAH25 −3424/−4144 K Na⁺ fixes form+$\tau$ | equilibrium-shell binding | OQ6 resolved (fixed $E_\infty$ caps strip); shape via $\kappa$ | OQ6, OQ7 |
| 12b | Cooling | $E_\text{elec}(N)$ (eV) | derived from row 12 and the production ladder | = B2 | electrostriction binding (collective excess) | **Sourced/Derived** ($-(\|S(N)\|-\sum_i D_0)$; **dominant** per GAH25 geometry) | — | marginal release → bath on shed (A8); $\partial\|S\|/\partial n{\approx}118$ cm⁻¹ ($n^*{=}21$, corr. 2026-06-21) ≈ $D_0(1)$ | K2, A8, OQ7, OQ8 |
| 13 | Budget | $f_\text{ret}\in[0,1]$ (dimensionless) | 0.1 (never swept; D0 GAP) | = B2 | S1 pickup binding-release retained fraction | **Bounded** (prior small) | Tier 2 size dist | 9/18 Å density contrast (feedback gain → density-dependence of terminal $n$) | — |
| 14 | Early | $f_\text{int}\in[0,1]$ (dimensionless) | **0.1221** (= 0.3297 / 2.70) | 0.15 (= 0.405 / 2.70) | S2 onset **sweep coordinate**, $f_\text{int}=E_\text{int}(0)/E_\text{avail}^\text{ion}$ | **Derived coordinate (2026-07-09, RQ1 adjudication; was Bounded)** — the physical variable is the absolute row-16 $E_\text{int}(0)$; scenario keying dissolved (floor identity: $0.235{\cdot}0.80\approx0.065{\cdot}2.70\approx\Sigma(21)$); budget transfer holds $E_\text{int}(0)$ fixed ⇒ $f$ rescales ×0.80/2.70 | Tier 2 size dist (via row 16) | historical scenario-keyed floors (2026-06-21, superseded as primary read): 0.065 (mix)/0.09–0.10 ($X_2$) @ 2.70 eV; 0.21–0.24/0.31–0.35 @ 0.80 eV | R1 regime axis, OQ2→RQ1 (adjudicated) |
| 15 | Early | $E_\text{avail}^\text{ion}$ (eV) | 2.70 (`coulomb_available_eV`; $R_0=2.666$ Å, `E_coulomb_scale` 1.0), consistent with the TDDFT reference | = B2 | per-ion Coulomb onset budget | **Sourced/fixed, scenario-keyed (2026-06-21)** — **0.80 eV** validation (½·14.40/9, $d{=}9$ Å) / **2.70 eV** production (½·14.40/2.666, $R_e$). The Hatherly 1994 0.8·E_C calibration (2.16 / 4.32 eV per channel) is discussed, not adopted: `THESIS_DISCUSSION.md` §9 | fixed reference ($t^*$-window) | stamped to scenario guard (§6.5); reversible to pair (1.60/5.40); D0 §19 sensitivity $S_k$ 0.389 eV/eV (h405-era) | OQ-none |
| 16 | Early | $E_\text{int}(0)$ (eV) | **0.3297**; every ion starts at exactly this value (row 28 inert) | 0.405 | initial internal energy (S2 onset, **absolute** — the physical variable) | **Free (B2): selected by the experiment at G5 jointly with $c$** — re-placed per $c$ to n̄ ≈ 3.9 under `post_pickup` (B1c → B2; D0 §15.9), inside the Bounded RQ1 band below. History: h405 **Free**, fitted at G4 0.405 eV (reclassified 2026-10-01); before that **Bounded (2026-07-09, RQ1 adjudication; was Derived)** — sourced band **0.2–0.5 eV** (floor ≈0.22 eV solvation, Na⁺ benchmark; admissible tail ~0.8–1.0 eV full-³P deposit, discounted); **budget-independent**; deposited at $t_0$, nothing subtracted from fragment mechanics | Tier 2 size dist + midHot (G5) | dn̄/dE₀ ≈ −31 /eV (D0 §15.9); VMI discriminator: fine-structure satellites (~0.8 eV down-shift), **not** a smooth deficit; provenance in the RQ1 NB register | RQ1 (closed at convention level) |
| 17 | Early | $t_\times$ (ps) | diagnostic (not set) | = B2 | self-bound crossing (gate-open time) | **Derived diagnostic** | — | GAH25 $t_0$ = 5.0/6.53 ps (±factor-2, Na⁺) | §6.11 |
| 17b | Early | $\Pi(t)=\lambda(n)f_\text{ret}\tau$ (dimensionless) | diagnostic (not set) | = B2 | pickup↔gate order parameter (NEW 2026-06-21) | **Derived diagnostic** (§6.11) | — | $\Pi{>}1$ shed / $\Pi{<}1$ freeze; regime-axis spine; $\Pi\to0$ at exit ⇒ termination guaranteed; reconstructable post-hoc | §6.11, R1 |
| 18 | Ladder | $D_0^{\mathrm{I^+}}(1)$: $X_2$ 106.9 / mix 74.4 cm⁻¹ | **163.7 cm⁻¹ (20.29 meV)**: the mixture 74.4 × 2.2 (`rq4graded` taper, rungs 1–3 × 2.2 / 1.5 / 1.3). The production rung is **not** the bare sourced value | = B2 | first dissociation rung (picture-dependent) | **Sourced base × experiment-derived taper** — three separate provenances: (i) **base Sourced**: mixture 74.4 cm⁻¹ (±3) from IHe05 (Buchachenko/Viehland JCP 122 (2005), EPAPS fit + exact J = 0 ZPE; equal-weight (106.9 + 62.4 + 54.0)/3). (ii) **Taper experiment-derived, not literature**: 2.2 / 1.5 / 1.3 are the abundance step ratios of the experimental detected I⁺(He)ₙ histogram (n = 1–4 weights 0.175 / 0.080 / 0.052 / 0.040), mapped to D₀ ratios by the model assumption abundance(n) ∝ D₀(n) under ε ≈ 0 and a smooth reservoir (RESEARCH_QUESTIONS RQ4 "sharpened target"). No publication or author communication supplies them; the RQ4 literature search found no analog with a decreasing-from-n = 1 first-shell ladder (every characterised cation: flat plateau → RQ4 outcome (ii)), so the external evidence is **against** the taper, not for it. (iii) **Family selected on the same data**: rq4graded chosen over floor1 on the experimental size-distribution reads (W₁_solv, χ²_med, n₁_solv) at I86 (S6, N = 500; finc1v725-era capped-cubic law, not B2; D0 §8). The same experimental distribution therefore sets the taper, selects the family, and later arbitrates B2 (G5): D₀(1) = 163.7 cm⁻¹ is **not externally constrained** beyond its base. Ladder fixes → outlook (THESIS_FINAL_GOALS status 11). History: **Sourced (pinned 2026-06-17)** — IHe05 EPAPS fit, exact $J{=}0$ ZPE $G(0){=}37$ cm⁻¹ ($D_e{=}143.9$ is **not** $D_0$); $\pm3$ cm⁻¹ | external ab initio | mobility/ZEKE-validated curve; mixture/$X_2$=0.70 | OQ1 (which curve), M1 |
| 19 | Ladder | $\kappa$ (Form U) / `ladder_steepness` | **not a runtime knob**: `dissociation_ladder="tabulated"` overrides it. κ = 1 enters only when the generator builds the rq4graded table (`LADDER_KAPPA_PIN`) | = B2 | sigmoid steepness, gradual↔cliff (1 knob) | **Inactive (reclassified 2026-10-01)** — baked into the rq4graded table at κ = 1. History: **Free** (prior large: 7.5× radial cliff) | Tier 2 size dist (broad vs magic-peak) | $n^*{=}21$ [I2-notes]; joint w/ picture + $s$; **staircase probe 2026-07-06: inverted + normalisation-capped for 21→14** (cliff centre $n^*{+}\tfrac12$ above the stripping range; Form U floors $D_0(21)$ ≈ 0.53·$D_0(1)$) | R3, A5, OQ4 |
| 20 | Ladder | `ladder_electronic_picture` | **not a runtime knob**: baked into the rq4graded table at `statistical_mixture` (`LADDER_PICTURE_PIN`) | = B2 | statistical-mixture (default) / $X_2$-only / cooling_relaxed | **Inactive (reclassified 2026-10-01)** — baked into the rq4graded table at `statistical_mixture`. History: **Free choice** (3 options, 1 selection; `cooling_relaxed` added 2026-06-21, between mix & $X_2$) | Tier 2 size dist | rung scale: mix 74.4 / $X_2$ 106.9 cm⁻¹ (=0.70); **cross-check 2026-07-06: near-degenerate in shed magnitude** ($x=D_0/\Sigma$ picture-invariant, ≤4 % under reduced $s_\text{eff}$); timing-only sensitivity via $\Sigma(21)\to t_\times$, degenerate there with $f_\text{int}$ | A10, OQ1 |
| 21 | Ladder | $\sum_i D_0^{\mathrm{I^+}}(i)$ (eV) | Σ(21) = **0.206** (the rq4graded table; above the Form-U mix band 0.17–0.19 because of the rung 1–3 taper) | = B2 | integrated ladder = self-bound gate threshold | **Derived** — $X_2$ **0.25–0.28** / mix 0.17–0.19 eV (Form U pure-$\kappa$ range, corr. 2026-06-21; ~11%, **nearly $\kappa$-indep**, cliff at $n^*{+}\tfrac12$) | — | drag binding (row 3) + crowding-reduced value are **separate cross-checks, not band ends**; $\|S\|{=}0.308$ collective UB | A10, OQ7 |
| 22 | Ladder | $n^*=21$ | 21 (every ion starts at $n_0=21$; row 27) | = B2 | first-shell equilibrium occupancy | **Sourced** ([I2-notes] $X_2$-only; **adopt 21 for all per-atom energetics**, 2026-06-21); = the Tier-0 extraction shell (G11.2) | external | GAH25 $R_e$-scaling → ~20 (corroborates to ±1–2); n* beyond 21 ± 1–2 not defensible (D0 §1.G9b) | OQ4, OQ8 |
| 23 | Ladder | $D_\text{floor}$ (cm⁻¹) | 4.97 (the table tail, rung 32) | = B2 | outer-shell rung floor (Form U) | **Sourced** ($\|\mu_\text{He}^\text{bulk}\|\approx4.97$ cm⁻¹, 7.15 K) | external (bulk superfluid) | picture-independent; sets sigmoid lower anchor | R3 |
| 24 | Detection | $t_\text{detect}=8.53$ µs ($8.53{\cdot}10^6$ ps) | 8.53 µs (`detection_time_ps` 8.53·10⁶) | = B2 | terminal-read detection time (TOF flight to detector) — where the Tier-2 size distribution is read | **Sourced (2026-07-07)** — experimental-setup publication, user-confirmed | Tier 2 size dist (the *detected* read is the arbitration observable) | log-sensitivity band = free report-side re-read of stored shed events (detection-stage design §3.3) | — |
| 25 | Ensemble | `birth_position_law` + `initial_position_margin_angstrom` (Å) | **`boltzmann`** at 0.4 K with a well of **313.2 K** (the DFT fit 26.99 meV, per-run override; config default 573.3 K). Margin **0, inactive** (guard-refused under `boltzmann`). Mean birth depth ≈ 36 Å, min ≈ 20 Å (B2 members s1–s5 + r6, `post_pickup` droplets) | = B2 law; mean birth depth ≈ 34.9 Å (s1, raw droplets) | molecule-centre birth law: `boltzmann` (delivered thermal sampler, default) vs `uniform_volume` ($p(r)\propto r^2$ on $[0, R-m]$, hard margin — the 1D twin's L1 law; Slice T7, 2026-07-16) | law = **arm** (production `boltzmann`); margin **Inactive (reclassified 2026-10-01)**. History: law = **arm, not knob** (twin-parity/capability lever — §I.11.0 NB); **audit 2026-07-26: the Boltzmann center-pin is an $R$-dependent artifact (erf width 14.33 Å absolute ⇒ allowed shell $r\lesssim R-35$ Å; D0 §15.4)**; **parent-document anchor (2026-07-26): the same law is the parent model (its Fig. 6.9), whose mean solvation depths 29 Å @ $R{=}34$ / 40 Å @ $R{=}68.3$ our sampler reproduces to 0.34 / 0.43 Å at the DFT-fit well. Port-fidelity reproduction, NOT independent validation; depths stay Derived.** The finc1v725-era `uniform_volume` + 3 Å pin is retired with G2; margin = **Bounded** — firm band {3, 4.67, 6} Å (H.2b D5), read only under `uniform_volume` | Tier 2 (T9 oracle-chain legs A′–D) | margin read only under `uniform_volume` (guard-refused otherwise); default byte-identical to the pre-T7 sampler (exact-draw regression) | — |
| 26 | Channel | `evaporation_shed_convention` | **`co_moving`** (not the default `cold`) | = B2 | shed velocity/momentum convention: `cold` (He at lab rest, $v\to\frac{m}{m-m_\text{He}}v$, negative reduced-mass ledger; byte-inert default) vs `co_moving` (He leaves co-moving, $v$ unchanged, positive $+\tfrac12 m_\text{He}\|v\|^2$ ledger; OQ-J build 2026-07-17) | **arm, not knob** — convention pair; **working convention adjudicated `co_moving`** (user 2026-07-17, findings §4n/I59: cold injects +0.9 eV / ×1.61 on $n_1$ KE); `cold` retained as the diagnostic bound arm | Tier 2 $(n,\text{KE})$ curve (histograms are convention-blind) | RQ10 brackets: cold = maximal backward kick, co-moving = $\varepsilon\to0$; every KE claim states its convention basis | OQ-J→RQ10 |
| 27 | Ensemble | `initial_shell_model` | `density_tied`, **structurally inert**: at the 3.5 Å surface $\hat\rho$ saturates beyond ≈ 7 Å depth and the shallowest B2 birth is ≈ 20 Å, so $n_0=21$ for every ion | = B2 arm, inert on the 14.2 Å surface (s1: 2000/2000) | t = 0 shell-dressing law: `full` (21-for-all; byte-inert default) vs `density_tied` ($n_{0,i}=\operatorname{round}(n^*\cdot\hat\rho(d_{\text{birth},i}))$ through the shared erf-complement surface; Slice T5, 2026-07-18) | **arm, not knob** — parameter-free; a first-order occupancy statement, no shell-restructuring dynamics. Bare / partly-solvated starts = outlook (G8, G11.5) | Tier 2 (T9 leg B) | biphasic-only by config-load guard; dressing↔pickup re-filling is a reported twin-divergence candidate; E_int(0) onset coupling is row 28 | — |
| 28 | Ensemble | `internal_energy_partition_law` | `sigma_proportional`, **structurally inert**: Σ(n₀)/Σ(n*) = 1 because $n_0=21$, so $E_\text{int}(0)=0.3297$ eV for every ion | = B2 arm, inert; $E_\text{int}(0)=0.405$ eV | how the S2 onset $E_\text{int}(0)$ couples to the T5 dressing $n_0$: `constant` ($f_\text{int}\,E_\text{avail}$; byte-inert default, $p=0$) vs `sigma_proportional` ($\times(\Sigma(n_{0})/\Sigma(n^*))^{1}$; Slice T6, 2026-07-18) | **arm, not knob** — the p-law is parameter-free ($p\in\{0,1\}$); the $p=1$ prior (§4j finding 1) is not hard-wired | Tier 2 (T9 leg C) | biphasic-only by config-load guard; structurally inert at $n_0=n^*$; touches the $E_\text{int}(0)$ onset **only** | — |
| 29 | Channel | `relaxation_dissipation` | `landau_gated_drag` is set but **unread**: E2 is skipped since 2026-09-30 (`relaxation_stage_enabled=False`). $v_L$ = 58 m/s = 0.58 Å/ps (sourced bulk value), read only by the hard-sphere collision cutoff `E_min_eV`. Stage-I drag has no Landau gate | = B2 | E2 post-ejection dissipation arm: `zero_gamma` (conservative closure; byte-inert default) vs `landau_gated_drag` (§I.11.2 item 2; build 2026-07-20): $\gamma=0$ for $v\le v_L$, the production drag above | **arm, not knob** — parameter-free given $v_L$ + the drag bundle; gate on **speed** (preserves the mass-agnostic drag contract) | Tier 2 — the marginal droplet-retained class only | coulomb-only by config-load guard; 5-term invariant closes | §I.11.2 item 2 |
| 30 | Ensemble | `droplet_size_prior` + `droplet_size_sampler_mode` | **`legacy`** + **`post_pickup`**: nozzle correlation at 40 bar / 14 K, pickup-selected (∝ N^(2/3), one pickup, pickup-induced evaporation) → ⟨N⟩ ≈ 16.7k, R̄ ≈ 54 Å (B2 members) | `legacy` + `raw`: ⟨N⟩ ≈ 12.8k, δ 0.625, R̄ ≈ 49.9 Å (s1) | droplet-size prior of the neutral-stage seed: `legacy` (fixed $N{=}2000$ or the ported pickup-cell MC) vs the analytic D4 family `kornilov_lognormal` / `pickup_weighted_lognormal` (Slice T8, 2026-07-20); sampler `raw` vs `post_pickup` | **arm (selection)** — `post_pickup` adopted for the thesis (THESIS_FINAL_GOALS G1: the experiment selects doped droplets; the legacy MATLAB production sampler), a stated modelling choice with measured sensitivity, not a validated property. History: `raw` (h405; basis a circumstantial quantile match, margin ≈ the N→R convention spread); the analytic ⟨N⟩ = 2000 pin is finc1v725-era, retired with G2 (**audit 2026-07-26:** inherited from the 9 Å TDDFT droplet, never re-derived; D0 §15.5) | Tier 2 (G5 A/B) | raw → `post_pickup`: n̄ +0.70 (law-independent, h405 +0.72), KE / shape flat on lq; chord map sampler-invariant (D0 §15.8–15.9); D4 marginal-partner safeguard tracked (no firing in the B2 pool) | — |

### Cross-cutting (not single parameters)

| Item | Nature | Status | Tier / check | Flag |
|---|---|---|---|---|
| Mass↔drag coefficient pairing | consistency, not a parameter | structural: production runs `allow_inconsistent_mass_pairing`; bundle stamp stays 19 He (202.954 amu). Defence: B2 passes Tier 0 unchanged at the 18 Å shell mass m = 21 He = production's $n_0$ at birth (G11.4); afterwards rests on the mass-agnostic force law — a **stated assumption** (G11.5) | G11.4 step 2; history: §6.6 mid-window defence | R6 |
| Total-stripping (Calvo24) limit | reachable far end of regime axis | secondary/sensitivity-run evaluation, not default | Tier 2 size dist (terminal $n\to$ small) | OQ6 ($E_\infty$ reach) |
| Terminal regime ($f_\text{int},\tau$, ladder depth) | reported output, not a knob | regime determination (shell-retaining ↔ total strip) | Tier 2 size dist | R1, §6.11 |
| Integrator↔mass-jump operator split | accuracy + conservation, not a parameter | **SQ1–SQ3** built (Tier 1a) and in production (variable-mass O-step, post-jump $m^+$) | implementation precondition; $dt$-convergence | A13 |
| Detector-stage pipeline (handover $t_h$ + Coulomb closure) | staging convention + exact closure maths, not a fitted parameter | **adopted 2026-09-30, unchanged for B2**: Stage I full physics to $t_h=500$ ps, E2 skipped (`relaxation_stage_enabled=False`), `detection_coulomb_closure="partner_aware"`, exact jump chain to $t_\text{detect}$. $t_h$ has no measurable effect (500 vs 502.84 ps: ≤ 4·10⁻⁶) | `TIER2_DetectorStageFix.md` §5.6–5.7 | limitations: point charges, fixed droplet, no mass change after $t_h$ |
| Retained-ion policy | scoring convention | production `exclude_all_coupled`; the marginal class is **0** at B2 (trap 0.049, all bound; h405: 0.0014), so the policy question is empirically empty | G5 winners table (`g5_winners_table.csv`) | D0 §14.2 |
| Suppressed ions | scoring convention | the n = 0 class, excluded from the solvated read (B2 0.111) | G5 Outcome; G6 figure | RQ8 |
| CE channel mixture + exit strip | candidate replacement for row 15's scalar budget | **off** in production (`ce_channel_mode="off"`, `exit_strip_mode="off"`); **SHELVED** 2026-07-30 after all three kills fired | D0 §20 | re-openable only on new external evidence |

---

## Tally — how few true knobs remain

**Production tally (B2, 2026-10-14).** Supersedes the h405 and pre-G4
tallies below, which are kept as the record.

- **TDDFT ridge (Tier 0, 18 Å, m = 21 He, 3.5 Å gate):** a one-parameter
  family $(a, E_\text{bind})(c)$; at B2 $a$ 18.007, $E_\text{bind}$
  0.11727 eV (rows 2–3, Derived on the ridge). Authority covers the
  18 Å fit band 2.54–3.02 Å/ps only; 9 Å = one-sided exit bound
  (c ≤ 5 pass; B2 margin +0.823 Å/ps at m = 21).
- **Sourced (9 + 1 composite):** $s_\rho$ 3.5 Å (gate + density),
  $\lambda_\text{attach}$ (λ₀ 0.9), $\nu$ 2.42, $E_\infty$ (|S| 0.308),
  $E_\text{avail}^\text{ion}$ 2.70, $n^*$ 21, $D_\text{floor}$ 4.97 cm⁻¹,
  $t_\text{detect}$ 8.53 µs, $\rho_\text{He}$ shape. **Composite:**
  $D_0(1)$ = Sourced base 74.4 cm⁻¹ (IHe05) × taper 2.2 → 163.7 cm⁻¹
  (row 18). The taper 2.2 / 1.5 / 1.3 is **experiment-derived** (step
  ratios of the detected I⁺(He)ₙ histogram, not literature; literature
  is against it), and the rq4graded family was selected on the same
  distribution (I86, finc1v725-era law) that later arbitrates B2.
- **Derived (4):** $f_\text{int}$ 0.1221, $t_\times$, $\sum_i D_0$
  (Σ(21) 0.206 eV), $\Pi$.
- **Bounded, pinned inside the band (3):** $s_\text{eff}$ 8,
  $f_\text{ret}$ 0.1 (never swept), the occupancy exponent $p$ 1. Noise
  $T_\text{eff}$ is off (Tier 3).
- **Bounded, selected in preliminary calibration, frozen at G5 (1):**
  τ 6.4 ps — chosen against the same experimental observables at
  c 2.0 / `raw` (G4.14 A5), inherited by B2 without re-scan. External
  band [2.6, 16.5] ps is Na⁺ only (experiment 2.6; GAH25 Table III
  first / second shell 7.3 / 16.5); 6.4 = ×0.88 of GAH25's 7.3. The
  often-quoted 6.55 is the band's geometric mean, not a GAH25 value
  (row 11). Strictly a third experiment-selected quantity, selected at
  a different law.
- **Free — selected by the experiment at G5 (2):** the ridge position
  $c$ 1.6 and $E_\text{int}(0)$ 0.3297 eV (inside RQ1 [0.2, 0.5]).
- **Inactive (5):** $v_c$ and $p_\text{tail}$ (B2 is uncapped), κ (baked
  in at 1), the electronic picture (baked in at `statistical_mixture`),
  the birth margin (0, guard-refused under `boltzmann`).
- **Arms (selections, not fitted):** rows 5a, 25–30 (sampler
  `post_pickup`, G1) and the cross-cutting pipeline / retained-policy
  rows. T5/T6 (rows 27/28) are structurally inert at production.

**Production headline (B2; supersedes the headlines below).** "TDDFT
constrains, the experiment selects" (THESIS_FINAL_GOALS G11.1): Tier 0
fixes a ridge, not a law; the Tier-2 experiment selects two knobs, $c$
on the ridge and $E_\text{int}(0)$. The selected $c$ depends on the
droplet description (raw → 2.0, `post_pickup` → 1.6), and a differently
shaped law (h405) lands comparably — neither TDDFT nor the detected
I⁺(He)ₙ data identify a unique drag function (G5). Remaining misfits
(n₁ / KE₁, n ≥ 15 tail, KE width at fixed n) are outlook
(`THESIS_DISCUSSION.md`). Effective variant B3 departs from the ridge
(row 3) to buy the high tail.

**h405 tally (2026-10-01 → 2026-10-08; kept as record, now the contrast law):**

- **Locked (3) — fitted to TDDFT (Tier-0 Method B, jointly), then
  fixed:** $b$ 2.5154, $a$ 0, $E_\text{bind}$ 0.1168 eV.
- **Sourced (8 + 1 composite):** $\lambda_\text{attach}$ (λ₀ 0.9),
  $\rho_\text{He}$, $\nu$ 2.42, $E_\infty$ (|S| 0.308),
  $E_\text{avail}^\text{ion}$ 2.70, $n^*$ 21, $D_\text{floor}$ 4.97 cm⁻¹,
  $t_\text{detect}$ 8.53 µs; composite $D_0(1)$ 163.7 cm⁻¹.
- **Derived (5):** $g(\text{depth})$ (14.2 Å reuse), $f_\text{int}$ 0.15,
  $t_\times$, $\sum_i D_0$, $\Pi$.
- **Bounded, pinned (3):** $s_\text{eff}$ 8, $f_\text{ret}$ 0.1, $p$ 1.
- **Free — fitted to the Tier-2 observables at G4 (4 + 1 selection):**
  $v_c$ 5.5 Å/ps, τ 4.4 ps, $E_\text{int}(0)$ 0.405 eV, the rq4graded
  ladder-family selection; $p_\text{tail}=-1$.
- **Inactive (3):** κ, the electronic picture, the birth margin.

*h405 headline:* the genuinely fitted set was $v_c$, τ and
$E_\text{int}(0)$ (plus $p_\text{tail}$ and the ladder family), arbitrated
at G4 — the set D0 §17 marks as effective (**E**).

**Historical tally (pre-G4; kept as record):**

- **Locked (3):** $b$, $a$, $E_\text{bind}$ (Tier-0 complete).
- **Sourced (9):** $\lambda_\text{attach}$, $\rho_\text{He}$, $\nu$, $E_\infty$,
  $E_\text{avail}^\text{ion}$, $D_0(1)$, $n^*$, $D_\text{floor}$,
  $t_\text{detect}$ (8.53 µs, added 2026-07-07).
- **Derived (5):** $g(\text{depth})$, $f_\text{int}$ (Derived coordinate
  since 2026-07-09; it swapped roles with $E_\text{int}(0)$ per rows
  14/16), $t_\times$,
  $\sum_i D_0$, $\Pi$ (pickup↔gate order parameter, §6.11). *(Was 6 — $s$
  promoted to Bounded 2026-07-06.)*
- **Bounded (6–7):** $\tau_\text{dissip}$ (RQ9 reclassification pending,
  2026-07-16), $f_\text{ret}$, $E_\text{int}(0)$,
  $s_\text{eff}$ (promoted 2026-07-06, band ≈[5, 20], staircase landing
  [8, 12]), $v_c$ (the realized `capped_cubic` cap, 2026-07-16; was the
  contingent $v_\text{ceiling}$ — Derived at the §I.10 T4 winner; **the
  re-pilot + v_c sensitivity ring (§4w/I81, 2026-07-21) locate the c1
  joint-landing basin at $v_c\approx7.25$–$7.5$** at τ3.2/E₀0.27;
  **superseded at the corrected geometry: the basin is 5.5–6.0**),
  the birth margin $m\in\{3, 4.67, 6\}$ Å (row 25, Slice T7 2026-07-16;
  active only under the `uniform_volume` twin-parity arm),
  ($T_\text{eff}$/noise). *Plus one Free choice:* the tail exponent
  $p_\text{tail}\in\{0,-1\}$ (row 4b; **pinned $p_\text{tail}=-1$ for the
  c1 finalist** — the E-2 twin scan §4u/I78 showed steeper
  $p\in\{-2,-3\}$ re-couples the KE↔histogram axes, 2026-07-21).
- **Free (2):** the ladder **steepness** $\kappa$ (Form U sigmoid, $D_0(n{>}1)$)
  and the **electronic picture** (3 options — mix/$X_2$/cooling_relaxed — but
  one selection) — both choices the Tier-2 size distribution
  arbiters, co-fit (not separable). *(2026-07-06: both shown **near-flat on the
  9 Å staircase** — κ inverted/capped, picture magnitude-degenerate — so the
  size distribution is their only realistic anchor; the Stage-1 campaign
  co-fit re-scopes to $s_\text{eff}$×τ, Phase F plan F2 NB.)* *Conditional:*
  the occupancy exponent $p$ (A12) is **tied to $\kappa$ by default** (0 net
  new free knobs); it becomes a 3rd free knob only if the size-distribution
  first-shell cutoff forces the $p\!\leftrightarrow\!\kappa$ split.

*Historical headline (pre-G4):* after sourcing/derivation/bounding, the
genuinely free calibration collapsed to **the ladder shape and the
electronic picture**, both discriminated by the Tier-2 size distribution,
with the early-window scalars bounded and cross-checked rather than fit.
*(2026-07-06 amendment: the dominant cascade-magnitude lever became the
Bounded $s_\text{eff}$, co-anchored with τ by the 9 Å staircase via
first-shed timing.)*

## Anchor coverage by tier

- **Tier 0:** the drag ridge $(a, E_\text{bind})(c)$ at the 3.5 Å gate,
  18 Å fit + 9 Å exit bound (done).
- **Tier 1 (retired, historical):** drag/mass-scenario sensitivity via the
  shell trajectory (21→19→14) — touched $\lambda_\text{attach}$, $g$,
  $v_\text{ceiling}$, the early-window tolerance. Results not used in the
  thesis (G11.3).
- **Tier 2 (load-bearing):** selects $c$ and $E_\text{int}(0)$ (G5); the
  sampler sensitivity and the pinned / sourced knobs ($\lambda_\text{attach}$,
  $f_\text{ret}$, $\nu$, $s$, τ, ladder family) are reported, not re-fit.
  **Single observable carrying the calibration** — its discriminating
  power is the project's central dependency.
- **Tier 3:** the noise model (open).
- **External / author-contact:** OQ1 (drag electronic state), OQ2
  ($KE_\text{shed}$ — **resolved 2026-07-10, ε≈0**), OQ3 ($E_\text{bind}(N)$ vs ladder — **RQ4 sharpened target: the $D_0(1{:}2{:}3)$ ratios**), OQ4 ($S_{\mathrm{I^+}}$/$n^*$
  provenance), OQ5 (cascade timing), OQ6 ($E_\infty$ stripping reach), OQ7
  ($S_{\mathrm{I^+}}$ energy reference — interprets the $|S|$-vs-$\sum D_0$ gap);
  M1 (a sourced I⁺–He D₀; gates the KE₁ deficit, `TIER2_MASS_SCENARIOS.md`).

### Flagged identifiability item — $E_\text{bind}$'s provenance spread exceeds the $\bar n$ gate (2026-08-10; h405-era record)

> **B2 note (2026-10-14).** Measured on h405 (cubic, 14.2 Å gate). Under
> B2 the well is no longer an independent co-extraction: the ridge ties
> $E_\text{bind}$ to $c$ (row 3; G4.12), and the n̄ gate at G5 is met by
> re-placing $E_\text{int}(0)$ per $c$. The numbers below are kept as the
> record of the cubic-era argument.

Measured by the atlas §6.5 E_bind study (twin Step 2 + MD Step 3; D0
§9.2/§9.3). $E_\text{bind}$ is **row 3** above — then classed **Locked**
(Tier-0, co-fit with the drag pair); D0 §17 carries the same quantity as
*Derived / P (paired)* with the note "pairing is the constraint, not the
value" (D0 now reads **Locked**, harmonized 2026-10-01). Both framings are too weak on one axis:

- The measured lever is $\partial\bar n/\partial E_\text{bind} \approx
  \mathbf{-7.3\ \text{eV}^{-1}}$, agreeing across **three independent
  systems** — MD $-7.25$, twin $-7.36$, and the §6.7 item-2 lq scan
  ($-0.50$ per the $0.0686$ eV step).
- The Tier-0 **co-extracted** well spread is $0.048$–$0.154$ eV
  ($0.106$ eV wide; every value is some Tier-0 form's own joint fit).
  That maps to $\Delta\bar n \approx \mathbf{0.77}$.
- The MD-side $\bar n$ hard gate is $[3.77, 4.37]$ — **$0.60$ wide**.

So **the extraction uncertainty on $E_\text{bind}$ alone spans $1.3\times$
the entire $\bar n$ gate band**. Concretely, at the corrected geometry
h405 measures $\bar n = 3.748$ at the bundle well — $0.022$ *below* the
gate floor (the §6.1 ring recorded the same shortfall) — and
$\bar n = 4.245$, comfortably mid-band, at the shallow Tier-0 well, with
trap dropping $0.063 \to 0.012$ at the same time.

**Consequence for the calibration map.** $E_\text{bind}$ is functionally
an $\bar n$ knob with a KE side-effect: per eV it moves $\bar n$ roughly
**30× more efficiently** than it moves $\text{KE}_1$ (the KE axis is
closed — the whole well span buys $\le 0.064$ eV, §9.2). Any $\bar n$-gated
arbitration is therefore conditional on *which* Tier-0 well its drag
partner was co-extracted with, not merely on the pairing being
self-consistent. This does **not** license moving $E_\text{bind}$ (§9.1
bounds the physically defensible variation at $\le 0.009$ eV); it means
an $\bar n$ verdict inside $\pm 0.4$ of a band edge should be read as
*not resolved by this observable* until the pairing is stated.

*Note:* Tier 2 carries a heavy load (8+ quantities on one observable). The
identifiability arguments are documented per-parameter (e.g. $f_\text{ret}$ via
9/18 Å contrast, $\tau$ via post-crossing tail, $s$ jointly with ladder shape) —
but whether the size distribution actually separates all of them was the
key open empirical question of the calibration campaign; the thesis answer
is that it does not identify the drag law (G5). *(2026-07-06: first
load shed — the 9 Å anchored staircase (TDDFT prior) co-anchors
$s_\text{eff}$×τ, separated by first-shed timing; κ and picture
near-flat on the staircase, so their weight rests entirely on the size
distribution.)*
