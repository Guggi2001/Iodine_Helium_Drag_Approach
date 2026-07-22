# Tier 2 — Research Questions (literature-research phase entry document)

> **Purpose.** Created 2026-07-09 (user decision, post-Wave-7 discussions)
> as the consolidated register of every open mechanism-level physics
> question in the biphasic Tier-2 model, stated so that a literature
> search / domain-expert consultation can attack each one directly. The
> probe program (Waves 1–7, `TIER2_STAIRCASE_PROBE_FINDINGS.md`) exhausted
> the calibration space; what remains open are **energy-bookkeeping
> conventions**, not knob values. This document is the working target of
> the cross-validation phase: findings and provenance get recorded here
> (per-RQ NBs), the user adjudicates, and any resulting model change gets
> its own design document (dimensional analysis, interchangeable enum arm)
> behind the `[PROCEED TO IMPLEMENTATION]` trigger.
>
> **Stance (program-wide, unchanged):** document-only; TDDFT is not ground
> truth; the arbitration observable is the experimental terminal I⁺Heₙ
> size distribution (`data/reference/integrated_i_he_abundance.csv`) read
> at the Sourced detection time t_detect = 8.53 µs.

---

## 1. The model in brief (for a domain reader)

I₂ inside a helium nanodroplet is ionized; the Coulomb explosion ejects
two I⁺, each dressed with a He shell (n ≈ 21 at the 9 Å condition). The
MD carries: TDDFT-calibrated drag and Poisson He pickup, both gated on the
local He density ρ_He (shared erf-complement bubble boundary); an internal
energy reservoir `E_int` per ion with

- **S2 onset:** `E_int(0) = f_int·E_avail` deposited once at ionization,
- **S1 pickup heating:** `+f_ret·D₀(n)` per captured He,
- **K1 evaporation:** RRK rate `k = ν·(1 − D₀(n)/E_int)^(s−1)` active in
  the band `D₀(n) < E_int < Σ(n)`, each shed draining exactly `D₀(n)`,
- **K2 Newton cooling:** `−E_int/τ`, density-gated on the production-
  intended arm (`cooling_spatial_gate = density_scaled`, i.e. zero in
  vacuum),

on a Form-U dissociation ladder `D₀(n)` (flat bottom ≈ 9.22 meV at
mixture/κ = 1; Σ(21) ≈ 0.188 eV), with `s` either the classical `3n−3` or
the promoted constant `s_eff` (Bounded, landing band arm-dependent ≈ 8–30).
An exact event-driven detection stage continues the post-ejection cascade
to the detector (8.53 µs).

**Where the probe program left the model** (findings §4b–§4e): the
post-ejection cascade is an exactly closed system (the self-bound margin
`G = E_int − Σ(n)` is shed-invariant), sits in the Klots
evaporative-ensemble regime at the detector (descent ≈ 0.4–0.9 He per time
decade), and every gated observable organizes along the race margin
Δ× = t_eject − t× between gate-open and droplet ejection. The kinetic and
race structure is internally validated; the remaining gaps are the
conventions below.

## 2. The problem, stated plainly

Three facts localize the open physics:

1. **The two largest experimental bins are structurally unreachable from
   the cascade side.** Bare I⁺ (43.5 %) requires an exactly-zero in-bubble
   leak (measure zero, by G-invariance); I⁺He (17.5 %) requires a sub-rung
   leak *and* sits behind the collapsing RRK bracket (absent in 10 200
   simulated detections).
2. **The gated cascade never self-terminates** (an eternal
   `time_exhausted` class at the detector) — because a shed drains only
   `D₀(n)`, evaporation is exactly self-sustaining by construction.
3. **The initial reservoir has no defensible provenance.** The
   `f_int·E_avail` parametrization is indefensible at the values that make
   the mechanism work (no actual energy partition is implemented) and
   inert at the values that are physically defensible as a literal
   Coulomb-energy fraction (~1 %).

Each maps onto a research question below. RQ1–RQ3 live in one
energy-bookkeeping layer (the ionization onset, each shed, the
self-unbound boundary); RQ4 is the ladder itself; RQ5 the flight; RQ6 a
cross-check of an already-adjudicated convention; RQ7 (added post-Wave-8)
is the one *kinematics* gap — the production budget currently changes
only the bookkeeping, not the mechanics.

## 3. The research questions

### RQ1 — Provenance and magnitude of E_int(0) (= CALIBRATION_MAP OQ2, fired 2026-07-09)

- **Model convention:** `E_int(0) = f_int·E_avail` at t = 0; nothing is
  subtracted from the fragment mechanics (`e_int_onset_eV`,
  `ion_initial_state.py`).
- **Problem:** at the working values (f_int ≈ 0.24–0.65) the model books
  0.2–1.35 eV of shell heat with no mechanical source; at the physically
  defensible literal coupling (~1 % → 8–27 meV) the mechanism goes inert.
  The row-14 scenario-keyed floors multiply out to the same absolute
  energy (≈ Σ(21) at both budgets) — the natural variable is an
  **absolute `E_int(0)` [eV]**, not a fraction.
- **Working hypothesis (NOT adjudicated):** E_int(0) ≈ 0.2–0.5 eV,
  essentially budget-independent, sourced from (i) solvation
  reorganization on vertical ionization (≤ ~Σ(21) ≈ 0.19 eV) and (ii)
  electronic / spin–orbit relaxation of nascent I⁺ partially degrading
  into the shell (fine-structure scale ~0.7–0.9 eV; ¹D ~1.7 eV — this
  would tie the electronic-**picture** knob to the E_int(0) provenance),
  with retained drag heating (tens of meV) and literal KER coupling
  (~1 %) as minor additive terms.
- **Literature needed:** I⁺ fine-structure level energies, populations
  after strong-field/photo-ionization, and relaxation pathways/timescales
  in a He environment; solvation reorganization energies on dopant
  ionization in He droplets; energy partition measurements for Coulomb
  explosion of dopants in droplets (KER deficits).
- **Discriminators:** VMI KER (an actual partition would show a KER
  deficit ∝ f_int); budget-(in)dependence of the bare peak; any measured
  picture ↔ E_int(0) coupling.
- **If confirmed:** f_int reclassifies fraction → absolute Bounded energy
  [eV]; row-14 scenario keying dissolves; the B.1(2) production
  prediction inverts (t× budget-independent; production changes via the
  earlier t_eject).

#### RQ1 findings — NB register (deep-research run, 2026-07-09)

> **Method/provenance of this register.** Multi-agent literature search
> (5 angles: I⁺ electronic-state populations; electronic relaxation of
> embedded ions in He; sudden-solvation/snowball reorganization energy;
> CE KER deficits in droplets; electronic-to-shell coupling). 22 sources
> fetched (all primary literature + NIST), 90 claims extracted, top 25
> adversarially verified by independent 3-voter panels (20 confirmed 3-0,
> 5 refuted 0-3, 0 unverified), synthesis pass merged to 11 findings.
> **Status: recorded, NOT adjudicated.** No source measures the exact
> target system (I₂ CE *inside* a droplet → two recoiling I⁺Heₙ); all
> values transfer by analogy from (a) gas-phase strong-field I₂ CE,
> (b) alkali cations in He droplets, (c) *neutral* alkyl-iodide
> photodissociation in droplets. Neutral-vs-ion is the sharpest transfer
> risk.

**Confirmed NBs (vote 3-0 unless noted):**

- **NB-RQ1-1 (electronic provenance — mechanism).** Strong-field CE of
  I₂ populates excited I⁺ fine-structure/¹D₂ states *at the expense of
  KER*: the (1,1) channel KER distribution shows two lower-KER shoulders
  whose splittings are consistent with I⁺ ³P/¹D₂ intervals. Forbes et
  al., *J. Phys. Chem. A* 2022 (PMC9706571; gas phase, 800 nm, VMI).
  Caveats: assignment explicitly non-conclusive (needs I₂²⁺ PECs); **no
  branching fractions quantified**; gas-phase, not droplet-embedded.
- **NB-RQ1-2 (level energies — anchored).** I⁺ (I II 5p⁴) levels above
  ³P₂ ground: **³P₀ = 0.799 eV, ³P₁ = 0.879 eV, ¹D₂ = 1.702 eV** (NIST
  ASD; independently via PES Rydberg limits 11.33/12.15 eV minus
  IE = 10.451 eV, PubMed 20815564). Confirms the working hypothesis'
  0.7–0.9 eV spin-orbit scale; note the *first* excited level is ³P₀.
- **NB-RQ1-3 (solvation reorganization — measured benchmark).** Sudden
  Na⁺ creation in a He droplet releases a measured dissipated solvation
  energy **E_disp(∞) ≈ 224 meV, droplet-size-independent**
  (3600–9000 He); He-DFT total ≈ 384 meV; single-ion energy lowering
  0.40 eV (García-Alfonso 2024). Albrechtsen et al. 2025
  (arXiv:2502.11783, companion to *Nature* 623, 319 (2023)). Caveat:
  Na⁺ created at the droplet *surface*, not an embedded halogen ion.
  **Amends the working hypothesis:** the solvation term's ceiling is
  ~0.22–0.40 eV, above the assumed ≤ Σ(21) ≈ 0.19 eV.
- **NB-RQ1-4 (dissipation law — external precedent for K2).** The
  released solvation energy leaves by He-atom ejection and **follows
  Newton's law of cooling for the first ~5 ps**; dissipation rate is
  roughly ion-independent across the alkali series; ~half the solvation
  energy is gone by ~4 ps (Albrechtsen 2025; time-resolved alkali TDDFT,
  *J. Chem. Phys.* 160, 164308 (2024)). Cross-link: direct literature
  precedent for the K2 Newton-cooling term and for reorganization energy
  *leaving* the shell rather than staying resident.
- **NB-RQ1-5 (terminal size distributions are evaporative).** HeₙX⁺
  snowball size distributions are evaporative-dissociation products of
  the ionization excess energy; abundance ≈ ∝ per-atom dissociation
  energy Dₙ (Klots/Hansen evaporative ensemble). Bartl/Scheier/Echt,
  *J. Phys. Chem. A* 2013 (10.1021/jp406540p); *Int. Rev. Phys. Chem.*
  2020 review. Cross-link: underwrites the Tier-2 arbitration observable
  itself and the RQ4 ladder read.
- **NB-RQ1-6 (first-shell size).** First He shells around heavy cations
  are large: 18–20 (Mg⁺), **20 (I₂⁺)**, 60/62 (C₆₀⁺/C₇₀⁺) — bracketing
  the n* ≈ 21 assumed for I⁺ (same JPCA 2013 source). Caveat: I₂⁺, not
  atomic I⁺.
- **NB-RQ1-7 (KER deficit + dressed escape — precedent).** Fragments
  born inside a droplet show speed/KE distributions "notably modified"
  vs gas phase and **escape still dressed in finite Heₙ** (IHeₙ,
  CH₃Heₙ). Braun & Drabbels, *J. Chem. Phys.* 127, 114303 & 114304
  (2007) (266 nm alkyl-iodide photodissociation, ion imaging). Caveat:
  *neutral* photofragments.
- **NB-RQ1-8 (against full budget-independence — kinematic coupling).**
  Terminal Heₙ size correlates with fragment speed ("dynamical
  adjustment of the solvation structure size to the relative speed") —
  shell size/energy is *kinetically coupled* to the velocity history
  (Braun & Drabbels paper II). Argues against treating E_int(0)/shell
  state as fully independent of the KER budget.
- **NB-RQ1-9 (energy-gated boil-off — precedent).** High-internal-energy
  fragments (C₂H₅, CF₃) exit with **no attached He** — internal energy
  rapidly destroys nascent complexes (same source). Supports the
  energy-gated evaporation picture; caveat: absence-based inference,
  neutral vibrational energy vs I⁺ electronic energy.
- **NB-RQ1-10 (microphysics of shell heating).** TDDFT of halogen-dimer
  photodissociation in He: two-regime **ICVF mechanism** — perfectly
  inelastic collision with the shell/cavity wall in the first
  ~0.08–0.16 ps, then viscous-flow friction (González group, *PCCP*
  2016, 10.1039/C6CP04315A; "probably general character"). NOTE: the
  quantitative sub-claim (effective He mass M = 10.1 amu ⇒ ~14 % of
  fragment KE deposited) was **refuted 0-3** — mechanism yes, numbers no.
- **NB-RQ1-11 (surface limiting case — KER coupling small).** Alkali
  dimers Coulomb-exploded on the droplet *surface*: fragment KER peak
  centers match gas-phase Coulomb energies from equilibrium distances
  (broadened, not shifted) — no large (≫100 meV) KER deficit in the
  surface geometry. *Phys. Rev. A* 107, 023104 (2023) +
  arXiv:2401.09211. NOTE: the companion "~1 % small-loss" TDDFT
  sub-claim was **refuted 0-3** — the small-coupling *magnitude* is not
  endorsed, only the surface-case absence of a large shift.

**Refuted claims (0-3, recorded so they are not re-imported):**

1. Per-He binding-scale summation bound (D∞ = 0.616 meV /
   He₁₂Ar⁺ D₁₂ = 17.8 meV ⇒ shell sum ≲ few 100 meV) — refuted as a
   bound for E_int(0).
2. "Heavy cations (Au⁺, I⁺, Cs⁺) show He-tagged signal exceeding bare
   ion" — refuted; do **not** cite toward the RQ3 bare-peak question.
3. "KE-loss fraction depends strongly on fragment mass" (as extracted) —
   refuted.
4. Inelastic effective-mass estimate M = 10.1 amu ⇒ ~14 % per-fragment
   KE into shell (~0.06 eV at 0.80 eV / ~0.19 eV at 2.70 eV) — refuted.
5. TDDFT "~1 %-scale KER-to-shell coupling" for surface alkalis —
   refuted.

Consequence of (4)+(5): **both** quantitative estimates of the
KER-coupled channel died in verification — its magnitude is genuinely
open in the literature.

**Unverified leads (extracted, ranked below the top-25 verification
cut — NOT verified, follow-up candidates for the relaxation-timescale
gap):**

- Droplet-induced electronic relaxation/spin quenching timescales
  600–1100 fs (*PCCP* 2022, 10.1039/D2CP03335F).
- ">1 eV of electronic energy dissipated into nuclear/nanofluid dof in
  < 1 ps" in He droplets (*Nat. Commun.* 10, 5735 (2019) family).
- Rb/Cs desorbing atoms populate lower electronic states via
  droplet-induced relaxation (arXiv:1406.4713).
- **Counter-lead:** Ba⁺ in liquid He remains optically active (LIF
  observed) — excited-cation electronic energy is *not* necessarily
  fully quenched non-radiatively into the bath (AIP Advances 8, 015328
  (2018)). Directly tempers the assumption that ³P₀/³P₁/¹D₂ energy must
  land in E_int; couples to RQ5 (radiative channels in flight).

**Synthesis (for adjudication, not adjudicated):**

- **Magnitude:** E_int(0) ≈ 0.2–0.5 eV is *defensible* — but as a
  partial capture from two distinct reservoirs, not from KER coupling:
  (i) electronic/spin-orbit (0.80/0.88 eV; ¹D₂ 1.70 eV) with
  droplet-embedded branching fractions unmeasured (NB-1/2, leads), and
  (ii) solvation reorganization ~0.22–0.40 eV (NB-3).
- **Budget-independence:** holds for the solvation and electronic
  components (NB-2/3); *contradicted* for the kinematic shell-heating
  component (NB-8), whose magnitude is unquantified (both numerical
  estimates refuted). Best-supported model form:
  `E_int(0) = E_solv + E_elec (budget-independent) + E_kin-coupling
  (budget-dependent, magnitude open)`.
- **Dominance ranking (as hypothesized, now sourced):** electronic
  channel ≳ solvation reorganization ≫ literal KER coupling — with the
  caveats that the electronic branching is unquantified in droplets and
  the Ba⁺ counter-lead allows radiative escape of part of it.
- **Supports the RQ1 "if confirmed" consequence:** reclassifying
  `f_int` fraction → absolute Bounded `E_int(0)` [eV] is the natural
  variable per NB-3 (absolute, size-independent) + NB-2 (absolute atomic
  scales).
- **Open (from the run):** (1) droplet-embedded I⁺ fine-structure/¹D₂
  branching fractions; (2) timescale/efficiency of electronic→shell
  degradation for I⁺ in He (vs radiative escape); (3) true magnitude of
  the KER-coupled term for a 127 amu fragment; (4) whether the embedded-
  I⁺ ~21-He-shell reorganization energy matches or exceeds the
  0.22–0.40 eV Na⁺ benchmark.

#### RQ1 adjudication — 2026-07-09 (user decision)

**Decision (document-level convention; no code touched in this phase):**

1. **Onset convention retained, no longer naive.** `E_int(0)` is
   deposited once at t = 0 and nothing is subtracted from the fragment
   mechanics — now *sourced* rather than assumed: the dominant channels
   (electronic relaxation, NB-RQ1-1/2; solvation reorganization,
   NB-RQ1-3) do not draw on the Coulomb/KER budget, so no-subtraction is
   the physically correct bookkeeping for them, not a simplification.
   The only channel that would justify a KER subtraction — the kinematic
   leak (NB-RQ1-8) — has genuinely open magnitude (both quantitative
   estimates refuted 0-3) and stays **un-modeled, absorbed by the
   sweep**.
2. **Sweep variable reinterpreted: absolute, not fractional.** The
   physical variable is the absolute
   `E_int(0) [eV] = f_int · E_avail`. `f_int` remains the config knob
   and sweep coordinate, but sweeps are *specified, reported, and
   transferred in absolute eV*. Budget transfer 0.80 → 2.70 eV holds
   `E_int(0)` [eV] fixed — i.e. `f_int` rescales by 0.80/2.70 ≈ 0.296 —
   **not** `f_int` fixed (holding the fraction would silently triple the
   reservoir at production, contradicting the sourced
   budget-independence of both dominant channels).
3. **Consequences (the RQ1 "if confirmed" branch is taken).** `f_int`
   reclassifies fraction → absolute **Bounded** energy [eV]; row-14
   scenario keying dissolves; the B.1(2) production prediction inverts
   (t× budget-independent; production behavior changes via the earlier
   t_eject). Follow-through edits done (2026-07-09): `CALIBRATION_MAP.md`
   rows 14/16 swapped (row 16 `E_int(0)` Derived → Bounded absolute,
   `f_int` → Derived coordinate; update note added) and the decision
   entry recorded in `drag_migration_log_tier2.md`.
4. **Sourced sweep band (prior, not a constraint).** Floor ≈ 0.22 eV
   (solvation alone, Na⁺ benchmark); working band **0.2–0.5 eV**; tail
   to ~0.8–1.0 eV (full ³P₀/³P₁ deposit) admissible, ~1.7 eV (¹D₂)
   fringe — the tail discounted by the two open unknowns
   (droplet-embedded branching fractions; deposit-vs-radiate efficiency,
   Ba⁺ counter-lead → RQ5).
5. **Timing idealization stated.** Single t = 0 deposit + K2 drain is
   coherent with the Na⁺ precedent (NB-RQ1-4: prompt deposit, ~ps
   Newton-type dissipation) for the solvation component. Lumping the
   electronic channel into the t = 0 deposit is the acknowledged
   idealization — its relaxation timescale in He is unmeasured; any
   in-flight remainder is RQ5's file.
6. **Discriminator refined.** Under this convention no smooth VMI KER
   deficit ∝ f_int is expected (nothing is subtracted). The electronic
   channel instead predicts *structured* fine-structure satellites — a
   sub-population down-shifted by ~0.8 eV (Forbes-type shoulders).
   Falsifiable both ways: smooth deficit observed → the no-partition
   convention is wrong; satellites observed → electronic channel
   confirmed, its weight directly readable from the VMI.

**Status: RQ1 closed at the convention level** (provenance understood;
variable reclassified fraction → absolute). The magnitude remains
Bounded and swept. Synthesis open items (1)–(4) stay live as literature
gaps; item (2)'s in-flight part migrates to RQ5.

> **NB (post-Wave-8, 2026-07-09):** consequence 3's "production behavior
> changes via the earlier t_eject" assumed the budget moves the
> mechanics; in the delivered model it does not (the budget's single
> physics reader is the S2 deposit — Wave 8, findings I28). The earlier
> production ejection is *physics still to be modeled*, not current model
> behavior → RQ7. The absolute-eV sweep convention itself is unaffected
> (indeed strengthened: the delivered model is exactly budget-invariant
> in absolute E_int(0)).

> **NB (Wave 9, 2026-07-09) — the histogram-implied p(E₀) lands in the RQ1
> band, no branching signature.** The Wave-9 E₀-mixture inversion (findings
> §4g, I29–I31) inverts the experimental I⁺Heₙ histogram into an implied
> p(E₀): it fits to Wasserstein-1 = 0.086 bins with mass **concentrated
> at/above E\* ≈ 0.46 eV** (43.6 % bare-class) and an in-band declining tail
> through ≈ 0.28 eV — squarely inside this section's sourced 0.2–0.5 eV band
> **without tuning** (the RQ1-required and histogram-required ranges
> coincide — the strongest quantitative support the biphasic reservoir
> magnitude has received). Two RQ1-relevant reads: (i) the implied density
> is **unimodal with no genuine second mode** → the fit does **not** demand
> the electronic fine-structure branching of consequence 6 (a branching
> signature would need multi-modal p(E₀) surviving basis de-collinearization;
> this one does not) — consistent with but not proof of the smooth
> reorganization + electronic-relaxation picture; (ii) a **≥ 9.6 %** ensemble
> tail must sit **below** the 0.22 eV solvation floor (deep-shell bins) — not
> an E_int(0)-provenance demand but the **droplet-radius axis's** job
> (larger K, deeper leak at in-band E₀; findings §4f discussion, W9-P2).
> Conditional on RQ3 spec (b) throughout.

> **NB (Wave 10 Steps 1+1.5, 2026-07-10/11) — the E₀-width demand dissolves;
> narrow E₀ is supported, conditional on RQ7.** The Wave-10 zero-MD analysis
> (findings §4h, I33–I36) shows the Wave-9 broad p(E₀) was a **reduction
> artifact of the pinned droplet/K axis**: a single sharp
> **E₀ = 0.28 eV** — one value inside this section's sourced solvation band —
> plus the physical Kornilov droplet-size distribution reproduces the
> experimental histogram untuned (bare 42.6 % vs 43.5 %, W₁ ≈ 1 bin) **at the
> production K-scale bracket (K₀ ≈ 0.49)**, while being excluded ~19× at the
> pinned 9 Å probe kinematics (W10-P3/P5). Consequences for RQ1: (i) the
> super-cliff (> 0.46 eV) and sub-floor (< 0.22 eV) E₀ mass Wave 9 demanded
> are re-carried by the small-K and large-K ends of the droplet axis — **no
> super-solvation or kinematic-coupling E₀ source is required on current
> evidence**; (ii) the Bounded E_int(0) sweep reduces from a distribution
> question to a **scalar calibration in the solvation band**, augmented by a
> physical **tens-of-meV solvation smear** (adjudicated useful, 2026-07-11 —
> sourced by the same NB-RQ1-3 reorganization physics; to be carried as a
> sensitivity leg of the droplet slice); (iii) with RQ7 resolved and the VMI
> speed check pinning K₀, **the bare fraction reads E₀ nearly directly** —
> RQ1's magnitude becomes measurable rather than swept. All conditional on
> the RQ7 kinematics measurement (K₂.₇₀); if K₂.₇₀ lands near the probe 0.9
> instead of ~0.5, this NB inverts and the super-solvation route reopens.

### RQ2 — Per-shed kinetic-energy release ε (= findings OQ-F)

- **Model convention:** a shed drains exactly `D₀(n)`
  (`internal_energy_budget.dE_int_shed_eV`); the evaporated He leaves with
  zero translational release from `E_int` (the cold-shed kick is
  mechanical bookkeeping only). Hence `G` is shed-invariant and the
  cascade is exactly self-sustaining.
- **Problem:** a statistical evaporation releases translational energy
  ε ~ (E_int − D₀)/s ≈ 5–20 meV per shed — over a ~17-shed cascade,
  ~0.1–0.3 eV, the same order as the entire Σ(21) budget. First-order,
  not a correction.
- **Literature needed:** evaporative-ensemble theory (Klots) KER
  conventions; measured/computed kinetic-energy-release distributions for
  atom evaporation from cold ionic (He) clusters.
- **Discriminators:** existence of the eternal `time_exhausted` class
  (with ε, cascades self-terminate and terminals converge); the endpoint
  distribution of fragmenting suppressed complexes (couples RQ3); the
  n = 1 bin.
- **NB (deep-research, 2026-07-10):** the ε literature register is joint
  with RQ3 (one coupled question) — see the **RQ2+RQ3 joint findings**
  block under RQ3. Headline: theory prescribes ε ≈ c·D/G (Klots–Hansen,
  c ≈ 1.5, G ≈ 23.5) ≈ D/16 ≈ 0.5–0.6 meV — nonzero (**ε = 0 is not the
  prescription**) but ~10–30× **below** the feared 5–20 meV; the ε ~ D scale
  is physically ruled out (ns-scale decay vs µs survival).
  **(ADJUDICATED 2026-07-10: ε ≈ 0 — small ε neglected; a shed drains exactly
  D₀(n). See the RQ2/RQ3 adjudication block under RQ3.)**

### RQ3 — Fate of a net self-unbound complex (= findings OQ-B, fragmentation hypothesis)

> **Status update (2026-07-11, two-channel adjudication — log entry same
> day): partially UN-ANCHORED.** The experimental bare-I⁺ KE (~3 eV,
> channel-distinct from the snowballs) removes the experimental anchor
> for "suppressed → bare sources the 43.5 % peak": an *evaporative* bare
> class would arrive slow (~0.2–0.3 eV after drag). The adopted
> sequential-shed physics stands (Stapelfeldt), but its **marginal
> self-termination-at-small-n sub-case** gains weight (a candidate n = 1
> feeder), and the bare bin is reinterpreted as channel branching →
> RQ8. Discriminator: a slow shoulder (or its absence) in the bare-KE
> distribution.

- **Model convention:** `E_int > Σ(n)` ⇒ evaporation suppressed; the
  complex rides to the detector intact at its handover n (gated dead arms:
  100 % `suppressed` at n = 21).
- **Problem/hypothesis:** a complex carrying more internal energy than its
  total shell binding plausibly fragments over 8.5 µs. If it does, **bare
  I⁺ is the crossed-after-ejection side of the Δ× race** and the
  experimental bimodality is the cliff structure itself. Quantitative
  structure: gateless sequential boil-off *without* ε goes exactly to bare
  (G > 0 invariant; the n = 1 direct channel has no barrier); *with* ε
  (RQ2), the one-sided G₀ distribution yields bare (bulk) + a decreasing
  small-n tail including n = 1 (fringe) — predicting a budget-dependent
  bare fraction. **RQ2 and RQ3 are one coupled question.**
- **Literature needed:** sequential evaporation vs (multi)fragmentation of
  clusters excited far above total binding; timescales for shell blow-off
  of ionic cores in He.
- **Discriminators:** the 43.5 % bare peak; the 17.5 % n = 1 bin; bare
  fraction vs budget (0.80 vs 2.70 eV — **RQ7-gated**: no in-model budget
  mechanism exists until production kinematics are modeled).
- **NB (post-Wave-8, 2026-07-09):** the weight-level side is measured
  (findings §4f, I26/I27): the suppressed fraction at a knob point is a
  step 0 → 1 at E* = 0.4612 eV (no 43.5 % split without injected ensemble
  heterogeneity), and the **sub-rung sliver below E* expresses n = 1**
  through the ordinary cascade (all ions at n_detect = 1, s_eff = 8) — a
  third n = 1 route besides the RQ2-ε fringe and the RQ4 deep rung,
  active within one rung of the cliff in leak units (= 22.6 meV wide in
  E_int(0), D₀(1)·e^K).

#### RQ2+RQ3 joint findings — NB register (deep-research run, 2026-07-10)

> **Method/provenance of this register.** Multi-agent literature search
> (5 angles: evaporative-ensemble/Klots KER prescription; measured/computed
> He KER distributions; He-nanodroplet snowball size distributions;
> sequential evaporation vs multifragmentation of over-energized clusters;
> incremental He binding around heavy cations). 17 sources fetched (all
> reported findings rest on primary peer-reviewed papers — Hansen, Klots,
> Märk/Echt/Scheier groups; unreliable hits filtered), 39 claims extracted,
> top 25 adversarially verified by independent 3-voter panels (20 confirmed,
> 5 refuted 0-3, 0 unverified), synthesis pass merged to 7 findings.
> **Status: recorded; ADJUDICATED 2026-07-10** (ε ≈ 0 — small ε neglected;
> spec (a) rejected → sequential shed-to-self-termination; RQ4 promoted — see
> the primary-source grounding + adjudication block below). RQ2 and RQ3 are
> researched jointly
> because ε (RQ2) sets the RQ3 fragmentation endpoint. **Transfer risk is the
> dominant caveat:** no source measures the target system (sub-Kelvin
> halogen-cation I⁺Heₙ); every value transfers by analogy from noble-gas
> cations (Ar/Kr/Xe), alkali snowballs (Cs⁺/Rb⁺/K⁺/Na⁺), metal anions (Na⁻),
> or pure Heₙ⁺.

**Confirmed NBs (vote 3-0 unless noted):**

- **NB-RQ23-1 (RQ2 ε magnitude — the headline).** Evaporative-ensemble
  (Klots/Hansen) theory sets the per-shed kinetic-energy release to
  `ε ≈ c·D/G` with c ≈ 1.5–2 and the Gspann parameter G ≈ 23.5 (±1) for a
  ~10 µs timescale → **ε ≈ D/12 to D/16 ≈ 0.5–0.6 meV for D ≈ 9 meV**. c = 2
  decomposes as 3/2 (mean KE of the escaping atom) + 1/2 (speed-biased
  outgoing flux); an ion-induced-dipole (Langevin r⁻⁴) attraction cancels
  the +1/2 term → c → 1.5. **Two consequences:** ε = 0 (drain D₀ only) is
  *not* the theory's prescription; but the model's feared ε ~ 5–20 meV (≈ D)
  is physically ruled out — it implies T ≈ 30–115 K → G ~ 1–4 → ns-scale
  evaporation, incompatible with 8.5 µs survival to the detector. Source:
  Hansen, "Do we know the value of the Gspann parameter?" (RG 222428878).
- **NB-RQ23-2 (RQ2 finite-heat-bath correction).** For a few-mode
  sub-Kelvin He shell (small heat capacity C_v ≈ (3N−7)k_B), a
  finite-heat-bath correction reduces ε *below* c·D/G and can drive it
  slightly negative: `ε = c(D/G − D/2C_v)`, daughter temperature
  `T_d = T_e − D/2C_v` (one reported O₂ datum is negative). The correction is
  largest exactly in the few-active-mode regime the target shell occupies →
  the correct ε is even smaller than the naive D/16. Same source.
- **NB-RQ23-3 (RQ2 foundation — ε ≠ 0 is baked into the method).** KER from
  unimolecular evaporation is *the* standard observable used to back out
  per-atom binding energies via finite-heat-bath / evaporative-ensemble
  theory — a framework that presupposes ε ≠ 0 statistically coupled to D,
  structurally undercutting the ε = 0 convention. Sources:
  Parajuli/Matt/Echt/Scheier/Märk, *Chem. Phys. Lett.* **352**, 288 (2002)
  (Ne, Kr cluster ions); *Int. J. Mass Spectrom.* (avg-KER → binding);
  Bartl/Scheier/Echt, *J. Phys. Chem. A* 2013 (10.1021/jp406540p, anchored
  in Klots 1988). Caveat: measured systems are hundreds-of-meV
  rare-gas/molecular clusters, not ~9 meV He shells — the He-specific ε
  magnitude is not directly established.
- **NB-RQ23-4 (RQ3 mechanism — sequential, not prompt).** Ionic rare-gas
  clusters excited above their binding budget decay by **sequential
  single-atom evaporation** on the µs flight timescale, not prompt
  multifragmentation: Ar₁₀⁺–Ar₂₅⁺ metastable decay proceeds exclusively by
  one-monomer loss (the Ar₃⁺*→Ar₂⁺*→Ar⁺ series). Source: Märk group,
  *Int. J. Mass Spectrom.* (0168117686850121). Vote 3-0 (single-atom
  conclusion); the specific ~26,000 s⁻¹ rate example is actually a two-atom
  channel (claim 19, 2-1). Caveat: Ar is more strongly bound than He —
  transfer to sub-Kelvin I⁺Heₙ is by analogy.
- **NB-RQ23-5 (RQ3 self-termination — kills spec (a)).** The unimolecular
  rate is extremely energy-sensitive (a 10 % excitation increase raises the
  rate 10–20×), self-selecting survivors near a characteristic T ≈ D/G: an
  over-energized cascade **burns down its internal energy and self-terminates**
  rather than riding intact. Every observed X⁺Heₙ has undergone ≥ 1
  evaporation; internal energy is dissipated by evaporative cooling *after*
  the complex leaves the droplet. **Spec (a) inert/ride-intact is directly
  contradicted.** Sources: Hansen (Gspann); alkali snowballs
  (arXiv:0902.4713); Na⁺/K⁺-in-He (PMC3350777); Bartl/Scheier/Echt.
- **NB-RQ23-6 (RQ3 endpoint shape — favors (c) over (b)).** The terminal
  small-n distribution is **structured by shell-closure magic numbers**, not
  a flat drain to bare: abundance steps at first-shell closure
  (He₁₂Ar⁺/He₁₂Kr⁺; alkali magic n ≈ 12–16), and pure Heₙ⁺ shows a dominant
  tightly-bound **He₂⁺ dimer "last survivor"** with magic n = 10, 14 and a
  dip at n = 12. Favors spec (c) (resolved small-n tail) over spec (b) (pure
  drain to bare). Sources: jp406540p; arXiv:0902.4713; "Protonated and
  Cationic Helium Clusters" (RG 339634192); Märk (rate anomalies ↔ magic
  numbers). **MEDIUM confidence** — all noble-gas/alkali/pure-He, not a
  weakly-bound I⁺ core; the strong incremental-binding numbers were *refuted
  for transfer* (below), so flat-D₀ vs shell-structured **for I⁺ remains
  open** (→ RQ4).
- **NB-RQ23-7 (methodology validation).** A quantitative evaporative-cascade
  analysis (Rabinovitch, Hansen & Kresin, *J. Phys. Chem. A* 2011,
  arXiv:1102.3476) reproduces observed terminal cluster abundances from a
  deposited-energy distribution — validating internal-energy-driven
  sequential-evaporation cascades as the correct tool for predicting terminal
  size distributions. **MEDIUM confidence** — validated for warm classical Na
  metal *anion* clusters (~0.8 eV binding, many modes), not He-tagged
  sub-Kelvin ionic clusters; supplies neither the per-shed ε nor an
  over-bound endpoint, only the generic cascade methodology.

**Refuted claims (0-3, recorded so they are not re-imported):**

1. "Abundance ∝ relative dissociation energy, so terminal magic numbers map
   directly onto Dₙ anomalies" — refuted as a clean readout of the D₀(n)
   ladder (jp406540p).
2. "The KER coefficient c is material-independent at ≈ 1.5 (fit coefficients
   a₁ = −1.892, a₂ = 0.311, a₃ = −0.054), transferable across ionic rare-gas
   systems" — refuted (Hansen Gspann). *(The c ≈ 1.5–2 range in NB-RQ23-1
   survives; the over-specific material-independent transferability did not.)*
3. "Each single evaporation cools the cluster enough to cut the rate ≥ 10×,
   self-terminating per event" — refuted (PMC3350777). *(The aggregate
   self-selection to T ≈ D/G, NB-RQ23-5, survived 3-0; the strong per-event
   version did not.)*
4. "Resolved small snowballs Na⁺He₁₋₃ directly detected, Ak⁺He₂/Ak⁺ ratio
   0.5–1.0 %, cores extend to ≥ Ak⁺He₄₁" — refuted (arXiv:0902.4713). *(This
   was the most direct empirical support for a resolved small-n tail; its loss
   weakens the spec-(c) n = 1 case specifically.)*
5. "He binding is strongly shell-dependent — 250 meV inner / ~6 meV second
   layer / ~0.6 meV outer" — refuted for transfer (RG 339634192). Leaves the
   flat-D₀ assumption **neither confirmed nor refuted** for I⁺ → RQ4.

**Open questions (carried to RQ4 / the fragmentation design):**

1. The actual per-shed ε for He evaporation from a weakly-bound (~9 meV) I⁺
   shell at sub-Kelvin T — does the classical Klots–Hansen `ε ≈ c·D/G`
   (G ≈ 23.5, hot-classical provenance) transfer to a quantum/superfluid
   few-mode He shell, or does the small-C_v correction dominate (ε → 0 or
   negative)? **The single largest uncertainty.**
2. Is D₀(n) for the I⁺ core genuinely flat at ~9 meV, or shell-structured
   with magic-number stabilization? (→ RQ4; the strong-snowball numbers were
   refuted for transfer, so the ladder shape is undetermined.)
3. Can spec (c) with the *small* theory ε (~0.5–1 meV) reproduce the
   experimental 43.5 % bare / 17.5 % n = 1 split — or does matching n = 1
   need a larger-than-theoretical ε **or** an explicit deep first rung (RQ4)?
   *Note the sign:* smaller ε → deeper strip before self-termination →
   *thinner* small-n tail, so a too-small ε may under-produce n = 1.
4. Over µs flight, is whole-shell collective blow-off ever a competing channel
   vs strict single-atom sequential evaporation for a sub-Kelvin weakly-bound
   n ≈ 21 shell (untested in this regime)?

**Synthesis (for adjudication, NOT adjudicated):**

- **Spec (a) inert is out** (NB-RQ23-4/5: over-energized ionic He complexes
  demonstrably shed sequentially in flight, never ride intact).
- **Spec (c) is favored over spec (b)** (NB-RQ23-6: structured, magic-number
  small-n endpoints, not a flat drain to bare) — *architecture* recommendation:
  a self-terminating ε > 0 boil-off producing bare bulk **plus** a decreasing
  small-n tail.
- **ε form:** the Klots–Hansen prescription `ε ≈ c·D(n)/G` with c ≈ 1.5
  (ion-induced-dipole limit) and G ≈ 23.5 (~µs) → **ε ≈ D/16 ≈ 0.5–0.6 meV**,
  ~an order of magnitude below the model's feared 5–20 meV, with an optional
  finite-heat-bath reduction (−D/2C_v). One nearly parameter-free term governs
  both RQ2 (per-shed release) and RQ3 (endpoint).
- **What the literature cannot give:** the He-specific ε magnitude for an I⁺
  core, and the n = 1 population — which is set *jointly* by ε and the ladder
  bottom (RQ4) and cannot be separated by this run. **Recommended posture:**
  adopt spec (c) + the Klots ε form as the architecture; keep the ε magnitude
  a **Bounded** knob (theory prior ~D/16, swept because the quantum-He transfer
  is untested); promote **RQ4 (ladder-bottom depth) to the next research
  target** since ε and D₀(1) jointly own the small-n tail and cannot be
  calibrated independently.

**Cross-links:** `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4e (OQ-F/OQ-B
derivations; the "top two experimental bins" spec (b)/(c) n = 1 analysis this
register grounds); RQ4 (coupled — ladder bottom); RQ5 (in-flight ε-like
drains). Deliverable stance mirrors the RQ1 register: recorded with full
provenance, the user adjudicates, any accepted change gets its own design doc
(dimensional analysis, interchangeable enum arm) behind `[PROCEED TO
IMPLEMENTATION]`.

**Primary-source grounding (Albrechtsen/Stapelfeldt et al., *Nature* 623, 319
(2023) — `Stapelfeld_Paper_Ion_Solvation_in_Helium_Droplets.pdf`, the model's
own foundational paper; read 2026-07-10).** The deep-research fan-out reached
the target only by analogy (noble-gas/alkali/anion/pure-He); the model's
*foundational* reference is itself a He shell around a **cation** (Na⁺Heₙ,
Coulomb-ejected, dissociating in flight, VMI-detected) — the same physics class,
and the paper our biphasic mechanism ports (D₀(N), E_bind = Σ D₀, "shed one He
when E_int > D₀(N), repeat until E_int < D₀(N−K)", Poisson pickup, Newton-type
dissipation). It speaks directly to RQ2/RQ3/RQ4:

- **RQ2 (ε ≈ 0) — explicit in the source.** "the He atom will have very low
  kinetic energy because the dissociation is the result of a statistical
  redistribution of the internal energy … just above the dissociation
  threshold … **by neglecting the kinetic energy of the dissociation products,
  E_int will be lowered by D₀(N)**" (Methods). The model's drain-exactly-D₀
  convention *is* the primary source's convention, justified as near-threshold
  statistically-cold dissociation. Independent magnitude datum: atoms emitted
  during solvation carry **~14 cm⁻¹ ≈ 1.7 meV each** (p.9) — small, consistent
  with Klots ε ≈ D/G ≈ 0.6 meV (NB-RQ23-1), far below the feared 5–20 meV.
  **Theory (small nonzero) and the primary source (neglected) agree ε is
  small.**
- **RQ3 (sequential shed to self-termination; kills spec (a)) — MD-verified in
  the source.** "the process will repeat until E_int is lower than D₀(N−K) …
  further dissociation is energetically forbidden"; classical MD of hot
  Na⁺He₆/He₁₄ + Coulomb ejection: "dissociation takes place almost exclusively
  during the first tens of picoseconds … the vast majority … finished after
  400 ps", "Na⁺He₁₄ fragments significantly even for low total kinetic
  energies", and "we measure the complex sizes as they are after any
  post-ejection dissociation." The over-energized complex **sheds sequentially
  and self-terminates over tens–hundreds of ps (≪ 8.5 µs flight); it never
  rides intact** — spec (a) contradicted by the model's own source directly,
  supplying the He-specific fate + timescale the web literature (NB-RQ23-4/5)
  gave only by analogy.
- **RQ4 (bottom-heavy ladder) — the source's Na⁺ ladder is *not* flat.** D₀ is
  **high and ~constant for N = 1–6 (≳ 140 cm⁻¹ ≈ ≳ 17 meV)**, drops through the
  middle (D₀(9) = 140, D₀(10) = 72 cm⁻¹ ≈ 8.9 meV), and flattens **shallow for
  N ≥ 13 (25 cm⁻¹ ≈ 3.1 meV)** (Methods + Extended Data Table 1). A sequential
  ε ≈ 0 cascade **stops at small n because the bottom rungs are deep** —
  n = 1–6 are physical "last survivors." **The small-n / n = 1 population is a
  ladder-bottom (RQ4) effect, not an ε tail** — the structured endpoint the
  deep-research spec (c) reached via ε, the source reaches via the ladder with
  ε ≈ 0. Our flat 9.22 meV bottom is the outlier.
- **Caveat:** Na⁺ (alkali) vs I⁺ (halogen) is the residual transfer gap — the
  deep I⁺ bottom is the *leading, source-motivated* hypothesis, not established
  (→ RQ4).

#### RQ2/RQ3 adjudication — 2026-07-10 (user decision)

Supersedes the harness "spec (c) with ε" lean in the synthesis above: the
primary source (Stapelfeldt) reaches the same structured small-n endpoint with
**ε ≈ 0 + a deep ladder bottom**, the more parsimonious and source-consistent
reading.

1. **RQ2 — ε ≈ 0, small ε neglected (ADJUDICATED; convention retained, now
   sourced).** The per-shed kinetic-energy release is **small and neglected**:
   a shed drains exactly D₀(n), the evaporated He leaving translationally cold.
   Sourced two ways — Stapelfeldt explicitly neglects the dissociation-product
   KER (near-threshold statistical dissociation), and Klots–Hansen gives only
   ε ≈ c·D/G ≈ D/16 ≈ 0.5–0.6 meV (finite-heat-bath reduction toward 0 for the
   few-mode shell); the feared 5–20 meV is physically ruled out (µs survival).
   The optional Klots ε(n) = c·D(n)/G refinement is **recorded but not adopted**
   (negligible vs. the ladder-shape uncertainty). *(Thesis-ready statement: the
   kinetic-energy release accompanying each He evaporation is small — of order
   D₀/G ≈ D₀/16 ~ 1 meV, and translationally cold at threshold — and is
   neglected, following Albrechtsen/Stapelfeldt (Nature 2023); each shed
   therefore removes exactly the dissociation energy D₀(n) from the internal
   reservoir.)* Resolves findings **OQ-F**.
2. **RQ3 — spec (a) inert REJECTED; sequential shed-to-self-termination ADOPTED
   (ADJUDICATED).** The "suppressed = rides to the detector intact at n = 21"
   convention is **not physical** and is rejected: an over-energized complex
   sheds sequentially, draining D₀ per atom and stopping when E_int < D₀(n),
   over tens–hundreds of ps (Stapelfeldt MD) — deep inside the 8.5 µs flight.
   With ε ≈ 0 the *endpoint* is set entirely by the ladder: a flat bottom →
   drain to bare (the deep-research spec (b) limit); the physical deep bottom
   (RQ4) → self-termination at small n. RQ3's endpoint is therefore **deferred
   to RQ4** — the mechanism (sequential shed) is fixed, the stopping rung is
   ladder-determined. Reframes findings **OQ-B**: suppressed → fragments; the
   43.5 % bare / 17.5 % n = 1 split is a joint (Δ×-race weight into the
   over-bound class) × (ladder-bottom stopping) readout.
3. **RQ4 promoted to the load-bearing next research target.** The small-n /
   n = 1 population is a ladder-bottom effect (decision 2); the primary source
   supplies a concrete template (deep ~17 meV rungs N = 1–6, shallow ~3 meV
   top). The flat 9.22 meV Form-U bottom is the outlier to fix. RQ4 (I⁺–He
   incremental binding; He-tagging / ab-initio) is the next deep-research pass,
   since ε (now fixed ≈ 0) no longer competes with it for the small-n tail.

**Consequences / follow-through.** Any implemented change — removing the
spec-(a) suppression in favour of a sequential shed-to-self-termination
continuation, and reshaping the ladder bottom per RQ4 — is a model change that
gets its own design document (strict dimensional analysis, interchangeable enum
arm) behind `[PROCEED TO IMPLEMENTATION]`; nothing is coded in this phase.
Cross-doc edits still pending (offered, not yet made): `CALIBRATION_MAP.md`
OQ-F (→ resolved ε ≈ 0) / OQ-B (→ reframed), and a `drag_migration_log_tier2.md`
decision entry.

**Status:** RQ2 **closed** (ε ≈ 0, sourced). RQ3 spec **adjudicated** (a
rejected; sequential shed adopted; endpoint deferred to RQ4). **RQ4 is the open
successor** and the next research target. Residual caveat: the deep I⁺ bottom is
the leading hypothesis (Na⁺ analogy), pending RQ4.

### RQ4 — Ladder-bottom depth (= findings OQ-G)

> **Status update (2026-07-11, two-channel adjudication): PROMOTED to
> the critical path** (was "parallel non-blocking", post-Wave-10 entry).
> Under the solvated-branch re-targeting the small-n shape is RQ4's
> direct observable, and the flat Form-U bottom is **falsified at ~2×**
> on that branch (solvated n₁ = 31.0 % vs model ≤ 17 % at every
> (E₀, prior); flat rungs ⇒ flat bins, drag-law-robust). The
> experimental n₁/n₂ = 2.18 sits on the pre-registered ratio target
> below. Competing/combinable steepness source: position heterogeneity
> (Addendum H, Wave 12) — separated by the (n, KE) curve (ladder moves
> bins, not speeds).
>
> **Amendment (2026-07-11, Addendum H.2b):** the F.5 graded-ladder
> **floor variant** (D₀(1) = 13.3 meV X₂, transition-at-n=1 — knob-free,
> [IHe05]-sourced) is exercisable *before* this RQ's external
> calculation returns; the calculation is re-cast as the **arbiter of
> the transition location** (H.2b outcome (b) = a standing prediction
> for it). Anti-RQ4 fingerprint on the (n, KE) curve: the W12b
> `n_eject(depth)` feeder puts *fast* fragments into n = 1–3, the
> ladder feeder slow ones.

- **Model convention:** Form-U bottom is flat — D₀(1..5) = 9.22 meV
  (mixture, κ = 1).
- **Problem/alternative:** the real I⁺–He first rung is plausibly much
  deeper (ion-induced dipole, compressed first shell), making n = 1 a
  thermodynamic "last survivor" — a competing (and combinable) explanation
  for the elevated n = 1/n = 2 step (2.2×) that is budget-*robust*, unlike
  the RQ2/RQ3 race explanation.
- **Literature needed:** I⁺Heₙ binding energies / incremental binding
  curves (He-tagging spectroscopy, snowball structures, ab initio I⁺–He
  potentials + many-body shell calculations).
- **Discriminators:** the small-n abundance tail shape — the **only
  identified observable that reads the ladder bottom** (every in-window
  quantity probes the top rungs near n = 21); budget dependence separates
  it from RQ2/RQ3.

#### RQ4 sharpened target — the incremental-binding *ratios*, not just "is it deep" (2026-07-10 discussion)

Post-RQ2/RQ3 adjudication (ε ≈ 0), the small-n tail is a **pure ladder-bottom
readout**, which lets us state RQ4's target far more precisely than
"deep vs flat." Under ε ≈ 0 and a smooth reservoir p(E₀), the per-bin abundance
is `abundance(n) ≈ p(E₀,n) · window(n)` with `window(n) = D₀(n)·e^K`, so with a
slowly-varying p the **tail shape ∝ D₀(n)**. The experimental tail

| bin | bare | n=1 | n=2 | n=3 | n=4 | n=5 |
|---|---|---|---|---|---|---|
| weight | 0.435 | 0.175 | 0.080 | 0.052 | 0.040 | 0.033 |
| step to next | — | **2.2×** | 1.5× | 1.3× | 1.2× | — |

therefore demands, under a smooth reservoir, a **monotonically *decreasing*
incremental binding from n = 1**: `D₀(1):D₀(2):D₀(3) ≈ 2.2 : 1.5 : 1.3`-scaled
(a first rung ~2× the second, falling thereafter). **The target is these
ratios, not the absolute depth** (the absolute scale is absorbed by e^K / where
E\* sits).

- **The Na⁺ analogy actively *warns* here.** Stapelfeldt's Na⁺ ladder is
  "**high and essentially constant for N = 1–6**" — a *flat first-shell
  plateau*, which under a smooth reservoir gives a **flat** small-n tail, not
  the observed 2.2× step. So a "deep bottom like Na⁺" is **insufficient**;
  I⁺ would have to depart from the alkali template with a genuinely
  *decreasing-from-n=1* incremental binding. Physically plausible for a cation
  (first He closest / most polarized, each subsequent one weaker and more
  screened) but **not** what the primary-source analogue does — it must be
  checked, not assumed.
- **Why the ratios must come from RQ4, not the fit (the degeneracy).** The
  observed decreasing tail can be produced by a *decreasing p(E₀)* on a flat
  ladder (Wave-9's tuned fit; L2 = 0.016 but interior weights are
  basis-collinear / underdetermined), by a *decreasing D₀(n)* under a smooth
  reservoir, or by any mix — **all degenerate on a single-budget histogram.**
  Deepening the ladder by *fitting* D₀(n) to the same bins only relocates the
  overfitting from the reservoir to the ladder. The parsimony/physicality gain
  is real **only if D₀(n) is fixed independently** (ab-initio / He-tagging),
  after which the *required* p(E₀) is whatever remains and can be checked for
  smoothness.
- **Two decisive outcomes.** (i) I⁺ incremental binding **decreasing** at
  ≈ 2.2 : 1.5 : 1.3 → the "deep/structured ladder + smooth reservoir" reading
  holds and is more physical than Wave-9's tuned p(E₀); (ii) I⁺ bottom is a
  **flat plateau** (Na⁺-like) → the ladder gives a flat tail and the decreasing
  structure must come back from the reservoir (RQ1 provenance) or the
  droplet-R/K axis — i.e. the small-n shape is *not* a pure ladder effect.
  Either way RQ4 decides it, because it is an **independent input**.
- **Sharpened literature target:** the **absolute and incremental** I⁺Heₙ
  binding energies at the *bottom* — specifically **D₀(1), D₀(2), D₀(3)** and
  their ratios (He-tagging spectroscopy of I⁺Heₙ; ab-initio I⁺–He pair
  potential + many-body first-shell calculations; snowball/electrostriction
  structure) — plus whether the first shell fills flat (alkali-like) or
  tapered (decreasing) around a heavy halogen *cation*.
- **Independent confirmation:** the production budget test (RQ7; R0 = 2.666 Å /
  2.70 eV) breaks the same degeneracy from the other side — the ladder ratios
  are budget-robust, the reservoir/race weight is budget-dependent (E\* =
  Σ·e^K). RQ4 (structure) and the production run (budget dependence) are the
  two independent handles on the flat-ladder-+-structured-reservoir vs
  deep-ladder-+-smooth-reservoir degeneracy.

#### RQ4 findings — NB register (deep-research run, 2026-07-10)

> **Method/provenance of this register.** Multi-agent literature search
> (5 angles: primary ab-initio I⁺–He anchor; n-resolved He-tagging /
> messenger spectroscopy; many-body first-shell snowball structure;
> alkali-cation D₀(n) comparison template; PIMC / He-DFT ladder
> computation). 19 sources fetched, 82 claims extracted, top 25
> adversarially verified by independent 3-voter panels (18 confirmed,
> 7 refuted, 0 unverified), synthesis pass merged to 8 findings.
> **Status: recorded, NOT adjudicated.** The dominant caveat is
> **species transfer**: *no source in the surviving evidence supplies a
> measured or computed I⁺Heₙ many-body incremental ladder*
> (D₀(1):D₀(2):D₀(3)). The only genuine I⁺–He datum is the two-body
> Buchachenko [IHe05] pair potential, which fixes D₀(1) **only**. Every
> quantitative many-body ladder in the literature comes from *closed-shell
> alkali* cations (Na⁺/K⁺/Li⁺/Rb⁺); heavy noble-gas cations (Ar⁺/Kr⁺/Xe⁺)
> are size-relevant but electronically distinct.
>
> **Headline (the load-bearing result).** The RQ4-sharpened hope — that a
> genuinely *decreasing-from-n=1* incremental ladder (ratios ≈ 2.2:1.5:1.3)
> makes the small-n tail a pure, budget-robust ladder-bottom readout — is
> **not corroborated by any analog**. Every well-characterized cation fills
> its *first shell* as a **FLAT plateau**, with the ladder downslope located
> at **shell closure (N ≈ 9–18)**, not at the bottom rungs. So either
> open-shell I⁺ departs qualitatively from all studied cations (unknown,
> untested), or the observed 2.2× n=1/n=2 step is **not** a pure ladder
> effect and must be carried by the reservoir shape p(E₀) or by a non-zero
> per-shed ε (which would reopen the ε ≈ 0 adjudication).

**Confirmed NBs (vote 3-0 unless noted):**

- **NB-RQ4-1 (the I⁺ anchor is two-body only — fixes D₀(1), not the
  taper).** The sole genuine ab-initio I⁺–He datum, Buchachenko et al.
  [IHe05] (UCCSD(T), relativistic small-core pseudopotential), is a
  **diatomic pair potential** ("RG–I" = one rare-gas atom + one iodine
  species; built for ion mobilities). It contains no multi-He I⁺Heₙ cluster
  and therefore **cannot resolve flat-vs-tapered first-shell filling** — it
  anchors D₀(1) and nothing above it. (JCP 122, 194311 (2005).) Directly
  confirms the §RQ4 premise that the ladder *shape* is an independent input
  the pair potential does not supply.
- **NB-RQ4-2 (neutral He–halogen is the wrong charge state).** HeI / HeBr
  ab-initio wells (RSC b501253h) are **dispersion-bound doublet** (²Σ⁺/²Π)
  *neutral* van der Waals complexes; I⁺ (³P, p⁴) instead gives
  **charge-induced-dipole triplet** (³Π/³Σ⁻) states (Buchachenko: ³Π
  D_e = 143.9 cm⁻¹, ³Σ⁻ D_e = 63.6 cm⁻¹). Records which neutral-halogen
  sources are off the correct charge state and must **not** be imported as
  I⁺ anchors.
- **NB-RQ4-3 (the alkali template is a FLAT first-shell plateau — decisive
  counter-template).** For the closed-shell alkali cations the *first-shell*
  incremental (evaporation) energy is **high and essentially constant**, not
  decreasing from n=1: Na⁺ D₀(N) "high and essentially constant for N = 1–6"
  ([Nat23]); PIMC E_evap(N) "remain around 31 meV up to N = 6 where a marked
  decrease is observed" (arXiv:2502.11783); K⁺ "first 12 He attachments have
  a similar E_evap before decreasing", Na⁺ "gradual reduction reaching a
  plateau at N=10", Li⁺ constant to N ≈ 6 (arXiv:2510.12330); TDDFT "binding
  of the first five He atoms occurs at a constant rate" (Poissonian
  attachment ~2.0 He/ps). *(3-0 across the group; one member 2-1.)*
  **A flat plateau produces a flat small-n abundance tail — so on the alkali
  template the ladder alone cannot generate the observed 2.2× n=1/n=2 step.**
- **NB-RQ4-4 (the alkali downslope sits at SHELL CLOSURE, not the bottom
  rungs).** Na⁺ D₀(N) is flat through the first shell and then **drops
  sharply at closure** — 140 cm⁻¹ at N=9 → 72 cm⁻¹ at N=10, first shell
  completing near N=12 (icosahedral) ([Nat23] Extended Data Table 1). The
  *same abundance-ratio-tail methodology RQ4 uses* (I_N/I_{N−1}) locates its
  decisive dip **at the full first-shell size** (Rb⁺: sudden drop at N=15,
  closure N₁ = 16–18.5; arXiv:0902.4713), **not** as a small-n taper. This
  directly tensions the RQ4 premise that the small-n tail is a
  bottom-rung readout — in every alkali case the informative feature is the
  shell-closure size, not D₀(1):D₀(2):D₀(3).
- **NB-RQ4-5 (shape is set by many-body screening, not pair-well depth).**
  The Na⁺–He pair well (ILJ ε = 43.0 meV at r_m = 2.31 Å, ~347 cm⁻¹) is
  **~2.4× deeper** than the I⁺–He ³Π pair well (143.9 cm⁻¹ ≈ 17.8 meV; even
  more so vs the SO-averaged ~9.2 meV) — yet Na⁺ still fills flat. The
  plateau mechanism is **charge screening**: first-shell He are "equivalently
  arranged … binding energy determined by the two-body cation–He potential
  alone"; only atoms *beyond* the first shell are weakened (arXiv:0902.4713,
  2502.11783). *(Well-depth arithmetic 3-0; plateau-mechanism claim 2-1.)*
  Consequence: the pair-well *magnitude* is non-transferable, and depth does
  **not** predict taper — a deep first rung does not by itself make a
  decreasing ladder.
- **NB-RQ4-6 (heavy noble-gas size analog — weak magic, larger shell,
  smooth tail; MEDIUM, 2-1).** Ar⁺/Kr⁺/Xe⁺ close the first He shell near
  **n=12** (icosahedral) with abrupt D_n drops, but for the closest I⁺
  *size* analog **Xe⁺** the He₁₂ anomaly is weak (~10% vs ~50% for Ar⁺)
  because the He–Xe⁺ bond (~3.4 Å) far exceeds the He–He spacing (2.97 Å) —
  "helium atoms too small to complete a solvation shell" (PMC4166691,
  Renzler/Scheier). Implies a comparably large heavy halogen I⁺ (R_e = 3.25 Å)
  has a first shell **likely larger than 12 and only weakly magic → a smooth,
  low-structure abundance tail**. *Caveat:* noble-gas cations form a covalent
  charge-resonance [He–Ng–He]⁺ core (magic at n=2) unlike I⁺
  electrostriction, so their small-n D_n is not a clean I⁺ proxy; this source
  gives magic numbers, not the D₀(1):D₀(2):D₀(3) ratios RQ4 needs. **Note the
  tension with the model's n\* ≈ 21** (this analog suggests a smaller,
  weakly-structured first shell).
- **NB-RQ4-7 (a real bottom feature exists — but as a *yield* magic, not a
  theory incremental step).** Both Na⁺Heₙ and K⁺Heₙ show an abrupt
  ion-yield drop at **n=2** (Na⁺ by >2×), signalling enhanced M⁺He₂
  stability (An der Lan et al., *Chem. Eur. J.* 2012, PMC3350777). **But
  many-body theory does not reproduce it as a smooth incremental-binding
  step** — calculated K⁺ D_n "gradually decline except for the drop at
  K⁺He₂." So a genuine n=1/n=2 abundance step is *observed* even for
  closed-shell cations, yet it is a discrete magic feature, not a monotone
  D₀(n) taper of the RQ4 form.
- **NB-RQ4-8 (methodology identity — the arbitration read is the Stapelfeldt
  read; ε ≈ 0 is their lower bound).** The project's energy-gated
  evaporation mechanism *is* the Stapelfeldt-group mechanism: a departing
  Ak⁺Heₙ "quickly gets rid of any internal energy by shedding the number of
  He atoms energetically available … below the dissociation energy of the
  now reduced-size ion complex", with ⟨E_disp⟩ reconstructed from the
  terminal size distribution weighted by summed incremental E_bind(n)
  (arXiv:2510.12330, Eq. 7). Confirms the arbitration methodology and that
  the **ε ≈ 0 convention corresponds to the lower bound of their bracketed
  ⟨E_disp⟩** — mechanism identity holds; only the species (alkali vs I⁺)
  differs.

**Refuted claims (recorded so they are not re-imported):**

1. "Alkali Na⁺/K⁺ ladders **decrease monotonically from n=1** (no plateau)"
   (arXiv:physics/0702169) — **refuted 0-3** against the stronger PIMC/[Nat23]
   plateau evidence. *(Its loss is load-bearing: the one source that would
   have handed RQ4 a decreasing-from-n=1 alkali template did not survive
   verification — the flat-plateau reading is the robust one.)*
2. "physics/0702169 computes exactly the incremental D₀(n) ladder RQ4 needs"
   — **refuted 1-2** (its ladder claim is the one killed in (1)).
3. "K⁺Heₙ evaporation energies **decline gradually to n=10** = a tapered
   ladder supporting decreasing-from-n=1" (PMC3350777) — **refuted 0-3**
   (the flat-through-shell reading dominates; the only sharp bottom feature
   is the n=2 magic, NB-RQ4-7).
4. "Heavy noble-gas incremental binding is **strongly tapered** via a
   charge-resonance [He–Ng–He]⁺ trimer core (first two He far more bound)"
   (PMC4166691) — **refuted 0-3** as a transfer to I⁺ (the covalent-core
   physics is noble-gas-specific, not I⁺ electrostriction).
5. "K⁺Heₙ binding **peaks at n=12** (magic), i.e. non-monotone filling"
   (S1093326321000814) — **refuted 1-2**.
6. "Ar⁺ abundance tail is governed by geometric 12/20/12 shell packing
   (magic 12, 32, 44), so a heavy-cation tail encodes packing not D₀(n)"
   (jp406540p) — **refuted 1-2** (recorded as a live *caveat direction*, not
   an established transfer).
7. "Xe⁺ bond length washes out first-shell geometric closure → no sharp I⁺
   first-shell structure" (jp406540p) — **refuted 1-2** as stated; the
   surviving, more careful version is **NB-RQ4-6** (weak-magic/larger-shell,
   MEDIUM).

*(Net of the refutations: the two claims that would have *supported* a
decreasing-from-n=1 ladder — the "monotone alkali" readings (1)/(3) — both
died 0-3, while the flat-plateau NBs survived 3-0. The literature's verdict
is robustly **flat first shell, downslope at closure**, for every studied
cation.)*

**Unverified leads (extracted below the top-25 cut — NOT verified;
follow-up candidates):**

- **Cs₂⁺-in-He (heavy cation, near-flat):** He binding "within 2.4–3.1 meV
  per He for up to 12 helium atoms" — a *near-plateau* incremental filling
  for a heavy dimer cation (supports the flat-shell reading, but Cs₂⁺, not
  atomic I⁺).
- **H₂O⁺/H₃O⁺ He-tagging (molecular cation, tapered):** stepwise He addition
  in a 5 K trap shows incremental He binding **decreasing** with n across the
  first shell, driven by noncooperative three-body induction (PCCP 2022,
  d2cp01192a) — the one lead toward a *tapered* first shell, but for a
  molecular cation with directional H-bonding sites, not a monatomic heavy
  halogen. A possible mechanism by which an *anisotropic* charge (I⁺ ³P is
  open-shell/anisotropic) could taper — flagged, unverified.

**Synthesis (for adjudication, NOT adjudicated):**

- **The tapered-ladder hypothesis is uncorroborated.** No measured or
  computed analog supplies a 2.2:1.5:1.3 *decreasing-from-n=1* first-shell
  ladder. Every quantitative cation ladder is a **flat first-shell plateau
  with the downslope at shell closure** (NB-RQ4-3/4), and the two sources
  that claimed a monotone-decreasing alkali ladder were refuted 0-3.
- **So the RQ4 "two decisive outcomes" are resolved toward outcome (ii).**
  The sharpened target posed: (i) I⁺ decreasing ≈ 2.2:1.5:1.3 → deep-ladder
  + smooth-reservoir holds; (ii) I⁺ flat plateau (Na⁺-like) → the decreasing
  small-n structure must come back from the **reservoir p(E₀)** (RQ1
  provenance) or the **droplet-R/K axis**, i.e. the small-n shape is *not* a
  pure ladder effect. **The literature points to (ii)** for every studied
  cation — the analog first shells are flat, so the observed 2.2× step is
  most likely a reservoir/race feature (Wave-9's tuned p(E₀); findings §4f/§4g)
  rather than a ladder-bottom feature.
- **The one escape hatch is genuinely open, not closed.** Open-shell I⁺ (³P,
  anisotropic charge) *could* fill its first shell tapered where all studied
  *closed-shell* cations plateau (the H₂O⁺/H₃O⁺ lead shows anisotropic
  induction *can* taper) — but **no source demonstrates this either way for
  I⁺**. Establishing it requires a dedicated many-body / PIMC / ⁴He-DFT
  I⁺Heₙ calculation on the Buchachenko (or an updated) pair surface; the web
  literature cannot supply it.
- **First-rung magnitude stands, with an SO caveat.** The Buchachenko
  D₀(1) anchor (NB-RQ4-1) is confirmed as the *only* genuine I⁺–He value; no
  post-2005 higher-level recomputation of the He–I⁺ well depths surfaced. The
  X₂-only (106.9 cm⁻¹) vs SO-mixture (74.4 cm⁻¹) split remains the model's
  own `ladder_electronic_picture` fork, unresolved by this run.
- **Cross-cutting tension with n\*.** The heavy-cation size analog (Xe⁺,
  NB-RQ4-6) suggests a first shell that is **larger than 12 yet only weakly
  magic** — a smooth tail — which both undercuts a sharp ladder-bottom
  feature *and* sits against the model's n\* ≈ 21 assumption; worth a
  cross-check against the [I2-notes] I⁺@He₂₀₀₀ n\* = 21 anchor.

**Open questions (carried forward):**

1. **Does open-shell I⁺ (³P) fill its first He shell with a genuinely
   tapered incremental ladder from n=1, unlike all studied closed-shell
   cations?** No source answers — needs a dedicated many-body/PIMC or
   ⁴He-DFT I⁺Heₙ calculation. **The single load-bearing gap.**
2. If the analogs are flat and heavy shells close at n ≈ 12–18 (below the
   assumed n\* ≈ 21), can the observed 2.2× n=1/n=2 step be reproduced by
   the ladder alone, or must it be carried by p(E₀) or by non-zero ε
   (reopening the ε ≈ 0 adjudication)?
3. What is the correct first-shell size n\* for I⁺? The Xe⁺ analogy points
   to >12 and weakly magic — tensioned with n\* ≈ 21.
4. Does SO averaging of the I⁺ (³Π/³Σ⁻) pair curves materially change
   D₀(1) and the many-body screening onset (106.9 vs 74.4 cm⁻¹), and has any
   post-2005 higher-level ab-initio calculation updated the Buchachenko
   He–I⁺ well depths?

**Cross-links:** `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (Form U
ladder, `ladder_steepness` κ, `ladder_electronic_picture` fork, [IHe05]
first-rung provenance); `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4f/§4g (Wave-9
p(E₀) inversion — the reservoir-shape alternative this run points back
toward); RQ1 (reservoir provenance — the outcome-(ii) fallback);
RQ2/RQ3 (ε ≈ 0 adjudication that open question 2 would reopen); RQ7
(budget-dependence handle that independently breaks the
flat-ladder-+-structured-reservoir vs deep-ladder-+-smooth-reservoir
degeneracy). Deliverable stance mirrors the RQ1 and RQ2/RQ3 registers:
recorded with full provenance, the user adjudicates, any accepted change
gets its own design document (strict dimensional analysis, interchangeable
enum arm) behind `[PROCEED TO IMPLEMENTATION]`. Nothing coded in this phase.

### RQ5 — µs-flight channels (= findings OQ-E, design doc §4)

- **Model convention:** over the 8.53 µs continuation the only active
  channel is RRK evaporation — no radiative cooling, no electronic
  relaxation in flight, no residual-gas collisions.
- **Literature needed:** radiative (vibrational/recurrent-fluorescence)
  cooling rates for small ionic clusters on the µs scale; metastable
  electronic states of I⁺ and their in-flight relaxation; typical
  background-collision rates at the experiment's vacuum.
- **Discriminators:** validity of the detected read as-is; any additional
  E_int drain in flight acts like a weak, slow ε (couples RQ2's
  conclusions at the detector).

### RQ6 — Effective RRK dof (s_eff) — resolved in-model, literature cross-check outstanding

- **Status:** the classical `s = 3n−3` (60 modes at n = 21) was falsified
  by Waves 1–2 (freeze at n ≈ 20); constant s_eff promoted to Bounded
  (landing ungated ≈ 8–12, gated ≈ 30 at the staircase; detector-read
  nearly s-blind, §4e).
- **Literature needed:** effective RRK/RRKM dof counts for quantum
  clusters (mode freezing / weak coupling in cold He shells); He-cluster
  heat capacities at ~0.4–few K; evaporative-ensemble treatments with
  reduced dof.
- **Purpose:** cross-validate the Bounded band and the physical reading
  (I5: bath ~an order of magnitude smaller than classical), and inform
  the parked n-dependent alternative `s = α·(3n−3)` (findings OQ-C).

### RQ7 — Production Coulomb kinematics (= findings OQ-H, fired by Wave 8, 2026-07-09)

- **Model convention:** `coulomb_available_eV` has exactly one physics
  reader — the S2 onset deposit `E_int(0) = f_int·E_avail`. The
  Coulomb-explosion mechanics (initial separation, charge state, fragment
  speeds), ejection time, and cooling exposure are **budget-blind**:
  Wave 8 measured K₂₇₀ ≡ K₀₈₀ = 0.898297 exactly, and the suppression
  step is budget-invariant in absolute E_int(0) (E* = 0.4612 eV at both
  budgets, MD-verified).
- **Problem:** physically the 2.70 eV channel means faster fragments,
  earlier ejection, shorter exposure — none modeled. Every
  budget-dependent prediction is void in-model: B.1(2) in both versions
  (t×-scaling and the RQ1-inverted earlier-t_eject route), and the
  RQ3-vs-RQ4 bare-vs-budget discriminator. What survives: the whole gated
  map transfers across budgets verbatim in absolute E_int(0).
- **Resolution shape:** model the production channel's kinematics
  (initial-separation / charge-state arm behind an enum, strict
  dimensional analysis, `[PROCEED TO IMPLEMENTATION]`), or justify the
  9 Å kinematics as shared by both scenarios. Domain input: what CE
  channel produces the 2.70 eV KER (I⁺+I⁺ at ~2.7 Å? I⁺+I²⁺?) and at
  what fragment speeds.
- **Discriminators:** fragment-speed / VMI comparison at production
  conditions; any observed bare-fraction budget dependence (unusable as a
  model test until this RQ is resolved).

#### RQ7 status — PROMOTED to the active next target (post-Wave-10, 2026-07-11)

Wave 10 (findings §4h) made RQ7 **load-bearing and co-requisite**: the
narrow-E₀ reading of the experimental histogram works only at the production
K-scale (bracket K₀ ≈ 0.49) and is excluded ~19× at the pinned 9 Å
kinematics — so the production exposure **K₂.₇₀ is now the single unknown
that decides between "narrow E₀ + faster ejection" (Wave-10 landing) and
"E₀ mass above 0.46 eV" (super-solvation, back to RQ1). Resolution path
(Move 1 of the post-Wave-10 program order, log 2026-07-11):

- **Measure, don't model-from-scratch:** the production channel is fragments
  born at the ground-state separation R₀ = 2.666 Å (14.4 eV·Å/2.666 Å =
  5.4 eV shared = 2.70 eV/fragment; the 9 Å condition is the same formula at
  1.6 eV shared). A **deterministic Tier-0-style run at R₀ = 2.666 Å** with
  the locked drag bundle, read through the Wave-10-validated exposure
  integral, measures K₂.₇₀ directly and replaces the ballistic ×0.545
  bracket. Expectation (pre-registered): drag eats part of the extra speed →
  **K₂.₇₀ ∈ ≈ [0.5, 0.7]**, above ballistic, below probe.
- **Independent discriminator:** final fragment speeds vs
  `data/reference/vmi_iplus_he.csv` — pins K₀ with no evaporation physics.
- **Validity caveat to argue in the design doc:** the pure-cubic
  γ = g·b·v² was calibrated in the 9/18 Å TDDFT windows; ~1.8× speed may
  exit the calibrated band. Include the n\*(R₀) opening check already named
  in the log (2026-07-10 RQ2/RQ3 entry).
- **Design doc delivered (2026-07-11):** `TIER2_STAIRCASE_PROBE_PLAN.md`
  **Addendum G (Wave 11)** — run matrix, dimensional analysis, the two
  validity decompositions (beyond-band + overlap exposure shares, both
  biasing K high → measured K₂.₇₀ is an upper bound), predictions
  W11-P1–P5, outcome shapes. MD leg gated behind
  `[PROCEED TO IMPLEMENTATION]`.
- **MEASURED (2026-07-11, Wave 11 — findings §4i, I37–I40):**
  **K₂.₇₀ = 0.746** (upper edge; ballistic 0.49 = hard lower bracket —
  41.7 % of the exposure accrues beyond the calibrated drag speed band,
  and the VMI speed comparison undershoots ~2×, both in the over-drag
  direction). The branch question softens into a calibration statement:
  the bare-crossing scalar E₀ tracks E\*(K₀) = Σ(21)·e^(K₀) ≈
  0.30–0.41 eV over the whole K bracket — in-RQ1-band everywhere — with
  the measured anchor E\*(0.746) = 0.396 eV. Neither the ~0.5-branch nor
  the ~0.9-branch fires as pre-registered; the landing sits between,
  at E₀ ≈ 0.38–0.41 eV. **Pending user adjudication** (supersede the
  0.28 eV working value? fire the drag-recalibration item?).
- **Law-conditionality NB (2026-07-11, two-channel adjudication):**
  K₂.₇₀ = 0.746 is a *pure-cubic-law* measurement. The endorsed Wave-13
  saturated-drag arm (Addendum H), once calibrated on the n = 1-gated
  VMI peak, **re-measures K₂.₇₀** with the same exposure machinery; the
  drag-recalibration item has **fired** (adjudicated 2026-07-11) and is
  carried by Addendum H, not by a TDDFT request (TDDFT breaks at 2.666 Å
  kinematics — the program's premise).

### RQ8 — Bare-I⁺ channel provenance & branching (opened 2026-07-11, two-channel adjudication)

- **Experimental input:** mean fragment KE falls 3.706 eV (n = 0) →
  0.066 eV (n = 17); bare I⁺ is KE-distinct from every snowball → the
  experimental reading is a **different ionization channel** for bare I⁺
  vs I⁺Heₙ ejection. (The n = 1 gate `vmi_iplus_he.csv` peaks at
  10.1 Å/ps ↔ 0.69 eV at mass 131.)
- **NB (2026-07-11, H.2b design freeze; committed-reference update
  2026-07-18):** the full per-n **mean-KE table** is **committed** at
  `data/reference/ihe_ked/IHe_KED_reference.csv` (per-n `meanKE_eV` =
  first moment of the 3-D P(E) — the value an MD forward model compares
  its own ⟨E⟩ against; `COLUMNS.md` + `IHe_KED_reference.provenance.json`;
  full n0/n1 spectra + trusted 3-D curves n0–n4). **The export
  prerequisite is SATISFIED** — the earlier "provenance-pending / the
  export prerequisite stands" language is retired. Committed values
  supersede the provisional quotes: **bare (n = 0) meanKE = 3.706 eV**
  (mode 4.758 eV; the stale ~2.9 eV is retired), **n = 1 meanKE =
  1.302 eV** (the stale 0.974 eV retired), decreasing monotonically to
  0.066 eV at n = 17; `noiseLimited = 0` on all 18 fragments (every value
  a genuine measurement). Error model (the scoring tolerance): per-point
  √(statErr² + sysErr²) + the calib (4 %) and condition (6 %) fractional
  bands as two *correlated* whole-curve shifts. **Kinematic corroboration
  of the two-channel reading (firmer on committed data):** bare at
  3.706 eV sits ≈ 1.4× above the ≈ 2.71 eV per-fragment ballistic ceiling
  of the 2.70 eV channel — bare cannot come from the production channel
  under *any* drag law, independently favouring candidates (a)/(b).
- **Question:** what channel produces fast bare I⁺, and what is the
  branching vs the biphasic (solvated) channel? Candidates:
  (a) vertical ionization near the inner turning point (R ≈ 2.4–2.55 Å,
  singly charged → 2.8–3.0 eV/fragment, ≈ undissipated prompt escape);
  (b) I⁺+I²⁺ (10.8 eV shared → 5.4 eV born, ~45 % dissipated — but then
  why zero He retention, and where is its snowball tail?);
  (c) any surface-sited sub-population (couples to the Wave-12 position
  axis).
- **Discriminators:** I²⁺ at m/z ≈ 63.5 in the mass spectrum (+ its KE);
  the bare-KE distribution *shape* (narrow/prompt vs broad with a slow
  evaporative shoulder ~0.2–0.3 eV — the RQ3 discriminator); any
  droplet-size correlation of bare KE; gas-phase reference comparison.
- **Consequences:** the bare bin (43.5 %) is a **channel-branching
  quantity, not a model target** — the biphasic arbitration surface is
  the solvated branch (n ≥ 1 renormalized + the (n, mean-KE) curve,
  Addendum H §H.2). The Wave-9/10/11 "bare lands untuned" results are
  reinterpreted as coincidental under this reading. Model-side bare
  production is *reported* (an upper-bound consistency check), pending
  this RQ. **Leg B (findings §4p / I62–I64, 2026-07-18) is the newest
  instance:** the dressed suppressed/bare class brackets 43.5 %
  *coincidentally*, and its detected KE resolves to fragmentation
  bookkeeping (the co-moving-vs-momentum-conserving break-up read), not
  the experimental bare channel. **Standing scoring convention (decided
  post-leg-B):** every remaining leg (C, D) and the re-pilot **renormalize
  n = 0 out** and score the **solvated branch** (n ≥ 1 histogram + the
  per-n mean-KE curve) against the committed `IHe_KED_reference.csv` and
  its error model — see `TIER2_STAIRCASE_PROBE_PLAN.md` §I.11 T9 and the
  migration-log decision entry (2026-07-18).
- **Method:** deep-research pass (CE channels of I₂ ionized in/on He
  droplets; vertical vs sequential double ionization; KER distributions)
  + experimental cross-checks above + the (n, mean-KE) reference export
  (data-contract prerequisite, Addendum H §H.2).

### RQ9 — Physicality of the fitted cooling time τ ≈ 4.1 ps vs the GAH25 pin (opened 2026-07-15, Addendum I Step 1b; parked by user decision same day)

- **Question:** the Addendum-I joint (v_c, τ) closure (findings §4k
  Step 1b, I49) selects τ ≈ 4.1 ps — ×0.63 of the GAH25-sourced 6.55 ps
  probe pin. Is that defensible physics? Specifically: (a) what is the
  actual uncertainty / regime of validity of the GAH25 cooling time
  (thermal impurity vs a hot, fast-moving ion in a collapsing bubble);
  (b) does the effective Newton-cooling constant plausibly shorten for
  a translationally fast ion (denser local environment, shockfront
  contact, ripplon/roton emission channels); (c) are there independent
  literature constraints bracketing τ for cationic dopants in He
  droplets at these energies?
- **Status:** **parked (user decision 2026-07-15)** — kept open, not
  researched now. The τ re-classification event is pre-registered in
  plan §I.8: at the Step-2 build the pin becomes a fitted knob
  (Bounded → Free/Derived) in `CALIBRATION_MAP.md` unless this RQ
  re-anchors it.
- **Consequences:** if τ ≈ 4.1 ps is indefensible, the Step-1b joint
  closure loses its arm-(a) mechanism and OQ-I arms (b)/(c) (per-shed
  ε / fate-map form / in-band revisit) reopen; if GAH25 carries ±40 %
  room, the closure stands on sourced ground and the calibration-map
  reclassification is cosmetic.
- **Method (when fired):** focused literature pass on He-droplet
  cooling rates (GAH25 primary source + successors; impurity
  translational vs internal cooling; bubble-collapse timescales), plus
  the model-side sensitivity already measured (the joint basin in τ —
  Step-1c scan).
- **NB (Step 1c, same day):** the basin refinement centers the fitted
  clock at **τ ≈ 3.8–4.0 ps (×0.58–0.61 of the pin)**, with joint
  closures spanning τ ∈ [3.0, 4.8] and full-house cells τ ∈ [3.4, 4.2]
  — so the question is whether GAH25 admits a ≈ ×0.6 cooling time, with
  ≈ ±0.4 ps of model-side slack (findings §4k Step-1c, I51).

### RQ10 — Shed-frame / momentum convention of the evaporation channel (= findings OQ-J, fired 2026-07-17 by the I.11.2 item-1 re-read)

- **Question:** what momentum does an evaporated He carry away from a
  *translationally fast* I⁺Heₙ complex? The delivered Tier-2 evaporation
  channel uses the **cold-shed** operator (He left at rest in the lab
  frame; the complex keeps its full momentum, so KE × m/m′ per shed —
  +0.9 eV / ×1.61 cumulative on a full 21→1 post-ejection strip at
  production speeds). The physically expected convention for thermal
  evaporation is **co-moving** (He leaves with the complex velocity plus
  an isotropic thermal recoil ε ~ meV — the Tier-1a-adjudicated
  continuous-velocity path). [Nat23]'s "cold He (≈ 0 KE)" grounding is an
  *at-rest* Na⁺ result where the two conventions coincide.
- **Why it matters (measured, findings §4n / I58–I60):** first-order on
  the (n, mean-KE) observable — the co-moving counterfactual reproduces
  the 1D twin's per-bin KE to ≤ 2 %, so the entire twin−MD KE divergence
  (I57) is this convention; the experimental n₁ mean (1.302 eV) sits
  between the two conventions (cold 2.42–2.48, co-moving 0.93–0.96 eV,
  capped-tail configs); the §4k (v_c, τ) KE closure is co-moving-based
  and does not transfer to a cold-shed MD at small n; the RQ3 bare-bin
  mean KE is a fragmentation-convention discriminator (1.18 vs 3.3 eV).
- **Couples with:** RQ2/OQ-F (ε is the *energy* the He carries; RQ10 is
  the *momentum* — one coupled "what leaves with the He" resolution),
  RQ3 (fragment speeds), the leg-B/T5 KE pre-registration basis.
- **Candidate resolution:** interchangeable shed-convention enum
  (cold / continuous_velocity — both operators already exist in
  `physics/mass_jump.py`) behind its own trigger; user adjudication
  required before the T9 endgame.
- **Status: working convention ADJUDICATED (user, 2026-07-17) —
  co-moving for the twin-parity legs (T5 onward);** cold retained as
  the interchangeable bound arm (`cold` stays the byte-inert config
  default, legs stamp `co_moving` — the T7 precedent; **enum build
  DELIVERED 2026-07-17**, CALIBRATION_MAP row 26). Grounds and consequences: findings §7 OQ-J adjudication
  NB. **The RQ itself stays open** as the physical-resolution question:
  the true convention is co-moving + isotropic thermal recoil, i.e.
  the momentum side of RQ2's ε — how much recoil, correlated how, is
  the literature/domain-expert item; the two operators bracket it
  (ε → 0 vs maximal backward kick). Histograms are essentially
  convention-blind (mass events never enter the RRK rate/gate/cooling;
  small-n occupants shed post-exit) — the switch moves only the KE
  axis (n₁ 2.42–2.48 → ≈ 0.93–0.96 eV on the capped-tail configs).

## 4. Coupling map — what each answer changes

| RQ | Primary observable consequence | Couples with |
|---|---|---|
| RQ1 E_int(0) | whether *anything* sheds; t×, the race margin, production behavior at 2.70 eV | RQ3 (G₀ at ejection), picture knob, VMI KER |
| RQ2 ε per shed | eternal vs converged terminals; deep-strip reach; fragmentation endpoint | RQ3 (endpoint), RQ5 (any in-flight drain), n = 1 |
| RQ3 suppressed fate | the 43.5 % bare peak; bimodality; the n = 1 tail | RQ1 (G₀), RQ2 (endpoint), budget dependence |
| RQ4 ladder bottom | n = 1/n = 2 ratio; small-n tail shape | RQ2/RQ3 (competing n = 1 explanation) |
| RQ5 flight channels | detected-read validity | RQ2 (acts as slow ε) |
| RQ6 s_eff | clock/arrival-state only (detector nearly s-blind) | staircase prior; RQ2 (ε ∝ 1/s) |
| RQ7 production kinematics | any budget-dependent observable; validity of 2.70 eV predictions | RQ1 (E_avail provenance), RQ3 (bare-vs-budget discriminator), F5 |
| RQ8 bare-channel provenance | the interpretation of the 43.5 % bare bin (branching vs model); the solvated-branch renormalization | RQ3 (slow-shoulder discriminator), RQ1 (KER/E_avail per channel), Wave-12 position axis |
| RQ9 τ physicality (parked) | whether the Step-1b joint (v_c, τ) closure stands on sourced ground; the CALIBRATION_MAP class of τ at the Step-2 build | RQ4 (the closure is rq4graded-conditional), OQ-I arms (b)/(c) (reopen if τ = 4.1 indefensible), W13 Step 2 |
| RQ10 shed-frame convention | every small-n mean-KE read (×1.6 at n₁); the twin↔MD KE comparability; transfer of the (v_c, τ) KE closure to MD | RQ2 (one "what leaves with the He" resolution), RQ3 (fragment speeds, bare-bin KE discriminator), leg-B pre-registration basis |

The experimental distribution's information, as currently understood:
**bare fraction ↔ the Δ× race + RQ3; small-n tail (n = 1–3) ↔ RQ2/RQ3
fringe vs RQ4 depth; mid-shell weight ↔ the in-bubble leak (f_int/Δ×);
n ≳ 19 absence ↔ consistent with any opened gate.**

> **Superseded reading (2026-07-11, two-channel adjudication):** under
> RQ8 the first clause is retired — **bare fraction ↔ channel branching
> (RQ8), not the Δ× race**; the model's arbitration surface is the
> solvated branch: **n = 1 steepness ↔ RQ4 depth vs Wave-12 position
> geometry (separated by the (n, KE) curve); the (n, KE) scale ↔ the
> Wave-13 drag arm (v_c); mid-shell weight ↔ in-bubble leak,
> unchanged.**

## 5. Method and cross-links

Method for this phase: per-RQ literature search → record findings with
full provenance as NBs under the RQ (source, values, applicability
caveats) → user adjudicates each convention → any accepted model change
gets a design document with strict dimensional analysis and an
interchangeable enum arm, behind `[PROCEED TO IMPLEMENTATION]`. Nothing
in this phase touches code, checkpoints, or the delivered probe
artifacts.

### RQ11 — Deep-bin KE undershoot: why do deeply-solvated survivors arrive too cold? (= findings I87, fired 2026-07-21 by the S6 N = 500 read)

- **Question:** at the adjudicated standing production point
  (`finc1v725`: capped_cubic v_c 7.25 / τ3.2 / E₀0.27 / rq4graded,
  N = 500) the detected mean KE at n ≥ 12 sits **30–60 % below** the
  experimental IHe_KED reference (n = 16: sim 0.028 vs ref 0.069 eV),
  carrying ≈ 72 % of χ²_med, while n = 1–11 agrees under the committed
  error model. Which mechanism owns the cold tail: (i) E2
  exposure/over-cooling of the slow deeply-solvated class, (ii) the
  ladder tail (rungs above the tapered bottom ride Form U verbatim —
  a freeze-out-ordering effect on the KE axis), or (iii) a missing
  recoil/relaxation channel (per-shed ε recoil — RQ2-coupled) whose
  cumulative effect is largest on deep survivors?
- **Why it matters (measured, findings §4z / I85–I87):** the sharpest
  experimental constraint the program has produced — the sole reason no
  S6 finalist lands the frozen joint acceptance (χ²_med 68.6–202.5 vs
  ≤ 30) while both histogram axes land; invisible at N = 50
  (sim-SE-dominated sigmas, the I85 unmasking mechanism), i.e. a
  genuinely new N = 500 observable.
- **Couples with:** RQ2/RQ10 (what leaves with each shed He — recoil
  and momentum-convention effects accumulate over the deep cascade);
  RQ5 (µs-flight channels). The E2 Landau arm is **excluded as the
  owner**: §4y measured the scored KE bins bit-flat across
  v_L ∈ {0.30, 0.40, 0.58} (the deep bins were already cold under
  `zero_gamma`).
- **Candidate resolution:** onset-n anatomy + velocity-class
  decomposition of the deep survivors at the standing point (read-only);
  the RQ2 per-shed ε literature number; a ladder-tail counterfactual in
  the twin (the KE axis is chord-cached — cheap). Any accepted model
  change enters as an interchangeable enum behind
  `[PROCEED TO IMPLEMENTATION]`.
- **Status (2026-07-22, findings §4bb / I89–I90):** the read-only
  decomposition ran at the standing point. **Candidate (i) is refuted**
  — KE is frozen at handover on every solvated bin (E2 `E_dissip` gain
  0.0000; the Landau-gated E2 drag touches only the in-droplet retained
  class; detection-stage sheds are n ≤ 15 and KE-negligible), so the E2
  cap is causally disconnected from the KE curve. The deficit accrues
  **inside the 30 ps MD window**, in the in-band pure-cubic segment
  ≈ 2–14 ps at v ≈ 2–7 Å/ps; fate is birth-dressing/droplet-size
  ordered and the post-window n-mapping is near-diagonal. Remaining:
  (ii) reframed as the in-window (KE, n) exit correlation, (iii) sized
  at ~4–8 meV recoil per shed over the ≈ 5-shed deep cascade, and a new
  (iv) beam-frame-dependent retained-population channel (RQ5-coupled;
  the retained class sits on the Landau floor at 0.0025 eV in the MD
  frame).

**Cross-links:** `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4c–§4e (the
derivations and numbers behind every RQ; insight register I13–I25);
`TIER2_DETECTION_STAGE_DESIGN.md` (§4 scope caveats → RQ5);
`CALIBRATION_MAP.md` (row 14 → RQ1; row 10 → RQ6; OQ2 register);
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (S1/S2/K1/K2
definitions); `drag_migration_log_tier2.md` (decision history);
`TIER2_STAIRCASE_PROBE_PLAN.md` (probe program contracts).
