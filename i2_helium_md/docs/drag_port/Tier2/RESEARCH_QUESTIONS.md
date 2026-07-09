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
cross-check of an already-adjudicated convention.

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

### RQ3 — Fate of a net self-unbound complex (= findings OQ-B, fragmentation hypothesis)

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
  fraction vs budget (0.80 vs 2.70 eV).

### RQ4 — Ladder-bottom depth (= findings OQ-G)

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

## 4. Coupling map — what each answer changes

| RQ | Primary observable consequence | Couples with |
|---|---|---|
| RQ1 E_int(0) | whether *anything* sheds; t×, the race margin, production behavior at 2.70 eV | RQ3 (G₀ at ejection), picture knob, VMI KER |
| RQ2 ε per shed | eternal vs converged terminals; deep-strip reach; fragmentation endpoint | RQ3 (endpoint), RQ5 (any in-flight drain), n = 1 |
| RQ3 suppressed fate | the 43.5 % bare peak; bimodality; the n = 1 tail | RQ1 (G₀), RQ2 (endpoint), budget dependence |
| RQ4 ladder bottom | n = 1/n = 2 ratio; small-n tail shape | RQ2/RQ3 (competing n = 1 explanation) |
| RQ5 flight channels | detected-read validity | RQ2 (acts as slow ε) |
| RQ6 s_eff | clock/arrival-state only (detector nearly s-blind) | staircase prior; RQ2 (ε ∝ 1/s) |

The experimental distribution's information, as currently understood:
**bare fraction ↔ the Δ× race + RQ3; small-n tail (n = 1–3) ↔ RQ2/RQ3
fringe vs RQ4 depth; mid-shell weight ↔ the in-bubble leak (f_int/Δ×);
n ≳ 19 absence ↔ consistent with any opened gate.**

## 5. Method and cross-links

Method for this phase: per-RQ literature search → record findings with
full provenance as NBs under the RQ (source, values, applicability
caveats) → user adjudicates each convention → any accepted model change
gets a design document with strict dimensional analysis and an
interchangeable enum arm, behind `[PROCEED TO IMPLEMENTATION]`. Nothing
in this phase touches code, checkpoints, or the delivered probe
artifacts.

**Cross-links:** `TIER2_STAIRCASE_PROBE_FINDINGS.md` §4c–§4e (the
derivations and numbers behind every RQ; insight register I13–I25);
`TIER2_DETECTION_STAGE_DESIGN.md` (§4 scope caveats → RQ5);
`CALIBRATION_MAP.md` (row 14 → RQ1; row 10 → RQ6; OQ2 register);
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (S1/S2/K1/K2
definitions); `drag_migration_log_tier2.md` (decision history);
`TIER2_STAIRCASE_PROBE_PLAN.md` (probe program contracts).
