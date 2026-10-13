# Thesis discussion points (final discussion round, 2026-10-02)

Material for the discussion chapter of the thesis. These are physical
arguments and open questions, **not** goals: the Q = 3 / channel-mixture
line is out of scope for the thesis (user decision, 2026-10-02) and appears
there at most as discussion. Measured results and projections are labelled
as such. No runs were made for this document.

Sources read: D0 `TIER2_PARAMETER_INFLUENCE.md` §14.3 (2026-10-02 note),
§15.8, §20; `drag_migration_log_tier2.md` 2026-10-02 entry;
`data/runs/h2b_forward_model/detfix_postpickup_table.csv` and
`detfix_postpickup_chord.csv`; `TIER2_DetectorStageFix.md` §5.7;
`CALIBRATION_MAP.md` row 30; `RESEARCH_QUESTIONS.md` (NBs of 2026-07-29);
`data/reference/README.md`, `ihe_ked/README.md`, `paper_cov/README.md`.
Later additions: §0 storyline, §8 mass influence and §9 Coulomb budget (2026-10-13, from
`THESIS_FINAL_GOALS.md` G5 / G11); §1 is to be refreshed by G6 (chord map
at B2) and G7 (KE₁ geometry-invariance caveat).

---

## 0. Storyline — "TDDFT constrains, the experiment selects" (G11.1; added 2026-10-13)

Spine of the thesis argument; numbers from `THESIS_FINAL_GOALS.md` G5
Outcome and G11.4 status. Tier 1 does not appear (retired, G11.3).

1. **Tier 0 — TDDFT constrains a ridge, not a law.** The 18 Å TDDFT
   trace, at the true shell mass m = I⁺ + 21 He (the shell stays at 21
   throughout), pins the drag force near v ≈ 3 Å/ps and the well / drag
   split at low energy (E_bind ≈ 0.117 eV), not the high-speed shape. The
   output is a one-parameter family F = a·v + c·v²: all 12 ridge points
   c 0–7.2 amu/Å pass 18 Å (RMSE 0.113–0.119 Å/ps, band 0.19). The 9 Å
   case (non-radial, transverse motion) is a one-sided exit bound only:
   it cuts c > 5 (c 7.2 fails). The cross-case cubic preference (G4.11)
   is contested by the radial 9 Å runs → caveat, not a foundation.
2. **Tier 2 — the experiment selects two knobs.** Only c on the ridge
   and E₀; τ = 6.4 ps ≈ the sourced 6.55 ps; physical droplet ensemble
   (`post_pickup`, §2). Production law **B2**: c 1.6, a 18.01 amu/ps,
   E_bind 0.1173 eV, E₀ 0.3297 eV (pooled N = 6000): n̄ 3.90 (gate
   3.83–4.31), midHot 1.002, W₁ 0.967, n₁ 0.211 (exp 0.310), KE₁ 0.609 eV
   (exp 1.302). B2 stays Tier-0-consistent at m = 21 (RMSE 0.118, 9 Å
   margin +0.82 Å/ps).
3. **Identifiability — the non-uniqueness is a result.** (a) The c the
   experiment selects depends on the droplet description (raw → 2.0,
   `post_pickup` → 1.6). (b) A differently shaped law, h405 (capped
   cubic, v_c 5.5), lands comparably (W₁ 0.755, KE₁ 0.640, midHot 0.948)
   but fails the 9 Å bound (−1.21 Å/ps). Neither TDDFT nor the detected
   I⁺(He)ₙ data identify a unique drag function in He nanodroplets.
   Caveat: h405 is not a controlled A/B partner (gate 14.2 Å, τ 4.4 ps,
   bundle fitted at the old gate).
4. **Effective variant B3 — an explicit departure.** A shallow well
   (E_bind 0.0482 eV, a 19.8, E₀ 0.3389) buys the high tail (n ≥ 15
   0.0075 vs B2 0.0008, deep-bin KE 0.77 of experiment, W₁ 0.795) at the
   price of leaving the TDDFT trace (18 Å RMSE 0.58 > 0.19). Presented as
   effective, not TDDFT-based.
5. **Geometry only reweights; the map sets KE** (§1; G6 figure, G7
   caveat). KE₁ is geometry-invariant within the deterministic model.
6. **Mass** (§8): the extraction mass is no longer an assumption; the
   transfer to small n is.
7. **Outlook — outside the deterministic model.** n₁ share and KE₁ (≈ half
   the experiment for every law; §3, §5, G8), the n ≥ 15 deficit (§7),
   the KE width at fixed n (Tier-3 collision noise), a sourced I⁺–He D₀
   (M1), a drag law for small n (§8), the calibrated Coulomb budget with
   the Q = 2 / Q = 3 mixture (§9; production runs E_C, ≈ 25 % above the
   source's 0.8·E_C — the generous bound for KE₁).

Fewer tuned knobs than in the earlier sweeps (v_c, E₀, τ); the
non-uniqueness is reported as a finding.

## 1. The chord → (n, KE) map (MEASURED; refreshed at B2 by G6, 2026-10-13)

Figure: `scripts/post_processing/tier2_g6_chord_map.py` →
`data/figures/thesis_tier2/g6_chord_map.{png,pdf}` (not committed): (a) scored n
vs chord (suppressed → n = 0, scorer convention), (b) KE vs chord, B2 per ion (N = 6000 pool, 12 000 fragments)
with the old-bundle map (h405 `post_pickup`, binned) overlaid, (c) chord
distribution per sampler. Numbers: `g5_postpickup_chord.csv` (binned,
`lqppb2`; the script anchors on it to 4 dp), `detfix_postpickup_chord.csv`
(h405), `g6_chord_stats.csv` (n = 1 window, KE₁ spread).

- The detected (n, KE) is a near-deterministic function of the
  **emission-axis chord**: the straight path from the birth point along the
  explosion axis (partner → fragment) to the droplet surface. The map does
  not depend on the size sampler (h405: `raw` vs `post_pickup` bin by bin
  within ≈ 0.01 eV and 0.1 in n; A5 likewise, G5 step A).
- Structure of the B2 map:
  - chord < 34 Å: the **suppressed** class (11 % of scored ions) — the
    complex is ejected net-self-unbound (E_int > Σ(n)) with its full shell.
    The code carries it at n_detected = 21; the scorer counts it as
    **n = 0** (bare candidate, RQ3; excluded from the solvated read, RQ8),
    and the figure follows the scorer. Physically such a complex
    fragments; how (and so its final KE) is open, so panel (b) shows the
    intact-complex KE;
  - 33–35 Å: the edge between suppressed and n = 1 (mixed bin: n = 1
    share 0.42, mean KE 0.90 eV);
  - **n = 1 window 34.5–39.9 Å** (q05–q95 of the n = 1 ions, median
    37.2 Å, width 5.3 Å); KE₁ 0.609 eV;
  - beyond ≈ 40 Å: n rises again (2 → 13) — **re-pickup after
    stripping**, not partial evaporation;
  - > 80 Å: mostly retained (0.89).
- **Old vs new bundle** (h405 → B2, panels a–b): the suppressed / n = 1 edge and the
  n = 1 window sit at the same chord (h405 window ≈ 35–39 Å) with nearly
  the same KE₁ (0.67 vs 0.63 eV in the 35–37 Å bin); the bundle change acts
  on the long-chord branch: under B2 KE falls more slowly with chord
  (60 Å: 0.26 vs 0.10 eV) and re-pickup is slower (mean n at 60 Å: 6.9 vs
  10.3).
- **Correction to a tempting reading:** the ion does not "choose the
  closest path to exit". The I–I axis at ionization fixes the direction and
  the partner's Coulomb push drives the fragment along it; the relevant
  length is the chord along that axis, not the radial distance to the
  surface. This is why radial birth depth is the weaker proxy.
- **Consequences.**
  - Geometry (droplet size, sampler, birth depth) only reweights the chord
    distribution (panel c: B2 and A5 `post_pickup` have bit-identical
    chords — the law does not enter; `post_pickup` shifts weight from
    < 40 Å to long chords); drag, τ, E₀ and the ladder set the map. The
    histogram is the chord distribution pushed through one map.
  - **No geometry change can move KE₁.** The n = 1 window position on the
    chord axis fixes KE₁; only the map can move it: Coulomb budget, drag,
    or flight mass (`TIER2_MASS_SCENARIOS.md`, M1 gate).
  - The n = 1 KE width in the model (SD 0.038 eV, the "needle") is the
    window width × the KE slope; the experimental width needs scatter in
    the map itself (Tier-3 noise or a channel mixture), not geometry. Scope
    of the claim: G7 caveat below.
  - Projection (not measured): on the B2 map KE ≈ 1.0 eV sits at a chord of
    ≈ 32–34 Å, at the suppressed / n = 1 edge where most ions are still
    suppressed.
    The ≈ 0.5–0.7 eV deficit (exp KE₁ median 1.13 / mean 1.30 eV) is
    equivalent to stripping completing ≈ 3–5 Å of chord earlier, or to KE
    falling less steeply along the chord.
- Caveats: the chord ignores deflection (partner Coulomb, drag) and the soft
  erf edge; h405 is binned only (its run directories are no longer on
  disk); the arms are not CRN, so the old-vs-new agreement is ensemble-level.

### 1.1 G7 caveat: scope of "KE₁ is geometry-invariant" (2026-10-13; zero MD)

Numbers: `g6_chord_stats.csv`, `g7_ke_width.csv` (both from
`tier2_g6_chord_map.py`, B2 pool).

**Re-check on the B2 map (MEASURED).**
- The n = 1 window is narrow: 5.3 Å of chord (q05–q95, 34.5–39.9 Å).
- The map is near-deterministic: KE₁ SD 0.038 eV, of which 0.036 eV is
  the KE slope along the chord (−0.021 eV/Å) times the chord spread inside
  the window; the residual at fixed chord is 0.012 eV. n at fixed chord is
  equally sharp (panel a: discrete bands, adjacent n overlap over ≈ 1–2 Å).
- KE width at fixed n, model vs experiment (IQR / median):

  | n | model (B2) | experiment | ratio |
  |---|---|---|---|
  | 1 | 0.101 | 0.804 | 8.0 |
  | 2 | 0.086 | 0.702 | 8.1 |
  | 3 | 0.092 | 0.706 | 7.7 |
  | 4 | 0.090 | 0.691 | 7.7 |

  The model is ≈ 8× too narrow in relative width at every n (consistent
  with §7.4, measured there at the earlier law). In absolute terms the gap
  at n = 1 is larger (σ 0.70 vs SD 0.038 eV, ≈ 18×) because the model
  median is about half the experimental one. Experimental quartiles are
  read from the unsmoothed 3-D P(E) (`IHe_KED_curves_n*.csv`); they match
  the reference median / mean columns to ≈ 4–5 % only (cuts not re-done),
  so only the relative width is taken from them.

**Statement for the thesis.**
1. "KE₁ is geometry-invariant" holds **within the deterministic model**:
   geometry moves ions along the chord axis, the n = 1 window sits at a
   fixed chord, and the map through it is near-deterministic, so KE₁
   cannot change.
2. The experimental KE width at fixed n shows the real map **scatters**
   (≈ 8× wider). With scatter the n = 1 window gets blurred edges and a
   chord-window selection can bias KE₁. Size (PROJECTION, not measured):
   at most the KE slope × a shift of the window's mean chord by up to half
   its width, −0.021 eV/Å × 2.7 Å ≈ 0.06 eV — small against the
   0.5–0.7 eV KE₁ deficit, so the deficit is not a geometry artefact of the
   deterministic map.
3. The near-determinism is partly **by construction**: the model has no
   fluctuation channel (zero noise; Langevin channel deferred to Tier 3),
   so "identical map in both sampler arms" is expected (same physics
   parameters). Map scatter (collision shot noise, §7.4) is the outlook,
   not modelled here.

## 2. Size sampler: `post_pickup` is production (MEASURED, D0 §15.8–15.9; updated 2026-10-13)

- **Production uses `post_pickup`** (user decision 2026-10-09, G5): the
  pickup-conditioned droplet ensemble, as described in thesis Ch. 3 and
  as used by legacy MATLAB production. The thesis production law B2 was
  arbitrated under it (lq ridge c 1.6, E₀ 0.3297 eV; pooled N = 6000).
  `raw` (the plain source log-normal at 40 bar / 14 K) was the sampler of
  the earlier standing points (h405, A5) and is kept only as the measured
  sensitivity. Ch. 3 and production now agree.
- **Sampler effect** (`post_pickup` − `raw`, same law, nothing re-tuned):
  - h405 (pooled N = 5000 each): ⟨N⟩ +29 %, R̄ +4.5 Å, n̄ +0.72, n₁ −0.033,
    W₁ +0.24, trap +0.057, supp −0.059, KE₁ +0.007 eV;
  - lq A5 (G5 step A): n̄ +0.70, n₁ −0.04; W₁, midHot, KE₁, S flat.
  The larger droplets only move ions from short chords (the n = 1 window)
  to long ones; the chord → (n, KE) map is sampler-invariant (§1).
- **Re-arbitration under `post_pickup`** (G5): E₀ or τ alone restore n̄
  only at midHot ≈ 0.85; lowering c on the TDDFT ridge (2.0 → 1.6) with
  E₀ re-placed holds every observable of A5 and improves W₁ (−0.13) and
  S (−0.25). The selected c therefore depends on the droplet description
  (raw → 2.0, `post_pickup` → 1.6) — part of the identifiability result
  (§0).
- The thesis headline claims (KE₁ deficit, chord-map argument) are
  sampler-invariant.

## 3. Q = 2 / Q = 3 channel mixture as an explanation of KE₁ (DISCUSSION ONLY)

### 3.1 The hypothesis

Both Q = 2 (I⁺ + I⁺, ≈ 2.16 eV per fragment at 0.8·E_C) and Q = 3
(I⁺ + I²⁺, ≈ 4.32 eV) feed n = 0, 1 and possibly 2; the high-n bins are
Q = 2. The model simulates Q = 2 only, which would explain why KE₁ and the
histogram cannot be landed together. Q = 3 would supply the high-KE part
of n = 1 (and 2) and pull the mean up, without retuning the Q = 2 point.

### 3.2 What supports it

- The (C) probe (D0 §20, 2026-07-29, E2-era pipeline, h405 pins,
  N = 1000): the channel mixture alone placed **KE₁ at 0.920** and broke
  the needle (KE₁ SD 0.038 → 0.286). This is the width mechanism the chord
  map requires.
- Mean arithmetic (projection): with a Q = 3 share w ≈ 0.2 of n = 1 and
  Q = 2 n = 1 at 0.64 eV, a 1.00 eV mean needs Q = 3 n = 1 ions at
  ≈ 2.4 eV, inside the 1.5–2.7 eV dressed-Q3 exit estimate (RQ NB).
- Reframing of "Q = 3 sits deeper": both channels share one chord
  distribution (same molecule, same droplet). Q = 3 has a **different
  map**: with a 4.32 eV start and ≈ 0.043 eV/Å capped force, its n = 0–2
  window should lie at much longer chords (the 40–65 Å bulk where Q = 2
  gives n = 2–10). Not measured.

### 3.3 What argues against it, or limits it

- **Q = 2 owns most of n = 1.** Event covariance (thesis Fig. 6.2): Q = 3
  is ≈ 7 % of I⁺He at 1.47×10¹⁴ W/cm² and ≈ 20 % at 2.94×10¹⁴ W/cm². So
  ≥ 80 % of n = 1 is Q = 2, and censoring n = 1 from a Q = 2 scoring
  discards a real constraint.
- **Shape tension.** The experimental n = 1 KED has 48.7 % of its weight
  above 1.15 eV, against a ≤ 20 % Q = 3 share. A Q = 2 needle at 0.64 eV
  plus a ≤ 20 % Q = 3 bump cannot produce that; roughly a third of the
  Q = 2 n = 1 ions would also need to sit above 1.15 eV. The deficit
  therefore remains a Q = 2 statement, which keeps the flight-mass arm
  (M1) as the open explanation. The two explanations are not exclusive.
- **Saturation (user's reading).** If Q = 3 is already saturated at
  300 mW, then ≈ 20 % is the ceiling of its n = 1 share at the KE-reference
  condition, which strengthens the point above.
- **KED power cross-check** (`ihe_ked/README.md`): between the 300 and
  600 mW series, n = 1 mean KE is invariant (+0.02 %), n ≥ 2 shifts
  +5–7 %, n = 0 shifts −6.3 %. Read loosely: Q = 3 mainly feeds n = 0 and
  contributes little to the n = 1 mean. Caveats: w at 600 mW is unmeasured,
  saturation would also explain the invariance, and n = 0 is entangled with
  count-rate effects.
- **The §20 histogram failure** (trap 0.465, n̄ 6.69, n₁ 0.080) came
  mainly from the single-ionization channel (E ≈ 0.53 eV cannot exit), not
  from Q = 3. A Q = 2 + Q = 3 arm with w_single = 0 was never run; it was
  argued only arithmetically, before the chord map and the detector-stage
  fix.
- The Q = 3 partner is not propagated as I²⁺ (CE design doc l. 331); any
  "no retune" Q = 3 prediction rests on an emulated partner.

### 3.4 Extracting the Q = 3 share per n from two powers

- Idea: unmix H(P) = (1 − w_P)·H_Q2 + w_P·H_Q3 from histograms at two
  powers, giving H_Q3 = [(1 − a)·H₂ − (1 − b)·H₁] / (b − a).
- Limits:
  - ill-conditioned: with a = 0.07, b = 0.20 the per-bin errors grow
    ≈ 6–8×, and inter-power systematics are already 5–7 % (KED
    cross-check);
  - the 7 % / 20 % values are shares of the I⁺He row, not global channel
    weights, so the unmixing is circular without another weight source;
  - channel shapes are not power-independent (He⁺ₘ covariance grows
    ≈ 6 % → 28 % with intensity; possible Q = 4; bare-I⁺ count-rate
    effects).
- The direct measurement is the covariance of each I⁺Heₙ row (n = 0, 1,
  2, …) with I²⁺ / I²⁺Heₙ at both powers. Two-power differencing is at most
  a consistency check.

## 4. Reference conditions (OPEN, data provenance)

- **Power ↔ intensity mapping corrected by the user (2026-10-02):**
  **150 mW ↔ 1.47×10¹⁴, 300 mW ↔ 2.94×10¹⁴ W/cm²** (150 mW = the
  "160 mW" series of `paper_v2`/`paper_v4`). This supersedes the
  300 / 600 mW mapping adopted on 2026-07-29 (RQ8 NB follow-up, later
  recorded as user-confirmed, CE design OQ-K), which is off by a factor
  of 2. Corrected in place in `RESEARCH_QUESTIONS.md` and the CE design
  doc; dated records (log, atlas plan, atlas findings) carry
  `[Corrected 2026-10-02]` markers. Consequence: the ~20 % Q3 share pairs
  with 300 mW; the 600 mW `ihe_ked` reference has no measured share.
- The gas-phase VMI band read (atlas findings, §3.5l read (iii)) still
  rises 300 → 600 mW (Q3/(Q2+Q3) 16.8 → 23.1 %). That bears on the
  saturation reading of §3.3: in the gas phase at least, the fast channel
  does not look saturated at 300 mW. Band reads carry Q2-tail leakage, so
  this is suggestive only.
- Conditions as the repo records them versus the user's statement:

  | Reference | Repo | User |
  |---|---|---|
  | `integrated_i_he_abundance.csv` (histogram arbiter) | no provenance ("to be completed") | 150 mW |
  | `ihe_ked/IHe_KED_reference.csv` (KE per n) | 600 mW, 17.10.24; 300 mW same-day cross-check | 300 mW |
  | `paper_cov` (covariance) | 300 mW, 03.12.24 | — |

  The KE entry is unresolved. Either way the histogram and the KE reference
  come from different conditions (and probably different days); scoring
  one model point against both assumes a power-independent ensemble.

## 5. Should n = 1 carry weight in the arbitration? (DISCUSSION)

- **Against using n = 1 as a target:** mixed channels (≥ 7–20 % Q = 3);
  a ≈ 4 Å chord window that is hypersensitive to the map; the reference's
  most background-sensitive bin (+29 % bg-off shift at 600 mW); the
  condition mismatch of §4. Currently n₁ is a hard gate clause
  (`n1_solv ∈ [0.19, 0.30]`).
- **Against dropping it:** n = 1 is the cleanest probe of where the map
  places the stripping threshold, and ≥ 80 % of it is Q = 2. It is most
  useful as a held-out prediction (share and KE) of a point fitted on the
  bins Q = 2 dominates.
- Arithmetic (projection; normalization equivalence of `_solv` and the
  covariance row share not verified): a Q = 2-only n₁ target is
  0.243 × (0.80–0.93) ≈ 0.19–0.23. The gate floor (0.19) and h405 (0.21)
  already sit in that range, so a Q = 3 correction of the share changes
  little. The KE side is where the residual lives.
- Framing for the thesis: KE₁ is reported as a finding ("the Q = 2 n = 1 KE
  is ≈ 0.36 eV low; a ≤ 20 % Q = 3 admixture cannot close the shape"), not
  as a fit residual to be weighted away.

## 6. Items raised in this round that are out of thesis scope

Recorded so the discussion chapter can name them as future work:

- zero-MD read of the existing `…cebonly` run (v8 `ce_channel` field):
  per-channel chord → (n, KE) map and the Q = 3 share per n;
- a Q = 2 + Q = 3 mixture arm (w_single = 0, Q = 3 pre-registered as
  "no retune") scored against the full histogram;
- per-n Q = 3 share from I²⁺ covariance per I⁺Heₙ row at 150 and 300 mW;
- interpretation of the 300/600 mW KED cross-check as a Q = 3 test;
- an arbitration objective that moves n = 1 from gate clause to held-out
  prediction (relevant to the planned re-arbitration);
- the high-n tail as a Tier-3 (fluctuation) problem, §7.

## 7. The high-n tail: a ceiling and a gap (DISCUSSION, added 2026-10-08)

No runs were made for this section. Sources: D0 §1.G4.14 (escape cliff,
well arms), §1.G9, §1.G9b; `THESIS_FINAL_GOALS.md` status item 11;
`data/reference/integrated_i_he_abundance.csv`;
`data/runs/h2b_forward_model/g9_columns.csv`, `g415_table.csv`,
`detfix_postpickup_table.csv`.

### 7.1 What the experiment shows (MEASURED, reference CSV)

Share of solvated ions (n ≥ 1):

| n | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|
| share | 0.020 | 0.019 | 0.016 | 0.014 | 0.009 | 0.007 | 0.006 | 0.005 |

- n 15–16 = 0.031, n 17–20 = 0.026 (n ≥ 15 = 0.057). The mass window ends
  at n = 20; "max n > 20" means only that the window ends there.
- Smooth decay, ≈ × 0.8 per step, with one step at 16 → 17 (× 0.60) — a
  candidate ladder edge (weaker binding above n ≈ 16), an input for M1.
- **No pile-up toward n = 20–21.**
- KE at high n is low but not near zero (n = 15: 0.075 eV; model 0.019).

### 7.2 Why the model stops at n ≈ 16–19 (MEASURED mechanism, estimated extrapolation)

- Inside the droplet pickup holds n = 21; the shell cools as
  E_int,exit = E₀·e^(−X/τ), X = ∫ρ̂ dt; after exit about one He is shed per
  ≈ 9 meV. Final n and KE are both functions of X: high n = cold exit =
  long exposure = slow ion (D0 §1.G9b).
- Measured escaper ceilings: X/τ ≈ 1.86 (lq A5) and 2.09 (h405);
  n ≥ 15 sits at X/τ ≈ 1.73 (E_exit ≈ 55 meV). Max n 16–17.
- Extrapolation (estimate, from E_exit = E₀·e^(−X/τ) and ≈ 9 meV per
  rung): n ≈ 17 needs X/τ ≈ 2.0, n ≈ 20 needs X/τ ≈ 2.7–3.0
  (E_exit ≈ 15–20 meV). The ceiling must rise ≈ 1.5× in X/τ.
- Levers that raise the ceiling, all measured:
  - faster cooling τ (lq would need τ ≈ 4.3 ps; inside the Bounded band
    [2.6, 16.5]) — trades the tail against KE₁ and n₁ (G4.14);
  - cheaper escape for slow ions (shallow well): well-only h405-type ring
    at τ 4.8 reaches n ≥ 15 = 0.046, max n 19; well-only lq A5 at τ 6.4
    reaches 0.007, max n 17. The 18 Å trace pins only the total escape
    cost, so a Tier-0-consistent shallower well is compensated by
    stronger drag and the gain vanishes (G4.15);
  - colder start E₀ — fills n 10–14 to experiment, not n ≥ 15 (G9).
- These get the model to n ≈ 19 at best. They do not give the shape.

### 7.3 The gap is the real obstacle (MEASURED gap, argued consequence)

- The model's exposure distribution has a **gap**: escapers end at
  X/τ ≤ 1.86, retained ions sit at X/τ ≥ 5 (cold, intact n = 21), nothing
  in between (D0 §1.G9b).
- Releasing the retained class does not fit: it would pile up at
  n ≈ 20–21, and its size (trap 0.041 lq / 0.087 h405) exceeds the
  experimental n 17–20 total (0.026). The data show a falling top end.
- The experiment instead needs a **continuum** of exit exposures over
  X/τ ≈ 2–3.5. A static well with deterministic dynamics makes a cliff by
  construction: at the surface an ion either has the KE to leave or never
  leaves. A smooth tail needs the escape threshold blurred.

### 7.4 Missing ingredient: collision shot noise (OUTLOOK, Tier 3)

- The drag law is the mean of discrete He collisions. Replacing the
  hard-sphere model by continuous drag removed the fluctuation channel
  (CLAUDE.md: "loss of the Langevin fluctuation channel").
- The relevant noise is collision shot noise at 2–10 Å/ps, not thermal
  noise (k_B T at 0.37 K is negligible against a ≈ 0.1 eV well). Its size
  follows from the collision kinematics (effective temperature of order
  m_He·v²), so it is tied to the drag and adds no free parameter.
- Order of magnitude (estimate): each collision removes ≈ 4–6 % of the
  ion's energy (mass ratio), so ≈ 20–40 collisions to stop; relative
  energy straggling ≈ 1/√30 ≈ 0.2.
- Near the threshold this gives a continuum of surface dwell times, hence
  exit exposures, hence a smooth n tail. The same noise widens KE at fixed
  n — the model is ≈ 8× too narrow there (IQR/median 0.09–0.11 vs
  0.71–0.82). Whether shot noise alone reaches the experimental width is
  open; the straggling estimate suggests it covers part of it.
- **Thesis statement (draft).** The size histogram is reproduced to
  n ≈ 15 and the n = 2–8 KE curve is matched. The tail beyond n ≈ 16 needs
  two things: a higher exposure ceiling for escapers (reachable within the
  model, at a KE₁ / Tier-0 cost) and fluctuations at the escape threshold
  (absent by construction in a mean-field drag model). The n ≥ 17 tail and
  the too-narrow KE width at fixed n are the same missing ingredient:
  Tier-3 noise.

### 7.5 Considered and rejected in this discussion

- **Droplet breakup weakening E_bind late** (rejected 2026-10-08). The
  deposit (≈ 3–5 eV per molecule) is below the binding of a typical
  production droplet (N̄ 12.8k, ≈ 8 eV; `detfix_postpickup_table.csv`), so
  droplets shrink rather than break. The well is local (1/r⁴), so
  shrinking changes it by ≈ 0.004 eV (50 → 35 Å, continuum estimate), and
  the retained ions sit in the large (long-chord) droplets that survive.
- **Stronger cooling only for long-residence ions.** No physical basis
  found; dropped.
- **Per-ion E₀ / f_int scatter and the Q = 3 channel** — tried
  (G4.14 projection, G9, D0 §20); see `THESIS_FINAL_GOALS.md` status 9–11.
- **Spin-orbit-state-dependent ladder.** Ruled out on energy: a fast
  escaper leaves with E_int ≈ 0.19 eV and would need ≈ 30 meV per He to
  stay at n ≥ 15, about three times the ladder.

## 8. Mass influence (DISCUSSION; added 2026-10-13, `THESIS_FINAL_GOALS.md` G11.5; zero MD)

Discussion, not validation: with the 18 Å-only extraction at the true
shell mass nothing mass-related is left to validate (Tier 1 retired,
G11.3); what remains is model scope.

| level | question | result (MEASURED) | source |
|---|---|---|---|
| extraction | sensitivity of the fitted law to the fit mass | 19 → 21 He: a −1.4 to −2.3 % per c, E_bind +0.0004 eV, 9 Å margin −0.00 to −0.04 Å/ps; B2 / A5 pass unchanged | G11.4 step 2, D0 §1.G4.12 |
| transfer | fixed-mass law under production's variable mass | exact at birth (n₀ = 21 = fit mass); afterwards rests on the mass-agnostic force law — **stated assumption:** the drag force does not depend on the shell size n (untestable: no small-n TDDFT trace) | G11.2 (3) |
| initial shell | ions born bare (n₀ = 0) | below (outlook) | G8, D0 §1.G8 |

**Bare start (OUTLOOK).** The drag law is extracted for one complex
only, I⁺ + 21 He. G8 applies it unchanged to ions born bare (n₀ = 0,
E_int(0) = 0; N = 1000, s1; lq A5 cell = the ridge point next to B2,
h405 as contrast):
- bare ions do not evaporate; n = pickup count (peak n = 2–3);
- KE₁: lq 0.491 eV, below the intact 0.621 and far below the
  experimental median 1.128 — the lighter ion's larger Coulomb speed is
  taken back by the ridge law's v² term;
- h405 (capped high-speed drag) keeps the gain (KE₁ 0.865), still below
  the experiment.

Reading: with a law extracted for the 21-He complex, a bare start cannot
supply the fast n = 1 ions. This is **consistent with** the law
over-dissipating on an ion it was not extracted for (a bare I⁺ has a
different solvation structure), **not shown**: G8 cannot separate the
drag from pickup and the lighter flight mass. Outlook: a drag law for
small n (TDDFT for I⁺Heₙ, n ≪ 21) before a bare / partly-solvated source
channel can be tested; ties to the n = 1 outlook (§3, §5; RQ1 / RQ8, M1).

## 9. Coulomb budget: production at E_C, the source calibration at 0.8·E_C (DISCUSSION; added 2026-10-13; zero MD)

Sources: D0 `TIER2_PARAMETER_INFLUENCE.md` §19; `CALIBRATION_MAP.md`
row 15; `RESEARCH_QUESTIONS.md` RQ7 / RQ8 NBs 2026-07-29; log 2026-07-29.

**What production uses (MEASURED, cfg.json).** Every Tier-2 production
run (h405, A5, B2, B3) runs `E_coulomb_scale` = 1.0 at R₀ = 2.666 Å: the
point-charge Coulomb energy, 2.70 eV per I⁺ fragment for the Q = 2
channel. It is set explicitly by `tier2_common.py`; the legacy
droplet-distribution preset's 0.8 is deliberately not inherited.

**The source calibration (SOURCED).** Hatherly et al., J. Phys. B 27
(1994) 2993 (doi:10.1088/0953-4075/27/14/032), gas-phase I₂ in intense
fields: channel KERs are a channel-independent fraction of the
point-charge E_C (0.75·E_C at 90 fs, 0.65·E_C at 200 fs), because
ionization is sequential during the separation (stretched states). Not
helium screening. The legacy `E_coulomb_scale = 0.8` (gas and droplet
presets) is this calibration; it gives 2.16 eV (Q = 2) / 4.32 eV (Q = 3)
per fragment, in line with the measured gas-phase peaks 2.26 / 4.11 eV
(`vmi_iplus_gas.csv`, ≈ 5 %). Against it the production budget is
**≈ 25 % high**.

**What changing it does (MEASURED, D0 §19 probe 2026-07-29; at the h405
pins, E₀ pinned, N = 1000, CRN-paired — not re-measured at B2).**
- KE₁ responds at S_k = dKE₁/dbudget = 0.389 eV/eV: budget 2.26 / 2.70 /
  4.11 eV → KE₁ 0.486 / 0.637 / 1.200 eV.
- A uniform re-anchor to the calibrated ≈ 2.2 eV breaks the landing:
  n̄ 5.33, n₁ 0.126, trapping doubled (0.179), W₁ 1.365 — a full
  re-arbitration (E₀, c) would be needed.
- A single budget cannot widen the KE at fixed n (KE₁ SD 0.034–0.047):
  the needle (§1.1) is not a budget effect.

**Why production stays at 1.0 (decision state).** The re-anchor was left
as a user adjudication in 2026-07 and never taken; the thesis keeps 1.0.
Arguments for the thesis text:
1. **It is the generous bound for KE₁.** With ≈ 25 % more Coulomb energy
   than the source calibration, the model still reaches only about half
   the experimental KE₁ (0.609 vs 1.302 eV at B2). With the calibrated
   budget the deficit would grow (projection from S_k: ≈ −0.17 eV at
   2.26 eV). The KE₁ deficit is therefore not a too-small Coulomb budget.
2. **Re-anchoring is not a one-knob change**: it moves the whole landing
   and would re-open the (c, E₀) arbitration of G5.
3. **Caveat to state:** the production budget is a known overestimate
   relative to the source's own calibration; the selected (c, E₀) are
   conditional on it, as they are on the sampler (§2).

**Outlook.** The calibrated budget belongs together with the Q = 2 /
Q = 3 channel mixture (§3): the 0.8·E_C calibration sets both channel
energies (2.16 / 4.32 eV), and a re-arbitration at the calibrated budget
is only meaningful with both channels present.
