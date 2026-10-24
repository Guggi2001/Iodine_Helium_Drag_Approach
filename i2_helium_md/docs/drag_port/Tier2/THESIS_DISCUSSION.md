# Thesis discussion points (final discussion round, 2026-10-02)

Material for the discussion chapter of the thesis. These are physical
arguments and open questions, **not** goals: the Q = 3 / channel-mixture
line is out of scope for the thesis (user decision, 2026-10-02) and appears
there at most as discussion. Measured results and projections are labelled
as such. No runs were made for this document.

Sources read: D0 `TIER2_PARAMETER_INFLUENCE_h405_finc1v725.md` §14.3 (2026-10-02 note),
§15.8, §20; `drag_migration_log_tier2.md` 2026-10-02 entry;
`data/runs/h2b_forward_model/detfix_postpickup_table.csv` and
`detfix_postpickup_chord.csv`; `TIER2_DetectorStageFix.md` §5.7;
`CALIBRATION_MAP.md` row 30; `RESEARCH_QUESTIONS.md` (NBs of 2026-07-29);
`data/reference/README.md`, `ihe_ked/README.md`, `paper_cov/README.md`.
Later additions: §0 storyline, §8 mass influence and §9 Coulomb budget (2026-10-13, from
`THESIS_FINAL_GOALS.md` G5 / G11); §1 is to be refreshed by G6 (chord map
at B2) and G7 (KE₁ geometry-invariance caveat). §10 (dissociation
timescale and the RRK s_eff, 2026-10-14) is the one section backed by new
runs made for it (s_eff scan + two N = 6000 pools at B2). §12 (exit vs
detection readout, 2026-10-16) is a design decision for future work; no
runs. §13 (over-alignment: model β = 2 vs experiment ≈ 0.3, the 2-D /
3-D comparison and the corrected VMI inversion, 2026-10-23;
`THESIS_FINAL_GOALS.md` G12) is zero-MD.

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
with the old-bundle map (h405 `post_pickup`, binned) overlaid, (c) B2 chord
distribution (the raw-sampler A5 overlay was dropped 2026-10-24: its runs
are no longer on disk). Numbers: `g5_postpickup_chord.csv` (binned,
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
    distribution (B2 and A5 `post_pickup` had bit-identical chords,
    checked 2026-10-13 — the law does not enter; `post_pickup` shifts weight from
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
- KE width at fixed n, model vs experiment (IQR / median; corrected KED
  reference v2, G12; refreshed 2026-10-24 and extended to n = 1–6, the
  B2-3 KED panels):

  | n | model (B2) | experiment (v2) | ratio |
  |---|---|---|---|
  | 1 | 0.101 | 0.847 | 8.4 |
  | 2 | 0.086 | 0.730 | 8.5 |
  | 3 | 0.092 | 0.733 | 8.0 |
  | 4 | 0.090 | 0.733 | 8.2 |
  | 5 | 0.116 | 0.773 | 6.6 |
  | 6 | 0.152 | 0.789 | 5.2 |

  The model is ≈ 8× too narrow in relative width at n = 1–4, falling to
  ≈ 5× by n = 6 (the model width grows with n, the experimental one
  barely); consistent with §7.4, measured there at the earlier law. (v1
  values, 2026-10-13: experiment 0.69–0.80, ratio 7.7–8.1 at n = 1–4; the
  model side is unchanged.) In absolute terms the gap at n = 1 is larger
  (σ 0.67 vs SD 0.038 eV, ≈ 18×) because the model median is about half
  the experimental one. Experimental quartiles are read from the
  unsmoothed 3-D P(E) (`IHe_KED_curves_n*.csv`, v2); their medians sit
  4–11 % above the reference median column (cuts not re-done), so only
  the relative width is taken from them.
  Caption sentence (B2-3): the simulated distributions are 5–8× narrower
  than the experiment (IQR / median 0.09–0.15 vs 0.73–0.85, n = 1–6;
  ≈ 8× for n ≤ 4).

**Statement for the thesis.**
1. "KE₁ is geometry-invariant" holds **within the deterministic model**:
   geometry moves ions along the chord axis, the n = 1 window sits at a
   fixed chord, and the map through it is near-deterministic, so KE₁
   cannot change.
2. The experimental KE width at fixed n shows the real map **scatters**
   (≈ 5–8× wider). With scatter the n = 1 window gets blurred edges and a
   chord-window selection can bias KE₁. Size (PROJECTION, not measured):
   at most the KE slope × a shift of the window's mean chord by up to half
   its width, −0.021 eV/Å × 2.7 Å ≈ 0.06 eV — small against the
   0.45–0.6 eV KE₁ deficit (v2: exp median 1.05 / mean 1.20 vs model
   0.61 eV), so the deficit is not a geometry artefact of the
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
- the high-n tail as a Tier-3 (fluctuation) problem, §7;
- the exit-versus-detection comparison (post-exit evaporation), readout
  defined in §12.

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

Sources: D0 `TIER2_PARAMETER_INFLUENCE_h405_finc1v725.md` §19; `CALIBRATION_MAP.md`
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

## 10. Over what time does dissociation happen? The RRK s_eff at B2 (MEASURED; added 2026-10-14)

Sources: runs made for this section (code `scripts/gen_tier2_b2_seff_scan.py`,
scorers `scripts/post_processing/tier2_b2_seff_table.py` and
`tier2_b2_seff_timing.py`; tables `data/runs/h2b_forward_model/b2_seff_table.csv`
and `b2_seff_timing.csv`; commits cef04ee, b47778e, d39da86, 52b519b).
Reference: Albrechtsen, Schouder, Viñas Muñoz, Christensen, Engelbrecht
Petersen, Pi, Barranco, Stapelfeldt, *Nature* 623, 319 (2023),
doi:10.1038/s41586-023-06593-5 (Methods, "Molecular dynamics simulations of
the dissociation time of Na⁺He_N complexes"; Extended Data Fig. 4).
Earlier in-model record: findings I10, I11, I13, I20, I32
(`TIER2_STAIRCASE_PROBE_FINDINGS.md`), RQ6.

### 10.1 The question

The He shell sheds by energy-gated RRK evaporation,
k = ν·(1 − D₀(n)/E_int)^(s−1), ν = 2.42 ps⁻¹ (`physics/evaporation.py`).
Production uses a constant effective dof s_eff = 8 (`evap_rrk_dof`,
CALIBRATION_MAP row 10). The earlier reading (I13, I20) was that s sets the
**clock** of the cascade, not its endpoint, and that the detector is nearly
s-blind. That was measured before the lq ridge, the `post_pickup` sampler
and the absolute-E₀ arm; it had never been checked at B2. Two questions:
(1) how long does dissociation take in the model, compared with the
literature statement that it is essentially over within tens of
picoseconds; (2) does s_eff change the detected observables at B2.

### 10.2 Literature: dissociation is fast (SOURCED)

Albrechtsen et al. 2023 justify that their detector records final sizes by
classical MD (TINKER, AMOEBA-09) of isolated Na⁺He₆ and Na⁺He₁₄, heated to
20 / 50 / 100 K, 1 ns NVE, 1000 runs each: "dissociation takes place
almost exclusively during the first tens of picoseconds"; the timescale is
"mostly unaffected" by size and energy; a slow decrease of ⟨N⟩ near
threshold ("especially evident for total kinetic energies close to the
binding energy") is attributed to rare events; "in a conservative estimate,
the vast majority of dissociation events are finished after 400 ps".

Transfer caveats: Na⁺, not I⁺ (different D₀ ladder and shell); classical
MD keeps every mode, whereas a low s_eff stands for quantum mode freezing,
so the paper contains no number that maps onto s; their complexes start
well above threshold, where the rate is fast for any s. The paper
constrains the **timescale**, not s directly. Its slow near-threshold tail
is the same effect as the RRK slowdown below.

### 10.3 Mechanism: why s can move the endpoint at finite flight time (DERIVED)

After ejection the cooling is density-gated off, so E_int falls only by
one D₀ per shed. Near threshold the bracket (1 − D₀/E_int) is small and
the exponent s − 1 collapses the rate. An ion effectively stops shedding
once 1/k exceeds the flight time t_detect = 8.53 µs = 8.53·10⁶ ps, i.e.
k < k_min ≈ 1.2·10⁻⁷ ps⁻¹. The effective freeze threshold is

E*/D₀ = 1 / [1 − (k_min/ν)^(1/(s−1))]:

| s | ion stops shedding below E_int ≈ |
|---|---|
| 5 | 1.02 D₀ |
| 8 | 1.10 D₀ |
| 12 | 1.28 D₀ |
| 20 | 1.70 D₀ |
| 30 | 2.27 D₀ |

For s ≲ 8 this is essentially the energetic floor E_int = D₀: the cascade
runs to its endpoint and s only sets the clock (the earlier reading holds).
For large s an ion can carry more than a full extra rung and still never
shed again before the detector: s then moves the detected n. A second,
smaller channel acts inside the droplet: a slower cascade gives the
Newton cooling more time to drain E_int before the shed (the race, I13),
and it lowers the number of pickups because n stays high under the
occupancy cap. The 1/k > t_detect criterion is an order-of-magnitude
threshold (shedding is a stochastic waiting time).

### 10.4 s_eff scan at B2 (MEASURED; seed s1, N = 1000 per cell)

The B2 s1 member with only `evap_rrk_dof` moved (cfg diff = that key;
s_eff = 8 is the committed B2 s1 run; scorer anchor = its committed
`g5_postpickup_table.csv` row, 4 dp).

| s_eff | n̄ | n₁_solv | W₁ | midHot | χ²_med | KE₁ [eV] | frozen | still shedding at detector |
|---|---|---|---|---|---|---|---|---|
| 5 | 3.74 | 0.233 | 1.01 | 0.957 | 250 | 0.606 | 0.82 | 0.01 |
| **8 (B2)** | **3.94** | **0.214** | **0.93** | **0.994** | **292** | **0.609** | **0.77** | **0.07** |
| 12 | 4.21 | 0.183 | 0.90 | 1.05 | 406 | 0.615 | 0.65 | 0.20 |
| 20 | 4.68 | 0.121 | 0.98 | 1.14 | 612 | 0.629 | 0.37 | 0.49 |
| 30 | 5.28 | 0.030 | 1.32 | 1.25 | 998 | 0.650 | 0.04 | 0.82 |

"Still shedding" = `time_exhausted` (next shed scheduled after 8.53 µs);
fractions over all ions (suppressed ≈ 0.11 at every s).

**He shed-event timing** (every shed of every ion over the whole flight:
MD window 0–500 ps from the stored `n_shell` history, detection stage
500 ps – 8.53 µs from the event record; re-runs with `ion.npz` kept,
detection bit-identical to the scan):

| s_eff | sheds per ion | median shed time | sheds after 50 ps | after 100 ps | after 400 ps |
|---|---|---|---|---|---|
| 5 | 14.8 | 7.8 ps | 2.2 % | 1.5 % | 1.0 % |
| **8 (B2)** | **14.5** | **9.0 ps** | **5.2 %** | **3.7 %** | **2.3 %** |
| 12 | 14.2 | 10.7 ps | 10 % | 7.1 % | 4.4 % |
| 20 | 13.7 | 17 ps | 23 % | 15 % | 9.1 % |
| 30 | 13.1 | 32 ps | 40 % | 28 % | 16 % |

Pickups inside the MD window fall with s (1557 → 1194 → 876 → 548 → 333
per 2000 ions). Undercount disclosed: a pickup and a shed inside the same
0.16 ps storage interval cancel. (The fraction of *ions* with at least one
shed after 400 ps is larger — 14 / 34 / 61 / 84 / 86 % — because one late
near-threshold shed per ion suffices; the paper's statement is about
events.)

**Reading.**
1. **Dissociation in the model is fast at production s:** at s_eff = 8
   the median shed happens at 9 ps and 97.7 % of all sheds are over by
   400 ps — consistent with the supplementary classical MD of
   Albrechtsen et al. ("first tens of picoseconds"; "vast majority …
   finished after 400 ps"; Methods + Extended Data Fig. 4, Na⁺He_N —
   **not** an experimental timescale; their experiment, Yₙ(t) vs Pₙ(t),
   only shows whether dissociation occurs, not when). s_eff = 5–12
   pass both statements; s_eff = 20 is borderline (23 % after 50 ps);
   s_eff = 30 fails (40 % after 50 ps, 16 % after 400 ps).
2. **"s sets the clock" holds for ≈ 98 % of the sheds**: the shed count
   per ion barely moves (14.8 → 13.1 over s 5 → 30). What changes the
   detected n is the small tail of last, near-threshold sheds (§10.3):
   pending at 8.53 µs, or late enough to suppress pickups in the droplet.
3. **The detector is not s-blind at B2 above s ≈ 8**: n₁ collapses
   (0.214 → 0.030) and χ² triples from s 8 to 30. The experiment rejects
   s ≥ 20 (n₁ 0.12 / 0.03); s = 12 is disfavoured, not excluded (n₁ 0.183
   vs the gate edge 0.19 ≈ 0.7 single-seed sd; χ² +40 %; one seed).
4. **KE₁ is insensitive** (0.606–0.650 eV over the whole range): s_eff is
   not a route to the KE₁ deficit (flight mass, `TIER2_MASS_SCENARIOS.md`).

**Explored range, not an established bound:** only s_eff = 5, 8, 12, 20,
30 were run at B2; **s < 5 and 8–12 are untested**. The size distribution
excludes s_eff ≥ 20 and disfavours 12; the cross-system MD timescale
(plausibility check only) fails s = 30. 5–8 is the plausible explored
range; production s_eff = 8 lies in it and is not uniquely selected.

### 10.5 What a switch to s_eff = 5 achieves (MEASURED; pooled N = 6000, s1–s5 + r6)

Three pools at the B2 cell, identical seeds, production pipeline; each
cfg differs from B2 in only the listed keys (oracle). The E₀ step
restoring n̄ came from the committed B1 E₀ column (dn̄/dE₀ ≈ −32 /eV):
ΔE₀ = −0.182/32 = −0.0057 eV → 0.3240 eV.

| pool | s_eff | E₀ [eV] | n̄ | n₁_solv | midHot | χ²_med | KE₁ [eV] | deepKE | W₁ | n ≥ 15 | suppressed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| B2 | 8 | 0.3297 | 3.903 | 0.211 | 1.002 | 527 | 0.609 | 0.477 | 0.967 | 0.0008 | 0.111 |
| s5 | 5 | 0.3297 | 3.721 | 0.226 | 0.960 | 461 | 0.607 | 0.512 | 1.044 | 0.0002 | 0.111 |
| s5, E₀ compensated | 5 | 0.3240 | 3.909 | 0.210 | 0.994 | 442 | 0.628 | 0.556 | 1.066 | 0.0002 | 0.084 |

Pooled containers (figures in each `figures/`):
`…N6000_detfix_conf270_g5lqppc160e3297seff05pooled_th500` and
`…g5lqppc160e3297seff05e3240pooled_th500` (build targets
`g5_lqppb2seff05` / `g5_lqppb2seff05e3240`; builder verified by rebuilding
the committed B2 pool array-for-array).

**Reading.**
1. **At matched n̄, n₁ and midHot, s_eff = 5 is not an E₀ shift.** The
   KE side improves — χ²_med −16 % (442 vs 527), KE₁ +0.019 eV, deepKE
   0.556 vs 0.477 — while the size histogram worsens: W₁ 1.066 vs 0.967
   and the n ≥ 15 tail nearly empties (0.0002 vs 0.0008). The suppressed
   fraction drops 0.111 → 0.084; that comes from the lower E₀, not from
   s (s = 5 alone leaves it unchanged).
2. **The same trade-off as τ** (G4.14: τ trades the large-n tail against
   KE at every c): s_eff is a second lever on that one KE ↔ tail axis,
   not a way to fit both. The experiment does not choose between s_eff 5
   and 8; the preference depends on weighting χ² / KE against W₁ / tail.
3. **KE₁ moves by 0.02 eV**, against a deficit of ≈ 0.7 eV to experiment
   (1.302 ± 0.071 eV): the KE₁ question stays with the flight mass.
4. **Thesis decision state:** production stays at s_eff = 8 (B2 unchanged;
   no re-arbitration). s_eff is reported as a convention within the
   explored range (5–8 consistent; ≥ 20 excluded by the size distribution;
   timescale consistent with cross-system classical MD; s < 5 and 8–12
   untested).

Noise: the W₁ and KE₁ differences (≈ 0.10 and ≈ 0.02) are several pooled
seed errors (member sd 0.036 / 0.003 → pooled ≈ 0.015 / 0.001), but those
member sd were measured at c = 2.0 (G4.14), not at B2. The s scan in §10.4
is single-seed.

### 10.6 Text for the discussion chapter (suggested)

- In the model, He dissociation is essentially complete within the first
  tens of picoseconds (median ≈ 9 ps, ≈ 98 % of sheds before 400 ps at
  production), in line with the supplementary classical MD estimate of
  Albrechtsen et al. (2023) for Na⁺He_N (a cross-system plausibility
  check, not an experimental timescale). The remaining late sheds are
  near-threshold events whose RRK rate collapses as (1 − D₀/E_int)^(s−1).
- The effective RRK dof sets the clock of the cascade, but at a finite
  flight time (8.53 µs) it also sets how far above threshold an ion can
  stay without shedding again; for s_eff ≥ 20 this moves the detected
  size distribution measurably and the experiment rejects it (s_eff = 12
  disfavoured).
- Within the explored range (s_eff 5–8 consistent; s < 5 and 8–12 not
  tested), s_eff acts on the same
  KE ↔ large-n-tail trade-off as the cooling time τ; it neither resolves
  the KE₁ deficit nor the missing n ≥ 17 tail.
- Outlook: a sourced s_eff for I⁺Heₙ (quantum mode counts / He-cluster
  heat capacities; RQ6) and an I⁺-specific dissociation-time simulation
  would turn the explored range into a bound or a value.

## 11. h405 at the production settings: a like-for-like contrast to B2 (MEASURED; added 2026-10-14)

Sources: `scripts/gen_tier2_h405_pp35.py`, scorer
`scripts/post_processing/tier2_h405_pp35_table.py`, table
`data/runs/h2b_forward_model/h405_pp35_table.csv` (commits 9110806,
0379a1d); pooled containers `…shared_pure_cubic_N6000_detfix_conf270_h405pp35e4310pooled_th500`
(E₀ only) and `…h405pp35v540e4240pooled_th500` (re-landed), figures in each
`figures/`.

The contrast law h405 was quoted with the `raw` sampler and the 14.2 Å
helium-boundary width, so B2 vs h405 differed in law, sampler and gate at
once. h405 has now been run at **both production settings**
(`post_pickup`, s_ρ 3.5 Å as `h405g35`) and re-landed with its own knobs:
E₀ alone restores n̄ but leaves midHot at 0.852; lowering the cap v_c
5.5 → 5.40 Å/ps (τ unchanged at 4.4 ps; τ 4.0 lifts n̄ by ≈ 1 He) and
E₀ 0.405 → 0.424 eV land midHot as well. h405 therefore **also lands under
the production settings**, so B2 and h405 now differ in the drag law only.

| pooled N = 6000, `post_pickup`, 3.5 Å | n̄ | n₁_solv | W₁ | midHot | χ²_med | KE₁ [eV] | deepKE | n ≥ 15 | trapped | suppressed |
|---|---|---|---|---|---|---|---|---|---|---|
| B2 (lq c 1.6, τ 6.4, E₀ 0.3297) | 3.903 | 0.211 | 0.967 | 1.002 | 527 | 0.609 | 0.477 | 0.0008 | 0.049 | 0.111 |
| h405 re-landed (v_c 5.40, τ 4.4, E₀ 0.424) | 3.941 | 0.189 | 0.934 | 1.007 | 704 | 0.636 | 0.537 | 0.0040 | 0.128 | 0.220 |
| *h405 as quoted (`raw`, 14.2 Å; v_c 5.5, E₀ 0.405)* | 3.890 | 0.209 | 0.755 | 0.948 | 456 | 0.640 | 0.571 | — | 0.080 | 0.192 |

**Reading.** At matched n̄ and midHot the two laws land comparably, each
with its own trade-off: h405 is slightly hotter (KE₁ +0.027 eV, deepKE)
with a marginally better W₁ and more n ≥ 15 tail; B2 has the better n₁
(h405 just below the 0.19 gate edge) and χ², and loses half as many ions
to trapping and suppression (h405 traps 13 % vs 5 %). h405's earlier W₁
advantage (0.755) came from the `raw` / 14.2 Å settings, not from the law.
This sharpens the identifiability argument (§0, G5): two differently shaped
drag laws reproduce the detected I⁺(He)ₙ observables comparably **under
identical settings**, and both leave KE₁ at about half the experiment.
The thesis contrast row is the re-landed h405, not the `raw` / 14.2 Å
one. Caveat: h405's 14.2 Å-extracted bundle under the 3.5 Å gate is a
production setting, not a Tier-0 refit (as for `h405g35`).

## 12. Exit versus detection: the readout definition (DESIGN DECISION, future work; added 2026-10-16; no runs)

Decision (user, 2026-10-16): "droplet exit" for the exit-versus-detection
comparison is the state at an **actual outward crossing of the sharp
droplet radius, r − R = 0**, not the state at the fixed handover time
t_h = 500 ps. Not executed for the thesis; this section fixes the
definition so that the discussion chapter can state it and future work
can implement it unchanged.

### 12.1 Why not the handover state

The two readouts are different comparisons:

- **Evaporation continues outside the droplet.** Only pickup, drag and
  cooling are density-gated (the P1–P3 handover guard,
  `simulation/detection_stage.py` module docstring); the RRK shedding is
  not. Sheds between the physical crossing and t_h would be booked as
  "in-droplet" by a handover readout.
- The handover state is a pipeline boundary (Stage I end), not a physical
  event; it would change with t_h.

### 12.2 Definition

- **Boundary:** depth dᵢ(t) = |rᵢ(t)| − Rᵢ = 0, with Rᵢ the per-atom
  droplet radius (`droplet_radii_angstrom`): the same depth and the same
  sharp radius as `b_ion_outside` (`simulation/ion.py`
  `_write_final_state`). On the 3.5 Å erf profile this is the half-density
  point, ρ̂(0) = 0.5 (`physics/_gates.py` `_erf_complement`): the sharp
  radius is a stated convention, not the point where the helium coupling
  vanishes (drag, pickup and cooling are still active just outside it).
- **Exit time t_exit,i:** the **last** outward crossing (dᵢ ≤ 0 → dᵢ > 0)
  in [0, t_h], so an ion that leaves and falls back is read at its final
  departure.
- **Exit state:** (n, E_int, m, v, KE) at the first stored column after
  that crossing; time resolution = the ion stage's stored stride.
- **Unit:** each iodine fragment separately (2N rows), not the
  per-molecule OR rule of `b_ion_outside`.

### 12.3 Which fragments enter both samples

- **Both samples are the same ions:** the detected ensemble,
  `DetectionResult.detected_mask` (state `frozen`, `suppressed`,
  `time_exhausted`). Every such ion has left the droplet before t_h (the
  P1–P3 guard and the closure's well-plateau check require it).
- **Excluded from both:** `droplet_retained` and
  `droplet_retained_marginal` (`RETAINED_REASONS`); their fractions are
  reported separately (B2: trapped 0.049, no D4 firing).
- `suppressed` ions are included in both; they keep their exit n only if
  no shed occurs after exit.

### 12.4 What is compared

- **Primary:** n at exit vs n at detection. Their difference is the
  post-exit evaporation: sheds in (t_exit, t_h] in the ion stage plus the
  detection-stage events (all at t > t_h).
- **Secondary, not an evaporation measure:** KE at exit vs detection. The
  residual pair Coulomb repulsion (booked at detection by the
  partner-aware closure) and the remaining well climb also change KE after
  exit, so a KE difference mixes Coulomb work with evaporation.

### 12.5 What it costs

- The B2 member directories hold `cfg.json`, `neutral.npz`,
  `detection.npz` only (`ion.npz` dropped after detection verified), so
  the crossing read needs the six members (s1–s5 + r6) re-run with
  `ion.npz` kept (deterministic seeds, ≤ 2 concurrent) and a small
  post-processing read of the crossing (code: `[PROCEED TO
  IMPLEMENTATION]`).
- Zero-MD by-product (not the exit comparison): n and E_int at **t_h** are
  exact from `detection.npz` alone, n(t_h) = n_det + number of events and
  E_int(t_h) = E_int,det + Σ D₀(n_pre) over the events.

---

## 13. Over-alignment: the model keeps the laser alignment, the experiment does not (MEASURED + DISCUSSION; added 2026-10-23; zero MD)

Origin: the 2-D KED figure B2-4 (`tier2_b2_thesis_figure.py`,
`b2_4_ked_shapes_2d.png`) puts every simulated 2-D peak to the right of
the experiment, while the 3-D overlay (B2-3) looks fine. Tracing this led
to the angular distribution. Goal record and plan: `THESIS_FINAL_GOALS.md`
G12. All numbers below are zero-MD reads of the B2 pool
(`…g5lqppc160e3297pooled_th500`) and of the VMI fragment scan of
17.10.24; scratch probes, not yet committed. Experimental numbers use the
**corrected** inversion (§13.2, G12); where the committed (v1) reference
differs this is said.

### 13.1 Geometry and how alignment is measured

- **Polarisation = image X (horizontal) = droplet-beam axis.** Evidence:
  (i) gas-phase I⁺ from free I₂ (same laser) has lobes along X, up to
  8 : 1 (strong-field ionisation of I₂ selects molecules along the field,
  ions fly along the polarisation); (ii) the per-fragment image centres
  shift along X with fragment mass (droplet-beam velocity × flight time),
  so X is the beam axis, and in the standard crossed-beam geometry the
  in-plane polarisation is along the beam; (iii) cues in the
  predecessor's thesis (user, 2026-10-23); (iv) the rotated inversion
  reproduces the measured 2-D curves slightly better (n 0–6).
- **Anisotropy parameter.** I(θ) ∝ 1 + β·P₂(cos θ), θ from the
  polarisation; β = 2 is pure cos² (all ions along the field), β = 0
  isotropic, β = −1 pure sin². Model β: measured directly from the final
  velocities (5⟨P₂⟩ about lab x). Experimental β: pyabel `rIbeta` on the
  inverted image, per speed bin, P(v)-weighted over the ⟨E⟩ signal mask
  (`IHe_KED_beta.csv`, `IHe_KED_beta_vs_v.csv`).
- **Two pitfalls found on the way** (both matter for reading old numbers):
  - pyabel's symmetry axis is the image *vertical*. For weak anisotropy
    an inversion about the wrong in-plane axis returns β ≈ −β_true
    (synthetic test, exact to ±0.01). The β sign alone therefore cannot
    tell which axis is right.
  - an image-plane fit a + b·cos²(φ − φ₀) with free φ₀ measures only the
    *strength*, not the sign, of the anisotropy (φ₀ absorbs a 90° turn).

### 13.2 Side result: the 3-D reference was inverted about the wrong axis and centre (MEASURED)

The same investigation found that the committed `ihe_ked` 3-D curves
and ⟨E⟩ were inverted (i) about the wrong symmetry axis (vertical,
polarisation horizontal) and (ii) about the array centre instead of the
per-fragment centre (`recenter_VMI` relabels the X/Y grids only; pyabel
receives the bare matrix; offset 3–5 grid px). Synthetic tests with the
exporter's call: (i) ⟨E⟩ ≈ +5 % at β 0.3, σ ≈ −10 %; (ii) ⟨E⟩ +1 % (n 1)
to +5 % (n 17), σ +2 to +9 %. The corrected export lowers ⟨E⟩ by 3–8 %
coherently for every n (n 1 −8.0 %); the 2-D columns are identical (never
inverted). **Adopted as reference v2 (2026-10-23).** Re-scored
(`g5_winners_table_ihekedv2.csv`): B2 midHot 1.002 → 1.050, deepKE
0.477 → 0.499, χ² 527 → 350, KE₁ model / exp 0.609 / 1.198 (deficit
53 % → 49 %); B3 midHot 1.016 → 1.065, deepKE 0.601 → 0.628. Sizes
unchanged. KE numbers elsewhere in this document that predate
2026-10-23 are v1. Table: `THESIS_FINAL_GOALS.md` G12.3.

### 13.3 What is measured (MEASURED)

**Model vs experiment, per fragment** (β along the polarisation):

| n | model β (B2) | experiment β (corrected) |
|---|---|---|
| 0 (bare I⁺) | 1.88 | 0.05 (core 0.16) |
| 1 | 2.04 | 0.30 |
| 2 | 1.97 | 0.25 |
| 3 | 1.98 | 0.25 |
| 4 | 1.98 | 0.33 |
| 5 | 1.83 | 0.29 |
| 6 | 1.93 | 0.31 |
| 7 | 2.04 | 0.34 |
| 8–17 | — | 0.22–0.34 |

- The model sits at the cos² limit for every n (β 1.83–2.04; image-plane
  fit 1.89–2.02). The experiment is weakly aligned, **flat in n** from
  n = 1 to 17. Experimental n = 0 is droplet-born I⁺ only: the
  gas-phase background (run 43635, no He) is subtracted; caveat — an
  over/under-subtraction of the strongly aligned gas I⁺ would move it.
- **Gas phase vs droplet** (raw images, horizontal / vertical intensity,
  crude ring sectors): gas I⁺ 8.3 / 5.2 / 3.3 at 160 / 300 / 600 mW
  (alignment falls with power: saturation); droplet-beam I⁺ before
  background subtraction 5.3 / 3.1 (160 / 600 mW; carries the gas
  component); I⁺He 1.28 / 1.38; I⁺He₂ ≈ 1.2. The free molecule is
  strongly aligned; everything that comes out of the droplet is not.
- **β inside a fragment depends on speed** (corrected export, P(v)
  terciles slow / mid / fast): n 1 +0.33 / +0.29 / +0.28; n 4 +0.64 /
  +0.40 / −0.01; n 6 +0.83 / +0.31 / −0.06; n 7 +0.62 / +0.42 / +0.04;
  n 0 −0.05 / +0.11 / +0.09. For n ≥ 3 the slow ions keep a clear
  alignment and the fast ones are isotropic. **Not an inversion
  artefact** (synthetic check 2026-10-23, exporter's pyabel call):
  a constant β = 0.3 on log-normal shells at the n = 2 / 4 / 7 radii and
  widths, 0.5 px residual offset, comes back flat — noise-free
  +0.30 / +0.30 / +0.31; at the real counts (N_eff, 5 seeds)
  e.g. n 4 +0.28 ± 0.03 / +0.27 ± 0.07 / +0.31 ± 0.04, n 7 +0.27 ± 0.11 /
  +0.27 ± 0.10 / +0.31 ± 0.06. The measured slow–fast differences (n 4
  0.65, n 7 0.58) are ≈ 10σ and ≈ 4–5σ. Not excluded: structure in
  the background subtraction (global background 43655, scaled in a
  4–8 eV window), which the synthetic test does not contain.
  Across fragments the same pattern appears: the fastest fragment (bare
  I⁺, ≈ 3.6 eV) is the least aligned.

### 13.4 Why the model is at β = 2 (MECHANISM, by construction)

1. **Initial orientation.** `sampling/orientations.py` draws the I–I axis
   with probability ∝ cos²(angle to the polarisation, lab x) — the
   one-photon/strong-field selection law (legacy recipe, PRA 1993). This
   alone gives β = 2 for axial recoil. Direction correct (§13.1).
2. **Axial recoil.** The Coulomb explosion pushes both ions along the
   I–I axis.
3. **No deflection afterwards.** The drag force is antiparallel to the
   velocity (F = −γ(v)·v), so it slows an ion without turning it; the
   fluctuation channel is off (noise zero). Nothing in the flight changes
   the direction, so the initial cos² arrives at the detector unchanged.
4. The legacy hard-sphere model used the same sampler; that the legacy
   runs showed the same 2-D offset is consistent with this.

### 13.5 Consequences

**(a) Energy / speed observables are unaffected (ARGUED, exact in the
model).** All forces depend on relative geometry, the droplet is
spherical and positions are drawn isotropically, so the model is
rotation-invariant: the orientation distribution changes only where each
ion points, never its speed, n or fate. Every Tier-2 score (size
histogram, KE per n, midHot, deepKE, χ²) is independent of the β
mismatch. **The 3-D comparison against the corrected reference is clean
and is the one the thesis presents.**

**(b) The 2-D comparison mixes speed and angle (MEASURED decomposition).**
A projection onto the detector depends on the speed distribution *and*
the angular distribution. Exact forward projection (shell kernel, line of
sight z, polarisation in-plane; validated: the kernel at β = 2
reproduces the model's direct √(vx² + vy²) histogram to ≤ 0.1 Å/ps).
Medians in Å/ps, n = 4–7:

| n | exp 2-D | model as is | + exp β | + exp β + exp width | left over | 3-D: model / exp (corr.) |
|---|---|---|---|---|---|---|
| 4 | 5.61 | 6.85 | 6.45 | 5.95 | +0.34 (+6 %) | 7.35 / 7.07 (+4 %) |
| 5 | 5.05 | 6.25 | 5.95 | 5.45 | +0.40 (+8 %) | 6.85 / 6.25 (+10 %) |
| 6 | 4.51 | 5.75 | 5.35 | 5.05 | +0.54 (+12 %) | 6.25 / 5.69 (+10 %) |
| 7 | 4.21 | 5.15 | 4.85 | 4.45 | +0.24 (+6 %) | 5.65 / 5.15 (+10 %) |

("exp width": the model 3-D P(v) Gaussian-broadened to the experimental
spread — a stand-in, not physics; experimental medians quantised at
0.25 Å/ps.) Roughly one third each: over-alignment, the 3-D width deficit
(model 5–12× too narrow in speed SD; a narrow shell projects to an edge
at its own speed whatever β is, so only the width moves the 2-D *peak*),
and a genuine speed offset that the corrected 3-D comparison shows as
well. "3-D fits, 2-D does not" came from comparing a narrow model peak
with the peak of a broad, skewed experiment, plus a 3–8 % hot reference.
For n = 1–2 the sign flips (model too slow in 3-D: the KE₁ deficit) and
2-D agrees with that.

**(c) Direction of the comparison.** Inverting the experiment (3-D) and
projecting the model (2-D) are both valid. 3-D tests speeds only and is
the score; 2-D tests speeds and angles together and is a diagnostic. 2-D
becomes a clean speed/width test only if the model is projected with the
measured β (post-processing re-orientation — not implemented).

### 13.6 Why is the experiment only weakly aligned? (OPEN)

The final β ≈ 0.3 cannot separate "weak from the start" from "lost on the
way": 15 % cos²-selected + 85 % isotropic gives β = 2·0.15 = 0.3, and so
does full cos² followed by strong randomisation. Candidates:

**A. Weak alignment already at ionisation (in the droplet).**
- A second, orientation-blind ionisation route next to direct
  strong-field ionisation — e.g. ionisation of the dopant via the helium
  (charge transfer to a He⁺ hole, Penning-type processes), which does not
  select the molecular axis. With axial recoil, β ≈ 2·f for a selected
  fraction f, so f ≈ 0.15.
- For: β is flat in n (an initial property is n-independent); the free
  molecule under the same laser is strongly aligned, so the selection law
  itself works — what changes is the environment.
- Against / open: an He-mediated channel should grow with intensity
  (clusters *less* aligned at high power); the crude raw-image ratio
  says I⁺He is slightly *more* aligned at 600 than at 160 mW
  (b/a 0.44 vs 0.30) — weak evidence, not inverted. The speed trend
  (§13.3) needs the two routes to give different speeds.

**B. Alignment lost in flight (deflection by the helium).**
- The physical picture behind the drag is many He collisions; the mean
  force has no transverse part, but each collision does.
- **Kinematic limit (ESTIMATE).** A heavy ion scattering off a light atom
  at rest is deflected by at most sin θ_max = m_He/m_ion ≈ 4/127 → θ_max
  ≈ 1.8° (less for I⁺Heₙ). With ≈ 20–40 collisions to stop (§7.4) the
  random-walk angle is ≈ 7–10°. Reducing β from 2 to 0.3 needs
  ⟨P₂⟩ × 0.15, i.e. an rms angular spread of ≈ 60° (small-angle
  diffusion: ⟨P₂⟩ factor ≈ exp(−1.5 ⟨δ²⟩)). Binary elastic deflection by
  single He atoms is short by almost an order of magnitude.
- So B needs a stronger momentum exchange than single elastic He
  collisions: collective / many-body He response, the ion being steered
  by the droplet's own asymmetry around an off-centre ionisation site, or
  transverse kicks in the early, dense phase (bubble formation). None of
  these is in the model, and none is quantified here.
- The speed trend argues against the simplest form of B: the slowest ions
  in a fragment are the most braked and should be the least aligned, but
  they are the most aligned.

**C. Two populations with different speeds (HINT; §13.3's trend survives
the inversion-artefact check, background structure not yet excluded).** A slower population that keeps the laser alignment and a
faster isotropic one. Fast and isotropic is what an orientation-blind or
higher-charge route would look like (the faster Q = 3 channel of §3 is
one candidate, but its ≤ 20 % share of n = 1 is too small to make all
fast ions isotropic). Speculative until the artefact check is done.

### 13.7 What the thesis can state

- The initial orientation (cos² along the polarisation) is correct in
  direction; the gas-phase images confirm a strong alignment for the free
  molecule.
- The model transports that alignment unchanged (β = 2); the experiment
  shows β ≈ 0.3 for every He-attached fragment and ≈ 0 for droplet-born
  bare I⁺. **The model over-aligns.**
- This does not affect any energy, speed or size observable (rotation
  invariance); the KE comparison is done in 3-D against the corrected
  reference.
- The 2-D comparison tests speed and angular distribution together. With
  the stochastic channel switched off the model cannot lose alignment;
  whether the strict alignment argument already fails at ionisation in
  the droplet or the alignment is lost in flight is not resolved here.
  An elementary kinematic estimate shows that single elastic He
  collisions cannot isotropise a heavy ion, so either the ionisation step
  or a stronger, collective momentum exchange must be responsible.
- Together with the too-narrow KE width at fixed n (§1.1) and the n ≥ 17
  tail (§7), this is the third observable pointing at the same missing
  physics: momentum exchange with the helium beyond the mean drag.

**Suggested text (draft).** "In the simulation the I⁺(He)ₙ fragments
keep the initial cos² alignment of the molecular axis (β = 2), because
the Coulomb recoil is axial and the mean drag force only decelerates the
ions. The experiment shows β ≈ 0.3 for all He-attached fragments. Since
the model is rotationally invariant, this does not affect the
kinetic-energy comparison, which we therefore perform on the
Abel-inverted (3-D) distributions; the detector projections (2-D) depend
on the angular distribution as well and serve as a diagnostic. Whether
the alignment is already weak at ionisation in the droplet or is lost
during the flight cannot be decided from the final distributions; single
elastic He collisions are kinematically too weak to isotropise a heavy
ion."

### 13.8 Open checks and outlook

1. **β(v) artefact test — DONE 2026-10-23:** pyabel does not create the
   trend (§13.3). Remaining check: background-subtraction structure
   (β(v) with the global-background scale varied, or bg-off).
2. **Power dependence, properly inverted:** β of the clusters at 160 /
   300 / 600 mW (corrected inversion). An He-mediated ionisation route
   (A) predicts falling β with power; in-flight loss (B) predicts none.
3. **Background check for n = 0:** β of droplet-born I⁺ vs the gas-phase
   scale factor (over/under-subtraction sensitivity).
4. **Model side (outlook, behind `SimConfig` enums, not implemented):**
   an initial-orientation option (cos² / channel mixture with fraction
   f) and a Tier-3 stochastic channel with a transverse component; and,
   for figures only, a post-processing re-orientation that projects the
   model with the measured β so that 2-D tests speed and width alone.

## 14. B2 thesis figures: caption drafts (added 2026-10-24; zero MD)

Figures: `scripts/post_processing/tier2_b2_thesis_figure.py` →
`data/figures/thesis_tier2/b2_{1,2,3,4}_*.png`, full text width
(16.6 cm), blue = experiment, red = B2. Every quoted number is checked
against the B2 row of `g5_winners_table_ihekedv2.csv` (reference v2,
G12.3) before drawing. Pipeline vocabulary that used to sit on the
figures (scored / detected read, median-anchor provenance, χ², smoothing
windows) is in these captions. Layout decisions (user, 2026-10-24):
sim / ref ratio replaces the profiled residual z in B2-2 (c) (the
profiled fit carries large uncertainties; the z values stay in G12.3);
B2-3 shows n = 1–6, so n = 0 (suppressed class vs bare I⁺, different
populations) is not in the results grid; B2-4 belongs to §13.

**B2-1 — size distribution.** "Terminal I⁺Heₙ size distribution for the
production law B2 (6000 molecules, 12 000 ions). Ions that stay bound
to the droplet (4.9 %) are excluded; 11 415 ions are scored. (a) Scored
distribution; the hatched n = 0 bar is the suppressed class (11.0 %).
Dashed line: mean size n̄ = 3.90; grey band: the n̄ landing band
[3.83, 4.31]. (b) Solvated branch (n ≥ 1, renormalised) compared with
the experimental abundance. n₁: share of n = 1 among solvated ions;
W₁: Wasserstein-1 distance between the two histograms in units of n;
n ≥ 15: solvated share in the high-n tail."

**B2-2 — mean kinetic energy.** "Mean kinetic energy per fragment size,
simulation (B2) vs the 3-D experimental reference (corrected inversion).
(a) ⟨E⟩ vs n on a logarithmic axis. Experiment: ⟨E⟩
with statistical ⊕ systematic errors; shaded: correlated calibration
(4 %) and condition (6 %) bands; open diamond: experimental median at
n = 1. Simulation: mean ± standard error of the detected ions; open
squares: sizes with fewer than two ions. (b) Ratio ⟨E⟩_sim / ⟨E⟩_exp per
size (error: simulation SE ⊕ experimental point error); line and band:
experiment with its correlated bands. Dashed: midHot = 1.050, geometric
mean of the ratio over n = 2–8; dotted: deepKE = 0.499, arithmetic mean
over n = 10–17."
(Optional sentence in the text, not on the figure: profiled
χ²_med = 349.7 over 15 sizes with ≥ 2 ions, n = 1 anchored to the
median.)

**B2-3 — speed distributions (3-D).** "Speed distributions of I⁺Heₙ,
n = 1–6. Experiment: Abel-inverted (3-D) VMI distribution, raw (dots)
and moving mean over 10 points (solid). Simulation (B2): |v| of the
detected ions, 0.04 Å/ps bins, moving mean over 15 bins (dashed); N =
ions per size. All curves are normalised to their maximum. Dotted
vertical lines: speed of the mean kinetic energy, v(⟨E⟩). The simulated
distributions are much narrower than the measured ones (§1.1, G7)."

**B2-4 — detector projection (2-D, §13).** As B2-3 with the 2-D
detector projection: experiment = projected VMI distribution (not
inverted), simulation = √(v_x² + v_y²). Caption must state that the
projection depends on the angular distribution: the model keeps β = 2,
the experiment shows β ≈ 0.3 (§13), so this comparison is a diagnostic,
not a score.
