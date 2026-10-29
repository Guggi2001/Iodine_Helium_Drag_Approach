# Thesis discussion points (final discussion round, 2026-10-02)

> **Status: COMPLETED (user, 2026-10-26).** The last model change
> (THESIS_FINAL_GOALS G13: self-bound gate lifted outside the droplet,
> `self_bound_gate_release = decoupled`, the thesis production setting) is
> folded in, and §15 holds the final observable table of B2 / B3 / h405
> against the v2 experimental reference. Any change after this date is a
> correction, not new material. Exception (user request, 2026-10-27): §16
> (shed timing and exit vs detection, zero MD on the §15 runs), which
> supersedes the §10.4 timing table and executes §12.

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
`THESIS_FINAL_GOALS.md` G12) is zero-MD. §14 (improved landing over B2:
the c lever on the v2 reference, the per-n shape of every measured lever,
untested candidates for a local n = 5–8 correction, 2026-10-25) is
zero-MD. §15 (final observable table, 2026-10-26) reads the G13
`decoupled` pools (`data/runs/{B2,B3,h405}_N6000_pooled`).

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
   but fails the 9 Å bound (−1.15 Å/ps re-landed: v_c 5.40, 3.5 Å
   gate, m = 21 He; −1.21 as fitted: v_c 5.5, 14.2 Å, m = 19). Neither TDDFT nor the detected
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
`data/figures/thesis_tier2/g6_chord_map_b2_only.jpg` and
`g6_chord_map_with_h405.jpg` (not committed): (a) scored n
vs chord (formerly suppressed ions boiled off → n = 0, G13; the panels read the `decoupled` re-runs), (b) KE vs chord, B2 per ion (N = 6000 pool, 12 000 fragments),
(c) chord distribution; the `with_h405` variant overlays the **re-landed**
h405 (§11; `post_pickup`, 3.5 Å, v_c 5.40, τ 4.4, E₀ 0.424; per ion from its
s1–s5 + r6 members, B2 chord bins) in all three panels (changed 2026-10-25;
before, the old h405 at 14.2 Å was overlaid binned from
`detfix_postpickup_chord.csv`). The raw-sampler A5 overlay was dropped
2026-10-24 (its runs are no longer on disk). Numbers: `g5_postpickup_chord.csv` (binned,
`lqppb2`; the script anchors on it to 4 dp), `h405_pp35_table.csv` (the
re-landed h405 pooled row; anchored to 4 dp), `g6_chord_stats.csv` (n = 1
window, KE₁ spread).

- The detected (n, KE) is a near-deterministic function of the
  **emission-axis chord**: the straight path from the birth point along the
  explosion axis (partner → fragment) to the droplet surface. The map does
  not depend on the size sampler (h405: `raw` vs `post_pickup` bin by bin
  within ≈ 0.01 eV and 0.1 in n; A5 likewise, G5 step A).
- Structure of the B2 map:
  - chord < 34 Å: the formerly **suppressed** class (11 % of scored
    ions) — the complex is ejected net-self-unbound (E_int > Σ(n)) with its
    full shell. Inside the droplet the self-bound gate holds the shell
    (confinement); once the complex is outside the droplet surface the
    confinement ends and the shell boils off (G13,
    `self_bound_gate_release = decoupled`): ε = 0, so the cascade reaches
    **bare I⁺** in ≲ 0.3 ns ≪ 8.5 µs (THESIS_FINAL_GOALS G13.6), co-moving
    sheds. The scorer already counted these ions at n = 0 (RQ8: excluded
    from the solvated read); re-running B2 with the release moves every
    score by ≤ 3 × 10⁻⁴ (`g13_decoupled_table.csv`). Panel (b) shows
    their bare-I⁺ KE;
  - 33–35 Å: the edge between the bare (boiled-off) class and n = 1
    (mixed bin: n = 1 share 0.42; mean KE 0.90 eV on the committed
    chord CSV, the gated read). Across it the speed is continuous and
    the KE step is only m(I⁺)/m(I⁺He) ≈ 0.97;
  - **n = 1 window 34.5–39.9 Å** (q05–q95 of the n = 1 ions, median
    37.2 Å, width 5.3 Å); KE₁ 0.609 eV;
  - beyond ≈ 40 Å: n rises again (2 → 13) — **re-pickup after
    stripping**, not partial evaporation;
  - > 80 Å: mostly retained (0.89).
- **Old vs new bundle (historical, 2026-10-13).** These numbers compare
  B2 with the *old* h405 at its own settings (`raw` sampler, 14.2 Å gate,
  gated `never` read), the overlay used before 2026-10-25. The figure now
  overlays the re-landed h405 instead, and its chord numbers (n = 1
  window, KE per bin) are not in `g6_chord_stats.csv`, which holds only
  B2. Do not quote these as the re-landed comparison. Old reading: the
  suppressed / n = 1 edge and the n = 1 window sit at the same chord
  (h405 window ≈ 35–39 Å) with nearly the same KE₁ (0.67 vs 0.63 eV in
  the 35–37 Å bin); the bundle change acts on the long-chord branch:
  under B2 KE falls more slowly with chord (60 Å: 0.26 vs 0.10 eV) and
  re-pickup is slower (mean n at 60 Å: 6.9 vs 10.3).
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
    ≈ 32–34 Å, at the bare / n = 1 edge where most ions leave
    net self-unbound and boil off to bare I⁺.
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

**Full-shell start (JUSTIFICATION; added 2026-10-27).** Ionization
occurs at the ground-state geometry, where HeDFT shows a well-defined
first solvation shell around the neutral (García-Alfonso et al., I₂
working notes). The single-atom shell counts, ≈20 He for I and 21 for
the relaxed I⁺, motivate starting each ion at the full I⁺ first shell,
n₀ = 21. How the molecular shell divides between the two fragments is
not computed, so n₀ = 21 is an upper estimate. Caveats:
- single pulse, R₀ = 2.666 Å: no neutral flight before ionization, so
  the ion inherits the equilibrium neutral shell (the TDDFT ion run
  instead ionizes after A-state flight to 9 Å, effectively pump–probe);
- the match is in count, not structure: the neutral shell is looser
  (first density peak ≈0.055 Å⁻³ for I vs ≈0.165 Å⁻³ for I⁺) and is
  compressed after ionization;
- at 2.666 Å the shell surrounds the molecule, not each atom; 21 per
  ion = 42 He at birth likely exceeds the uncounted I₂ first shell.
The extraction row above bounds the consequence: the drag fit is
insensitive to the 19 ↔ 21 fit mass (G11.4).

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

> **Superseded for the thesis by §16 (2026-10-27):** the timing table
> below is one seed under the old `never` release. §16 reads the pooled
> thesis runs (`decoupled`, s1–s6, scored ions) and gives 9.1 ps and
> 4.8 / 3.4 / 2.1 % at s_eff = 8. Kept as the s_eff comparison.

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

> **Executed 2026-10-27 in §16** on the thesis runs (ion.npz kept, no
> re-runs needed), with this definition unchanged. §12.3 (suppressed
> class) and §12.5 (cost) describe the state before G13 and are
> historical.

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
  no shed occurs after exit. (These reads are on the committed gated
  runs; under the G13 release, `self_bound_gate_release = decoupled`, the
  same ions boil off to bare I⁺ once outside the droplet.)

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

## 14. Improved landing over B2 (MEASURED re-scores + DISCUSSION; added 2026-10-25; zero MD)

Question: on the corrected KED reference (v2, G12) B2 runs ≈ 5 % hot
over n 2–8 (midHot 1.050), and the overshoot sits in n = 5–8 (per-n
ratio 1.15–1.21). Is there a lever that cools n = 5–8 and leaves the
other sizes less affected?

Numbers: zero-MD re-scores of existing run directories against v2 with
the G12 read path (`tier2_detfix_battery_table._row`, per-n
`ke_by_n(min_count = 1)`). Oracle: the c 1.6 / c 1.8 s1 runs reproduce
their committed v1 rows in `g5_postpickup_table.csv` to 4 dp, and the
B2 pool reproduces its committed v2 row in
`g5_winners_table_ihekedv2.csv`. Scratch scripts, **not yet committed
tables**; single-seed rows are N = 1000, pooled rows N = 6000.

### 14.1 The c lever at re-placed E₀ (B1c pair, s1, v2)

The G5 B1c step ran c 1.6 (E₀ 0.3297) and c 1.8 (E₀ 0.3331) on the same
seed (identical births; paired), E₀ placed per c from the twin line to
hold n̄. The selection of c 1.6 was made against v1.

| | c 1.6 s1 | c 1.8 s1 | B2 pool |
|---|---|---|---|
| n̄ / n₁ / W₁ (reference-independent) | 3.94 / 0.214 / 0.930 | 3.90 / 0.217 / 0.926 | 3.90 / 0.211 / 0.967 |
| suppressed / trapped | 0.114 / 0.050 | 0.122 / 0.060 | 0.111 / 0.049 |
| KE₁ mean [eV] (exp 1.198) | 0.609 | 0.577 | 0.609 |
| n 2–4 (geometric) | 0.918 | 0.861 | 0.916 |
| **n 5–8 (geometric)** | **1.146** | **1.061** | 1.164 |
| midHot n 2–8, v1 → v2 | 0.994 → 1.042 | 0.925 → 0.970 | 1.002 → 1.050 |
| n 9–12 (geometric) | 0.795 | 0.714 | 0.824 |
| deepKE n 10–17 (arithmetic, v2) | 0.593 | 0.597 | 0.499 |
| χ²_med v1 → v2 | 292 → 186 | 249 → 154 | 527 → 350 |
| S (W₁ + midHot) v1 → v2 | 1.67 → 1.93 | 2.18 → 1.84 | 1.71 → 2.04 |

- The sizes do not discriminate (E₀ re-placed per c absorbs them).
- On v1, c 1.6 won on midHot (0.994 vs 0.925); χ²_med already favoured
  c 1.8. On v2 the two cells sit about symmetrically around midHot = 1
  (1.042 / 0.970) and both χ² and S favour c 1.8. Interpolated (not
  run): c ≈ 1.7 gives midHot ≈ 1.0 on v2. The v2 correction moves the
  experiment's pick along the TDDFT ridge from c 1.6 to roughly
  c 1.7–1.8 — a strengthening of the ridge-freedom argument (§0), not
  a change of the production law (G12.2: no re-selection).
- c is a **global KE scale**, not a shape lever: c 1.6 → 1.8 multiplies
  ⟨E⟩ by 0.91–0.95 at every n = 1–10 (0.87–0.88 at n 11–12). It removes
  the n = 5–8 overshoot (1.146 → 1.061) but cools the sizes that are
  already too cold by the same factor (n 2–4 0.918 → 0.861, n 9–12
  0.795 → 0.714) and widens the KE₁ deficit (−0.032 eV, −5 %).

### 14.2 The shape of every measured lever (v2, per-n variant / baseline)

| lever (pair) | n 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | pattern |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| c 1.6 → 1.8, E₀ re-placed (s1) | .95 | .94 | .94 | .93 | .93 | .93 | .93 | .91 | .95 | .90 | .87 | .88 | uniform −7 % |
| E₀ 0.3155 → 0.3355 (c 2.0, s1) | .89 | .89 | .89 | .87 | .88 | .89 | .87 | .84 | .92 | .86 | .74 | .94 | uniform −12 %; n̄ 4.56 → 3.91 |
| τ 6.4 → 6.8 ps (c 2.0, s1) | .95 | .94 | .92 | .91 | .89 | .89 | .87 | .80 | .88 | .80 | .71 | 1.00 | uniform, slight tilt; n̄ 4.56 → 4.05 |
| s_eff 8 → 5, n̄ matched (pooled) | 1.03 | 1.02 | 1.02 | 1.02 | 1.00 | .98 | .97 | .94 | .91 | .88 | .89 | .92 | tilt, pivot ≈ n 5 |
| B2 → B3, E_bind 0.0482 + a 19.8 (pooled) | .97 | .98 | .98 | .99 | 1.00 | 1.02 | 1.04 | 1.09 | 1.12 | 1.17 | 1.25 | 1.28 | opposite tilt |

(The E₀ and τ rows are at the A5 cell under `post_pickup`, baseline
`g5lqpps1`; n̄ is not matched there.)

- What the misfit needs: the B2 ratio is 0.51 / 0.78 / 0.94 / 1.05 at
  n = 1–4, 1.15–1.21 at n = 5–8 and 1.03 → 0.63 at n = 9–12 — a hump.
  The correction is a **U** (heat n ≤ 4, cool n 5–8 by ≈ 13 %, heat
  n ≥ 9), the same diagnosis as `TIER2_MASS_SCENARIOS.md` §5 ("the
  required correction is a U, not a boost").
- c, E₀ and τ scale the hump; s_eff and the shallow well tilt it. No
  measured lever makes a local dip. s_eff 5 comes closest: n ≤ 4 +2–3 %
  (KE₁ +0.019 eV), n 7–8 −3 to −6 %, but n ≥ 9 cools most (−9 to −12 %,
  sizes already too cold); W₁ +0.10; χ²_med 350 → 241 (the scorer
  weights the well-populated small sizes). Combinations multiply tilts
  (e.g. s_eff 5 × B3 well mainly repairs the deep tail, not n 5–8;
  c ↑ × B3 well cools n ≤ 6 uniformly) and the shallow well is not
  Tier-0-consistent (G10.5).
- Why: the detected n is set by the emission chord (§1; n 5–8 ⇔ chord
  ≈ 45–55 Å), and every measured map lever acts smoothly along the whole
  chord axis — as a level (drag strength, evaporation) or as a tilt
  (s_eff, well + drag). A local correction needs structure at
  intermediate chords or sizes.

**Statement for the thesis.** The remaining KE-per-n misfit is a shape —
too cold at small and large n, too hot at n = 5–8 — and none of the
measured levers (c on the TDDFT ridge, E₀, τ, s_eff, well depth) changes
that shape locally; they scale or tilt it. The KE misfit is therefore not
a drag-calibration problem. On the corrected reference the experiment
would place c at ≈ 1.7–1.8 on the ridge (interpolated) instead of 1.6;
the identifiability argument is unchanged.

### 14.3 Untested candidates for a local correction (future investigation; no runs)

1. **Binding ladder D₀(n) at rungs 5–8.** The ladder sets at which chord
   an ion ends at a given n. Weaker binding at n ≈ 5–8 would let those
   ions shed one or two more He, so the n 5–8 bins fill from longer,
   slower chords and come out cooler. Unmeasured; it also moves the size
   histogram, and the neighbouring bins receive faster ions. The current
   ladder (rq4graded) only modifies rungs 1–3; the ladder shape is one of
   the two genuinely free knobs (Tier 2). First test: a ladder replay on
   frozen trajectories as in G9b (those are the c 2.0 `raw` trajectories,
   not B2).
2. **Flight mass / shell stripping** (`TIER2_MASS_SCENARIOS.md`). The
   existing explanation for a U-shaped correction; gated by M1 (a
   sourced I⁺–He D₀), open.
3. **Pickup rate λ₀ (+ Langmuir p).** Acts at shallow chords (small n);
   could lift n = 1–4 and ease the hump from below. Never swept at lq
   (D0 §7, §15.9).
4. **Combined c ≈ 1.7 pool on v2** (interpolation check, not a shape
   fix). Would turn the c ≈ 1.7 projection into a measured candidate:
   6 × N = 1000 MD runs with E₀ from the twin-line rule. Only if the
   G12.2 no-re-selection decision is revisited.

## 15. Final observable table: B2, B3, h405 vs experiment (MEASURED; added 2026-10-26; the thesis numbers)

The final state of all three laws, after the last model change
(THESIS_FINAL_GOALS G13: self-bound gate lifted outside the droplet,
`self_bound_gate_release = decoupled`). Every simulated value is a pooled
N = 6000 run (members s1–s6), `post_pickup` sampler, 3.5 Å density gate,
Stage I to 500 ps, partner-aware Coulomb closure:

| law | role | drag law | E_bind [eV] | E₀ [eV] | τ [ps] | pool |
|---|---|---|---|---|---|---|
| B2 | thesis production law | `linear_quadratic` a 18.007, c 1.6 | 0.1173 | 0.3297 | 6.4 | `data/runs/B2_N6000_pooled` |
| B3 | effective high-tail variant (not Tier-0-consistent) | `linear_quadratic` a 19.8, c 1.6 | 0.0482 | 0.3389 | 6.4 | `data/runs/B3_N6000_pooled` |
| h405 | re-landed contrast law | `capped_cubic` b 2.515, v_c 5.40 | 0.1168 | 0.424 | 4.4 | `data/runs/h405_N6000_pooled` |

Experiment = the v2 references the scorer reads: the size distribution
`integrated_i_he_abundance.csv` and the kinetic energies `IHe_KED_reference.csv`
(v2, 2026-10-23, G12). Sources: §15.1 = `g13_decoupled_table.csv` (the
decoupled rows, v2 scorer `tier2_g12_ihekedv2_table._score`), except the
rows marked † and the experimental KE medians / modes. Those, and all of
§15.2, come from `scripts/post_processing/tier2_thesis_observable_table.py`
→ `data/runs/h2b_forward_model/thesis_observables_extra.csv` and
`thesis_observables_per_n.csv`. That script reads the same pooled
`detection.npz` with the scorer's reader (`read_confirmation_detection`,
`solvated_renormalized`), and its oracle requires f(1) and KE(1) of every
law to reproduce `g13_decoupled_table.csv` to 4 dp. "—" means there is no
like-for-like experimental value.

### 15.1 Scalar observables

| group | observable | definition | experiment (v2) | B2 | B3 | h405 |
|---|---|---|---|---|---|---|
| sample | scored ions | ions reaching the detector, of 12 000 | — | 11 415 | 11 668 | 10 462 |
| sample | trapped fraction | droplet-retained (bound + marginal), not scored | — | 0.0488 | 0.0277 | 0.1282 |
| sample | — bound | conservatively bound (physics) | — | 0.0488 | 0.0277 | 0.1243 |
| sample | — marginal | still-coupled unbound (modelling exclusion) | — | 0 | 0 | 0.0039 |
| sample | suppressed share | intact net self-unbound complexes | — | 0 | 0 | 0 |
| sample | bare I⁺ share † | n = 0 among scored ions | 0.4352 ‡ | 0.1105 | 0.1514 | 0.2204 |
| size | n̄ (scored) | mean n over scored ions incl. bare | — (pre-registered band 3.77–4.37) | 3.9028 | 3.9559 | 3.9404 |
| size | n̄ solvated † | mean n over n ≥ 1 | 4.8886 | 4.3875 | 4.6615 | 5.0544 |
| size | n₁ solvated | n = 1 share of n ≥ 1 | 0.3103 | 0.2110 | 0.2111 | 0.1888 |
| size | W₁ solvated | Wasserstein-1 distance to experiment, n ≥ 1 | 0 | 0.9672 | 0.7952 | 0.9355 |
| size | n ≥ 15 mass | share of n ≥ 15 among n ≥ 1 | 0.0569 | 0.0008 | 0.0075 | 0.0040 |
| size | n_max | largest scored n | ≥ 20 (reference ends at 20) | 17 | 18 | 16 |
| KE | KE₁ mean [eV] | mean KE at n = 1 | 1.1975 | 0.6088 | 0.5898 | 0.6356 |
| KE | KE₁ median [eV] | | 1.0527 | 0.6086 | 0.5898 | 0.6367 |
| KE | KE₁ mode [eV] § | | 0.8516 | 0.5688 | 0.5688 | 0.6563 |
| KE | KE₂ mean [eV] | mean KE at n = 2 | 0.6687 | 0.5212 | 0.5085 | 0.5541 |
| KE | KE₂ median [eV] | | 0.6202 | 0.5212 | 0.5085 | 0.5538 |
| KE | KE₂ mode [eV] § | | 0.5701 | 0.4813 | 0.4813 | 0.5688 |
| KE | midHot | mid-n KE ratio sim / experiment (bins) | 1 | 1.0502 (7) | 1.0652 (7) | 1.0562 (7) |
| KE | deepKE | deep-bin KE ratio sim / experiment (bins) | 1 | 0.4986 (8) | 0.6276 (8) | 0.5603 (7) |
| KE | KE n = 10–14 / ref | mean of per-n KE ratios, bins with > 3 ions | 1 | 0.6705 | 0.8063 | 0.5638 |
| fit | χ²_med | profiled KE-curve χ² vs experiment (points) | — | 349.69 (15) | 365.98 (17) | 422.25 (16) |
| fit | S | joint score on W₁, midHot, deepKE (provisional norms) | 0 | 2.0440 | 1.8445 | 2.0297 |
| fit | MD gate | n₁ in 0.19–0.30 and n̄ in 3.77–4.37 | — | pass | pass | fail (n₁ 0.1888) |

† From `thesis_observables_extra.csv` (see above); not a column of
`g13_decoupled_table.csv`.
‡ Not like-for-like: the experimental n = 0 bin carries the
channel-external bare peak (RQ8), which the model does not produce. That
is why every size score uses the solvated (n ≥ 1) distribution.
§ The modes differ in method: model = midpoint of the tallest of 40 bins
on 0–3.5 eV; experiment = the reference's own `mode_KE_eV`. They are a
display diagnostic, not a scored quantity.

### 15.2 Per-n profile (solvated size distribution and mean KE)

f(n) = share of n among n ≥ 1. KE(n) = mean kinetic energy at that n [eV],
with the number of scored ions in brackets where it is below 30; "—" = fewer
than 2 ions (the scorer's minimum bin count). The experimental KE ends at
n = 17.

| n | f exp | f B2 | f B3 | f h405 | KE exp | KE B2 | KE B3 | KE h405 |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.3103 | 0.2110 | 0.2111 | 0.1888 | 1.1975 | 0.6088 | 0.5898 | 0.6356 |
| 2 | 0.1422 | 0.1623 | 0.1545 | 0.1420 | 0.6687 | 0.5212 | 0.5085 | 0.5541 |
| 3 | 0.0918 | 0.1345 | 0.1277 | 0.1157 | 0.4811 | 0.4512 | 0.4444 | 0.4825 |
| 4 | 0.0705 | 0.0892 | 0.0868 | 0.0837 | 0.3778 | 0.3971 | 0.3944 | 0.4206 |
| 5 | 0.0576 | 0.0850 | 0.0806 | 0.0760 | 0.3058 | 0.3505 | 0.3505 | 0.3635 |
| 6 | 0.0494 | 0.0730 | 0.0695 | 0.0720 | 0.2522 | 0.3041 | 0.3105 | 0.3027 |
| 7 | 0.0388 | 0.0603 | 0.0567 | 0.0671 | 0.2142 | 0.2529 | 0.2629 | 0.2350 |
| 8 | 0.0362 | 0.0511 | 0.0483 | 0.0549 | 0.1823 | 0.2047 | 0.2232 | 0.1845 |
| 9 | 0.0299 | 0.0422 | 0.0438 | 0.0505 | 0.1587 | 0.1638 | 0.1830 | 0.1371 |
| 10 | 0.0295 | 0.0344 | 0.0351 | 0.0413 | 0.1345 | 0.1241 | 0.1453 | 0.1016 |
| 11 | 0.0262 | 0.0273 | 0.0295 | 0.0378 | 0.1174 | 0.0900 | 0.1126 | 0.0751 |
| 12 | 0.0222 | 0.0175 | 0.0238 | 0.0345 | 0.1033 | 0.0650 | 0.0835 | 0.0560 |
| 13 | 0.0197 | 0.0086 | 0.0163 | 0.0222 | 0.0916 | 0.0531 | 0.0590 | 0.0420 |
| 14 | 0.0187 | 0.0029 | 0.0090 | 0.0094 | 0.0793 | 0.0360 (29) | 0.0428 | 0.0335 |
| 15 | 0.0162 | 0.0006 | 0.0049 | 0.0034 | 0.0717 | 0.0175 (6) | 0.0302 | 0.0392 (28) |
| 16 | 0.0144 | 0.0001 | 0.0019 | 0.0006 | 0.0666 | — (1) | 0.0196 (19) | 0.0371 (5) |
| 17 | 0.0087 | 0.0001 | 0.0005 | 0 | 0.0623 | — (1) | 0.0171 (5) | — (0) |
| 18 | 0.0072 | 0 | 0.0001 | 0 | — | — (0) | — (1) | — (0) |
| 19 | 0.0058 | 0 | 0 | 0 | — | — (0) | — (0) | — (0) |
| 20 | 0.0047 | 0 | 0 | 0 | — | — (0) | — (0) | — (0) |

### 15.3 Reading

- **Size.** All three laws hit the scored mean n̄ ≈ 3.9 by construction:
  E₀ was placed on it. Against the experiment they share the same shape
  misfit. There are too few n = 1 (0.19–0.21 vs 0.31). The bulk at
  n = 3–10 is over-populated in all three. The tail at n ≥ 14 is depleted
  (n ≥ 15 mass 0.001–0.008 vs 0.057). B3 is closest in W₁ and in the tail,
  the reason it is kept as the effective high-tail variant.
- **KE.** midHot (mid-n KE) lands within 7 % for all three laws. KE₁ and
  KE₂ are low: KE₁ is about half of the experimental value, which is the
  KE₁ deficit that `TIER2_MASS_SCENARIOS.md` attributes to the flight mass
  (gated by M1). The per-n curve is too flat. From n = 3–4 to n ≈ 9 the
  model is faster than the experiment; above n ≈ 10–11 it is slower, so
  deepKE and the n = 10–14 ratio lie well below 1.
- **G13 effect.** Against the committed `never` rows the release changed
  the suppressed share (to 0). It moved every other observable by ≤ 0.003
  and χ²_med by ≤ +0.16 (`g13_decoupled_table.csv`; THESIS_FINAL_GOALS
  G13.7). The misfits above are therefore not caused by the release.
- **Ranking.** B2 and B3 pass the MD gate; h405 fails it on n₁ (0.1888 <
  0.19). The ranking by S (B3 1.84 < h405 2.03 ≈ B2 2.04) and by χ²_med
  (B2 < B3 < h405) is unchanged from the `never` rows.

---

## 16. When and where the ion sheds He (MEASURED; added 2026-10-27; zero MD)

Sources: `scripts/post_processing/tier2_g13_shed_timing.py` →
`data/runs/h2b_forward_model/g13_shed_timing.csv`,
`g13_exit_vs_detection.csv`, `g13_exit_detection_hist.csv`,
`g13_shed_chord.csv`; figures `data/figures/thesis_tier2/g13_shed_timing.png`
and `g13_exit_vs_detection.png` (B2). Reads the thesis production members
(`{B2,B3,h405}_N1000_s1` … `_s6`, `decoupled`, ion.npz kept), i.e. the same
ions as §15, plus their Stage I event log `stage1_events.npz` (added
2026-10-27, see 16.1). Anchors (fail loudly): per member, the event log
replays the stored `n_shell` at every column, its exit falls in the
stored-depth exit interval, and n(t_h) − #detection events = n_det for
every scored ion; pooled scored count, n̄ and bare share reproduce §15 to
4 dp. This section **executes the §12 exit-vs-detection
readout** and **replaces the single-seed `never` shed-timing table of
§10.4** for the thesis.

### 16.1 Method

- **Shed events over the whole flight**, scored ions only (the §15
  sample). Both stages now have a complete event record:
  - *Stage I (0–500 ps): `stage1_events.npz`.* The ion checkpoint stores
    the state only every 16th step (dt 0.01 ps → 0.16 ps grid, 3126
    columns), so the 18 thesis members were re-run with an RNG-free event
    log (`simulation/stage1_events.py`, written by `run_ion_propagation(…,
    record_events=True)`; `scripts/regen_tier2_stage1_events.py`). Every
    re-run reproduced cfg.json and every array of neutral.npz, ion.npz and
    detection.npz bit for bit, so the §15 numbers are untouched. One row
    per shed and per pickup at the internal step: ion, step, time, kind,
    n before and after, and the position at which the event is applied.
    The biphasic channels fire at the pre-step seam of step k (time
    t_{k−1}, positions x(t_{k−1}), which the gate and pickup density read);
    a same-step shed + pickup gives two rows, shed first (the frozen draw
    order). Event times are therefore exact on the 0.01 ps grid.
  - *Stage II (500 ps – 8.53 µs).* One row per shed in the detection-stage
    event record, with a continuous (Gillespie) time.
  - *What the old frame reconstruction missed* (measured, `frac_stage1_*
    _hidden_by_frames`). Reading events from the net change of `n_shell`
    between stored columns cancels a pickup and a shed in the same
    0.16 ps interval. At B2 this hid **0.25 % of the Stage I sheds** and
    **19 % of the Stage I pickups** (B3 17 %, h405 20 %; the earlier
    estimate of ≈ 10 % was too low). Timing shares changed by ≤ 0.1
    percentage points; the frame reconstruction is kept only for this
    count.
- **Sign convention.** §16 and its figures use the signed radial
  coordinate **r − R** (positive outside; r from the droplet centre, R =
  the per-atom static radius `droplet_radii_angstrom`, the half-density
  point of the erf profile). It is written r − R, not d, because the
  thesis defines the depth as d = R − r (positive inside). (§12.2 and the
  code's `_depth` use the opposite sign under the name "depth"; the
  CSV column names `shed_depth_A`, `frac_sheds_depth_le_0`, … keep the
  code sign, i.e. they are r − R.)
- **Exit** as fixed in §12.2: last outward crossing of r − R = 0
  (r − R = 0 counts as inside), taken from the log's crossing rows (every
  sign change of r − R at the internal step, time linearly interpolated
  inside the step). Every scored ion of all three laws crosses **exactly
  once**, now confirmed at the 0.01 ps step (no re-entry hidden between
  stored columns). **n_exit is the occupancy right after the exit step.**
  The log's exit lies inside the 0.16 ps interval of the old stored-column
  exit for every scored ion (anchor); the old t_exit was 0–0.16 ps late
  (median 0.08 ps).
- **Booking before / after exit.** A Stage I shed is *before exit* if its
  step is ≤ the exit step: it fires at the pre-step seam, before the
  motion across the surface. Then n_exit − n_det = sheds after exit −
  pickups after exit exactly. The geometric booking (r − R ≤ 0 where the
  shed is applied) gives **the same number** for every class and law,
  because each ion crosses once. The former bracket 15.5–16.7 % is
  replaced by **16.3 %** (B2).
- **Position of an event.** *Stage I:* r − R at the position where the
  event is applied (the pre-step position, exact). *Stage II:* the
  event record has no positions and positions are not advanced; no
  handover position is reused and no ballistic position reconstructed.
  The events are booked as outside every threshold (r − R = +∞).
  Justification (checked, B2, all 3859 scored ions with a Stage II
  event): at handover these ions sit ≥ 230 Å outside (1st percentile
  802 Å, median 3557 Å) with outward radial velocity (min +0.56 Å/ps), so
  any straight-line continuation stays beyond every threshold used (0,
  3.5, 10.5, 150 Å). Figure `g13_exit_vs_detection` (c) draws only the
  finite Stage I values (cut at 150 Å) but normalises the bars to all
  sheds, so the Stage II 2.0 % is not drawn (see the caption, §16.6).
- **Classes** by detected n: bare (n = 0), n = 1, n ≥ 2. "Formerly
  suppressed" = the ions that the committed `never` members of the same
  seeds flagged `suppressed` (B2, h405).

### 16.2 Results at B2 (pooled N = 6000, 11 415 scored ions)

| class | share | n at exit | n at detection | exit time (median) | sheds per ion | before exit | median shed | sheds after 50 / 100 / 400 ps | ions with a shed after 400 ps | last shed (median) |
|---|---|---|---|---|---|---|---|---|---|---|
| **all** | 1 | **18.4** | **3.90** | **4.4 ps** | 17.3 | **16.3 %** | **9.0 ps** | **4.8 / 3.3 / 2.1 %** | **36 %** | 143 ps |
| bare | 0.110 | 21.0 | 0 | 2.7 ps | 21.0 | 0 % | 9.1 ps | 0.3 / 0 / 0 % | 0 % | 24 ps |
| n = 1 | 0.188 | 20.4 | 1 | 3.3 ps | 20.0 | 2.8 % | 9.4 ps | 3.9 / 2.7 / 1.8 % | 36 % | 123 ps |
| n ≥ 2 | 0.702 | 17.4 | 5.29 | 5.1 ps | 16.0 | 24.1 % | 8.8 ps | 6.0 / 4.3 / 2.6 % | 42 % | 228 ps |

From the Stage I event log (§16.1): exact event times and positions,
interpolated exit. The geometric booking (`frac_sheds_depth_le_0`)
equals "before exit" in every class.

Further numbers (all scored ions):

- Shed-time quantiles: q10 3.9 ps, q90 24 ps, q99 5.3 ns. 2.0 % of the
  sheds fall in the detection stage (their median 5.8 ns).
- Per ion the slow tail is long: 27 / 15 / 7.3 / 2.9 % of the ions still
  shed after 1 ns / 10 ns / 100 ns / 1 µs.
- Exit: 10–90 % at 3.0–7.7 ps, born at r − R = −35 Å (median; thesis depth d = 35 Å).
  Outside sheds come a median 6.2 ps after exit.
- Position at shed: 79.7 % of all sheds lie more than 3.5 Å (one erf gate
  width) outside the sharp radius, 72.2 % more than 10.5 Å.
- n changes after exit for all but 3 of the 11 415 ions. Per ion, 14.5
  sheds and 0.023 pickups come after exit, against 2.8 sheds and 0.20
  pickups before it (mean n_exit − n_det = 14.5; n_exit ≥ 9 for every
  ion).
- Chord (`g13_shed_chord.csv`): the in-droplet share of the sheds rises
  with the emission-axis chord, 0 % (< 30 Å), 1.8 % (30–40), 12 %
  (40–50), 25 % (50–60), 45 % (60–70), 65 % (> 70 Å); the exit time rises
  from 2.5 to 9.4 ps along it. Bare ions come only from chords < 40 Å.
- **Formerly suppressed = bare.** 1260 of the 1261 B2 bare ions are the
  ions the `never` run flagged `suppressed` (the remaining formerly
  suppressed ion ends at n = 1, as logged for G13). They leave at 2.7 ps
  with the full shell (n_exit 21.0, no shed inside) and boil off
  completely outside (last shed: median 24 ps); none sheds after 100 ps.

Contrast laws (all scored ions): B3: median 9.1 ps, 4.6 / 3.2 / 2.0 %
after 50 / 100 / 400 ps, 34 % of ions after 400 ps, 16.3 % before exit,
n_exit 18.4. h405: median 9.1 ps, 4.3 / 3.1 / 1.9 %, 32 %, 11.5 % before
exit, n_exit 19.1 (2301 of its 2306 bare ions formerly suppressed). The
timing and the place of shedding are a property of the RRK cascade, not
of the drag law.

Against §10.4 (s1 only, `never`, all ions incl. suppressed, frame
reconstruction): 9.0 ps and 5.2 / 3.7 / 2.3 % there, 9.0 ps and
4.8 / 3.3 / 2.1 % here. The G13 release and the pooling change the
timing by less than its single-seed spread.

### 16.3 Reading

1. **Fast, as in the literature, per event.** The median shed at 9 ps,
   90 % by 24 ps and 97.9 % of all sheds before 400 ps match the classical
   MD of Albrechtsen et al. 2023 ("almost exclusively during the first
   tens of picoseconds"; "the vast majority of dissociation events are
   finished after 400 ps"; Na⁺He_N, a simulation, not an experimental
   timescale, §10.2).
2. **Per ion, the near-threshold tail is not small.** 36 % of the scored
   ions shed once more after 400 ps, 2.9 % still after 1 µs. These are
   single last sheds just above threshold (§10.3), the counterpart of the
   paper's "slow decrease … close to the binding energy". For these ions
   the detected n depends on the flight time, which is why a too-large
   s_eff shows in n₁ (§10.4).
3. **In the model the shell boils off outside the droplet.** The ion
   leaves within ~3–8 ps, faster than the cascade (median 9 ps), carrying
   18 of its 21 He on average. 83.7 % of all sheds happen after the exit
   and 80 % more than one gate width outside; the internal energy loaded in the
   droplet drives the cascade in vacuum. Only the long-chord ions (> 60 Å),
   which stay inside for 7–10 ps, lose a large part of their shell inside.
4. **The bare class is a post-exit product.** The bare ions are the
   short-chord, formerly suppressed ions; they exit first (2.7 ps) with
   the full shell and lose it within tens of ps outside. No bare ion is
   made inside the droplet.
5. **Consequence for the readout (§12):** the detected size distribution
   is not a damped copy of the size distribution at exit. The ions that
   leave with the largest shells (n_exit 20–21, short chords, earliest
   exit) end smallest (bare / n = 1; figure `g13_exit_vs_detection`,
   panel b). A size comparison with the experiment is a comparison of the
   end of the vacuum cascade, i.e. of E_int at exit and of the D₀ ladder.

**Caveats.** (i) That most shedding happens outside follows from the
model's construction: RRK evaporation is not density-gated, while
pickup, drag and cooling are (§12.1), so an ion that leaves hot keeps
shedding and stops cooling. No observable here tests whether the real
complex sheds inside or outside. (ii) "Outside" refers to the sharp
radius (half-density point); the r − R numbers show it is not an
edge-of-profile artefact (72 % beyond three gate widths). (iii) Stage I
events, positions and the exit come from the complete event log at the
0.01 ps step; Stage II depths are bounded, not measured; see §16.1.
(iv) The literature comparison carries the transfer caveats of §10.2
(Na⁺ vs I⁺, classical MD, complexes well above threshold).

### 16.4 Text for the discussion chapter (suggested)

> In the model the He shell is not lost in the droplet: the ion leaves the
> droplet within 3–8 ps, faster than its evaporation cascade, and sheds
> about 84 % of its He in vacuum, driven by the internal energy it picked up in
> the droplet. Dissociation is fast per event (median 9 ps, 98 % of all
> sheds before 400 ps), consistent with the classical MD estimate of
> Albrechtsen et al. 2023, but about a third of the ions shed a last He
> near threshold after 400 ps, some after microseconds. The detected size
> therefore reflects the end of the vacuum cascade, not the state at
> droplet exit, and the bare I⁺ ions are a product of this cascade: they
> are the ions that leave first, with the full shell, and lose it outside.

### 16.5 Kinetic energy at exit, handover, after the closure and at detection (MEASURED; added 2026-10-27)

Source: same script → `g13_ke_readouts.csv`. KE at four points per ion:
exit (stored state at the first stored column after the exit, 0–0.16 ps
late: the log gives the exit time, not the state; `n_exit_col_mean` is n
there), handover t_h (last stored column,
= just before the partner-aware Coulomb closure, which acts before the
Stage II event loop), post-closure (½ m_h |v_det|²: the closure changes v
at fixed m, and the `co_moving` Stage II sheds keep v), detection. Anchors
(fail loudly, all three laws, scored ions, ≤ 10⁻¹² eV; measured ≤ 3·10⁻¹⁶):
stored KE = ½mv² at exit, handover and detection; closure ΔKE =
−(E_pot,det − E_pot,h − Σ fold); Stage II ΔKE = −Σ `event_dE_mass_transfer`;
pooled KE₁ mean = §15 to 4 dp. n̄ at handover is read directly:
**n̄_h = 4.241** (B2; = n̄_det + 0.338 Stage II sheds per ion).

| B2 class | KE exit | KE handover | KE after closure | KE detection [eV] | ln-split: mass | ln-split: speed | v_det / v_exit |
|---|---|---|---|---|---|---|---|
| all | 0.633 | 0.439 | 0.440 | 0.437 | −0.342 | −0.025 | 0.992 |
| bare | 1.108 | 0.713 | 0.714 | 0.714 | −0.505 | +0.066 | 1.034 |
| n = 1 | 0.916 | 0.613 | 0.615 | 0.609 | −0.464 | +0.056 | 1.028 |
| n ≥ 2 | 0.482 | 0.349 | 0.351 | 0.347 | −0.284 | −0.061 | 0.975 |

The split is exact: ln(KE_det/KE_exit) = ln(m_det/m_exit) + 2 ln(v_det/v_exit),
order-independent. The speed term lumps all forces after exit (pair
Coulomb, well, drag tail, pickup recoil); it is **not** split by force.

Reading: (1) the KE drop after exit is mass loss, ≈ 93 % of
ln(KE_det/KE_exit) for all ions: co-moving sheds carry their share of the
KE away; (2) the closure adds only ≈ +0.001 eV, so by 500 ps the pair
Coulomb work is practically done; (3) the speed change after exit is a
few per cent and its sign depends on the class and the law (B2 n ≥ 2
−2.5 %, light classes +3 %; B3 positive in every class, h405 n ≥ 2
−4.6 %). Without a force split this sign is not attributed to a channel.

### 16.6 Figure captions (suggested)

**`g13_shed_timing`.** He evaporation events of the scored ions at the
production law B2 (pooled N = 6000, 11 415 ions). (a) Event times over
the whole flight on logarithmic bins (59 bins from 0.1 ps to 9 µs; ≈ 0.13
decade each). The height is a density per decade, dF/dlog₁₀ t: the
fraction of all events in the bin divided by the bin width in decades, so
the area over log₁₀ t is 1 (and heights can exceed 1). The stack is split
by readout, not by position: Stage I events before the droplet exit
(16.3 %), Stage I events after the exit (81.8 %), and Stage II events
(500 ps – 8.53 µs, 2.0 %). Both stages are complete event records: Stage
I times on the 0.01 ps MD step, Stage II times continuous.
(b) Cumulative fraction of events versus time, for all scored ions and by
detected occupancy n_det = 0, 1, ≥ 2 (median 9.0 ps; 2.1 % of the events
after 400 ps). (c) Cumulative fraction of ions versus the time of their
last recorded evaporation event: 36.3 % of the ions have their last event
after 400 ps and 2.9 % after 1 µs. In (b) and (c), the grey band and the
dotted line mark the statements of Albrechtsen et al. 2023 (classical MD
of Na⁺He_N) that dissociation happens "within the first tens of
picoseconds" and that the vast majority of events is "finished after
400 ps". Dashed: handover from Stage I to Stage II (500 ps); dash-dotted:
detector arrival (8.53 µs).

**`g13_exit_vs_detection`.** Shell size at droplet exit versus at
detection, B2 (pooled N = 6000, 11 415 scored ions). (a) Distribution of
n at the droplet exit (light; right after the last outward crossing of
r = R, median 4.4 ps) and at detection (dark), both
as fractions of the scored ions; the readout gives the mean n at exit,
at the 500 ps handover and at detection. (b) Per-ion map of n at exit against n at
detection (log colour scale; dotted line n_det = n_exit). (c) Radial
position r − R of each shed (positive outside the droplet; R = sharp
radius at half density); solid line r − R = 0, dotted lines at one and
three density-gate widths (3.5, 10.5 Å). The histogram is normalised to
all sheds but drawn only up to r − R = 150 Å: the 10.8 % of sheds
beyond 150 Å are not shown (annotated in the panel). These include all detection-stage sheds
(2.0 %), which carry no stored position and lie beyond 230 Å at handover.
Ion-stage positions are those at which each shed is applied (MD event
log, 0.01 ps step).
