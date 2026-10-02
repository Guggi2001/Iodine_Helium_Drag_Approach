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

---

## 1. The chord → (n, KE) map (MEASURED, D0 §15.8 reads 5–6)

- The detected (n, KE) is a near-deterministic function of the
  **emission-axis chord**: the straight path from the birth point along the
  explosion axis (partner → fragment) to the droplet surface. The map is the
  same under the `raw` and `post_pickup` size samplers (bin by bin within
  ≈ 0.01 eV and 0.1 in n).
- Structure of the map (`detfix_postpickup_chord.csv`):
  - chord < 33 Å: exits intact at n = 21 (supp class; matched on fraction
    only, per-ion identity with the scorer's suppressed class unchecked);
  - 35–39 Å: n = 1 (KE 0.61–0.67 eV);
  - beyond 39 Å: n rises again (2 → 14) — **re-pickup after stripping**,
    not partial evaporation; stripping itself is a step 21 → 1 within
    ≈ 33–37 Å;
  - > 65 Å: mostly retained.
  - KE falls ≈ 0.03–0.04 eV per Å of chord, consistent with the capped-tail
    force b·v_c³ ≈ 0.043 eV/Å.
- **Correction to a tempting reading:** the ion does not "choose the
  closest path to exit". The I–I axis at ionization fixes the direction and
  the partner's Coulomb push drives the fragment along it; the relevant
  length is the chord along that axis, not the radial distance to the
  surface. This is why radial birth depth is the weaker proxy.
- **Consequences.**
  - Geometry (droplet size, sampler, birth depth) only reweights the chord
    distribution; (v_c, τ, E₀, ladder) set the map. The histogram is the
    chord distribution pushed through one map.
  - **No geometry change can move KE₁.** The n = 1 window position on the
    chord axis fixes KE₁; only the map can move it: Coulomb budget, drag,
    or flight mass (`TIER2_MASS_SCENARIOS.md`, M1 gate).
  - The n = 1 KE width in the model (SD 0.038, the "needle") is the window
    width × the KE slope. The experimental width needs scatter in the map
    itself (Tier-3 noise or a channel mixture), not geometry.
  - Projection (not measured): on the current map KE ≈ 1.0 eV sits at a
    chord of ≈ 33–35 Å, where n ≈ 18. The 0.36 eV deficit is equivalent to
    stripping completing ≈ 3–5 Å of chord too late, or to KE falling too
    steeply along the chord.
- Caveats: the chord ignores deflection (partner Coulomb, drag) and the soft
  erf edge; the arms are not CRN, so the agreement is ensemble-level.

## 2. Size sampler: `raw` vs `post_pickup` (MEASURED, D0 §15.8 reads 1–4)

- Production (h405, N = 6000 pool) uses `raw`, the plain source log-normal
  at 40 bar / 14 K. Thesis Ch. 3 describes the pickup-conditioned ensemble
  (`post_pickup`), which legacy MATLAB production used.
- `post_pickup` − `raw` (pooled N = 5000 each): ⟨N⟩ +29 %, R̄ +4.5 Å,
  n̄ +0.72, n₁ −0.033, W₁ +0.24, trap +0.057, supp −0.059, KE₁ +0.007 eV.
  Both gate clauses fail under `post_pickup`.
- Discussion points:
  - The standing point is **conditional on the sampler**. `post_pickup` is
    not ruled out; it would need its own (v_c, τ, E₀) arbitration
    (projection: ≈ +0.03 eV in E₀ to recover the n̄ shift, untested).
  - The thesis headline claims (KE₁ deficit, chord-map argument) are
    sampler-invariant.
  - Ch. 3 and production must agree: either `raw` stated as a modelling
    assumption with §15.8 as its measured sensitivity, or `post_pickup` at
    the planned re-arbitration. The user intends to redo the arbitration
    with further improvements; the sampler choice belongs to that redo.

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
  prediction (relevant to the planned re-arbitration).
