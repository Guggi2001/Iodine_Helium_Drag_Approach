# Tier 2 — Parameter Influence at the Thesis Production Law (B2)

> **What this is (2026-10-14).** The per-lever influence reference for the
> thesis: what each lever does to the detected observables, at which law
> it was measured, and where the evidence lives. It replaces the atlas D0,
> now frozen as `TIER2_PARAMETER_INFLUENCE_h405_finc1v725.md` ("D0 §x"
> everywhere refers to that archive). **No new analysis** — every number
> is a quote of a committed measurement.
>
> **Not here.** Parameter values and classes: `CALIBRATION_MAP.md`.
> Narrative (storyline, chord map, sampler, mass, Coulomb budget):
> `THESIS_DISCUSSION.md`. Run records: `drag_migration_log_tier2.md`,
> `THESIS_FINAL_GOALS.md`.

## 0. Reference point and how to read the table

**B2** (CALIBRATION_MAP production column): `post_pickup` sampler, lq
ridge c 1.6 / a 18.007 / E_bind 0.11727 eV, 3.5 Å density gate, τ 6.4 ps,
E₀ 0.3297 eV; pooled N = 6000. Observables at B2 (THESIS_FINAL_GOALS G5
Outcome): n̄ 3.903, n₁ 0.211 (exp 0.310), W₁ 0.967, midHot 1.002, KE₁
0.609 eV (exp 1.302), deepKE 0.477, n ≥ 15 0.0008 (exp 0.057), max n 17,
trapped 0.049, suppressed 0.111.

**Observables.** n̄ detected mean size; n₁ solvated n = 1 share; W₁
Wasserstein distance of the solvated size histogram to experiment; midHot
/ deepKE mean-KE ratio model / experiment over n 2–8 / n 10–17; KE₁ mean
n = 1 KE; trapped = droplet-retained fraction; suppressed = n = 0 class
(excluded from the solvated read, RQ8). All read at t_detect = 8.53 µs.

**Basis flag — read this before quoting any number.**

| flag | meaning | how far it carries to B2 |
|---|---|---|
| **B2** | measured on B2 itself (pool or its members) | directly |
| **lq** | measured on the same law family next to B2: lq ridge, 3.5 Å gate, τ 6.4, mostly the A5 cell (c 2.0, E₀ 0.3155) under `raw` or `post_pickup` | sign and order of magnitude; magnitudes shift with c and sampler |
| **T** | **transferred**: measured at h405 or finc1v725 (cubic / capped cubic, 14.2 Å gate, `raw`, other τ / E₀) | mechanism and sign only — **magnitude not established at B2** |
| **—** | never measured | stated gap |

## 1. The organising picture: the chord map (B2)

The detected (n, KE) is a near-deterministic function of the
emission-axis **chord** (birth point → droplet surface along the
explosion axis) (`THESIS_DISCUSSION.md` §1, §1.1; G6). At B2: suppressed
below ≈ 34 Å, n = 1 window 34.5–39.9 Å, re-pickup beyond ≈ 40 Å, mostly
retained beyond 80 Å; KE slope along the chord −0.021 eV/Å.

This splits the levers into two classes:

- **Geometry levers** (droplet size, sampler, birth depth) only
  **reweight the chord distribution**; the map is unchanged (B2 and A5
  `post_pickup` have bit-identical chords). They move the size histogram
  but cannot move KE at fixed n — in particular not KE₁.
- **Map levers** (drag c / a, well E_bind, E₀, τ, ladder, Coulomb budget,
  flight mass) change (n, KE) per chord. Only these can move KE₁.

## 2. Lever table

| lever | B2 value | influence (measured) | basis | source |
|---|---|---|---|---|
| **Map levers** | | | | |
| drag shape c on the TDDFT ridge (a, E_bind follow) | c 1.6 | lower c = faster mid-n ions: at matched n̄ ≈ 3.9 under `post_pickup` (E₀ re-placed per c) midHot 0.87 (c 2.0) → 0.925 (1.8) → 0.994 (1.6). The E₀ that holds n̄ falls with c (≈ 0.017 eV per amu/Å: c 1.8 0.3331, c 1.6 0.3297). On the long-chord branch B2 vs h405: KE at 60 Å 0.26 vs 0.10 eV, slower re-pickup (mean n 6.9 vs 10.3). Experiment selects c (raw → 2.0, `post_pickup` → 1.6) | lq / B2 | G5 (B1, B1c); D0 §15.9; `THESIS_DISCUSSION.md` §1 |
| drag law family (ridge vs capped cubic) | lq | a different shape lands comparably (h405: W₁ 0.755, KE₁ 0.640, midHot 0.948 vs B2 0.967 / 0.609 / 1.002); KE₁ ≈ half the experiment for every law | B2 + T | G5 Outcome, identifiability |
| well depth E_bind | 0.11727 eV (ridge-tied) | well-only cut to 0.0482 at A5: midHot 0.995 → 1.130, KE₁ 0.617 → 0.653, n ≥ 15 ×10, trapped 0.035 → 0.013, n̄ +0.16. With drag retuned to midHot ≈ 1 (B3 at B2: a 19.8): n ≥ 15 ×9, W₁ −0.17, KE₁ −0.019. **Leaves the TDDFT band** (18 Å RMSE 0.58 vs band 0.19; no a rescues 0.0482, best 0.226; a refit still passes at 0.080 eV, band edge between the two). Slope dn̄/dE_bind ≈ −7.3 /eV | lq / B2; slope T + lq | G10, G10.5; G5 (B3); D0 §9.2, §6.7 item 2 |
| E₀ = E_int(0) | 0.3297 eV | sharpest n̄ lever: dn̄/dE₀ ≈ −31 /eV. At A5 `post_pickup` E₀ 0.3355 / 0.3455 / 0.3555 → n̄ 3.91 / 3.60 / 3.34, midHot 0.87 / 0.82 / 0.78, KE₁ 0.55 / 0.52 / 0.49 (extra evaporation drops slow heavy clusters into low n). Cold minority (≤ 0.05 eV below E₀) makes no tail; the tail needs E₀ ≲ 0.18 eV | lq | D0 §15.9; G9 |
| τ (cooling clock) | 6.4 ps | at A5 `post_pickup` τ 6.8 / 7.2 / 7.6 → n̄ 4.05 / 3.59 / 3.20, midHot 0.89 / 0.79 / 0.70, KE₁ 0.59 / 0.55 / 0.52 (dn̄/dτ ≈ −1.1 /ps). Under `raw` at every c 1.8–2.5: low τ keeps the large-n tail (W₁ ↓) but costs KE₁; high τ recovers KE₁, empties n ≥ 10 | lq | D0 §15.9; G4.14 |
| ladder D₀(n), n* | rq4graded (flat 74.4 cm⁻¹ base, rungs 1–3 × 2.2 / 1.5 / 1.3: empirical low-n strengthening inspired by the abundance steps, not a histogram inversion), n* 21 | on frozen lq trajectories (c 2.0) at matched n̄: every ladder / n* variant keeps n ≥ 15 ≤ 0.005 — high n stays on slow ions. Ladder family (rq4graded vs floor1) owns the n₁ / small-n shape | lq (replay); family choice T | G9b; D0 §8 |
| Coulomb budget | 2.70 eV / fragment | KE₁ slope S_k = 0.389 eV/eV; a uniform re-anchor to the calibrated ≈ 2.2 eV breaks the landing (n̄ 5.33, trapped ×2); does not widen KE at fixed n | T (h405) | `THESIS_DISCUSSION.md` §9; D0 §19 |
| RRK s_eff | 8 | sets the shedding clock; moves the detected n once near-threshold sheds outlast the 8.53 µs flight (s ≳ 12): s 5 / 8 / 12 / 20 / 30 → n̄ 3.74 / 3.94 / 4.21 / 4.68 / 5.28, n₁ 0.233 / 0.214 / 0.183 / 0.121 / 0.030, KE₁ flat (0.606–0.650). At matched n̄ (s 5, E₀ 0.3240): χ² −16 %, KE₁ +0.019, W₁ +0.10, n ≥ 15 → 0.0002 — same KE ↔ tail trade as τ. Explored s 5–30: 5–8 consistent, ≥ 20 excluded; s < 5 and 8–12 untested (§2.1) | B2 | §2.1; `THESIS_DISCUSSION.md` §10 |
| flight mass / Tier-0 fit mass | n₀ 21, fit 19 → 21 He | fit mass 19 → 21: a −1.4 to −2.3 % per c, E_bind +0.0004 eV, B2 passes unchanged. The flight mass, not the drag, carries the KE₁ deficit (open, gated by M1) | B2 (Tier 0); M-question open | G11.4, G11.5; `TIER2_MASS_SCENARIOS.md` |
| initial shell (bare start) | n₀ = 21 for all | bare-born ions (n₀ 0, E_int 0) land at n = 2–4 by pickup; KE₁ 0.49 (lq) / 0.87 (h405) vs exp 1.13 median — no fast n = 1 source | lq (A5) + T | G8; `THESIS_DISCUSSION.md` §8 |
| shed momentum convention | `co_moving` | `cold` injects ×1.61 KE₁ over a full strip; histogram convention-blind | T | D0 §12 |
| **Geometry levers** | | | | |
| size sampler (`raw` → `post_pickup`) | `post_pickup`, ⟨N⟩ 16.7k | n̄ +0.70, n₁ −0.04; W₁, midHot, KE₁ flat at A5 (h405: n̄ +0.72, W₁ +0.24, trapped +0.057). Moves ions from short chords (n = 1 window) to long ones; map invariant | lq + T | D0 §15.8–15.9; `THESIS_DISCUSSION.md` §2 |
| droplet size R | R̄ ≈ 54 Å | acts only through the chord distribution (sampler row). R is a selection knob, birth depth the physics knob; pinning R̄ alone instead of the distribution worsens W₁ | B2 (chord) + T | `THESIS_DISCUSSION.md` §1; D0 §14 |
| birth depth / chord | mean depth ≈ 36 Å | sets n and KE per ion through the map (n = 1 ⇔ chord 34.5–39.9 Å). Radial depth is the weaker proxy; ≈ 1 He/Å depth slope at the finc1v725 geometry | B2 (map) + T (slope) | `THESIS_DISCUSSION.md` §1; D0 §14.1, §14.3 |
| density width s_ρ | 3.5 Å | at h405 3.5 vs 14.2 Å: W₁ +0.06, KE₁ +0.003; width scan null on KE₁ (+0.002 eV over 4.5×). Enters mainly through the Tier-0 refit; birth dressing saturates (n₀ = 21) | T | G4.14 (`h405g35`); D0 §9.6.5; G2 |
| **Mechanism levers not measured on the lq family** | | | | |
| pickup λ₀ (+ Langmuir p) | 0.9 /ps (prescribed transferred estimate within the GAH25 heavy-alkali range ≈ 0.5–1.1 /ps; CALIBRATION_MAP row 7), p 1 | re-filling ≈ 7 ions / 100 at shallow chords; p = 1 lifts n₁ KE toward experiment. Named as the untested pickup-side n̄ lever at lq | T | D0 §7; D0 §15.9 |
| RRK prefactor ν | 2.42 /ps | cascade rate prefactor (Sourced from the I⁺–He V''(R_e)); never swept. s_eff moved to the map levers (measured at B2, §2.1) | — | D0 §0, §16 |
| retention f_ret | 0.1 | never swept | — | D0 §7 |
| noise (Tier 3) | off | not modelled; KE width at fixed n ≈ 8× too narrow at B2 | B2 (width) | `THESIS_DISCUSSION.md` §1.1, §7.4 |

### 2.1 RRK s_eff at B2 (measured 2026-10-14; basis B2)

Full record: `THESIS_DISCUSSION.md` §10; tables
`data/runs/h2b_forward_model/b2_seff_table.csv` (scan, seed s1, N = 1000)
and `b2_seff_timing.csv` (shed timing); pools `…g5lqppc160e3297seff05pooled_th500`
and `…seff05e3240pooled_th500` (N = 6000). Supersedes the transferred
reading "detector-compressed except the n = 1 bin" (I20 / I32).

- **Mechanism.** k = ν(1 − D₀/E_int)^(s−1). s sets the clock: sheds per
  ion barely change (14.8 → 13.1 over s 5 → 30), median shed time
  7.8 → 32 ps. At finite flight time an ion stops shedding below
  E* ≈ D₀ / [1 − (k_min/ν)^(1/(s−1))] (k_min = 1/8.53 µs): 1.02 / 1.10 /
  1.28 / 1.70 / 2.27 D₀ at s 5 / 8 / 12 / 20 / 30. For s ≲ 8 the cascade
  reaches its floor (clock only); for s ≳ 12 the near-threshold tail is
  pending at the detector (`time_exhausted` 0.01 / 0.07 / 0.20 / 0.49 /
  0.82) and fewer pickups occur in the droplet → n̄ up, n₁ down.
- **Sensitivity (single seed).** n̄ 3.74 / 3.94 / 4.21 / 4.68 / 5.28;
  n₁ 0.233 / 0.214 / 0.183 / 0.121 / 0.030; χ² 250 / 292 / 406 / 612 / 998;
  W₁ 1.01 / 0.93 / 0.90 / 0.98 / 1.32; midHot 0.957 → 1.25; KE₁ 0.606 →
  0.650 (insensitive).
- **Timing vs Albrechtsen et al. 2023** — their **supplementary classical
  MD**, not experiment (Methods + Extended Data Fig. 4: isolated
  Na⁺He₆ / Na⁺He₁₄; dissociation within tens of ps, vast majority done by
  400 ps): sheds after 400 ps 1.0 / 2.3 / 4.4 / 9.1 / 16 %, after 50 ps
  2.2 / 5.2 / 10 / 23 / 40 %. B2 (s 8) consistent; s 30 fails, s 20
  borderline. A cross-system plausibility check only.
- **Explored range, not a bound.** Tested only s 5 / 8 / 12 / 20 / 30;
  s < 5 and 8–12 untested. Size distribution: ≥ 20 excluded, 12
  disfavoured (n₁ 0.183 vs gate edge 0.19, ≈ 0.7 single-seed sd; χ²
  +40 %), 5–8 consistent. 8 is a convention within the explored range.
- **s_eff 5 vs 8 at matched n̄ (pooled N = 6000, E₀ 0.3240 vs 0.3297).**
  n̄ 3.909 / 3.903, n₁ 0.210 / 0.211, midHot 0.994 / 1.002; χ² 442 / 527,
  KE₁ 0.628 / 0.609, deepKE 0.556 / 0.477; W₁ 1.066 / 0.967, n ≥ 15
  0.0002 / 0.0008. Not an E₀ shift: a second lever on the τ-like KE ↔
  large-n-tail trade-off. Production stays at 8.

## 3. What this means for the thesis

- **Two knobs are selected by the experiment (c, E₀)**, both measured at
  B2 / lq. τ, s_eff, the ladder and the pickup are pinned or sourced;
  their influence is reported, not fitted (CALIBRATION_MAP tally). τ and
  s_eff (explored 5–30, 5–8 consistent, §2.1) both act on the KE ↔ large-n-tail trade-off.
- **KE₁ is a map property.** No geometry lever moves it; the measured map
  levers either cost the landing (E₀, τ, Coulomb re-anchor) or the TDDFT
  consistency (shallow well) — the deficit is an outlook item (flight
  mass M1, channel mixture, Tier-3 noise).
- **The high-n tail** needs a raised exposure ceiling *and* a blurred
  escape threshold; no measured lever supplies both
  (`THESIS_DISCUSSION.md` §7).
- **Transferred rows (T)** are qualitative statements at B2. Quote their
  magnitudes only with the law they were measured at.
