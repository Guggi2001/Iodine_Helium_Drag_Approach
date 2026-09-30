# TIER2_DetectorStageFix — retire the E2 relaxation stage?

**Status:** OPEN — Step 0 (zero-MD read) done 2026-09-29 on h405, extended
2026-09-30 to b031 + linclones s1–s3 (§2d), Coulomb closure tested against
E2 (§2e), retained split explained (§2f); decisions (1)/(2) confirmed, (3)
revised (§3a). Validation MD (h405 + b031, skip path, t_h = 502.84 ps) run
2026-09-30: **PASS** (§3d). Decisions for production + the new thesis
reference battery recorded in **§5 (fresh-chat handoff)**. Next: settle the
implementation design (§5.3), then `[PROCEED TO IMPLEMENTATION]`.

## 1. The problem

The h405 pipeline is three stages: **Stage I** ion MD (30 ps, full physics),
**Stage II / E2** relaxation (8000 ps, `relaxation_stage.py`), **Stage III**
Gillespie detection stage (to 8.53 µs, `detection_stage.py`).

E2 runs `biphasic_step` under a *relaxation view*
`replace(cfg, pickup_rate_coefficient=0.0, dt_ion=dt_relax)` — **pickup is
forced off**, while density-gated cooling, RRK evaporation, Landau-gated drag
and the full conservative translation (Coulomb + droplet well) stay live.

History (from `drag_migration_log_tier2.md` l. 212–217 and the
`relaxation_stage.py` docstring):

- E2 was decided 2026-06-29 — **before** the detection stage existed
  (~2026-07-07) — with the mission "reach the terminal size distribution by
  matched-time integration", as the "free-flight, evaporation-only limit" with
  pickup **and** drag off. Premise: all ions are already ejected at the end of
  Stage I (true at the old R ≈ 27 Å geometry).
- Forcing λ₀ = 0 was **redundant where the premise holds** (pickup is
  ρ̂-gated, so it is dead outside the droplet anyway) and **wrong where it
  fails** (ions still in helium).
- It also made the freeze early-exit valid (pickup off ⇒ E_int monotone ⇒
  frozen is absorbing), but in h405 the early exit never fires anyway
  (`landau_gated_drag` additionally requires every ion |v| ≤ v_L = 58 m/s,
  which an ejected ion never satisfies) ⇒ E2 always runs the full 8 ns.
- When the detection stage arrived, E2 was re-framed as a "handover bridge",
  but λ₀ = 0 was never revisited; `landau_gated_drag` (2026-07-20) later
  re-enabled drag in E2, **not** pickup ⇒ asymmetric channel set.
- At the corrected geometry (R ≈ 50 Å), ⅓–½ of ions are still in the droplet
  during E2 (`drag_migration_log_tier2.md` l. 10602–10615) — the premise fails.

E2 is by far the most expensive stage (800 000 steps vs 3000 for Stage I, same
per-step cost). No compute argument for λ₀ = 0 exists: the pickup RNG draw is
consumed every E2 step regardless.

Concern raised by the user: pickup-off could have biased the calibrated result
(E₀, τ were fit with this E2 in place).

## 2. Step 0 — zero-MD read on the committed h405 run (2026-09-29)

Run: `data/runs/9A_drag_shared_pure_cubic_N1000_tier2atlas_conf270_g4fh405`
(`ion.npz`, `relaxation.npz` 4147 stored cols, stride 1.93 ps, 29.99 → 8029.99
ps; `detection.npz`).

**Anchor (reproduced exactly vs committed
`data/runs/h2b_forward_model/atlas_g4finals_table.csv` h405 row):**
num_scored 1826, nbar_det 3.9551, n1_solv 0.2078, supp 0.1829 (scorer
convention: `suppressed` ions score in bin 0).

### 2a. Missed pickups during E2

Per ion, along its stored (pickup-off) E2 trajectory, with the package's own
`lambda_attach` / `rho_he_ratio` / `drag_gate_steepness` at the h405 pins
(λ₀ 0.9/ps, Langmuir p = 1, `density_only`, steepness 14.2 Å), trapezoid over
the stored columns:

$$N_\text{miss,i}=\int_{E2}\lambda_0\,\hat\rho(d_i(t))\,(1-n_i/n^*)_+\,dt$$

| class | N | mean N_miss | frac > 0.1 |
|---|---|---|---|
| **scored (detected)** | 1826 | **0.0003** | 0.05 % (1 ion, max 0.109) |
| retained (excluded) | 174 | ~1270 | 98 % |

- **Scored ensemble total: 0.55 expected missed pickups over 1826 ions** —
  ≲ one ion off by one He. 98 % of scored ions already have ρ̂ < 10⁻³ at E2
  start; the scored exposure accrues 50 / 90 / 99 / 99.9 % by
  36 / 57 / 117 / 163 ps.
- The single exposed scored ion lands at n = 13; **no low-n bin** (where the W₁
  residual lives) holds an exposed ion.
- Retained ions sit in helium for the full 8 ns (pickup-off is grossly wrong for
  them) but are excluded by `exclude_all_coupled`, and pickup (mass gain,
  at-rest capture) could only slow them further — it cannot promote one into
  the scored set.

**Verdict: pickup-off in E2 is physically unjustified scaffolding but
numerically negligible for the h405 scored result.** No MD A/B needed for the
h405 number (predicted Δ ≈ 1/1826).

### 2b. Would an earlier handover pass the detection guard?

Replaying the P1–P3 guard (`EPS_DRAIN` = 1e-6; exposure ρ̂·max(λ₀,1/τ)·(t_detect
− t_h); cooling ρ̂·(t_detect − t_h)/τ on non-frozen ions) on the stored E2 state:

| t_h (ps) | scored violators | retained violators | scored non-frozen |
|---|---|---|---|
| 49 | 64 / 1826 | 174 / 174 | 1606 |
| 100 | 6 | 174 | 1412 |
| 200 | 1 | 174 | 1223 |
| 300 | 1 | 174 | 1136 |
| **500** | **0** | 174 | 1038 |
| 1000 – 8030 | 0 | 174 | 924 → 722 |

At t_h ≈ 500 ps the scored / excluded partition is **identical** to the
committed one (all 174 retained are the violators, no scored ion is). Ions
still non-frozen at handover continue in the Gillespie stage, which is the
*exact* solver of the same Markov process E2 approximates with fixed dt
(cooling is dead outside the droplet under `density_scaled`).

### 2c. Residual Coulomb not propagated by Stage III

Half pair-Coulomb per scored fragment (14.3996/|r₁−r₂| / 2 eV):

| t (ps) | median | 95 % | max | median / detected KE |
|---|---|---|---|---|
| 200 | 0.0025 | 0.0048 | 0.0157 | 0.52 % |
| 500 | 0.0010 | 0.0019 | 0.0055 | 0.21 % |
| 1000 | 0.0005 | 0.0010 | 0.0025 | 0.10 % |
| 8030 | 0.0001 | 0.0001 | 0.0003 | 0.01 % |

At 500 ps the un-propagated Coulomb is ~0.2 % of KE — small; could be closed
analytically (asymptotic-KE helper already in `detection_stage.py`,
half-credit split).

Estimates caveat: 2b/2c are read along the pickup-off / Landau-drag E2
trajectory, not along a hypothetical extended Stage I.

### 2d. Generalisation beyond h405 (2026-09-30)

§2a–2c are h405-specific. Replayed on the committed runs that bracket the axes
controlling helium exposure — **drag strength** (b031: `capped_cubic`,
v_c = 6.0, τ 6.4, E0 0.31, trap 0.203 — the slowest-exit finalist) and **drag
form** (linclones s1–s3: `pure_linear`, N = 500, three seeds). The v_c = 5.5
τ/E0 ladder (f1–f3, a037, h345–h415) was not replayed: same drag and geometry
as h405, so the same exit kinematics. All runs are at the corrected geometry
(`legacy` prior, `raw` sampler), E2 `landau_gated_drag`, λ₀ 0.9/ps from cfg.

**Anchors (exact vs committed rows; scorer convention: `suppressed` → bin 0):**
h405 1826 / 3.9551 / 0.2078; b031 1594 / 4.2629 / 0.1771
(`atlas_g4finals_table.csv`); s1 945 / 3.9122 / 0.1937, s2 932 / 4.1041 /
0.1734, s3 937 / 3.9648 / 0.2095 (`atlas_linclone_table.csv`) — num_scored /
nbar_det / n1_solv. The h405 replay reproduces §2a/§2b exactly.

**Missed pickups (scored ions, §2a formula):**

| run | scored | Σ N_miss | ions > 0.1 | max | exposed ions' detected n | 99 % accrued by |
|---|---|---|---|---|---|---|
| h405 | 1826 | 0.55 | 1 | 0.11 | 13 | 117 ps |
| **b031** | 1594 | **1.81** | 3 | 0.26 | 10, 11, 13 | 109 ps |
| lin s1 | 945 | 0.62 | 3 | 0.25 | 13, 14, 14 | 116 ps |
| lin s2 | 932 | 0.66 | 2 | 0.26 | 13, 15 | 120 ps |
| lin s3 | 937 | 0.90 | 1 | 0.48 | 14 | 80 ps |

Retained (excluded) ions: 95–99 % exposed in every run (mean N_miss
~1400–1700) — pickup-off is grossly wrong for them everywhere, and irrelevant
to the score everywhere.

**Handover guard (§2b formula), scored violators at t_h = 100 / 200 / 300 /
500 ps:** b031 7 / 2 / 1 / **0**; s1 4 / 1 / 0 / **0**; s2 5 / 2 / 1 / **0**;
s3 5 / 0 / 0 / **0**. Retained violators = all retained (406 / 55 / 68 / 63)
at every t_h ⇒ at 500 ps the partition equals the committed one in all five
runs.

**Residual half pair-Coulomb (§2c formula), median / median detected KE:**

| run | median KE_det | 500 ps | 1 ns |
|---|---|---|---|
| h405 | 0.493 eV | 0.0010 eV (0.2 %) | 0.1 % |
| **b031** | 0.233 eV | 0.0015 eV (**0.6 %**) | 0.3 % |
| lin s1–s3 | 0.42–0.43 eV | 0.0010 eV (0.23 %) | 0.12 % |

The absolute residual is geometry-set (~1 meV at 500 ps everywhere); the
*relative* bias scales as 1/KE, so slow cells are hit harder — a one-sided,
cell-dependent downward KE bias if left open.

**Verdict:** pickup-off in E2 is score-inert across drag strength and drag
form, not only at h405 ⇒ the calibration (E0, τ, v_c ranking) is not biased
by it. A 0.5 ns handover reproduces the scored/excluded partition in all five
runs. The residual Coulomb should be closed analytically (§3a (3)), not
accepted.

### 2e. Coulomb closure — candidates tested against E2 ground truth (2026-09-30)

Outside helium E2 is conservative (pair Coulomb + flat well; the drag gate
g(depth) ≈ 0), so for a scored ion the committed E2 state at 8030 ps —
projected to infinity with the exact two-body closure (residual ~6·10⁻⁵ eV) —
is the ground truth for any closure applied at t_h ≈ 500 ps. Error = per-ion
KE-equivalent deviation from that truth; three classes because the right
closure depends on the partner's fate:

| class (h405 / b031 / lin s1) | none | ½ each | full to this ion | exact two-body |
|---|---|---|---|---|
| **A** both scored, no evaporation after t_h (1124 / 794 / 604) | −0.19 / −0.60 / −0.24 % | median 0, **p95 0.8–2.6 %, max 36–107 %** | +0.15 / +0.52 / +0.20 % | **0.0000 %** |
| **B** partner retained in the droplet (148 / 282 / 43) | −0.52 / −1.54 / −0.73 % | −0.26 / −0.77 / −0.37 % | **p95 ≤ 0.06 %** | p95 0.17–0.32 % |
| **C** scored, evaporates after t_h (316 / 320 / 175) | −0.31 / −0.95 / −0.30 % | p95 1.3–2.7 % | +0.22 to +0.48 % | **p95 0.03–0.20 %, median ≈ −0.005 %** |

- **½ each is unbiased in the median but wrong per ion.** The per-fragment
  share of the pair work is fixed by momentum conservation (fragment masses
  = n, and the pair's centre-of-mass motion), not ½/½ — so it would bias the
  **per-n KE observables** (midHot, deepKE) even though the ensemble median
  is right. Rejected.
- **Both fragments free (A, C) → exact two-body asymptote.** CM moves
  uniformly; relative motion is the repulsive Kepler hyperbola:
  $w_\infty^2 = |\mathbf v_{rel}|^2 + 2\alpha/r$, $\alpha = k_e/\mu$
  [Å³/ps²] with $k_e$ = 14.39964548 eV·Å converted to amu·Å³/ps²; asymptotic
  direction from the Laplace–Runge–Lenz vector
  $\mathbf e = -(\mathbf v_{rel}\times\mathbf h)/\alpha - \hat{\mathbf r}$,
  $\hat{\mathbf u} = (-\mathbf e + c\,\mathbf e\times\hat{\mathbf h})/(1+c^2)$,
  $c = w_\infty|\mathbf h|/\alpha$; then
  $\mathbf v_{1,\infty} = \mathbf V + (m_2/M)\,w_\infty\hat{\mathbf u}$,
  $\mathbf v_{2,\infty} = \mathbf V - (m_1/M)\,w_\infty\hat{\mathbf u}$.
  Exact against E2 to 10⁻⁴ %; residual error only from mass change after t_h
  (class C, p95 ≤ 0.2 %).
- **Partner retained (B) → full credit to the free ion.** The droplet is
  fixed at the origin (infinite mass) and holds the partner, so the partner
  does not recoil: the free ion takes the full k/r (fixed-centre repulsion).
  p95 ≤ 0.06 %.
- **Well climb:** U(depth) − U_∞ is 0 to machine precision for every scored
  ion at 500 ps — no well term needed (assert it, don't model it).
- **Both fragments retained:** no closure (excluded).
- **Post-hoc equivalence:** under `co_moving` the detection stage never
  changes velocity (verified: detected v ≡ handover v for all 1826 h405
  scored ions, 287 of which changed mass there). So the closure applied to
  the detected velocity with the detected mass is **identical** to applying
  it at handover — it can be evaluated on any finished run without code in
  the stage.

What this validates: the closure reproduces what E2 would have integrated. It
does not validate the E2 physics itself (point charges, fixed droplet, no
polarisation of the droplet).

### 2f. What the bound / marginal split actually measures (2026-09-30)

Classification (`detection_stage.py`): guard violators → `escape_energetics`
compares E_tot = KE + U(depth) + **E_c(pair) credited in full to both
fragments** (deliberate conservative over-estimate) against the outward
effective-potential barrier (well + centrifugal). Below → `droplet_retained`
(bound); above → `droplet_retained_marginal` (under `exclude_all_coupled`).

| | h405 | b031 | lin s1 |
|---|---|---|---|
| bound / marginal at t_h = 500 ps | 172 / 2 | 354 / 52 | 53 / 2 |
| bound / marginal at 8 ns (committed) | 170 / 4 | 354 / 52 | 53 / 2 |
| marginal ions whose partner is also marginal | **4 / 4** | **52 / 52** | **2 / 2** |
| bound ions whose partner is also retained | 0 | 4 | 0 |

- **Marginal = both fragments of one molecule trapped** (a doubly charged
  droplet): r_sep ≈ 110–130 Å (opposite sides, 10–24 Å below the surface),
  E_c ≈ 110–130 meV double-credited to each ion lifts them +2 to +65 meV
  over the barrier, while both sit below v_L and oscillate in place for the
  full 8 ns. The "unbound" verdict is an artefact of the full-credit rule.
- **Bound = one fragment trapped, partner escaped:** E_c ≈ 0.2–5 meV, ~110 meV
  below the barrier (h405 / b031; ~45 meV lin s1).
- The split therefore counts **molecules with one vs two fragments
  retained**, mislabelled by the conservative rule — not binding depth.
- Handover-time sensitivity: h405 one both-retained molecule (2 ions,
  margin ≈ 0) flips bound (500 ps) → marginal (8 ns); b031 and lin s1
  identical. The total retained count is identical in all three.
- Doc drift: the `escape_energetics` / `_conservatively_bound` docstrings
  still say "E2 is zero-gamma" — false since `landau_gated_drag`, and
  meaningless on the skip path. The criterion remains conservative (drag only
  makes escape harder); the docstring should be corrected when the stage
  change is implemented.

## 3. Proposed fix (to discuss)

**Stage I to ~0.5 ns (full physics, pickup live) → skip E2
(`relaxation_stage_enabled=False`, the detection stage then seeds from
`ion.npz` — an already-supported path) → Stage III.** Cost ≈ 50 000 steps vs
803 000 today (~16× less). Stage I auto-strides its checkpoint under the byte
budget (`ion.py` `_decide_stride_ion`), so storage is not a blocker.

Differences an extended Stage I introduces vs today's E2 (all expected small
for the scored set, per §2):

1. pickup live 30–500 ps (§2a: ≤ 0.55 events);
2. Stage-I drag is `capped_cubic` with **no Landau gate**; E2 used the
   Landau-gated **same `capped_cubic` bundle** (`_make_relaxation_gamma_fn`
   wraps `drag_gamma(cfg.drag_coefficients)`; `[Corrected 2026-09-30]` — this
   item previously said "pure cubic") — so the only drag difference is the
   gate, which only matters in helium, i.e. essentially for the excluded
   retained class; decided in §3a (1);
3. 0.5–8 ns evaporation moves from fixed-dt E2 to exact Gillespie;
4. residual Coulomb after 0.5 ns not propagated (§2c);
5. the retained class is classified at 0.5 ns instead of 8 ns — bound vs
   marginal split (`trap_bound` / `trap_marg`) will shift; both are excluded, so
   the score is unaffected, but the reported split changes.

RNG note `[Corrected 2026-09-30]`: Stage I and E2 consume the same *pattern*
of draws per step, but **not the same stream** — Stage I continues
`default_rng(cfg.seed)`, while E2 spawns its own
`stage_stream_rng(cfg.seed, RELAXATION_STREAM_KEY)` (`relaxation_stage.py`
l. 128, 494). An extended Stage I is therefore draw-aligned with the committed
run only for 0–30 ps; every stochastic event after 30 ps is an independent
re-draw (see §3b for the consequence). The earlier claim "draw-aligned until
the first physics divergence — a near-CRN comparison" was wrong.

**Validation proposal:** one MD run = h405 with `ion_simulation_time` 500 ps,
relaxation disabled, detection as-is; compare the full observable vector
(num_scored, trap, supp, nbar_det, n1_solv, w1_solv, midHot, deepKE, chi2)
against the committed h405 row; acceptance = within the h405 seed-to-seed
spread. Config-only (generator constants), but it is a pipeline change →
decide in the re-discussion before running. *(Superseded by the paired design
in §3b.)*

### 3a. Discussion of the open decisions (2026-09-30)

Recommendations; none is adopted until the user confirms.

**(1) Landau gate in Stage I — recommend: no (keep Stage I as calibrated).**
Scored-inert (≥ 98 % of scored ions have ρ̂ < 10⁻³ at 30 ps; sub-v_L cubic
friction is ∝ v³ ≈ 0 anyway). Leaving Stage I untouched keeps 0–30 ps
**bit-identical** to the committed `ion.npz` — a free hard anchor for §3b —
and keeps the pipeline change to one variable. A time-switched gate (on only
after 30 ps) would be scaffolding. Consequence to record: with E2 skipped,
`relaxation_dissipation` has no production reader, and v_L loses its *drag*
reader → the v_L "P" row in `TIER2_PARAMETER_INFLUENCE.md` §17 must say so.
`[Corrected 2026-09-30]` `v_limit_m_per_s = 58` itself **stays read** in
production: it also sets `cfg.E_min_eV`, the Landau cutoff of the hard-sphere
collision step in the neutral stage (`propagation_step.py`) — so the field
cannot be retired, only its E2 reader (`docs/POST_THESIS_CLEANUP.md` C1, C4).
A physics case for a global Landau gate on the drag is a separate enum change
with its own validation.

*A1–A4 decided 2026-09-30:* **A1** fixed handover ≈ 0.5 ns **plus a
safeguard**: fail loudly if any `droplet_retained_marginal` ion has a partner
that is not retained (the signature of an escaper cut off by too early a
handover; 0 occurrences in 88 finished runs, whose 1867 marginal ions all
have a retained partner — known blind spot: both fragments of one molecule
slow escapers at t_h). **A2** closure in the detection stage at handover,
behind a new enum, default = today's behaviour (note: the closure converts
pair Coulomb *potential* into KE — `E_pot_detected` must drop by the same
amount to keep the 5-term ledger closed; check how the pair term is split
per ion first). **A3** keep E2 as legacy (banner), removal registered in
`docs/POST_THESIS_CLEANUP.md` (C1). **A4** pipeline-version stamp in
`cfg.json`; no mixing of E2-era and skip-path rows in one CRN comparison.

**(2) Validation MD — recommend: one run, judged by paired (CRN) checks, not
seed spread.** The predicted effect (~1/1826) is far below the h405 seed
spread, so a seed-spread criterion cannot fail and proves nothing. Design and
pre-registered predictions in §3b.

*User decisions 2026-09-30:* (1) **confirmed** — no Landau gate in Stage I
for now (may be revisited later as its own change); (2) **confirmed** —
paired/statistical validation, not seed spread (design under discussion,
§3b); (3a) closure design **replaced** after critical review (below, §2e);
(3b) detailed in §2f.

**(3a) Residual Coulomb — close with the partner-aware exact closure, do not
accept.** One-sided (KE low), cell-dependent (0.2 % h405 → 0.6 % b031, §2d),
and in the direction of the KE₁ deficit under investigation
(`TIER2_MASS_SCENARIOS.md`). ~~Credit each fragment ½·E_coul at fixed
direction~~ — **rejected** (§2e: median-unbiased but p95 0.8–2.6 % per ion,
biases the per-n KE observables). Replacement, per molecule at handover:
- both fragments scored → **exact two-body asymptote** (CM + repulsive-Kepler
  / LRL; magnitude *and* direction; momentum-conserving split);
- one fragment retained → **full k/r to the free fragment** (fixed-centre;
  the droplet absorbs the partner's recoil);
- both retained → nothing (excluded).
Validated zero-MD against E2 (§2e): exact for non-evaporating pairs, p95 ≤
0.2 % for pairs that evaporate after t_h, ≤ 0.06 % for retained partners; the
result becomes t_h-independent and more complete than today's 8 ns run
(0.01 % open). Under `co_moving` it is exactly equivalent post hoc (§2e). A
code change behind its own enum (default = today's behaviour) — requires
`[PROCEED TO IMPLEMENTATION]`; placement (detection stage at handover vs the
scorer) to be decided then.

**(3b) Bound/marginal split — recommend: accept; report it as "one fragment
retained" vs "both fragments retained".** Score-inert (both excluded). §2f
shows the split is not binding depth but the number of fragments a molecule
leaves in the droplet, mislabelled by the full-credit rule; its handover
sensitivity is one molecule (h405) or none (b031, lin s1), and the total
retained count is unchanged. Footnote `trap_bound` / `trap_marg` as
pipeline-dependent; stamp a pipeline version and never mix E2-era and
Stage-I-0.5 ns rows in one CRN comparison. A pair-aware escape criterion
(molecule total energy) is possible but not needed for the score.

### 3b. Validation MD — design (DRAFT, under discussion 2026-09-30; not run)

*Revised 2026-09-30:* the first draft assumed draw alignment with E2 after
30 ps and pre-registered per-ion equality of (n, E_int, r, v) at t_h. That is
void (RNG note, §3): every stochastic event after 30 ps is re-drawn. Only 13
scored ions per run (h405 and b031 alike) have no stochastic event after
30 ps (frozen at 30 ps with a frozen partner); 1149 / 1189 scored ions
evaporate between 30 and 500 ps. The n-observables must therefore be judged
statistically; the kinematics stay near-paired (co-moving sheds do not change
velocity, so a re-drawn shed perturbs an ion's trajectory only through the
mass-dependence of the residual Coulomb acceleration).

Run(s): generator config, `ION_TIME_PS` ≈ 500 ps set to an exact stored E2
column time (equal-time comparison), relaxation disabled (skip path),
detection as-is (closure applied post hoc, exact under `co_moving`, §2e), same
seed, new run tag. Proposed: **h405** (reference) **and b031** (stress case:
406 retained, 52 marginal, 0.6 % Coulomb) — two concurrent runs, ~50 000
steps each.

Predictions, in order:

1. **Bit-identity 0–30 ps** at common stored columns vs committed `ion.npz`.
   Fail ⇒ stop.
2. **Frozen-pair ions (13 per run):** state at t_h equals the committed E2
   column to round-off (≲ 10⁻¹⁰ relative).
3. **Kinematics, all scored ions:** |Δv|/|v| at t_h vs the committed E2 column
   small — estimate ≲ 10⁻³ (residual Coulomb speed gain after 30 ps ~1.5 % ×
   re-drawn mass difference ~2–3 %); exact tolerance to fix before running.
4. **Partition:** same excluded ion ids (174 / 406), zero scored violators;
   bound/marginal per §2f (h405 172/2, b031 354/52 predicted).
5. **n-observables:** nbar_det within 2σ of committed, σ = Poisson upper
   bound on the post-30 ps resampling noise (√(2·⟨events⟩/N): 0.035 h405,
   0.041 b031); physics-predicted shift ≈ 1/N (pickup). n1_solv, w1_solv,
   supp: tolerances to fix before running.
6. **KE observables (after closure):** per-ion KE agrees with the committed
   run closed the same way to the item-3 level; per-n-bin KE differences only
   from ions re-binned by different n.

Seed spread is reported only as the "does it change any conclusion" context.

*User decision 2026-09-30:* runs h405 + b031 together; **sharp checks 1–4
decide**, 5–6 are reported; `[PROCEED TO IMPLEMENTATION]` for the generator
and the runs.

### 3d. Validation MD — results (2026-09-30)

Generator `scripts/gen_detstagefix_validation.py` (tests
`tests/test_gen_detstagefix_validation.py`, 8 passed): the finals `build_cell`
verbatim + `ion_simulation_time = 502.845`, `relaxation_stage_enabled =
False`. Pre-MD guards: base cell ≡ committed `cfg.json` (empty diff);
validation diff = exactly those two keys; **t_h = 502.84 ps** = committed E2
column 245 (step 50284); storage stride 13 divides 50284, so the detection
seed column is the final state at t_h. Run dirs
`9A_drag_shared_pure_cubic_N1000_detfix_conf270_{h405,b031}_th503`.
**Wall-clock ≈ 2.5 min per run** (committed h405: ≈ 42 min, of which E2
≈ 41 min) — ≈ 17× cheaper.

**Sharp checks (decide):**

| check | h405 | b031 |
|---|---|---|
| 1 `neutral.npz` identical; `ion.npz` 231 common columns (t ≤ 29.90 ps) bit-identical | **pass** | **pass** |
| 2 no-event scored pairs bit-identical at t_h (revised criterion, below) | **184 / 188**; 4 explained | **54 / 54** |
| 3 kinematics \|Δv\|/\|v\| at t_h, median / p99 / max | 5.1·10⁻⁵ / 1.6·10⁻³ / 2.4·10⁻² | 2.0·10⁻⁴ / 4.2·10⁻³ / 2.4·10⁻² |
| 4 excluded ion ids identical; bound / marginal (old → new) | **pass** (174); 170/4 → 170/4 | **pass** (406); 354/52 → 352/54 |

- **Check 2 criterion revised after the run (disclosed):** the pre-registered
  "frozen at 30 ps (E_int < D₀), frozen partner" does **not** exclude events —
  9 / 13 (h405) and 11 / 13 (b031) of those pairs had mass events in the new
  run (pickup is live). The correct criterion is "no mass change after 30 ps
  for the ion *and* its partner in *either* run" (from `mass_history_kg`).
  Under it, every pair is bit-identical after ~47 000 steps except four h405
  ions, all explained by the **Landau gate** (the one intended physics
  difference): ion 800 is a slow escaper that drops to 0.63 v_L inside the
  surface tail (ρ̂ ≤ 2·10⁻²) — ungated cubic drag in Stage I, frictionless in
  E2 (|Δv|/|v| 9.7·10⁻³); 1800 is its scored partner (10⁻⁶); 261 and 1394
  have sub-v_L droplet-retained partners at ρ̂ ≈ 1 whose motion differs by
  the gate, felt through the pair Coulomb (≤ 2.5·10⁻⁶).
- **Check 3:** the pre-registered ≲ 10⁻³ holds for the median, not the tail.
  The tail comes from pairs with re-drawn events shortly after 30 ps, while
  the residual Coulomb is still strong. After the closure (which maps each
  state to its conserved asymptote) the per-ion KE of ions with the same
  detected mass agrees to **median 10⁻¹² (h405) / 1.3·10⁻⁵ (b031), p99
  0.3 % / 0.7 %** — i.e. most of the raw Δv is Coulomb phase, not a
  different asymptote.
- **Check 4:** b031's bound/marginal moves by one both-retained molecule —
  exactly the §2f mechanism (margin ≈ 0 oscillation).
- **Pickups after 30 ps in the new run** (n increases between stored
  columns): scored ions **1 (h405), 2 (b031)** vs the §2a predictions of
  0.55 and 1.81 expected events. Retained ions: 171 / 405 pick up, as
  expected (excluded).

**Observable vector (reported, committed scorer; anchor: re-scoring the
committed `detection.npz` reproduces the committed rows exactly):**

| | nbar_det | n1_solv | w1_solv | midHot | deepKE | chi2_med |
|---|---|---|---|---|---|---|
| h405 committed | 3.9551 | 0.2078 | 0.7089 | 0.9567 | 0.6026 | 349.63 |
| h405 committed + closure | 3.9551 | 0.2078 | 0.7089 | 0.9569 | 0.6039 | 348.81 |
| h405 new, raw | 3.9584 | 0.2058 | 0.7076 | 0.9529 | 0.5832 | 360.95 |
| **h405 new + closure** | **3.9584** | **0.2058** | **0.7076** | **0.9572** | **0.6057** | **346.57** |
| b031 committed | 4.2629 | 0.1771 | 1.1134 | 0.5198 | 0.7316 | 68.26 |
| b031 committed + closure | 4.2629 | 0.1771 | 1.1134 | 0.5203 | 0.7341 | 68.01 |
| b031 new, raw | 4.2629 | 0.1784 | 1.1160 | 0.5117 | 0.7073 | 72.94 |
| **b031 new + closure** | **4.2629** | **0.1784** | **1.1160** | **0.5195** | **0.7570** | **68.31** |

num_scored, trap, supp unchanged in both. Noise scale: σ(nbar) ≤ 0.035 /
0.041 (Poisson bound), σ(n1_solv) ≈ 0.015; deepKE sampling SE ≈ 0.05 (h405)
/ 0.09 (b031) (half-sample estimate).

- **n-observables** differ only through re-drawn *late* evaporations: 30
  (h405) / 26 (b031) scored ions end ±1 He apart, all flipping frozen ↔
  `time_exhausted` near t_detect (net +6 He h405, 0 b031). All shifts are
  ≪ 1σ.
- **The closure is necessary:** without it, KE observables are biased low
  (h405 deepKE −0.019, chi2 +11; b031 midHot −0.008, chi2 +4.7); with it,
  midHot agrees to ≤ 0.0008 and chi2 to ≤ 2.3.
- **b031 deepKE +0.023:** decomposed — new kinematics on the committed n
  assignment gives 0.7228 (−0.011), so the shift is mainly the 26 re-binned
  ions in its 5 sparse deep bins; ≪ its sampling SE (0.09).
- **Side finding:** the committed E2-era runs are themselves slightly
  KE-low: closing them at 8 ns moves scored KE +0.014 % (h405) / +0.05 %
  (b031) median, deepKE +0.0013 / +0.0025 — the "0.01 % open" in §2c
  under-stated b031. Below every decision margin.

**Verdict: PASS.** Sharp checks 1–4 pass on both cells (check 2 under the
disclosed corrected criterion; every non-identical ion traced to a redrawn
event or the Landau gate). The reported observables agree with the committed
runs within the post-30 ps resampling noise once the §2e closure is applied.
The proposed pipeline (Stage I to ~0.5 ns, E2 skipped, partner-aware Coulomb
closure) reproduces the committed h405 and b031 results at ≈ 1/17 of the cost.

### 3c. Next steps

1. ~~User confirms §3a (1)–(3b).~~ (1), (2) confirmed 2026-09-30; (3a)
   revised to the partner-aware closure; (3b) detailed (§2f).
2. ~~Settle the §3b design, then run it.~~ Done 2026-09-30 — **PASS** (§3d).
3. Next, on the user's trigger: `[PROCEED TO IMPLEMENTATION]` for (a) generator constants
   (`ION_TIME_PS`, relaxation off) with a pipeline-version stamp, and (b) the
   Coulomb-closure enum + focused tests.
4. Ledger / doc updates: v_L row in `TIER2_PARAMETER_INFLUENCE.md` §17; E2
   retirement entry in `drag_migration_log_tier2.md`; decide whether
   `relaxation_stage.py` is kept as a non-production enum or retired (rule 2).
5. Thesis: E2-era atlas results stand (§2d verdict); state the stage change
   and the §2d numbers as the justification.

## 5. Fresh-chat handoff (2026-09-30)

**Entry point for the next chat.** §1–§3d are the evidence; this section is
the decision record and the open design agenda. Nothing below §5.2 is
implemented yet. Working rules: CLAUDE.md (no code before
`[PROCEED TO IMPLEMENTATION]`); commit on `drag_implementation`; max two
concurrent MD runs.

### 5.1 Where things stand

- E2 audit (§2a–2f) and the paired validation MD (§3d) are **done — PASS**.
- Repo files from this phase: this doc; `docs/POST_THESIS_CLEANUP.md` (C1–C4);
  `scripts/gen_detstagefix_validation.py` +
  `tests/test_gen_detstagefix_validation.py` (the validation generator,
  t_h = 502.84 ps — a validation device, **not** the production value).
- Validation run dirs (gitignored):
  `data/runs/9A_drag_shared_pure_cubic_N1000_detfix_conf270_{h405,b031}_th503`.
- The closure maths as working analysis code (not yet in the package):
  scratchpad `vmd_checks.py` (`two_body_asymptote`, `closure_dv`,
  `closed_detection`) — see §4 for the path. §2e has the formulas; the
  package implementation must be written fresh from §2e, not copied blindly.

### 5.2 Decisions (user, 2026-09-30)

| id | decision |
|---|---|
| Landau gate | **Not** in Stage I for now (may be revisited later as its own change). |
| Validation | Paired design; sharp checks decide — done, PASS (§3d). |
| A1 handover | **Fixed t_h = 500 ps** + **safeguard**: fail loudly if any `droplet_retained_marginal` ion has a partner that is not retained (signature of an escaper cut off by too early a handover; 0 of 1867 marginal ions in 88 finished runs; blind spot: both fragments slow escapers at t_h). |
| A2 closure | Partner-aware exact closure (§2e) **in the detection stage at handover**, behind a **new config enum, default = today's behaviour** (byte-inert for existing runs). |
| A3 E2 | Keep as **legacy** (docstring banner + config comments); removal registered in `docs/POST_THESIS_CLEANUP.md` C1. |
| A4 stamp | Pipeline-version stamp; never mix E2-era and skip-path rows in one CRN comparison (see §5.3 D5 — possibly derivable from existing cfg fields). |
| C t_h value | **Round 500 ps** for production (502.84 ps was only the validation device). Handover exactly at 500.00 ps ⇒ `ion_simulation_time = 500.005` (50001 steps, last id 50000) and storage stride **16** (smallest divisor of 50000 ≥ the default 13; ≈ 3126 columns, ≈ 750 MB), set via `run_ion_propagation(max_bytes=…)` in the generator — Stage-I code untouched; guard: last stored column = final state at 500.00 ps. |
| Reference battery | **h405**, N = 1000 × **the same five seeds as the E2-era battery** (20260730–20260734, `gen_tier2atlas_g4step2_battery.py`, committed `atlas_g4step2_battery.csv`, run dirs `…g4s2h405s1…s5`, pooled `…g4s2h405pooled`) — each member CRN-paired with its E2-era twin; pooled N = 5000 = the new thesis reference. |
| Run 6 | **Include**: h405, seed 20260729, t_h = 500 ps — pairs with the validation run (502.84 ps) and the committed finals h405; answers "does 500 vs 502.84 matter" (predicted: no). Not part of the reference bundle. |
| Order | Implement (closure, safeguard, stamp, production generator, banner, docstrings, tests) **before** the runs, so the battery is a production artifact. |

### 5.3 Implementation design — agenda (settled 2026-09-30, see §5.4)

**D1 Closure enum.** Name / values (proposal:
`detection_coulomb_closure: Literal["none", "partner_aware"] = "none"`),
`_KNOWN_…` reject guard (enum-surface convention, CLAUDE.md). Placement in
`run_detection_stage`: **after** the handover classification (needs the
retained / scored flags per ion) and **before** the Gillespie loop. Per
molecule: both scored → two-body asymptote (magnitude + direction); one
retained → fixed-centre (free fragment takes full k/r; implement as the
two-body limit m_partner → ∞, v_partner = 0, or in closed form); both retained
→ nothing. Masses = handover masses. Under `co_moving` the result is identical
to post-hoc application (§2e); state whether the enum is allowed with the
`cold` shed convention (there a later shed changes v — closure-at-handover is
then an approximation; guard or document).

**D2 Energy ledger (the real design question).** Per-atom `E_pot` carries
**half** the pair Coulomb (`ion_propagation_step.py` l. 299, `_E_pot_per_atom`);
`E_pot_detected_eV` is initialised from the seed's `E_pot_eV`
(`detection_stage.py` l. 564). The closure converts the pair's E_c into KE:
- pair level: Σ ΔKE = E_c exactly (both-free case) — the 5-term invariant
  closes **per pair** if each fragment's E_pot drops by E_c/2;
- per ion: ΔKE_i ≠ E_c/2 (momentum-fixed share) ⇒ per-ion
  E_kin + E_pot shifts by an inter-fragment transfer — the same thing the MD
  pair force does every step, so check whether any test/assertion checks the
  invariant **per ion** (detection-stage docstring l. ~65: "closes across the
  stage boundary: per fire") and decide the convention;
- fixed-centre case: the free ion gains the full E_c while its E_pot holds
  only E_c/2; the other half sits on the retained (excluded) partner — decide
  where that half is booked.

**D3 Well-term assertion.** Assert |U(depth) − U_∞| below a tolerance for
every closed scored ion (it was 0 to machine precision in all reads);
tolerance value and failure message.

**D4 Safeguard.** Location (detection stage, after classification), exact
condition (`droplet_retained_marginal` & partner not in the retained
vocabulary), behaviour (raise with remedy "raise ion_simulation_time"),
scope: all paths/policies (it would have fired 0× historically ⇒ byte-inert),
or only when E2 is off. Unit test on a synthetic handover state.

**D5 Version stamp.** A new `SimConfig` field (cfg.json only; no checkpoint
schema change) vs deriving the pipeline identity from existing fields
(`relaxation_stage_enabled`, the closure enum) — the cfg may already be
self-describing, making a stamp redundant (rule 2).

**D6 Production generator.** New script (proposal
`scripts/gen_tier2_detfix_battery.py`): `build_cell(h405)` + seed swap (the
battery pattern) + `ion_simulation_time = 500.005`,
`relaxation_stage_enabled = False`, closure on; stride via `max_bytes`;
guards: cfg diff vs the E2-era battery member's `cfg.json` = exactly the
pre-registered key set; alignment guard (last stored column = final step at
500.00 ps); members s1–s5 + run 6 (seed 20260729); run-dir naming distinct
from the E2-era namespace; concurrency 2.

**D7 Legacy banner + C2 docstrings.** `relaxation_stage.py` banner;
`relaxation_*` config comments; fix the "E2 is zero-gamma" text in
`escape_energetics` / `_conservatively_bound` / `EscapeEnergetics`.

**D8 Tests.** Closure (two-body exactness vs a direct ODE integration on a
synthetic pair; fixed-centre limit; both-retained no-op; default `"none"`
byte-inert on an existing detection path), ledger convention (D2), safeguard
(D4), generator guards (D6).

**D9 Scoring + checks for the battery.** Per member: 0–30 ps bit-identity vs
its E2-era twin (sharp); excluded ids; pooled full observable vector (commit
a CSV with **all** columns incl. KE, per member + pooled, memory
"scorer CSVs carry full vector"); paired old/new table, with the E2-era
battery also closed post hoc for a like-for-like KE comparison. Run 6:
bit-identity of its `ion.npz` vs the validation run's at all common stored
times ≤ 500 ps; closed per-ion KE agreement.

### 5.4 Implementation design — settled (user adopted all recommendations, 2026-09-30)

Not implemented yet; code waits for `[PROCEED TO IMPLEMENTATION]`, runs for
the user's go.

**D1 Closure enum.**
- `detection_coulomb_closure: Literal["none", "partner_aware"] = "none"`;
  `_KNOWN_DETECTION_COULOMB_CLOSURES` rejected via `_reject_unknown_enum` at
  the top of `check_detection_config` (the retained-policy pattern).
- Config-load refusals: `partner_aware` with `evaporation_shed_convention =
  "cold"` (later sheds rescale v ⇒ closing at handover is unvalidated there;
  production is `co_moving`); `partner_aware` with the detection stage
  disabled (silently inert — the `config.py` l. 1876 precedent). Allowed on
  the E2 path (closing an 8 ns seed is correct physics, §3d side finding).
- Placement: `run_detection_stage`, after the retained classification, before
  the event loop (the loop starts from closed velocities). Partner of ion i =
  (i + N) mod 2N. "Free" = in neither retained class. Both free → exact
  two-body asymptote (§2e). One free → fixed-centre closure = the
  m_partner → ∞ limit in closed form (μ = m_free, V = 0, v_rel = v_free,
  centre = partner's handover position): speed identical to the validated
  "full k/r", plus the Kepler deflection. Both retained → no-op. Handover
  masses and positions.
- The maths lives once, in a pure physics module (proposal
  `physics/coulomb_closure.py`), called by the stage **and** by the post-hoc
  scorer for the E2-era twins (rule 1).

**D2 Energy ledger — per-ion exact booking.** For every closed ion
ΔE_pot,i = −ΔKE_i (the Coulomb work it received). Consequences:
- the per-ion 5-term invariant stays exact (the existing
  `test_detection_stage.py:853` / `test_evaporation_shed_convention.py:290`
  contract generalises unchanged);
- the pair sum is physically exact: ΣΔE_pot = −E_c (both free), and in the
  fixed-centre case the free ion's E_pot drops by the full E_c while the
  retained partner stays **verbatim** (V0-2 convention intact; the pair sum
  U₁ + U₂ + E_c − E_c is again correct);
- no new field, no `detection.npz` schema change: the closure amount is
  recoverable per ion as E_pot_detected − E_pot_seed − Σ fold;
- per-ion E_pot after closure means "handover E_pot minus the Coulomb work
  received" — the ½/½ split was a convention for a pair term anyway, and the
  per-ion sum was never conserved in Stage I (the pair force does unequal work).
- Rejected: "physical" E_pot (−E_c/2 each) — breaks the per-ion invariant and
  forces either editing the retained partner or a new ledger field.

**D3 Well-term assertion.** Only when the closure is on: for every closed free
ion, |U(depth) − E_bind| ≤ 1·10⁻⁶ eV (U_∞ = `binding_energy_I_ion_eV`,
`droplet_potential` saturates; reads were 0 to machine precision; 10⁻⁶ eV is
~10⁻⁶ of the KE, far below the closure's own class-C error). Failure:
`ValueError` naming the count, the worst ion id, its depth and gap, remedy
"raise ion_simulation_time".

**D4 Safeguard.** In `run_detection_stage` right after classification,
**unconditional** (all paths, closure on or off; can only fire under
`exclude_all_coupled`, where the marginal class exists; fired 0× in 88
finished runs ⇒ inert for every committed run): any
`droplet_retained_marginal` ion whose partner is in neither retained class ⇒
`ValueError` with the ion ids and remedy "raise ion_simulation_time (an
unbound escaper was cut off by the handover)". Known blind spot (both
fragments slow escapers) stated in the docstring.

**D5 Version stamp — none.** cfg.json is already self-describing: the pipeline
identity is (`relaxation_stage_enabled`, `detection_coulomb_closure`). A stamp
field would duplicate it (rule 2). The no-mixing rule is enforced by the
generator's cfg-diff oracle and the scorer's pairing oracle (D9).

**D6 Production generator** `scripts/gen_tier2_detfix_battery.py`.
- Members s1–s5 = seeds 20260730–20260734, plus `r6` = seed 20260729.
- cfg = `replace(build_cell(h405), seed, ion_simulation_time=500.005,
  relaxation_stage_enabled=False, detection_coulomb_closure="partner_aware")`.
- Oracles before any MD: diff vs the E2-era twin's committed `cfg.json`
  (s1–s5: `…g4s2h405sK`; r6: `g4fh405`) = exactly {`ion_simulation_time`,
  `relaxation_stage_enabled`, `detection_coulomb_closure`} (the pre-closure
  cfg.json loads with the default `"none"` — verify `cfg_diff_vs_reference`
  reports it); alignment: 50001 steps, last id 50000 × dt = 500.00 ps;
  `max_bytes` chosen so `_decide_stride_ion` gives stride 16 (asserted), and
  16 | 50000 ⇒ the seed column is the final state.
- Run dirs `…detfix_conf270_h405{s1..s5,r6}_th500` (distinct from the atlas
  `tier2atlas_` namespace); required artifacts without `relaxation.npz`;
  skip-complete / rebuild-incomplete (battery precedent); **concurrency 2**.

**D7 Legacy banner + docstrings.** `relaxation_stage.py` module banner
("LEGACY — not in the production pipeline since 2026-09-30; see this doc;
removal `POST_THESIS_CLEANUP.md` C1"); `relaxation_*` config comments point
there; fix the "E2 is zero-gamma" text in `detection_stage.py`
(`EscapeEnergetics`, `_conservatively_bound`, `escape_energetics`): E2 is the
Landau-gated `capped_cubic` bundle, and on the skip path the criterion is the
conservative energy verdict at handover — still conservative because drag
only makes escape harder.

**D8 Tests.**
- Closure maths: two-body asymptote vs a direct high-resolution ODE
  integration of a synthetic pair (magnitude and direction); momentum
  conservation; ΣΔKE = E_c; fixed-centre speed² = v² + 2k/(m r) and direction
  vs a fixed-centre ODE; both-retained no-op.
- Stage: default `"none"` leaves every output bit-identical (existing suite
  green + explicit check); per-ion 5-term closure with the closure on
  (1e-12 eV); pair ΣΔE_pot = −E_c; retained partner verbatim.
- Config: unknown value, `partner_aware + cold`, `partner_aware` with the stage
  disabled — each refused.
- D3 assertion fires on a synthetic near-surface free ion; D4 safeguard fires
  on a synthetic marginal ion with a free partner and stays silent otherwise.
- Generator: cfg-diff key set, alignment (stride 16, last column 500.00 ps),
  run-dir names, seeds.

**D9 Scorer** `scripts/post_processing/tier2_detfix_battery_table.py` (pure;
observable vector from the committed battery scorer — rule 1).
- Oracles first: rescoring the E2-era members reproduces the committed
  `atlas_g4step2_battery.csv` rows; the post-hoc closure of an E2-era run
  reproduces the §3d "committed + closure" h405 row.
- Sharp, per member vs its E2-era twin: `neutral.npz` identical; `ion.npz`
  bit-identical at common stored columns ≤ 30 ps (strides differ, so common
  steps are the common multiples); excluded ion ids identical (§3d
  precedent).
- Committed CSV: every column of the committed battery header (incl. KE₁ /
  KE₂ / S) for new s1–s5 + pooled, E2-era raw, and E2-era closed post hoc;
  plus a paired old-closed vs new table.
- Run 6: `ion.npz` bit-identical vs the th503 validation run at common stored
  times ≤ 500 ps (Stage I shares the `default_rng(seed)` stream there); closed
  per-ion KE agreement reported; prediction: no observable moves beyond the
  §3d resampling noise. Not part of the reference bundle.

### 5.5 Implementation (2026-09-30)

Delivered per §5.4. Files: `physics/coulomb_closure.py` (new),
`simulation/detection_stage.py` (closure, D3, D4, `_pair_coulomb_eV` shared
with `escape_energetics`, D7 docstrings), `simulation/ion_propagation_step.py`
(`_eV_to_amu_ang2_ps2`, the inverse of the existing conversion), `config.py`
(enum, refusals, legacy comments), `simulation/relaxation_stage.py` (banner),
`scripts/gen_tier2_detfix_battery.py`,
`scripts/post_processing/tier2_detfix_battery_table.py`; tests
`test_detection_coulomb_closure.py`, `test_gen_tier2_detfix_battery.py`.

Deviations / findings while implementing:

- **D6 oracle is two-step.** `cfg_diff_vs_reference` refuses a non-default
  value in a field an older `cfg.json` cannot witness, so the member cfg is
  diffed against the twin with the closure reset to `"none"` (exactly
  {`ion_simulation_time`, `relaxation_stage_enabled`}), and the closure is
  asserted to be the only further in-memory difference.
- **D4 changed three existing synthetic fixtures**
  (`TestExcludeAllCoupledPolicy` in `test_detection_stage.py`): they built a
  marginal ion with an ejected partner — exactly the shape D4 refuses. The
  fixtures now keep the marginal ion's partner inside (a trapped pair), the
  asserted fractions follow; the old shape is now the D4 refusal test.
- **k = 0 pass-through.** The two-body routine returns zero-coupling pairs
  verbatim (the CM recomposition otherwise leaves round-off).
- **Maths vs ODE:** the asymptote is invariant along DOP853-integrated orbits
  to ~10⁻¹¹ Å/ps (two-body and fixed-centre); the integrated velocity
  approaches it to the residual k/r.
- **Zero-MD anchor:** re-running the detection stage with the closure on the
  §3d th503 runs reproduces the §3d "new + closure" rows exactly (h405:
  midHot 0.9572, deepKE 0.6057, chi2 346.57; b031: 0.5195 / 0.7570 / 68.31)
  with `n_detected`, `state_reason` and event times bit-identical to the
  committed `detection.npz` — the independent scratch implementation and the
  package agree. D3 and D4 did not fire on real data.
- **Post-hoc closure of E2-era runs = re-running the detection stage** with
  the closure on the committed seed checkpoint (stage-private RNG, jump chain
  velocity-independent under `co_moving`); the scorer asserts the event
  record is identical, so no second closure implementation exists.

### 5.6 Reference battery — results (2026-09-30)

`gen_tier2_detfix_battery.py` (6 members, 2 concurrent, ≈ 3 min each);
scorer `tier2_detfix_battery_table.py`; committed CSVs
`data/runs/h2b_forward_model/detfix_battery_{table,paired}.csv` (full vector,
per member × {new, E2-era raw, E2-era closed post hoc} + pooled, r6 vs th503).
Full `pytest`: 3116 passed, 3 skipped (pre-existing), + 25 generator tests.

**Oracles (all pass):** each E2-era member rescored = its committed battery
row (4 dp), r6's twin = the committed finals h405 row, pooled E2-era = the
committed pooled row; pipeline identity new = `skip+closure:partner_aware`,
twins = `E2+closure:none`; post-hoc closure reproduces every twin's event
record bit-exactly; §3d anchors (committed + closure, th503 + closure)
reproduced.

**Sharp checks (all 6 members pass):** `neutral.npz` identical; `ion.npz`
bit-identical at all 188 common stored columns ≤ 30 ps; excluded ion ids
identical (173 / 136 / 145 / 152 / 174 / 174). D3 and D4 never fired.

**Run 6 (500 vs 502.84 ps):** `ion.npz` bit-identical at all 241 common
columns; identical n assignment and excluded set; closed per-ion KE median
1.8·10⁻¹³, p99 3.4·10⁻⁶; every observable equal to ≤ 4·10⁻⁶ (chi2 −0.0009).
**The handover time does not matter** (prediction confirmed).

**Paired new − E2-era closed** (observables per §3d; kinematics per ion for
ions with the same detected n: median ≈ 1.1·10⁻⁴, p99 ≈ 0.3 %):

| | Δnbar | Δn1_solv | Δw1_solv | ΔmidHot | ΔdeepKE | Δchi2 | ΔKE1 | ΔS |
|---|---|---|---|---|---|---|---|---|
| s1 | 0.000 | −0.0007 | −0.0027 | −0.0001 | +0.0014 | −0.9 | +0.0002 | −0.004 |
| s2 | −0.002 | +0.0033 | −0.0040 | −0.0008 | −0.0015 | −2.3 | −0.0008 | −0.001 |
| s3 | −0.003 | 0.0000 | 0.0000 | −0.0002 | −0.0007 | +2.4 | +0.0001 | +0.002 |
| s4 | +0.002 | +0.0014 | +0.0007 | +0.0004 | +0.0036 | +2.0 | −0.0003 | −0.002 |
| s5 | +0.002 | −0.0013 | +0.0027 | −0.0023 | **+0.0907** | **−124** | +0.0003 | +0.022 |
| **pooled** | −0.0003 | +0.0005 | −0.0007 | −0.0007 | +0.0154 | +0.2 | −0.0001 | +0.004 |

- **s5 deepKE is one ion in a 2-ion bin:** with the new kinematics on the
  E2-era n assignment deepKE is 0.4989 (vs 0.4988); with E2-era kinematics on
  the new n assignment 0.5872 — the shift is re-binning only. The n = 16 bin
  holds 2 ions in both runs but a different one (mean KE 0.008 → 0.057 eV);
  the arithmetic mean of 7 bin ratios at `min_count = 1` amplifies it. Same
  mechanism as b031 in §3d; deepKE's sparse-tail sensitivity is a scorer
  property, not a pipeline effect.
- Everything else moves within the post-30 ps resampling noise (§3d scales).

**Pooled reference (new pipeline, N = 5000, 9220 scored):** nbar 3.876,
n1_solv 0.2100, w1_solv 0.7646, midHot 0.946, deepKE 0.520, chi2 460,
KE1 0.6405, KE2 0.5491, S 1.736; trap 0.078 (bound 0.0767 / marginal 0.0013).

**Verdict:** the production pipeline reproduces the E2-era battery member by
member; the new pooled battery is the thesis reference.

**After implementation + runs:** records (TIER2_PARAMETER_INFLUENCE §17 v_L /
closure / t_h rows; migration-log entry; ≤ 2 CLAUDE.md pointer lines) and
thesis framing (§3c item 5; discussed 2026-09-30: pipeline, E2 rationale,
validation, closure, limitations incl. no sub-Landau regime in the drag).

## 4. Artifacts

Scratchpad (session-specific, may be cleaned; the reads are fully specified
above and re-derivable from the committed run dir):

`C:\Users\user\AppData\Local\Temp\claude\T--github-synchronized-Iodine-Helium-Drag-Approach-i2-helium-md\65ef168f-b9cf-45ec-9f11-20adc391b9e8\scratchpad\`

- `step0_missed_pickup.py` — §2a (+ anchor)
- `step0b_guard_vs_handover.py` — §2b
- §2c was an inline one-off (formula above).

§2d (2026-09-30, session scratchpad
`...\ee333bdc-ea9d-40f7-8194-43ac0efd7c5d\scratchpad\`):

- `step0_multi.py` — §2a + §2b + §2c on h405, b031, linclones s1–s3, each
  anchored on its committed row; λ₀ read from cfg; pairs = `[:N]`/`[N:]`
  halves.
- `closure_and_split.py` — §2e closure test (h405, b031, lin s1) + §2f split.
- `pair_check.py` — §2f partner fates of bound / marginal ions.
- `covel.py` — §2e post-hoc equivalence (detected v ≡ handover v).
- `vmd_design.py` — §3b frozen-pair counts and Poisson noise bound.
- `vmd_checks.py` — §3d checks 1–6 (closure, anchored scorer).
- `diag1.py`, `diag2.py`, `diag3.py` — §3d per-ion divergence, no-event
  pairs / Landau-gate attribution, deepKE decomposition + sampling SE.
- Repo files: `scripts/gen_detstagefix_validation.py`,
  `tests/test_gen_detstagefix_validation.py`.

Also found in passing: `drag_migration_log_tier2.md` l. 10700 states
"production runs `cooling_spatial_gate = "none"`" — **false for h405 and all
atlas/finals generators** (`density_scaled`); corrected in place 2026-09-29
(marked `[Corrected 2026-09-29]`).
