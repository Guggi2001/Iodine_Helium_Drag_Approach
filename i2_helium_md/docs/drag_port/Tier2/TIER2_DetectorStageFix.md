# TIER2_DetectorStageFix — retire the E2 relaxation stage?

**Status:** OPEN — Step 0 (zero-MD read) done 2026-09-29; no code changed, no MD
run. Entry document for a fresh-chat re-discussion.

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
   Landau-gated pure cubic — only matters in helium, i.e. essentially for the
   excluded retained class; the Landau-gate question (should it apply in
   Stage I too? v_L is a "P" ledger row) must be decided explicitly;
3. 0.5–8 ns evaporation moves from fixed-dt E2 to exact Gillespie;
4. residual Coulomb after 0.5 ns not propagated (§2c);
5. the retained class is classified at 0.5 ns instead of 8 ns — bound vs
   marginal split (`trap_bound` / `trap_marg`) will shift; both are excluded, so
   the score is unaffected, but the reported split changes.

RNG note: with noise off (T_eff = 0), Stage I and E2 consume the same
per-step draw stream (evaporation + pickup draw every step), so an extended
Stage-I run is draw-aligned with the committed h405 run until the first
physics divergence — a near-CRN comparison.

**Validation proposal:** one MD run = h405 with `ion_simulation_time` 500 ps,
relaxation disabled, detection as-is; compare the full observable vector
(num_scored, trap, supp, nbar_det, n1_solv, w1_solv, midHot, deepKE, chi2)
against the committed h405 row; acceptance = within the h405 seed-to-seed
spread. Config-only (generator constants), but it is a pipeline change →
decide in the re-discussion before running.

## 4. Artifacts

Scratchpad (session-specific, may be cleaned; the reads are fully specified
above and re-derivable from the committed run dir):

`C:\Users\user\AppData\Local\Temp\claude\T--github-synchronized-Iodine-Helium-Drag-Approach-i2-helium-md\65ef168f-b9cf-45ec-9f11-20adc391b9e8\scratchpad\`

- `step0_missed_pickup.py` — §2a (+ anchor)
- `step0b_guard_vs_handover.py` — §2b
- §2c was an inline one-off (formula above).

Also found in passing: `drag_migration_log_tier2.md` l. 10700 states
"production runs `cooling_spatial_gate = "none"`" — **false for h405 and all
atlas/finals generators** (`density_scaled`); corrected in place 2026-09-29
(marked `[Corrected 2026-09-29]`).
