# Tier 2 — Phase D Bridge Findings (Slice Z)

**Date:** 2026-07-03. **Run:** `data/runs/9A_drag_shared_pure_cubic_N50_tier2_bridge_biphasic`
(9 Å / 0.80 eV `biphasic`, N = 50 → 100 ions, 30 ps, Tier-0 dt/seed/drag bundle,
`allow_inconsistent_mass_pairing=True`). **Pinned priored point (plan §0/§4):**
λ₀ = 0.9/ps, f_int = 0.5, f_ret = 0.1, τ = 6.55 ps (the config default — the
geometric mid √(2.6·16.5); the plan's "6.5" is rounded shorthand), κ = 1.0,
picture = `statistical_mixture`, ν = 2.42/ps, s = 3n−3. Numbers below are from
`bridge_summary.txt` in the run dir (report script
`scripts/post_processing/tier2_bridge_report.py`).

**Verdict (editorial, per the plan §7 reporting stance):** the composed
generative driver is **wiring-clean** — every sharp oracle passes — but the
pinned priored point **does not reproduce the anchored 21→19→14 staircase**:
the energy-gated cascade sheds ~0.7 He on average and freezes near n ≈ 20,
not 7 He to n = 14. The miss is *localized* with its lever attached (below):
it sits in the **evaporation-side priors** (RRK kinetics vs. Newton cooling),
not in the mechanism's composition, the drag law, or the pickup channel. This
is exactly the de-risk output the bridge exists to produce before Phase E/F.

---

## 1. Oracle outcomes

| Oracle | Expectation | Result | Outcome |
|---|---|---|---|
| `t×` (sharp, §2.2) | closed form `τ·ln(f_int·0.80/Σ(21))`; all ions agree | closed form 4.951 ps; reconstructed **4.960 ps** (+0.009 = within one stored dt); **all 100 ions agree exactly** (spread 0.0000) | **pass** |
| `t×` bands | [1, 15] ps sanity; GAH25 5–6.5 ps ± factor 2 | inside sanity; 4.96 ps sits just below the literal 5.0 edge (Σ(21) = 0.1878 eV at κ = 1 → 4.95; the plan's own band was 4.8–5.6), well inside ± factor 2 | **pass** |
| `Π(t)` regime (corrected §2.2) | freeze side: Π = 0 at gate-open, Π < 1 throughout, Π → 0 at exit | Π(0) = 0 exactly; max mean Π = **0.0054** (max per-ion 0.027); terminal Π = 0 | **pass** (deep freeze side) |
| 5-term invariant | closes to Verlet drift | max residual 2.23·10⁻⁵ eV = **0.0011 %** of E_system(0) | **pass** |
| `R(t)` / `\|v(t)\|` vs 9 Å TDDFT | qualitative reproduction (Tier-0 quality) | R: RMSE 17.75 Å, mean ratio 1.300; \|v\|: I1 0.876 / I2 1.087 — **indistinguishable from the fixed null** (17.43 Å / 1.299 / 0.868 / 1.073) **and the anchored t★=5.0 run** (17.54 Å / 1.297 / 0.872 / 1.079), same full-window `compare_*` calls | **pass at baseline** (the known Tier-0 9 Å non-radial character; no bridge regression) |
| mean `n(t)` staircase | overlap the anchored 21→19→14 family (loose early-window) | 21.00 → 20.97 (5 ps) → 20.34 (10 ps) → **flat 20.33 to 30 ps**; per-ion terminal range [18, 21]; anchored terminal 14 | **miss** (~0.7 sheds vs 7) |

## 2. The miss, quantified, with its lever (plan §3 framing)

The decline that does occur happens where the mechanism says it should — in a
burst after gate-open (t× ≈ 5 ps → frozen by ~8 ps), driven by the initial
`E_int` cascade, with pickup contributing nothing (Π ≤ 0.005 ≪ 1: no
pickup-sustained shedding, as expected at this condition). The *size* of the
burst is the problem:

- At gate-open, `E_int = Σ(21) = 0.1878 eV` **by construction** (the crossing
  definition) — note this is independent of `f_int`, which only sets *when*
  the gate opens, not the cascade budget.
- The RRK rate there is `k = ν·(1 − D₀(21)/Σ(21))^(s−1) = 2.42·(1 − 0.0060/0.1878)^59
  ≈ 0.36 /ps`, and it collapses as K2 cooling (τ = 6.55 ps) drains `E_int`
  away from the threshold. The time-integral of the rate along the cooling
  trajectory gives **≈ 0.7 expected sheds** — precisely the observed mean
  Δn = 0.67. The run is the mechanism's genuine prediction at this point, not
  a wiring artifact (the 1e-12 isolated-event closures and the exact `t×`
  agreement rule wiring out).
- Energetics are *not* the bottleneck: stripping 21→14 costs
  `Σ(21) − Σ(14) ≈ 0.06 eV`, well inside the 0.188 eV budget. The bottleneck
  is **kinetic** — the s − 1 = 59 RRK exponent keeps `k ≪ ν` near threshold,
  and cooling wins.

**Levers (for Phase F, in priority order):**

1. **κ + electronic picture** — the two genuinely-free co-fit knobs. They set
   `D₀(n)/Σ(n)` (the RRK suppression factor) and the rung sizes; a sharper
   cliff (larger κ) shrinks the near-`n*` rungs and raises `k` steeply.
2. **τ within its [2.6, 16.5] band** — slower cooling holds `E_int` near the
   threshold longer (τ → 16.5 roughly triples the shed integral; τ alone does
   not reach 7 sheds on the anchored timescale, but it moves terminal n).
3. **The `s = 3n−3` dof convention** (Derived) and **ν** (Sourced, pinned) are
   *not* free knobs; if Phase F's staged calibration cannot land the staircase
   inside the κ/picture/τ/f_int bands, that is a genuine mechanism-level
   finding to surface (the RRK dof-counting assumption becomes the suspect —
   an OQ-class item, not a silent retune).
4. **λ₀ is not the lever here** — pickup is structurally dead through the
   decline window (Π ≈ 0), exactly as the corrected §2.2 oracle predicted.

## 3. Interpretation boundaries (plan §8, restated for the record)

- **The bridge is an existence probe, not an arbiter.** This result says the
  *single pinned point* does not land the TDDFT-anchored staircase; it does
  **not** falsify the mechanism (the priored bands are wide and the staircase
  is itself a prior — TDDFT is not ground truth; experiment arbitrates at
  Phase E/F via the I⁺Heₙ size distribution).
- The kinematic observables `R(t)`/`|v(t)|` are **insensitive to the mass
  scenario** at this scale (fixed ≈ anchored ≈ biphasic to three figures) —
  consistent with Tier-1a; the size distribution really is the only
  discriminating observable, which is the Tier-2 premise.
- **Do not carry the freeze-side Π reading to 2.70 eV** (condition-specific).
- A stale-artifact note validating pre-build decision #3: the delivered
  Tier-1a anchored run dirs on this machine carry pre-Phase-B `cfg.json`
  fields (`mass_rate_*`) and no longer load through `RunDirectory` — the
  schedule-family overlay (`build_shell_schedule`, zero artifact dependence)
  was the right call; the anchored comparator above was scored by loading
  `ion.npz` directly (v6→v7 shim).

## 4. Artifacts

- Run dir: `data/runs/9A_drag_shared_pure_cubic_N50_tier2_bridge_biphasic/`
  (`cfg.json`, `neutral.npz`, `ion.npz` v7, `bridge_summary.txt`,
  `figures/bridge_{mean_n_overlay,kinematics,regime_pi}.png`) — machine-local,
  not committed (repo data policy).
- Code: `i2_helium_md/postprocess/bridge_diagnostics.py` (helpers; Phase-E D2
  generalizes; renamed `postprocess/derived_diagnostics.py` at Phase-E E5a),
  `scripts/tier2_common.py` (Phase-F F1 extends),
  `scripts/gen_tier2_bridge_run.py`,
  `scripts/post_processing/tier2_bridge_report.py`.
- Tests: `tests/test_bridge_diagnostics.py` (19; incl. the §2.2 analytic
  oracle and the few-step driver smoke with 5-term closure; renamed
  `tests/test_derived_diagnostics.py` at Phase-E E5a).
- Delivery record: `docs/drag_port/Tier2/drag_migration_log_tier2.md`
  (Phase D — Slice Z DELIVERED entry).
