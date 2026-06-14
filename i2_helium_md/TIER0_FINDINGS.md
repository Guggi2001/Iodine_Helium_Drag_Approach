# Tier-0 TDDFT Drag Comparison — Findings

**Status:** Settled findings from the §6.4 Tier-0 comparison — what the
deterministic, fixed-`m_eff`, `linear_cubic` drag model reproduces against the
TDDFT references in-window, and why the two cases diverge. This doc carries the
**present-state findings only**; the chronological diagnosis (the withdrawn
"frame systematic" / "windowing + bubble-mode" / "over-damping" readings and why
each was withdrawn) lives in `drag_migration_log.md` → "Withdrawn Tier-0
readings". The final section folds in the **drag-form fit RMSE summary**
(formerly `data/reference/drag/shared/trajectory_matching/rmse_summary.md`).

Date of record: 2026-06-08 (Tier-0 findings); fit-summary tables generated
2026-06-14. Branch: `drag_implementation`.

## Setup

- Model: `drag_form=linear_cubic`, `mass_scenario=fixed` at `m_eff=202.954`,
  `noise=none` — the Tier-0 envelope. Coefficients/window from
  `data/reference/drag/<case>/linear_and_cubic/fit_parameters.json`.
- Runs: `single_pulse_N2000{,_18Angst}_drag`, reduced **N=50**, full 20 ps,
  `dt_ion=0.01`, scored in-window only via the `window=` parameter on
  `compare_distance` / `compare_velocity_magnitude`.
- Two comparison modes: **from-onset** (integrate from explosion onset, score
  `[t*, t_end]`; conflates the uncalibrated pre-`t*` transient with the form) and
  **t\*-seeded** (seed the ion at `t*` with the reference state and integrate
  forward; isolates the form). Scoring is also reported **same-smoothed** (MD
  vs the reference run through the extraction's own CEEMDAN+SG denoising) *and*
  raw, never one instead of the other.

## Finding 1 — 18 Å is a clean pass

`linear_cubic` reproduces the 18 Å trace in-window: t\*-seeded distance 0.30 Å,
mean |v| RMSE ~0.09–0.16 Å/ps, ratios ~1.0; from-onset nearly as good. The
extraction window `[4.54, 8.0]` stays in the clean radial regime throughout. The
same-smoothed residual is zero net bias with a pure-shape ~0.09 (the textbook
constant-`m_eff` §2 signature) — the cleaned instrument and the machine are
sound. This is the **committed regression floor**
(`tests/test_tier0_drag_comparison.py`), with named thresholds: distance RMSE
≤ **3.0 Å** (run 2.512), mean(I1,I2) |v| RMSE ≤ **0.25 Å/ps** (run 0.163). N=50
vs N=2000 in-window RMSE is essentially identical, so the reduced-N reference is
a faithful stand-in.

## Finding 2 — the 9 Å reference is genuinely NON-RADIAL (first-class)

Real 3D per-atom velocities **and** positions show atom 2 (the clean extraction
atom) carries a sustained transverse drift of order **~4 Å/ps across the whole
window**, in a **radial↔transverse oscillation** (a y-velocity peak follows each
radial z-peak, repeating ~twice). At `t*` the split is ~99.6 % transverse:

| case | \|v2\| | radial (v·R̂) | transverse | reading |
|---|---|---|---|---|
| 9 Å  | 4.90 | −0.29 | 4.89 | ~99.6 % transverse — atom 2 genuinely **co-translates** |
| 18 Å | 3.12 | −3.12 | 0.17 | ~99.9 % radial — clean |

(`|r1−r2| == R_distance` holds exactly for both cases.)

**The radial-projection convention.** The extraction takes the speed magnitude
`|v2|`, fits a smooth γ to it, and the MD projects that law **radially**. This is
a **defined modelling convention** — "drag depends on total speed, applied along
the radial direction the MD evolves" — *not* a measurement of a purely-radial
drag law. The 9 Å drag law is labelled as such: `|v|`-based, projected radially,
applied to a reference that is genuinely non-radial.

**The residual is a model-dimensionality statement, not a drag error.** The
surviving same-smoothed residual (~0.39 Å/ps t\*-seeded, vs ~0.09 for 18 Å) is the
central-force MD (Coulomb + radial droplet + radial drag) being structurally
unable to carry the real transverse co-translation — the cost of the missing
dimension, recorded and moved past, not a defect to patch by distorting γ. It is
**not** a harness artifact, **not** a drag-form error, and **not** the withdrawn
frame story (see the migration log).

## Finding 3 — correct drag traps the ions; effective binding is calibrated jointly (first-class)

With the in-window-correct drag, the ions arrive at the droplet boundary with
**radial kinetic energy below the static solvation barrier**
`binding_energy_I_ion_eV = 0.308 eV`, so the confining potential pulls them back
— `R(t)` reverses, no ejection. **This is not a bug.** The drag correctly
delivers TDDFT-like (low) surface KE; a static well that deep traps an ion with
that energy. The predecessor's TD-HeDFT analysis is explicit: ~77 % of the
Coulomb energy dissipates within 5 ps, the ions slow to ~4.5 Å/ps, and they
*"have less than the solvation energy ... but are still able to escape ... ion
ejection cannot be predicted with static solvation potential values alone ...
dynamical effects play a major role."* Real ejection is **dynamical** (the He
reorganizes); the MD has no dynamical bubble, so it imposes the full static
barrier the real dynamics bypass. The old hard-sphere model only escaped because
its collisions over-accelerated the ions past the barrier — the same root cause
as the 9 Å over-damping seen from the other side.

**Decision — effective binding depth, calibrated jointly with drag against the
VMI observable:**

- **Do NOT reduce the drag to force escape** as a separate after-the-fact knob —
  that detunes the anchored in-window velocity to mask the binding (the canonical
  Method-B trap, made concrete). Drag and binding move *together* against the data.
- **`binding_energy_I_ion_eV` becomes an effective, calibrated parameter** — a
  pragmatic stand-in for absent dynamical-barrier physics, *not* a re-measured
  solvation energy (the static 0.308 eV is the upper-bound starting point).
- **Joint extraction over the full window `[2.67, 14 ps]`** against the VMI
  final-velocity distribution (target), with the TDDFT escape energy as a sanity
  cross-check (ions escape with < 0.308 eV). The coupling is only visible through
  the forward-integrated trajectory to ejection, which is why it surfaced here.
- **Degeneracy caveat:** drag + binding are under-determined on a single
  trajectory (stronger drag + shallower binding mimics the reverse for in-window
  velocity); the drag is anchored by the in-window velocity *shape* and the
  **held-out VMI distribution** breaks the degeneracy.
- **Coupled pair:** the effective binding is **stamped alongside the drag
  coefficients**; the §6.5 consistency machinery extends to refuse an unvalidated
  drag↔binding pairing (`DRAG_PORT_DESIGN_DECISIONS.md` §6.5.1).
- **Dynamical-barrier structure** (a surface-weakened or velocity-dependent
  depth) is the principled long-term fix, deferred until the effective-static
  depth is tested against VMI.

This **vindicates Method B's mandatory held-out validation**: a trajectory-matched
drag can reproduce the in-window velocity beautifully and *still* trap the ions /
fail the VMI observable. In-window match ≠ correct production behaviour.

## Consequences — extraction method and Tier-0 role

- **Extraction A → B.** The non-radial reality plus the hand-tuning observation
  motivate **trajectory-matching calibration**
  (`METHOD_B_trajectory_matching_extraction.md`): fit the coefficients (jointly
  with the effective binding) by minimizing the forward-integrated trajectory
  RMSE against the same-smoothed reference, over the full post-dynamic-start
  window `[2.67, 14 ps]` (final velocity is the production-relevant quantity).
- **Accepted 9 Å cost.** The full window re-includes the post-6 ps non-radial
  region, so the 9 Å coefficients can absorb some unrepresentable transverse
  drift — a risk **deliberately accepted** (no truncation both reaches final
  velocity and excludes the non-radial region). 18 Å is clean-radial throughout.
- **Tier 0 repurposed: consistency → held-out generalization.** Under B the
  trajectory match *is* the fit objective, so the consistency-check framing is
  circular and retired. The infrastructure (the `window=` parameter, the
  harnesses, the gate) survives but now scores **held-out** data: the cross-case
  shared-form check (the transport-physics signal) and the downstream VMI
  observable. This held-out role is mandatory for 9 Å specifically — it is the
  only thing standing between the accepted transverse contamination and a
  corrupted coefficient set. See METHOD_B §3.5, §4, §5.
- `EXTRACTION_FRAME_FIX_milestone.md` is demoted to a contingency (the He-field
  relative-velocity route, revived only if the radial-projection convention
  cannot be made to generalize).

## Infrastructure delivered (independent of the physics outcome)

- `compare_trajectories.py`: additive `window=(t_start, t_end)` parameter
  (`window=None` bit-identical) + windowing tests.
- `scripts/post_processing/tier0_drag_comparison.py` (from-onset),
  `tier0_tstar_seeded_comparison.py` (t\*-seeded, real `|v|`-radial seed + real
  radial/transverse split diagnostic + `|r1−r2| == R_distance` cross-check),
  `tier0_same_smoothed_comparison.py` (`compare_speed_to_reference` +
  `load_smoothed_speed_reference`).
- `scripts/gen_tier0_runs.py`: reproducible run generation (`data/runs/`
  gitignored).
- `data/reference/drag/{9A,18A}/tier0/md_mean_trajectory_N50.csv`: committed
  tiny ensemble-mean references.
- `tests/test_tier0_drag_comparison.py`: the 18 Å-anchored regression gate +
  named thresholds (9 Å recorded finite-only).

## Drag-form fit results — RMSE summary (Method-B §9 + §10)

One-glance comparison of every trajectory-matching fit: the shared-form joint
fits and Stage-1 analogs, plus the §10.8 per-case single-curve diagnostic fits.
Pulled from the `fit_parameters.json` / `stage1_*.json` artifacts in
`data/reference/drag/shared/trajectory_matching/` and
`data/reference/drag/<case>/trajectory_matching/<variant>/`. **Generated
2026-06-14.**

- RMSE / objective in **Å/ps**; objective = equal-weight mean of the two
  per-case in-window `|v2|` RMSEs (+ escape penalty, 0 for every row — all fits
  fully escape). **Escape = yes** ⇒ `escape_fraction = 1.0`.
- Parameter units: `a` [amu/ps], `b` [amu·ps/Å²], `c` [amu/Å],
  `C` [amu·Å^(1−n)·ps^(n−2)], `n` dimensionless, `E_bind` [eV].
- Gating bands (§9.4, reused by §10): 18 Å RMSE ≤ 0.19, 9 Å RMSE ≤ 0.45,
  escape = 1.0 both cases. Incumbent (`shared_pure_cubic`) objective =
  **0.1292649398514104**; `Δ = variant − incumbent`.

### Stage-2 — shared joint fits (both cases in the objective; gating)

| Variant (try) | Form | Parameters | E_bind | 18 Å RMSE | 9 Å RMSE | Escape (18/9) | Objective | §9.4 | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| `shared_pure_cubic` *(§9, PRODUCTION)* | linear_cubic | a = 0, b = 2.5154 | 0.1168 | 0.1344 | 0.1241 | yes / yes | **0.129265** | PASS | incumbent / production candidate |
| `shared_3param` *(§9)* | linear_cubic | a = 0.00025, b = 2.5159 | 0.1168 | 0.1345 | 0.1240 | yes / yes | 0.129265 | PASS | `a → 0`; `T_a0`-equivalent to incumbent |
| `pl_shared_3param` *(§10)* | power_law | C = 2.8352, n = 2.9269 | 0.1126 | 0.1305 | 0.1278 | yes / yes | **0.129171** | PASS | **EQUIVALENT** (Δ = −0.0001); n̂ ≈ 3 confirms pure-cubic |
| `lq_shared_3param` *(§10)* | linear_quadratic | a = 0.0001, c = 12.7922 | 0.0482 | 0.1308 | 0.1955 | yes / yes | 0.163152 | PASS | **WORSE** (Δ = +0.0339); rejected-by-objective |
| `lq_shared_pure_quadratic` *(§10)* | linear_quadratic | a = 0, c = 12.7973 | 0.0483 | 0.1313 | 0.1950 | yes / yes | 0.163150 | PASS | **WORSE** (Δ = +0.0339); rejected-by-objective |

Every variant passes the §9.4 bands; the **objective** (not the bands) separates
the forms.

### Stage-1 analog — 18 Å-only fit, predict held-out 9 Å (recorded, NON-gating)

Each family is fit on 18 Å alone; the untouched 9 Å trajectory is scored as a
prediction (band ≤ 0.45 Å/ps).

| Family | Fit parameters (18 Å only) | E_bind | 9 Å pred. RMSE | Escape | Band ≤ 0.45 |
|---|---|---|---|---|---|
| pure-cubic *(§9, reused §8 18 Å bundle)* | a = 1.456, b = 3.316 | 0.0709 | **0.2685** | yes | PASS |
| `power_law` *(§10)* | C = 1.0529, n = 4.000 | 0.0825 | 0.4110 | yes | PASS |
| `linear_quadratic` *(§10)* | a = 16.6254, c = 6.1313 | 0.0532 | 0.6987 | yes | **FAIL** |

**Form discrimination.** `power_law` (free n) recovers the pure-cubic exponent:
joint fit lands at **n̂ = 2.927 ≈ 3** with the lowest objective (0.129171, inside
the ±0.005 equivalence band) and a passing held-out 9 Å prediction — *not* the
Method-A power-law exponent (n ≈ 2.06). Forced `v²` (`linear_quadratic`) is
measurably worse: both lq variants collapse to the pure-quadratic corner
(`a → 0`), sit at objective 0.1632 (Δ = +0.034, outside the band), and the lq
9 Å prediction (0.699) fails. **Conclusion:** the trajectory objective
discriminates the exponent and picks **n = 3** → `shared_pure_cubic` confirmed,
no escalation, presets unchanged. (Full record: METHOD_B §10.7;
`drag_migration_log.md`.)

### Per-case single-curve fits (§8 + §10.8 — DIAGNOSTIC, never preset-wired)

Each row fits **one curve only** — the per-case diagonal across all three forms:
`linear_cubic` (§8) at `.../trajectory_matching/fit_parameters.json`, and the
`linear_quadratic` / `power_law` variants (§10.8) at
`.../trajectory_matching/<variant>/fit_parameters.json`. These are
**calibrated-not-validated**: each fit trivially matches its own curve, so the
RMSE is the *self-fit* (**not** held-out, not comparable to Stage-1 predictions).
**9 Å is the genuinely non-radial / transverse-contaminated case.** The
production law is the §9 *shared* `shared_pure_cubic` joint fit — **no per-case
row is a production candidate**.

Half-width conventions differ by provenance: `linear_cubic` (§8) bands are
**seed-sweep std + RMSE sensitivity**; §10.8 form bands are **sensitivity-only**.
For `power_law`, `C_err` is the fixed-`n` partial band `C_err = γ_ref_err /
v_ref**(n−1)` (not a marginal uncertainty).

| Case | Variant | Form | Parameters (± half-width) | E_bind | Self RMSE | Escape | Diagnostic note |
|---|---|---|---|---|---|---|---|
| 18 Å | `linear_cubic` *(§8)* | linear_cubic | a = 1.46 ± 1.46, b = 3.32 ± 0.16 | 0.0709 ± 0.0068 | 0.0968 | yes | `a_err = a` exactly → `a` consistent with 0: §8 weak-`a` ridge that motivated the shared refit |
| 18 Å | `lq_shared_3param` *(§10.8)* | linear_quadratic | a = 16.63 ± 1.59, c = 6.13 ± 0.59 | 0.0532 ± 0.0071 | 0.0994 | yes | `a` *not* collapsed — but degenerate with the `a=0` corner (next row, equal RMSE): §10.4.1 weak-`a` ridge |
| 18 Å | `lq_shared_pure_quadratic` *(§10.8)* | linear_quadratic | a ≡ 0, c = 11.27 ± 0.55 | 0.0598 ± 0.0080 | 0.0993 | yes | same RMSE as 3-param → `a` carries no in-curve information |
| 18 Å | `pl_shared_3param` *(§10.8)* | power_law | C = 1.05 ± 0.10, **n = 4.000 ± 1.333** | 0.0825 ± 0.0079 | 0.0929 | yes | **exponent unidentified** — railed to n = 4 bound; 18 Å's narrow speed range has no exponent leverage |
| 9 Å | `linear_cubic` *(§8)* | linear_cubic | a = 7.50 ± 0.72, b = 1.75 ± 0.03 | 0.1543 ± 0.0031 | 0.0409 | yes | `a` resolved away from 0 here (unlike 18 Å) — but jointly the §9 shared fit still drives `a → 0` |
| 9 Å | `lq_shared_3param` *(§10.8)* | linear_quadratic | **a → 0**, c = 10.32 ± 0.10 | 0.1303 ± 0.0044 | 0.0446 | yes | **collapses to pure-quadratic** independently |
| 9 Å | `lq_shared_pure_quadratic` *(§10.8)* | linear_quadratic | a ≡ 0, c = 10.32 ± 0.10 | 0.1303 ± 0.0064 | 0.0446 | yes | identical to the 3-param fit (already at `a = 0`) |
| 9 Å | `pl_shared_3param` *(§10.8)* | power_law | C = 3.59 ± 0.04, **n = 2.651 ± 0.026** | 0.1566 ± 0.0031 | 0.0406 | yes | **exponent identified tightly** — 9 Å's broad speed range pins `n` within the single curve |

*The 18 Å rows are the same fits as the 18 Å-only Stage-1 analogs (the §8
`linear_cubic` row = the "pure-cubic, reused §8 18 Å bundle" Stage-1 row; the two
§10.8 rows = the lq / pl Stage-1 analogs). The RMSE here is their 18 Å self-fit;
Stage-1 reports their held-out 9 Å prediction instead.*

**Reading the per-case table.** Per-case exponent leverage is
**case-asymmetric**, refining the §10.4.1 "exponent unidentified within a single
curve" expectation: 9 Å alone pins **n̂ = 2.651 ± 0.026** and collapses `a → 0`,
while 18 Å alone **rails to n = 4** (half-width 1.333) and its `a` is the flat
weak-`a` ridge. The shared joint fit's **n̂ ≈ 2.927 ≈ 3** is the *cross-case*
combination of the two speed scales — neither single case lands at 3 alone (9 Å
biased low and contaminated; 18 Å unidentified). The §8 `linear_cubic` rows show
the same per-case-vs-cross-case pattern in the *linear* term: 18 Å's `a` is
consistent with 0 (`a_err = a`) while 9 Å's is resolved (`a = 7.50`), yet the §9
shared fit collapses `a → 0` jointly. Diagnostic only; the incumbent is
unaffected. (Full record: METHOD_B §8 / §10.8; `drag_migration_log.md`.)
