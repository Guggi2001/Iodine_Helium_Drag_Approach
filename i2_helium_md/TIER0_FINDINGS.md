# Tier-0 TDDFT Drag Comparison — Findings

**Status:** Settled findings from the §6.4 Tier-0 comparison — what the
deterministic, fixed-`m_eff`, `linear_cubic` drag model reproduces against the
TDDFT references in-window, and why the two cases diverge. This doc carries the
**present-state findings only**; the chronological diagnosis (the withdrawn
"frame systematic" / "windowing + bubble-mode" / "over-damping" readings and why
each was withdrawn) lives in `drag_migration_log.md` → "Withdrawn Tier-0
readings".

Date of record: 2026-06-08. Branch: `drag_implementation`.

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
