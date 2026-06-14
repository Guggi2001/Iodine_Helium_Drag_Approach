# Method B — Trajectory-Matching Drag Extraction

**Status:** Implemented and run; this is the **active extraction method**. It
fits the drag coefficients (jointly with the effective droplet binding) by
minimizing the *forward-integrated trajectory* RMSE against the smoothed
reference — as opposed to Method A, which regresses the force-balance residual
`F_drag` against speed directly. §1–§7 are the live method specification;
§8–§10.7 record the outcomes (per-case fit → shared-form joint refit →
alternative-form discrimination). The production law is `shared_pure_cubic`
(`γ = g·b·v²`); full decision/delivery history is in `drag_migration_log.md`.

**Why it exists.** The Tier-0 investigation showed that hand-adjusting the fit
parameters while comparing the forward-integrated MD against the smoothed
reference can drive the in-window difference very small. That is end-to-end
calibration: fit the parameters so the *downstream observable you care about*
(the in-MD trajectory) matches, rather than fitting an intermediate quantity
(`F_drag` vs `v`). It is a legitimate and arguably superior extraction method —
**but it changes what the validation must be**, which is the load-bearing point
of this document (§4).

---

## 1. Method A vs. Method B

**Method A — direct force-balance regression (current pipeline,
`Drag_extraction_code.md`).** Compute `F_drag(t) = m_eff·a(t) − F_C(t)`
point-by-point, then regress `|F_drag|` against the speed `v` (log-linear for
power-law, `curve_fit` for linear+cubic). The fit never sees the
forward-integrated trajectory; it fits an *intermediate* quantity. Tier 0 then
*independently* asked "does forward-integrating this γ reproduce the reference?"
— a real question with a possible "no."

**Method B — trajectory-matching calibration (this doc).** Choose `{a,b}` to
minimize the **forward-integrated, radial-projected, in-window trajectory RMSE**
against the smoothed reference:
$$\{a,b\}^\star = \arg\min_{a,b}\; \text{RMSE}_{[t^\star,\,t_\text{end}]}
\big(\,|v_2|^\text{MD}(a,b)\;,\; |v_2|^\text{smoothed-ref}\,\big),$$
where `|v2|^MD(a,b)` is produced by the *actual BAOAB ion driver* run with those
coefficients (the same forward integration Tier 0 used). The fit objective **is**
the trajectory match.

The key structural difference: **under B, the trajectory match is the fit
objective, so it can no longer also be the validation.** This retires Tier 0's
*consistency-check* role (§3) and forces validation to move to *held-out* data
(§4).

---

## 2. Why B is attractive (and legitimate)

- **It optimizes the observable that matters.** The drag law's job is to make the
  MD trajectory (and downstream, the VMI distribution) match reality. B fits
  exactly that, rather than an intermediate `F_drag(v)` regression whose good fit
  does not guarantee a good *trajectory* (the integration can compound small
  per-step force errors).
- **It absorbs the integration faithfully.** A is blind to how the BAOAB O-step,
  the spatial gate, and the per-step mass treatment transform `γ` into a
  trajectory; B sees the whole pipeline and fits through it.
- **It matches the project's own finding.** The hand-tuning observation that
  motivated this — small parameter changes producing near-perfect in-window
  agreement — is exactly B done manually. Formalizing it removes the "by hand"
  unrigour (no provenance, no uncertainty) and replaces it with a defined
  objective and the existing seed-sweep uncertainty machinery.

---

## 3. What B retires — Tier 0's consistency-check role

Under A, Tier 0 = "forward-integrate the independently-extracted γ; does it
reproduce the reference?" Under B, the parameters were *chosen* to minimize
exactly that trajectory RMSE, so re-running it reads back the fit objective:
**"did I tune the parameters so the answer matches? — yes."** That is training on
the test set; it proves nothing and must not be committed as a "pass."

**So the consistency-check framing of Tier 0 is obsolete under B — correctly
retired, not patched.** What is *not* obsolete is the broader role Tier 0's
infrastructure can play: checking **generalization** beyond the fit (§4). The
distinction is the entire methodological crux:

> **Litmus test.** Is there any data the parameters were **not** fit against, on
> which the agreement is then checked? **Yes** → calibration + validation
> (legitimate, even superior). **No** — every point used to judge agreement was
> also used to choose the parameters → circular, proves nothing.

---

## 3.5 The calibration window — full post-dynamic-start window `[2.67, 14 ps]` (2026-06-09)

**Decision: calibrate against the full window from the dynamic start to the end
of the dynamics, `[t* ≈ 2.67, ~14 ps]`, for BOTH cases**, scored against the
extended smoothed reference `cleaned_data_long.csv` (the **same** CEEMDAN+SG
pipeline as before, extended in span to 14 ps — same anti-laundering discipline).

**Rationale — final velocity is the production-relevant quantity.** The drag
law's job is the VMI **final-velocity** distribution. A window truncated at ~6 ps
fits the early/mid trajectory but says nothing about whether the ion reaches the
*correct terminal speed* — and terminal speed is exactly what the VMI observable
measures and what the drag+binding pair must reproduce (§5.5). Calibrating on the
full window puts the final velocity inside the fit objective.

**Cost, recorded as an accepted risk (not solved) — the 9 Å transverse
re-entry.** The ~6 ps truncation existed because 9 Å's atom 2 develops a
directional **transverse drift after ~6 ps** that a central-force MD structurally
**cannot** reproduce (the §5 non-radial finding). Extending the window to 14 ps
**re-includes that region in the 9 Å fit.** Under B's optimizer this means the
9 Å `{a,b}` (and the coupled effective binding, §5.5) will **absorb some
unrepresentable transverse drift** — a fit that is numerically good and partly
physically corrupted. **This risk is deliberately accepted** (9 Å makes
transverse motion regardless; there is no truncation that both captures final
velocity *and* excludes the non-radial region). It is **accepted, not
eliminated** — and it changes the validation (§4):

- The **held-out-window axis is forfeited** — the whole window is now fit, so
  there is no in-window sub-region left to validate on. Accepting the full-window
  risk therefore *raises* the stakes on the remaining held-out axes, it does not
  lower them.
- The 9 Å coefficients carry an expanded standing flag:
  `|v|`-projected-radially, **full-window**, **transverse-contaminated** —
  trustworthy *only* insofar as they survive the cross-case and VMI held-out
  checks (§4).

18 Å is unaffected by the cost: it is clean-radial across the whole `[4.54, 8.0]`
(now extended) window, so the full-window fit for 18 Å is purely beneficial
(captures final velocity, no transverse contamination).

---

## 4. What B requires — held-out validation (mandatory, not ceremony)

B is *training*. It needs a *test set* the training never saw. With the
full-window calibration (§3.5) the **held-out-window axis is gone**, so the
remaining two axes are now *load-bearing*, not optional extras:

1. **Held-out window — FORFEITED under full-window calibration (§3.5).** Retained
   here only as the record of why it no longer applies: the calibration now spans
   `[t*, 14 ps]`, leaving no untouched in-window sub-region. (It would return only
   if a future calibration truncated the window again.)
2. **Held-out case — the shared-form / cross-case test (now the PRIMARY held-out
   axis; strongest "is this real transport physics" signal, §3.6).** A form fit to
   one case must also fit the *other* case with its own coefficients, **without**
   independently trajectory-tuning both. If you trajectory-match 9 Å *and* 18 Å
   separately, each matches trivially and the cross-case signal is destroyed. The
   honest version: derive the form (and ideally the coefficient relationship)
   constrained so it generalizes across both cases, or fit one and check it
   predicts the other. **This is the main thing standing between the accepted 9 Å
   transverse contamination (§3.5) and a corrupted coefficient set** — a 9 Å fit
   that absorbed transverse drift will not match the form that fits clean-radial
   18 Å.
3. **Held-out observable — the downstream tiers (the ultimate arbiter).** The VMI
   final-velocity distribution (Tier 2/3) and ensemble second moments are
   *different observables entirely* from the calibration trajectory, so they are
   inherently held out from a trajectory fit. With the full-window choice fitting
   final velocity on the *single* calibration trajectory, the VMI **distribution**
   is the test that the fit generalizes across the ensemble. These are the
   ultimate test that B did not overfit — and the only check that can catch a 9 Å
   fit whose transverse contamination still happens to pass cross-case.

**Without at least one held-out axis, a B fit is unfalsifiable.** With the
held-out-window axis forfeited (§3.5), **at least one of cross-case / VMI is now
mandatory, not a nice-to-have.** The existing Tier-0 infrastructure (`window=`
parameter, the harnesses, the gate machinery) is **reused** to score these
held-out comparisons — it survives; only the *data it scores* changes from in-fit
to held-out (§3 of the repurposed Tier-0 docs).

---

## 5. The 9 Å danger — fitting a radial model to non-radial data

This is the warning that makes held-out validation **non-optional specifically
for 9 Å.** The 9 Å reference is genuinely **non-radial**: a transverse (y) drift
of order ~4 Å/ps, in a radial↔transverse oscillation (a y-velocity peak follows
each radial (z) peak, repeating ~twice across the window). The MD is a
**central-force** model (Coulomb + radial droplet + radial drag) and evolves
**only radial dynamics.**

Trajectory-matching `{a,b}` to force the *radial-projected* MD onto a
*non-radial* reference means **the drag coefficients absorb the transverse
discrepancy** — `a,b` become a fudge factor compensating for the dimension the MD
lacks. B will do this silently and report an excellent fit. The resulting law:

- is **physically corrupted** (it encodes "whatever radial drag makes the radial
  projection match a non-radial trajectory"), and
- will be **wrong on any trajectory with a different radial/transverse balance** —
  i.e. the production ensemble, and 18 Å.

**B cannot detect this on its own** — the fit looks perfect *because* the
coefficients ate the error. Only the held-out checks catch it:

- 9 Å's B-fit coefficients, if corrupted, will **fail the held-out-case check**
  (won't match the form that fits genuinely-radial 18 Å) and likely **fail the
  VMI observable.**
- 18 Å, being genuinely radial, can be trajectory-matched **safely** — its
  coefficients are not absorbing a missing dimension.

**Consequence for method policy:** 18 Å is the trustworthy B-calibration case;
9 Å's B fit must be treated as *suspect until it passes a held-out check*, and the
9 Å drag law carries an explicit flag that it is a `|v|`-projected-radially
convention applied to a non-radial reference (see `TIER0_FINDINGS.md` → 9 Å
non-radial finding).

---

## 5.5 Drag and droplet binding are a jointly-calibrated coupled pair

A second, independent reason the **held-out observable (VMI) is the real test** —
discovered when the validated drag was first run to ejection
(`TIER0_FINDINGS.md` → "Correct drag traps the ions"):

A trajectory-matched (in-window-correct) drag delivers TDDFT-like *low* kinetic
energy at the droplet surface. The MD's droplet potential is a **static** well of
depth `binding_energy_I_ion_eV = 0.308 eV` (the computed solvation energy), and a
sub-barrier ion is **trapped** — `R(t)` reverses, the ion never ejects. Yet the
TD-HeDFT ions *do* escape with less than 0.308 eV, because real ejection is
**dynamical** (the He reorganizes; the static barrier is bypassed). The old
hard-sphere model only escaped because it over-accelerated the ions past the
static barrier by brute force; the correct drag removes that excess and exposes
that **a static barrier this deep is incompatible with correct surface kinetic
energies.**

**Therefore the binding depth is not a fixed constant but a calibrated parameter
coupled to the drag:**

- **Do not reduce the drag to force escape** — that detunes the validated
  quantity to mask a different component's error (the canonical trap this whole
  document warns about, now concrete).
- **`binding_energy_I_ion_eV` is an effective parameter**, calibrated **jointly
  with the drag** (not drag-fixed-then-binding-patched) over the full window
  `[2.67, 14 ps]` against the **VMI final-velocity distribution** (target (b)),
  with the TDDFT escape energy as a sanity cross-check (c). It is a pragmatic
  stand-in for absent dynamical-barrier physics — **not** a re-measured solvation
  energy (the static 0.308 eV is the upper-bound starting point). *Degeneracy
  caveat:* drag and binding can trade off for the in-window velocity, so the joint
  fit is under-determined on the calibration trajectory alone; the drag is
  anchored by the in-window velocity shape and the **held-out VMI distribution**
  breaks the degeneracy among near-equivalent pairs.
- **The effective binding is STAMPED alongside the drag coefficients** in the
  same bundle/`fit_parameters.json` that carries `extraction_mass_amu` etc.
  (field e.g. `effective_binding_energy_I_ion_eV` + its calibration provenance).
  A drag bundle thus records the binding it was validated with; the §6.5
  config-load consistency guard extends to **refuse a drag↔binding pairing that
  was not jointly validated** (parallel to the mass↔coefficient guard). Swapping
  drag coefficients now *requires* re-checking the binding-permits-escape — a
  detectable inconsistency, not a silent one.

**Why this belongs in the Method-B document:** the binding↔drag coupling is only
visible *through the forward-integrated trajectory to ejection* — exactly the
end-to-end view B adopts. A Method-A `F_drag`-vs-`v` fit never sees the droplet
well at all, so it cannot surface this. And the failure mode is the canonical
B-overfit signature: an in-window-perfect drag that fails the held-out VMI
observable because it was calibrated against the wrong target. The dynamical-
barrier structure (a surface-weakened or velocity-dependent depth) is the
principled long-term alternative, deferred until the effective-static depth is
tested against VMI (`TIER0_FINDINGS.md`).

---

## 6. Procedure (what a B extraction must record)

To have the same rigour as Method A, a B extraction produces:

- **Objective:** full-window forward-integrated trajectory RMSE against the
  *same-smoothed* reference `cleaned_data_long.csv` (the extraction's own
  CEEMDAN+SG extended to 14 ps, not a hand-tuned filter — same anti-laundering
  discipline as the same-smoothed comparison).
- **Fit window** `[t* ≈ 2.67, ~14 ps]`, full post-dynamic-start window for both
  cases (§3.5), recorded — including the accepted 9 Å transverse-contamination
  risk and the forfeited held-out-window axis.
- **Coefficients `{a,b}`** with an **uncertainty band** from the existing
  seed-sweep machinery applied to the *trajectory-RMSE* objective (not the
  `F_drag`-vs-`v` objective).
- **Provenance stamp** carried into `fit_parameters.json` exactly as Method A's
  bundle is: `extraction_method = trajectory_matching`, plus `m_eff`,
  `extraction_mass_model`, window, date/branch. The Slice-3 loader's
  provenance/consistency guards apply unchanged — B coefficients are a coefficient
  swap behind the same interchangeable surface; `drag.py`/`baoab.py` consume them
  unchanged.
- **Effective droplet binding stamped alongside (§5.5):**
  `effective_binding_energy_I_ion_eV` + its calibration provenance, recording the
  binding the drag was jointly validated with against VMI. The §6.5 guard refuses
  a drag↔binding pairing that was not jointly validated.
- **Held-out validation result** (§4) recorded alongside, with the 9 Å flag (§5).

---

## 7. Relationship to the other docs

- **Supersedes** the Method-A *default* for extraction going forward, **but does
  not delete it** — A remains the independent cross-reference (an A-fit and a
  B-fit that agree is itself reassuring; divergence is informative).
- **Retires** the consistency-check framing of Tier 0 (§3); the Tier-0 docs are
  rewritten to the **held-out generalization** role (§4), keeping the
  infrastructure.
- **Demotes further** `EXTRACTION_FRAME_FIX_milestone.md` — the 9 Å non-radial
  reality is now handled by (a) the explicit radial-projection convention flag and
  (b) the held-out validation that guards against the dimensionality fudge, not by
  a relative-velocity re-extraction. The He-field relative-velocity route remains a
  *contingency* only if held-out validation shows the radial-projection convention
  cannot be made to generalize.
- **Feeds** Tiers 1–3: B produces the coefficients; the tiers (now including the
  repurposed held-out Tier 0) validate generalization.

---

## 8. Per-case Method-B fit — outcome (superseded by §9; delivery record in `drag_migration_log.md`)

B was implemented as a named extraction method
(`i2_helium_md/extraction/trajectory_matching.py`,
`scripts/extraction/method_b_extraction.py`) producing `{a, b}` jointly with
`E_bind` (§5.5) by minimizing the forward-integrated trajectory RMSE against the
same-smoothed reference. The cross-case held-out axis ran; VMI stamped `pending`.

**Result: both per-case fits are in-window excellent but FAIL the cross-case
held-out bands** — and stand flagged **not usable**, not wired into presets:

- **18 Å:** in-window 0.097 Å/ps, full escape, `E_bind` 0.071 eV.
- **9 Å:** in-window 0.041 Å/ps (*suspiciously* better than clean-radial 18 Å —
  the §5 signature), full escape, `E_bind` 0.154 eV.

The dominant finding is that **the linear coefficient `a` is weakly identified by
the full-window trajectory objective** (the cubic term dominates over both
windows): 18 Å's `a` pinned at the optimizer bound, 9 Å showed a live
`a`↔`E_bind` degeneracy ridge, and the cross-case `E_bind` disagreement (0.083 eV)
is entangled with it. This motivates the cross-case **shared-form joint refit**
(§9), which resolves the ridge.

---

## 9. Shared-form joint refit — spec (resolves the §8 weak-`a` ridge)

The chosen response to §8 is the **cross-case shared-form joint refit**: fit one
size-independent drag law plus one effective binding across both cases. A
successful shared fit is the strongest available statement that the drag law is a
*general working principle* (one size-independent transport law), not a per-case
fudge. The pure-cubic reduced form (`a ≡ 0`) is folded in as an explicit variant.

### 9.1 The VMI axis is unreachable until Tier 1 — gate policy

The VMI channels are **mass-selected** and the final velocity depends on the
mass history, so a fixed-`m_eff` ensemble cannot be honestly scored against
`vmi_iplus_he.csv`. The designated ultimate arbiter (§4 axis 3) therefore
requires Tier 1 mass evolution. Gate policy: if the shared-form refit passes its
pre-registered bands (§9.4), **Tier 1 ungates** and the **VMI tier becomes the
post-Tier-1 final arbiter** (stamped `pending` until then). Cross-case is the
only falsification axis reachable now, so the design preserves held-out content
(Stage 1) *before* the joint fit consumes the axis (Stage 2).

### 9.2 Hypothesis and stages

**Hypothesis:** one size-independent drag law `{a, b}` plus one effective binding
`E_bind` carries both cases (the bulk He density the ion traverses is the same
liquid in both droplets; the spatial gate handles the surface). Fully shared =
**3 parameters against 2 full trajectories — over-constrained, hence
falsifiable**, restoring the falsifiability the per-case (6-parameter) fits lack.

- **Stage 1 (held-out, non-gating)** — score the untouched 9 Å prediction from
  the §8 18 Å single-case bundle (`a = 1.456`, `b = 3.316`, `E_bind = 0.0709`).
  Pure prediction: one 9 Å forward integration + scoring, no optimizer, no
  variants. 9 Å is strictly held-out (zero refit). If the law generalizes, the
  residual should land near the Tier-0 ~0.39 Å/ps dimensionality residual.
- **Stage 2 (gating)** — joint fit over both cases. Objective: per-case
  (RMSE + escape penalty) first, then the **equal-weight mean** (prevents
  long-window dominance). Variants:
  - `shared_3param` — `a` free, lower bound 0 (§9.5);
  - `shared_pure_cubic` — `a ≡ 0`, fit `{b, E_bind}`. If the two variants are
    indistinguishable (Δ ≤ `T_a0`, §9.4), the pure-cubic reduced form is the
    **empirical conclusion** on `a`'s identifiability;
  - `diagnostic_4param` — shared `{a, b}`, per-case `E_bind`; run **only if the
    primary fails its bands** (localizes drag-vs-binding failure).
- **Stage-1↔Stage-2 parameter-shift check** — qualitative, recorded-only (no
  band): a large shift when 9 Å enters the objective would flag 9 Å dragging the
  law toward its transverse contamination (the §5 warning, made measurable).

### 9.3 Normalization and anchors

Optimization uses the extractor's normalized bounded Nelder-Mead pattern with
multi-start (E_bind pre-scan + ridge-probing starts); the conditioning anchors
`{a0, b0}` come from the **18 Å Method-A bundle only** (a convention choice, kept
for the record).

### 9.4 Pre-registered thresholds (first-runs rule — fixed before any run, not re-tuned after)

| name | check | band |
|---|---|---|
| `S1_PRED_RMSE_MAX` | Stage-1 9 Å prediction same-smoothed \|v2\| RMSE, full 9 Å window | ≤ 0.45 Å/ps |
| `S1_PRED_ESCAPE` | Stage-1 9 Å prediction escape fraction at the 18 Å-fit `E_bind` | = 1.0 |
| `S2_RMSE_18A_MAX` | Stage-2 joint per-case 18 Å RMSE | ≤ 0.19 Å/ps |
| `S2_RMSE_9A_MAX` | Stage-2 joint per-case 9 Å RMSE | ≤ 0.45 Å/ps |
| `S2_ESCAPE` | Stage-2 escape fraction, both cases | = 1.0 |
| `T_a0` | pure-cubic equivalence: Δobjective = objective(`a≡0`) − objective(`a`-free) | ≤ 0.005 Å/ps |

Anchoring (recorded so a later re-derivation can re-judge the same numbers):
0.45 = the Tier-0 9 Å dimensionality residual (~0.39 Å/ps) + margin — the
residual expected *if the law generalizes and only the known transverse
deficit remains*; 0.19 = 2× the 18 Å single-case B-fit RMSE (0.097); the 9 Å
band is deliberately **not** 2× its single-case 0.041 (that value is the §5
suspiciously-good signature, not a clean anchor); 0.005 Å/ps ≈ the optimizer
`fatol` scale. The `a`-free optimum can never be worse than the nested
`a ≡ 0` optimum (up to optimizer noise), so the Δ above is ≥ 0 by
construction.

**Gating policy:** only the Stage-2 bands (`S2_RMSE_18A_MAX`, `S2_RMSE_9A_MAX`,
`S2_ESCAPE`) carry verdict power; the Stage-1 rows are recorded for the held-out
audit trail but are non-gating. `T_a0` is a classification threshold (which
conclusion to record on `a`'s identifiability), not a pass/fail gate. These §9.4
bands supersede the §8 provisional cross-case bands (the 18 Å A↔B `a`-ratio check
is dead under the weak-`a` finding and not re-applied).

**Verdict mapping.** Stage-2 primary passes its bands → the shared bundle is the
production candidate, Tier 1 ungates (§9.1), VMI = post-Tier-1 final arbiter.
Primary fails → run `diagnostic_4param` to localize; per-case bundles stay
flagged not-yet-usable; `EXTRACTION_FRAME_FIX_milestone.md` remains the
contingency.

### 9.5 `a = 0` and the §3.3 guard

`a = 0` with `b > 0` is still strictly dissipative
(`γ = g·(a + b·v²) ≥ 0`, vanishing only at `v = 0` where no energy can be
added). Decision: realize the pure-cubic variant by **relaxing the §3.3
guard from `a > 0` to `a ≥ 0` when `b > 0`** for `linear_cubic`, rather than
adding a `pure_cubic` form tag — no loader/enum churn; the variant remains
`linear_cubic` with `a = 0`. (The new-form-tag alternative is noted here in
case a later phase wants the explicit tag.)

### 9.6 Artifacts and provenance (outcome in §9.7)

`data/reference/drag/shared/trajectory_matching/`: `fit_parameters.json`
(`extraction_method = "trajectory_matching"` — unchanged loader contract — plus
`calibration_cases`, `stage`, `variant` provenance fields),
`stage1_prediction.json`, and the verdict records; the per-case
`held_out_validation.json` files updated with the shared-refit verdicts under the
§9.4 bands. Both variant bundles (`shared_3param`, `shared_pure_cubic`) are
recorded when they differ beyond `T_a0` (the `variant` field disambiguates). The
production candidate is `shared_pure_cubic` and the presets are wired to it
(§10).

---

## 9.7 Outcome — RAN 2026-06-11, verdict PASS (delivery record in `drag_migration_log.md`)

Implemented and run the same day the pre-run clarifications were locked
(N=50, seed 20260604, 20 ps @ 0.01 ps; anchors `a0=14.5556, b0=2.0534` from
the 18 Å Method-A bundle; joint objective bitwise-deterministic, verified
by DRY_RUN).

- **Stage 1 (held-out 9 Å prediction from the §8 18 Å bundle, recorded
  non-gating): PASS.** RMSE **0.2685 Å/ps** ≤ 0.45 band — comfortably
  *below* the ~0.39 Tier-0 dimensionality residual the band was anchored
  on — with full escape at the 18 Å-fit `E_bind`. The 18 Å-calibrated law
  generalizes to the untouched 9 Å case.
- **Stage 2 primary (`shared_3param`): ALL BANDS PASS.** `a` ran to the 0
  bound (0.0002 amu/ps), `b = 2.5159 amu·ps/Å²`, `E_bind = 0.1168 eV`;
  per-case RMSE 18 Å **0.1345** (≤ 0.19), 9 Å **0.1240** (≤ 0.45), escape
  1.0/1.0; all three multi-starts converged to the same point — the §8
  `a`↔`E_bind` ridge is resolved at `a → 0`.
- **`shared_pure_cubic`: ALL BANDS PASS**, `b = 2.5154`, `E_bind = 0.1168`,
  indistinguishable per-case RMSEs.
- **`T_a0` = −0.000000 Å/ps ≤ 0.005 → EQUIVALENT.** The pure-cubic reduced
  form is recorded as the **empirical conclusion on `a`'s identifiability**
  (§9.2): the linear term carries no information over these windows; the
  transport law is effectively `γ = g·b·v²`.
- **Stage-1↔Stage-2 shift (qualitative, recorded):** Δa = −1.456 (to the
  bound), Δb = −0.800, ΔE = +0.046 eV. The shared `E_bind` 0.117 eV sits
  between the per-case 0.071 (18 Å) and 0.154 (9 Å).
- **Verdict (per the §9.4 mapping): PASS → Tier 1 UNGATES; VMI becomes the
  post-Tier-1 final arbiter.** The diagnostic variant was not needed.
- **Artifacts:** both variant bundles under
  `data/reference/drag/shared/trajectory_matching/{shared_3param,shared_pure_cubic}/fit_parameters.json`
  (each loads through `load_drag_coefficients`; sensitivity-only uncertainty
  band, seed sweep omitted per the §8 seed-insensitivity finding), plus
  `stage1_prediction.json` and `verdict.json`; per-case
  `held_out_validation.json` files updated. Production candidate =
  `shared_pure_cubic`; presets are wired to it (§10). Regression coverage:
  `tests/test_extraction_trajectory_matching.py::TestSharedBundleArtifacts`.

---

## 10. Production candidate + alternative-form discrimination phase

Defines the production candidate and the alternative-form discrimination phase
(now complete — outcome §10.7). Decision/delivery records in
`drag_migration_log.md`.

### 10.1 Decisions

1. **Production candidate = `shared_pure_cubic`.** The two §9.7 bundles are
   `T_a0`-equivalent; `shared_pure_cubic` (`a = 0` exactly) is the honest
   encoding of the §9.7 identifiability conclusion (`shared_3param`'s
   `a = 0.0002` is optimizer noise at a bound).
2. **Presets are wired to `shared_pure_cubic`.** Both drag presets load the
   shared bundle and set `binding_energy_I_ion_eV` from its stamped
   `effective_binding_energy_I_ion_eV` (the §6.5.1 identity holds by
   construction); the transitional `allow_unvalidated_binding_pairing` hatch is
   removed and the guard runs at full strength.
3. **Alternative drag-form discrimination phase:** realize the reserved
   `power_law` and `linear_quadratic` forms (the latter including its
   pure-quadratic `a ≡ 0` variant) behind the interchangeable surface, and run
   each family through the same shared-form joint-refit machinery (§9) against
   the pure-cubic incumbent. `threshold` stays reserved (out of scope).
4. **Ordering:** the form phase precedes Tier 1 — the transport form is settled
   before mass dynamics builds on it. VMI remains the post-Tier-1 final arbiter.

### 10.2 Motivation — the exponent question

§9.7's empirical conclusion is a **pure-cubic force law**:
`F_drag = g·b·v³` (`γ = g·b·v²`, exponent `n = 3`). But the Method-A
power-law export on the same force-balance data independently found
**`n ≈ +2`** (`F ∝ v²`, the inertial/form-drag wing — design doc §3
finding note). Two extraction routes therefore point at **different
effective exponents (3 vs 2)**. This phase asks two questions, and both
answers are informative:

- **Discrimination:** can the full-window trajectory objective tell the
  exponents apart at all? The §8/§9 weak-`a` finding showed it cannot see
  the *linear* term; whether it can resolve `v²` vs `v³` force scaling over
  these windows is an open empirical question.
- **Generalization:** does any alternative family carry **both** cases
  (the shared, over-constrained fit) as well as or better than pure-cubic?

If the objective discriminates → a sharper law (and possibly a new
production candidate). If it is degenerate across exponents → that
degeneracy is itself the recorded finding (the trajectory objective fixes
the *magnitude scale* of the drag but not its exponent), pure-cubic stays
the candidate, and the exponent question passes to the post-Tier-1 VMI
arbiter.

### 10.3 Forms and dimensional analysis (working-method rule)

All forms keep the unified friction convention: `γ(v)` is a **force
coefficient** [amu/ps], `F_drag = γ(v)·v`, mass never enters the drag
module. `g = g(depth)` is the dimensionless spatial gate. `γ` is always
exposed via its **closed form**, never via `|F|/v` (the Slice-1 rule).

- **`linear_quadratic`** (coefficients `{a, c}`, design doc §3.8):
  $$\gamma(v) = g\,(a + c\,|v|), \qquad
    F_\text{drag} = g\,(a\,v + c\,|v|\,v).$$
  Units: `[a] = amu/ps`; `[c·|v|] = (amu/Å)·(Å/ps) = amu/ps` ✓. Force
  `[γ·v] = amu·Å/ps²` ✓. **Pure-quadratic variant:** `a ≡ 0`,
  `F = g·c·v²` — the Method-A `n ≈ +2` hypothesis as a closed form.
- **`power_law`** (coefficients `{C, n}`, design doc §3.8):
  $$F_\text{drag} = g\,C\,|v|^{n}, \qquad
    \gamma(v) = g\,C\,|v|^{\,n-1}.$$
  Units: `[C] = amu·Å^(1−n)·ps^(n−2)`, `n` dimensionless. Check:
  `C·|v|^(n−1) → amu·Å^(1−n)·ps^(n−2)·Å^(n−1)·ps^(1−n) = amu/ps` ✓.
- **Nesting identities (cross-check obligations for the implementation):**
  `power_law(n=2, C=c) ≡ pure-quadratic(c)` and
  `power_law(n=3, C=b) ≡ pure-cubic(b)` exactly. The free-`n` fit
  therefore **nests both incumbents**: the fitted `n̂` (with its
  sensitivity half-width) is the direct measurement of exponent
  identifiability — the sharpest single number this phase produces.
- **Dissipativity (§3.3 guard arms, mirroring the §9.5 relaxation):**
  `linear_quadratic`: `a ≥ 0`, `c ≥ 0`, `a + c > 0` (strictly dissipative
  for `v > 0`, no turnover; pure-quadratic = `a = 0, c > 0`).
  `power_law`: `C > 0`; for `n ≥ 1`, `γ` is finite at `v = 0`. A fitted
  `n < 1` would make `γ` diverge at rest — the fit bounds (§10.4) keep
  `n` above 1, so `drag_low_v_floor` stays inert; if a future re-fit ever
  releases that bound, the floor obligation (design doc §3.8) goes live.

### 10.4 Design — same machinery, per family

- **Stage-1 analog per family (recorded, non-gating, honestly weakened):**
  fit the family on 18 Å only, score the untouched 9 Å prediction. Unlike
  §9's Stage 1 this is **no longer strictly held-out** — the 9 Å data has
  been seen repeatedly — so it carries audit-trail value (does the
  18 Å-calibrated family generalize the way pure-cubic did at
  0.2685 Å/ps?), not verdict power.
- **Stage-2 shared joint fit per family**, same objective (per-case
  RMSE + escape penalty, then equal-weight mean), same windows, same
  N/seed/anchoring discipline. **Anchor sourcing (user decision
  2026-06-12, Method-A future-role record in `drag_migration_log.md`):
  the anchors are PRE-REGISTERED CONSTANTS** — values read off the
  Method-A bundles once (18 Å `linear_and_cubic` for the `{a, c}`-family
  scale; the `power/` export may inform the `{C, n}` scale) and locked
  with provenance in the script's USER SETTINGS at the pre-run session,
  **not** a runtime read of any bundle (the preset-derived `setup.a0/b0`
  reads `a0 = 0` from the shared bundle since the re-wiring and must not
  be used). Family-specific anchor mapping (values **LOCKED** in §10.4.1;
  carried as named constants in the driver's USER SETTINGS at
  implementation):
  - `lq_shared_3param` — `{a, c, E_bind}`, `a` lower bound 0;
  - `lq_shared_pure_quadratic` — `a ≡ 0`, fit `{c, E_bind}`; `T_a0`-analog
    equivalence classification between the two, as in §9.2;
  - `pl_shared_3param` — `{C, n, E_bind}`, `n` free within the locked
    bounds `[1, 4]` (§10.4.1); the optimizer works in the **pivot
    parameterization** `(γ_ref, n)`, `γ_ref = C·v_ref^(n−1)` at the
    locked pivot speed `V_REF_APS` (§10.4.1) — it axis-aligns the
    `log C ≈ const − n·log v̄` matching ridge so the fitted `n̂`'s
    sensitivity half-width is meaningful; the stamped bundle records raw
    `{C, n}` (the §10.3 closed form — loader contract unchanged).
- **Bands: the §9.4 Stage-2 bands are reused unchanged**
  (`S2_RMSE_18A_MAX` 0.19, `S2_RMSE_9A_MAX` 0.45, escape 1.0). They are
  form-agnostic statements about trajectory reproduction; reusing them
  avoids any post-hoc tuning. **The new pre-registered numbers are
  LOCKED (2026-06-12 pre-run session, §10.4.1):** the form-equivalence
  classification `T_form` scores
  Δ = (a family's **best-variant** objective) − (the pure-cubic
  incumbent's 0.1293 Å/ps) under a **two-threshold scheme** —
  equivalent within ±`T_FORM_EQUIV_APS` (0.005 Å/ps, the
  optimizer-noise/`T_a0` scale); "genuinely better → escalate" only
  beyond `T_FORM_BETTER_APS` (Δ ≤ −0.013 Å/ps, the 10%-sensitivity
  scale of the incumbent objective); the zone between is recorded as
  "marginally better, not escalation-worthy"; worse beyond +0.005 →
  rejected. *Note:* unlike the §9 nested `a ≡ 0` case, the alternative
  families are **not** nested in the incumbent, so Δobjective can be
  negative (a genuinely better form) — and the asymmetric escalation
  bar keeps an improvement inside the objective's own flatness scale
  from entering the record stamped "genuinely better".
- **Verdict mapping (two-threshold, §10.4.1):** family passes the §9.4
  bands AND beats the incumbent beyond `T_FORM_BETTER_APS` → competing
  production candidate, decision escalated to the user. Within
  ±`T_FORM_EQUIV_APS` → **exponent degeneracy recorded** as the finding;
  pure-cubic stays candidate; VMI post-Tier-1 arbitrates. In the
  marginal zone (−0.013 < Δ < −0.005) → recorded "marginally better,
  not escalation-worthy"; pure-cubic stays. Fails bands or worse beyond
  +`T_FORM_EQUIV_APS` → incumbent confirmed, family bundle recorded as
  rejected-by-trajectory-objective.
- **Methodological status of the axes (recorded honestly):** the §9
  cross-case axis was *spent* by the Stage-2 joint fit, and these
  comparisons re-use the same two trajectories — this is **model selection
  on seen data**, legitimate for ranking forms under a pre-registered
  protocol but not fresh held-out validation. No new validation claim is
  made; the winner's external test remains VMI after Tier 1.

### 10.4.1 Pre-run constants — LOCKED (2026-06-12 pre-run session)

First-runs rule satisfied: fixed before any form-phase run, not re-tuned
after. Decision record in `drag_migration_log.md`. The driver scripts
carry these as named constants with this provenance in their USER
SETTINGS blocks.

| constant | value | provenance / anchoring |
|---|---|---|
| `T_FORM_EQUIV_APS` | 0.005 Å/ps | equivalence-zone half-width; the optimizer-`fatol`/`T_a0` scale (§9.4) |
| `T_FORM_BETTER_APS` | 0.013 Å/ps | escalation bar (Δ ≤ −0.013 = genuinely better); the 10%-sensitivity scale of the incumbent objective 0.1293 Å/ps |
| `power_law` `n` bounds | `[1, 4]` | nesting points 2 and 3 interior; `n ≥ 1` keeps `γ` finite at `v = 0` (`drag_low_v_floor` stays inert, §10.3) |
| `V_REF_APS` | 3.0 Å/ps | `power_law` pivot speed — inside both cases' in-window speed ranges (18 Å 2.54–3.02, 9 Å 2.83–4.95 Å/ps) |
| `a0` | 14.555626399148123 amu/ps | `data/reference/drag/18A/linear_and_cubic/fit_parameters.json` |
| `b0` | 2.0534044239692157 amu·ps/Å² | same bundle |
| `c0` | 11.016050300970692 amu/Å | `data/reference/drag/18A/quadratic/fit_parameters.json` (Method-A quadratic export, commit `d1ca029`) |
| `C0` | 10.36139380949775 amu·Å^(1−n)·ps^(n−2) | `data/reference/drag/18A/power/fit_parameters.json` |
| `n0` | 2.0557526931077588 (dimensionless) | same bundle |

Decisions bound to these numbers (user, 2026-06-12):

- **Two-threshold `T_form`** (supersedes the original single-0.005
  candidate): Δ = best-variant objective − 0.1293 is classified
  better-beyond-`T_FORM_BETTER_APS` (escalate) / marginally-better
  (−0.013 < Δ < −0.005; recorded, no escalation) / equivalent within
  ±0.005 / worse beyond +0.005 (rejected). Rationale: 0.005 is the
  right scale for "indistinguishable from optimizer noise" but
  hair-trigger for "genuinely better" — the asymmetric bar keeps the
  verdict record honest (a Δ = −0.007 must not enter the record stamped
  "genuinely better" when it sits inside the objective's own flatness
  scale).
- **`power_law` pivot parameterization** `(γ_ref, n)` with
  `γ_ref = C·v_ref^(n−1)` at `V_REF_APS` — optimizer-internal only; the
  bundle stamps raw `{C, n}`.
- **Anchors stay 18 Å-only** — now a convention choice (the 9 Å
  Method-A artifact is restored, migration log 2026-06-12), kept for
  consistency with the §9.3 record; anchors are numerical conditioning
  only.
- **Stage-1 analog convention:** each family's 18 Å-only fit uses the
  family's full variant (`lq_shared_3param` / `pl_shared_3param`),
  scored against the same 0.45 Å/ps number for comparability with
  pure-cubic's 0.2685 (recorded, non-gating per §10.4).
- **Interpretive context (recorded for reading the outcome):** the
  in-window reference speed ranges are narrow (18 Å 2.54–3.02 Å/ps,
  ratio 1.19; 9 Å 2.83–4.95, ratio 1.75), so exponent leverage comes
  mainly from the from-onset transient and the **cross-case
  speed-scale difference**, not the in-window shape — a large `n̂`
  half-width (the degeneracy outcome) is a live expectation, and §10.2
  already treats it as a finding, not a failure.

### 10.5 Implementation surface (delivered — see §10.7)

- `physics/drag.py` — realize `linear_quadratic` and `power_law` behind
  the existing dispatch (a branch, not a signature change — the Slice-1
  promise); closed-form `γ` per form; nesting-identity unit tests.
- `config.py` — §3.3 guard arms per §10.3; `load_drag_coefficients`
  content validation for the new coefficient keys (`{a, c}`, `{C, n}`).
- `i2_helium_md/extraction/trajectory_matching.py` — form-generic
  parameter mapping in the joint objective / shared fit (currently
  `linear_cubic`-specific), writer support for the new families.
- `scripts/extraction/` — extend `method_b_shared_refit.py` (or sibling
  scripts per family, USER SETTINGS + DRY_RUN pattern), `T_form` and `n`
  bounds as named pre-registered constants.
- Tests — per-form dissipativity/guard coverage, nesting identities,
  extraction fit recovery on stub laws, artifact checks.
- Artifacts — per-family bundles under
  `data/reference/drag/shared/trajectory_matching/<variant>/` with the
  `variant` provenance field; a comparison verdict record alongside
  `verdict.json`.

### 10.6 Definition of done — complete

All items done (delivery records in `drag_migration_log.md`): production
candidate `shared_pure_cubic` chosen and presets re-wired (transitional §6.5.1
hatch removed); the §10.4.1 pre-run constants locked; `linear_quadratic` and
`power_law` realized behind the dispatch with guard arms, form-generic loader,
and nesting-identity/dissipativity/extraction-recovery/artifact tests (full suite
705 passed / 0 failed); both families run through the §9 machinery; incumbent
`shared_pure_cubic` **confirmed** (§10.7, no escalation). Tier-1 start is now on
the table.

### 10.7 Outcome — RAN 2026-06-14, incumbent CONFIRMED (delivery record in `drag_migration_log.md`)

Both families ran through the §9 machinery (N=50, seed 20260604, 20 ps @ 0.01 ps;
anchors the §10.4.1 LOCKED 18 Å-only constants passed explicitly — NOT the
re-wired preset's `a0 = 0`). Scored against the reused §9.4 Stage-2 bands plus
the two-threshold `T_form` vs the exact incumbent objective
`0.1292649398514104 Å/ps`. **Model selection on seen data** (the §9 cross-case
axis was spent) — a pre-registered ranking, not fresh held-out validation; VMI
after Tier 1 remains the external arbiter.

- **`power_law` (`pl_shared_3param`, free `n`): EQUIVALENT → pure-cubic
  confirmed.** All four starts converged to **`n̂ = 2.927`** (C ≈ 2.835,
  E_bind ≈ 0.113 eV), objective **0.129171** → Δ = **−0.0001 Å/ps**, inside
  ±`T_FORM_EQUIV` (0.005). Bands PASS (18 Å 0.1305, 9 Å 0.1278, escape 1.0/1.0).
  The free-exponent fit **independently recovers the pure-cubic exponent**
  (`n ≈ 3`), *not* the Method-A `n ≈ 2.06` — resolving the §10.2 tension in
  favour of `n = 3`. `n_err` half-width **0.279** (a genuine measurement, not a
  wild degeneracy: the trajectory objective *does* see the exponent here).
  Stage-1 analog (18 Å→9 Å predict, non-gating): 0.411 ≤ 0.45 PASS.
- **`linear_quadratic` (forced `n = 2`): WORSE → rejected-by-objective.** Both
  variants pass the §9.4 bands (18 Å 0.131, 9 Å 0.195/0.195, escape 1.0) and
  collapse to the **pure-quadratic** corner (`a → 0`, `T_a0`-analog
  −1e-6 EQUIVALENT, `c ≈ 12.8`, E_bind ≈ 0.048 eV), but the family objective
  **0.16315** → Δ = **+0.0339 Å/ps**, well past +`T_FORM_EQUIV` → recorded
  **rejected-by-trajectory-objective**. Stage-1 analog 9 Å prediction **0.699 >
  0.45 FAIL** (non-gating) — the `n = 2` family does **not** generalize from
  18 Å the way pure-cubic did (0.2685).
- **Joint conclusion:** the full-window trajectory objective **discriminates the
  exponent** (contra the §10.2 worry that it might be exponent-degenerate the way
  it was `a`-degenerate): forced `v²` is measurably worse, and the free-`n` fit
  lands at `v³`. **No alternative family beats `shared_pure_cubic`; no
  escalation.** The incumbent stands; the Method-A `n ≈ 2` is not supported by
  trajectory matching. The exponent question does **not** pass to VMI undecided —
  it is settled here at `n ≈ 3`, with VMI (post-Tier-1) the final external check.
- **Artifacts** under `data/reference/drag/shared/trajectory_matching/`:
  `{lq_shared_3param,lq_shared_pure_quadratic,pl_shared_3param}/fit_parameters.json`
  (each loads through `load_drag_coefficients`; sensitivity-only band, seed sweep
  omitted per §8), `stage1_analog_{linear_quadratic,power_law}.json`,
  `verdict_{linear_quadratic,power_law}.json`, and the combined
  `form_comparison_verdict.json`. **Presets remain on `shared_pure_cubic`
  (NOT re-wired by this phase).**
