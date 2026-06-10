# Method B — Trajectory-Matching Drag Extraction

**Status:** Method specification. Defines an **alternative drag-extraction
method** that fits the drag coefficients by minimizing the *forward-integrated
in-window trajectory* RMSE against the smoothed reference — as opposed to the
current Method A, which regresses the force-balance residual `F_drag` against
speed directly. No implementation until `[PROCEED TO IMPLEMENTATION]`.

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

## 8. Definition of done

- B implemented as a named extraction method producing `{a,b}` by minimizing the
  forward-integrated in-window trajectory RMSE against the same-smoothed
  reference, with seed-sweep uncertainty and a provenance stamp
  (`extraction_method = trajectory_matching`).
- At least one **held-out axis** exercised and recorded (held-out window the
  minimum; held-out case strongly preferred as the transport-physics signal).
- 18 Å B-fit produced and held-out-validated (the trustworthy case).
- 9 Å B-fit produced **with the radial-to-non-radial flag** and **held-out
  validated**; if it fails held-out, the 9 Å law is recorded as not-yet-usable for
  the radial MD (the honest finding), not silently committed.
- The Tier-0 docs rewritten from consistency to held-out generalization (separate
  doc pass).
- A one-line verdict per case: does the B-extracted `linear_cubic` generalize off
  its fit data, and does the 9 Å fit survive the held-out check or stand flagged?
