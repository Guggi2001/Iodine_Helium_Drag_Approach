# Tier 2 — Phase D Implementation Plan (Generative-vs-Anchored Bridge)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python, no pseudo-code until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger. Equations are the *locked* formulations from
> `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`; module descriptions are *interface
> contracts*, not implementations.
>
> **Entry docs:** `TIER2_IMPLEMENTATION_PLAN.md` §3/§4 (Phase D = Slice Z in the
> dependency graph); `TIER2_PHASE_C_IMPLEMENTATION_PLAN.md` (the `biphasic` driver +
> 5-term invariant this slice *runs*); `MASS_DYNAMICS_LOCKED_*.md` §6.11 (the two early
> scalars `f_int`/`τ`, the crossing `t×`, the `Π` order parameter, the regime axis);
> `CALIBRATION_MAP.md` (the priored bands). **Predecessor / kinematic targets:**
> `docs/drag_port/Tier1/TIER1A_IMPLEMENTATION_PLAN.md` (the anchored 21→19→14 run) and the
> frozen 9 Å TDDFT reference `data/reference/9A_All_Data.csv`.

---

## 0. Status and intent

Tier 0 / Tier 1a are delivered; Tier-2 Phases A (energetics), B (channels), and C (schema
+ generative driver) are delivered *as plans*. **Phase D is the bridge**: a single
*validation-first* slice (**Z**) that drives the Phase-C `biphasic` integrator at the
**9 Å / 0.80 eV** condition and checks that the *emergent* mass dynamics reproduce the
**known kinematics** — the Tier-1a anchored shell decline and the 9 Å TDDFT trace —
**before** the experimental size-distribution tier (Phases E/F).

Phase D asks **one** question: *with the priored knobs (λ₀, f_int, f_ret, τ) set to a
single representative in-band point, does the generative mechanism's emergent mean `n(t)`
reproduce the anchored 21→19→14 staircase and the 9 Å TDDFT loss, with `t×` near the
GAH25 prior and the 5-term invariant closing?* This is the **full-driver smoke** validation
step (item 7 of the scientific-code-caution order) *gated by the invariant* — the first
time the composed driver meets a real physical target.

**Central engineering fact.** Phase D builds **no new physics and no new modules**. It is a
**run + compare + report** slice: it composes the delivered Phase-C driver, the Tier-1a
run-orchestration scaffolding (`tier1a_common` / `gen_tier1a_runs`), and the existing
`compare_*` post-process helpers, plus **two tiny reconstruction helpers** (`t×`, `Π`). The
correctness bar is **editorial** (a reported diagnostic), *not* an automated pass/fail —
mirroring the Tier-1a reporting-gate stance. Its job is to catch a mechanism/parameter bug
cheaply, with the actionable lever attached, before the costly experimental tier.

**Locked decisions (2026-06-29, this session).**
- **Minimalistic existence demonstration.** A **single representative priored parameter
  point** (not a sweep, not a mean-field layer, not a calibration) — the goal is to show a
  priored point *exists* that lands the 9 Å kinematics. Parameter *sweeps* and the
  generative-vs-anchored *quantitative* arbitration stay in Phase F.
- **Ensemble = the ions in one generative run.** The driver already propagates many
  independent ion trajectories per run; the mean `n(t)` is averaged over those ions. No
  multi-run seeding.
- **Minimal reconstruction self-contained in D.** Slice Z builds only the thin pieces it
  needs (mean `n(t)` from the v7 `n_shell` array; tiny `t×` and `Π(t)` helpers; reuse
  `compare_distance` / `compare_velocity_magnitude` for `R(t)`/`|v(t)|`). **Phase-E Slices
  H / W / D2 are not pulled forward** — terminal-n extraction, the abundance loader, and
  the Wasserstein comparison stay in Phase E.
- **Targets are loaded, not regenerated.** The Tier-1a anchored `n(t)` artifact and the
  frozen 9 Å TDDFT curve are read as fixed comparison targets (no anchored-artifact churn).

**Locked decisions (2026-07-02, pre-build refinement — supersede where they overlap).**
- **Staircase target = the schedule family, not a run artifact.** The anchored `n(t)` is a
  *deterministic schedule*, numerically identical to the delivered run artifacts' `n_shell`
  (which are gitignored / machine-local). The bridge evaluates
  `build_shell_schedule(t★)` for **t★ ∈ {0.5, 5.0, 9.0} ps** as an overlay *family* (band),
  with **t★ = 5.0 primary** (the generative `t×` at the pinned point lands ≈ 5.2 ps, aligning
  with it). This honors "loaded, not regenerated" with zero artifact fragility.
- **Representative priored point pinned (§4):** λ₀ = 0.9/ps, τ = 6.5 ps (geometric mid of
  [2.6, 16.5]; coincides with GAH25's shell-1 `t₀`), **f_int = 0.5**, f_ret = 0.1,
  picture = `statistical_mixture`, **κ = 1.0** (the config default and the Slice-L oracle
  sweep center); ν = 2.42/ps and s = 3n−3 are Sourced/Derived (no choice). `f_int` is chosen
  to land `t×` in the GAH25 window (see §2.2's closed form) — "just above the 0.21–0.24
  mixture floor" (f_int ≈ 0.3) would give `t×` ≈ 1.9 ps, a factor ~3 below the prior; the
  existence demonstration deliberately takes the point that lands the prior.
- **`t×` is per-ion and has a sharp analytic oracle** (§2.2): the helper returns the per-ion
  first-crossing array; in the bridge run all ions must agree (pre-crossing dynamics are
  deterministic), and the value has a closed form — the check is *sharp*, not editorial.
- **Π oracle corrected** (§2.2): at the 9 Å / 0.80 eV priored central point the expected
  regime is **Π < 1 throughout** (freeze side; gentle shell-retaining regime), **Π = 0 at
  gate-open** (full-shell Langmuir cap), **Π → 0 at exit**. The earlier "Π > 1 during dense
  traversal" expectation belongs to the strip end of the regime axis and is *wrong* here.
- **Π helper re-derives `ρ_ratio(t)`** from checkpoint positions + `droplet_radii_angstrom`
  + the cfg density-gate steepness — the reconstruction step is part of the helper contract,
  not an unspecified input.
- **Run parameters inherit Tier-1a exactly** (N = 50 molecules → 100 ions, 30 ps, Tier-0
  dt/seed/drag bundle, `allow_inconsistent_mass_pairing=True`); only `mass_scenario` and the
  knob values differ.
- **Deliverable shape:** helpers in a new `postprocess/bridge_diagnostics.py` (Phase-E D2
  generalizes it); orchestration script mirroring `gen_tier1a_runs.py`; a report script
  emitting overlay figures into the run dir + the numeric summary; the editorial verdict in
  `TIER2_PHASE_D_BRIDGE_FINDINGS.md` + the tier2 log (TIER0_FINDINGS precedent).
- **Carry pickup — already closed upstream.** The Phase-C follow-up #4 (config→module
  threading test for `gate_onset_override_eV` / ν) was found to be **already delivered**
  (`test_gate_onset_override_threads_to_evaporation` + `test_nu_prefactor_threads_to_evaporation`
  in `tests/test_biphasic_step.py`, landed in the Slice-G re-review pass). Phase D adds no
  duplicate; the follow-up is retired in the log.

---

## 1. Scope lock — what Phase D does and does not do

**Does:**
- run **one** generative `biphasic` trajectory set at the **9 Å / 0.80 eV** validation
  condition (same conditions as the Tier-1a anchored run; `mass_scenario` swapped to
  `biphasic`) at a single representative priored parameter point, into a self-describing
  run dir;
- extract the **emergent mean `n(t)`** (average of the v7 `n_shell` array over the run's
  ions), and reconstruct **`t×`** (first `E_int(t) < Σ(n)`) and **`Π(t) = λ(n)·f_ret·τ`**
  with two tiny helpers;
- **compare** against the *fixed* targets — the anchored 21→19→14 staircase family
  (schedule-evaluated, §2.1), the 9 Å TDDFT loss curve, the GAH25 `t×` prior + the §2.2
  closed form, and the corrected freeze-side `Π` expectation (§2.2) — reusing
  `compare_distance` / `compare_velocity_magnitude` for `R(t)`/`|v(t)|`;
- confirm the **5-term invariant closes** (reuse `ion_ledger_closure`, now 5-term from
  Phase C);
- emit a **written diagnostic report** with the actionable lever attached to any miss.

**Does NOT (fails review if it leaks in):**
- any **parameter sweep / calibration** (Phase F) or the **experimental** size-distribution
  comparison / Wasserstein (Phase E W) — the bridge touches *neither* the abundance CSV nor
  the VMI velocity refs;
- the **mean-field ODE** cross-check (§6.11 flow) — explicitly out of the minimalistic scope;
- **Phase-E machinery** (terminal-n extraction H, t×/Π *generalized* diagnostics D2 — D
  builds only the minimal in-line versions);
- any **new physics**, the **noise model** (Tier 3), changes to the **Tier-0 drag law**, the
  **RNG draw order** (locked at Phase C), or **regeneration of the anchored artifact**.

---

## 2. Locked physics — the bridge targets and reconstructions (the comparison oracles)

These are *comparison targets*, not new derivations. The numbers are the priors the report
is read against.

### 2.1 Kinematic targets (loaded)

| Target | Value / source | Read against |
|---|---|---|
| anchored shell decline | **21 → 19 → 14** He over the 9 Å traversal — the deterministic anchored schedule family `build_shell_schedule(t★)`, **t★ ∈ {0.5, 5.0, 9.0} ps** overlaid as a band, **t★ = 5.0 primary** (numerically identical to the delivered run artifacts' `n_shell`, which are gitignored/machine-local) | emergent mean `n(t)` |
| 9 Å TDDFT trace | `R(t)`, `|v(t)|` from `data/reference/9A_All_Data.csv` (the frozen reference) | `compare_distance` / `compare_velocity_magnitude` |
| `t×` prior (GAH25 `t₀`) | **5.0–6.5 ps**, treated as a **±factor-2** prior (Na⁺ number; I⁺ is Rb⁺-like, so a factor-2 offset is unsurprising; a factor-10 miss is the genuine flag) | reconstructed `t×` |
| budget | `E_avail^ion = 0.80 eV` (`d = 9 Å`, the budget the drag + shell refs were generated under) | scenario-stamped config |

### 2.2 Reconstructions (two tiny helpers, built in D)

**Crossing time `t×` (MASS §6.11 — this is literally GAH25's `t₀`):**
$$
t_\times = \min\{t : E_\text{int}(t) < \Sigma(n(t))\},\qquad \Sigma(n)=\textstyle\sum_{i\le n}D_0(i),
$$
reconstructed post-hoc from the v7 `E_int_eV (2N,T)` array + the Phase-A ladder cumulative
`Σ(n)` (Slice L). Zero schema cost. *Sanity band:* flag only if absurd (`t× < 1 ps` or
`> 15 ps`).

**Closed-form expected value (2026-07-02 refinement — the check is sharp, not editorial).**
Before the gate opens there is *no stochasticity at all*: evaporation is suppressed
(`E_int > Σ`), and pickup is dead because `n(0) = 21 = n*` zeroes the Langmuir cap. So
`E_int` evolves under pure K2 Newton cooling, identically for every ion, and
$$
t_\times \;=\; \tau \,\ln\!\frac{f_\text{int}\,E_\text{avail}^\text{ion}}{\Sigma(21)}
$$
exactly (up to `dt` discreteness and the cooling-before-draw step ordering). At the pinned
point (f_int = 0.5, τ = 6.5 ps, mixture Σ(21) ≈ 0.17–0.19 eV): **t× ≈ 4.8–5.6 ps ≈ 5.2 ps**
— inside the GAH25 window by construction. Consequences: (a) the helper returns the
**per-ion** first-crossing array, and in the bridge run all ions must agree — disagreement
is a wiring flag, not physics; (b) the bridge smoke checks the reconstructed `t×` against
this analytic value, a *sharp* oracle upgraded from the editorial band; (c) `f_int` and the
GAH25 prior are **coupled choices** — the pinned point deliberately lands the prior.

**Regime order parameter `Π` (MASS §6.11):**
$$
\Pi(n) \equiv \lambda(n)\,f_\text{ret}\,\tau\quad(\text{dimensionless}),\qquad
\Pi>1\ \text{shedding persists},\ \Pi<1\ \text{freeze},
$$
with `λ(n)` the pickup rate (Phase B). **Corrected qualitative expectation (2026-07-02;
supersedes the earlier "Π > 1 during dense traversal" wording, which belongs to the *strip*
end of the regime axis and is wrong at this condition):** at the 9 Å / 0.80 eV priored
central point the expected regime is the **freeze side** —
- `Π = 0` at gate-open (`n = 21 = n*` → the Langmuir factor is 0);
- `Π < 1` throughout (even at `n = 14` the occupancy factor is 1/3, giving
  `Π ≈ 0.9 · ⅓ · 0.1 · 6.5 ≈ 0.2` at the pinned point — the 21→14 decline is driven by the
  *initial `E_int` cascade* after gate-open, not by pickup-sustained quasi-steady shedding,
  consistent with MASS §6.11's "self-binds early → retains a shell" gentle-9 Å regime);
- `Π → 0` at exit (`ρ_He → 0`, termination on every trajectory).

A *computed* `Π > 1` at this point would itself be a flag (it requires the upper corner of
every band simultaneously: `f_ret ≳ 0.4`, `τ → 16.5`, `λ₀ → 1.1`). The `Π(t)` trace is
reported either way — it remains the regime-axis diagnostic Phase E generalizes.

*Reconstruction contract:* `λ(n)` needs `ρ_ratio(t)` at each ion's position; the helper
**re-derives it** from the checkpoint positions + `droplet_radii_angstrom` + the cfg
density-gate steepness (same profile the run used), then evaluates `lambda_attach`
per ion. Reported as mean ± envelope over the ions.

### 2.3 Invariant (reused, Phase C)

The **5-term** closure `E_kin + E_pot + E_dissip + E_mass_transfer + E_int ≈ const` must
hold to Verlet drift over the bridge run (reuse `ion_ledger_closure`). A divergent residual
is a **miswire**, surfaced *before* any kinematic interpretation.

---

## 3. New-work slice

### Slice Z — Bridge self-consistency check *(run + compare + report; composes only delivered modules)*

**Modules.** A `scripts/` orchestration entry (mirroring `gen_tier1a_runs`) + the two tiny
reconstruction helpers in a **new `postprocess/bridge_diagnostics.py`** (settled 2026-07-02;
Phase-E Slice D2 generalizes this module rather than re-implementing) + a report script under
`scripts/post_processing/` emitting the overlay figures into the run dir and the numeric
summary. **No new physics module.**

**Purpose.** Demonstrate that the generative mechanism *can* reproduce the 9 Å kinematics at
a single priored point — the de-risk gate before the experimental tier.

**Interface (settled 2026-07-02).**
- *Orchestration:* a script that builds the **9 Å / 0.80 eV `biphasic`** config (reuse the
  Tier-1a conditions + `tier1a_common` scaffolding; swap `mass_scenario="biphasic"`; set the
  pinned priored knobs), runs the Phase-C driver into a self-describing run dir. Run
  parameters inherit Tier-1a exactly: N = 50 molecules (100 ions — both CE fragments are
  ions, so "mean over ions" = mean over all 2N rows), 30 ps, Tier-0 dt/seed/drag bundle,
  `allow_inconsistent_mass_pairing=True`.
- *Reconstruction helpers (`postprocess/bridge_diagnostics.py`):*
  - `crossing_time_ps(E_int_eV, n_shell, time_ps, *, picture, kappa) -> per-ion t× array`
    (first index with `E_int < Σ(n)`, consuming Slice L's `ladder_cumsum`; NaN where the gate
    never opens in-window; all-ions-agree is a bridge wiring check, §2.2);
  - `regime_parameter(ckpt, cfg) -> Π(t) per ion` — re-derives `ρ_ratio(t)` from positions +
    `droplet_radii_angstrom` + the cfg density-gate steepness, then `lambda_attach(...)·f_ret·τ`
    (§2.2 reconstruction contract); reported mean ± envelope.
- *Comparison + report:* mean `n(t) = mean(n_shell, over the 2N ion rows)` overlaid on the
  `build_shell_schedule(t★)` family (t★ ∈ {0.5, 5.0, 9.0}, t★ = 5.0 primary);
  `compare_distance` / `compare_velocity_magnitude` for `R(t)`/`|v(t)|` vs the 9 Å CSV;
  `ion_ledger_closure` for the invariant; assemble the written diagnostic (overlay numbers +
  the `t×`/`Π`/residual summary + the actionable-lever note) into
  `TIER2_PHASE_D_BRIDGE_FINDINGS.md` + the tier2 log entry.

**Encoded form.** §2.2 (`t×`, `Π`); the run physics is the Phase-C driver.

**Knobs (config §4):** consumes the full A/B/C surface at the **pinned representative
priored point** (2026-07-02) — `coulomb_available_eV = 0.80` (scenario-stamped),
**λ₀ = 0.9/ps** (central of 0.7–1.1), **τ = 6.5 ps** (geometric mid of [2.6, 16.5]),
**f_int = 0.5** (well above the 0.21–0.24 mixture floor; chosen so the closed-form `t×`
lands the GAH25 window, §2.2), **f_ret = 0.1** (prior small, nonzero to exercise S1),
picture = **`statistical_mixture`**, **κ = 1.0** (config default / Slice-L oracle center);
ν = 2.42/ps and s = 3n−3 are Sourced/Derived. **No knob is fitted** — the point is chosen,
not calibrated, to demonstrate existence.

**Oracle / acceptance values (reported, not auto-verdict — except where marked sharp).**
- emergent mean `n(t)` overlaps the anchored 21→19→14 staircase family within a **loose
  early-window tolerance** (a qualitative staircase overlay, not a tight metric);
- `t×`: **sharp** — per-ion values all agree, and match the §2.2 closed form
  `τ·ln(f_int·0.80/Σ(21)) ≈ 5.2 ps` (up to dt discreteness); also inside the [1, 15] ps
  sanity band and the GAH25 5–6.5 ps prior (±factor-2) by construction of the pinned point;
- `Π(t)`: **corrected expectation (§2.2)** — `Π = 0` at gate-open, `Π < 1` throughout
  (freeze side, ≈ 0.2 max at the pinned point), `Π → 0` at exit; a computed `Π > 1` is a flag;
- 5-term invariant closes to Verlet drift.
- A **miss is framed with its lever:** mean `n(t)` decays too fast → λ₀ too low or
  evaporation too aggressive; `t×` too late → f_int too high or τ too long. The bridge
  *localizes* a bug; it does not adjudicate fidelity.

**Independence.** Composes only **delivered** modules (Phase-C driver, Tier-1a scaffolding,
`compare_*`). The two helpers depend only on the v7 arrays + Slice L/B. No dependence on
Phases E/F.

**Test spec (`tests/test_bridge_diagnostics.py` + a tiny driver smoke).**
- `crossing_time_ps` on a synthetic `E_int(t)` ramp returns the exact first-crossing index
  per ion; flag-NaNs if the gate never opens within the window.
- `crossing_time_ps` against the **§2.2 analytic oracle** on a tiny cooling-only synthetic
  (pure Newton decay at fixed `n = 21`): the reconstructed `t×` matches
  `τ·ln(E_int(0)/Σ(21))` to dt discreteness.
- `regime_parameter` reproduces `λ·f_ret·τ` on known inputs; `→ 0` as `ρ_ratio → 0`; the
  `ρ_ratio` re-derivation matches a direct `helium_density` call on hand-built positions.
- mean-`n(t)` reduction matches a hand-average on a tiny synthetic `n_shell` array.
- config→module threading (Phase-C follow-up #4): **already covered upstream** by
  `test_gate_onset_override_threads_to_evaporation` / `test_nu_prefactor_threads_to_evaporation`
  in `tests/test_biphasic_step.py` — do not duplicate; verify they pass in the pre-build
  suite run.
- a **few-step** generative smoke run closes the 5-term invariant (reuse the Phase-C gate);
  **no production-sized run, no figures** in pytest.

**Acceptance.** The report is produced from a real (out-of-pytest) bridge run; the helpers
pass their unit tests; the smoke invariant closes. The kinematic "pass" is editorial.

---

## 4. Config contract (Phase D)

No new `SimConfig` fields. Phase D **selects values** on the existing A/B/C surface:
- `mass_scenario = "biphasic"`; `coulomb_available_eV = 0.80` (scenario-stamped → trips the
  §6.5 pairing guard structurally → runs under `allow_inconsistent_mass_pairing=True`, R6);
- the pinned representative priored point (2026-07-02): λ₀ = 0.9/ps, τ = 6.5 ps,
  f_int = 0.5, f_ret = 0.1, κ = 1.0, picture = `statistical_mixture` — chosen, not fitted
  (rationale in §0/§3);
- noise stays inert (Tier 3); the `s≥1` and pairing guards are the Phase-B/C ones, unchanged.

The two reconstruction helpers read config (picture, κ, λ₀, f_ret, τ) but **add no field**.

---

## 5. Dependency / build order (within Phase D)

```
(delivered) Phase C biphasic driver + 5-term invariant
        │
        ├─ helper: crossing_time_ps(E_int, n_shell, t; L.ladder_cumsum) ── unit-tested first
        ├─ helper: regime_parameter(ckpt, cfg; rho-ratio re-derivation + B.lambda_attach)
        │
        ▼
Slice Z: 9A/0.80 eV biphasic run  ->  mean n(t) + per-ion t_x + Pi(t)
         compare vs anchored schedule family (t* in {0.5,5,9}, 5.0 primary)
         + 9A TDDFT (compare_*) + closed-form t_x (sharp) + GAH25 prior
         -> written diagnostic report (reported, not gated)
```

The two helpers are built and unit-tested first (pure, tiny). The run + comparison + report
is the slice body. Nothing in Phase D depends on Phase E/F; Phase D depends only on the
delivered Phase-C driver and the Tier-1a scaffolding.

---

## 6. Testing methodology

- **Helpers:** closed-form/identity asserts on tiny synthetic arrays (first-crossing index;
  `Π` formula; mean-`n` average) — tight tolerances.
- **Driver smoke:** a **few-step** generative run closes the 5-term invariant (reuse the
  Phase-C `ion_ledger_closure`); this is the only "live driver" check in pytest.
- **The bridge run itself** (the production-sized 9 Å / 0.80 eV trajectory + the kinematic
  overlay + the report) is **out-of-pytest** orchestration — like the Tier-1a runs. **No
  figures, no production-sized checkpoints, no reference-data reads in pytest.**

Suites: `tests/test_bridge_diagnostics.py` (+ the smoke in the existing driver test). Run the
narrowest first, then the full suite. Interpreter (Python may not be on PATH):

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q
```

---

## 7. Acceptance criteria

**Slice Z (Phase D overall):**
- the two reconstruction helpers pass their unit tests; the few-step driver smoke closes the
  5-term invariant;
- a real bridge run at 9 Å / 0.80 eV produces the **written diagnostic report**
  (`TIER2_PHASE_D_BRIDGE_FINDINGS.md` + tier2 log entry): mean `n(t)` overlay vs the
  anchored staircase family (t★ ∈ {0.5, 5.0, 9.0}, t★ = 5.0 primary) + 9 Å TDDFT, the
  reconstructed per-ion `t×` vs the §2.2 closed form and the GAH25 prior, the `Π(t)` regime
  trace (corrected freeze-side expectation), and the invariant residual;
- the report carries the **actionable-lever** framing for any miss.
- The kinematic "pass" is **editorial, not a pytest gate** (mirrors the Tier-1a reporting
  stance). A miss is a *flag* that localizes a mechanism/parameter bug before Phase E/F.

**Out-of-scope guard (fails review if it leaks in):** any Phase-D code path that runs a
parameter sweep / calibration (Phase F), reads the experimental abundance/VMI references or
computes Wasserstein (Phase E W), builds the Phase-E terminal-n extraction (H) or the
generalized t×/Π diagnostics (D2), integrates the mean-field ODE, reads a noise amplitude,
changes the Tier-0 drag law or the locked RNG draw order, or regenerates the anchored
artifact.

---

## 8. Risks / notes carried

- **The bridge is an existence demonstration, not an arbiter.** A single priored point
  landing the 9 Å kinematics shows the mechanism *can* reproduce them; it does **not** pin the
  parameters (that is Phase F) and it does **not** validate the experimental observable (that
  is Phase E). Do not over-read a single-point match as calibration.
- **TDDFT is not ground truth — experiment arbitrates at Tier 2.** The 9 Å staircase and the
  GAH25 `t×` are *priors with wide bars* (the predictive shell-timing "1b" variant was
  rejected for exactly this reason). The bridge lands *near* them for the right reasons; a
  factor-2 `t×` offset is expected (Na⁺→Rb⁺-like I⁺), a factor-10 miss is the flag.
- **`E_int` is a constructed reservoir (R2).** If the bridge misses, the §6.11 levers (`t×`
  via f_int/τ; mean-`n` decay via λ₀/evaporation; Π via the regime axis) localize *which*
  construction knob is off — the report must carry this mapping so the miss is diagnosable.
- **Minimal reconstruction now, generalized in Phase E.** The `t×`/`Π` helpers built here are
  the thin versions; Slice D2 (Phase E) generalizes them into the full regime-determination
  diagnostic. Keep them small and reusable so D2 extends rather than re-implements.
- **Mechanism locked, values open** — the bridge *selects* priored values to demonstrate
  existence; it fits nothing. Calibration is Phase F.
- **`f_int = 0.5` is a demonstration choice, not a finding.** It is picked so the closed-form
  `t×` lands the GAH25 prior — and it sits *inside* the derived 0.80 eV window: above the
  mixture floor 0.21–0.24, below the MASS A7 **scenario-invariant velocity-consistency
  ceiling ~0.6** (`E_trans/E_avail ≈ 40%` at both budgets). The older "soft upper ~0.2"
  advisory (CALIBRATION_MAP row 14 wording) is vacuous at 0.80 eV — the floor sits above
  it (the known A7 margin compression); the ~0.6 ceiling is the operative upper edge,
  provisional pending OQ2. Phase F owns the real value.
- **The freeze-side Π expectation is condition-specific.** At 2.70 eV production (Phase E/F)
  the regime question re-opens; do not carry "Π < 1" forward as a general oracle — it is the
  9 Å / 0.80 eV priored-central-point expectation only.
