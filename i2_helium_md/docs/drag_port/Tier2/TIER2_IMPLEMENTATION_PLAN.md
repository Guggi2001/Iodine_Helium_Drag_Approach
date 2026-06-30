# Tier 2 — Implementation Plan (Generative Mass Dynamics + Size-Distribution Comparison)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python, no LAMMPS, no pseudo-code until
> the explicit `[PROCEED TO IMPLEMENTATION]` trigger. Equations are the *locked*
> formulations from `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`; module
> descriptions are *interface contracts*, not implementations.
>
> **Entry docs:** `DRAG_PORT_DESIGN_DECISIONS.md` §0/§6.4 (Tier 2);
> `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` (mechanism §4, budget §6,
> schema §7, config §11); `CALIBRATION_MAP.md` (parameter classes + Tier anchors).
> Predecessor: `docs/drag_port/Tier1/TIER1A_IMPLEMENTATION_PLAN.md`.

---

## 0. Status and intent

Tier 0 (drag form locked: `shared_pure_cubic`, `γ=g·b·v²`) and Tier 1a (anchored
kinematic mass-dynamics validation) are complete. Tier 2 replaces the *anchored*
`n(t)` schedule with the **generative** `biphasic_energy_gated` mass mechanism and
arbitrates it against the experimental I⁺Heₙ **size distribution** — the only
observable that separates the mass scenarios and the two genuinely-free knobs
(ladder steepness κ + electronic picture).

The mechanism is *locked* (mechanism, energy structure, 5-term invariant) in the
MASS doc; only calibration values are open. This plan breaks that mechanism + its
Tier-2 comparison layer into small, independently-testable, deliverable slices,
following the Tier-1a methodology: pure functions first, mockable RNG, the energy
invariant as the correctness gate, validate in tier order.

**Confirmed scope decisions (2026-06-29):**
1. **Full program** — generative mechanism + comparison layer + calibration campaign.
2. **Validation-first** — first generative runs at the **0.80 eV** (d=9 Å) budget,
   cross-checked against the TDDFT 21→19→14 shell decline and the Tier-1a anchored
   results; switch to the **2.70 eV** production onset only after that lands.
3. **Bridge slice included** — a generative-vs-anchored self-consistency check to
   de-risk the mechanism against known kinematics before facing the experimental
   histogram.

---

## 1. Scope lock — what Tier 2 does and does not do

**Does:** build the generative two-channel mechanism (Poisson pickup + energy-gated
RRK evaporation), the `E_int` reservoir (schema bump + 5-term invariant), the
`D_0(n)` Form-U ladder + electronic picture, K2 Newton cooling, the He density
gating, the generative driver, the size-distribution + per-fragment velocity
comparison layer, and the κ/picture co-fit calibration campaign at both budgets.

**Does NOT (fails review if it leaks in):**
- the **noise model** (Tier 3 — stays inert behind its enum);
- changes to neutral propagation, the **Tier-0-locked drag law**, the checkpoint
  RNG **draw order** (once locked), or default scope beyond the schema bump;
- mechanism re-litigation — only calibration *values* move, via config knobs.

---

## 2. What already exists (reuse) vs. what is new

**Reference data — already present (no export slice needed):**
- `data/reference/integrated_i_he_abundance.csv` — experimental I⁺Heₙ size
  distribution, n=0..20 (`ionCounts`, `ionPercent`). **The Tier-2 arbiter.**
- `data/reference/vmi_summary/vmi_iplus_he.csv` (+ `_high_snr`, `vmi_iplus_gas`) —
  per-fragment velocity references.
- `data/reference/9A_All_Data.csv` — TDDFT shell/trajectory cross-check.

**Code spine — delivered by Tier 1a (reuse, do not rebuild):**
- `physics/mass_jump.py` — cold-shed / continuous reset primitives, reduced-mass
  `dE_mass_transfer`, `kick_factor`. Tier 2 adds a **+He capture** counterpart.
- `physics/baoab.py` — SQ1 O-step (`e^{−γ·dt/m}`, frozen γ, exact dissipation),
  rebuilt every step so `m` can vary. **Reused unchanged.**
- `simulation/ion_propagation_step.py::shed_step` — the jump-then-O seam (SQ2/SQ3).
- `simulation/checkpoint.py` — v6 schema + v5 shim; 4-term ledger arrays.
- `postprocess/energy_balance.py::ion_ledger_closure` — 4-term closure gate.
- `postprocess/velocity_distribution.py` — mass-selected final-velocity histograms.
- `scripts/tier1a_common.py`, `gen_tier1a_runs.py` — run-orchestration scaffolding.
- `config.py` — `mass_scenario ∈ {fixed, biphasic, anchored_discrete}` (`biphasic`
  currently inert / guard-only); `validation_histogram_metric="wasserstein"` declared.

**New for Tier 2 (this plan's slices):** the generative channels (Poisson pickup,
RRK evaporation, self-bound gate), the **E_int reservoir** (new per-ion state →
schema v6→v7, **5-term** invariant), the **D_0(n) Form-U ladder** + electronic
picture, **K2 Newton cooling** of E_solv.struct, the He **density profile** gating,
the **generative driver**, and the **size-distribution comparison + calibration**
layer.

---

## 3. Dependency graph and build order

```
Phase A — energetics primitives (pure, no state, no integrator)
   L (ladder D0(n), ΣD0 gate)   K (cooling K2, E∞(N))   U (E_int budget rules)
        │                          │                       │
Phase B — stochastic channels (stateful primitives; mock RNG; reuse mass_jump reset)
   ρ (He density profile) ──► P (pickup: Poisson + Langmuir + capture reset + S1)
                               Q (evaporation: gate + RRK + cold-shed reset + K1)
        │
Phase C — schema + generative integrator (compose A+B)
   X (checkpoint v6→v7 + 5-term invariant gate)
   G (generative driver: per-step pickup+evap+cooling+E_int; one event/step)
        │
Phase D — bridge (validation-first de-risk, 0.80 eV)
   Z (generative mean n(t) reproduces anchored 21→19→14; t× ~ GAH25; invariant closes)
        │
Phase E — Tier-2 observable (comparison layer)
   H (terminal-n + per-fragment velocity extraction)
   W (abundance/VMI loader + Wasserstein comparison)
   D2 (derived diagnostics: t×, Π(t), regime determination)
        │
Phase F — calibration campaign + production switch
   R2 (run matrix + κ/picture co-fit + bounded sweep; report; then 2.70 eV)
```

**Independent / parallelizable:** L, K, U (Phase A); ρ (Phase B). P and Q each
compose A-primitives but mock the integrator. Everything downstream of G composes
only **accepted** modules. The **5-term invariant closure** is the cross-cutting
correctness gate from X onward.

---

## 4. New-work slices

### Phase A — Energetics primitives (pure functions)

Maps to validation steps 1–2 (direct formula + shape/unit). No state, no RNG, no
integrator. Tested against closed-form oracle values from the MASS doc.

#### Slice L — Dissociation ladder `D_0(n)` + integrated gate threshold
- **Purpose.** The Form-U sigmoid ladder and its cumulative sum (the self-bound
  gate threshold). No physics state.
- **Encoded form (MASS §8 R3 / §11):**
  `D_0(n)=D_floor+(D_0(1)−D_floor)(1−σ(n))/(1−σ(1))`,
  `σ(n)=[1+e^{−κ(n−n*−½)}]^{−1}`; cumulative `Σ(n)=Σ_{i≤n} D_0(i)`.
- **Knobs:** `ladder_electronic_picture ∈ {statistical_mixture(default), x2_only,
  cooling_relaxed}` sets `D_0(1)` (74.4 / 106.9 / between, cm⁻¹); `ladder_steepness`
  κ (Free); `D_floor=4.97 cm⁻¹` (sourced); `n*=21` (sourced).
- **Oracle:** first rung 106.9/74.4 cm⁻¹; floor 4.97 cm⁻¹; `Σ_{i≤21} D_0` = X₂
  0.25–0.28 eV / mix 0.17–0.19 eV (≈κ-independent, ~11% over κ∈[0.3,5]); cliff at
  n*+½; tabulated-ladder fallback.
- **Acceptance.** Rungs/sum match oracle to 4 figures across κ∈[0.3,5] and all three
  pictures; gate threshold ≈κ-independent (decoupling check).

#### Slice K — Newton cooling K2 + occupancy-resolved asymptote
- **Purpose.** The cooling driver for `E_solv.struct` and its binding split.
- **Encoded form (MASS §6 K2):** `dE_solv.struct/dt = −(E_solv.struct − E∞(N))/τ`,
  `E∞(N)=−|S(N)|`, `|S(N)|=|S|·Σ_{i≤N}D_0/Σ_{i≤n*}D_0` (occupancy-resolved),
  `E_elec(N)=−(|S(N)|−Σ_{i≤N}D_0)≤0` (electrostriction, dominant term).
- **Knobs:** `internal_energy_cooling_tau_ps` τ (Bounded sweep [2.6,16.5] ps);
  `solv_struct_asymptote` |S|=0.308 eV; `electrostriction_binding` (derived).
- **Oracle:** E∞→0 as N→0 (OQ6 resolved, strip reachable); eV/ps; cold-shed
  neutrality ΔE_solv.struct=0 up to the marginal-electrostriction bath booking (A8);
  E_int^eq=0 exactly under the split (R12 resolved).
- **Acceptance.** Newton step + split match oracle; neutrality identity holds.

#### Slice U — `E_int` budget bookkeeping rules
- **Purpose.** Pure per-channel ΔE_int helpers + reconstruction; the bricks the
  driver and ledger consume.
- **Encoded rules (MASS §6):** S2 onset `E_int(0)=f_int·E_avail^ion`; S1 pickup
  heating `+f_ret·D_0(n+1)` (remainder to bath); K1 shed drain `−D_0(n)`;
  reconstruction `E_int=E_solv.struct−E_bind^pair(N)−E_elec(N)` (valid only
  **post-t×**, A9).
- **Knobs:** `internal_energy_partition_fraction` f_int (Bounded; scenario-keyed
  floor); `internal_energy_retained_fraction` f_ret (Bounded); `coulomb_available_eV`
  E_avail (0.80 validation / 2.70 production).
- **Oracle:** f_int floor = ΣD_0/E_avail (0.21–0.35 @ 0.80 eV; 0.065–0.10 @ 2.70 eV);
  S1/K1 sign + magnitude; jump-consistency (variable swap leaves S1/K1 unchanged).
- **Acceptance.** Each rule reproduces the MASS §6 identities; reconstruction guarded
  against pre-t× use.

### Phase B — Stochastic channels (stateful; mock RNG)

Maps to validation step 6 (statistical) at the unit level: distribution/moment checks
with a seeded/mock RNG, neighbours mocked. Reuse `mass_jump` reset primitives.

#### Slice ρ — Helium density profile `ρ_He(depth)/ρ_bulk`
- **Purpose.** The density gating field for pickup. Sourced from the baseline/TDDFT
  density (CALIBRATION_MAP row 8). Small, pure; kept separate as a distinct sourced
  quantity feeding P.
- **Acceptance.** Profile →0 outside the droplet (exit-driven termination), unit ratio
  in bulk; shape sane vs the confining-potential radius (14.2 Å).

#### Slice P — Pickup channel (Poisson → Bernoulli)
- **Purpose.** The gain channel: one independent Bernoulli draw per ion per step.
- **Encoded form (MASS §4):** `P_attach=1−e^{−λ dt}`,
  `λ_attach=λ_0(ρ_He/ρ_bulk)(1−n/n*)_+^p` (Langmuir cap, A12). On fire: `n→n+1`,
  `m→m+m_He`, **capture momentum reset** `v⁺=m/(m+m_He)·v⁻` (He at rest), S1 heating
  to E_int via U. Reuses a **+He** counterpart of the `mass_jump` reset.
- **Knobs:** `pickup_rate_coefficient` λ_0 (Sourced+Bounded, ~0.7–1.1/ps prior);
  `pickup_rate_form ∈ {density_only(default), sweeping, dwell_time}`;
  `pickup_occupancy_cap ∈ {langmuir(default), none}`; `pickup_occupancy_exponent` p
  (default tied to κ).
- **Test/oracle:** Bernoulli stats reproduce Poisson P_n(t) under a seeded RNG (A1);
  capture defect = reduced-mass form to machine precision; cap →0 at n=n*; density-only
  ≡ langmuir for ejection (cap inert).
- **Acceptance.** Momentum/defect exact; Poisson moments within sample-size band.

#### Slice Q — Evaporation channel (energy-gated + RRK)
- **Purpose.** The loss channel: self-bound gate then saturating RRK shed.
- **Encoded form (MASS §4):** suppress while `E_int>Σ_{i≤n}D_0` (gate); else
  `k=ν(1−D_0(n)/E_int)^{s−1}`, `P_shed=1−e^{−k dt}`, `s=3n−3` (n≥2), **n=1 direct
  dissociation k=ν**; on fire `n→n−1`, `E_int−=D_0(n)`, **cold-shed reset** (reuse
  `mass_jump.cold_shed`), K1 drain via U.
- **Knobs:** `evap_rate_prefactor_per_ps` ν (Sourced, pinned 2.42 ps⁻¹);
  `evap_rrk_dof` s (Derived 3n−3; override **guarded s≥1** at config-load).
- **Test/oracle:** k bounded in [0,ν) for s>1; k=ν at s=1/n=1; gate suppresses all
  sheds while net self-unbound; no gate-open avalanche (per-step ≤ ~ν·dt); config-load
  guard rejects s<1.
- **Acceptance.** Rate matches oracle; gate + avalanche guarantees hold; s≥1 guard fires.

### Phase C — Schema + generative integrator

#### Slice X — Checkpoint v6→v7 + 5-term invariant gate
- **Purpose.** Persist the new `E_int` reservoir and extend closure to 5 terms.
- **Delta (MASS §7):** add per-step `E_int_eV (2N,T)`; scenario metadata
  `mass_scenario=biphasic`; keep the v6 `n_shell`/`E_mass_transfer_eV` fields;
  **back-compat shim** for v6 (E_int absent → synthesize/flag). Extend `energy_balance`
  to the 5-term closure `E_kin+E_pot+E_dissip+E_mass_transfer+E_int≈const`.
- **Lock the RNG draw order** in the spec (forbidden-list item): pickup vs shed draw
  sequence fixed and documented so results are reproducible.
- **Acceptance.** v7 round-trips; v6/v5 still load; 5-term closure gate reuses
  `ion_energy_totals`; relabel-fault still caught.

#### Slice G — Generative driver (compose A + B into the per-step loop)
- **Purpose.** The production `biphasic_energy_gated` integrator: per step run
  conservative kicks → cooling (K, K2) → **at most one mass event/step** (shed via Q,
  else pickup via P; shed-then-pickup if both draw) **jump-then-O** (reuse the Tier-1a
  seam + SQ1 O-step at post-jump m⁺/SQ3) → update E_int (S1/S2/K1/K2 via U).
- **Reuse:** `baoab.py` O-step, `shed_step` seam, `mass_jump` resets; **only** the
  channel composition + E_int evolution are new.
- **Config guards:** `mass_scenario=biphasic` enters the §6.5 `time_resolved` arm →
  trips the mass↔coefficient pairing guard structurally → runs under
  `allow_inconsistent_mass_pairing=True` (R6 mid-window defence); `s≥1` guard; no-noise
  (Tier 3 stays inert).
- **Acceptance.** Constant-mass / fixed regression unaffected; force-free run conserves
  the 5-term invariant; one-event-per-step + draw-order honored; E_int trajectory finite
  and gate behavior correct.

### Phase D — Generative-vs-anchored bridge (validation-first, 0.80 eV)

#### Slice Z — Bridge self-consistency check
- **Purpose.** Before any experimental comparison, show the generative mechanism *can*
  reproduce known kinematics. With λ_0/f_int/f_ret/τ in their priored bands at **0.80
  eV**, the **mean** `n(t)` should track the anchored 21→19→14 decline and `t×` should
  land near GAH25 5–6.5 ps (±factor-2). 5-term invariant closes.
- **Reuse:** Tier-1a anchored `n(t)` and the 9 Å TDDFT curve as comparison targets;
  `compare_*` helpers for R(t)/|v(t)|.
- **Acceptance (reported, not auto-verdict):** emergent mean n(t) overlaps the anchored
  staircase within the loose early-window tolerance; t× in [1,15] ps and near the GAH25
  prior; Π(t) regime axis behaves (Π>1 during dense traversal). A miss here flags a
  mechanism/parameter bug **before** the experimental tier.

### Phase E — Tier-2 observable (comparison layer)

> **Detailed plan: `PHASE_E_IMPLEMENTATION_PLAN.md`** (slices E1–E5).

Re-sliced H/W/D2 → **E1–E5** with three scope decisions locked (user, 2026-06-29):
(1) build a post-ejection **relaxation stage** (E2) for the R5 truncation — the
locked Q/K/U channels with pickup + drag off, no new physics; (2) **size
distribution only** — velocity comparison cut (the `vmi_summary` refs are aggregate,
not per-fragment); (3) **pure functions only** — scripts deferred to Phase F (R2).

- **E1** — terminal integer-n I⁺Heₙ size-distribution extractor
  (`postprocess/size_distribution.py`; reuses `complex_mass_amu`, `n_shell`).
- **E2** — post-ejection relaxation stage (`simulation/relaxation_stage.py`):
  free-flight evaporation-only cascade to the experimental timescale; reuses the
  Slice-G loop + `mass_jump.cold_shed`; 5-term invariant is its correctness gate.
- **E3** — abundance reference loader mirroring the `hedft_loader` contract
  (`postprocess/abundance_loader.py`).
- **E4** — pure-numpy integer-support Wasserstein `W₁=Σ|F_sim−F_ref|`
  (`postprocess/distribution_compare.py`; activates `validation_histogram_metric`);
  R5 matched-time + sim-end upper-bound both reported.
- **E5** — derived diagnostics `t×` / `Π(t)` / regime / total-strip reachability
  from the v7 arrays (`postprocess/derived_diagnostics.py`).

**Prereqs:** Slice X (v7 `E_int`) gates E2/E5; L+P+K gate E5; G+Q+K+U gate E2.
E1/E3/E4 are buildable now on synthetic v7 checkpoints; real runs wire in at R2.

### Phase F — Calibration campaign + production switch

> **Detailed plan: `PHASE_F_IMPLEMENTATION_PLAN.md`** (slices F1–F6).

R2 expanded into the campaign that absorbs everything postponed to this stage
(κ+picture co-fit, the f_int/f_ret/τ bounded scalars, the 0.80→2.70 eV switch,
identifiability + regime reporting, total-strip runs, and all scripts/figures Phase
E deferred). Composes **only accepted A–E modules** — no new physics, no mechanism
re-litigation; only calibration *values* move via config knobs. Scope decisions
locked (user, 2026-06-29): (1) **staged** clean-fallback calibration — co-fit
κ+picture, pin f_int @ fixed τ, then a τ-insensitivity sweep, f_ret via the 9/18 Å
contrast (full grid is a documented escalation); (2) **larger-N single-seed** runs
(N≈500) for a stable histogram; (3) **build the total-strip path now**, R6
re-extraction + p↔κ split stay documented triggers.

- **F1** — campaign harness `scripts/tier2_common.py` (`build_biphasic_cfg` +
  run-dir convention, mirroring `tier1a_common`).
- **F2** — `scripts/gen_tier2_runs.py` staged run-matrix generator (0.80 eV first,
  9 Å + 18 Å, ion→relaxation pipeline; mirrors `gen_tier1a_runs.py`).
- **F3** — scoreboard `scripts/post_processing/tier2_size_distribution_table.py`
  composing E1/E2/E4/E5 → per-run W₁ (matched + upper-bound) + diagnostics; reported,
  not adjudicated (mirrors `tier1a_rmse_table.py`).
- **F4** — `tier2_identifiability_report.py`: κ/picture landscape, f_ret 9/18 Å
  contrast, τ tail-insensitivity, s↔κ, regime determination.
- **F5** — 2.70 eV production switch (gated on 0.80 eV landing) + total-strip
  secondary runs; R6 / p↔κ as documented conditional triggers.
- **F6** — overlay figures (size-dist overlay, Π/t× regime, W₁-vs-κ), `SHOW_FIGURE`
  gated; no figures in pytest.

**Prereqs:** all of A–E + the Phase-D bridge (the 0.80 eV cross-check that gates the
production switch). No new `SimConfig` physics fields — uses the A–E knobs.

---

## 5. Config contract (new `SimConfig` fields — MASS §11)

All new unless noted. Each lands with its owning slice; declared-but-unread fields follow
the rule-2 exception convention until activated.

- **Ladder (L):** `dissociation_ladder` (Form U), `ladder_steepness` κ,
  `ladder_electronic_picture ∈ {statistical_mixture, x2_only, cooling_relaxed}`,
  `solv_struct_asymptote` |S(N)|, `electrostriction_binding`.
- **Cooling (K):** `internal_energy_cooling_tau_ps` τ.
- **E_int (U):** `internal_energy_partition_fraction` f_int,
  `internal_energy_retained_fraction` f_ret, `internal_energy_initial_eV`
  (retired→optional override). `coulomb_available_eV` (exists; scenario-keyed).
- **Pickup (P/ρ):** `pickup_rate_coefficient` λ_0 (rename/repurpose `mass_rate_coefficient`),
  `pickup_rate_form` (exists as `mass_rate_form`), `pickup_occupancy_cap`,
  `pickup_occupancy_exponent` p, `helium_density_profile`.
- **Evaporation (Q):** `evap_rate_prefactor_per_ps` ν, `evap_rrk_dof` s (guarded s≥1),
  `evap_gate_onset_eV` (retired→optional override).
- **Integrator (G/X):** `mass_jump_velocity_reset=momentum_conserving` (forbid
  `label_only` in production), `he_capture_velocity=at_rest`, `one_mass_event_per_step=true`,
  `jump_o_step_ordering=jump_then_O`.
- **Comparison (W):** `validation_histogram_metric=wasserstein` (exists).
- **Naming reconcile:** MASS §11 names the production value `biphasic_energy_gated`;
  `config.py` currently has `biphasic`. Pick one at build (recommend keeping `biphasic`,
  alias in docs) — record in `drag_migration_log_tier2.md`.

---

## 6. Cross-cutting gates and out-of-scope guards

- **5-term invariant** is the correctness detector from X onward (R4 double-count guard:
  drag→E_dissip; binding release→split E_int/E_dissip).
- **RNG draw order is locked** and documented (forbidden-list: do not change draw order
  without explicit approval).
- **Config-load guards:** `s≥1`; mass↔coefficient pairing (production trips it →
  `allow_inconsistent_mass_pairing=True`, R6); scenario-stamped E_avail.
- **Out of scope (fails review if it leaks in):** the **noise model** (Tier 3 — stays
  inert behind its enum); any change to neutral propagation, the drag law (Tier-0 locked),
  checkpoint draw order, or default scope beyond the schema bump.
- **OQ1 (electronic picture provenance)** stays open → carry **both** pictures as the
  co-fit free knob; do **not** cite E_bind as independent support for the mixture. Update
  if the author contact resolves it.

---

## 7. Testing methodology (per slice)

Mirror Tier-1a: pure functions (L, K, U, ρ) against closed-form oracles; stochastic units
(P, Q) with seeded/mock RNG and distribution/moment checks; composed slices (G) on analytic
limits + the invariant before any real run; comparison slices (H, W, D2) on tiny synthetic
checkpoints — **no figures or production-sized checkpoints in pytest**. Justify tolerances
(tight analytical; sample-size-based Monte Carlo; FD-consistent where relevant). Run the
narrowest relevant test first, then broaden; interpreter:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q`.

---

## 8. Acceptance criteria

**Per slice:** the slice's own suite green with all other slices mocked.

**Generative driver (G):** constant-mass / `fixed` regression unaffected; force-free run
conserves the 5-term invariant; one-event-per-step + locked draw order; E_int gate behavior
correct.

**Bridge (Z):** emergent mean n(t) reproduces the anchored 21→19→14 staircase within the
loose early-window tolerance and t× near the GAH25 prior at 0.80 eV — **reported**, not an
automated verdict.

**Tier-2 observable (H/W/D2):** size-distribution + per-fragment histograms extracted on
integer-n support; Wasserstein vs `integrated_i_he_abundance.csv` and `vmi_summary` computed
with the R5 matched-time caveat carried in the report.

**Campaign (R2):** the κ/picture co-fit + bounded sweep runs to completion at **0.80 eV**
then **2.70 eV**, emitting a Wasserstein scoreboard + regime-determination report with
identifiability notes. No automated pass/fail on fidelity — adjudication is reported and
left to the user.

**Out-of-scope guard:** any code path reading a noise amplitude, changing the Tier-0 drag
law, or changing the locked RNG draw order **fails review**.

---

## 9. Risks / notes

- **Doc drift — RESOLVED (2026-06-29).** `CLAUDE.md` now records **Tier 1a — DELIVERED**
  and **Tier 2 — ACTIVE (current goal)** (Project Snapshot + validation hierarchy),
  consistent with the Tier-1a log + plan (all slices delivered, suite all-green). No
  reconcile outstanding.
- **Calibration load (HIGH).** One observable (size distribution) carries 8+ quantities
  (CALIBRATION_MAP). Whether it *separates* them is the key open empirical question — the
  campaign (R2) must report identifiability, not assume it.
- **R5 truncation.** Terminal n at 20 ps is an upper bound; the comparison must be
  matched-time or post-relaxation. Carried explicitly in W.
- **Mechanism is locked, values are not.** No mechanism re-litigation in slices; only
  calibration values move, and only via config knobs the slices expose.
