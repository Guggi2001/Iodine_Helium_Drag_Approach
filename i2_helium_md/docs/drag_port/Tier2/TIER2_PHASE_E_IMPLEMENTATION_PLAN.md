# Tier 2 — Phase E Implementation Plan (Tier-2 Observable / Comparison Layer)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger. Equations are the *locked* formulations
> from `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`; the module descriptions
> are *interface contracts*, not implementations.
>
> **Parent:** `TIER2_IMPLEMENTATION_PLAN.md` §3–§4 (Phase E = slices E1–E5,
> re-sliced from the original H/W/D2 headline).
> **Entry docs:** MASS doc §2.1–2.2 (the observable), §4 (channels), §6.11
> (regime / `Π` / `t_×`), §R5 (truncation), §7 (schema), §11 (config);
> `CALIBRATION_MAP.md` (parameter classes); `TIER2_PHASE_D_BRIDGE_FINDINGS.md`
> (the delivered 9 Å / 0.80 eV pinned-point behaviour Phase E's oracles must
> respect).

---

## 0. Status and intent

Phase E is the layer that turns a finished generative (`biphasic_energy_gated`)
run into the **one observable that arbitrates Tier 2**: the terminal **I⁺Heₙ
integer-n size distribution**, scored against the experimental reference
`data/reference/integrated_i_he_abundance.csv` with a Wasserstein metric. It also
reconstructs the derived diagnostics (`t_×`, `Π(t)`, regime label) that explain
*why* a given distribution came out.

Phase E ships and is tested **entirely on synthetic checkpoints + the real
reference CSV**; production runs are wired in only at Phase F (R2). It composes
only **accepted** upstream modules.

**Confirmed scope decisions (2026-06-29):**
1. **R5 truncation → build a post-ejection relaxation stage.** Terminal `n` at
   20 ps is only an *upper bound* (MASS §R5: the cascade continues to hundreds of
   ps; loosely-bound outer He keep evaporating). Phase E adds a reduced relaxation
   driver that propagates the energy-gated cascade forward to the experimental
   timescale before the size distribution is read off.
2. **Size distribution only.** The experimental velocity references
   (`data/reference/vmi_summary/*.csv`) are *aggregate* I⁺He / I⁺-gas
   distributions, not per-fragment per-`n`. Phase E therefore drops velocity
   comparison entirely and commits to the integer-n size distribution as the sole
   observable. (The original Slice-H per-fragment velocity histograms are cut.)
3. **Pure functions only.** No CLI / figure / report scripts; all orchestration
   and figures are deferred to the Phase F campaign (R2). (The one sanctioned
   non-code artifact is the E3 README data-contract entry, decision R3 below.)

The locked mechanism is **not** re-litigated: the relaxation stage is the locked
Q/K/U channels composed with pickup and drag switched off — exactly the
post-ejection regime MASS §R5 describes.

**Pre-build refinement pass (2026-07-03).** The plan below was refined against
the delivered Phases A–D at HEAD (full code-surface audit + MASS / CALIBRATION_MAP
cross-reference + the Phase-D bridge findings; record:
`drag_migration_log_tier2.md`, "Phase E — pre-build refinement decisions").
Three decisions locked:

- **R1 (user). E2 composes the delivered `biphasic_step` verbatim** under a
  λ₀ = 0 "relaxation view" of the config — **not** a separate shed-only reduced
  loop. Pickup is structurally inert at λ₀ = 0 but its RNG draw is still
  consumed, so the Slice-X frozen two-draw stream contract is preserved
  byte-identically (`config.py`'s λ₀ = 0 advisory already names this the
  relaxation stage's "close cousin"). Drag drops out via a zero-γ BAOAB closure.
  The earlier "draws **only** shed Bernoullis" wording is **superseded** — it
  contradicted the delivered unconditional-draw channel contract.
- **R2 (user). E5 lives in the generalized Phase-D module:**
  `git mv postprocess/bridge_diagnostics.py → postprocess/derived_diagnostics.py`
  (module named for its concern, Quality Principle 6; matches the parent / Phase-F
  plan naming; rule 1 — no second t_×/Π implementation). The delivered helpers
  (`mean_shell_count`, `crossing_time_ps`, `regime_parameter`) move unchanged;
  the Phase-D importers (`scripts/post_processing/tier2_bridge_report.py`,
  `tests/test_bridge_diagnostics.py` → `test_derived_diagnostics.py`) are updated
  mechanically.
- **R3 (recommendation applied; user idle on the ask — Phase-D precedent).**
  `integrated_i_he_abundance.csv` has **no documented provenance** (no exporter
  under `data/reference/scripts/`, no `data/reference/README.md` entry; the only
  repo reference is the consumer `plotting_histogram.py`). E3 adds a README
  **data-contract entry** (columns, units, normalization, consumer) with
  provenance marked **"to be completed"**, and the gap is carried as an open
  item in the log. Does not block E3.

---

## 1. Hard prerequisites and test strategy

Phase E binds to these **accepted** upstream modules — **all delivered at HEAD**
(Phases A–D complete, full suite 1840 green as of 2026-07-03), so every slice is
buildable now; the "gated on Phase C" phasing in the original plan is obsolete.

| Prereq | Provides | Needed by |
|---|---|---|
| **X** (Phase C) — checkpoint v7 + 5-term invariant | per-step `E_int_eV (2N,T)`; extended `ion_ledger_closure` | E2, E5 |
| **L** (Phase A) — ladder `D_0(n)` + cumulative `Σ_{i≤n}D_0` | the self-bound gate threshold | E2, E5 |
| **K** (Phase A) — Newton cooling K2 + occupancy-resolved `E_∞(N)` | cooling drain; total-strip reachability; `e_bind_pair_eV` fold | E2, E5 |
| **U** (Phase B/A) — `E_int` budget rules (S1/K1) | K1 drain on shed | E2 |
| **Q** (Phase B) — energy-gated RRK evaporation | the cascade rate engine | E2 |
| **P/ρ** (Phase B) — `λ_attach(rho_ratio, n, ...)` | the `Π` order parameter (delivered in `regime_parameter`) | E5 |
| **G** (Phase C) — `biphasic_step` generative seam | the per-step mass subsystem E2 reuses verbatim | E2 |
| **Z** (Phase D) — `bridge_diagnostics.py` helpers | the E5 generalization seed (R2) | E5 |

**Test strategy (mirrors Tier-1a).** Pure functions against closed-form oracles
(E1, E3, E4); stochastic units against seeded RNG with moment checks + real-PCG64
stream-consumption parity (E2 — the pickup `TestRNGConsumptionParity` precedent);
the **5-term invariant** as the cross-cutting correctness gate (E2). **No figures
or production-sized checkpoints in pytest** — tiny synthetic checkpoints only
(fixture precedents: `tests/test_checkpoint.py::_make_ion_checkpoint` for
hand-built v7 arrays; `tests/test_biphasic_step.py::_tiny_neutral` /
`_biphasic_driver_cfg` for driver-produced checkpoints).

---

## 2. Slices E1–E5

### E1 — Terminal I⁺Heₙ size-distribution extractor (integer-n)

- **Purpose.** From a finished generative run (or a relaxation result, E2), build
  the simulated size distribution on integer-n support — the direct counterpart of
  the experimental abundance reference.
- **Encoded form.** Per ion atom `i`, terminal shell count
  $n_i = \texttt{n\_shell}[i,-1]$ (sim-end) **or** the relaxed terminal `n` from E2.
  Histogram on integer support $n = 0,\dots,n^*$ ($n^*=21$):
  $$
  \mathrm{counts}(n)=\sum_i \mathbf 1[n_i=n],\qquad
  \mathrm{fraction}(n)=\frac{\mathrm{counts}(n)}{\sum_m \mathrm{counts}(m)} .
  $$
  Label each rung by $m(n)=\texttt{complex\_mass\_amu}(n)$ to cross-check against
  the reference mass windows. The native Poisson spread (MASS §2.2) is preserved —
  **no binning of a continuous endpoint**.
- **Dtype note (audit).** `n_shell` is stored as an **int-valued float** array
  (v7 contract). E1 validates integer-valuedness (`rint` residual = 0, fail-loud
  otherwise — a fractional entry means upstream corruption) and casts to `int`.
- **Support note (audit).** Simulated support is $0\dots21$: runs start at
  $n_0 = n^* = 21$ and the Langmuir cap forbids $n > n^*$, so **n = 21 is a
  legal simulated outcome** (the bridge run produced per-ion terminal
  $n \in [18,21]$). The experimental reference stops at $n = 20$; the mismatch is
  handled by E4's zero-filled union support, not by clipping in E1.
- **Reuse (exact paths).**
  - `i2_helium_md/physics/shell_schedule.py :: complex_mass_amu(n)`
    ($m(n)=126.90+4.0026\,n$ amu) — rung labelling / reference-window match.
  - `i2_helium_md/simulation/run_directory.py :: RunDirectory.load_ion/load_cfg`;
    `i2_helium_md/simulation/checkpoint.py :: IonCheckpoint.n_shell (2N,T)`,
    `mass_final_kg`.
- **New surface.** `postprocess/size_distribution.py ::
  compute_terminal_shell_distribution(source, *, n_max=21) -> ShellDistribution`
  (frozen: `n_values (Nn,)` int, `counts (Nn,)`, `fraction (Nn,)`, `source`).
  Accepts **either** a raw `IonCheckpoint` (sim-end = R5 upper bound) **or** a
  `RelaxationResult` (E2, matched/relaxed time — duck-typed on the terminal
  `n_shell` column, so E1's sim-end mode lands first and the relaxed mode is a
  type admission, not new logic).
- **Oracle / tests.** Hand-counted tiny checkpoint (e.g. `n=[0,0,1,2,2,2]`
  → `fraction` matches by hand); integer support exact; a terminal `n = 21` entry
  is counted (not clipped); int-valued-float cast + fractional-entry rejection;
  monotone-falling envelope reproduced on a synthetic stream; empty/degenerate
  ensemble raises loudly.
- **Acceptance.** Matches hand-count; integer support incl. n = 21; both input
  modes accepted; `fraction` sums to 1; fractional `n_shell` fails loudly.

### E2 — Post-ejection relaxation stage (R5 mitigation)

- **Purpose.** Propagate the **mass subsystem** forward past the 20 ps ion stage
  to the experimental timescale, so the size distribution is read at *matched
  time*, not at the truncated upper bound (MASS §R5).
- **Composition (decision R1 — reuse, not re-composition).** Per relaxation step
  the stage calls the **delivered `biphasic_step` verbatim** under a
  **relaxation view** of the run config:
  `replace(cfg, pickup_rate_coefficient=0.0, dt_ion=dt_relax)` (constructed
  internally; not re-`validate()`d — the λ₀ = 0 advisory is the sanctioned
  evaporation-only limit, `config.py::check_biphasic_config`). Then:
  1. `biphasic_step` runs the locked K2 → gate → RRK → K1 sequence: **K2** drains
     `E_int` toward `E_∞(N)` (Slice K, inside the step); the self-bound gate
     $E_\text{int}>\sum_{i\le n}D_0(i)$ suppresses; else RRK
     $k=\nu\,(1-D_0(n)/E_\text{int})^{\,s-1}$ ($n{=}1$ direct dissociation
     $k=\nu$); on fire: $n\to n-1$, K1 $E_\text{int}\mathrel{-}=D_0(n)$,
     **cold-shed momentum reset** (`mass_jump.cold_shed` via
     `evaporation_step_components`). The pickup channel runs too but is
     **structurally inert** (λ₀ = 0 ⇒ `P_attach ≡ 0`); its draw is still
     consumed.
  2. Translation: rebuild the BAOAB closure at $m^+$ with a **zero-γ**
     `gamma_fn` (decay = 1, `dE_dissip ≡ 0` — drag off with the delivered
     bookkeeping intact) and run `baoab_propagation_step`
     (`relaxation_forces="coulomb"`, default: the existing
     `make_ion_accel_fn` conservative field — residual inter-ion Coulomb +
     droplet solvation, the latter vanishing outside the droplet), **or** skip
     the force evaluation entirely (`"free_flight"`: ballistic positions,
     `E_pot` held between sheds).
  3. The driver applies the same **`e_bind_pair(n)` E_pot fold** as the ion
     stage's biphasic arm (rule 1 — same call, `solvation_cooling.e_bind_pair_eV`),
     so a shed shifts `E_pot` by exactly $+D_0(n)$ and the 5-term closure holds
     in both force arms.
- **Decoupling fact (audit).** With λ₀ = 0 and γ = 0 the mass subsystem
  (K2 → gate → RRK → K1) is **independent of translation** — positions enter only
  through the dead density gate. Translation is retained solely for ledger and
  asymptotic-state fidelity; this is what makes the `free_flight` arm exact for
  the observable, not an approximation of it.
- **Seeding.** Initial state = `ion_state_from_checkpoint_column(ion, -1)` (the
  final stored column carries positions, velocities, `mass_kg`, `E_int_eV`,
  `n_shell`, and all ledger columns).
- **RNG (stream extension, not change).** The relaxation stage uses its **own
  generator**, seeded `SeedSequence((cfg.seed, RELAXATION_STREAM_KEY))` with a
  fixed module-level key — the ion-stage stream is untouched (the Slice-X freeze
  is *extended by a new stage*, never re-ordered). Within the stage the frozen
  two-draw order holds verbatim (evaporation draw, then the inert pickup draw,
  both unconditional, both consumed every step until termination).
- **Termination.** Run to `relaxation_time_ps` **or** stop early when every
  fragment is **frozen**: $n_i = 0$ or $E_{\text{int},i} < D_0(n_i)$ (then
  $k \equiv 0$ forever, since `E_int` is monotone non-increasing under K2 + K1
  with pickup off). Record `freeze_flags (2N,)` and the time actually reached.
- **Artifact (no schema bump — reuse the checkpoint machinery).** The stage emits
  a **bona fide v7 `IonCheckpoint` covering the relaxation window** (all required
  fields are available: trajectories via the ion-stage **auto-stride** convention
  (`ion.py::_decide_stride_ion` size budget), `n_shell`/`E_int_eV`/ledger columns
  from the loop, pass-throughs for `droplet_radii_angstrom`/`b_ion_outside`,
  zero `number_of_collisions`, `mass_scenario="biphasic"`), saved as
  **`relaxation.npz`** beside `ion.npz` (self-describing run directory), via the
  existing `save_ion_checkpoint`/`load_ion_checkpoint` — zero new I/O code, and
  `ion_ledger_closure` applies to it unchanged.
  `run_relaxation_stage(ion, cfg) -> RelaxationResult` returns a thin frozen
  wrapper: `checkpoint` (the relaxation-window `IonCheckpoint`),
  `freeze_flags (2N,)`, `time_relaxed_ps`, `terminal n` view for E1.
- **Config (new; all opt-in; default off keeps default scope unchanged).**
  - `relaxation_stage_enabled: bool = False`.
  - `relaxation_time_ps: Optional[float] = None` — **required-when-enabled, no
    default** (fail-loud): MASS §R5 gives no sourced $t_\text{exp}$, only
    "hundreds of ps"; the concrete value is a Phase-F campaign choice. The
    freeze early-exit bounds the cost of generous values.
  - `relaxation_dt_ps: Optional[float] = None` → defaults to `dt_ion` (0.01 ps).
    Guard: refuse $\nu\cdot dt_\text{relax} > 0.1$ — since $k \le \nu$, this
    bounds the one-event-per-step bias at $(k\,dt)^2/2 \lesssim 0.5\%$ (MASS
    authorizes no larger step; at the default, $\nu\,dt \approx 0.024$).
  - `relaxation_forces ∈ {"coulomb" (default), "free_flight"}` — enum reject-arm
    guard for anything else.
  - `check_relaxation_config` wired into `validate()` (the `check_biphasic_config`
    pattern): no-op when disabled; when enabled require `mass_scenario ==
    "biphasic"`, `relaxation_time_ps` set and > 0, the dt guard, the enum guard.
- **Oracle / tests (corrected against the delivered bookings).**
  - **5-term invariant closes** on the relaxation window (the cross-cutting gate).
  - **Per-channel sharp oracle (supersedes the wrong "E_dissip constant"):** with
    γ = 0 and λ₀ = 0, $\Delta E_\text{dissip}(t)$ **equals the cumulative K2
    cooling drain exactly** (the drag and pickup-bath contributions are
    identically zero; the K2 drain books to `E_dissip` at
    `ion_propagation_step.py` — delivered Slice-G booking).
  - A fully self-bound cold input ($E_\text{int}<D_0(n)\ \forall$) relaxes with
    **zero sheds** (`n` unchanged) while **both draws are still consumed** —
    real-PCG64 post-state parity (the pickup `TestRNGConsumptionParity`
    precedent).
  - A hot synthetic input sheds monotonically until freeze; no avalanche
    (per-step shed $\lesssim \nu\,dt$). *(Synthetic-hot only: the bridge showed
    the real 9 Å pinned point is a ~0.7-shed burst frozen by ~8 ps — at that
    point the relaxation stage adds ≈ 0 sheds, a characterization to record,
    not a target.)*
  - Matched-time terminal `n` distinct from (and $\le$) the sim-end upper bound
    on the hot synthetic case; `coulomb` vs `free_flight` arms give the
    **identical shed sequence** under the same seed (the decoupling fact made
    testable).
  - Seeded RNG, tiny synthetic input throughout.
- **Acceptance.** Reduced driver conserves the 5-term invariant; ΔE_dissip ≡
  cumulative K2 drain; gate + no-avalanche + draw-consumption parity hold;
  freeze-termination reached; matched-time terminal `n` ≤ sim-end on the hot
  case; the artifact round-trips through `load_ion_checkpoint`.

### E3 — Experimental abundance reference loader

- **Purpose.** Load `integrated_i_he_abundance.csv` (the Tier-2 arbiter) under the
  same validate-early contract as the HeDFT loader.
- **Reuse (mirror the contract).** `postprocess/hedft_loader.py ::
  load_hedft_trajectory` — frozen-dataclass return, **order-independent
  membership** header match (missing *and* extra columns raise; mirror this
  variant, not the module's stricter exact-tuple sibling
  `load_smoothed_speed_reference`), loud `FileNotFoundError`/`ValueError`,
  `source_path = p.resolve()` retained.
- **Data contract (verified against the real file, 2026-07-03).** Columns
  `n, label, massCenter_u_per_e, massWindowLower_u_per_e,
  massWindowUpper_u_per_e, ionCounts, ionPercent`; 21 rows, `n` contiguous
  $0\dots20$ ($n=0$ = bare I⁺, `label` `I^+`…`I^+He_20`); `massCenter` step =
  4.0026 u/e (He); `ionPercent` sums to 100.0000.
- **New surface.** `postprocess/abundance_loader.py ::
  load_he_abundance_reference(path) -> HeAbundanceReference` — **function name
  pinned**: the Phase-F plan (F3 reuse row) already binds to
  `load_he_abundance_reference`. Frozen: `n (int,Nn,)`, `label`,
  `mass_center_u (Nn,)`, `ion_counts (Nn,)`, `ion_fraction (Nn,)` =
  `ionPercent` normalized **by its own sum** (sums to exactly 1.0),
  `source_path`.
- **Validation.** Columns present (order-independent, extras refused); `n`
  integer, contiguous, ascending from 0; `ionCounts ≥ 0`; `ionPercent` sums to
  100 within a stated tolerance (loud otherwise); raise loudly on all failures.
- **Provenance (decision R3).** E3 adds a `data/reference/README.md` entry for
  the file: columns/units/normalization, the consumer script, and **provenance
  marked "to be completed"** (no exporter or measurement record exists in-repo).
  Carried as an open item in the log; the user supplies measurement IDs / the
  producing script when available.
- **Oracle / tests.** Loads the real CSV (values spot-checked: n=0 → 43.52 %,
  n=20 → 0.2645 %); validation fires on corrupted synthetic CSVs (missing column,
  extra column, non-contiguous n, negative counts, bad percent sum); normalized
  `ion_fraction` sums to 1.
- **Acceptance.** Real CSV loads; every validation arm fires on a synthetic
  corruption; `ion_fraction` sums to 1; README entry present.

### E4 — Wasserstein size-distribution comparison (integer-n, pure numpy)

- **Purpose.** Score simulated (E1) vs experimental (E3) size distributions with
  the config-selected metric on integer-n support.
- **Encoded form.** On unit-spaced integer support the 1-D $W_1$ reduces to the
  CDF-gap sum:
  $$
  W_1=\sum_n \big|\,F_\text{sim}(n)-F_\text{ref}(n)\,\big|,
  $$
  with $F$ the cumulative fraction over the **union** support (zero-filled where a
  side is absent — this is where the sim's legal $n=21$ meets the reference's
  $0\dots20$ support). Units: He atoms (unit rung spacing). **Pure numpy** — a
  1-D integer-support $W_1$ is exact and trivial; note this is a *style/precision*
  choice, not a dependency constraint (scipy is already a hard package dependency
  via `erf`/`curve_fit`). The optional cross-check against
  `scipy.stats.wasserstein_distance` lives **in tests only**.
- **Config activation (retires the rule-2 carry).** `config.py ::
  SimConfig.validation_histogram_metric` is the three-arm Literal
  `{"wasserstein","chi2","ks"}`, default `"wasserstein"`, currently
  declared-but-unread. E4 adds the dispatching entry point — the field's first
  **physics-live** reader — and gives the unbuilt `chi2`/`ks` arms
  **point-of-use `NotImplementedError` refusals** (the established
  unbuilt-enum-arm convention). Removed from the rule-2 exception table at the
  E4 build.
- **New surface.** `postprocess/distribution_compare.py ::`
  - `wasserstein_integer_support(sim: ShellDistribution, ref:
    HeAbundanceReference) -> float` — aligns union supports; raises loudly on
    disjoint/empty support;
  - `compare_size_distributions(sim, ref, *, metric: str) -> float` — the
    config-dispatch activation point (`metric = cfg.validation_histogram_metric`
    at the caller; `chi2`/`ks` refuse).
- **R5 caveat (carried in every comparison record).** Score at **matched time**
  (post-E2 relaxation) **and** report the **sim-end upper-bound** value; never
  over-read absolute `n`. E4 computes both scores; the *pairing record*
  (`W1_matched` + `W1_simend_upper` per run) is assembled by the Phase-F F3
  scoreboard — the E4 API must make both trivially computable (it does: call it
  on the E1 extraction of each source).
- **Oracle / tests.** Hand-computed micro-case (two 3-point distributions → known
  $W_1$); identical distributions → 0; shifting all mass by one $n$ → exactly the
  shift magnitude; sim mass at $n=21$ vs zero-filled reference handled; disjoint
  support raises; `chi2`/`ks` refuse; scipy cross-check (import-guarded, tests
  only).
- **Acceptance.** Matches hand-computed $W_1$; metric selected via the config
  field (carry retired); union support incl. n = 21 exact; refusal arms fire.

### E5 — Derived diagnostics (`t_×`, `Π(t)`, regime label)

- **Purpose.** Post-hoc reconstruct, at **zero schema cost** from the v7 arrays,
  the quantities that *explain* the size distribution (MASS §6.11).
- **Module (decision R2).** `git mv postprocess/bridge_diagnostics.py →
  postprocess/derived_diagnostics.py`; update the two Phase-D importers
  (`tier2_bridge_report.py`; `test_bridge_diagnostics.py` →
  `test_derived_diagnostics.py`) mechanically; delivered helpers
  (`mean_shell_count`, `crossing_time_ps`, `regime_parameter`) move **unchanged**
  (their contracts — per-ion `min{t: E_int < Σ(n(t))}` via the driver's own
  `is_self_bound` incl. the `gate_onset_eV` passthrough; Π via `lambda_attach`
  with ρ re-derived through the shared steepness resolver — are already the E5
  formulations).
- **Absorbed deferrals (recorded owners: Phase-E D2/E5, Slice-Z review pass).**
  1. Promote the private `simulation.ion::_drag_gate_steepness` import to a
     **public helper** (public name on `simulation.ion`, single source; the
     module and any future caller import it without the underscore contract).
  2. Add the **`postprocess/__init__` exports** for the Phase-E surfaces
     (E1/E3/E4/E5 public functions + dataclasses).
  3. Add **shape / monotonic-`time_ps` guards** to the thin delivered helpers.
- **New encoded forms (added to the module).**
  - **Crossing time $t_\times$:** delivered (`crossing_time_ps`), plus an
    ensemble `t_cross_summary` (median / spread / all-agree flag — the
    deterministic pre-crossing dynamics make all-agree the wiring check).
    Sanity band: flag only if $t_\times<1$ ps or $>15$ ps; GAH25
    $t_0\approx5$–6.5 ps is a **±factor-2 prior, not a threshold** (a factor-10
    miss is the genuine flag; Na⁺ number, I⁺ is Rb⁺-like).
  - **Order parameter $\Pi(t)=\lambda(n)\,f_\text{ret}\,\tau$:** delivered
    (`regime_parameter`). Interpretation per MASS §6.11: $\Pi>1$ shedding
    persists, $\Pi<1$ freeze, $\Pi\to0$ at exit guarantees termination.
    **Condition-specific expectation (bridge finding):** at 9 Å / 0.80 eV the
    run sits deep on the freeze side — Π = 0 at gate-open, max mean Π ≈ 0.005 —
    so a computed $\Pi>1$ *there* is itself a wiring/units flag; the reading
    must **not** be carried to 2.70 eV.
  - **Regime label:** shell-retaining (self-binds early → moderate terminal `n`)
    ↔ total-strip (never self-binds in window → Calvo24 limit). Concrete
    reporting rule (a documented convention, **not** a physics knob — no config
    field): `total_strip` iff the majority of ions have no finite $t_\times$ in
    the window **or** the median terminal $n \le 1$; else `shell_retaining`.
    Expected at the 9 Å pinned point: `shell_retaining` (bridge finding).
  - **Total-strip reachability:** $E_\infty(N)\to0$ as $N\to0$
    (occupancy-resolved; OQ6 resolved) — evaluate `e_infinity_eV(0)` under the
    run's picture/κ. *Annotation (K/U precedent): this is structurally true for
    the delivered Slice-K form — a consistency confirmation, not an independent
    anchor.*
- **New surface.** `postprocess/derived_diagnostics.py ::
  reconstruct_diagnostics(ion, cfg) -> Diagnostics` (frozen: `t_cross_ps (2N,)`,
  `t_cross_summary`, `Pi_t (2N,T)`, `regime_label`, `total_strip_reachable:
  bool`, `sanity_flags`) — **field names pinned**: the Phase-F F3/F4 scoreboard
  binds to `t_×`, `Π`, `regime_label`, `total_strip_reachable`. Composes the
  three delivered helpers; re-derives nothing.
- **Reuse (exact paths).** Slice L cumulative ladder (through
  `evaporation.is_self_bound` — the delivered gate-unification fix); Slice P
  `lambda_attach(rho_ratio, n, *, lambda0, ...)` (through `regime_parameter`);
  Slice K `e_infinity_eV` (reachability);
  `simulation/checkpoint.py :: IonCheckpoint.E_int_eV (2N,T)`, `n_shell (2N,T)`;
  `postprocess/energy_balance.py :: ion_ledger_closure` (5-term) for invariant
  context.
- **Oracle / tests.** Delivered suites move with the rename (regression-locked);
  new: `t_cross_summary` on an all-agree synthetic + a constructed disagreement;
  regime label flips across a constructed early-bind vs never-bind pair;
  reachability true under all three pictures; sanity flags fire on absurd inputs
  ($t_\times$ out of band; Π > 1 on a freeze-side construction); the new
  shape/monotonic-`time_ps` guards fire.
- **Acceptance.** All four diagnostics reconstructable with zero schema cost;
  rename lands with Phase-D report/tests green; deferral items 1–3 closed;
  sanity flags fire; regime label emitted with the documented rule.

---

## 3. Dependency graph and build order

All upstream prereqs (Phases A–D) are delivered — every slice is buildable now.
Internal order:

```
E5a (git mv bridge_diagnostics → derived_diagnostics; deferral items)   ── first (pure rename, keeps later diffs clean)
E3 (loader)  ──┐
E4 (W1)      ──┼──►  E1 (sim-end mode)  ──►  E2 (relaxation)  ──►  E1 relaxed-input admission
E5b (reconstruct_diagnostics + summary/label/reachability/flags) ──►  (any time after E5a)
```

- **E3 / E4 / E1(sim-end) are parallelizable** (pure functions on synthetic v7
  checkpoints + the real CSV).
- **E2** composes only delivered code; its `RelaxationResult` unlocks E1's
  relaxed-input admission (a type admission, not new logic).
- **E5a before E5b** (rename first so additions land in the final module name).
- Real runs are wired in **Phase F (R2)** — Phase E is tested entirely on
  synthetic checkpoints + the real reference CSV.

---

## 4. Config contract (new / activated)

| Slice | Field | Default | Status |
|---|---|---|---|
| E2 | `relaxation_stage_enabled` (bool) | `False` | new, opt-in |
| E2 | `relaxation_time_ps` (`Optional[float]`) | `None` — **required-when-enabled** | new (no sourced $t_\text{exp}$; Phase-F picks the value) |
| E2 | `relaxation_dt_ps` (`Optional[float]`) | `None` → `dt_ion` | new; guard $\nu\cdot dt \le 0.1$ |
| E2 | `relaxation_forces ∈ {coulomb, free_flight}` | `coulomb` | new; enum reject-arm |
| E4 | `validation_histogram_metric` (`{"wasserstein","chi2","ks"}`) | `"wasserstein"` | **activated** (first physics-live reader; `chi2`/`ks` = point-of-use refusals; rule-2 carry retired) |

`check_relaxation_config` (new, wired into `validate()`): no-op when
`relaxation_stage_enabled=False`; when enabled requires
`mass_scenario == "biphasic"`, `relaxation_time_ps` set and positive, the
$\nu\cdot dt$ bound, and the forces enum. Default
`relaxation_stage_enabled=False` keeps the default simulation scope unchanged
(forbidden-list compliant). New fields follow the rule-2 declared-but-unread
convention until E2 activates them (they land *with* E2, so they are born live).

---

## 5. Cross-cutting gates

- **5-term invariant** is the correctness detector for E2, with the corrected
  per-channel form: γ = 0 and λ₀ = 0 ⇒ $\Delta E_\text{dissip}$ ≡ the cumulative
  K2 cooling drain (drag and pickup-bath contributions identically zero); the
  shed bookings close through the `e_bind_pair` E_pot fold +
  `E_mass_transfer` + K1, exactly as in the ion stage.
- **RNG:** the Slice-X frozen two-draw stream is **reused verbatim inside the
  stage** (pickup draw consumed though inert) and the stage runs on its **own
  seeded generator** (`SeedSequence((cfg.seed, RELAXATION_STREAM_KEY))`) — the
  locked ion-stage draw order is extended by a new stage, never re-ordered. The
  superseded "shed-only draws" wording is retired (refinement R1).
- **Integer-n support** everywhere; **Wasserstein pure-numpy** (style choice;
  scipy cross-check tests only).
- **R5 matched-time caveat** carried in every comparison record (E4): report both
  the relaxed/matched-time score and the sim-end upper bound (paired by F3).
- **Phase-F interface names are pinned:** `load_he_abundance_reference` (E3);
  matched + upper-bound W₁ computable per run (E4); `t_cross_ps` / `Pi_t` /
  `regime_label` / `total_strip_reachable` / `sanity_flags` (E5);
  `relaxation_stage_enabled` (E2).
- **Out of scope (fails review if it leaks in):** velocity comparison (cut per
  the scope decision); any noise-amplitude read (Tier 3 stays inert); changes to
  the Tier-0-locked drag law, neutral propagation, ion-stage RNG draw order, or
  default scope beyond the opt-in relaxation flag; CLI/figure scripts (deferred
  to R2; the E3 README data-contract entry is the one sanctioned doc artifact);
  any new mass-mechanism physics (E2 composes `biphasic_step`, it does not
  modify it).

---

## 6. Testing methodology and verification

- **Per slice:** narrowest pytest first, all other slices mocked; tolerances
  justified (tight analytical for E1/E3/E4; sample-size-based Monte-Carlo for E2's
  stochastic cascade; invariant-closure for E2; real-PCG64 stream-parity for the
  E2 draw-consumption lock). Interpreter:
  `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q`.
- **Fixture reuse:** hand-built v7 checkpoints via the
  `test_checkpoint.py::_make_ion_checkpoint` pattern; driver-produced tiny runs
  via `test_biphasic_step.py::_tiny_neutral` / `_biphasic_driver_cfg` (N=2,
  0.2 ps — the delivered bridge-smoke precedent).
- **Integration smoke (synthetic only, no figures):** tiny synthetic v7
  checkpoint → E2 relaxation → E1 size distribution → E4 Wasserstein vs the real
  abundance CSV → E5 diagnostics; assert the 5-term invariant closes on both
  windows and the metric is finite.
- **Suite:** confirm the baseline (1840 green at HEAD) stays green; any new red
  is a real regression. The E5a rename must keep the Phase-D report script and
  its test suite green under the new module/test names.

---

## 7. Risks / notes

- **R5 upper-bound (carried).** Even with the relaxation stage, the cascade may
  not fully complete at $t_\text{exp}$; E4/F3 always report both the matched-time
  and the sim-end values and never over-read absolute `n`.
- **Relaxation-stage cost (downgraded to LOW by R1/audit).** The mass subsystem
  decouples from translation, so the `free_flight` arm is exact for the
  observable and nearly free per step; even the default `coulomb` arm is O(N)
  FD-potential evals. At `dt = 0.01` ps, 500 ps ≈ 50 k steps on a tiny state —
  bounded further by the freeze early-exit and the auto-stride storage budget.
  No larger `dt` is authorized beyond the $\nu\,dt \le 0.1$ guard.
- **Pinned-point expectation (bridge finding, recorded).** At the 9 Å / 0.80 eV
  representative point the ensemble freezes by ~8 ps ($E_\text{int}(30\,\text{ps})
  \approx 0.004$ eV $< D_0$), so the relaxation stage adds ≈ 0 sheds there —
  matched-time ≈ sim-end. The stage earns its keep at hotter points (2.70 eV
  production, larger κ) where the cascade is still live at 20 ps. E2's dynamic
  tests therefore use synthetic hot inputs.
- **Calibration load (HIGH, carried into R2).** One observable (size
  distribution) carries 8+ quantities; whether it *separates* them is the open
  empirical question the Phase F campaign (R2) must report, not assume.
- **Velocity data asymmetry (recorded).** The `vmi_summary` references are
  aggregate, not per-`n`; this is *why* Phase E commits to the size distribution
  only. If per-fragment velocity references are later exported, a
  velocity-comparison slice can be reopened.
- **Provenance gap (open item, R3).** `integrated_i_he_abundance.csv` ships
  without a documented producer; the E3 README entry records the verified data
  contract and flags provenance "to be completed" (user to supply measurement
  IDs / source script).

**Cross-links:** `TIER2_IMPLEMENTATION_PLAN.md` §4 (Phase E headline);
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` §2.1–2.2 / §4 / §6.11 / §R5 /
§7 / §11; `TIER2_PHASE_D_BRIDGE_FINDINGS.md` (pinned-point oracles);
`TIER2_PHASE_F_IMPLEMENTATION_PLAN.md` (the consumer of every Phase-E surface);
`drag_migration_log_tier2.md` (the Phase E decision + delivery record).
