# Tier 2 — Phase E Implementation Plan (Tier-2 Observable / Comparison Layer)

> **Boundary.** This is a *plan*, not code. The strict Physics-Definition /
> Software-Implementation boundary holds: no Python until the explicit
> `[PROCEED TO IMPLEMENTATION]` trigger. Equations are the *locked* formulations
> from `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`; the module descriptions
> are *interface contracts*, not implementations.
>
> **Parent:** `TIER2_IMPLEMENTATION_PLAN.md` §3–§4 (Phase E = slices H, W, D2).
> This doc expands that headline into an executable slice plan.
> **Entry docs:** MASS doc §2.1–2.2 (the observable), §4 (channels), §6.11
> (regime / `Π` / `t_×`), §R5 (truncation), §7 (schema), §11 (config);
> `CALIBRATION_MAP.md` (parameter classes).

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
   and figures are deferred to the Phase F campaign (R2).

The locked mechanism is **not** re-litigated: the relaxation stage is the locked
Q/K/U channels composed with pickup and drag switched off — exactly the
post-ejection regime MASS §R5 describes.

---

## 1. Hard prerequisites and test strategy

Phase E binds to these **accepted** upstream modules:

| Prereq | Provides | Needed by |
|---|---|---|
| **X** (Phase C) — checkpoint v6→v7 + 5-term invariant | per-step `E_int_eV (2N,T)`; extended `ion_ledger_closure` | E2, E5 |
| **L** (Phase A) — ladder `D_0(n)` + cumulative `Σ_{i≤n}D_0` | the self-bound gate threshold | E2, E5 |
| **K** (Phase A) — Newton cooling K2 + occupancy-resolved `E_∞(N)` | cooling drain; total-strip reachability | E2, E5 |
| **U** (Phase A) — `E_int` budget rules (S1/K1) | K1 drain on shed | E2 |
| **Q** (Phase B) — energy-gated RRK evaporation | the cascade rate engine | E2 |
| **P/ρ** (Phase B) — `λ_attach(n,ρ)` | the `Π` order parameter | E5 |
| **G** (Phase C) — generative driver | the per-step loop E2 reuses (pickup+drag off) | E2 |

**Parallelizable now (against v7-shaped synthetic checkpoints):** E1 (sim-end
mode), E3, E4 — none of these need `E_int`. **Gated on Phase C:** E2 (then E1
gains its relaxed-input mode). **Gated on L+P+K+X:** E5.

**Test strategy (mirrors Tier-1a).** Pure functions against closed-form oracles
(E1, E3, E4); stochastic units against seeded/mock RNG with moment checks (E2);
the **5-term invariant** as the cross-cutting correctness gate (E2). **No figures
or production-sized checkpoints in pytest** — tiny synthetic checkpoints only.

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
- **Reuse (exact paths).**
  - `i2_helium_md/physics/shell_schedule.py :: complex_mass_amu(n)`
    ($m(n)=126.90+4.0026\,n$ amu) — rung labelling / reference-window match.
  - `i2_helium_md/simulation/run_directory.py :: RunDirectory.load_ion/load_cfg`;
    `i2_helium_md/simulation/checkpoint.py :: IonCheckpoint.n_shell (2N,T)`,
    `mass_final_kg`.
- **New surface.** `postprocess/size_distribution.py ::
  compute_terminal_shell_distribution(source, *, n_max=21) -> ShellDistribution`
  (frozen: `n_values (Nn,)`, `counts (Nn,)`, `fraction (Nn,)`, `source`). Accepts
  **either** a raw `IonCheckpoint` (sim-end = R5 upper bound) **or** a
  `RelaxationResult` (E2, matched/relaxed time).
- **Oracle / tests.** Hand-counted tiny checkpoint (e.g. `n=[0,0,1,2,2,2]`
  → `fraction` matches by hand); integer support exact; monotone-falling envelope
  reproduced on a synthetic stream; empty/degenerate ensemble raises loudly.
- **Acceptance.** Matches hand-count; integer support; both input modes accepted;
  `fraction` sums to 1.

### E2 — Post-ejection relaxation stage (R5 mitigation)

- **Purpose.** Propagate the **mass subsystem** forward past the 20 ps ion stage
  to the experimental timescale, so the size distribution is read at *matched
  time*, not at the truncated upper bound (MASS §R5).
- **Encoded form (reduced dynamics — locked channels, no new physics).** Outside
  the droplet $\rho_\text{He}\to0$, so **pickup is off** ($\lambda_\text{attach}=0$)
  and the **He drag is off** (no bubble). The only active channels are **Q**
  (energy-gated RRK evaporation) and **K2** (Newton cooling of
  $E_\text{solv.struct}$), with $E_\text{int}$ drained by **K1** (U) on each shed.
  Translation is free-flight (residual inter-ion Coulomb retained for invariant
  fidelity; negligible at large $R$). Per relaxation step:
  1. **K2** drains $E_\text{int}$ toward $E_\infty(N)$ (Slice K);
  2. at most one **shed** via Q — the self-bound gate
     $E_\text{int}>\sum_{i\le n}D_0(i)$ suppresses; else RRK
     $k=\nu\,(1-D_0(n)/E_\text{int})^{\,s-1}$ ($n{=}1$ direct dissociation $k=\nu$);
  3. on fire: $n\to n-1$, $E_\text{int}\mathrel{-}=D_0(n)$, **cold-shed
     momentum reset** (reuse `physics/mass_jump.cold_shed`), free-flight position
     update.
  Run to `relaxation_time_ps` ($t_\text{exp}$) **or** until every fragment freezes
  ($E_\text{int}<D_0(n)\ \forall$ atoms → no further sheds; $\Pi\to0$).
- **Reuse (exact paths).** The Slice-G per-step composition and the
  `simulation/ion_propagation_step.py :: shed_step` jump-then-O seam (SQ2/SQ3) —
  **with pickup and drag disabled** (the relaxation stage is G in its free-flight,
  evaporation-only limit); `physics/mass_jump.py :: cold_shed`; Slice K cooling;
  Slice Q rate; Slice U K1 drain.
- **New surface.** `simulation/relaxation_stage.py ::
  run_relaxation_stage(ion, cfg) -> RelaxationResult` (frozen: terminal `n` per
  atom at matched time, `E_int_final (2N,)`, `time_relaxed_ps`, `freeze_flags`,
  5-term ledger arrays for the closure gate). **No checkpoint schema bump** — emit
  a separate lightweight relaxation artifact; reuse v7 `E_int`.
- **Config (new, all opt-in; default off keeps default scope unchanged).**
  `relaxation_stage_enabled=False`, `relaxation_time_ps` ($t_\text{exp}$),
  `relaxation_dt_ps` (default `dt_ion`; larger allowed only with a justified
  $k\,dt\ll1$), `relaxation_forces ∈ {coulomb(default), free_flight}`.
- **RNG.** Extends the **locked draw order** (Slice X): the relaxation stage draws
  **only** shed Bernoullis (pickup off). Documented in the draw-order spec.
- **Oracle / tests.** Pickup-off invariant — in free-flight with drag frozen the
  **5-term invariant still closes** ($E_\text{dissip}$ constant; the $E_\text{int}$
  drain closes through bath + $E_\text{pot}/E_\text{mass\_transfer}$). A fully
  self-bound input ($E_\text{int}<D_0(n)\ \forall$) relaxes with **zero sheds**
  ($n$ unchanged); a hot input sheds monotonically until freeze; no avalanche
  (per-step shed $\lesssim \nu\,dt$). Seeded RNG, tiny synthetic input.
- **Acceptance.** Reduced driver conserves the 5-term invariant; gate +
  no-avalanche hold; freeze-termination reached; matched-time terminal `n`
  distinct from (and $\le$) the sim-end upper bound on a hot synthetic case.

### E3 — Experimental abundance reference loader

- **Purpose.** Load `integrated_i_he_abundance.csv` (the Tier-2 arbiter) under the
  same validate-early contract as the HeDFT loader.
- **Reuse (mirror the contract).** `postprocess/hedft_loader.py ::
  load_hedft_trajectory` — frozen-dataclass return, order-independent header match,
  loud `FileNotFoundError`/`ValueError` on missing/extra columns, `source_path`
  retained.
- **Data contract (verified columns).** `n, label, massCenter_u_per_e,
  massWindowLower_u_per_e, massWindowUpper_u_per_e, ionCounts, ionPercent`
  (rows $n=0\dots20$; $n=0$ = bare I⁺; fraction column `ionPercent`).
- **New surface.** `postprocess/abundance_loader.py ::
  load_he_abundance_reference(path) -> HeAbundanceReference` (frozen: `n (int,Nn,)`,
  `label`, `mass_center_u (Nn,)`, `ion_counts (Nn,)`, `ion_fraction (Nn,)` =
  `ionPercent/100` normalized, `source_path`).
- **Validation.** Columns present; `n` monotone $0\dots20$; `ionPercent` sums
  $\approx100$ (tolerance-checked); raise loudly otherwise.
- **Acceptance.** Loads the real CSV; validation fires on a corrupted synthetic
  CSV; normalized `ion_fraction` sums to 1.

### E4 — Wasserstein size-distribution comparison (integer-n, pure numpy)

- **Purpose.** Score simulated (E1) vs experimental (E3) size distributions with
  the config-default Wasserstein metric on integer-n support.
- **Encoded form.** On unit-spaced integer support the 1-D $W_1$ reduces to the
  CDF-gap sum:
  $$
  W_1=\sum_n \big|\,F_\text{sim}(n)-F_\text{ref}(n)\,\big|,
  $$
  with $F$ the cumulative fraction over the **union** support (zero-filled where a
  side is absent). **Pure numpy — no scipy dependency** (a 1-D integer-support
  $W_1$ is exact and trivial; matches the package's minimal-deps style). An optional
  cross-check against `scipy.stats.wasserstein_distance` lives **in tests only**,
  guarded on scipy being importable.
- **Reuse.** `config.py :: SimConfig.validation_histogram_metric` (default
  `"wasserstein"`, currently declared-but-unread → **activated here**, removed from
  the rule-2 exception table).
- **New surface.** `postprocess/distribution_compare.py ::
  wasserstein_integer_support(sim: ShellDistribution, ref: HeAbundanceReference)
  -> float`; aligns supports; raises loudly on disjoint/empty support.
- **R5 caveat (carried in every comparison record).** Score at **matched time**
  (post-E2 relaxation) **and** report the **sim-end upper-bound** value; never
  over-read absolute `n`.
- **Oracle / tests.** Hand-computed micro-case (two 3-point distributions → known
  $W_1$); identical distributions → 0; shifting all mass by one $n$ → exactly the
  shift magnitude; disjoint support raises.
- **Acceptance.** Matches hand-computed $W_1$; metric selected via config;
  matched-time + upper-bound both reported.

### E5 — Derived diagnostics (`t_×`, `Π(t)`, regime label)

- **Purpose.** Post-hoc reconstruct, at **zero schema cost** from the v7 arrays,
  the quantities that *explain* the size distribution (MASS §6.11).
- **Encoded forms.**
  - **Crossing time $t_\times$:** first time
    $E_\text{int}(t)<\sum_{i\le n(t)}D_0(i)$ (gate opens), per atom + ensemble
    statistic. Sanity band: flag only if $t_\times<1$ ps or $>15$ ps; GAH25
    $t_0\approx5$–6.5 ps is a **±factor-2 prior, not a threshold**.
  - **Order parameter $\Pi(t)=\lambda(n)\,f_\text{ret}\,\tau$:** $\Pi>1$ shedding
    persists, $\Pi<1$ freeze, $\Pi\to0$ at exit guarantees termination.
  - **Regime label:** shell-retaining (self-binds early → moderate terminal `n`)
    ↔ total-strip (never self-binds in window → Calvo24 limit), from $t_\times$ +
    the terminal-n envelope.
  - **Total-strip reachability:** $E_\infty(N)\to0$ as $N\to0$
    (occupancy-resolved; OQ6 resolved) — confirm the Calvo24 limit is dynamically
    reachable.
- **Reuse (exact paths).** Slice L cumulative ladder $\sum_{i\le n}D_0$ (gate);
  Slice P `λ_attach(n,ρ)` (for $\Pi$); Slice K `E_∞(N)` (reachability);
  `simulation/checkpoint.py :: IonCheckpoint.E_int_eV (2N,T)`, `n_shell (2N,T)`;
  `postprocess/energy_balance.py :: ion_ledger_closure` (5-term, extended in X) for
  invariant context.
- **New surface.** `postprocess/derived_diagnostics.py ::
  reconstruct_diagnostics(ion, cfg) -> Diagnostics` (frozen: `t_cross_ps (2N,)`,
  `t_cross_summary`, `Pi_t (2N,T)`, `regime_label`, `total_strip_reachable: bool`,
  `sanity_flags`).
- **Oracle / tests.** $t_\times$ reconstructable from a synthetic
  `E_int_eV`/`n_shell` pair with a known crossing; $\Pi$ matches
  $\lambda f_\text{ret}\tau$ on synthetic inputs; regime label flips correctly
  across a constructed early-bind vs never-bind pair; reachability true when
  $E_\infty(N)\to0$.
- **Acceptance.** All four diagnostics reconstructable with zero schema cost;
  sanity flags fire on absurd inputs; regime label emitted.

---

## 3. Dependency graph and build order

```
prereqs:  X (v7) ── L ── K ── U ── Q ── P ── G        (Phases A–C, accepted)
                │
E3 (loader) ──┐ │
E4 (W1)     ──┼─┼──►  [E1 sim-end mode]   ── start now (synthetic v7 checkpoints)
E1 (size)   ──┘ │
                ▼
E2 (relaxation; needs G+Q+K+U+X) ──►  E1 relaxed / matched-time mode
                │
E5 (diagnostics; needs L+P+K+X) ───────────────────────────────────────►
```

- **Parallel now:** E3, E4, E1 (sim-end mode).
- **After Phase C accepted:** E2, then E1 gains its relaxed-input mode.
- **After L+P+K+X accepted:** E5.
- Real runs are wired in **Phase F (R2)** — Phase E is tested entirely on synthetic
  checkpoints + the real reference CSV.

---

## 4. Config contract (new / activated)

| Slice | Field | Default | Status |
|---|---|---|---|
| E2 | `relaxation_stage_enabled` (bool) | `False` | new, opt-in |
| E2 | `relaxation_time_ps` ($t_\text{exp}$) | — | new |
| E2 | `relaxation_dt_ps` | `dt_ion` | new |
| E2 | `relaxation_forces ∈ {coulomb, free_flight}` | `coulomb` | new |
| E4 | `validation_histogram_metric` (`"wasserstein"`) | `"wasserstein"` | activated (was declared-but-unread) |

Default `relaxation_stage_enabled=False` keeps the default simulation scope
unchanged (forbidden-list compliant). The relaxation draw-order extension is added
to the Slice-X locked RNG draw-order spec; new fields follow the rule-2
declared-but-unread convention until their owning slice activates them.

---

## 5. Cross-cutting gates

- **5-term invariant** is the correctness detector for E2 (pickup off, drag frozen
  → $E_\text{dissip}$ constant; the drain closes through bath +
  $E_\text{pot}/E_\text{mass\_transfer}$).
- **Integer-n support** everywhere; **Wasserstein pure-numpy** (no scipy dep).
- **R5 matched-time caveat** carried in every comparison record (E4): report both
  the relaxed/matched-time score and the sim-end upper bound.
- **Locked RNG draw order** is *extended* (not changed) for the relaxation stage.
- **Out of scope (fails review if it leaks in):** velocity comparison (cut per the
  scope decision); any noise-amplitude read (Tier 3 stays inert); changes to the
  Tier-0-locked drag law, neutral propagation, or default scope beyond the opt-in
  relaxation flag; CLI/figure scripts (deferred to R2).

---

## 6. Testing methodology and verification

- **Per slice:** narrowest pytest first, all other slices mocked; tolerances
  justified (tight analytical for E1/E3/E4; sample-size-based Monte-Carlo for E2's
  stochastic cascade; invariant-closure for E2). Interpreter:
  `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q`.
- **Integration smoke (synthetic only, no figures):** tiny synthetic v7 checkpoint
  → E2 relaxation → E1 size distribution → E4 Wasserstein vs the real abundance CSV
  → E5 diagnostics; assert the 5-term invariant closes and the metric is finite.
- **Suite:** confirm the baseline (all-green on `drag_implementation`) stays green;
  any new red is a real regression.

---

## 7. Risks / notes

- **R5 upper-bound (carried).** Even with the relaxation stage, the cascade may not
  fully complete at $t_\text{exp}$; E4 always reports both the matched-time and the
  sim-end values and never over-reads absolute `n`.
- **Relaxation-stage cost / `dt` (MEDIUM).** 20 ps → hundreds of ps is a large step
  count at `dt_ion`; `relaxation_dt_ps` may be raised only with a justified
  $k\,dt\ll1$ (the mass subsystem is slow, but the cold-shed resets must stay
  resolved). The free-flight reduction (no force/drag evaluation) keeps per-step
  cost low.
- **Calibration load (HIGH, carried into R2).** One observable (size distribution)
  carries 8+ quantities; whether it *separates* them is the open empirical
  question the Phase F campaign (R2) must report, not assume.
- **Velocity data asymmetry (recorded).** The `vmi_summary` references are
  aggregate, not per-`n`; this is *why* Phase E commits to the size distribution
  only. If per-fragment velocity references are later exported, a velocity-comparison
  slice can be reopened.

**Cross-links:** `TIER2_IMPLEMENTATION_PLAN.md` §4 (Phase E headline);
`MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md` §2.1–2.2 / §4 / §6.11 / §R5 / §7
/ §11; `drag_migration_log_tier2.md` (record the Phase E decision + delivery here).
