# Tier 2 — G4 Step 2 (h405 battery + W₁ anatomy) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:executing-plans`
> (inline) or `superpowers:subagent-driven-development` to implement this plan
> task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and execute atlas plan §3.5f — the pooled 5 × N = 1000
verification battery at h405 (Block V), the signed per-bin CDF-gap
decomposition of the W₁ residual (Block D), and the pre-registered mechanism
fingerprints (Block F), so the user's W₁-floor call (chase a mechanism knob
vs record the honest residual) is made on a measured residual anatomy.

**Architecture:** Three small additions, all reusing committed machinery.
The battery generator reuses `gen_tier2atlas_g4finals.build_cell` on the
frozen h405 cell with only the seed swapped (zero pin duplication). The W₁
decomposition is a new shared primitive `cdf_gap_profile` in
`postprocess/distribution_compare.py` that the committed
`wasserstein_integer_support` delegates to (rule 1 — one CDF-alignment
implementation, bit-exact delegation guarded by the existing test suite).
Pooling is a new pure function `pool_confirmation_reads` in
`postprocess/tier2_confirmation.py` (no pooled container dir is built —
the standing container was a hand-assembled figures container, not a
pattern to repeat for scoring). Two new pure-scorer reports sit in
`scripts/post_processing/`.

**Tech Stack:** Python 3.14, numpy, pytest. Interpreter:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe'`.

## Global Constraints

- **Physics is frozen for this stage.** No change to drag form, b, the
  committed scorer conventions (`postprocess/tier2_confirmation` numeric
  behavior), checkpoint schema, RNG draw order, constants, or `SimConfig`
  defaults. The only MD is Block V: 5 × N = 1000 at h405.
- **h405 cell = the frozen finals cell:** v_c 5.5 / τ 4.4 / E₀ 0.405,
  eb1168 bundle well, corrected geometry (`legacy`+`raw`, Boltzmann
  313.2 K, margin 0), `rq4graded` ladder, p_tail −1, co-moving shed,
  Landau-gated drag 0.58 (58 m/s), budget 2.70, `exclude_all_coupled`.
  All of it comes from `gen_tier2atlas_g4finals._SPEC_BY_LABEL["h405"]` +
  its `build_cell` — this plan adds **no new pin values**, only seeds.
- **Battery seeds: 20260730, 20260731, 20260732, 20260733, 20260734**
  (five fresh; verified unused by any committed generator before launch).
- **Oracles run before any new number is read** (§1.4): the h405 frozen
  twin row string-exact; each member cfg diffs against the committed
  `g4fh405` run cfg in **exactly** `{"seed"}`; the pooled-battery
  scorer-drift oracle; the committed h405/h410/h415 finals rows rescored
  equal to `atlas_g4finals_table.csv`; and the decomposition-identity
  oracle Σ|gap| == `w1_solvated` per distribution.
- **Block-F-before-Block-D ordering is hard:** the fingerprints section
  must be committed to the findings doc before the anatomy scorer's
  output is read. Building the scorer earlier is fine; running it and
  reading numbers is not.
- **Nothing adopts.** `finc1v725` stands until the user's G4
  adjudications. GV-P1..P3 make them *fireable*, not fired.
- **Do not read generator launch prints as n̄** (the `n_detect_mean` trap
  — launch print lacks the §4r suppressed→bin-0 convention).

## File Structure

- Modify: `i2_helium_md/postprocess/distribution_compare.py` — add
  `CdfGapProfile` + `cdf_gap_profile`; `wasserstein_integer_support`
  delegates (bit-exact).
- Modify: `i2_helium_md/postprocess/tier2_confirmation.py` — add
  `pool_confirmation_reads`.
- Modify: `scripts/post_processing/tier2atlas_geometry_table.py` — extract
  `observable_columns(read, ...)` from `observable_row` (pure refactor;
  scorer-drift oracle guards it).
- Create: `scripts/gen_tier2atlas_g4step2_battery.py` — Block V generator.
- Create: `scripts/post_processing/tier2atlas_g4step2_battery_table.py` —
  Block V scorer (per-member + pooled rows, GV-P1..P3).
- Create: `scripts/post_processing/tier2atlas_g4step2_w1_anatomy.py` —
  Block D scorer (signed gap profiles, ≥70 % bin set, E₀-lever signature).
- Tests: `tests/test_distribution_compare.py` (extend),
  `tests/test_tier2_confirmation.py` (extend),
  `tests/test_tier2atlas_g4step2_anatomy.py` (new).
- Docs: findings "G4 Step 2" sections (Block F first), migration-log
  entry, D0 refinements (E₀-lever bin signature under the E₀ knob).
- Committed artifacts: `data/runs/h2b_forward_model/atlas_g4step2_battery.csv`,
  `data/runs/h2b_forward_model/atlas_g4step2_w1_anatomy.csv`.

---

### Task 1: `cdf_gap_profile` — the shared W₁ decomposition primitive

**Files:**
- Modify: `i2_helium_md/postprocess/distribution_compare.py`
- Test: `tests/test_distribution_compare.py`

**Interfaces:**
- Produces: `cdf_gap_profile(sim, ref) -> CdfGapProfile` with fields
  `support: np.ndarray[int]`, `f_sim`, `f_ref`, `gaps` (signed
  `F_sim − F_ref` per support point) and property `w1: float`.
  `wasserstein_integer_support(sim, ref)` behavior unchanged (bit-exact).
- Consumes: the module's existing `_validated_distribution`.

- [ ] **Step 1: Write the failing tests** (append to
  `tests/test_distribution_compare.py`; reuse its `_sim`/`_ref` helpers):

```python
class TestCdfGapProfile:
    def test_gap_sum_equals_w1_bit_exact(self):
        sim = _sim([0, 1, 2, 3], [0.1, 0.4, 0.3, 0.2])
        ref = _ref([1, 2, 4], [0.5, 0.3, 0.2])
        prof = cdf_gap_profile(sim, ref)
        assert float(np.sum(np.abs(prof.gaps))) == \
            wasserstein_integer_support(sim, ref)
        assert prof.w1 == wasserstein_integer_support(sim, ref)

    def test_signed_gaps_hand_case(self):
        # sim all at 1, ref all at 2: F_sim - F_ref = +1 at n = 1, 0 at 2.
        prof = cdf_gap_profile(_sim([1], [1.0]), _ref([1, 2], [0.0, 1.0]))
        assert list(prof.support) == [1, 2]
        assert prof.gaps[0] == pytest.approx(1.0)
        assert prof.gaps[1] == pytest.approx(0.0)

    def test_union_support_is_contiguous(self):
        prof = cdf_gap_profile(_sim([0, 1], [0.5, 0.5]),
                               _ref([3, 4], [0.5, 0.5]))
        # (disjoint support raises before this — use overlapping instead)

    def test_validation_still_fires_through_profile(self):
        with pytest.raises(ValueError):
            cdf_gap_profile(_sim([0, 1], [0.6, 0.6]), _ref([1], [1.0]))
```

  (Replace `test_union_support_is_contiguous` body with an overlapping
  case, e.g. sim `[0, 1]`, ref `[1, 4]`, and assert
  `list(prof.support) == [0, 1, 2, 3, 4]`.)

- [ ] **Step 2: Run tests to verify they fail**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_distribution_compare.py -q`
Expected: FAIL with `ImportError: cannot import name 'cdf_gap_profile'`.

- [ ] **Step 3: Implement** — in `distribution_compare.py`, directly above
  `wasserstein_integer_support`:

```python
@dataclass(frozen=True)
class CdfGapProfile:
    """Signed per-bin CDF-gap profile of sim vs ref on the union support.

    Attributes
    ----------
    support : np.ndarray, shape (K,), int
        Contiguous union support (unit rung spacing).
    f_sim, f_ref : np.ndarray, shape (K,)
        Zero-filled fraction vectors on ``support``.
    gaps : np.ndarray, shape (K,)
        ``F_sim(n) - F_ref(n)`` per support point (signed).
    """

    support: np.ndarray
    f_sim: np.ndarray
    f_ref: np.ndarray
    gaps: np.ndarray

    @property
    def w1(self) -> float:
        """``sum |gaps|`` — exactly the committed W1 (unit spacing)."""
        return float(np.sum(np.abs(self.gaps)))


def cdf_gap_profile(sim: ShellDistribution,
                    ref: HeAbundanceReference) -> CdfGapProfile:
    """Signed CDF-gap decomposition of the committed W1 (same validation,
    same union-support alignment; ``profile.w1`` is bit-exact W1)."""
    sim_n, sim_f = _validated_distribution(sim.n_values, sim.fraction, side="sim")
    ref_n, ref_f = _validated_distribution(ref.n, ref.ion_fraction, side="ref")
    if not np.intersect1d(sim_n, ref_n).size:
        raise ValueError(
            f"disjoint integer support: sim n in [{sim_n.min()}, {sim_n.max()}], "
            f"ref n in [{ref_n.min()}, {ref_n.max()}] share no shell count -- "
            f"likely a units/labelling bug."
        )
    n_lo = int(min(sim_n.min(), ref_n.min()))
    n_hi = int(max(sim_n.max(), ref_n.max()))
    support = np.arange(n_lo, n_hi + 1)
    f_sim = np.zeros(support.shape, dtype=float)
    f_ref = np.zeros(support.shape, dtype=float)
    f_sim[sim_n - n_lo] = sim_f
    f_ref[ref_n - n_lo] = ref_f
    return CdfGapProfile(
        support=support, f_sim=f_sim, f_ref=f_ref,
        gaps=np.cumsum(f_sim) - np.cumsum(f_ref),
    )
```

  Then rewrite the body of `wasserstein_integer_support` (docstring
  unchanged) as:

```python
    return cdf_gap_profile(sim, ref).w1
```

- [ ] **Step 4: Run the full existing distribution-compare suite**

Run: `& '...python.exe' -m pytest tests/test_distribution_compare.py -q`
Expected: ALL PASS (the pre-existing W1 hand cases are the bit-exactness
guard for the delegation).

- [ ] **Step 5: Commit**

```bash
git add i2_helium_md/postprocess/distribution_compare.py tests/test_distribution_compare.py
git commit -m "tier2 atlas G4 Step 2 Task 1: cdf_gap_profile primitive; committed W1 delegates bit-exact"
```

---

### Task 2: `pool_confirmation_reads` + `observable_columns` extraction

**Files:**
- Modify: `i2_helium_md/postprocess/tier2_confirmation.py`
- Modify: `scripts/post_processing/tier2atlas_geometry_table.py:153-191`
- Test: `tests/test_tier2_confirmation.py`

**Interfaces:**
- Produces: `pool_confirmation_reads(reads, *, label) ->
  ConfirmationDetectionRead` (pure; exact integer recount of the state
  classes via `int(round(frac * denominator))`).
- Produces: `observable_columns(read, abundance_ref, ked_ref) -> dict`
  in `tier2atlas_geometry_table.py`; `observable_row` becomes
  `load_confirmation_run` + `geometry_columns` (optional) +
  `observable_columns`, output dict unchanged key-for-key.

- [ ] **Step 1: Write the failing tests** (append to
  `tests/test_tier2_confirmation.py`, following its existing synthetic
  `ConfirmationDetectionRead` construction style — build two small reads
  directly or via `read_confirmation_detection` on synthetic
  `DetectionResult`s, whichever the file already does):

```python
class TestPoolConfirmationReads:
    def _read(self, n_scored, ke, *, num_ions, trapped, bound, marginal,
              suppressed, n_max=25):
        counts = np.bincount(n_scored, minlength=n_max + 1).astype(float)
        return ConfirmationDetectionRead(
            label="m", num_ions=num_ions, num_scored=len(n_scored),
            trapped_frac=trapped / num_ions,
            trap_bound_frac=bound / num_ions,
            trap_marginal_frac=marginal / num_ions,
            suppressed_frac=suppressed / len(n_scored),
            n_scored=np.asarray(n_scored), ke_scored_eV=np.asarray(ke),
            n_values=np.arange(n_max + 1), fraction=counts / counts.sum(),
            n_mean=float(np.mean(n_scored)),
            n1_frac=float(counts[1] / counts.sum()),
        )

    def test_pool_two_members_equals_concatenation(self):
        a = self._read([1, 2, 2], [0.1, 0.2, 0.3], num_ions=4,
                       trapped=1, bound=1, marginal=0, suppressed=0)
        b = self._read([0, 4], [0.4, 0.5], num_ions=3,
                       trapped=1, bound=0, marginal=1, suppressed=1)
        pooled = pool_confirmation_reads([a, b], label="pool")
        assert pooled.num_ions == 7 and pooled.num_scored == 5
        assert pooled.trapped_frac == pytest.approx(2 / 7)
        assert pooled.trap_marginal_frac == pytest.approx(1 / 7)
        assert pooled.suppressed_frac == pytest.approx(1 / 5)
        assert pooled.fraction[2] == pytest.approx(2 / 5)
        assert pooled.n_mean == pytest.approx(np.mean([1, 2, 2, 0, 4]))

    def test_pool_self_preserves_fractions(self):
        a = self._read([1, 1, 3], [0.1, 0.1, 0.2], num_ions=3,
                       trapped=0, bound=0, marginal=0, suppressed=0)
        pooled = pool_confirmation_reads([a, a], label="pool")
        assert np.allclose(pooled.fraction, a.fraction)

    def test_empty_and_mismatched_support_raise(self):
        with pytest.raises(ValueError):
            pool_confirmation_reads([], label="pool")
        a = self._read([1], [0.1], num_ions=1, trapped=0, bound=0,
                       marginal=0, suppressed=0, n_max=25)
        b = self._read([1], [0.1], num_ions=1, trapped=0, bound=0,
                       marginal=0, suppressed=0, n_max=30)
        with pytest.raises(ValueError):
            pool_confirmation_reads([a, b], label="pool")
```

- [ ] **Step 2: Run to verify failure** (ImportError on
  `pool_confirmation_reads`).

- [ ] **Step 3: Implement** in `tier2_confirmation.py` (below
  `load_confirmation_run`; add to `__all__`):

```python
def pool_confirmation_reads(
    reads: Sequence[ConfirmationDetectionRead], *, label: str
) -> ConfirmationDetectionRead:
    """Pool member reads into one read (the §4cc battery convention).

    Equivalent to reducing the concatenated detection arrays: scored ions
    concatenate; the retained/suppressed classes recount exactly via
    ``int(round(frac * denominator))`` (each frac was ``count / denom``).
    No pooled container dir is involved.

    Raises
    ------
    ValueError
        On an empty member list or mismatched histogram supports.
    """
    if not reads:
        raise ValueError("pool_confirmation_reads: no member reads.")
    base = reads[0].n_values
    for r in reads[1:]:
        if r.n_values.shape != base.shape or not np.array_equal(r.n_values, base):
            raise ValueError(
                f"pool {label!r}: mismatched histogram support "
                f"({r.label!r} vs {reads[0].label!r})."
            )
    num_ions = sum(r.num_ions for r in reads)
    n_scored = np.concatenate([r.n_scored for r in reads])
    ke = np.concatenate([r.ke_scored_eV for r in reads])
    trapped = sum(int(round(r.trapped_frac * r.num_ions)) for r in reads)
    bound = sum(int(round(r.trap_bound_frac * r.num_ions)) for r in reads)
    marginal = sum(int(round(r.trap_marginal_frac * r.num_ions)) for r in reads)
    suppressed = sum(int(round(r.suppressed_frac * r.num_scored)) for r in reads)
    counts = np.bincount(n_scored, minlength=base.size).astype(float)
    fraction = counts / counts.sum()
    return ConfirmationDetectionRead(
        label=label, num_ions=num_ions, num_scored=int(n_scored.size),
        trapped_frac=trapped / num_ions,
        trap_bound_frac=bound / num_ions,
        trap_marginal_frac=marginal / num_ions,
        suppressed_frac=suppressed / n_scored.size,
        n_scored=n_scored, ke_scored_eV=ke,
        n_values=base.copy(), fraction=fraction,
        n_mean=float(n_scored.mean()), n1_frac=float(fraction[1]),
    )
```

  Then in `tier2atlas_geometry_table.py` split `observable_row`: move the
  body from `hist = score_histogram_vs_reference(...)` through the
  `row.update({...})` dict into

```python
def observable_columns(read, abundance_ref, ked_ref) -> dict[str, Any]:
    """The committed observable vector for one reduced read (pure)."""
```

  returning that dict (keys `num_scored` … `ke_npts` unchanged), and have
  `observable_row` call `load_confirmation_run` + optional
  `geometry_columns` + `observable_columns`. No key or numeric change.

- [ ] **Step 4: Run the tests**

Run: `& '...python.exe' -m pytest tests/test_tier2_confirmation.py tests/test_distribution_compare.py -q`
Expected: ALL PASS.

- [ ] **Step 5: Refactor guard — rescore the standing pooled battery**

Run: `& '...python.exe' scripts/post_processing/tier2atlas_geometry_table.py`
Expected: "scorer-drift oracle" section passes (trap 0.067 / supp 0.187 /
nbar 4.068 / n1 0.243 / w1 0.571 / midHot 1.0110 / deepKE 0.631 within
tol 0.002) — proves `observable_columns` extraction changed nothing.

- [ ] **Step 6: Commit**

```bash
git add i2_helium_md/postprocess/tier2_confirmation.py scripts/post_processing/tier2atlas_geometry_table.py tests/test_tier2_confirmation.py
git commit -m "tier2 atlas G4 Step 2 Task 2: pool_confirmation_reads + observable_columns extraction (oracle-guarded)"
```

---

### Task 3: Block V generator — 5 fresh-seed h405 members

**Files:**
- Create: `scripts/gen_tier2atlas_g4step2_battery.py`

**Interfaces:**
- Consumes: `gen_tier2atlas_g4finals._SPEC_BY_LABEL`, `.build_cell`,
  `.verify_corrected_geometry`, `.finals_run_dir_name`, `.TWIN_ROWS`,
  `.TWIN_SCAN_CSV`; `scripts.tier0_common.run_dir_name`;
  `scripts.tier2_common.cfg_diff_vs_reference`.
- Produces: run dirs
  `9A_drag_shared_pure_cubic_N1000_tier2atlas_conf270_g4s2h405s{1..5}`
  (exact names via `run_dir_name(CASE, VARIANT, 1000,
  run_tag="tier2atlas_conf270_g4s2h405s{i}")`), each with
  `cfg.json/neutral.npz/ion.npz/relaxation.npz/detection.npz`.
  Member labels `s1..s5` are the scorer's contract.

- [ ] **Step 1: Verify the seeds are fresh**

Run: `grep -rn "2026073" scripts/*.py scripts/post_processing/*.py`
Expected: no generator uses 20260730–20260734 (committed generators end
at 20260729). If any hit, shift the block to the next free range and
record the change in the plan file.

- [ ] **Step 2: Write the generator.** Same skeleton as
  `gen_tier2atlas_g4finals.py` (thread pinning, UTF-8, argparse with
  `--dry-run`/`--concurrency`/member labels, `SKIP_COMPLETED_RUNS = True`,
  process pool). The heart — zero pin duplication:

```python
from scripts.gen_tier2atlas_g4finals import (
    TWIN_ROWS, TWIN_SCAN_CSV, _SPEC_BY_LABEL,
    build_cell, finals_run_dir_name, verify_corrected_geometry,
)

H405 = _SPEC_BY_LABEL["h405"]
MEMBER_SEEDS = {"s1": 20260730, "s2": 20260731, "s3": 20260732,
                "s4": 20260733, "s5": 20260734}
COMMITTED_H405_RUN = finals_run_dir_name("h405")   # the Block-3 run

def build_member(member: str) -> SimConfig:
    """h405's frozen cell, seed swapped — nothing else may move."""
    cfg = dataclasses.replace(build_cell(H405), seed=MEMBER_SEEDS[member])
    cfg.validate()
    verify_corrected_geometry(cfg, H405)
    return cfg

def verify_member_vs_committed_h405(cfg: SimConfig, member: str) -> None:
    """Oracle: the member cfg diffs vs the committed g4fh405 cfg.json in
    exactly {'seed'} — same cell, fresh seed, no drifted pin."""
    ref = PROJECT_ROOT / "data" / "runs" / COMMITTED_H405_RUN / "cfg.json"
    diff = cfg_diff_vs_reference(cfg, ref, context=f"g4s2 {member}")
    if set(diff) != {"seed"}:
        raise AssertionError(
            f"[{member}] cfg vs committed h405 differs in {sorted(diff)}; "
            f"pre-registered exactly ['seed']."
        )

def member_run_dir_name(member: str) -> str:
    name = run_dir_name(CASE, VARIANT, N,
                        run_tag=f"tier2atlas_conf270_g4s2h405{member}")
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"atlas namespace violated by {name!r}")
    return name
```

  Launch oracle (before any MD, and under `--dry-run`): re-run the h405
  slice of the finals twin oracle — read `TWIN_SCAN_CSV`, assert the row
  at h405's key equals `TWIN_ROWS["h405"]` string-exact (same code shape
  as `verify_twin_preregistration`, restricted to h405). `_run_one`
  mirrors the finals generator (neutral → ion → relaxation → detection
  into a `RunDirectory`), printing the standard stage lines. Docstring
  carries the §3.5f Block V registration (GV-P1..P3 text) and the
  launch-print n̄ warning.

- [ ] **Step 3: Dry-run**

Run: `& '...python.exe' scripts/gen_tier2atlas_g4step2_battery.py --dry-run`
Expected: twin-row oracle PASSED; five members print
`cfg OK ... diff vs committed h405: ['seed']`; namespace guard silent.

- [ ] **Step 4: Commit the generator, then launch in the background**

```bash
git add scripts/gen_tier2atlas_g4step2_battery.py
git commit -m "tier2 atlas G4 Step 2 Task 3: Block V battery generator (h405 x 5 fresh seeds, seed-only cfg oracle)"
```

Run (background, ~60–90 min at concurrency 3):
`& '...python.exe' scripts/gen_tier2atlas_g4step2_battery.py --concurrency 3`
Proceed to Task 4 while it runs.

---

### Task 4: Block F — the fingerprint pre-registration (doc work, while the battery runs)

**Files:**
- Modify: `docs/drag_port/Tier2/TIER2_SENSITIVITY_ATLAS_FINDINGS.md` —
  new section `## G4 Step 2 — Block F: mechanism fingerprints (FROZEN
  before any Block-D read) (2026-07-28)`.

**Interfaces:**
- Consumes (read-only): `MASS_DYNAMICS_LOCKED_energy_gated_evaporation.md`
  (mechanism equations), `docs/drag_port/Tier2/RESEARCH_QUESTIONS.md`
  (RQ2 per-shed ε, RQ4 ladder), `TIER2_PARAMETER_INFLUENCE.md` (c1 /
  λ0 / ε knob entries + §17 ledger), the biphasic mechanism source
  (`i2_helium_md/physics/` evaporation/pickup modules) as needed.
- Produces: the frozen fingerprint table the Block-D read is judged
  against (plan §3.5f match rule).

- [ ] **Step 1: Write the section.** One subsection per knob — **ladder
  shape (c1 grading)**, **pickup (λ0)**, **per-shed ε (RQ2)** — each
  containing exactly these fields:
  1. *Mechanism route* — 2–4 sentences deriving how the knob enters the
     cascade (equation references, not vibes).
  2. *Sign on n₁_solv* (+/−/0, with the one-line reason).
  3. *Sign on n̄_det* (+/−/0).
  4. *Bin-pattern* — which n-bins of the detected size distribution the
     knob reshapes, stated as a signed CDF-gap direction per bin range
     (e.g. "raises F in n = 2–4, lowers it in n ≥ 6").
  5. *The §3.5f constraint answer* — can this knob add solvated-shoulder
     mass **without** paying it out of n₁? (yes/no/conditional, derived.)
  Close the section with the frozen match rule and decision table copied
  verbatim from plan §3.5f (so the findings section is self-contained)
  and the sentence: "Frozen before any Block-D number was read."

- [ ] **Step 2: Commit** (this commit timestamp is the freeze evidence —
  it must predate the anatomy scorer's first run):

```bash
git add docs/drag_port/Tier2/TIER2_SENSITIVITY_ATLAS_FINDINGS.md
git commit -m "tier2 atlas G4 Step 2 Block F: mechanism fingerprints FROZEN (pre-Block-D)"
```

---

### Task 5: Block V scorer — battery table + GV-P1..P3

**Files:**
- Create: `scripts/post_processing/tier2atlas_g4step2_battery_table.py`
- Test: none beyond reuse (pure scorer built entirely from tested parts;
  the finals-table precedent).

**Interfaces:**
- Consumes: `load_confirmation_run`, `pool_confirmation_reads` (Task 2),
  `observable_columns` / `check_oracle` / `ORACLE_RUN` / reference CSVs
  from `tier2atlas_geometry_table.py`, `joint_score` / `PROVISIONAL_NORM`
  from `tier2atlas_g4_transfer.py`, member dir names from Task 3's
  `member_run_dir_name`, the committed `atlas_g4finals_table.csv`.
- Produces: `data/runs/h2b_forward_model/atlas_g4step2_battery.csv` —
  rows s1..s5 + `pooled`; the GV verdict block.

- [ ] **Step 1: Write the scorer.** Order of operations (§1.4 — oracles
  first, hard-fail before any new number prints):
  1. Scorer-drift oracle on the standing pooled battery (`observable_row`
     on `ORACLE_RUN` + `check_oracle`) — abort on drift.
  2. Committed-row oracle: rescore the committed `g4fh405` run via
     `load_confirmation_run` + `observable_columns`; assert nbar_det /
     n1_solv / w1_solv / midHot / deepKE equal the committed
     `atlas_g4finals_table.csv` h405 row to 4 decimals — abort on drift.
  3. Load s1..s5 (`load_confirmation_run` each), print per-member rows
     (`observable_columns` + `md_gate` + S on the licensed axes
     `("w1", "midhot")`, exactly the finals `_s` convention).
  4. `pooled = pool_confirmation_reads(members, label="h405pooled")` →
     the pooled row the same way.
  5. Seed statistics: per-member W₁ mean ± SD next to the Block-0
     normalizers (0.5791 ± 0.0954 at the standing point).
  6. GV verdicts, each printed CONFIRMED/REFUTED with its number:
     - **GV-P1:** pooled n1_solv ∈ [0.19, 0.30] ∧ pooled nbar_det ∈
       [3.77, 4.37].
     - **GV-P2:** pooled w1_solv ∈ [0.64, 0.78]; per-seed SD printed
       next to 0.095 (informational, not a hard clause).
     - **GV-P3:** pooled S < 1.683 (a037's N = 1000 S); distance to the
       Block-3 single-seed 1.559 printed.
     - Closing line: "GV pass makes the G4 adjudications FIREABLE
       (user); nothing adopted here."
  7. `write_rows_csv` to the artifact path.

- [ ] **Step 2: Run it** (battery must be finished; if members are
  missing it must print them and stop, the finals-table pattern).

Run: `& '...python.exe' scripts/post_processing/tier2atlas_g4step2_battery_table.py`
Expected: both oracles pass, six rows print, GV verdicts print.

- [ ] **Step 3: Commit scorer + artifact**

```bash
git add scripts/post_processing/tier2atlas_g4step2_battery_table.py data/runs/h2b_forward_model/atlas_g4step2_battery.csv
git commit -m "tier2 atlas G4 Step 2 Block V: battery table + GV verdicts"
```

---

### Task 6: Block D scorer — the W₁ anatomy

**Files:**
- Create: `scripts/post_processing/tier2atlas_g4step2_w1_anatomy.py`
- Test: `tests/test_tier2atlas_g4step2_anatomy.py`

**Interfaces:**
- Consumes: `cdf_gap_profile` (Task 1), `pool_confirmation_reads`
  (Task 2), `solvated_renormalized` / `load_confirmation_run` /
  `score_histogram_vs_reference` from `tier2_confirmation`,
  `load_he_abundance_reference`, member dirs (Task 3), the committed
  finals run dirs `g4fh405`/`g4fh410`/`g4fh415`, the standing pooled
  container `ORACLE_RUN`, `atlas_g4finals_table.csv`.
- Produces: `data/runs/h2b_forward_model/atlas_g4step2_w1_anatomy.csv`
  (one row per (distribution, n): signed gap, |gap|, share, cum-share,
  in-top-70 flag, f_sim, f_ref) + two helpers the tests cover:

```python
def solvated_gap_profile(read, abundance_ref) -> CdfGapProfile:
    """Gap profile on the RQ8 solvated branch (n >= 1 renormalized) —
    EXACTLY the branch w1_solvated is scored on."""
    sim_n, sim_f = solvated_renormalized(read.n_values, read.fraction)
    ref_n, ref_f = solvated_renormalized(abundance_ref.n,
                                         abundance_ref.ion_fraction)
    return cdf_gap_profile(
        SimpleNamespace(n_values=sim_n, fraction=sim_f),
        SimpleNamespace(n=ref_n, ion_fraction=ref_f),
    )

def top_share_bins(profile: CdfGapProfile, share: float = 0.70
                   ) -> tuple[np.ndarray, float]:
    """Minimal bin set carrying >= share of W1: sort by (-|gap|, n),
    take until the cumulative |gap| fraction first reaches share.
    Returns (bins ascending, exact share reached). share in (0, 1]."""
```

- [ ] **Step 1: Write the failing tests** (synthetic profiles only — no
  run dirs, no figures):

```python
class TestTopShareBins:
    def _profile(self, gaps):
        g = np.asarray(gaps, dtype=float)
        z = np.zeros_like(g)
        return CdfGapProfile(support=np.arange(1, g.size + 1),
                             f_sim=z, f_ref=z, gaps=g)

    def test_minimal_set_and_share(self):
        prof = self._profile([0.5, -0.3, 0.1, -0.1])   # W1 = 1.0
        bins, reached = top_share_bins(prof, 0.70)
        assert list(bins) == [1, 2] and reached == pytest.approx(0.8)

    def test_tie_broken_by_n(self):
        prof = self._profile([0.4, -0.4, 0.2])
        bins, _ = top_share_bins(prof, 0.5)
        assert list(bins) == [1, 2]   # equal |gap|: lower n first

    def test_share_one_returns_all_nonzero(self):
        prof = self._profile([0.4, 0.0, -0.6])
        bins, reached = top_share_bins(prof, 1.0)
        assert list(bins) == [1, 3] and reached == pytest.approx(1.0)

    def test_bad_share_raises(self):
        with pytest.raises(ValueError):
            top_share_bins(self._profile([1.0]), 0.0)

class TestSolvatedGapIdentity:
    def test_matches_committed_w1_on_synthetic(self):
        # synthetic read + reference; the identity oracle in unit form
        ...  # build tiny ConfirmationDetectionRead + HeAbundanceReference,
             # assert solvated_gap_profile(read, ref).w1 ==
             #   pytest.approx(score_histogram_vs_reference(read, ref)
             #                 .w1_solvated, abs=1e-12)
```

  (Fill `TestSolvatedGapIdentity` with the same synthetic-read
  constructor as Task 2's tests plus a two-bin reference; the assertion
  shown is the whole test.)

- [ ] **Step 2: Run to verify failure**, then **Step 3: implement** the
  two helpers in the scorer module (`top_share_bins` raising
  `ValueError` unless `0 < share <= 1`) and the report body:
  1. Oracles first: scorer-drift (standing pooled battery); identity
     oracle `solvated_gap_profile(read).w1 == w1_solvated` to 1e-12 for
     **every** distribution scored in the report; committed-row oracle:
     rescored W₁ of `g4fh405`/`g4fh410`/`g4fh415` equal
     `atlas_g4finals_table.csv` to 4 decimals.
  2. Profiles + tables printed for: **pooled h405** (Task 3 members),
     **standing pooled finc1v725** (`ORACLE_RUN` container), and the
     three arm cells. Per profile: signed gap per bin, share, cum share,
     `top_share_bins(profile, 0.70)`.
  3. **E₀-lever signature:** per-bin Δgap between (h405→h410) and
     (h410→h415) scaled to +0.005 eV — the measured knob-signature
     template the Block-F fingerprints are compared against.
  4. Closing print: "Block F is frozen in the findings doc; apply the
     §3.5f match rule there. This scorer reports data, not verdicts."
  5. CSV artifact.

- [ ] **Step 4: Run the tests + full twin/postprocess suites**

Run: `& '...python.exe' -m pytest tests/test_tier2atlas_g4step2_anatomy.py tests/test_distribution_compare.py tests/test_tier2_confirmation.py -q`
Expected: ALL PASS.

- [ ] **Step 5: Commit the scorer + tests** (running it on real data is
  Task 7 — do NOT run the report before the Block-F commit exists):

```bash
git add scripts/post_processing/tier2atlas_g4step2_w1_anatomy.py tests/test_tier2atlas_g4step2_anatomy.py
git commit -m "tier2 atlas G4 Step 2 Block D: W1 anatomy scorer (identity-oracle-gated); not yet run"
```

---

### Task 7: Execute + record

- [ ] **Step 1: Preconditions** — battery finished (5 complete run
  dirs); Task 4's Block-F commit exists (`git log --oneline -- docs/...`
  shows it). Both must hold before any report runs.

- [ ] **Step 2: Run the Block V table** (Task 5 Step 2 if not yet run) —
  record GV-P1..P3.

- [ ] **Step 3: Run the anatomy report**

Run: `& '...python.exe' scripts/post_processing/tier2atlas_g4step2_w1_anatomy.py`
All oracles must pass; then read the profiles **against the frozen
Block-F fingerprints** and write the match verdict per the §3.5f rule
(live candidate(s) / none).

- [ ] **Step 4: Docs, in this order** (memory rule: D0 first):
  1. `TIER2_PARAMETER_INFLUENCE.md` — the E₀ knob entry gains the
     measured bin-level signature; add/refresh mechanism-knob GAP
     markers if the match verdict names a candidate.
  2. Findings: `## G4 Step 2 — Block V` (battery + GV) and
     `## G4 Step 2 — Block D` (anatomy + fingerprint-match verdict)
     sections.
  3. Plan §3.5f: EXECUTED stamp with the one-paragraph outcome.
  4. Migration log: one EXECUTED entry (chronology + artifacts + the
     decision-table outcome).
- [ ] **Step 5: Full test suite** —
  `& '...python.exe' -m pytest -q` (report exact counts honestly).
- [ ] **Step 6: Commit** docs + CSVs on `drag_implementation`.
- [ ] **Step 7: Report to the user** — GV verdicts, the residual
  anatomy, the fingerprint match outcome, and the adjudication menu
  (successor / retained policy / ledger / the W₁-floor call), which are
  **user decisions** — stop there.

## Self-review notes

- Spec coverage: §3.5f Block V → Tasks 3+5; Block D → Tasks 1+6;
  Block F → Task 4; ordering discipline → Global Constraints + Task 6
  Step 5 + Task 7 Step 1; adjudication boundary → Task 7 Step 7.
- The delegation refactor of `wasserstein_integer_support` is the one
  committed-code numeric surface touched; it is guarded three ways
  (hand-case tests bit-exact, scorer-drift oracle, committed-row
  oracles).
- `observable_columns` extraction is the second refactor; guarded by
  Task 2 Step 5's oracle rescore.
- No figures, no production checkpoints in tests; all new tests are
  synthetic-only.
