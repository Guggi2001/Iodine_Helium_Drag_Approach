# Tier 2 — G4 Step 1 (fine ridge sweep) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: use `superpowers:executing-plans`
> (inline) or `superpowers:subagent-driven-development` to implement this plan
> task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the machinery for atlas plan §3.5e — measure the twin's ranking
authority from the 14 paired G3-ring cells (Block 0), then scan the fine
(v_c, τ, E₀) ridge at the corrected geometry (Block 1), so the G4 successor
point is chosen on resolved axes rather than a grid coarser than the basin.

**Architecture:** Two additions, both reusing committed machinery.
Block 0 is a **pure scorer report** in `scripts/post_processing/` that joins the
committed MD ring table to the committed twin scan CSV and measures the
transfer — it runs nothing and mutates nothing. Block 1 is a **new stage** in
`scripts/tier2_h2b_forward_model.py` that reuses `_g3scan_chord_family`
(npz-cached, so the 12 existing chords cost zero) and `g3_score` unchanged,
with its own fine grids and its own gate constants.

**Tech Stack:** Python 3.14, numpy, scipy (hard dependency; `scipy.stats.spearmanr`),
pytest. Interpreter:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe'`.

## Global Constraints

- **Physics is frozen for this stage.** No change to drag form, b, the
  committed scorer (`postprocess/tier2_confirmation`), checkpoint schema, RNG
  draw order, constants, or `SimConfig` defaults. Blocks 0–1 are **zero MD**.
- **Scorer reuse is mandatory** (§1.4). Block 0 imports `observable_row` /
  `check_oracle` from `tier2atlas_geometry_table.py`; Block 1 calls the
  existing `g3_score`. No new histogram, KE-band, or W₁ code.
- **Oracles run before any new number is read** (§1.4). Block 0:
  `verify_twin_preregistration()` + the pooled-battery scorer-drift oracle.
  Block 1: the S6 machinery oracle + the committed corrected-row landmark
  oracle + G4-P1 (the 24 Step-2 gated cells reproduce bit-exact on the shared
  sub-lattice). An oracle failure aborts the stage — no partial reads.
- **Standing pins for every Block-1 cell:** corrected geometry (`legacy`+`raw`,
  Boltzmann 313.2 K, margin 0), `rq4graded` ladder, p = 1, p_tail = −1,
  co-moving shed, Landau 0.58, budget 2.70, m = 20000, `G3_SEED = 20260727`.
- **Nothing adopts.** `finc1v725` stands until the user's G4 adjudication.
- Run-directory / artifact namespace: `h2b_g4*` under
  `data/runs/h2b_forward_model/`. No literal `_tier2_` substring in any new
  run-dir tag.

## File Structure

| file | responsibility |
|---|---|
| `scripts/post_processing/tier2atlas_g4_transfer.py` (create) | Block 0. Twin↔MD transfer measurement, seed-SD normalizers, joint score S. Pure report. |
| `tests/test_tier2atlas_g4_transfer.py` (create) | Block 0 unit tests on synthetic inputs — no run dirs, no figures. |
| `scripts/tier2_h2b_forward_model.py` (modify) | Block 1. Adds `G4SCAN_*` constants, `stage_g4scan`, and the `g4scan` CLI mode. |
| `tests/test_tier2_h2b_forward_model.py` (modify) | Block 1 grid/gate unit tests, matching the existing g3scan test style. |
| `data/runs/h2b_forward_model/h2b_g4_transfer.csv` (generated) | Block 0 output: per-cell paired rows + the fitted bias model + normalizers. |
| `data/runs/h2b_forward_model/h2b_g4scan_*.csv` (generated) | Block 1 output: full scan, gated subset, gated KE rows, ridge summary. |

---

### Task 1: Block 0 — twin ranking authority, seed-SD normalizers, joint score

**Files:**
- Create: `scripts/post_processing/tier2atlas_g4_transfer.py`
- Create: `tests/test_tier2atlas_g4_transfer.py`
- Read-only inputs: `data/runs/h2b_forward_model/atlas_g3ring_table.csv`,
  `data/runs/h2b_forward_model/h2b_g3scan_predictions.csv`,
  `data/runs/h2b_forward_model/h2b_g3scan_chords.csv`, the five battery run
  dirs `9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s{1..5}`

**Interfaces:**

- Consumes: `observable_row(label, run_dir, abundance_ref, ked_ref,
  with_geometry=False) -> dict`, `check_oracle(row) -> list[str]`,
  `RUNS_ROOT`, `ORACLE_RUN`, `ABUNDANCE_REFERENCE_CSV`,
  `IHE_KED_REFERENCE_CSV` (all from
  `scripts.post_processing.tier2atlas_geometry_table`);
  `verify_twin_preregistration()` from `scripts.gen_tier2atlas_g3ring`;
  `write_rows_csv`, `format_table` from
  `scripts.post_processing.tier2_confirmation_score`.
- Produces (imported by Task 2's report and by later blocks):
  - `RANK_LICENSE_RHO = 0.7`
  - `joint_score(w1: float, midhot: float, deepke: float, norm: dict[str, float]) -> float`
  - `nbar_bias_model(rows: list[dict]) -> dict` returning keys
    `{"a": float, "b": float, "resid_sd": float, "r2": float, "n": int}`
    for the linear fit `Δn̄ = a + b · twin_nbar` (Δn̄ = MD − twin, so
    `Δn̄ < 0` means twin-hot)
  - `predict_md_nbar(twin_nbar: float, model: dict) -> float`
  - CSV `h2b_g4_transfer.csv`

- [ ] **Step 1: Write the failing tests**

```python
"""Block 0 (atlas plan §3.5e) unit tests — pure functions, synthetic input."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.post_processing.tier2atlas_g4_transfer import (  # noqa: E402
    joint_score,
    nbar_bias_model,
    predict_md_nbar,
    rank_license,
)


def test_joint_score_reproduces_the_preregistered_standing_value():
    """finc1v725 pooled (W1 0.571, midHot 1.011, deepKE 0.631) -> S = 4.37.

    The §3.5e table is the pre-registration; with the design's provisional
    normalizers (W1 ref 0.571, ln(1.15) for both KE axes) it must reproduce.
    """
    norm = {"w1_ref": 0.571, "midhot_ln": np.log(1.15), "deepke_ln": np.log(1.15)}
    s = joint_score(0.571, 1.011, 0.631, norm)
    assert s == pytest.approx(4.37, abs=0.01)


def test_joint_score_is_zero_at_a_perfect_cell():
    norm = {"w1_ref": 0.571, "midhot_ln": np.log(1.15), "deepke_ln": np.log(1.15)}
    assert joint_score(0.0, 1.0, 1.0, norm) == pytest.approx(0.0)


def test_joint_score_is_symmetric_in_the_ke_ratio():
    """A ratio and its reciprocal are equally wrong — the log form guarantees it."""
    norm = {"w1_ref": 0.571, "midhot_ln": np.log(1.15), "deepke_ln": np.log(1.15)}
    assert joint_score(0.0, 1.25, 1.0, norm) == pytest.approx(
        joint_score(0.0, 1 / 1.25, 1.0, norm)
    )


def test_nbar_bias_model_recovers_a_planted_line():
    rows = [{"twin_nbar": x, "d_nbar": -(0.1 + 0.2 * x)} for x in
            (4.0, 5.0, 6.0, 7.0, 17.0)]
    model = nbar_bias_model(rows)
    assert model["a"] == pytest.approx(-0.1, abs=1e-9)
    assert model["b"] == pytest.approx(-0.2, abs=1e-9)
    assert model["resid_sd"] == pytest.approx(0.0, abs=1e-9)
    assert model["n"] == 5
    # the predictor inverts the fit: MD = twin + (a + b*twin)
    assert predict_md_nbar(5.0, model) == pytest.approx(5.0 - 1.1, abs=1e-9)


def test_rank_license_is_a_hard_threshold_at_0p7():
    twin = [1.0, 2.0, 3.0, 4.0, 5.0]
    assert rank_license(twin, [1.0, 2.0, 3.0, 4.0, 5.0])["licensed"] is True
    assert rank_license(twin, [5.0, 4.0, 3.0, 2.0, 1.0])["licensed"] is False
    lic = rank_license(twin, [1.0, 2.0, 3.0, 4.0, 5.0])
    assert lic["rho"] == pytest.approx(1.0)
    assert lic["n"] == 5


def test_rank_license_ignores_nan_pairs():
    twin = [1.0, 2.0, np.nan, 4.0, 5.0]
    md = [1.0, 2.0, 3.0, np.nan, 5.0]
    assert rank_license(twin, md)["n"] == 3
```

- [ ] **Step 2: Run the tests to verify they fail**

Run:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_tier2atlas_g4_transfer.py -q`
Expected: FAIL — `ModuleNotFoundError` / `ImportError` on
`scripts.post_processing.tier2atlas_g4_transfer`.

- [ ] **Step 3: Implement the module**

Structure (module docstring must state: plan §3.5e Block 0, pure report,
zero MD, oracle-first):

```python
RANK_LICENSE_RHO = 0.7          # pre-registered permission gate (§3.5e)
RANK_OBSERVABLES = ("w1_solv", "midHot", "deepKE")
TWIN_COLUMN = {"w1_solv": "w1_solv", "midHot": "midhot_geo", "deepKE": "deepke"}
BATTERY_RUNS = tuple(
    f"9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s{i}"
    for i in range(1, 6)
)
PROVISIONAL_NORM = {"w1_ref": 0.571, "midhot_ln": float(np.log(1.15)),
                    "deepke_ln": float(np.log(1.15))}


def joint_score(w1, midhot, deepke, norm):
    """§3.5e joint score. Lower is better; 0 = perfect. Convention, not truth."""
    return (w1 / norm["w1_ref"]
            + abs(np.log(midhot)) / norm["midhot_ln"]
            + abs(np.log(deepke)) / norm["deepke_ln"])


def rank_license(twin_vals, md_vals):
    """Spearman rho + a bootstrap CI over the paired cells (NaN pairs dropped)."""
    # -> {"rho", "ci_lo", "ci_hi", "n", "licensed"}; licensed = rho >= RANK_LICENSE_RHO


def nbar_bias_model(rows):
    """OLS fit of d_nbar = a + b * twin_nbar over the paired cells."""


def predict_md_nbar(twin_nbar, model):
    return twin_nbar + model["a"] + model["b"] * twin_nbar
```

`main()` block order, printed with headers, **oracles first**:

1. `verify_twin_preregistration()`, then `check_oracle(observable_row(...ORACLE_RUN...))`.
   On drift: print `*** SCORER-DRIFT ORACLE FAILED — do not read the transfer ***`
   and `return` without reading anything.
2. Load `atlas_g3ring_table.csv` (MD side) and join to
   `h2b_g3scan_predictions.csv` on `(v_c, E_bind_tag, tau_ps, E0_eV)` —
   float-compared with `math.isclose(rel_tol=0, abs_tol=1e-9)` after
   `float()`, matching `well` → `E_bind_tag`. Raise `AssertionError` naming
   the cell if any of the 14 fails to join (a join miss means the grids
   drifted, which invalidates the transfer).
   Cross-check: the joined `twin_nbar`/`twin_n1` must equal the ring table's
   frozen `twin_nbar`/`twin_n1` columns to 1e-6 — this is the G4-P1 half that
   applies to Block 0.
3. Join `h2b_g3scan_chords.csv` on `(v_c, E_bind_tag)` for `t_exit_q50`
   (the residence column the §3.5d bias reading is scaled by) — reported
   beside each cell, not fitted (14 points do not support two predictors).
4. Per observable in `RANK_OBSERVABLES`: `rank_license(twin, md)`; print
   `rho`, CI, `n`, and `LICENSED` / `NOT LICENSED (gate-only)`.
5. `nbar_bias_model(rows)` on all 14 cells, and again excluding `f725`
   (leverage check — f725 sits at twin n̄ 17 while every other cell is 4.4–6.7).
   Print both, and `resid_sd`, which is the band Block 1's corrected gate
   must carry.
6. Seed-SD normalizers: `observable_row` on each of `BATTERY_RUNS`; print the
   per-seed W₁ / midHot / deepKE / χ²_med and their SDs. Emit the measured
   normalizer dict `{"w1_ref": <pooled W1>, "midhot_ln": ..., "deepke_ln": ...}`
   using `ln(1 + SD/mean)` for the KE axes, next to `PROVISIONAL_NORM`.
   **Explicitly print whether the measured deep-bin SD reconciles with the
   §1.2 record "0.0603 ± 0.0033"** — the §3.5e open item.
7. Joint-score table: S for all 14 ring cells + the pooled battery row, under
   both normalizer sets; flag which ring cells beat the pooled battery's S.
8. `write_rows_csv(RUNS_ROOT / "h2b_forward_model" / "h2b_g4_transfer.csv", rows)`.

- [ ] **Step 4: Run the tests to verify they pass**

Run:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_tier2atlas_g4_transfer.py -q`
Expected: PASS (6 tests).

- [ ] **Step 5: Run the report and read the verdicts**

Run:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' scripts/post_processing/tier2atlas_g4_transfer.py`
Expected: oracles PASS; 14/14 cells join; per-observable ρ with CI; the bias
model; the measured normalizers; `h2b_g4_transfer.csv` written.

**Do not proceed to Task 2 without recording which observables are licensed** —
that result parameterizes Block 1's ranking.

- [ ] **Step 6: Commit**

```bash
git add scripts/post_processing/tier2atlas_g4_transfer.py \
        tests/test_tier2atlas_g4_transfer.py \
        data/runs/h2b_forward_model/h2b_g4_transfer.csv
git commit -m "tier2 atlas G4 Block 0: twin ranking authority + score normalizers (zero MD)"
```

---

### Task 2: Block 1 — the fine ridge scan stage

**Files:**
- Modify: `scripts/tier2_h2b_forward_model.py` (add constants after the
  `G3SCAN_*` block near line 2518; add `stage_g4scan` after `stage_g3scan`;
  register the `g4scan` mode in `main()`)
- Modify: `tests/test_tier2_h2b_forward_model.py` (append the G4 grid tests)

**Interfaces:**
- Consumes: `_g3scan_chord_family(v_c, eb_tag, e_bind_ev, m, ens, force_rebuild)`,
  `_g3_corrected_ensemble(m, ref_ke, solv_exp)`, `_g3_s6_oracle`, `g3_score`,
  `fate_map`, `g3_ref_mean_ke`, `load_experiment`, `G3_STANDING`, `TAU_PS`,
  `OUT` — all unchanged.
- Produces: `h2b_g4scan_predictions.csv` (every cell),
  `h2b_g4scan_gated_predictions.csv`, `h2b_g4scan_gated_ke.csv`,
  `h2b_g4scan_ridge.csv` (the gated set with S and the Pareto flag).

- [ ] **Step 1: Write the failing tests**

```python
def test_g4_grids_refine_the_g3_grids_as_exact_sublattices():
    """G4-P1 rests on this: every Step-2 grid point survives in the fine grid."""
    from scripts.tier2_h2b_forward_model import (
        G3SCAN_E0_GRID, G3SCAN_TAU_PS, G4SCAN_E0_GRID, G4SCAN_TAU_PS,
        G4SCAN_VC,
    )
    assert set(G4SCAN_VC) >= {5.0, 5.5, 6.0, 6.5}
    assert set(G4SCAN_TAU_PS) >= {4.8, 6.4}
    fine = {round(float(x), 3) for x in G4SCAN_E0_GRID}
    coarse = {round(float(x), 3) for x in G3SCAN_E0_GRID
              if 0.26 <= float(x) <= 0.42}
    assert coarse <= fine


def test_g4_vc_grid_is_quarter_stepped_inside_the_basin():
    from scripts.tier2_h2b_forward_model import G4SCAN_VC
    vc = sorted(float(v) for v in G4SCAN_VC)
    assert vc == [5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5]


def test_g4_tau_grid_is_0p4_stepped():
    from scripts.tier2_h2b_forward_model import G4SCAN_TAU_PS
    tau = sorted(float(t) for t in G4SCAN_TAU_PS)
    assert tau == [4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8]
    assert max(tau) > 6.55  # the flagged cell must exist and be flagged


def test_g4_gate_uses_the_bias_corrected_md_band_not_the_g3_bracket():
    from scripts.tier2_h2b_forward_model import (
        G3SCAN_GATE_NBAR, G4SCAN_GATE_NBAR, G4SCAN_GATE_N1SOLV,
    )
    assert G4SCAN_GATE_N1SOLV == (0.19, 0.30)      # unchanged, MD-grade transfer
    assert G4SCAN_GATE_NBAR != G3SCAN_GATE_NBAR    # the crude bracket is retired
    assert G4SCAN_GATE_NBAR[0] >= 3.77 - 1e-9      # anchored on the MD band
```

- [ ] **Step 2: Run the tests to verify they fail**

Run:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_tier2_h2b_forward_model.py -q -k g4`
Expected: FAIL — `ImportError: cannot import name 'G4SCAN_VC'`.

- [ ] **Step 3: Add the constants**

```python
# --- G4 Step 1 (plan §3.5e): the fine ridge scan -------------------------
# The G3 Step-2 grids are exact sub-lattices of these (G4-P1 oracle).
G4SCAN_VC = (5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5)
G4SCAN_TAU_PS = (4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4, 6.8)
G4SCAN_E0_GRID = np.round(np.arange(0.26, 0.4201, 0.005), 3)   # 33 values
G4SCAN_EBIND = G3SCAN_EBIND                                    # unchanged
# n1 transfer is MD-grade (<= 0.036 measured), so the n1 band is the MD band.
G4SCAN_GATE_N1SOLV = (0.19, 0.30)
# nbar: the Block-0 bias model replaces the crude [4.4, 7.1] bracket. The gate
# is applied to the PREDICTED MD nbar, widened by the fit's residual SD.
G4SCAN_GATE_NBAR = (3.77, 4.37)
G4SCAN_NBAR_RESID_SD = None   # set from Block 0's h2b_g4_transfer.csv at load
# Block-2 trigger (plan §3.5e): the KE tension is "broken" only if both hold.
G4SCAN_KE_MIDHOT_BAND = (0.85, 1.15)
G4SCAN_KE_DEEP_MIN = 0.60
```

- [ ] **Step 4: Implement `stage_g4scan`**

Block order, mirroring `stage_g3scan`:

1. **Oracles.** `_g3_s6_oracle(m, ref_ke, solv_exp)`; the committed
   `h2b_g3_corrected_row.csv` + `_ke.csv` re-derived string-identically
   (copy the `stage_g3scan` Block-0b code path verbatim — same convention,
   same failure message prefix changed to `g4scan`).
2. **G4-P1.** Re-score the shared sub-lattice and compare to
   `h2b_g3scan_predictions.csv`: for every row of the committed scan whose
   `(v_c, E_bind_tag, tau_ps, E0_eV)` lies in the G4 grids, the freshly scored
   `n1_solv`, `nbar_det`, `w1_solv`, `midhot_geo`, `deepke` must match the
   committed value **exactly as written** (compare the rounded strings, the
   convention the other oracles use). Any mismatch ⇒ `AssertionError` naming
   the cell; the stage aborts before writing anything.
3. **Chord families.** Loop `G4SCAN_VC × G4SCAN_EBIND` through
   `_g3scan_chord_family` — the 12 (v_c ∈ {5.0, 5.5, 6.0, 6.5}) families load
   from their existing npz caches; only the 9 new (5.25, 5.75, 6.25) families
   integrate. Print per family: cache hit/miss, trap, `K655_q50`, elapsed.
4. **Nested scoring.** For each family × `G4SCAN_TAU_PS` × `G4SCAN_E0_GRID`:
   `K = K655 * (TAU_PS / tau)`, `fate_map(ne, K, e0, 1, sig)`, `g3_score(...)`.
   Same `except ValueError` un-scoreable path as `stage_g3scan` (row of NaN,
   never gated). Per row add `tau_flag = int(tau > G3SCAN_TAU_SOURCED)` and
   `ebind_exception`, exactly as Step 2 did.
5. **Gate + score.** `gate_n1` on `n1_solv`; `gate_nbar` on the **predicted MD
   n̄** `predict_md_nbar(nbar_det, model)` where `model` is read from
   `h2b_g4_transfer.csv` (fail loudly with a clear message if that file is
   absent — Block 0 is a hard prerequisite), band widened by
   `G4SCAN_NBAR_RESID_SD`. Record both raw and predicted n̄ in the row.
   For gated rows compute `S = joint_score(...)` with the Block-0 normalizers,
   and the `licensed` flags so an unlicensed axis can be excluded from S.
6. **Ridge summary.** Write `h2b_g4scan_ridge.csv`: the gated cells with S,
   the Pareto-front flag over (W₁, |ln midHot|, |ln deepKE|), and the caveat
   stamps (`tau_flag`, `ebind_exception`). Print:
   - **G4-P2**: whether the gated set at `eb1168` connects (5.5, 4.8) to
     (6.0, 6.4) through 4-neighbour steps on the (v_c, τ) lattice — print
     `CONNECTED RIDGE` / `TWO ISLANDS` with the path or the gap;
   - **G4-P3**: whether any gated cell holds
     `midHot ∈ G4SCAN_KE_MIDHOT_BAND ∧ deepKE ≥ G4SCAN_KE_DEEP_MIN`;
     if none → print `BLOCK 2 TRIGGER FIRED (p_tail scan is authorized)`;
   - **G4-P4**: the best S vs the pooled battery's S (target: S < 4.37 under
     the provisional normalizers, or the Block-0 measured equivalent).
7. Register `g4scan` in `main()`'s mode dispatch alongside `g3scan`, with
   `--force-rebuild` honoured (passed through to `_g3scan_chord_family`).

- [ ] **Step 5: Run the tests to verify they pass**

Run:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest tests/test_tier2_h2b_forward_model.py -q`
Expected: PASS (the existing suite plus the four new G4 tests).

- [ ] **Step 6: Run the stage**

Run:
`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' scripts/tier2_h2b_forward_model.py g4scan`
Expected: oracles PASS, G4-P1 bit-exact, 9 new chord integrations (~15 min),
≈ 5544 cells scored, the four CSVs written, G4-P2/P3/P4 verdicts printed.

- [ ] **Step 7: Commit**

```bash
git add scripts/tier2_h2b_forward_model.py tests/test_tier2_h2b_forward_model.py \
        data/runs/h2b_forward_model/h2b_g4scan_*.csv
git commit -m "tier2 atlas G4 Block 1: fine ridge scan at the corrected geometry (zero MD)"
```

---

### Task 3 (conditional): Block 2 — the p_tail scan

**Trigger:** fires **only** if Task 2 Step 6 printed
`BLOCK 2 TRIGGER FIRED` (no gated cell holds midHot ∈ [0.85, 1.15] ∧
deepKE ≥ 0.60). If the trigger did not fire, this task is **skipped** and the
outcome recorded as "the ridge broke the KE tension without a new knob".

**Files:** Modify `scripts/tier2_h2b_forward_model.py` (`G4SCAN_PTAIL`,
extend `stage_g4scan` behind a `--ptail` flag), modify
`tests/test_tier2_h2b_forward_model.py`.

The chord cache key `_g3scan_chord_tag(v_c, eb_tag)` does **not** carry
`p_tail` — it is currently pinned at `G3_STANDING[2] = −1.0`. This task must
therefore extend the tag to `f"vc{...}_{eb_tag}_pt{...}"` for non-standing
values **while leaving the existing `p_tail = −1` cache files reachable under
their current names** (keep the old tag when `p_tail == −1.0`), so no committed
artifact is invalidated and no chord is re-integrated. A test must assert both:
the legacy tag is unchanged at −1, and distinct p_tail values get distinct tags.

Grid: `p_tail ∈ (−0.5, −1.0, −1.5, −2.0)` at the top two (v_c, E_bind) chords
from `h2b_g4scan_ridge.csv` — 6 new integrations (the −1.0 pair is cached).
Same gate, same score, same G4-P3 test re-evaluated. Output
`h2b_g4scan_ptail.csv`.

---

### Task 4: Block 3 — the MD finalists generator

**Prerequisite:** Task 2 (and Task 3 if it fired) complete. **The six cells are
chosen from `h2b_g4scan_ridge.csv` and cannot be written before it exists** —
this is a data dependency, not an unspecified step. The composition is fixed
now; only the parameter values come from the ridge:

| slot | cell | source |
|---|---|---|
| 1–3 | top-ranked ridge cells by S, Pareto-front members, distinct (v_c, τ) | `h2b_g4scan_ridge.csv` (or 2 + one p_tail cell if Task 3 fired) |
| 4 | `a037` replicate | frozen (v_c 5.5, eb1168, τ 4.8, E₀ 0.37) |
| 5 | `b031` replicate | frozen (v_c 6.0, eb1168, τ 6.4, E₀ 0.31) |
| 6 | off-ridge fail control, pre-registered to miss | nearest non-gated ridge neighbour |

**Files:** Create `scripts/gen_tier2atlas_g4finals.py` (model:
`scripts/gen_tier2atlas_g3ring.py` — same `RING_MATRIX`/`TWIN_ROWS`/
`verify_twin_preregistration` pattern, same safe-relaunch semantics), create
`scripts/post_processing/tier2atlas_g4finals_table.py` (model:
`tier2atlas_g3ring_table.py`).

**Fixed:** N = 1000 per cell, one fresh shared seed (CRN-paired), corrected
geometry, `exclude_all_coupled` interim retained policy with the bound/marginal
decomposition recorded, all standing pins otherwise.

Frozen twin rows for all six cells must be committed **before launch** and
re-verified string-exact by the scorer (the GR-P1 pattern). Predictions
G4-P5 and the acceptance gate are as recorded in plan §3.5e.

This task carries its own `[PROCEED TO IMPLEMENTATION]` — it is the first MD
spend of the stage and must not launch on Task 2's trigger.

---

## Self-review notes

- **Spec coverage:** §3.5e Block 0 → Task 1; Block 1 → Task 2; Block 2 →
  Task 3 (with its trigger); Block 3 → Task 4; Block 4 (pooled battery + the
  G4 adjudications) is deliberately **not** planned here — it is user
  adjudication plus an existing battery pattern, and it starts only after
  Task 4 reports.
- **Open item carried:** the deep-bin KE seed-SD ("0.0603 ± 0.0033" vs pooled
  deepKE 0.631) is resolved by Task 1 Step 3 item 6; until then the score's
  deepKE normalizer is the provisional `ln(1.15)` and every S is printed under
  both normalizer sets.
- **Risk carried from §3.5e:** if Task 1 licenses no observable, Task 2's
  ranking degrades to gate-only and the S column is reported as twin-level
  only — the stage still delivers the ridge map and G4-P2, and the MD budget
  shifts toward the two-stage pattern.
