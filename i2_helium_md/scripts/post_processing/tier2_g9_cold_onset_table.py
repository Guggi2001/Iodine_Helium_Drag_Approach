"""G9 scorer — cold-onset minority mixtures at the lq A5 cell (zero MD).

THESIS_FINAL_GOALS.md G9 (frozen 2026-10-08). Columns: the seed-s1 E0 grid
of ``gen_tier2_g9_cold_onset.py`` (k = -4 .. +2, step 0.0125 eV; k = 0 =
the committed A5 s1 member). All columns share births (CRN).

Mixture (exact pooling): molecule fractions ``p_k`` over the grid; both
fragments of a molecule come from one column (P1). Molecules are assigned
by a dedicated permutation (``np.random.default_rng(ASSIGN_SEED)``, not the
simulation RNG; largest-remainder counts); the mixed read is the pooled
``include_mask`` reads of the columns.

* main peak = blend of adjacent nodes (j, j+1), fraction f on j+1, chosen so
  the **expected** n̄ (sum of column n-sums / sum of column scored counts,
  weighted by p) = 3.89;
* shape P: fraction w at main - 0.05 eV (the blend shifted 4 nodes down);
* shape U: fraction w uniform over [main - 0.05, main] (trapezoid over the 5
  shifted nodes, weights 1/8, 1/4, 1/4, 1/4, 1/8);
* w in {0, 0.05, ..., 0.50}. A row whose nodes leave the grid is reported
  ``off_grid`` and not scored.

Anchors: the A5 s1 column rescored = committed ``c200t64confs1`` (4 dp);
the pure-k=0 mixture through the mask path reproduces the same row.

Marked per shape (reported, nothing adopted): min W1 and max n >= 15 among
rows with midHot in [0.85, 1.15] and KE1 >= (w = 0 row) - 0.01.

Outputs (committed): ``data/runs/h2b_forward_model/g9_columns.csv``,
``g9_mixtures.csv``.

    python scripts/post_processing/tier2_g9_cold_onset_table.py
"""

from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path
from typing import Any, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

# =============================================================================
# Pre-registered (G9; FIXED)
# =============================================================================
NBAR_TARGET = 3.89
COLD_SHIFT_NODES = 4                              # 0.05 eV / 0.0125 eV
TRAPEZOID = (1 / 8, 1 / 4, 1 / 4, 1 / 4, 1 / 8)   # shifts 0, -1, -2, -3, -4 nodes
W_GRID = tuple(round(0.05 * i, 2) for i in range(11))
SHAPES = ("P", "U")
ASSIGN_SEED = 20261008
MIDHOT_BAND = (0.85, 1.15)
KE1_TOL = 0.01
DEEP_N = (10, 14)
MIN_BIN_IONS = 5


def cold_weights(shape: str) -> dict[int, float]:
    """Node-shift -> weight of the cold part (shift <= 0, in grid nodes)."""
    if shape == "P":
        return {-COLD_SHIFT_NODES: 1.0}
    if shape == "U":
        return {-s: wt for s, wt in enumerate(TRAPEZOID)}
    raise ValueError(f"unknown shape {shape!r}; expected one of {SHAPES}")


def mixture_weights(j: int, f: float, w: float, shape: str) -> dict[int, float]:
    """Molecule fractions per grid node k for main blend (j, j+1; f on j+1),
    cold fraction ``w`` of ``shape``. Weights sum to 1."""
    main = {j: 1.0 - f, j + 1: f}
    p: dict[int, float] = {}
    for k, pk in main.items():
        p[k] = p.get(k, 0.0) + (1.0 - w) * pk
        for s, ws in cold_weights(shape).items():
            p[k + s] = p.get(k + s, 0.0) + w * pk * ws
    return {k: v for k, v in p.items() if v > 0.0}


def expected_nbar(p: dict[int, float], n_sum: dict[int, float],
                  n_scored: dict[int, float]) -> float:
    """Expected mixture n̄ = sum p_k S_k / sum p_k C_k (S = column n-sum over
    scored ions, C = column scored count, both per N molecules)."""
    return (sum(pk * n_sum[k] for k, pk in p.items())
            / sum(pk * n_scored[k] for k, pk in p.items()))


def solve_main_blend(w: float, shape: str, nodes: list[int], n_sum: dict[int, float],
                     n_scored: dict[int, float], target: float = NBAR_TARGET
                     ) -> Optional[tuple[int, float]]:
    """(j, f) putting the expected n̄ at ``target``, or None if no segment
    with every node on the grid brackets it. n̄(f) is linear-fractional on a
    segment, so f solves exactly."""
    lo, hi = min(nodes), max(nodes)
    for j in range(lo, hi):
        p0 = mixture_weights(j, 0.0, w, shape)
        p1 = mixture_weights(j, 1.0, w, shape)
        if any(k < lo or k > hi for k in list(p0) + list(p1)):
            continue
        a = sum(pk * n_sum[k] for k, pk in p0.items())
        d = sum(pk * n_scored[k] for k, pk in p0.items())
        b = sum(pk * n_sum[k] for k, pk in p1.items()) - a
        e = sum(pk * n_scored[k] for k, pk in p1.items()) - d
        denom = b - target * e
        if denom == 0.0:
            continue
        f = (target * d - a) / denom
        if -1e-12 <= f <= 1.0 + 1e-12:
            return j, float(min(1.0, max(0.0, f)))
    return None


def assign_molecules(p: dict[int, float], num_molecules: int,
                     seed: int = ASSIGN_SEED) -> dict[int, np.ndarray]:
    """Molecule index arrays per node: largest-remainder counts of
    ``p * num_molecules`` over one fixed permutation (common to every row)."""
    keys = sorted(p)
    raw = np.array([p[k] * num_molecules for k in keys])
    counts = np.floor(raw).astype(int)
    short = num_molecules - int(counts.sum())
    order = np.argsort(-(raw - counts), kind="stable")
    counts[order[:short]] += 1
    perm = np.random.default_rng(seed).permutation(num_molecules)
    out, start = {}, 0
    for k, c in zip(keys, counts):
        out[k] = np.sort(perm[start:start + c])
        start += c
    return out


def fragment_mask(molecules: np.ndarray, num_molecules: int) -> np.ndarray:
    """Per-ion mask (2N,) selecting both fragments (i and i + N) of each molecule."""
    mask = np.zeros(2 * num_molecules, dtype=bool)
    mask[molecules] = True
    mask[molecules + num_molecules] = True
    return mask


def _write(rows: list[dict[str, Any]], path: Path) -> None:
    keys: list[str] = []
    for r in rows:
        keys.extend(k for k in r if k not in keys)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference
    from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference
    from i2_helium_md.postprocess.tier2_confirmation import (
        pool_confirmation_reads, read_confirmation_detection, solvated_renormalized,
    )
    from i2_helium_md.simulation.detection_stage import load_detection_result
    from i2_helium_md.simulation.run_directory import RunDirectory
    from scripts.gen_tier2_g9_cold_onset import columns, run_dir_for
    from scripts.post_processing.tier2_detfix_battery_table import (
        _check_4dp, _csv_row, _pipeline, _row,
    )
    from scripts.post_processing.tier2_g413_table import ANCHOR_COLS
    from scripts.post_processing.tier2_g414_table import band_columns
    from scripts.post_processing.tier2atlas_geometry_table import (
        ABUNDANCE_REFERENCE_CSV, IHE_KED_REFERENCE_CSV, RUNS_ROOT,
    )

    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warnings
    out_dir = RUNS_ROOT / "h2b_forward_model"
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    exp_n, exp_f = solvated_renormalized(ab.n, ab.ion_fraction, n_min=1)
    exp_ke = dict(zip(np.asarray(ked.n, dtype=int).tolist(),
                      np.asarray(ked.mean_KE_eV, dtype=float).tolist()))

    cols = columns()
    missing = [c.label for c in cols if not (run_dir_for(c) / "detection.npz").exists()]
    if missing:
        print(f"columns not yet run: {missing}")
        return 1
    dets, cfg0 = {}, None
    for c in cols:
        rd = run_dir_for(c)
        cfg = RunDirectory(rd).load_cfg()
        if abs(cfg.internal_energy_partition_fraction * 2.70 - c.e0_eV) > 1e-12:
            raise AssertionError(f"[{c.label}] stored E0 != grid E0")
        dets[c.k] = load_detection_result(rd / "detection.npz")
        cfg0 = cfg0 or cfg
    num_mol = int(dets[0].num_molecules)

    def extra(read) -> dict[str, Any]:
        n, ke = np.asarray(read.n_scored), np.asarray(read.ke_scored_eV)
        out: dict[str, Any] = {"max_n": int(n.max())}
        ratios = [ke[n == m].mean() / exp_ke[m] for m in range(DEEP_N[0], DEEP_N[1] + 1)
                  if (n == m).sum() >= MIN_BIN_IONS and m in exp_ke]
        out["KE_n10_14_over_exp"] = float(np.mean(ratios)) if ratios else float("nan")
        for m in (1, 2, 3, 4):
            sel = ke[n == m]
            q1, med, q3 = np.percentile(sel, [25, 50, 75]) if sel.size else (np.nan,) * 3
            out[f"KE_iqr_over_med_n{m}"] = float((q3 - q1) / med) if sel.size else np.nan
        return out

    def score(label: str, read) -> dict[str, Any]:
        sn, sf = solvated_renormalized(read.n_values, read.fraction, n_min=1)
        row = {**_row(label, "s1", _pipeline(cfg0), cfg0.seed, read, ab, ked),
               **band_columns(sn, sf, exp_n, exp_f), **extra(read)}
        row["KE2_over_KE1"] = float(row["KE2_mean"]) / float(row["KE1_mean"])
        return row

    # 1. columns + anchors
    committed = _csv_row(out_dir / "g414_table.csv", "c200t64confs1")
    col_rows, n_sum, n_scored = [], {}, {}
    for c in cols:
        r = read_confirmation_detection(dets[c.k], label=c.label)
        row = {**score(c.label, r), "k": c.k, "E0_eV": c.e0_eV}
        if c.k == 0:
            _check_4dp("[anchor] A5 s1 column", row, committed, ANCHOR_COLS)
            print("ANCHOR PASSED: A5 s1 column = committed c200t64confs1 (4 dp)")
        n_sum[c.k] = float(np.asarray(r.n_scored).sum())
        n_scored[c.k] = float(r.num_scored)
        col_rows.append(row)
        print(f"[{c.label}] E0 {c.e0_eV:.4f}: nbar {row['nbar_det']:.3f} n>=15 "
              f"{row['mass_n15plus']:.4f} KE1 {row['KE1_mean']:.4f} W1 {row['w1_solv']:.3f}")
    nb = [r["nbar_det"] for r in col_rows]
    print("cross-check n̄ monotone (non-increasing in E0):",
          all(a >= b - 0.1 for a, b in zip(nb, nb[1:])))

    def mixed_read(p: dict[int, float], label: str):
        groups = assign_molecules(p, num_mol)
        reads = [read_confirmation_detection(dets[k], label=f"{label}_k{k}",
                                             include_mask=fragment_mask(m, num_mol))
                 for k, m in groups.items() if m.size]
        return pool_confirmation_reads(reads, label=label)

    pure = score("mask_k0", mixed_read({0: 1.0}, "mask_k0"))
    _check_4dp("[anchor] mask path pure k=0", pure, committed, ANCHOR_COLS)
    print("ANCHOR PASSED: mask-path pure k=0 mixture = committed c200t64confs1 (4 dp)")

    # 2. mixtures
    nodes = [c.k for c in cols]
    e0_of = {c.k: c.e0_eV for c in cols}
    mix_rows = []
    for shape in SHAPES:
        for w in W_GRID:
            sol = solve_main_blend(w, shape, nodes, n_sum, n_scored)
            base = {"shape": shape, "w": w}
            if sol is None:
                mix_rows.append({**base, "status": "off_grid"})
                continue
            j, f = sol
            p = mixture_weights(j, f, w, shape)
            label = f"g9_{shape}_w{round(w * 100):02d}"
            row = score(label, mixed_read(p, label))
            row.update({**base, "status": "ok",
                        "main_E0_eV": (1 - f) * e0_of[j] + f * e0_of[j + 1],
                        "nbar_expected": expected_nbar(p, n_sum, n_scored),
                        "weights": ";".join(f"{e0_of[k]:.4f}:{v:.4f}" for k, v in sorted(p.items()))})
            mix_rows.append(row)
    for shape in SHAPES:
        ok = [r for r in mix_rows if r["shape"] == shape and r["status"] == "ok"]
        ref = next((r for r in ok if r["w"] == 0.0), None)
        if ref is None:
            continue
        good = [r for r in ok if MIDHOT_BAND[0] <= r["midHot"] <= MIDHOT_BAND[1]
                and r["KE1_mean"] >= ref["KE1_mean"] - KE1_TOL]
        for r in ok:
            r["eligible"] = int(r in good)
            r["mark"] = ""
        if good:
            min(good, key=lambda r: r["w1_solv"])["mark"] += "minW1 "
            max(good, key=lambda r: r["mass_n15plus"])["mark"] += "maxN15 "
        for r in ok:
            print(f"{shape} w {r['w']:.2f}: main E0 {r['main_E0_eV']:.4f} nbar "
                  f"{r['nbar_det']:.3f} n>=15 {r['mass_n15plus']:.4f} max n {r['max_n']} "
                  f"W1 {r['w1_solv']:.3f} deepKE {r['deepKE']:.3f} KE1 {r['KE1_mean']:.4f} "
                  f"midHot {r['midHot']:.3f} {r['mark']}")
    _write(col_rows, out_dir / "g9_columns.csv")
    _write(mix_rows, out_dir / "g9_mixtures.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
