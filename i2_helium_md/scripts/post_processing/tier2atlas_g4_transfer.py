"""Atlas G4 Step 1 — Block 0: the twin's ranking authority (plan §3.5e).

Pure scorer report — runs nothing, mutates nothing, spends **zero MD**. It
reads three committed artifacts and the five pooled-battery run directories:

1. ``atlas_g3ring_table.csv``      — the 14 MD ring cells (G3 Step 3),
2. ``h2b_g3scan_predictions.csv``  — their frozen twin rows (G3 Step 2),
3. ``h2b_g3scan_chords.csv``       — the per-chord residence column,

and answers the three questions plan §3.5e made prerequisites of the fine
ridge scan:

* **Ranking authority.** §14.4 licenses the twin for n1/nbar/trap only — W1 is
  "bias-loaded" and KE "direction-only", yet those are exactly the observables
  the G4 joint score ranks on. The 14 paired (twin, MD) cells measure the
  transfer for the first time. Pre-registered permission gate: Spearman
  ``rho >= RANK_LICENSE_RHO`` licenses *ranking* on that observable; below it
  the observable is gate-only and MD must rank it. Never used as a correction
  to a twin value.
* **The nbar bias model.** ``d_nbar = a + b * twin_nbar`` replaces the crude
  [+0.3, +3] bracket that forced the Step-2 gate to nbar in [4.4, 7.1]; its
  residual SD is the band the Block-1 gate carries.
* **The score normalizers.** The five N = 1000 battery members give the
  seed-SDs of W1 / midHot / deepKE / chi2_med, resolving the plan §1.2 open
  item (the recorded deep-bin read "0.0603 +/- 0.0033" against the pooled
  deepKE 0.631) and replacing the design's provisional ln(1.15).

Oracles first (§1.4): the frozen twin pre-registration, then the
pooled-battery scorer-drift oracle. Either failing aborts the report before a
single new number is read.

Invocation:

    python scripts/post_processing/tier2atlas_g4_transfer.py
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.gen_tier2atlas_g3ring import (  # noqa: E402
    verify_twin_preregistration,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUNS_ROOT,
    check_oracle,
    observable_row,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

# Pre-registered permission gate (plan §3.5e Block 0).
RANK_LICENSE_RHO = 0.7

# Ranking observables: MD column -> twin column. chi2_med is deliberately
# absent — the twin does not compute a chi-square, and §3.5e reports it
# without scoring it.
TWIN_COLUMN = {
    "w1_solv": "w1_solv",
    "midHot": "midhot_geo",
    "deepKE": "deepke",
}

FORWARD_MODEL_DIR = RUNS_ROOT / "h2b_forward_model"
RING_TABLE_CSV = FORWARD_MODEL_DIR / "atlas_g3ring_table.csv"
TWIN_SCAN_CSV = FORWARD_MODEL_DIR / "h2b_g3scan_predictions.csv"
TWIN_CHORDS_CSV = FORWARD_MODEL_DIR / "h2b_g3scan_chords.csv"
SAVE_CSV_PATH = FORWARD_MODEL_DIR / "h2b_g4_transfer.csv"

# The five members of the pooled N = 5000 battery (findings §4cc).
BATTERY_RUNS = tuple(
    f"9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s{i}"
    for i in range(1, 6)
)

# The design's provisional normalizers (plan §3.5e). Every S is printed under
# these AND under the measured set until the §1.2 deep-bin item is settled.
PROVISIONAL_NORM = {
    "w1_ref": 0.571,
    "midhot_ln": float(np.log(1.15)),
    "deepke_ln": float(np.log(1.15)),
}

# The plan §1.2 record this report checks the measured deep-bin SD against.
RECORDED_DEEP_READ = (0.0603, 0.0033)

# Join tolerance for the frozen twin columns already carried by the ring table.
TWIN_CROSSCHECK_TOL = 1e-6

# ---------------------------------------------------------------------------


def joint_score(
    w1: float,
    midhot: float,
    deepke: float,
    norm: dict[str, float],
    axes: Sequence[str] = ("w1", "midhot", "deepke"),
) -> float:
    """The plan §3.5e joint score. Lower is better; 0 is a perfect cell.

    ``S = W1/w1_ref + |ln midHot|/midhot_ln + |ln deepKE|/deepke_ln``

    Equal weights are a **declared convention, not a truth**: each term is the
    observable's distance from the experiment divided by a scale that makes the
    three commensurable. ``axes`` drops any term Block 0 did not license, so an
    unlicensed observable never silently contributes to a ranking.

    Units: ``w1`` is a Wasserstein distance in He-number bins; ``midhot`` and
    ``deepke`` are dimensionless sim/ref KE ratios (1.0 = experiment). Returns
    NaN if any *scored* axis is NaN.
    """
    total = 0.0
    if "w1" in axes:
        total += w1 / norm["w1_ref"]
    if "midhot" in axes:
        total += abs(math.log(midhot)) / norm["midhot_ln"] if midhot > 0 \
            else float("nan")
    if "deepke" in axes:
        total += abs(math.log(deepke)) / norm["deepke_ln"] if deepke > 0 \
            else float("nan")
    if any(np.isnan(v) for v in (
        w1 if "w1" in axes else 0.0,
        midhot if "midhot" in axes else 1.0,
        deepke if "deepke" in axes else 1.0,
    )):
        return float("nan")
    return float(total)


def _paired(twin_vals: Iterable[float], md_vals: Iterable[float]):
    t = np.asarray(list(twin_vals), dtype=float)
    m = np.asarray(list(md_vals), dtype=float)
    if t.size != m.size:
        raise ValueError("twin/MD vectors must be the same length")
    ok = np.isfinite(t) & np.isfinite(m)
    return t[ok], m[ok]


def rank_license(
    twin_vals: Iterable[float],
    md_vals: Iterable[float],
    n_boot: int = 2000,
    seed: int = 20260728,
) -> dict[str, Any]:
    """Spearman rho over the paired cells + a bootstrap CI and the permission bit.

    NaN pairs are dropped (a cell can be un-scoreable on one axis and fine on
    another). ``licensed`` is the §3.5e gate ``rho >= RANK_LICENSE_RHO`` — a
    *permission* to rank on this observable in the fine scan, never a licence
    to correct a twin value. With n ~ 14 the CI is wide by construction and
    must be reported with the verdict.
    """
    from scipy.stats import spearmanr

    t, m = _paired(twin_vals, md_vals)
    n = int(t.size)
    if n < 3:
        return {"rho": float("nan"), "ci_lo": float("nan"),
                "ci_hi": float("nan"), "n": n, "licensed": False}
    rho = float(spearmanr(t, m).statistic)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if np.unique(t[idx]).size < 3 or np.unique(m[idx]).size < 3:
            continue
        boots.append(float(spearmanr(t[idx], m[idx]).statistic))
    if boots:
        ci_lo, ci_hi = (float(x) for x in np.nanpercentile(boots, [2.5, 97.5]))
    else:
        ci_lo = ci_hi = float("nan")
    return {"rho": rho, "ci_lo": ci_lo, "ci_hi": ci_hi, "n": n,
            "licensed": bool(rho >= RANK_LICENSE_RHO)}


def nbar_bias_model(rows: Sequence[dict[str, Any]]) -> dict[str, float]:
    """OLS fit of ``d_nbar = a + b * twin_nbar`` (d_nbar = MD - twin).

    d_nbar < 0 means the twin runs hot, which is the measured direction at
    14/14 ring cells. Returns ``a``, ``b``, the residual SD (the band the
    Block-1 gate must carry), ``r2`` and the point count ``n``.
    """
    x = np.asarray([float(r["twin_nbar"]) for r in rows], dtype=float)
    y = np.asarray([float(r["d_nbar"]) for r in rows], dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        raise ValueError(
            f"nbar bias model needs at least 3 finite points, got {x.size}"
        )
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    # ddof=2: two fitted parameters.
    resid_sd = float(np.sqrt((resid ** 2).sum() / max(x.size - 2, 1)))
    return {"a": float(a), "b": float(b), "resid_sd": resid_sd,
            "r2": float(r2), "n": int(x.size)}


def predict_md_nbar(twin_nbar: float, model: dict[str, float]) -> float:
    """Predicted MD nbar from a twin nbar under the fitted bias model."""
    return float(twin_nbar + model["a"] + model["b"] * twin_nbar)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required committed artifact missing: {path}")
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _close(a: str | float, b: str | float, tol: float = 1e-9) -> bool:
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)


def join_twin_rows(
    ring_rows: Sequence[dict[str, str]],
    twin_rows: Sequence[dict[str, str]],
    chord_rows: Sequence[dict[str, str]],
) -> list[dict[str, Any]]:
    """Join each MD ring cell to its twin scan row and its chord row.

    Keyed on (v_c, E_bind_tag, tau_ps, E0_eV) for the scan and (v_c,
    E_bind_tag) for the chord. A miss is fatal: it means the committed grids
    drifted apart, which invalidates the whole transfer measurement.
    """
    out: list[dict[str, Any]] = []
    for ring in ring_rows:
        key = (ring["v_c"], ring["well"], ring["tau"], ring["E0"])
        match = [
            t for t in twin_rows
            if t["E_bind_tag"] == key[1]
            and _close(t["v_c"], key[0])
            and _close(t["tau_ps"], key[2])
            and _close(t["E0_eV"], key[3])
        ]
        if len(match) != 1:
            raise AssertionError(
                f"twin join failed for ring cell {ring['label']!r} at "
                f"(v_c {key[0]}, {key[1]}, tau {key[2]}, E0 {key[3]}): "
                f"{len(match)} matching scan rows (expected exactly 1)"
            )
        twin = match[0]
        chords = [
            c for c in chord_rows
            if c["E_bind_tag"] == key[1] and _close(c["v_c"], key[0])
        ]
        if len(chords) != 1:
            raise AssertionError(
                f"chord join failed for ring cell {ring['label']!r}: "
                f"{len(chords)} rows (expected exactly 1)"
            )
        # The ring table carries frozen twin_nbar/twin_n1 columns; they must
        # agree with the scan CSV or one of the two artifacts has moved.
        for ring_col, twin_col in (("twin_nbar", "nbar_det"),
                                   ("twin_n1", "n1_solv")):
            if not _close(ring[ring_col], twin[twin_col], TWIN_CROSSCHECK_TOL):
                raise AssertionError(
                    f"G4-P1 (Block 0 half) FAILED at {ring['label']}: "
                    f"ring {ring_col} {ring[ring_col]} != scan {twin_col} "
                    f"{twin[twin_col]}"
                )
        row: dict[str, Any] = {
            "label": ring["label"],
            "v_c": float(ring["v_c"]),
            "well": ring["well"],
            "tau": float(ring["tau"]),
            "E0": float(ring["E0"]),
            "md_gate": int(ring["md_gate"]),
            "t_exit_q50": float(chords[0]["t_exit_q50"]),
            "twin_nbar": float(twin["nbar_det"]),
            "md_nbar": float(ring["nbar_det"]),
            "twin_n1": float(twin["n1_solv"]),
            "md_n1": float(ring["n1_solv"]),
        }
        row["d_nbar"] = row["md_nbar"] - row["twin_nbar"]
        for md_col, twin_col in TWIN_COLUMN.items():
            row[f"twin_{md_col}"] = float(twin[twin_col])
            row[f"md_{md_col}"] = float(ring[md_col])
        row["md_chi2_med"] = float(ring["chi2_med"])
        out.append(row)
    return out


def battery_normalizers(abundance_ref, ked_ref) -> tuple[dict[str, float],
                                                         list[dict[str, Any]]]:
    """Seed-SDs of the ranking observables from the five N = 1000 members."""
    rows = [
        observable_row(run.split("bigc1v725")[-1], RUNS_ROOT / run,
                       abundance_ref, ked_ref, with_geometry=False)
        for run in BATTERY_RUNS
    ]
    keep = ("label", "num_scored", "w1_solv", "midHot", "deepKE", "chi2_med")
    trimmed = [{k: r[k] for k in keep} for r in rows]
    stats: dict[str, float] = {}
    for col in ("w1_solv", "midHot", "deepKE", "chi2_med"):
        vals = np.asarray([float(r[col]) for r in rows], dtype=float)
        stats[f"{col}_mean"] = float(vals.mean())
        stats[f"{col}_sd"] = float(vals.std(ddof=1))
    norm = {
        "w1_ref": stats["w1_solv_mean"],
        "midhot_ln": float(np.log1p(stats["midHot_sd"] / stats["midHot_mean"])),
        "deepke_ln": float(np.log1p(stats["deepKE_sd"] / stats["deepKE_mean"])),
    }
    return {**norm, **stats}, trimmed


def main() -> None:
    abundance_ref = None
    ked_ref = None

    print("=== G4 Step 1 / Block 0 — twin ranking authority (plan §3.5e) ===")
    print("zero MD; oracles first.\n")

    # ---- Oracles ---------------------------------------------------------
    verify_twin_preregistration()
    from i2_helium_md.postprocess.abundance_loader import (
        load_he_abundance_reference,
    )
    from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference

    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    oracle = observable_row("pooled(oracle)", RUNS_ROOT / ORACLE_RUN,
                            abundance_ref, ked_ref, with_geometry=False)
    drift = check_oracle(oracle)
    if drift:
        print("*** SCORER-DRIFT ORACLE FAILED — do not read the transfer ***")
        for line in drift:
            print(f"    {line}")
        return
    print("scorer-drift oracle OK (pooled battery reproduces).\n")

    # ---- Join ------------------------------------------------------------
    rows = join_twin_rows(_read_csv(RING_TABLE_CSV), _read_csv(TWIN_SCAN_CSV),
                          _read_csv(TWIN_CHORDS_CSV))
    print(f"joined {len(rows)}/{len(rows)} ring cells to their frozen twin "
          f"rows (G4-P1 Block-0 half PASSED: ring twin columns agree with "
          f"the scan CSV).\n")

    # ---- 1. Ranking authority -------------------------------------------
    print("=== Ranking authority (Spearman twin vs MD, 95% bootstrap CI) ===")
    licenses: dict[str, dict[str, Any]] = {}
    for md_col in TWIN_COLUMN:
        lic = rank_license([r[f"twin_{md_col}"] for r in rows],
                           [r[f"md_{md_col}"] for r in rows])
        licenses[md_col] = lic
        verdict = "LICENSED (twin may rank)" if lic["licensed"] \
            else "NOT LICENSED (gate-only; MD must rank this axis)"
        print(f"  {md_col:9s} rho = {lic['rho']:+.3f}  "
              f"CI [{lic['ci_lo']:+.3f}, {lic['ci_hi']:+.3f}]  "
              f"n = {lic['n']:2d}   {verdict}")
    print(f"  (gate: rho >= {RANK_LICENSE_RHO})")

    # ---- 2. nbar bias model ---------------------------------------------
    print("\n=== nbar bias model  d_nbar = a + b * twin_nbar  (MD - twin) ===")
    full = nbar_bias_model(rows)
    no_f725 = nbar_bias_model([r for r in rows if r["label"] != "f725"])
    for name, mdl in (("all 14 cells", full), ("without f725", no_f725)):
        print(f"  {name:14s} a = {mdl['a']:+.4f}  b = {mdl['b']:+.4f}  "
              f"R2 = {mdl['r2']:.3f}  resid SD = {mdl['resid_sd']:.3f} He  "
              f"n = {mdl['n']}")
    print("  (f725 sits at twin nbar ~17 while every other cell is 4.4-6.7 — "
          "the second fit is the leverage check.)")
    print(f"  Block-1 gate band: MD nbar in [3.77, 4.37] widened by the "
          f"residual SD {no_f725['resid_sd']:.3f} He.")
    print("  Residence column reported beside each cell (not fitted — 14 "
          "points do not support a second predictor).")

    # ---- 3. Score normalizers -------------------------------------------
    print("\n=== Seed-SD normalizers from the five N = 1000 battery members ===")
    measured, battery_rows = battery_normalizers(abundance_ref, ked_ref)
    print(format_table(battery_rows))
    for col in ("w1_solv", "midHot", "deepKE", "chi2_med"):
        print(f"  {col:9s} mean {measured[f'{col}_mean']:.4f}   "
              f"seed-SD {measured[f'{col}_sd']:.4f}")
    lo, hi = RECORDED_DEEP_READ
    deep_mean, deep_sd = measured["deepKE_mean"], measured["deepKE_sd"]
    reconciles = abs(deep_sd - hi) <= 0.5 * hi and abs(deep_mean - lo) <= 0.5 * lo
    print(f"\n  plan §1.2 open item: recorded deep-bin read "
          f"{lo} +/- {hi} vs measured deepKE {deep_mean:.4f} +/- "
          f"{deep_sd:.4f} -> "
          f"{'RECONCILES' if reconciles else 'DOES NOT RECONCILE'}")
    if not reconciles:
        print("  => the §1.2 figure is NOT the deep-bin KE ratio's mean/SD; "
              "the measured seed-SD above is what the G4 score normalizes by, "
              "and §1.2 needs a provenance note (recorded, not silently "
              "overwritten).")
    measured_norm = {k: measured[k] for k in
                     ("w1_ref", "midhot_ln", "deepke_ln")}
    print(f"\n  provisional norm: {PROVISIONAL_NORM}")
    print(f"  measured    norm: {measured_norm}")

    # ---- 4. Joint score --------------------------------------------------
    axes = tuple(
        short for short, md_col in (("w1", "w1_solv"), ("midhot", "midHot"),
                                    ("deepke", "deepKE"))
        if licenses[md_col]["licensed"]
    )
    print(f"\n=== Joint score S (lower better; licensed axes: "
          f"{axes or '(none)'}) ===")
    battery_ref = {
        "w1_solv": measured["w1_solv_mean"], "midHot": measured["midHot_mean"],
        "deepKE": measured["deepKE_mean"],
    }
    score_rows: list[dict[str, Any]] = []
    for r in rows + [{"label": "pooled battery", "md_gate": "-",
                      "md_w1_solv": battery_ref["w1_solv"],
                      "md_midHot": battery_ref["midHot"],
                      "md_deepKE": battery_ref["deepKE"]}]:
        s_prov = joint_score(r["md_w1_solv"], r["md_midHot"], r["md_deepKE"],
                             PROVISIONAL_NORM)
        s_meas = joint_score(r["md_w1_solv"], r["md_midHot"], r["md_deepKE"],
                             measured_norm)
        s_lic = joint_score(r["md_w1_solv"], r["md_midHot"], r["md_deepKE"],
                            measured_norm, axes=axes) if axes else float("nan")
        score_rows.append({
            "label": r["label"], "md_gate": r["md_gate"],
            "W1": round(float(r["md_w1_solv"]), 4),
            "midHot": round(float(r["md_midHot"]), 4),
            "deepKE": round(float(r["md_deepKE"]), 4),
            "S_provisional": round(s_prov, 3),
            "S_measured": round(s_meas, 3),
            "S_licensed_axes": round(s_lic, 3) if axes else float("nan"),
        })
    print(format_table(score_rows))
    incumbent = score_rows[-1]
    better = [r["label"] for r in score_rows[:-1]
              if r["S_provisional"] < incumbent["S_provisional"]]
    print(f"\n  pooled-battery S_provisional = "
          f"{incumbent['S_provisional']:.3f} (the §3.5e 4.37 target)")
    print(f"  ring cells beating it: {better or 'NONE'}")

    # ---- 5. CSV ----------------------------------------------------------
    for r in rows:
        r["S_provisional"] = round(
            joint_score(r["md_w1_solv"], r["md_midHot"], r["md_deepKE"],
                        PROVISIONAL_NORM), 4)
        r["S_measured"] = round(
            joint_score(r["md_w1_solv"], r["md_midHot"], r["md_deepKE"],
                        measured_norm), 4)
        r["nbar_bias_a"] = round(no_f725["a"], 6)
        r["nbar_bias_b"] = round(no_f725["b"], 6)
        r["nbar_bias_resid_sd"] = round(no_f725["resid_sd"], 6)
        r["pred_md_nbar"] = round(predict_md_nbar(r["twin_nbar"], no_f725), 4)
        for md_col in TWIN_COLUMN:
            r[f"rho_{md_col}"] = round(licenses[md_col]["rho"], 4)
            r[f"licensed_{md_col}"] = int(licenses[md_col]["licensed"])
        r["norm_w1_ref"] = round(measured_norm["w1_ref"], 6)
        r["norm_midhot_ln"] = round(measured_norm["midhot_ln"], 6)
        r["norm_deepke_ln"] = round(measured_norm["deepke_ln"], 6)
    print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
