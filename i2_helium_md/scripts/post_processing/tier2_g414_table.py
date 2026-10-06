"""G4.14 scorer — ridge-law MD tau lines: surface fit, residual criterion, selection.

THESIS_FINAL_GOALS.md G4.14 (frozen 7726e43). Pure scorer: runs no MD. Same
observable vector as the production battery scorer (``_row`` from
``tier2_detfix_battery_table.py``).

Order (oracles first; any failure stops before a new number):

1. **Anchors** — h405 detfix ``s1`` rescored = committed ``s1_new`` (4 dp);
   the G4.13 ``h405g35`` s1 run rescored = its committed ``g413_table.csv``
   row (4 dp).
2. **h405g35 pool** (when s2-s5 + r6 exist) — pooled N = 6000 with s1: the
   3.5 A reference row.
3. **Step 1 rows** — every planned run (``gen_tier2_g414_tauline.plan``).
4. **Fit** (per c, once all its planned runs exist) — least squares in
   tau and dE = (E0 - E0_3.85(tau)) / 0.01 eV, E0_3.85 from the committed
   Step-0 line: cubic in tau (quadratic for the 10-run c) + dE + tau*dE;
   n1 additionally + dE^2. **Residual criterion:** residual sd (ddof = n - p)
   <= 1.5 x per-member seed sd (committed ``g413_transfer.csv`` row
   ``seed_std_t64e310``) for n1, nbar, KE1, W1, midHot. A failing
   observable with tau structure (nested F-test for one more tau order,
   p < 0.05) gets that order once; without tau structure the c needs a
   second independent-seed set; still failing after the extra order ->
   *unresolved* (no Step 2). chi2 / deepKE / trap are fitted, reported only.
   c = 2.0 out-of-sample: the fit's mean residual on the 16 G4.13 Stage-A s1
   cells (an estimate of s1's offset; never fit data).
5. **Selection** (per resolved c) — along the fitted n̄ = 3.89 curve,
   tau within the c's MD range at step 0.01, dE inside the sampled range:
   feasible iff predicted n1 >= 0.200 and midHot in [0.85, 1.15]
   (n̄ = 3.89 lies in [3.83, 4.31] by construction). tau_KE = max KE1,
   tau_W = min W1, merged (to tau_KE) if <= 0.2 ps apart -> Step-2
   candidates. Step 2 itself is NOT run from here.

Outputs (committed; every column) in ``data/runs/h2b_forward_model/``:
``g414_table.csv``, ``g414_fit.csv``, ``g414_tradeoff.csv``,
``g414_candidates.csv``.

    python scripts/post_processing/tier2_g414_table.py
"""

from __future__ import annotations

import csv
import math
import sys
import warnings
from pathlib import Path
from typing import Any, NamedTuple, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

# =============================================================================
# Pre-registered (G4.14; FIXED)
# =============================================================================
FIT_OBS = ("n1_solv", "nbar_det", "KE1_mean", "w1_solv", "midHot")
REPORT_OBS = ("chi2_med", "deepKE", "trap")
SEED_SD_COLS = {"n1_solv": "seed_n1_solv", "nbar_det": "seed_nbar_det",
                "KE1_mean": "seed_KE1_mean", "w1_solv": "seed_w1_solv",
                "midHot": "seed_midHot"}
RESID_FACTOR = 1.5
F_ALPHA = 0.05
TAU_CENTER_PS = 5.6
DE_UNIT_EV = 0.01
FULL_LINE_ORDER = 3           # cubic in tau (c = 2.0, 24 runs)
SHORT_LINE_ORDER = 2          # quadratic (10-run c)
SHORT_LINE_MAX_TAUS = 5
N1_MIN = 0.200
NBAR_LANDING = (3.83, 4.31)
MIDHOT_BAND = (0.85, 1.15)
NBAR_MATCH = 3.89
SELECT_TAU_STEP_PS = 0.01
MERGE_PS = 0.2


# =============================================================================
# Pure fit / selection helpers (unit-tested in tests/test_tier2_g414.py)
# =============================================================================
class SurfaceFit(NamedTuple):
    beta: np.ndarray
    tau_order: int
    with_de2: bool
    n: int
    p: int
    rss: float
    resid_sd: float
    resid: np.ndarray


def design(tau, de, tau_order: int, with_de2: bool) -> np.ndarray:
    """Columns: 1, t, ..., t^order, dE, t*dE[, dE^2]; t = tau - 5.6 ps."""
    t = np.asarray(tau, float) - TAU_CENTER_PS
    de = np.asarray(de, float)
    cols = [t ** k for k in range(tau_order + 1)] + [de, t * de]
    if with_de2:
        cols.append(de * de)
    return np.column_stack(cols)


def fit_surface(tau, de, y, tau_order: int, with_de2: bool) -> SurfaceFit:
    """Ordinary least squares; raises if there are no residual dof."""
    X = design(tau, de, tau_order, with_de2)
    y = np.asarray(y, float)
    n, p = X.shape
    if n <= p:
        raise ValueError(f"{n} runs cannot fit {p} parameters")
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    rss = float(resid @ resid)
    return SurfaceFit(beta, tau_order, with_de2, n, p, rss, math.sqrt(rss / (n - p)), resid)


def predict(fit: SurfaceFit, tau, de) -> np.ndarray:
    return design(tau, de, fit.tau_order, fit.with_de2) @ fit.beta


def tau_structure_p(tau, de, y, fit: SurfaceFit) -> float:
    """Nested F-test p-value for adding one tau order to ``fit``."""
    from scipy.stats import f as f_dist
    big = fit_surface(tau, de, y, fit.tau_order + 1, fit.with_de2)
    df2 = big.n - big.p
    if big.rss <= 0.0:
        return 0.0
    F = (fit.rss - big.rss) / (big.rss / df2)
    return float(f_dist.sf(F, 1, df2))


def fit_observable(tau, de, y, base_order: int, with_de2: bool,
                   seed_sd: float) -> dict[str, Any]:
    """The frozen residual-criterion ladder for one observable.

    Returns ``{"fit", "status", "threshold", "p_tau", "order"}``; status is
    ``pass`` | ``pass_after_order`` | ``second_seed`` | ``unresolved``.
    """
    thr = RESID_FACTOR * seed_sd
    fit = fit_surface(tau, de, y, base_order, with_de2)
    if fit.resid_sd <= thr:
        return {"fit": fit, "status": "pass", "threshold": thr, "p_tau": float("nan"),
                "order": base_order}
    p_tau = tau_structure_p(tau, de, y, fit)
    if p_tau >= F_ALPHA:
        return {"fit": fit, "status": "second_seed", "threshold": thr, "p_tau": p_tau,
                "order": base_order}
    fit2 = fit_surface(tau, de, y, base_order + 1, with_de2)
    status = "pass_after_order" if fit2.resid_sd <= thr else "unresolved"
    return {"fit": fit2, "status": status, "threshold": thr, "p_tau": p_tau,
            "order": base_order + 1}


def c_status(statuses: list[str]) -> str:
    """Per-c outcome over the gated observables."""
    if "unresolved" in statuses:
        return "unresolved"
    if "second_seed" in statuses:
        return "second_seed"
    return "resolved"


def de_at_nbar(nbar_fit: SurfaceFit, tau, target: float = NBAR_MATCH) -> np.ndarray:
    """dE where the fitted nbar equals ``target`` (nbar is linear in dE)."""
    if nbar_fit.with_de2:
        raise ValueError("the nbar model carries no dE^2 term")
    tau = np.asarray(tau, float)
    X0 = design(tau, np.zeros_like(tau), nbar_fit.tau_order, False)
    base = X0 @ nbar_fit.beta
    t = tau - TAU_CENTER_PS
    k = nbar_fit.tau_order + 1
    slope = nbar_fit.beta[k] + nbar_fit.beta[k + 1] * t
    return (target - base) / slope


def tradeoff_curve(fits: dict[str, SurfaceFit], tau_range, de_range) -> list[dict]:
    """Predicted observables along the nbar = 3.89 curve, with feasibility."""
    taus = np.round(np.arange(tau_range[0], tau_range[1] + 1e-9, SELECT_TAU_STEP_PS), 2)
    de = de_at_nbar(fits["nbar_det"], taus)
    rows = []
    for t, d in zip(taus, de):
        r = {"tau_ps": float(t), "dE": float(d),
             "in_sampled_dE": int(de_range[0] <= d <= de_range[1])}
        for obs, f in fits.items():
            r[f"pred_{obs}"] = float(predict(f, [t], [d])[0])
        r["feasible"] = int(bool(
            r["in_sampled_dE"] and r["pred_n1_solv"] >= N1_MIN
            and NBAR_LANDING[0] <= r["pred_nbar_det"] <= NBAR_LANDING[1]
            and MIDHOT_BAND[0] <= r["pred_midHot"] <= MIDHOT_BAND[1]))
        rows.append(r)
    return rows


def select_candidates(curve: list[dict]) -> list[dict]:
    """tau_KE (max KE1) and tau_W (min W1) over feasible rows; merged to
    tau_KE if <= 0.2 ps apart. Empty if no feasible row."""
    feas = [r for r in curve if r["feasible"]]
    if not feas:
        return []
    ke = max(feas, key=lambda r: r["pred_KE1_mean"])
    w = min(feas, key=lambda r: r["pred_w1_solv"])
    if abs(ke["tau_ps"] - w["tau_ps"]) <= MERGE_PS + 1e-9:
        return [{**ke, "role": "KE+W"}]
    return [{**ke, "role": "KE"}, {**w, "role": "W"}]


def e0_385_interp(line_rows: list[dict], c: float):
    """tau -> Step-0 E0 at predicted MD nbar 3.85 (linear in tau)."""
    pts = sorted((float(r["tau_ps"]), float(r["E0_eV"])) for r in line_rows
                 if float(r["c"]) == c and float(r["nbar_target_md"]) == 3.85
                 and math.isfinite(float(r["E0_eV"])))
    t, e = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])

    def f(tau):
        tau = np.asarray(tau, float)
        if np.any(tau < t[0]) or np.any(tau > t[-1]):
            raise ValueError("tau outside the Step-0 line")
        return np.interp(tau, t, e)
    return f


# =============================================================================
# Driver
# =============================================================================
def _write(rows: list[dict[str, Any]], path: Path) -> None:
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {path}")


def main() -> int:
    from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference
    from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference
    from i2_helium_md.postprocess.tier2_confirmation import (
        pool_confirmation_reads, read_confirmation_detection,
    )
    from i2_helium_md.simulation.detection_stage import load_detection_result
    from i2_helium_md.simulation.run_directory import RunDirectory
    from scripts.gen_tier2_detfix_battery import member_run_dir_name
    from scripts.gen_tier2_g413_ridge_battery import (
        H405_GATE_LABEL, STAGE_A_MEMBER, STAGE_B_MEMBERS, run_dir_for as g413_run_dir_for,
    )
    from scripts.gen_tier2_g414_tauline import LINE_CSV, load_plan, run_dir_for
    from scripts.post_processing.tier2_detfix_battery_table import (
        TABLE_CSV as DETFIX_TABLE_CSV, _check_4dp, _csv_row, _pipeline, _row,
    )
    from scripts.post_processing.tier2_g413_table import (
        ANCHOR_COLS, TABLE_CSV as G413_TABLE_CSV, TRANSFER_CSV,
    )
    from scripts.post_processing.tier2atlas_geometry_table import (
        ABUNDANCE_REFERENCE_CSV, IHE_KED_REFERENCE_CSV, RUNS_ROOT,
    )

    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warnings
    out_dir = RUNS_ROOT / "h2b_forward_model"
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    def read(rd: Path, label: str):
        return read_confirmation_detection(load_detection_result(rd / "detection.npz"),
                                           label=label)

    # 1. anchors
    rows: list[dict[str, Any]] = []
    h405_s1 = RUNS_ROOT / member_run_dir_name("s1")
    cfg = RunDirectory(h405_s1).load_cfg()
    anchor = _row("h405_s1", "s1", _pipeline(cfg), cfg.seed, read(h405_s1, "h405_s1"), ab, ked)
    _check_4dp("[anchor] h405 detfix s1", anchor, _csv_row(DETFIX_TABLE_CSV, "s1_new"),
               ANCHOR_COLS)
    g35_dir = g413_run_dir_for(H405_GATE_LABEL, STAGE_A_MEMBER)
    g35_reads = [read(g35_dir, "h405g35_s1")]
    c35 = RunDirectory(g35_dir).load_cfg()
    g35 = _row("h405g35_s1", "s1", _pipeline(c35), c35.seed, g35_reads[0], ab, ked)
    _check_4dp("[anchor] h405g35 s1", g35, _csv_row(G413_TABLE_CSV, "h405g35_s1"), ANCHOR_COLS)
    print("ANCHORS PASSED: h405 s1 = s1_new; h405g35 s1 = committed g413 row (4 dp)")
    ref6000 = _csv_row(DETFIX_TABLE_CSV, "pooled_new6000")
    rows.append({**anchor, "kind": "ref"})
    rows.append({**{k: ref6000[k] for k in ref6000}, "kind": "ref_pooled_14p2"})
    rows.append({**g35, "kind": "ref_member_3p5"})

    # 2. h405g35 pool
    g35_members = [g413_run_dir_for(H405_GATE_LABEL, m) for m in STAGE_B_MEMBERS]
    g35_pooled = None
    if all((d / "detection.npz").exists() for d in g35_members):
        for m, d in zip(STAGE_B_MEMBERS, g35_members):
            c_ = RunDirectory(d).load_cfg()
            r_ = read(d, f"h405g35_{m}")
            g35_reads.append(r_)
            rows.append({**_row(f"h405g35_{m}", m, _pipeline(c_), c_.seed, r_, ab, ked),
                         "kind": "ref_member_3p5"})
        g35_pooled = _row("h405g35_pooled6000", "pooled6000", "skip+closure:partner_aware",
                          "-", pool_confirmation_reads(g35_reads, label="h405g35_pooled6000"),
                          ab, ked)
        rows.append({**g35_pooled, "kind": "ref_pooled_3p5"})
        print(f"h405g35 pooled N=6000: n1 {g35_pooled['n1_solv']:.4f} nbar "
              f"{g35_pooled['nbar_det']:.3f} KE1 {g35_pooled['KE1_mean']:.4f} "
              f"W1 {g35_pooled['w1_solv']:.4f} midHot {g35_pooled['midHot']:.4f}")
    else:
        print("h405g35 pool: members not yet run")

    # 3. Step-1 rows
    plan = load_plan()
    with open(LINE_CSV, newline="", encoding="utf-8") as fh:
        line_rows = list(csv.DictReader(fh))
    step1: dict[float, list[dict]] = {}
    for cell in plan:
        rd = run_dir_for(cell)
        if not (rd / "detection.npz").exists():
            continue
        c_ = RunDirectory(rd).load_cfg()
        if _pipeline(c_) != "skip+closure:partner_aware" or c_.seed != cell.seed:
            raise AssertionError(f"[{cell.label}] pipeline/seed mismatch")
        r = _row(cell.label, "g414", _pipeline(c_), c_.seed, read(rd, cell.label), ab, ked)
        e385 = e0_385_interp(line_rows, cell.c)
        r.update({"kind": "step1", "c": cell.c, "tau_ps": cell.tau_ps, "E0_eV": cell.e0_eV,
                  "nbar_target": cell.nbar_target,
                  "dE": (cell.e0_eV - float(e385(cell.tau_ps))) / DE_UNIT_EV})
        rows.append(r)
        step1.setdefault(cell.c, []).append(r)
    planned = {}
    for cell in plan:
        planned[cell.c] = planned.get(cell.c, 0) + 1
    for c, n in planned.items():
        print(f"Step 1 c={c}: {len(step1.get(c, []))}/{n} runs scored")

    # 4. fits
    with open(TRANSFER_CSV, newline="", encoding="utf-8") as fh:
        seed_row = next(r for r in csv.DictReader(fh) if r["cell"] == "seed_std_t64e310")
    seed_sd = {o: float(seed_row[col]) for o, col in SEED_SD_COLS.items()}
    with open(G413_TABLE_CSV, newline="", encoding="utf-8") as fh:
        stage_a = [r for r in csv.DictReader(fh)
                   if r["stage"] == "A" and r["label"].startswith("t")]
    fit_rows, curve_rows, cand_rows = [], [], []
    for c in sorted(planned):
        runs = step1.get(c, [])
        if len(runs) != planned[c]:
            continue
        tau = np.array([r["tau_ps"] for r in runs])
        de = np.array([r["dE"] for r in runs])
        order = (SHORT_LINE_ORDER if len(set(tau.tolist())) <= SHORT_LINE_MAX_TAUS
                 else FULL_LINE_ORDER)
        fits, statuses = {}, []
        e385 = e0_385_interp(line_rows, c)
        for obs in FIT_OBS + REPORT_OBS:
            y = np.array([float(r[obs]) for r in runs])
            de2 = obs == "n1_solv"
            if obs in FIT_OBS:
                res = fit_observable(tau, de, y, order, de2, seed_sd[obs])
                statuses.append(res["status"])
            else:
                f0 = fit_surface(tau, de, y, order, de2)
                res = {"fit": f0, "status": "reported", "threshold": float("nan"),
                       "p_tau": float("nan"), "order": order}
            fits[obs] = res["fit"]
            fr = {"c": c, "obs": obs, "n_runs": res["fit"].n, "n_params": res["fit"].p,
                  "tau_order": res["order"], "with_dE2": int(de2),
                  "resid_sd": res["fit"].resid_sd, "threshold": res["threshold"],
                  "status": res["status"], "p_tau": res["p_tau"],
                  **{f"beta{i}": b for i, b in enumerate(res["fit"].beta)}}
            if c == 2.0 and obs in FIT_OBS:
                sa_tau = np.array([float(r["tau_ps"]) for r in stage_a])
                sa_de = (np.array([float(r["E0_eV"]) for r in stage_a])
                         - e385(sa_tau)) / DE_UNIT_EV
                sa_res = np.array([float(r[obs]) for r in stage_a]) - predict(
                    res["fit"], sa_tau, sa_de)
                fr["oos_stageA_s1_mean_resid"] = float(np.mean(sa_res))
                fr["oos_stageA_s1_sd_resid"] = float(np.std(sa_res, ddof=1))
            fit_rows.append(fr)
        status = c_status(statuses)
        print(f"fit c={c}: order {order}, status {status} "
              + ", ".join(f"{o} {fits[o].resid_sd:.4g}" for o in FIT_OBS))
        if status != "resolved":
            continue
        curve = tradeoff_curve({o: fits[o] for o in FIT_OBS + REPORT_OBS},
                               (float(tau.min()), float(tau.max())),
                               (float(de.min()), float(de.max())))
        for r in curve:
            r.update({"c": c, "E0_eV": float(e385(r["tau_ps"])) + r["dE"] * DE_UNIT_EV})
        curve_rows.extend(curve)
        cands = select_candidates(curve)
        cand_rows.extend({**r, "c": c} for r in cands)
        print(f"  c={c}: {sum(r['feasible'] for r in curve)}/{len(curve)} feasible tau; "
              f"candidates {[(r['role'], r['tau_ps'], round(r['E0_eV'], 4)) for r in cands]}")

    _write(rows, out_dir / "g414_table.csv")
    if fit_rows:
        _write(fit_rows, out_dir / "g414_fit.csv")
    if curve_rows:
        _write(curve_rows, out_dir / "g414_tradeoff.csv")
        _write(cand_rows or [{"c": "", "role": "none_feasible"}],
               out_dir / "g414_candidates.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
