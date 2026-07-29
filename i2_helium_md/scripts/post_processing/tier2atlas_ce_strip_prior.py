"""Atlas (C) pre-step P1 — strip-form prior calibration (zero MD, detection-only).

Re-runs the §3.5j exit-stripping counterfactual instrument on the committed
``g4fh405`` checkpoints with the ACTUAL (C)-design strip form
(`TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md` §3.3)

    P_knock(j) = P0(v_x) * G(j)
    P0(v)      = min(1, (v / v_strip)^a)          [dimensionless]
    G(j)       = 1 / (1 + exp(-(j - j0) / w_j))   [dimensionless, j = 1 innermost]

instead of the §3.5j CF-3 sharp/eta variants — mapping (a, j0, w_j) ->
(n1 gain, KE1, bare weight, W1, midHot) at zero MD, and computing the
retention-twin anchor read (P(all stripped) vs P(one survivor) at the
n = 1-producing exit speed ~10.2 A/ps). Pins the prior box before probe
registration; family fixed (design OQ-D) — values may move inside the box only.

CF conventions (§3.5j, unchanged — ONLY the CF-3 knockout rule is replaced):

- CF-1: strip at the FIRST outbound surface crossing; the post-exit cascade is
  not re-run. Normally-detected ions: terminal n_cf = min(n_detected, n_strip).
  Suppressed ions (mechanical shell n = 21): n_cf = n_strip, re-scored as
  normally detected (the scaffolding retires in the counterfactual; survivor
  counts at high n are ceiling reads — the real cascade would shrink them).
- CF-2: the velocity history is unchanged; KE_cf = 1/2 * m(n_cf) * v_det^2.
  The strip toll Delta_E_j (design §3.3 energetics) is deliberately NOT
  applied — same limitation as §3.5j; the MD probe measures it.
- Independent Bernoulli draws per rung implement the retention-twin
  stochasticity (fixed named seed; replica mean +/- SD reported).

Oracles (§1.4 — hard-fail before any new number is read):

- O1 scorer drift: the pooled standing battery reproduces its recorded row.
- O2 committed-row: the finals h405 row rescored to 4 decimals.
- O3 kinematic conventions: stored E_kin_detected_eV reproduced bit-tight from
  (mass model, stored velocity); mass model = MASS_I_ION_AMU + n*MASS_HE_AMU.
- O4 crossing read: the suppressed-class crossing speeds reproduce the
  committed §3.5j p10-p90 band 10.66-11.69 A/ps.
- O5 §3.5j variant reproduction: the CF-3 sharp-knockout row and the
  supp -> n = 1 bracket row reproduce the committed §3.5j variants table
  (n_bar / n1_solv / KE1 / bare / W1) — this pins the crossing extraction
  and the CF-1/CF-2 conventions to the committed record.

Pure scorer — runs nothing, mutates nothing.

Invocation:

    python scripts/post_processing/tier2atlas_ce_strip_prior.py
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import json

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.physics.constants import (  # noqa: E402
    EV,
    MASS_HE_AMU,
    MASS_I_ION_AMU,
    U,
)
from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.size_distribution import N_STAR  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
)
from scripts.gen_tier2atlas_g4finals import finals_run_dir_name  # noqa: E402
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
    ORACLE_ROW_COLS,
    committed_h405_row,
    lowke_columns,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUNS_ROOT,
    check_oracle,
    observable_columns,
    observable_row,
)
from scripts.post_processing.tier2atlas_sn_table import sn_ke_columns  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS — the P1 grid (design §3.3 prior box + one margin ring)
# ---------------------------------------------------------------------------

RUN_DIR = RUNS_ROOT / finals_run_dir_name("h405")
SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_ce_strip_p1.csv"

V_STRIP_APS = 9.9        # Sourced (design §3.3): 1/2 m_He v^2 = D0(1) threshold
V_TWIN_APS = 10.2        # the n=1-producing exit speed (retention-twin read)

GRID_A = (1.0, 1.5, 2.0, 2.5, 3.0)          # prior a ~ 2 (ram scaling)
GRID_J0 = (1.0, 1.5, 2.0, 2.5, 3.0)         # prior [1, 3]
GRID_WJ = (0.5, 1.0, 1.5, 2.0)              # prior ~ 1

N_REPLICAS = 48
RNG_SEED = 20260729

# Committed §3.5j anchors (findings §3.5j — the O4/O5 oracle values).
O4_SUPP_CROSS_P10_P90 = (10.66, 11.69)      # A/ps
O4_TOL = 0.10
O5_SHARP = {"nbar_det": 2.901, "n1_solv": 0.445, "KE1_mean": 0.484,
            "bare": 0.329, "w1_solv": 0.787}
O5_BRACKET = {"nbar_det": 4.138, "n1_solv": 0.353, "KE1_mean": 0.708,
              "bare": 0.000, "w1_solv": 0.751}
O5_TOL = 0.005            # committed values quoted to 3 decimals

# Reference targets printed next to the grid (read-only context, not gates):
REF_N1_SOLV = 0.31        # reference solvated n = 1 share (CP-4 direction)
TWIN_SAME_ORDER = (0.2, 1.2)   # P(0)/P(1) "same-order" display band

# ---------------------------------------------------------------------------

_SUPPRESSED = "suppressed"
_DETECTED_STATE = "frozen"           # the state a de-suppressed ion re-scores as
_KE_COEF = 0.5 * U * 1e4 / EV        # eV per (amu * (A/ps)^2)


def _mass_amu(n: np.ndarray | float) -> np.ndarray | float:
    """Complex mass [amu] at shell count n (Tier-1a shell-schedule model)."""
    return MASS_I_ION_AMU + np.asarray(n, dtype=float) * MASS_HE_AMU


def _ke_eV(n: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """KE [eV] of a complex with shell n at squared speed v2 [(A/ps)^2]."""
    return _KE_COEF * np.asarray(_mass_amu(n)) * v2


def _stage_crossings(path: Path, need: np.ndarray, out_v: np.ndarray,
                     out_n: np.ndarray) -> np.ndarray:
    """Fill first-outbound-crossing (|v|, n_shell) from one trajectory stage.

    Crossing = first stored step with r >= R_droplet (per-fragment radius).
    Only fragments flagged in ``need`` (not yet crossed) are filled; returns
    the updated still-needed mask.
    """
    with np.load(path, allow_pickle=False) as d:
        R2 = d["droplet_radii_angstrom"][:, None] ** 2
        r2 = d["positions_x"] ** 2 + d["positions_y"] ** 2 + d["positions_z"] ** 2
        outside = r2 >= R2
        has = outside.any(axis=1) & need
        idx = outside.argmax(axis=1)
        sel = np.flatnonzero(has)
        vx = d["velocities_x"][sel, idx[sel]]
        vy = d["velocities_y"][sel, idx[sel]]
        vz = d["velocities_z"][sel, idx[sel]]
        out_v[sel] = np.sqrt(vx * vx + vy * vy + vz * vz)
        out_n[sel] = np.round(d["n_shell"][sel, idx[sel]])
    return need & ~has


def _survivor_pmf_head(surv_probs: np.ndarray, kmax: int = 2) -> np.ndarray:
    """P(exactly k survivors) for k = 0..kmax (Poisson-binomial head, exact DP)."""
    pmf = np.zeros(kmax + 1)
    pmf[0] = 1.0
    for s in surv_probs:
        carry = pmf[:-1] * s
        pmf *= 1.0 - s
        pmf[1:] += carry
    return pmf


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # ---- O1: scorer drift.
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** O1 SCORER-DRIFT ORACLE FAILED — do not read P1 ***")
        for line in drift:
            print(f"    {line}")
        return
    print("O1 scorer-drift oracle OK (pooled standing battery reproduces).")

    # ---- Load the committed run once; all variants score from these arrays.
    with np.load(RUN_DIR / "detection.npz", allow_pickle=False) as det:
        state = det["state_reason"].copy()
        n_det = det["n_detected"].astype(float).copy()
        ke_det = det["E_kin_detected_eV"].copy()
        m_det_amu = (det["mass_detected_kg"] / U).copy()
        v2_det = (det["vx_detected"] ** 2 + det["vy_detected"] ** 2
                  + det["vz_detected"] ** 2).copy()
    cfg = json.loads((RUN_DIR / "cfg.json").read_text(encoding="utf-8"))
    rungs_eV = np.asarray(cfg["tabulated_ladder_rungs_eV"], dtype=float)
    n_rungs = rungs_eV.size

    def score(n_arr, ke_arr, state_arr, label) -> dict[str, Any]:
        read = read_confirmation_detection(
            SimpleNamespace(state_reason=state_arr, n_detected=n_arr,
                            E_kin_detected_eV=ke_arr),
            label=label, n_max=N_STAR,
        )
        obs = observable_columns(read, abundance_ref, ked_ref)
        row = {"label": label,
               "nbar_det": obs["nbar_det"], "supp": obs["supp"],
               "n1_solv": obs["n1_solv"], "w1_solv": obs["w1_solv"],
               "midHot": obs["midHot"], "deepKE": obs["deepKE"],
               # bin0 = the whole n = 0 bin; under the committed gate the
               # suppressed class scores there too (supp column separates it).
               "bin0": float(read.fraction[0])}
        row.update(lowke_columns(read))
        row.update(sn_ke_columns(read))
        return row

    # ---- O2: committed finals h405 row through this code path.
    committed = committed_h405_row()
    baseline = score(n_det, ke_det, state, "gate(h405)")
    bad = [
        f"{c}: rescored {float(baseline[c]):.4f} vs committed "
        f"{float(committed[c]):.4f}"
        for c in ORACLE_ROW_COLS
        if f"{float(baseline[c]):.4f}" != f"{float(committed[c]):.4f}"
    ]
    if bad:
        print("*** O2 COMMITTED-ROW ORACLE FAILED — do not read P1 ***")
        for line in bad:
            print(f"    {line}")
        return
    print("O2 committed-row oracle OK (finals h405 rescored to 4 decimals).")

    # ---- O3: kinematic conventions.
    ke_check = _KE_COEF * m_det_amu * v2_det
    ke_err = float(np.abs(ke_check - ke_det).max())
    scored_mask = ~np.isin(state, ("droplet_retained",
                                   "droplet_retained_marginal"))
    normal_mask = scored_mask & (state != _SUPPRESSED)
    mass_err = float(np.abs(m_det_amu[normal_mask]
                            - _mass_amu(n_det[normal_mask])).max())
    sup_mask = state == _SUPPRESSED
    sup_mass_err = float(np.abs(m_det_amu[sup_mask] - _mass_amu(21.0)).max())
    if ke_err > 1e-9 or mass_err > 1e-9 or sup_mass_err > 1e-9:
        print(f"*** O3 KINEMATICS ORACLE FAILED: ke_err {ke_err:.3e}, "
              f"mass_err {mass_err:.3e}, sup_mass_err {sup_mass_err:.3e} ***")
        return
    print("O3 kinematics oracle OK (stored KE = 1/2 m(n) v^2 bit-tight; "
          "suppressed carry the full n = 21 mechanical mass).")

    # ---- Crossing extraction (ion window first, then relaxation).
    num = state.size
    v_x = np.full(num, np.nan)
    n_x = np.full(num, np.nan)
    need = np.ones(num, dtype=bool)
    need = _stage_crossings(RUN_DIR / "ion.npz", need, v_x, n_x)
    need = _stage_crossings(RUN_DIR / "relaxation.npz", need, v_x, n_x)
    crossed = ~need
    n_x = np.clip(n_x, 0, n_rungs)

    det_uncrossed = int(np.count_nonzero(scored_mask & ~crossed))
    sup_v = v_x[sup_mask]
    p10, p90 = np.percentile(sup_v, [10, 90])
    n1_now = scored_mask & (n_det == 1) & (state != _SUPPRESSED)
    print(f"Crossings: {int(crossed.sum())}/{num} fragments crossed "
          f"(scored-but-uncrossed: {det_uncrossed}, strip not applied to "
          f"these); suppressed v_x p10-p90 {p10:.2f}-{p90:.2f} A/ps; "
          f"current n = 1 enders mean v_x "
          f"{float(v_x[n1_now].mean()):.2f} A/ps.")

    # ---- O4: the committed suppressed crossing-speed band.
    if (abs(p10 - O4_SUPP_CROSS_P10_P90[0]) > O4_TOL
            or abs(p90 - O4_SUPP_CROSS_P10_P90[1]) > O4_TOL):
        print(f"*** O4 CROSSING ORACLE FAILED: suppressed p10-p90 "
              f"{p10:.2f}-{p90:.2f} vs committed "
              f"{O4_SUPP_CROSS_P10_P90} +/- {O4_TOL} ***")
        return
    print("O4 crossing oracle OK (suppressed band matches §3.5j).")

    # ---- Counterfactual builders (CF-1/CF-2 shared by O5 and the grid).
    sup_crossed = sup_mask & crossed
    norm_crossed = normal_mask & crossed

    def apply_strip(n_strip_full: np.ndarray, label: str) -> dict[str, Any]:
        """CF-1/CF-2: build the counterfactual arrays and score them.

        ``n_strip_full`` is the per-fragment survivor count (only entries at
        crossed fragments are used).
        """
        n_cf = n_det.copy()
        state_cf = state.copy()
        n_cf[sup_crossed] = n_strip_full[sup_crossed]
        state_cf[sup_crossed] = _DETECTED_STATE
        n_cf[norm_crossed] = np.minimum(n_det[norm_crossed],
                                        n_strip_full[norm_crossed])
        ke_cf = ke_det.copy()
        # Every re-binned fragment gets the CF-2 mass-consistent KE; a
        # suppressed fragment was scored with FULL mechanical mass, so its
        # KE is recomputed even when n_cf happens to stay 21.
        changed = (n_cf != n_det) | sup_crossed
        ke_cf[changed] = _ke_eV(n_cf[changed], v2_det[changed])
        row = score(n_cf, ke_cf, state_cf, label)
        conv1 = sup_crossed & (n_cf == 1)
        row["KE1_conv"] = (float(ke_cf[conv1].mean())
                           if np.count_nonzero(conv1) >= 5 else float("nan"))
        row["conv_n0_frac"] = float(np.count_nonzero(sup_crossed & (n_cf == 0))
                                    / max(1, np.count_nonzero(sup_crossed)))
        return row

    # ---- O5a: the §3.5j CF-3 sharp-knockout variant.
    v_thr = np.sqrt(rungs_eV / (_KE_COEF * MASS_HE_AMU))   # per-rung [A/ps]
    j_idx = np.arange(1, n_rungs + 1)
    valid = j_idx[None, :] <= n_x[:, None]                  # rung present
    knocked_sharp = valid & (v_x[:, None] >= v_thr[None, :])
    n_sharp = (valid & ~knocked_sharp).sum(axis=1).astype(float)
    sharp = apply_strip(n_sharp, "sharp(CF-3)")
    bad5 = [
        f"{k}: {float(sharp[kk]):.3f} vs committed {ref:.3f}"
        for k, kk, ref in (
            ("nbar", "nbar_det", O5_SHARP["nbar_det"]),
            ("n1", "n1_solv", O5_SHARP["n1_solv"]),
            ("KE1", "KE1_mean", O5_SHARP["KE1_mean"]),
            ("bare", "bin0", O5_SHARP["bare"]),
            ("W1", "w1_solv", O5_SHARP["w1_solv"]),
        )
        if abs(float(sharp[kk]) - ref) > O5_TOL
    ]
    # ---- O5b: the §3.5j bracket (every suppressed -> n = 1 at existing v).
    bracket = apply_strip(np.where(sup_crossed, 1.0, n_det), "bracket(supp->1)")
    bad5 += [
        f"bracket {k}: {float(bracket[kk]):.3f} vs committed {ref:.3f}"
        for k, kk, ref in (
            ("nbar", "nbar_det", O5_BRACKET["nbar_det"]),
            ("n1", "n1_solv", O5_BRACKET["n1_solv"]),
            ("KE1", "KE1_mean", O5_BRACKET["KE1_mean"]),
            ("bare", "bin0", O5_BRACKET["bare"]),
            ("W1", "w1_solv", O5_BRACKET["w1_solv"]),
        )
        if abs(float(bracket[kk]) - ref) > O5_TOL
    ]
    if bad5:
        print("*** O5 §3.5j VARIANT ORACLE FAILED — do not read the grid ***")
        for line in bad5:
            print(f"    {line}")
        return
    print("O5 variant oracle OK (sharp CF-3 + supp->n=1 bracket reproduce "
          "the committed §3.5j table).\n")

    # ---- The P1 grid (new numbers — read only after O1..O5).
    print(f"=== P1 strip-form prior calibration on {RUN_DIR.name} ===")
    print(f"P_knock(j) = min(1, (v_x/{V_STRIP_APS})^a) / (1 + exp(-(j - j0)/w_j))"
          f"; {N_REPLICAS} replicas, seed {RNG_SEED}; CF-1/CF-2 conventions; "
          "toll NOT applied (§3.5j limitation carried).")
    print(f"baseline gate row: nbar {float(baseline['nbar_det']):.3f}, "
          f"n1 {float(baseline['n1_solv']):.3f}, "
          f"KE1 {float(baseline['KE1_mean']):.3f}, "
          f"W1 {float(baseline['w1_solv']):.3f}\n")

    rng = np.random.default_rng(RNG_SEED)
    rows: list[dict[str, Any]] = []
    rep_keys = ("nbar_det", "n1_solv", "bin0", "w1_solv", "midHot", "deepKE",
                "KE1_mean", "KE1_sd", "KE2_mean", "KE1_conv", "conv_n0_frac")
    cross_idx = np.flatnonzero(crossed)
    v_c = v_x[cross_idx]
    valid_c = valid[cross_idx]

    for a in GRID_A:
        p0_c = np.minimum(1.0, (v_c / V_STRIP_APS) ** a)
        p0_twin = min(1.0, (V_TWIN_APS / V_STRIP_APS) ** a)
        for j0 in GRID_J0:
            for wj in GRID_WJ:
                g = 1.0 / (1.0 + np.exp(-(j_idx - j0) / wj))
                surv_c = np.where(valid_c,
                                  1.0 - p0_c[:, None] * g[None, :], 0.0)
                # retention-twin read: 21-He shell at the n=1-producing speed
                pmf = _survivor_pmf_head(1.0 - p0_twin * g[:21], kmax=2)
                reps: dict[str, list[float]] = {k: [] for k in rep_keys}
                for _ in range(N_REPLICAS):
                    u = rng.random(surv_c.shape)
                    n_strip_c = (u < surv_c).sum(axis=1).astype(float)
                    n_strip_full = np.zeros(num)
                    n_strip_full[cross_idx] = n_strip_c
                    r = apply_strip(n_strip_full,
                                    f"a{a:g}_j{j0:g}_w{wj:g}")
                    for k in rep_keys:
                        reps[k].append(float(r[k]))
                row: dict[str, Any] = {
                    "a": a, "j0": j0, "w_j": wj,
                    "P0_twin": p0_twin,
                    "P_strip0": float(pmf[0]), "P_surv1": float(pmf[1]),
                    "P_surv2": float(pmf[2]),
                    "twin_ratio": float(pmf[0] / pmf[1]) if pmf[1] > 0
                    else float("inf"),
                }
                for k in rep_keys:
                    arr = np.asarray(reps[k], dtype=float)
                    ok = np.isfinite(arr)
                    row[k] = float(arr[ok].mean()) if ok.any() else float("nan")
                    row[f"{k}_repsd"] = (float(arr[ok].std(ddof=1))
                                         if ok.sum() > 1 else float("nan"))
                rows.append(row)

    # ---- Console: the w_j = 1.0 slice (full grid in the CSV).
    hdr = (f"{'a':>4} {'j0':>4} {'w_j':>4} | {'P(0)':>6} {'P(1)':>6} "
           f"{'P0/P1':>6} | {'n1':>6} {'KE1':>6} {'KE1sd':>6} {'bin0':>6} "
           f"{'nbar':>6} {'W1':>6} {'midHot':>7} {'KE1cnv':>7}")
    for wj_show in (1.0,):
        print(f"--- w_j = {wj_show} slice ---")
        print(hdr)
        for row in rows:
            if row["w_j"] != wj_show:
                continue
            print(f"{row['a']:>4g} {row['j0']:>4g} {row['w_j']:>4g} | "
                  f"{row['P_strip0']:>6.3f} {row['P_surv1']:>6.3f} "
                  f"{row['twin_ratio']:>6.2f} | {row['n1_solv']:>6.3f} "
                  f"{row['KE1_mean']:>6.3f} {row['KE1_sd']:>6.3f} "
                  f"{row['bin0']:>6.3f} {row['nbar_det']:>6.2f} "
                  f"{row['w1_solv']:>6.3f} {row['midHot']:>7.3f} "
                  f"{row['KE1_conv']:>7.3f}")
        print()

    # ---- The prior-box read (display summary; adjudication in the findings).
    inband = [r for r in rows
              if TWIN_SAME_ORDER[0] <= r["twin_ratio"] <= TWIN_SAME_ORDER[1]]
    print(f"Retention-twin same-order band P(0)/P(1) in {TWIN_SAME_ORDER}: "
          f"{len(inband)}/{len(rows)} cells.")
    if inband:
        best = sorted(inband, key=lambda r: abs(r["n1_solv"] - REF_N1_SOLV))
        print(f"closest n1 to the reference {REF_N1_SOLV} inside that band:")
        for r in best[:8]:
            print(f"    a {r['a']:g} j0 {r['j0']:g} w_j {r['w_j']:g}: "
                  f"n1 {r['n1_solv']:.3f}, KE1 {r['KE1_mean']:.3f} "
                  f"(sd {r['KE1_sd']:.3f}), bin0 {r['bin0']:.3f}, "
                  f"W1 {r['w1_solv']:.3f}, P0/P1 {r['twin_ratio']:.2f}")

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
