"""Atlas §3.5g — the low-n KE retro-scan over the committed surface (zero MD).

Retro-scores **KE₁/KE₂** (detected mean KE in bins n = 1 / n = 2) over every
committed corrected-geometry run plus the standing-geometry references — 52
rows, no new MD — to measure the in-surface (v_c, τ, E₀) freedom on the low-n
KE axis *before* any new axis (p_tail) is opened. Pre-registration:
`TIER2_SENSITIVITY_ATLAS_PLAN.md` §3.5g (committed before this scorer ran).

Convention (§3.5g amendment, user adjudication 2026-07-29): KE₁ is scored
against the experimental n = 1 **peak ≈ 1.00 eV** (the wide upper tail is
read as a second process and not chased); the old median anchor 1.128 stays
reported. KE₂ keeps the reference mean 0.706 eV.

Pure scorer — runs nothing, mutates nothing. Every physics convention is
imported from the committed modules (rule 1): reads via
`tier2_confirmation.load_confirmation_run` / `pool_confirmation_reads`,
observables via `tier2atlas_geometry_table.observable_columns`, KE columns
via `tier2atlas_g4step2_battery_table.lowke_columns`.

Order of operations (§1.4 — oracles first, hard-fail before any new number):

1. **O1 scorer-drift** — the standing pooled battery reproduces its recorded
   observable row (`check_oracle`).
2. **O2 committed-KE** — the five h405 battery members + pooled, rescored
   through this script's own code path, reproduce the committed
   `atlas_g4step2_battery.csv` KE₁/KE₂ columns to 1e-9 relative.
3. **O3 incumbent-KE** — the incumbent pooled battery lands on the recorded
   KE₁ 1.034 / KE₂ 0.754 (3-decimal comparison).

Then the pre-registered readings R1 (τ), R2 (v_c), R3 (E₀), R4 (geometry at
fixed chord) and R5 (the in-gate reachability verdict vs KE₁ ≥ 0.95).

Invocation:

    python scripts/post_processing/tier2atlas_ke_lown_scan.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any, Sequence

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    load_confirmation_run,
    pool_confirmation_reads,
)
from scripts.gen_tier2atlas_g3ring import (  # noqa: E402
    RING_MATRIX,
    ring_run_dir_name,
)
from scripts.gen_tier2atlas_g4finals import (  # noqa: E402
    FINAL_MATRIX,
    STANDING_E0_EV,
    STANDING_TAU_PS,
    STANDING_V_C,
    finals_run_dir_name,
)
from scripts.gen_tier2atlas_g4step2_battery import (  # noqa: E402
    MEMBER_SEEDS,
    member_run_dir_name,
)
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
    lowke_columns,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    CELL_LABELS,
    IHE_KED_REFERENCE_CSV,
    ORACLE_RUN,
    RUN_SELECTION,
    RUNS_ROOT,
    check_oracle,
    geometry_columns,
    observable_columns,
    observable_row,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_ke_lown_scan.csv"

# §3.5g scoring anchors (user adjudication 2026-07-29).
KE1_ANCHOR_EV = 1.00        # experimental n = 1 peak (supersedes the median
KE1_MEDIAN_OLD_EV = 1.128   # anchor for scoring; median stays reported)
KE2_REF_EV = 0.706          # reference n = 2 mean (unchanged)

# O2: the committed Block-V battery CSV, compared on the KE columns.
COMMITTED_BATTERY_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_g4step2_battery.csv"
KE_ORACLE_COLS = ("KE1_mean", "KE1_med", "KE1_mode",
                  "KE2_mean", "KE2_med", "KE2_mode")
KE_ORACLE_RTOL = 1e-9

# O3: the incumbent pooled battery's recorded low-n KE (Block-V constants).
INCUMBENT_KE_RECORDED = {"KE1_mean": 1.034, "KE2_mean": 0.754}

# The incumbent battery run dirs (standing geometry, standing chord).
INCUMBENT_MEMBER_RUNS = {
    f"s{i}": f"9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s{i}"
    for i in range(1, 6)
}
INCUMBENT_POOLED_RUN = ORACLE_RUN

# MD gate bands (§3.5d/e conventions, incumbent-anchored — see the §3.5f
# provenance correction) + the soft midHot band used by R5.
GATE_N1 = (0.19, 0.30)
GATE_NBAR = (3.77, 4.37)
MIDHOT_BAND = (0.85, 1.15)

# R5 target: KE₁ ≥ 1.00 − 1σ of the linearization (plan §3.5g).
R5_TARGET_EV = 0.95

# Pre-registered slope cell sets (finals labels; plan §3.5g R1–R3).
E0_ARM_T44 = ("f3", "h405", "h410", "h415")          # τ 4.4, v_c 5.5
E0_ARM_T48 = ("x345", "f1", "a037", "h375", "h380")  # τ 4.8, v_c 5.5
E0_ARM_T52 = ("f2", "h345", "h350", "h355")          # τ 5.2, v_c 5.5
TAU_GATED_CHORD = ("h405", "h375", "h345")           # gated ridge chord
VC_PAIR = ("v525", "f3")                             # 5.25 vs 5.5, τ 4.4 arm

# ---------------------------------------------------------------------------
# Pure helpers (unit-tested)
# ---------------------------------------------------------------------------


def fit_slope(points: Sequence[tuple[float, float]]) -> float:
    """Least-squares slope dy/dx through ``points`` [(x, y), ...].

    Requires >= 2 points with non-identical x; raises ValueError otherwise.
    """
    pts = [(float(x), float(y)) for x, y in points]
    if len(pts) < 2:
        raise ValueError(f"fit_slope needs >= 2 points, got {len(pts)}")
    xs = np.asarray([p[0] for p in pts])
    ys = np.asarray([p[1] for p in pts])
    if np.allclose(xs, xs[0]):
        raise ValueError("fit_slope: all x identical")
    return float(np.polyfit(xs, ys, 1)[0])


def max_ingate_move(value_by_obs: dict[str, float],
                    slope_by_obs: dict[str, float],
                    band_by_obs: dict[str, tuple[float, float]],
                    direction: float,
                    eps: float = 1e-12) -> float:
    """Largest |Δknob| along ``direction`` (±1) before any gated observable
    leaves its band, from current values and measured d(obs)/d(knob) slopes.

    Observables with |slope| < eps do not constrain the move. Returns +inf
    when nothing constrains it. Raises ValueError if a current value is
    already outside its band (headroom would be negative).
    """
    if direction not in (+1.0, -1.0):
        raise ValueError(f"direction must be ±1, got {direction}")
    limit = float("inf")
    for obs, slope in slope_by_obs.items():
        lo, hi = band_by_obs[obs]
        val = value_by_obs[obs]
        if not lo <= val <= hi:
            raise ValueError(f"{obs}={val} already outside band [{lo}, {hi}]")
        step = slope * direction
        if abs(step) < eps:
            continue
        headroom = (hi - val) if step > 0 else (val - lo)
        limit = min(limit, headroom / abs(step))
    return limit


def lowke_counts(read) -> dict[str, int]:
    """Per-bin populations behind KE₁/KE₂ (single-seed SD context)."""
    n = np.asarray(read.n_scored)
    return {"KE1_n": int((n == 1).sum()), "KE2_n": int((n == 2).sum())}


def ke_oracle_mismatches(rows_by_label: dict[str, dict[str, Any]],
                         committed_csv: Path,
                         cols: Sequence[str] = KE_ORACLE_COLS,
                         rtol: float = KE_ORACLE_RTOL) -> list[str]:
    """O2: compare rescored KE columns against the committed battery CSV.

    Returns human-readable mismatch lines (empty == oracle passes). Every
    committed label must be present in ``rows_by_label``.
    """
    problems: list[str] = []
    with open(committed_csv, newline="", encoding="utf-8") as fh:
        for committed in csv.DictReader(fh):
            label = committed["label"]
            if label not in rows_by_label:
                problems.append(f"{label}: missing from rescored rows")
                continue
            row = rows_by_label[label]
            for col in cols:
                want = float(committed[col])
                got = float(row[col])
                if not np.isclose(got, want, rtol=rtol, atol=0.0):
                    problems.append(
                        f"{label}.{col}: rescored {got!r} vs committed {want!r}"
                    )
    return problems


# ---------------------------------------------------------------------------


def _blank_geometry() -> dict[str, Any]:
    return {"R_A": "-", "depth_mean_A": "-"}


def scan_row(group: str, label: str, run_dir: Path,
             v_c: float, well: str, tau: float, e0: float,
             abundance_ref, ked_ref, *,
             read=None, geometry: dict[str, Any] | None = None) -> dict[str, Any]:
    """One uniform scan row (shared key set across all groups)."""
    if read is None:
        read = load_confirmation_run(run_dir, label=label)
    obs = observable_columns(read, abundance_ref, ked_ref)
    ke = lowke_columns(read)
    row: dict[str, Any] = {
        "group": group, "label": label,
        "v_c": v_c, "well": well, "tau": tau, "E0": e0,
    }
    row.update(geometry if geometry is not None else _blank_geometry())
    row.update({
        "num_scored": obs["num_scored"],
        "nbar_det": obs["nbar_det"], "n1_solv": obs["n1_solv"],
        "w1_solv": obs["w1_solv"], "midHot": obs["midHot"],
        "deepKE": obs["deepKE"], "supp": obs["supp"], "trap": obs["trap"],
    })
    row.update(ke)
    row.update(lowke_counts(read))
    ke1 = float(ke["KE1_mean"])
    ke2 = float(ke["KE2_mean"])
    row.update({
        "KE1_r100": ke1 / KE1_ANCHOR_EV,
        "KE1_rmed": ke1 / KE1_MEDIAN_OLD_EV,
        "KE2_r": ke2 / KE2_REF_EV,
        "md_gate": int(GATE_N1[0] <= float(obs["n1_solv"]) <= GATE_N1[1]
                       and GATE_NBAR[0] <= float(obs["nbar_det"]) <= GATE_NBAR[1]),
    })
    return row


def _print_group(title: str, rows: list[dict[str, Any]]) -> None:
    print(f"\n=== {title} ===")
    print(format_table(rows))


def _row(rows: list[dict[str, Any]], label: str) -> dict[str, Any]:
    for r in rows:
        if r["label"] == label:
            return r
    raise KeyError(f"no row {label!r}")


def _slopes_vs(rows: list[dict[str, Any]], labels: Sequence[str],
               knob_key: str, obs_keys: Sequence[str]) -> dict[str, float]:
    picked = [_row(rows, lb) for lb in labels]
    return {
        obs: fit_slope([(float(r[knob_key]), float(r[obs])) for r in picked])
        for obs in obs_keys
    }


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # ---- O1: scorer drift (before any new number).
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** O1 SCORER-DRIFT ORACLE FAILED — do not read the scan ***")
        for line in drift:
            print(f"    {line}")
        return
    print("O1 scorer-drift oracle OK (pooled standing battery reproduces).")

    # ---- h405 battery rows (needed by O2 before anything else is read).
    battery_rows: list[dict[str, Any]] = []
    battery_reads = []
    for member in MEMBER_SEEDS:
        run_dir = RUNS_ROOT / member_run_dir_name(member)
        read = load_confirmation_run(run_dir, label=member)
        battery_reads.append(read)
        battery_rows.append(scan_row(
            "h405bat", member, run_dir, 5.5, "eb1168", 4.4, 0.405,
            abundance_ref, ked_ref, read=read,
        ))
    pooled_read = pool_confirmation_reads(battery_reads, label="h405pooled")
    battery_rows.append(scan_row(
        "h405bat", "pooled", RUNS_ROOT, 5.5, "eb1168", 4.4, 0.405,
        abundance_ref, ked_ref, read=pooled_read,
    ))
    mismatches = ke_oracle_mismatches(
        {r["label"]: r for r in battery_rows}, COMMITTED_BATTERY_CSV,
    )
    if mismatches:
        print("*** O2 COMMITTED-KE ORACLE FAILED — do not read the scan ***")
        for line in mismatches:
            print(f"    {line}")
        return
    print("O2 committed-KE oracle OK (battery KE columns reproduce to 1e-9).")

    # ---- O3: incumbent pooled KE vs the recorded Block-V constants.
    incumbent_rows: list[dict[str, Any]] = []
    for member, run_name in INCUMBENT_MEMBER_RUNS.items():
        incumbent_rows.append(scan_row(
            "incumbent", member, RUNS_ROOT / run_name,
            STANDING_V_C, "eb1168", STANDING_TAU_PS, STANDING_E0_EV,
            abundance_ref, ked_ref,
        ))
    inc_pooled = scan_row(
        "incumbent", "pooled", RUNS_ROOT / INCUMBENT_POOLED_RUN,
        STANDING_V_C, "eb1168", STANDING_TAU_PS, STANDING_E0_EV,
        abundance_ref, ked_ref,
    )
    incumbent_rows.append(inc_pooled)
    o3_bad = [
        f"{col}: rescored {float(inc_pooled[col]):.3f} vs recorded {want:.3f}"
        for col, want in INCUMBENT_KE_RECORDED.items()
        if f"{float(inc_pooled[col]):.3f}" != f"{want:.3f}"
    ]
    if o3_bad:
        print("*** O3 INCUMBENT-KE ORACLE FAILED — do not read the scan ***")
        for line in o3_bad:
            print(f"    {line}")
        return
    print("O3 incumbent-KE oracle OK (pooled 1.034 / 0.754 reproduce).\n")

    # ---- The remaining groups (new numbers — read only after O1–O3).
    grid_rows = []
    for cell in CELL_LABELS:
        run_dir = RUNS_ROOT / RUN_SELECTION[cell]
        geo = geometry_columns(run_dir)
        grid_rows.append(scan_row(
            "geogrid", cell, run_dir,
            STANDING_V_C, "eb1168", STANDING_TAU_PS, STANDING_E0_EV,
            abundance_ref, ked_ref,
            geometry={"R_A": geo["R_A"], "depth_mean_A": geo["depth_mean_A"]},
        ))

    ring_rows = [
        scan_row("g3ring", c.label, RUNS_ROOT / ring_run_dir_name(c.label),
                 c.v_c, c.eb_tag, c.tau_ps, c.e0_eV, abundance_ref, ked_ref)
        for c in RING_MATRIX
    ]
    finals_rows = [
        scan_row("g4finals", c.label, RUNS_ROOT / finals_run_dir_name(c.label),
                 c.v_c, c.eb_tag, c.tau_ps, c.e0_eV, abundance_ref, ked_ref)
        for c in FINAL_MATRIX
    ]

    _print_group("Incumbent battery (standing geometry, standing chord)",
                 incumbent_rows)
    _print_group("Geometry grid (fixed geometry, standing chord)", grid_rows)
    _print_group("G3 ring (corrected ensemble, N = 500)", ring_rows)
    _print_group("G4 finals + ladder arm (corrected ensemble, N = 1000)",
                 finals_rows)
    _print_group("h405 battery (corrected ensemble, 5 × N = 1000)",
                 battery_rows)

    # ---- Pre-registered readings.
    pooled_row = _row(battery_rows, "pooled")
    obs_keys = ("KE1_mean", "KE2_mean", "n1_solv", "nbar_det", "midHot")

    print("\n=== R1 (τ): v_c 5.5 ===")
    tau_chord = _slopes_vs(finals_rows, TAU_GATED_CHORD, "tau", obs_keys)
    print(f"gated ridge chord {TAU_GATED_CHORD} (E₀ co-moves to hold the "
          f"gate): dKE1/dτ {tau_chord['KE1_mean']:+.4f} eV per ps, "
          f"dKE2/dτ {tau_chord['KE2_mean']:+.4f}")
    for name, arm in (("τ 4.4", E0_ARM_T44), ("τ 4.8", E0_ARM_T48),
                      ("τ 5.2", E0_ARM_T52)):
        vals = [float(_row(finals_rows, lb)["KE1_mean"]) for lb in arm]
        print(f"  arm {name}: KE1 {', '.join(f'{v:.3f}' for v in vals)}")

    print("\n=== R2 (v_c) ===")
    for label, rows_src in (("c50", ring_rows), ("v525", finals_rows),
                            ("f3", finals_rows), ("h405", finals_rows),
                            ("b031", finals_rows), ("c65", ring_rows),
                            ("f725", ring_rows)):
        r = _row(rows_src, label)
        print(f"  v_c {float(r['v_c']):.2f}  τ {float(r['tau']):.1f}  "
              f"E₀ {float(r['E0']):.3f}  [{r['group']}]  "
              f"KE1 {float(r['KE1_mean']):.3f}  KE2 {float(r['KE2_mean']):.3f}  "
              f"midHot {float(r['midHot']):.3f}  gate {r['md_gate']}")
    vc_slopes = _slopes_vs(finals_rows, VC_PAIR, "v_c", obs_keys)
    print(f"pair {VC_PAIR} (τ 4.4 arm, ΔE₀ 0.005): dKE1/dv_c "
          f"{vc_slopes['KE1_mean']:+.4f} eV per Å/ps, dKE2/dv_c "
          f"{vc_slopes['KE2_mean']:+.4f}, dmidHot/dv_c "
          f"{vc_slopes['midHot']:+.4f}")

    print("\n=== R3 (E₀) ===")
    e0_slopes_by_arm = {}
    for name, arm in (("t44", E0_ARM_T44), ("t48", E0_ARM_T48),
                      ("t52", E0_ARM_T52)):
        s = _slopes_vs(finals_rows, arm, "E0", obs_keys)
        e0_slopes_by_arm[name] = s
        print(f"  arm {name}: dKE1/dE₀ {s['KE1_mean']:+.3f} eV/eV, "
              f"dKE2/dE₀ {s['KE2_mean']:+.3f}, dn1/dE₀ {s['n1_solv']:+.3f}, "
              f"dn̄/dE₀ {s['nbar_det']:+.3f}, dmidHot/dE₀ {s['midHot']:+.3f}")

    print("\n=== R4 (geometry at fixed standing chord) ===")
    geo_sorted = sorted(grid_rows, key=lambda r: (float(r["R_A"]),
                                                  float(r["depth_mean_A"])))
    for r in geo_sorted:
        print(f"  R {float(r['R_A']):5.1f} Å  depth {float(r['depth_mean_A']):5.1f} Å  "
              f"KE1 {float(r['KE1_mean']):.3f} (n={r['KE1_n']})  "
              f"KE2 {float(r['KE2_mean']):.3f} (n={r['KE2_n']})  "
              f"n1_solv {float(r['n1_solv']):.3f}")
    print(f"  incumbent pooled (sampled ⟨R⟩ ≈ 26.6, depth ≈ 9): "
          f"KE1 {float(inc_pooled['KE1_mean']):.3f}  "
          f"KE2 {float(inc_pooled['KE2_mean']):.3f}")
    f725 = _row(ring_rows, "f725")
    print(f"  f725 (standing chord, corrected ensemble): "
          f"KE1 {float(f725['KE1_mean']):.3f}  KE2 {float(f725['KE2_mean']):.3f}")

    print("\n=== R5 (in-gate reachability vs KE₁ ≥ "
          f"{R5_TARGET_EV:.2f} eV) ===")
    current = {
        "n1_solv": float(pooled_row["n1_solv"]),
        "nbar_det": float(pooled_row["nbar_det"]),
        "midHot": float(pooled_row["midHot"]),
    }
    bands = {"n1_solv": GATE_N1, "nbar_det": GATE_NBAR, "midHot": MIDHOT_BAND}
    ke1_now = float(pooled_row["KE1_mean"])
    need = R5_TARGET_EV - ke1_now
    print(f"h405 pooled KE1 {ke1_now:.3f}; required ΔKE1 ≥ {need:+.3f} eV")
    total_best = 0.0
    for knob, slopes in (
        ("E0 [eV] (τ4.4 arm)", e0_slopes_by_arm["t44"]),
        ("tau [ps] (gated chord)", tau_chord),
        ("v_c [Å/ps] (v525-f3 pair)", vc_slopes),
    ):
        s_ke1 = slopes["KE1_mean"]
        gate_slopes = {k: slopes[k] for k in ("n1_solv", "nbar_det", "midHot")}
        direction = +1.0 if s_ke1 >= 0 else -1.0
        move = max_ingate_move(current, gate_slopes, bands, direction)
        dke1 = abs(s_ke1) * move if np.isfinite(move) else float("inf")
        total_best += 0.0 if not np.isfinite(dke1) else dke1
        move_txt = "unbounded" if not np.isfinite(move) else f"{move:.4f}"
        print(f"  {knob}: dKE1/dknob {s_ke1:+.3f}; max in-gate move "
              f"{move_txt} (direction {direction:+.0f}) -> ΔKE1 ≤ "
              f"{dke1 if np.isfinite(dke1) else float('nan'):+.4f} eV")
    print(f"  additive single-knob bound (generous — ignores joint gate "
          f"consumption): ΔKE1 ≤ {total_best:+.4f} eV")
    verdict = total_best >= need
    print(f"R5 VERDICT: in-surface freedom "
          f"{'MAY REACH' if verdict else 'CANNOT REACH'} KE1 "
          f"{R5_TARGET_EV:.2f} eV "
          f"({'identify the chord and register an MD arm' if verdict else 'the in-surface freedom is exhausted on the low-n KE axis; the p_tail axis is motivated with this record as its evidence'}).")

    all_rows = (incumbent_rows + grid_rows + ring_rows + finals_rows
                + battery_rows)
    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, all_rows)}")


if __name__ == "__main__":
    main()
