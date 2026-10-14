"""Shed timing of the B2 s_eff scan vs Albrechtsen et al. (Nature 623, 319 (2023)).

Reference statement (Methods, "Molecular dynamics simulations of the
dissociation time of Na+HeN complexes"; Extended Data Fig. 4): classical MD
of isolated Na+He6 / Na+He14 shows dissociation "almost exclusively during
the first tens of picoseconds"; "the vast majority of dissociation events are
finished after 400 ps". Transfer caveats: Na+ not I+, classical (all modes),
isolated complexes heated to 20-100 K.

Per s_eff cell (``gen_tier2_b2_seff_scan.py --keep-ion`` dirs, seed s1,
N = 1000): every He shed event, timed

* in the MD window (0-500 ps): from the stored ``n_shell`` history, one shed
  per unit decrease between stored columns (stride 16 x 0.01 ps = 0.16 ps;
  a pickup and a shed inside the same 0.16 ps interval cancel -- a small
  undercount of both, disclosed);
* in the detection stage (500 ps - 8.53 us): one shed per recorded event.

Columns: total sheds and pickups, fraction of all shed events after 50 /
100 / 400 / 500 ps, and the fraction of ions shedding at all after 400 ps.

Output (commit with ``git add -f``): ``data/runs/h2b_forward_model/
b2_seff_timing.csv``.

Invocation::

    python scripts/post_processing/tier2_b2_seff_timing.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2_b2_seff_scan import SEFF_GRID, SEFF_PRODUCTION, seff_run_dir  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table, write_rows_csv,
)
from scripts.post_processing.tier2atlas_geometry_table import RUNS_ROOT  # noqa: E402

OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "b2_seff_timing.csv"
#: Thresholds [ps] for the late-shed fractions (400 = the paper's bound).
LATE_PS = (50.0, 100.0, 400.0, 500.0)


def shed_times(run_dir: Path) -> tuple[np.ndarray, np.ndarray, int]:
    """All shed-event times [ps] and the owning ion index, plus the pickup
    count, for one run (window from n_shell history + detection events)."""
    rd = RunDirectory(run_dir)
    ion = rd.load_ion()
    n = np.asarray(ion.n_shell, dtype=float)
    t = np.asarray(ion.time_ps, dtype=float)
    if n.ndim != 2 or n.shape[1] != t.size:
        raise AssertionError(f"{run_dir.name}: n_shell shape {n.shape} vs time {t.size}")
    dn = np.diff(n, axis=1)                                   # (2N, T-1)
    k = np.clip(-dn, 0, None).astype(int)
    ion_idx, col = np.nonzero(k)
    reps = k[ion_idx, col]
    t_win = np.repeat(t[1:][col], reps)
    i_win = np.repeat(ion_idx, reps)
    pickups = int(np.clip(dn, 0, None).sum())

    det = load_detection_result(run_dir / "detection.npz")
    off = np.asarray(det.event_offsets)
    t_det = np.asarray(det.event_time_ps, dtype=float)
    i_det = np.repeat(np.arange(off.size - 1), np.diff(off))
    return (np.concatenate([t_win, t_det]), np.concatenate([i_win, i_det]), pickups)


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    rows: list[dict[str, Any]] = []
    missing: list[float] = []
    for s in sorted((*SEFF_GRID, SEFF_PRODUCTION)):
        rd = seff_run_dir(s, keep_ion=True)
        if not (rd / "ion.npz").exists() or not (rd / "detection.npz").exists():
            missing.append(s)
            continue
        times, ions, pickups = shed_times(rd)
        n_ions = 2 * RunDirectory(rd).load_cfg().num_molecules
        row: dict[str, Any] = {"s_eff": s, "sheds_total": int(times.size),
                               "sheds_per_ion": times.size / n_ions,
                               "pickups_window": pickups,
                               "t_shed_median_ps": float(np.median(times))}
        for thr in LATE_PS:
            row[f"frac_sheds_after_{thr:g}ps"] = float(np.mean(times > thr))
        row["frac_ions_shed_after_400ps"] = np.unique(ions[times > 400.0]).size / n_ions
        rows.append(row)

    show = tuple(rows[0]) if rows else ()
    print("\n=== B2 s_eff scan: shed timing (seed s1, N = 1000 per cell) ===")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    if missing:
        print(f"\nmissing: s_eff {missing}")
    print(f"\nCSV -> {write_rows_csv(OUT_CSV, rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
