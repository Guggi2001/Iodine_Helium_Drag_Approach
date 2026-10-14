"""RRK s_eff scan at B2 (scorer for ``gen_tier2_b2_seff_scan.py``).

Rows: one per s_eff cell at seed s1 (s_eff = 8 is the committed B2 s1 run),
with the standard G5 observable vector plus the arrival-state split
(``state_reason`` fractions over all ions) and the shedding clock (median /
p90 time of each ion's last detection-stage shed event, over ions with at
least one event) — the "clock vs endpoint" question.

Oracle (failure stops before a new number): the s_eff = 8 row re-scored here
equals the committed ``g5lqppc160e3297s1`` row of ``g5_postpickup_table.csv``
(4 dp, ``ANCHOR_COLS``).

Output (commit with ``git add -f``; every column):
``data/runs/h2b_forward_model/b2_seff_table.csv``.

Invocation::

    python scripts/post_processing/tier2_b2_seff_table.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference  # noqa: E402
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import read_confirmation_detection  # noqa: E402
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2_b2_seff_scan import (  # noqa: E402
    MEMBER, SEFF_GRID, SEFF_PRODUCTION, seff_run_dir,
)
from scripts.gen_tier2_detfix_postpickup import B2_C, B2_E0_EV, B2_TAU_PS, lq_label  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table, write_rows_csv,
)
from scripts.post_processing.tier2_detfix_battery_table import _check_4dp, _csv_row, _row  # noqa: E402
from scripts.post_processing.tier2_g413_table import ANCHOR_COLS  # noqa: E402
from scripts.post_processing.tier2_g415_table import tail_columns  # noqa: E402
from scripts.post_processing.tier2_g5_postpickup_table import OUT_CSV as G5_CSV  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV, IHE_KED_REFERENCE_CSV, RUNS_ROOT,
)

OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "b2_seff_table.csv"
STATES = ("frozen", "time_exhausted", "suppressed")


def clock_columns(det) -> dict[str, float]:
    """Arrival-state fractions over all 2N ions and the last-shed time [ps]
    (median / p90 over ions with >= 1 detection-stage event)."""
    reason = np.asarray(det.state_reason)
    out = {f"frac_{s}": float(np.mean(reason == s)) for s in STATES}
    off = np.asarray(det.event_offsets)
    t = np.asarray(det.event_time_ps)
    has = off[1:] > off[:-1]
    last = t[off[1:][has] - 1] if t.size else np.empty(0)
    out["frac_with_event"] = float(np.mean(has))
    out["t_last_shed_med_ps"] = float(np.median(last)) if last.size else float("nan")
    out["t_last_shed_p90_ps"] = float(np.percentile(last, 90)) if last.size else float("nan")
    return out


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    rows: list[dict[str, Any]] = []
    missing: list[float] = []
    for s in sorted((*SEFF_GRID, SEFF_PRODUCTION)):
        rd = seff_run_dir(s)
        if not (rd / "detection.npz").exists():
            missing.append(s)
            continue
        cfg = RunDirectory(rd).load_cfg()
        if cfg.evap_rrk_dof != s:
            raise AssertionError(f"{rd.name}: evap_rrk_dof {cfg.evap_rrk_dof} != {s}")
        det = load_detection_result(rd / "detection.npz")
        read = read_confirmation_detection(det, label=f"seff{s:g}")
        row = _row(f"b2_seff{s:g}_{MEMBER}", MEMBER, cfg.droplet_size_sampler_mode,
                   cfg.seed, read, ab, ked)
        row["s_eff"] = s
        row.update(tail_columns(read, ked))
        row.update(clock_columns(det))
        if s == SEFF_PRODUCTION:
            _check_4dp("[anchor] B2 s1", row,
                       _csv_row(G5_CSV, lq_label(MEMBER, B2_E0_EV, B2_TAU_PS, B2_C)),
                       ANCHOR_COLS)
            print("B2 s1 anchor PASSED (4 dp)", flush=True)
        rows.append(row)

    keys = list(dict.fromkeys(k for r in rows for k in r))
    rows = [{k: r.get(k, "") for k in keys} for r in rows]
    show = ("s_eff", "nbar_det", "n1_solv", "w1_solv", "midHot", "chi2_med", "KE1_mean",
            "supp", "mass_n15plus_solv", "frac_frozen", "frac_time_exhausted",
            "t_last_shed_med_ps", "t_last_shed_p90_ps")
    show = tuple(c for c in show if c in keys)
    print("\n=== B2 s_eff scan (seed s1, N = 1000 per cell) ===")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    if missing:
        print(f"\nmissing: s_eff {missing}")
    print(f"\nCSV -> {write_rows_csv(OUT_CSV, rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
