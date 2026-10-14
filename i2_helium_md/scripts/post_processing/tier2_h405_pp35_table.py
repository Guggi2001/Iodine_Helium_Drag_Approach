"""h405 at both production settings (``post_pickup`` + 3.5 A gate) -- scorer.

Scores ``gen_tier2_h405_pp35.py`` runs: the E0 column at s1 and, once run,
the pooled s1-s5 + r6 battery at ``POOL_E0_EV``; adds the B2 pool and the
committed h405 one-change-at-a-time pools as comparison rows.

Oracle (failure stops before a new number): the committed B2 s1 member
re-scored here equals its ``g5_postpickup_table.csv`` row (4 dp,
``ANCHOR_COLS``) -- same scorer path as the new rows.

Comparison rows copied from committed CSVs (not re-scored; their member run
dirs are no longer on disk): ``h405_pooled6000`` and ``h405g35_pooled6000``
(``g414_table.csv``), ``pooled_pp5000`` (``detfix_postpickup_table.csv``).

Output (commit with ``git add -f``): ``data/runs/h2b_forward_model/
h405_pp35_table.csv``.

Invocation::

    python scripts/post_processing/tier2_h405_pp35_table.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference  # noqa: E402
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    pool_confirmation_reads, read_confirmation_detection,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2_detfix_battery import MEMBER_SEEDS  # noqa: E402
from scripts.gen_tier2_detfix_postpickup import (  # noqa: E402
    B2_C, B2_E0_EV, B2_TAU_PS, lq_label, lq_run_dir,
)
from scripts.gen_tier2_g413_ridge_battery import _E0_KEY  # noqa: E402
from scripts.gen_tier2_h405_pp35 import (  # noqa: E402
    COLUMN_E0_EV, COLUMN_MEMBER, POOL_E0_EV, run_dir_for,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table, write_rows_csv,
)
from scripts.post_processing.tier2_detfix_battery_table import _check_4dp, _csv_row, _row  # noqa: E402
from scripts.post_processing.tier2_detfix_postpickup_table import OUT_CSV as H405_AB_CSV  # noqa: E402
from scripts.post_processing.tier2_g413_table import ANCHOR_COLS  # noqa: E402
from scripts.post_processing.tier2_g415_table import tail_columns  # noqa: E402
from scripts.post_processing.tier2_g5_postpickup_table import OUT_CSV as G5_CSV  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV, IHE_KED_REFERENCE_CSV, RUNS_ROOT,
)

OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "h405_pp35_table.csv"
G414_CSV = RUNS_ROOT / "h2b_forward_model" / "g414_table.csv"
COMMITTED_ROWS = ((G414_CSV, "h405_pooled6000", "h405 raw 14.2A (committed)"),
                  (G414_CSV, "h405g35_pooled6000", "h405 raw 3.5A (committed)"),
                  (H405_AB_CSV, "pooled_pp5000", "h405 pp 14.2A not re-landed (committed)"))


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    def scored(label: str, member: str, rd: Path):
        cfg = RunDirectory(rd).load_cfg()
        read = read_confirmation_detection(load_detection_result(rd / "detection.npz"),
                                           label=label)
        row = _row(label, member, cfg.droplet_size_sampler_mode, cfg.seed, read, ab, ked)
        row["E0_eV"] = float(getattr(cfg, _E0_KEY) * cfg.coulomb_available_eV)
        row["gate_A"] = float(cfg.drag_gate_steepness)
        row.update(tail_columns(read, ked))
        return row, read

    def pooled(label: str, reads, cfg_dir: Path) -> dict[str, Any]:
        cfg = RunDirectory(cfg_dir).load_cfg()
        pr = pool_confirmation_reads(reads, label=label)
        row = _row(label, f"pooled{1000 * len(reads)}", cfg.droplet_size_sampler_mode,
                   "-", pr, ab, ked)
        row["E0_eV"] = float(getattr(cfg, _E0_KEY) * cfg.coulomb_available_eV)
        row["gate_A"] = float(cfg.drag_gate_steepness)
        row.update(tail_columns(pr, ked))
        return row

    # Oracle: B2 s1 re-scored == committed row.
    b2 = {m: scored(lq_label(m, B2_E0_EV, B2_TAU_PS, B2_C), m,
                    lq_run_dir(m, B2_E0_EV, B2_TAU_PS, B2_C)) for m in MEMBER_SEEDS}
    _check_4dp("[anchor] B2 s1", b2["s1"][0], _csv_row(G5_CSV, b2["s1"][0]["label"]),
               ANCHOR_COLS)
    print("B2 s1 anchor PASSED (4 dp)", flush=True)

    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for e0 in COLUMN_E0_EV:
        rd = run_dir_for(COLUMN_MEMBER, e0)
        if not (rd / "detection.npz").exists():
            missing.append(rd.name)
            continue
        rows.append(scored(f"h405pp35_e{e0:.4f}_{COLUMN_MEMBER}", COLUMN_MEMBER, rd)[0])
    if POOL_E0_EV is not None:
        reads = []
        for m in MEMBER_SEEDS:
            rd = run_dir_for(m, POOL_E0_EV)
            if not (rd / "detection.npz").exists():
                missing.append(rd.name)
                continue
            row, read = scored(f"h405pp35_e{POOL_E0_EV:.4f}_{m}", m, rd)
            rows.append(row)
            reads.append(read)
        if reads:
            rows.append(pooled(f"h405pp35_pooled{1000 * len(reads)}", reads,
                               run_dir_for("s1", POOL_E0_EV)))
    rows.append(pooled("B2_pooled6000", [r for _, r in b2.values()],
                       lq_run_dir("s1", B2_E0_EV, B2_TAU_PS, B2_C)))
    for path, src, label in COMMITTED_ROWS:
        rows.append({**_csv_row(path, src), "label": label})

    keys = list(dict.fromkeys(k for r in rows for k in r))
    rows = [{k: r.get(k, "") for k in keys} for r in rows]
    show = ("label", "E0_eV", "gate_A", "num_scored", "nbar_det", "n1_solv", "w1_solv",
            "midHot", "chi2_med", "KE1_mean", "deepKE", "supp", "trap",
            "mass_n15plus_solv", "n_max")
    show = tuple(c for c in show if c in keys)
    print("\n=== h405 at post_pickup + 3.5 A (vs B2 and committed h405 pools) ===")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    if missing:
        print(f"\nmissing: {missing}")
    print(f"\nCSV -> {write_rows_csv(OUT_CSV, rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
