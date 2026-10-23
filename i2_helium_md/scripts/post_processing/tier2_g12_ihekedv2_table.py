"""THESIS_FINAL_GOALS G12 — B2 / B3 re-scored against the corrected ihe_ked reference (v2).

The ihe_ked 3-D reference was re-exported on 2026-10-23 with the inversion
symmetry axis on the laser polarisation and the inversion centre on the
per-fragment centre (v2; ``data/reference/ihe_ked/README.md``). This script
re-scores the two thesis laws only — **B2** (production) and **B3**
(effective high-tail variant) — with the unchanged G5 read path
(``tier2_g5_winners_table``: same ``_row`` / ``tail_columns`` / law columns,
same pooled containers). No other law is re-scored (G12.2).

Oracle (before any new number): scoring each container against the archived
v1 reference (``archive_v1_2026-07-14/IHe_KED_reference.csv``) reproduces its
committed ``g5_winners_table.csv`` row to 4 dp on the anchor columns.

Output (commit with ``git add -f``):
``data/runs/h2b_forward_model/g5_winners_table_ihekedv2.csv`` — full vector
per law plus the experiment row (v2 reference). The committed
``g5_winners_table.csv`` stays as the v1 record.

Invocation::

    python scripts/post_processing/tier2_g12_ihekedv2_table.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.post_processing.tier2_detfix_battery_table import (  # noqa: E402
    _check_4dp,
    _csv_row,
    _row,
)
from scripts.post_processing.tier2_g413_table import ANCHOR_COLS  # noqa: E402
from scripts.post_processing.tier2_g415_table import tail_columns  # noqa: E402
from scripts.post_processing.tier2_g5_winners_table import (  # noqa: E402
    LAWS,
    OUT_CSV as G5_WINNERS_CSV,
    _law_columns,
    experiment_row,
)
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    RUNS_ROOT,
)
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "g5_winners_table_ihekedv2.csv"
IHE_KED_V1_CSV = IHE_KED_REFERENCE_CSV.parent / "archive_v1_2026-07-14" / "IHe_KED_reference.csv"
LABELS: tuple[str, ...] = ("B2_pp_c160", "B3_pp_c160_w0482_a198")

# ---------------------------------------------------------------------------


def _score(label: str, read, cfg, container: str, role: str, ab, ked) -> dict[str, Any]:
    return {**_row(label, "pooled", cfg.droplet_size_sampler_mode, "-", read, ab, ked),
            "role": role, **tail_columns(read, ked), **_law_columns(cfg),
            "num_molecules": int(cfg.num_molecules), "container": container,
            "ihe_ked_reference": "v2_2026-10-23"}


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_v1 = load_ihe_ked_reference(IHE_KED_V1_CSV)
    ked_v2 = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    rows: list[dict[str, Any]] = []
    for label in LABELS:
        container, role, _, _ = LAWS[label]
        rd = RUNS_ROOT / container
        cfg = RunDirectory(rd).load_cfg()
        read = read_confirmation_detection(load_detection_result(rd / "detection.npz"),
                                           label=label)
        v1 = _score(label, read, cfg, container, role, ab, ked_v1)
        _check_4dp(f"[v1 oracle] {label}", v1, _csv_row(G5_WINNERS_CSV, label), ANCHOR_COLS)
        rows.append(_score(label, read, cfg, container, role, ab, ked_v2))
    print(f"v1 oracle PASSED: {len(rows)} laws reproduce g5_winners_table.csv (4 dp)")
    rows.append({**experiment_row(ab, ked_v2), "ihe_ked_reference": "v2_2026-10-23"})

    keys = list(dict.fromkeys(k for r in rows for k in r))
    rows = [{k: r.get(k, "") for k in keys} for r in rows]
    show = ("label", "nbar_det", "n1_solv", "w1_solv", "midHot", "KE1_mean", "KE2_mean",
            "deepKE", "KE10_14_over_ref", "mass_n15plus_solv", "S", "chi2_med", "ke_npts")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    print(f"\nCSV -> {write_rows_csv(OUT_CSV, rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
