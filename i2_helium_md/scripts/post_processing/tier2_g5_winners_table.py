"""THESIS_FINAL_GOALS G5 — one comparison table of the thesis laws (pooled containers).

Scores the committed pooled figure containers with one read path
(``read_confirmation_detection`` on each container's ``detection.npz``), so
every column is computed the same way for every law:

* A5 (lq c 2.0, raw sampler) — previous TDDFT-based production law;
* h405 (capped cubic, raw / post_pickup) — the identifiability contrast;
* G10 effective (A5 + E_bind 0.0482, a 19.5, raw sampler);
* A5 with post_pickup, parameters unchanged;
* **B2** (post_pickup, lq c 1.6) — TDDFT-based production law (2026-10-09);
* **B3** (B2 + E_bind 0.0482, a 19.8) — effective high-tail variant.

Oracle (before any new number): each container row reproduces its committed
pooled row (4 dp on the anchor columns) — A5 / B2 / B3 / A5-pp against
``g5_postpickup_table.csv``, h405 against ``detfix_battery_table.csv``
(``pooled_new6000``), G10 against ``g415_table.csv``.

Output (commit with ``git add -f``):
``data/runs/h2b_forward_model/g5_winners_table.csv`` — full vector per law
plus an experiment row (reference-derived columns only).

Invocation::

    python scripts/post_processing/tier2_g5_winners_table.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
    solvated_renormalized,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.post_processing.tier2_detfix_battery_table import (  # noqa: E402
    TABLE_CSV as DETFIX_TABLE_CSV,
    _check_4dp,
    _csv_row,
    _row,
)
from scripts.post_processing.tier2_g413_table import ANCHOR_COLS  # noqa: E402
from scripts.post_processing.tier2_g415_table import tail_columns  # noqa: E402
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

OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "g5_winners_table.csv"
G5_CSV = RUNS_ROOT / "h2b_forward_model" / "g5_postpickup_table.csv"
G415_CSV = RUNS_ROOT / "h2b_forward_model" / "g415_table.csv"
_PREFIX = "9A_drag_srho3p5A18A_linear_quadratic_N6000_detfix_conf270_"

#: label -> (container dir, role, anchor CSV, anchor row label or None).
LAWS: dict[str, tuple[str, str, Path, Optional[str]]] = {
    "A5_raw": (_PREFIX + "g414c200t64pooled_th500", "previous TDDFT-based production",
               G5_CSV, "pooled_lqraw_6000"),
    "h405_raw": ("9A_drag_shared_pure_cubic_N6000_detfix_conf270_h405pooled_th500",
                 "identifiability contrast", DETFIX_TABLE_CSV, "pooled_new6000"),
    "h405_pp5000": ("9A_drag_shared_pure_cubic_N5000_detfix_conf270_h405pppooled_th500",
                    "contrast, post_pickup (s1-s5)", RUNS_ROOT / "h2b_forward_model"
                    / "detfix_postpickup_table.csv", "pooled_pp5000"),
    "G10_raw_w0482_a195": (_PREFIX + "g415wo048a195pooled_th500",
                           "previous effective high-tail", G415_CSV, "g415wo048a195pooled"),
    "A5_pp": (_PREFIX + "g5lqpppooled_th500", "post_pickup, unchanged",
              G5_CSV, "pooled_lqpp_6000"),
    "B2_pp_c160": (_PREFIX + "g5lqppc160e3297pooled_th500",
                   "TDDFT-based production (2026-10-09)", G5_CSV, "pooled_lqppb2_6000"),
    "B3_pp_c160_w0482_a198": (_PREFIX + "g5lqppc160a198w0482e3389pooled_th500",
                              "effective high-tail variant (2026-10-09)",
                              G5_CSV, "pooled_lqppb3_6000"),
}

# ---------------------------------------------------------------------------


def _law_columns(cfg) -> dict[str, Any]:
    co = cfg.drag_coefficients.coefficients
    return {"drag_form": cfg.drag_form, "drag_coefficients": dict(co),
            "e_bind_eV": float(cfg.binding_energy_I_ion_eV),
            "E0_eV": float(cfg.internal_energy_partition_fraction * cfg.coulomb_available_eV),
            "tau_ps": float(cfg.internal_energy_cooling_tau_ps),
            "sampler": cfg.droplet_size_sampler_mode,
            "gate_A": float(cfg.drag_gate_steepness)}


def experiment_row(ab, ked) -> dict[str, Any]:
    """Reference-derived columns: solvated n1 share, n >= 15 mass, mean KE at
    n = 1 / 2 [eV]; midHot / deepKE / KE n10-14 ratios are 1 and W1 0 by
    construction."""
    sn, sf = solvated_renormalized(ab.n, ab.ion_fraction, n_min=1)
    sn, sf = np.asarray(sn), np.asarray(sf, dtype=float)
    ke = dict(zip(np.asarray(ked.n).astype(int), np.asarray(ked.mean_KE_eV, dtype=float)))
    return {"label": "experiment", "role": "reference",
            "n1_solv": float(sf[sn == 1].sum()), "w1_solv": 0.0, "midHot": 1.0,
            "deepKE": 1.0, "KE10_14_over_ref": 1.0,
            "mass_n15plus_solv": float(sf[sn >= 15].sum()),
            "KE1_mean": ke.get(1, float("nan")), "KE2_mean": ke.get(2, float("nan"))}


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    rows: list[dict[str, Any]] = []
    for label, (container, role, anchor_csv, anchor_label) in LAWS.items():
        rd = RUNS_ROOT / container
        cfg = RunDirectory(rd).load_cfg()
        read = read_confirmation_detection(load_detection_result(rd / "detection.npz"),
                                           label=label)
        row = {**_row(label, "pooled", cfg.droplet_size_sampler_mode, "-", read, ab, ked),
               "role": role, **tail_columns(read, ked), **_law_columns(cfg),
               "num_molecules": int(cfg.num_molecules), "container": container}
        _check_4dp(f"[anchor] {label}", row, _csv_row(anchor_csv, anchor_label), ANCHOR_COLS)
        rows.append(row)
    print(f"anchors PASSED: {len(rows)} containers reproduce their committed pooled rows (4 dp)")
    rows.append(experiment_row(ab, ked))

    keys = list(dict.fromkeys(k for r in rows for k in r))
    rows = [{k: r.get(k, "") for k in keys} for r in rows]
    show = ("label", "num_molecules", "nbar_det", "n1_solv", "w1_solv", "midHot", "KE1_mean",
            "KE2_mean", "deepKE", "KE10_14_over_ref", "mass_n15plus_solv", "n_max", "trap",
            "supp", "S", "chi2_med")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    print(f"\nCSV -> {write_rows_csv(OUT_CSV, rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
