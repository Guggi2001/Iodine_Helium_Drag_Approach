"""Score the fixed-well lq runs (``gen_tier2_g415_well.py``) -> ``g415_table.csv``.

Full production observable vector + histogram band columns (the G4.14 scorer's
convention), plus the high-n tail and deep-bin KE read of D0 §1.G4.14.

Anchor (before any new row is written): the paired committed A5 member
``g414c200t64confs1`` rescored here reproduces its ``g414_table.csv`` row (4 dp).

Usage::

    python scripts/post_processing/tier2_g415_table.py
"""

from __future__ import annotations

from pathlib import Path
import sys
import warnings

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference  # noqa: E402
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection, solvated_renormalized,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2_g414_tauline import confirm_cells, run_dir_for as a5_run_dir_for  # noqa: E402
from scripts.build_pooled_detection_container import TARGETS  # noqa: E402
from scripts.gen_tier2_g415_well import (  # noqa: E402
    E0_EV, MEMBER, POOL_MEMBERS, WELL_ONLY_EV, run_dir_for, well_point,
)
from scripts.post_processing.tier2_detfix_battery_table import (  # noqa: E402
    _check_4dp, _csv_row, _pipeline, _row,
)
from scripts.post_processing.tier2_g413_table import ANCHOR_COLS  # noqa: E402
from scripts.post_processing.tier2_g414_table import _write, band_columns  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV, IHE_KED_REFERENCE_CSV, RUNS_ROOT,
)

OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "g415_table.csv"
G414_TABLE_CSV = RUNS_ROOT / "h2b_forward_model" / "g414_table.csv"


def tail_columns(read, ked) -> dict[str, float]:
    """High-n tail and deep-bin KE (D0 §1.G4.14): solvated n >= 15 mass, max
    scored n, and mean KE / reference mean KE over n = 10-14 (bins with > 3 ions)."""
    ns, ke = np.asarray(read.n_scored), np.asarray(read.ke_scored_eV)
    solv = ns >= 1
    ref = dict(zip(np.asarray(ked.n).astype(int), np.asarray(ked.mean_KE_eV, dtype=float)))
    ratios = [ke[ns == k].mean() / ref[k] for k in range(10, 15) if (ns == k).sum() > 3]
    return {"mass_n15plus_solv": float((ns[solv] >= 15).mean()),
            "n_max": int(ns.max()),
            "KE10_14_over_ref": float(np.mean(ratios)) if ratios else float("nan")}


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warnings
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    exp_n, exp_f = solvated_renormalized(ab.n, ab.ion_fraction, n_min=1)

    def scored(label: str, run_dir: Path, member: str = MEMBER) -> dict:
        cfg = RunDirectory(run_dir).load_cfg()
        read = read_confirmation_detection(load_detection_result(run_dir / "detection.npz"),
                                           label=label)
        sn, sf = solvated_renormalized(read.n_values, read.fraction, n_min=1)
        return {**_row(label, member, _pipeline(cfg), cfg.seed, read, ab, ked),
                **band_columns(sn, sf, exp_n, exp_f), **tail_columns(read, ked),
                "e_bind_eV": float(cfg.binding_energy_I_ion_eV),
                "a_amu_per_ps": float(cfg.drag_coefficients.coefficients["a"]),
                "E0_eV": float(cfg.internal_energy_partition_fraction * cfg.coulomb_available_eV)}

    a5 = next(c for c in confirm_cells() if c.label.endswith(f"conf{MEMBER}"))
    anchor = scored(a5.label, a5_run_dir_for(a5))
    _check_4dp(f"[anchor] {a5.label}", anchor, _csv_row(G414_TABLE_CSV, a5.label), ANCHOR_COLS)
    print(f"anchor PASS: {a5.label} reproduces g414_table.csv (4 dp)")
    rows = [anchor]
    wp = well_point()
    print(f"fixed well: a {wp['a']:.3f}, E_bind {wp['e_bind_eV']}")
    for label in E0_EV:
        rows.append({**scored(f"g415{label}{MEMBER}", run_dir_for(label)),
                     "arm": "tier0_consistent"})
    for label in WELL_ONLY_EV:   # informative, NOT Tier-0-consistent (A5 law, well only)
        if (run_dir_for(label) / "detection.npz").exists():
            rows.append({**scored(f"g415{label}{MEMBER}", run_dir_for(label)),
                         "arm": "well_only_informative"})
    # Well-only pool (2026-10-08): the other production members, then the
    # pooled N = 6000 container (build_pooled_detection_container.py).
    for label in WELL_ONLY_EV:
        for member in POOL_MEMBERS:
            if member == MEMBER or not (run_dir_for(label, member) / "detection.npz").exists():
                continue
            rows.append({**scored(f"g415{label}{member}", run_dir_for(label, member), member),
                         "arm": "well_only_informative"})
        pooled_dir = RUNS_ROOT / TARGETS[f"g415_{label}"][2]
        if (pooled_dir / "detection.npz").exists():
            rows.append({**scored(f"g415{label}pooled", pooled_dir, "pooled"),
                         "arm": "well_only_informative_pooled"})
    rows[0]["arm"] = "a5_anchor"
    cols = ("label", "arm", "e_bind_eV", "a_amu_per_ps", "E0_eV", "nbar_det", "n1_solv", "w1_solv",
            "midHot", "chi2_med", "KE1_mean", "supp", "trap", "mass_n10_14",
            "mass_n15plus_solv", "n_max", "KE10_14_over_ref", "deepKE")
    for r in rows:
        print("  " + "  ".join(f"{c} {r[c]:.4f}" if isinstance(r[c], float) else f"{c} {r[c]}"
                               for c in cols if c in r))
    _write(rows, OUT_CSV)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
