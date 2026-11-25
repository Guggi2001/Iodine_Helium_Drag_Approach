"""THESIS_DISCUSSION §15 — per-n profile and extra observables, B2 / B3 / h405 vs experiment.

Zero-MD. Reads the thesis production pools (``data/runs/{B2,B3,h405}_N6000_pooled``,
``self_bound_gate_release = decoupled``, THESIS_FINAL_GOALS G13) with the
scorer's reader (``read_confirmation_detection``, ``solvated_renormalized``)
and the v2 references the scorer uses (``integrated_i_he_abundance.csv``,
``IHe_KED_reference.csv``). Writes the two blocks of THESIS_DISCUSSION §15
that are not columns of ``g13_decoupled_table.csv``:

* per n = 1-20: solvated share f(n) (share of n among n >= 1) and mean KE(n)
  [eV] with the scored-ion count; KE is blank below the scorer's minimum bin
  count ``MIN_KE_BIN_COUNT``;
* per law: bare I+ share of scored ions, solvated mean n (n >= 1), the
  KE1 / KE2 median (reference: also the mode), and the bare-ion (n = 0)
  mean / median KE [eV] -- a display quantity (n = 0 is not scored).
  Simulated medians: ``np.median`` of the scored-ion KE at that n.

Reported: where a law's f(1), KE(1) or KE1 / KE2 median differ from the
``n1_solv`` / ``KE1_mean`` / ``KE1_med`` / ``KE2_med`` of its row in
``g13_decoupled_table.csv`` at 4 dp.

Output (commit with ``git add -f``):
``data/runs/h2b_forward_model/thesis_observables_per_n.csv`` and
``data/runs/h2b_forward_model/thesis_observables_extra.csv``.

Invocation::

    python scripts/post_processing/tier2_thesis_observable_table.py
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
from i2_helium_md.postprocess.battery_score import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    MIN_KE_BIN_COUNT,
    write_rows_csv,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
    solvated_renormalized,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.thesis_runs import (  # noqa: E402
    TABLES_ROOT,
    pool_dir,
    require_run,
)
from scripts.post_processing.anchors import report, row_differences  # noqa: E402

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

LAWS = ("B2", "B3", "h405")
N_MAX = 20                          # the abundance reference ends at n = 20
SCORE_TABLE_CSV = TABLES_ROOT / "g13_decoupled_table.csv"
OUT_PER_N_CSV = TABLES_ROOT / "thesis_observables_per_n.csv"
OUT_EXTRA_CSV = TABLES_ROOT / "thesis_observables_extra.csv"

# ---------------------------------------------------------------------------


def _solvated(n_values, fraction) -> tuple[np.ndarray, np.ndarray]:
    sn, sf = solvated_renormalized(n_values, fraction, n_min=1)
    return np.asarray(sn), np.asarray(sf, dtype=float)


def scored_ke_median(read, n: int) -> float:
    """Median detected KE [eV] of the scored ions at ``n`` -- the
    ``lowke_columns`` recipe behind the committed ``KE{n}_med`` columns."""
    return float(np.median(np.asarray(read.ke_scored_eV)[np.asarray(read.n_scored) == n]))


def anchor_differences(per_n: list[dict[str, Any]], extra: list[dict[str, Any]]) -> list[str]:
    """f(1), KE(1) and the KE1 / KE2 medians that differ from each law's
    committed row of ``SCORE_TABLE_CSV`` at 4 dp."""
    n1 = per_n[0]
    medians = {r["label"]: r for r in extra}
    out = []
    for law in LAWS:
        out += row_differences(
            {"n1_solv": n1[f"f_{law}"], "KE1_mean": n1[f"KE_{law}_eV"],
             "KE1_med": medians[law]["KE1_median_eV"],
             "KE2_med": medians[law]["KE2_median_eV"]}, SCORE_TABLE_CSV, law)
    return out


def build_tables() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """The per-n rows and the per-law rows, from the thesis pools."""
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    exp_n, exp_f = _solvated(ab.n, ab.ion_fraction)
    ref = {name: {int(k): float(v) for k, v in zip(ked.n, getattr(ked, name))}
           for name in ("mean_KE_eV", "median_KE_eV", "mode_KE_eV")}
    reads = {law: read_confirmation_detection(
        load_detection_result(require_run(pool_dir(law) / "detection.npz", "gen_thesis_laws.py")), label=law)
        for law in LAWS}

    per_n: list[dict[str, Any]] = []
    for n in range(1, N_MAX + 1):
        row: dict[str, Any] = {"n": n, "f_exp": float(exp_f[exp_n == n].sum())}
        for law, r in reads.items():
            sn, sf = _solvated(r.n_values, r.fraction)
            row[f"f_{law}"] = float(sf[sn == n].sum())
        row["KE_exp_eV"] = ref["mean_KE_eV"].get(n, "")
        for law, r in reads.items():
            ke = np.asarray(r.ke_scored_eV)[np.asarray(r.n_scored) == n]
            row[f"KE_{law}_eV"] = float(ke.mean()) if ke.size >= MIN_KE_BIN_COUNT else ""
            row[f"count_{law}"] = int(ke.size)
        per_n.append(row)

    an, af = np.asarray(ab.n), np.asarray(ab.ion_fraction, dtype=float)
    extra: list[dict[str, Any]] = [{
        "label": "experiment",
        "bare_share": float(af[an == 0].sum()),
        "mean_n_solvated": float((exp_n * exp_f).sum()),
        "KE1_median_eV": ref["median_KE_eV"][1], "KE1_mode_eV": ref["mode_KE_eV"][1],
        "KE2_median_eV": ref["median_KE_eV"][2], "KE2_mode_eV": ref["mode_KE_eV"][2],
        "KE0_mean_eV": ref["mean_KE_eV"][0], "KE0_median_eV": ref["median_KE_eV"][0],
    }]
    for law, r in reads.items():
        ns = np.asarray(r.n_scored)
        sn, sf = _solvated(r.n_values, r.fraction)
        row = {"label": law, "bare_share": float((ns == 0).mean()),
               "mean_n_solvated": float((sn * sf).sum())}
        for n in (1, 2):
            row[f"KE{n}_median_eV"] = scored_ke_median(r, n)
        row["KE0_mean_eV"] = float(np.asarray(r.ke_scored_eV)[ns == 0].mean())
        row["KE0_median_eV"] = scored_ke_median(r, 0)
        extra.append(row)
    keys = list(dict.fromkeys(k for r in extra for k in r))
    return per_n, [{k: r.get(k, "") for k in keys} for r in extra]


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    per_n, extra = build_tables()
    for r in per_n:
        print(" ".join(f"{k}={v:.4f}" if isinstance(v, float) else f"{k}={v}"
                       for k, v in r.items()))
    for r in extra:
        print(r)
    report(anchor_differences(per_n, extra), f"anchor OK vs {SCORE_TABLE_CSV.name} (4 dp)")
    print(f"\nCSV -> {write_rows_csv(OUT_PER_N_CSV, per_n)}")
    print(f"CSV -> {write_rows_csv(OUT_EXTRA_CSV, extra)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
