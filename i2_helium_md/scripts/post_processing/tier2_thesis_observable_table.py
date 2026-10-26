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
* per law: bare I+ share of scored ions, solvated mean n (n >= 1), and the
  reference KE1 / KE2 median and mode.

Oracle (before anything is written): each law's f(1) and KE(1) reproduce the
``n1_solv`` / ``KE1_mean`` of its decoupled row in ``g13_decoupled_table.csv``
to 4 dp.

Output (commit with ``git add -f``):
``data/runs/h2b_forward_model/thesis_observables_per_n.csv`` and
``data/runs/h2b_forward_model/thesis_observables_extra.csv``.

Invocation::

    python scripts/post_processing/tier2_thesis_observable_table.py
"""

from __future__ import annotations

import csv
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference  # noqa: E402
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
    solvated_renormalized,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from scripts.gen_tier2_detfix_postpickup import thesis_pool_dir  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import write_rows_csv  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    MIN_KE_BIN_COUNT,
    RUNS_ROOT,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

LAWS = ("B2", "B3", "h405")
N_MAX = 20                          # the abundance reference ends at n = 20
SCORE_TABLE_CSV = RUNS_ROOT / "h2b_forward_model" / "g13_decoupled_table.csv"
OUT_PER_N_CSV = RUNS_ROOT / "h2b_forward_model" / "thesis_observables_per_n.csv"
OUT_EXTRA_CSV = RUNS_ROOT / "h2b_forward_model" / "thesis_observables_extra.csv"

# ---------------------------------------------------------------------------


def _solvated(n_values, fraction) -> tuple[np.ndarray, np.ndarray]:
    sn, sf = solvated_renormalized(n_values, fraction, n_min=1)
    return np.asarray(sn), np.asarray(sf, dtype=float)


def _check_oracle(per_n: list[dict[str, Any]]) -> None:
    with SCORE_TABLE_CSV.open(newline="") as f:
        scored = {r["label"]: r for r in csv.DictReader(f)}
    n1 = per_n[0]
    bad = [f"{law} {col}: {float(n1[key]):.4f} vs {float(scored[law][col]):.4f}"
           for law in LAWS
           for key, col in ((f"f_{law}", "n1_solv"), (f"KE_{law}_eV", "KE1_mean"))
           if round(float(n1[key]), 4) != round(float(scored[law][col]), 4)]
    if bad:
        raise RuntimeError("oracle vs g13_decoupled_table.csv FAILED: " + "; ".join(bad))


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    exp_n, exp_f = _solvated(ab.n, ab.ion_fraction)
    ref = {name: {int(k): float(v) for k, v in zip(ked.n, getattr(ked, name))}
           for name in ("mean_KE_eV", "median_KE_eV", "mode_KE_eV")}
    reads = {law: read_confirmation_detection(
        load_detection_result(thesis_pool_dir(law) / "detection.npz"), label=law)
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
    _check_oracle(per_n)

    an, af = np.asarray(ab.n), np.asarray(ab.ion_fraction, dtype=float)
    extra: list[dict[str, Any]] = [{
        "label": "experiment",
        "bare_share": float(af[an == 0].sum()),
        "mean_n_solvated": float((exp_n * exp_f).sum()),
        "KE1_median_eV": ref["median_KE_eV"][1], "KE1_mode_eV": ref["mode_KE_eV"][1],
        "KE2_median_eV": ref["median_KE_eV"][2], "KE2_mode_eV": ref["mode_KE_eV"][2],
    }]
    for law, r in reads.items():
        ns = np.asarray(r.n_scored)
        sn, sf = _solvated(r.n_values, r.fraction)
        extra.append({"label": law, "bare_share": float((ns == 0).mean()),
                      "mean_n_solvated": float((sn * sf).sum())})
    keys = list(dict.fromkeys(k for r in extra for k in r))
    extra = [{k: r.get(k, "") for k in keys} for r in extra]

    for r in per_n:
        print(" ".join(f"{k}={v:.4f}" if isinstance(v, float) else f"{k}={v}"
                       for k, v in r.items()))
    for r in extra:
        print(r)
    print(f"oracle PASSED: f(1) / KE(1) reproduce {SCORE_TABLE_CSV.name} (4 dp)")
    print(f"\nCSV -> {write_rows_csv(OUT_PER_N_CSV, per_n)}")
    print(f"CSV -> {write_rows_csv(OUT_EXTRA_CSV, extra)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
