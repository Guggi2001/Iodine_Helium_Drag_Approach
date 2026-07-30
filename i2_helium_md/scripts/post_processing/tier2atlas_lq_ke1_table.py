"""Atlas §6.7 item 3 — the lq KE₁ read (closure follow-up, 2026-07-30).

States the n = 1 / n = 2 KE observables (mean / median / mode / SD /
above-1.15 share) of the §6.7 lq sanity battery, paired-by-seed against
the standing cubic battery — the number the §6.7 item-1 record left
implicit (it reported the KED χ² but not KE₁ itself). Context: the (C)
closure discussion (findings "(C) adjudication CLOSED") asked whether a
lower-power drag form could raise the low-n KE; item 1 measured lq's KED
χ² worse on all five seeds, and this instrument turns that inference into
the KE₁ number directly. Pure scorer — runs nothing, mutates nothing.

The lq run dirs are local-only and were cleaned after §6.7; they are
regenerated seed-exactly by the committed ``gen_tier2atlas_lqbattery.py``
(frozen seeds 20260722–26, paired cfg-diff guard). Regeneration fidelity
is gated by oracle O3 below before any KE number is read.

Convention note (why O3 is column-restricted): the §6.7 findings numbers
came from the session-local ``atlas_spotcheck_score.py`` battery driver,
which is not repo-kept. Its W₁ / supp / trap / n̄ / n₁_solv conventions
match the committed `observable_columns` (the §4cc pooled oracle chain:
W₁ 0.5713 both), but its midHot (mean-anchored, 1.0139 pooled cubic) and
median-anchored χ² differ from the committed row (`midHot` 1.0110,
``chi2_med`` = profiled). O3 therefore anchors ONLY on the
convention-shared columns; midHot / χ²_med are reported, not oracled.

Order of operations (§1.4 — oracles first, hard-fail before any new
number):

1. **O1 scorer-drift** — the standing pooled battery reproduces its
   recorded observable row (`check_oracle`).
2. **O2 KE-read** — the committed ``g4fh405`` run rescored through this
   script's KE path equals the committed ``atlas_ce_probe.csv`` h405 row
   to 4 decimals on KE1_mean / KE1_med / KE1_sd / n1_solv.
3. **O3 committed §6.7 rows** — every cubic member reproduces its
   committed findings-table W₁ and supp to the printed 3 decimals
   (instrument consistency on the original dirs), and every REGENERATED
   lq member reproduces its committed W₁ and supp the same way, plus the
   pooled lq row (num_scored exact; trap / supp / n̄ / n₁_solv / W₁ at
   the printed precision) — the regeneration-fidelity gate.

Invocation:

    python scripts/post_processing/tier2atlas_lq_ke1_table.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any

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
from scripts.gen_tier2atlas_g4finals import finals_run_dir_name  # noqa: E402
from scripts.gen_tier2atlas_lqbattery import (  # noqa: E402
    BATTERY_MATRIX,
    atlas_run_dir_name,
)
from scripts.post_processing.tier2atlas_g4step2_battery_table import (  # noqa: E402
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
    format_table,
    write_rows_csv,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

SAVE_CSV_PATH = RUNS_ROOT / "h2b_forward_model" / "atlas_lq_ke1_table.csv"

CE_PROBE_CSV = RUNS_ROOT / "h2b_forward_model" / "atlas_ce_probe.csv"

ABOVE_EV = 1.15                 # the n = 1 KED upper-share threshold (§3.5j)

# Reference n = 1 / n = 2 KED anchors (committed: findings §3.5j read (ii),
# g4step2 LOWKE_REF; display context only, no verdict here).
REF_N1 = {"mean": 1.302, "median": 1.128, "mode": 0.891, "above115": 0.487}
REF_N2_MEAN = 0.706

# O2 columns compared to 4 decimals against the committed h405 probe row.
O2_COLS = ("KE1_mean", "KE1_med", "KE1_sd", "n1_solv")

# O3 — the committed §6.7 item-1 findings table, compared at PRINTED
# precision (format the regenerated value with the committed table's format
# and require string equality). A numeric half-band tolerance breaks at the
# printed-precision boundary: the regenerated pooled trap is exactly
# 365/10000 = 0.0365, which Python prints as "0.036" (binary sits under the
# half) — |0.0365 − 0.036| ≤ 0.0005 then fails by 4e-19. String-format
# comparison is the faithful encoding of "reproduces the printed value".
# Per-member (w1_solv, supp) for both arms; midHot/χ² excluded (convention).
COMMITTED_667_MEMBERS: dict[str, tuple[str, str]] = {
    "bigc1v725s1": ("0.728", "0.192"),
    "bigc1v725s2": ("0.606", "0.185"),
    "bigc1v725s3": ("0.519", "0.188"),
    "bigc1v725s4": ("0.560", "0.192"),
    "bigc1v725s5": ("0.482", "0.180"),
    "qccbigs1": ("0.692", "0.217"),
    "qccbigs2": ("0.615", "0.213"),
    "qccbigs3": ("0.504", "0.206"),
    "qccbigs4": ("0.481", "0.219"),
    "qccbigs5": ("0.487", "0.209"),
}

# O3 — the committed §6.7 pooled lq row (concat 5 × N = 1000), as
# (format spec, committed printed value).
COMMITTED_667_LQ_POOL: dict[str, tuple[str, str]] = {
    "num_scored": ("d", "9635"), "trap": (".3f", "0.036"),
    "supp": (".3f", "0.213"), "nbar_det": (".3f", "4.008"),
    "n1_solv": (".4f", "0.2447"), "w1_solv": (".4f", "0.5484"),
}


def ke_row_columns(read) -> dict[str, Any]:
    """The KE observables of one read: lowke + SD + the above-1.15 share."""
    out: dict[str, Any] = {}
    out.update(lowke_columns(read))
    out.update(sn_ke_columns(read))
    sel = np.asarray(read.ke_scored_eV)[np.asarray(read.n_scored) == 1]
    out["KE1_n"] = int(sel.size)
    out["KE1_above115"] = (
        float((sel > ABOVE_EV).mean()) if sel.size >= 5 else float("nan")
    )
    return out


def committed_h405_probe_row() -> dict[str, str]:
    """The committed (C)-probe h405 row (string fields)."""
    with open(CE_PROBE_CSV, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row["label"] == "h405":
                return row
    raise AssertionError(f"no h405 row in {CE_PROBE_CSV}")


def member_row(label: str, run_dir: Path, seed: int, form: str,
               abundance_ref, ked_ref):
    """One battery-member row + its read (for pooling)."""
    read = load_confirmation_run(run_dir, label=label)
    row: dict[str, Any] = {"label": label, "form": form, "seed": seed}
    row.update(observable_columns(read, abundance_ref, ked_ref))
    row.update(ke_row_columns(read))
    return row, read


def check_667_member(row: dict[str, Any]) -> list[str]:
    """O3 per-member drift vs the committed §6.7 (W₁, supp) pair."""
    w1_c, supp_c = COMMITTED_667_MEMBERS[str(row["label"])]
    return [
        f"{row['label']} {col}: got {float(row[col]):.4f} vs committed "
        f"{committed}"
        for col, committed in (("w1_solv", w1_c), ("supp", supp_c))
        if f"{float(row[col]):.3f}" != committed
    ]


def main() -> None:
    abundance_ref = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked_ref = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    # --- O1: scorer drift (standing pooled battery). -----------------------
    oracle_row = observable_row(
        "pooled(oracle)", RUNS_ROOT / ORACLE_RUN, abundance_ref, ked_ref,
        with_geometry=False,
    )
    drift = check_oracle(oracle_row)
    if drift:
        print("*** O1 SCORER-DRIFT ORACLE FAILED — do not read ***")
        for line in drift:
            print(f"    {line}")
        return
    print("O1 scorer-drift oracle OK (pooled standing battery reproduces).")

    # --- O2: the KE read path vs the committed h405 probe row. -------------
    committed = committed_h405_probe_row()
    h405_read = load_confirmation_run(
        RUNS_ROOT / finals_run_dir_name("h405"), label="g4fh405")
    rescored: dict[str, Any] = {}
    rescored.update(observable_columns(h405_read, abundance_ref, ked_ref))
    rescored.update(ke_row_columns(h405_read))
    bad = [
        f"{c}: rescored {float(rescored[c]):.4f} vs committed "
        f"{float(committed[c]):.4f}"
        for c in O2_COLS
        if f"{float(rescored[c]):.4f}" != f"{float(committed[c]):.4f}"
    ]
    if bad:
        print("*** O2 KE-READ ORACLE FAILED — do not read ***")
        for line in bad:
            print(f"    {line}")
        return
    print("O2 KE-read oracle OK (committed h405 probe row to 4 decimals).")

    # --- Score both arms (members must all exist). -------------------------
    missing = [
        name for spec in BATTERY_MATRIX
        for name in (spec.paired_cubic_run_dir, atlas_run_dir_name(spec.label))
        if not (RUNS_ROOT / name / "detection.npz").exists()
    ]
    if missing:
        print("battery incomplete — stopping before any KE number:")
        for name in missing:
            print(f"    missing {name}")
        return

    rows: list[dict[str, Any]] = []
    o3_bad: list[str] = []
    reads = {"cubic": [], "lq": []}
    for spec in BATTERY_MATRIX:
        cub_label = spec.paired_cubic_run_dir.rsplit("_", 1)[-1]
        for form, label, dir_name in (
            ("cubic", cub_label, spec.paired_cubic_run_dir),
            ("lq", spec.label, atlas_run_dir_name(spec.label)),
        ):
            row, read = member_row(label, RUNS_ROOT / dir_name, spec.seed,
                                   form, abundance_ref, ked_ref)
            o3_bad += check_667_member(row)
            rows.append(row)
            reads[form].append(read)

    # --- O3: committed §6.7 reproduction (incl. regeneration fidelity). ----
    lq_pool_read = pool_confirmation_reads(reads["lq"], label="lq_pool")
    lq_pool: dict[str, Any] = {"label": "lq_pool", "form": "lq", "seed": "-"}
    lq_pool.update(observable_columns(lq_pool_read, abundance_ref, ked_ref))
    lq_pool.update(ke_row_columns(lq_pool_read))
    o3_bad += [
        f"lq_pool {col}: got {lq_pool[col]} vs committed {committed}"
        for col, (spec, committed) in COMMITTED_667_LQ_POOL.items()
        if format(int(lq_pool[col]) if spec == "d" else float(lq_pool[col]),
                  spec) != committed
    ]
    if o3_bad:
        print("*** O3 §6.7 COMMITTED-ROW ORACLE FAILED — regeneration or "
              "instrument drift; do not read KE ***")
        for line in o3_bad:
            print(f"    {line}")
        return
    print("O3 §6.7 committed rows OK (10 members W₁+supp to 3 decimals; "
          "pooled lq row reproduced — regeneration is faithful).\n")

    cub_pool_read = pool_confirmation_reads(reads["cubic"], label="cubic_pool")
    cub_pool: dict[str, Any] = {"label": "cubic_pool", "form": "cubic",
                                "seed": "-"}
    cub_pool.update(observable_columns(cub_pool_read, abundance_ref, ked_ref))
    cub_pool.update(ke_row_columns(cub_pool_read))
    rows.append(cub_pool)
    rows.append(lq_pool)

    print("=== §6.7 item 3 — the lq KE₁ read (paired 5 × N = 1000) ===")
    print(format_table(rows))

    # --- Paired per-seed deltas (lq − cubic). ------------------------------
    print("\n=== Paired Δ (lq − cubic), per seed ===")
    deltas: dict[str, list[float]] = {}
    for spec in BATTERY_MATRIX:
        cub_label = spec.paired_cubic_run_dir.rsplit("_", 1)[-1]
        cub = next(r for r in rows if r["label"] == cub_label)
        lq = next(r for r in rows if r["label"] == spec.label)
        line = [f"seed {spec.seed}:"]
        for col in ("KE1_mean", "KE1_sd", "KE1_above115", "KE2_mean",
                    "n1_solv"):
            d = float(lq[col]) - float(cub[col])
            deltas.setdefault(col, []).append(d)
            line.append(f"Δ{col} {d:+.4f}")
        print("  " + "  ".join(line))
    print("\npaired mean ± sample SD (n = 5):")
    for col, vals in deltas.items():
        arr = np.asarray(vals)
        print(f"  Δ{col}: {arr.mean():+.4f} ± {arr.std(ddof=1):.4f}")

    # --- Context anchors (display only — no verdict fires here). -----------
    print(f"\nreference n = 1 KED: mean {REF_N1['mean']} / median "
          f"{REF_N1['median']} / mode {REF_N1['mode']} / above-1.15 "
          f"{REF_N1['above115']};  n = 2 mean {REF_N2_MEAN}")
    for pool in (cub_pool, lq_pool):
        print(f"{pool['label']}: KE1 {float(pool['KE1_mean']):.3f} ± SD "
              f"{float(pool['KE1_sd']):.3f} (med {float(pool['KE1_med']):.3f}"
              f", mode {float(pool['KE1_mode']):.3f}, above-1.15 "
              f"{float(pool['KE1_above115']):.3f}, n {int(pool['KE1_n'])}); "
              f"KE2 {float(pool['KE2_mean']):.3f} ± "
              f"{float(pool['KE2_sd']):.3f}")

    if SAVE_CSV_PATH is not None:
        print(f"\nCSV -> {write_rows_csv(SAVE_CSV_PATH, rows)}")


if __name__ == "__main__":
    main()
