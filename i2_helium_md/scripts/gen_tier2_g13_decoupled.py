"""THESIS_FINAL_GOALS G13 — re-run with the self-bound gate lifted at decoupling.

``self_bound_gate_release="decoupled"``: the self-bound gate (argued from
droplet confinement) is lifted wherever the ion is outside the nominal
droplet surface (depth > 0, memoryless) -- in Stage I and in the detection
stage. Full pipeline per member (neutral -> ion to t_h = 500 ps ->
detection, E2 skipped, partner-aware closure), **ion.npz kept**, for

* ``B2``   -- thesis production law (``gen_tier2_detfix_postpickup --b2``),
* ``B3``   -- effective high-tail variant (``--b3``),
* ``h405`` -- the re-landed contrast law (``gen_tier2_h405_pp35.SCAN_POOL_CELL``),

six members each (s1-s6; s6 = the committed member r6). Each member cfg is
the law's committed builder output (its pre-registered diff guard passes
first) with only ``self_bound_gate_release`` changed -- the thesis production
builders of the two generators (``build_thesis_lq`` / ``build_thesis``). Runs
go to ``data/runs/<law>_N1000_s<k>``; every other detail is in each run's
cfg.json. The committed ``never`` batteries are not touched.

Modes::

    python scripts/gen_tier2_g13_decoupled.py --oracle        # B2 s1 under "never"
    python scripts/gen_tier2_g13_decoupled.py --run B2 B3 h405
    python scripts/gen_tier2_g13_decoupled.py --pool

``--oracle`` (first): re-runs B2 s1 under ``never`` into
``B2_N1000_s1_never_oracle`` and requires its ``detection.npz`` to equal the
committed member's array for array -- the pipeline is deterministic, so every
later difference is the release alone. Prints the wall time.

``--pool``: pair-preserving pools (``<law>_N6000_pooled``), scored with the v2
read path (``tier2_g12_ihekedv2_table._score``) next to the intact pools under
the same scorer; the intact B2 / B3 pools must reproduce their committed v2
rows to 4 dp. Output (commit with ``git add -f``):
``data/runs/h2b_forward_model/g13_decoupled_table.csv``.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Callable

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import load_he_abundance_reference  # noqa: E402
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import read_confirmation_detection  # noqa: E402
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    DetectionResult,
    load_detection_result,
    save_detection_result,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.build_pooled_detection_container import TARGETS, pool_members  # noqa: E402
from scripts.gen_tier2_detfix_battery import MEMBER_SEEDS  # noqa: E402
from scripts.gen_tier2_detfix_postpickup import (  # noqa: E402
    THESIS_GATE_RELEASE,
    THESIS_LQ_JOBS,
    build_thesis_lq,
    lq_run_dir,
    run_thesis_lq,
    run_thesis_member,
    thesis_pool_dir,
    thesis_run_dir,
    verify_lq_member,
)
from scripts.gen_tier2_h405_pp35 import run_thesis as run_thesis_h405  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table,
    write_rows_csv,
)
from scripts.post_processing.tier2_g12_ihekedv2_table import (  # noqa: E402
    OUT_CSV as V2_TABLE_CSV,
    _score,
)
from scripts.post_processing.tier2_g5_winners_table import LAWS  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import (  # noqa: E402
    ABUNDANCE_REFERENCE_CSV,
    IHE_KED_REFERENCE_CSV,
    RUNS_ROOT,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

ORACLE_LAW, ORACLE_MEMBER = "B2", "s1"
MAX_CONCURRENCY = 2
OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "g13_decoupled_table.csv"

# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Law:
    """One law: its thesis production run (guard + MD), the committed intact
    pool, and its committed v2 row (None for h405)."""

    run: Callable[[str], str]
    intact_container: str
    v2_label: str | None


LAW_SPECS: dict[str, Law] = {
    "B2": Law(lambda m: run_thesis_lq(("B2", m)),
              LAWS["B2_pp_c160"][0], "B2_pp_c160"),
    "B3": Law(lambda m: run_thesis_lq(("B3", m)),
              LAWS["B3_pp_c160_w0482_a198"][0], "B3_pp_c160_w0482_a198"),
    "h405": Law(run_thesis_h405, TARGETS["h405_pp35_relanded"][2], None),
}


def oracle_run_dir() -> Path:
    """``data/runs/B2_N1000_s1_never_oracle``."""
    d = thesis_run_dir(ORACLE_LAW, ORACLE_MEMBER)
    return d.with_name(d.name + "_never_oracle")


def run_member(job: tuple[str, str]) -> str:
    """Thesis production member (guard, MD, ion.npz kept). ``job`` = (law, member)."""
    law_name, member = job
    return LAW_SPECS[law_name].run(member)


def _equal(a: DetectionResult, b: DetectionResult) -> list[str]:
    return [f.name for f in dataclasses.fields(DetectionResult)
            if not np.array_equal(np.asarray(getattr(a, f.name)),
                                  np.asarray(getattr(b, f.name)))]


def oracle() -> int:
    """B2 s1 re-run under "never" must reproduce the committed detection.npz."""
    t0 = time.perf_counter()
    job = THESIS_LQ_JOBS[ORACLE_LAW]
    verify_lq_member(ORACLE_MEMBER, *job)               # committed pre-registration
    cfg = dataclasses.replace(build_thesis_lq(ORACLE_LAW, ORACLE_MEMBER),
                              self_bound_gate_release="never")
    cfg.validate()
    run_thesis_member(cfg, oracle_run_dir(), f"{ORACLE_LAW} {ORACLE_MEMBER} never")
    committed = load_detection_result(lq_run_dir(ORACLE_MEMBER, *job) / "detection.npz")
    rerun = load_detection_result(oracle_run_dir() / "detection.npz")
    bad = _equal(rerun, committed)
    if bad:
        print(f"ORACLE FAILED: re-run differs from the committed detection.npz in {bad}")
        return 1
    print(f"oracle PASSED: {ORACLE_LAW} {ORACLE_MEMBER} re-run under 'never' reproduces "
          f"the committed detection.npz array for array ({time.perf_counter() - t0:.0f} s)")
    return 0


def run(laws: list[str], concurrency: int) -> int:
    jobs = [(law, m) for law in laws for m in MEMBER_SEEDS]
    workers = max(1, min(concurrency, len(jobs)))
    if workers == 1:
        for job in jobs:
            run_member(job)
        return 0
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=workers) as poolx:
        futures = {poolx.submit(run_member, j): j for j in jobs}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{futures[fut]}] FAILED: {exc!r}", flush=True)
    print(f"{len(jobs) - failures}/{len(jobs)} runs complete.", flush=True)
    return 1 if failures else 0


def _committed_v2_row(label: str) -> dict[str, str]:
    with V2_TABLE_CSV.open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["label"] == label]
    if len(rows) != 1:
        raise RuntimeError(f"expected one {label!r} row in {V2_TABLE_CSV.name}")
    return rows[0]


def _check_4dp(tag: str, got: dict[str, Any], want: dict[str, str]) -> None:
    bad = []
    for k, w in want.items():
        try:
            wf = float(w)
        except (TypeError, ValueError):
            continue
        if k in got and round(float(got[k]), 4) != round(wf, 4):
            bad.append(f"{k}: {float(got[k]):.4f} vs {wf:.4f}")
    if bad:
        raise RuntimeError(f"{tag}: " + "; ".join(bad))


def pool() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)
    rows: list[dict[str, Any]] = []
    for law_name, law in LAW_SPECS.items():
        intact_dir = RUNS_ROOT / law.intact_container
        intact_cfg = RunDirectory(intact_dir).load_cfg()
        label = law.v2_label or "h405pp35_relanded"
        role = LAWS[law.v2_label][1] if law.v2_label else "contrast law (re-landed)"
        intact_row = _score(label, read_confirmation_detection(
            load_detection_result(intact_dir / "detection.npz"), label=label),
            intact_cfg, law.intact_container, role, ab, ked)
        if law.v2_label:
            _check_4dp(f"intact {law_name} vs committed v2 row", intact_row,
                       _committed_v2_row(law.v2_label))
        dirs = [thesis_run_dir(law_name, m) for m in MEMBER_SEEDS]
        members = [load_detection_result(d / "detection.npz") for d in dirs]
        pooled = pool_members(members)
        out_dir = thesis_pool_dir(law_name)
        container = out_dir.name
        cfg = dataclasses.replace(intact_cfg, self_bound_gate_release=THESIS_GATE_RELEASE)
        RunDirectory(out_dir).save_cfg(cfg)
        save_detection_result(pooled, out_dir / "detection.npz")
        (out_dir / "README_POOLED.txt").write_text(
            "FIGURES CONTAINER, not an MD run.\n\n"
            f"{law_name} with self_bound_gate_release='decoupled' (THESIS_FINAL_GOALS "
            "G13): s1-s6 (ion.npz kept), pooled pair-preserving "
            "(build_pooled_detection_container.pool_members):\n"
            + "".join(f"  {d.name}\n" for d in dirs)
            + "Built by scripts/gen_tier2_g13_decoupled.py --pool.\n")
        dec_row = _score(law_name, read_confirmation_detection(pooled, label=label),
                         cfg, container, role, ab, ked)
        rows += [{**intact_row, "self_bound_gate_release": "never"},
                 {**dec_row, "self_bound_gate_release": THESIS_GATE_RELEASE}]
        print(f"[{law_name}] pooled {len(members)} members -> {container}", flush=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    rows = [{k: r.get(k, "") for k in keys} for r in rows]
    show = ("label", "num_scored", "trap", "supp", "nbar_det", "n1_solv", "w1_solv",
            "midHot", "KE1_mean", "deepKE", "chi2_med", "S")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    print(f"\nCSV -> {write_rows_csv(OUT_CSV, rows)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--oracle", action="store_true")
    mode.add_argument("--run", nargs="+", choices=tuple(LAW_SPECS))
    mode.add_argument("--pool", action="store_true")
    parser.add_argument("--concurrency", type=int, default=MAX_CONCURRENCY)
    args = parser.parse_args(argv)
    if args.oracle:
        return oracle()
    if args.run:
        return run(args.run, args.concurrency)
    return pool()


if __name__ == "__main__":
    raise SystemExit(main())
