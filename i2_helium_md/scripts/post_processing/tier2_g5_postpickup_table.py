"""THESIS_FINAL_GOALS G5 — droplet-size sampler A/B on the lq A5 cell (scorer).

Scores the ``gen_tier2_detfix_postpickup.py --cell lq`` runs against the
committed G4.14 A5 members (raw sampler). Same design as the h405 A/B
(``tier2_detfix_postpickup_table.py``): one field flipped, not CRN, so the
comparison is pooled arm vs pooled arm, s1–s5, N = 5000 each.

Rows:

* A — ``{m}_lqraw`` / ``g5lqpp{m}`` per member, ``pooled_lqraw_5000``,
  ``pooled_lqpp_5000``, ``delta_lqpp_minus_raw``;
* B1 — ``g5lqppe####s1`` per E0 of the column present (the A s1 row is the
  E0 = 0.3155 point); B1' — ``g5lqppt##s1`` per tau (A s1 = tau 6.4);
  B1c — ``g5lqppc###e####s1`` per ridge c (E0 re-placed per c);
* B2 — members, ``pooled_lqppb2_5000`` and ``delta_lqppb2_minus_raw`` once
  the B2 cell is set and all six members exist (N = 6000, vs the A5 pool);
* h405 contrast — ``h405_pooled_raw5000`` / ``h405_pooled_pp5000`` copied
  from ``detfix_postpickup_table.csv`` (the h405 A/B scorer's output).

A member absent from disk is listed as missing (not run, or the run raised —
e.g. the D4 marginal-partner safeguard); ``trap_marg`` tracks marginal ions
per row. Sizes and chords read ``neutral.npz`` (ion.npz is dropped).

Oracles (failure stops before a new number):

1. **A5 anchor** — the six committed A5 members re-pooled here equal the
   committed ``c200t64_pooled6000`` row of ``g414_table.csv`` (4 dp).
2. **Arm identity** — each lq post_pickup cfg differs from its A5 member's in
   exactly ``droplet_size_sampler_mode`` (raw → post_pickup), plus the E0 /
   tau key when re-placed.

Output (commit with ``git add -f``; every column):
``data/runs/h2b_forward_model/g5_postpickup_table.csv`` and
``g5_postpickup_chord.csv`` (chord -> (n, KE) map per lq arm).

Invocation::

    python scripts/post_processing/tier2_g5_postpickup_table.py
"""

from __future__ import annotations

import csv
import dataclasses
import sys
import warnings
from pathlib import Path
from typing import Any, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.postprocess.abundance_loader import (  # noqa: E402
    load_he_abundance_reference,
)
from i2_helium_md.postprocess.ihe_ked import load_ihe_ked_reference  # noqa: E402
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    pool_confirmation_reads,
    read_confirmation_detection,
)
from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2_detfix_battery import MEMBER_SEEDS  # noqa: E402
from scripts.gen_tier2_detfix_postpickup import (  # noqa: E402
    B1C_C, B1_E0_EV, B1_MEMBER, B1_TAU_PS, B2_C, B2_E0_EV, B2_TAU_PS, LAW_DIFF_KEYS, POOL_MEMBERS,
    a5_member_cell, b1c_e0, lq_label, lq_run_dir,
)
from scripts.gen_tier2_g413_ridge_battery import _E0_KEY, _TAU_KEY  # noqa: E402
from scripts.gen_tier2_g414_tauline import run_dir_for as a5_run_dir_for  # noqa: E402
from scripts.post_processing.tier2_detfix_battery_table import (  # noqa: E402
    OracleFailure,
    _check_4dp,
    _csv_row,
    _row,
)
from scripts.post_processing.tier2_detfix_postpickup_table import (  # noqa: E402
    OUT_CSV as H405_AB_CSV,
    chord_rows,
    size_columns,
)
from scripts.post_processing.tier2_g413_table import ANCHOR_COLS  # noqa: E402
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

OUT_CSV = RUNS_ROOT / "h2b_forward_model" / "g5_postpickup_table.csv"
CHORD_CSV = RUNS_ROOT / "h2b_forward_model" / "g5_postpickup_chord.csv"
G414_TABLE_CSV = RUNS_ROOT / "h2b_forward_model" / "g414_table.csv"
A5_POOLED_LABEL = "c200t64_pooled6000"
H405_ROWS = {"pooled_raw5000": "h405_pooled_raw5000", "pooled_pp5000": "h405_pooled_pp5000"}

# ---------------------------------------------------------------------------


def _a5_dir(member: str) -> Path:
    return a5_run_dir_for(a5_member_cell(member))


def _check_arm(member: str, e0_eV: Optional[float], tau_ps: Optional[float],
               c: Optional[float] = None) -> None:
    """Oracle 2: the pp run's cfg vs its A5 member's — sampler (+ E0 / tau /
    ridge law) only."""
    cfg_raw = RunDirectory(_a5_dir(member)).load_cfg()
    cfg_pp = RunDirectory(lq_run_dir(member, e0_eV, tau_ps, c)).load_cfg()
    diff = sorted(f.name for f in dataclasses.fields(cfg_raw)
                  if getattr(cfg_raw, f.name) != getattr(cfg_pp, f.name))
    want = sorted({"droplet_size_sampler_mode"} | ({_E0_KEY} if e0_eV is not None else set())
                  | ({_TAU_KEY} if tau_ps is not None else set())
                  | (set(LAW_DIFF_KEYS) if c is not None else set()))
    modes = (cfg_raw.droplet_size_sampler_mode, cfg_pp.droplet_size_sampler_mode)
    if diff != want or modes != ("raw", "post_pickup"):
        raise OracleFailure(f"[{lq_label(member, e0_eV, tau_ps, c)}] cfg diff {diff} / "
                            f"modes {modes}")


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    ab = load_he_abundance_reference(ABUNDANCE_REFERENCE_CSV)
    ked = load_ihe_ked_reference(IHE_KED_REFERENCE_CSV)

    def scored(label: str, member: str, run_dir: Path) -> tuple[dict[str, Any], Any]:
        cfg = RunDirectory(run_dir).load_cfg()
        read = read_confirmation_detection(load_detection_result(run_dir / "detection.npz"),
                                           label=label)
        row = _row(label, member, cfg.droplet_size_sampler_mode, cfg.seed, read, ab, ked)
        row["E0_eV"] = float(getattr(cfg, _E0_KEY) * cfg.coulomb_available_eV)
        row["tau_ps"] = float(getattr(cfg, _TAU_KEY))
        row["c_amu_per_A"] = float(cfg.drag_coefficients.coefficients["c"])
        row.update(size_columns([run_dir]))
        return row, read

    def pooled(label: str, mode: str, reads, run_dirs) -> dict[str, Any]:
        cfg = RunDirectory(run_dirs[0]).load_cfg()
        row = _row(label, "pooled5000", mode, "-",
                   pool_confirmation_reads(reads, label=label), ab, ked)
        row["E0_eV"] = float(getattr(cfg, _E0_KEY) * cfg.coulomb_available_eV)
        row["tau_ps"] = float(getattr(cfg, _TAU_KEY))
        row.update(size_columns(run_dirs))
        return row

    def delta(label: str, a: dict, b: dict) -> dict[str, Any]:
        out: dict[str, Any] = {"label": label, "member": "pooled5000",
                               "pipeline": "post_pickup - raw", "seed": "-"}
        for key, value in a.items():
            if key not in out:
                try:
                    out[key] = float(value) - float(b[key])
                except (TypeError, ValueError):
                    out[key] = ""
        return out

    # Oracle 1: the committed A5 pool reproduces (all six members).
    a5_reads = {m: scored(f"{m}_lqraw", m, _a5_dir(m)) for m in MEMBER_SEEDS}
    a5_pool = _row("a5_pooled6000", "pooled6000", "raw", "-",
                   pool_confirmation_reads([r for _, r in a5_reads.values()],
                                           label="a5_pooled6000"), ab, ked)
    _check_4dp("[anchor] A5 pooled6000", a5_pool, _csv_row(G414_TABLE_CSV, A5_POOLED_LABEL),
               ANCHOR_COLS)
    print(f"A5 anchor PASSED ({A5_POOLED_LABEL}, 4 dp)", flush=True)

    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    raw_dirs = [_a5_dir(m) for m in POOL_MEMBERS]
    for m in POOL_MEMBERS:
        rows.append(a5_reads[m][0])
    raw_pool = pooled("pooled_lqraw_5000", "raw", [a5_reads[m][1] for m in POOL_MEMBERS],
                      raw_dirs)
    rows.append(raw_pool)
    chord: list[dict[str, Any]] = chord_rows("lq_raw", raw_dirs)
    a5_6000 = pooled("pooled_lqraw_6000", "raw", [r for _, r in a5_reads.values()],
                     [_a5_dir(m) for m in MEMBER_SEEDS])
    raw_ref = {len(POOL_MEMBERS): raw_pool, len(MEMBER_SEEDS): a5_6000}

    def arm(tag: str, e0: Optional[float], tau: Optional[float], members,
            c: Optional[float] = None) -> None:
        reads, dirs = [], []
        for m in members:
            rd = lq_run_dir(m, e0, tau, c)
            if not (rd / "detection.npz").exists():
                missing.append(lq_label(m, e0, tau, c))
                continue
            _check_arm(m, e0, tau, c)
            row, read = scored(lq_label(m, e0, tau, c), m, rd)
            rows.append(row)
            reads.append(read)
            dirs.append(rd)
        if len(members) in raw_ref and len(reads) == len(members):
            size = 1000 * len(members)
            pool = pooled(f"pooled_{tag}_{size}", "post_pickup", reads, dirs)
            rows.append(pool)
            rows.append(delta(f"delta_{tag}_minus_raw_{size}", pool, raw_ref[len(members)]))
            chord.extend(chord_rows(tag, dirs))

    arm("lqpp", None, None, POOL_MEMBERS)                    # A
    # A at N = 6000 (s1-s5 + r6, user 2026-10-08) vs the committed A5 pool.
    r6_dir = lq_run_dir("r6")
    if (r6_dir / "detection.npz").exists():
        _check_arm("r6", None, None)
        rows.append(a5_reads["r6"][0])
        r6_row, r6_read = scored(lq_label("r6"), "r6", r6_dir)
        rows.append(r6_row)
        pp_reads = [r for r in (scored(lq_label(m), m, lq_run_dir(m))[1]
                                for m in POOL_MEMBERS)] + [r6_read]
        pp_6000 = pooled("pooled_lqpp_6000", "post_pickup", pp_reads,
                         [lq_run_dir(m) for m in MEMBER_SEEDS])
        rows += [a5_6000, pp_6000, delta("delta_lqpp_minus_raw_6000", pp_6000, a5_6000)]
    else:
        missing.append(lq_label("r6"))
    for e0 in B1_E0_EV:                                       # B1 (E0 column)
        arm("lqppb1", e0, None, (B1_MEMBER,))
    for tau in B1_TAU_PS:                                     # B1' (tau column)
        arm("lqppb1tau", None, tau, (B1_MEMBER,))
    for c in B1C_C:                                           # B1c (ridge c column)
        arm("lqppb1c", b1c_e0(c), None, (B1_MEMBER,), c)
    if B2_E0_EV is not None or B2_TAU_PS is not None or B2_C is not None:   # B2
        arm("lqppb2", B2_E0_EV, B2_TAU_PS, tuple(MEMBER_SEEDS), B2_C)

    if H405_AB_CSV.exists():                                  # h405 contrast
        for src, label in H405_ROWS.items():
            rows.append({**_csv_row(H405_AB_CSV, src), "label": label})
    else:
        print(f"h405 A/B table {H405_AB_CSV.name} missing: run "
              "tier2_detfix_postpickup_table.py for the contrast rows")

    keys = list(dict.fromkeys(k for r in rows for k in r))
    rows = [{k: r.get(k, "") for k in keys} for r in rows]
    show = ("label", "pipeline", "c_amu_per_A", "E0_eV", "tau_ps", "num_scored", "trap", "trap_marg", "supp",
            "nbar_det", "n1_solv", "w1_solv", "midHot", "deepKE", "chi2_med",
            "KE1_mean", "S", "N_mean")
    show = tuple(c for c in show if c in keys)
    print("\n=== G5: sampler A/B on the lq A5 cell — observable vector ===")
    print(format_table([{c: r[c] for c in show} for r in rows]))
    if missing:
        print(f"\nmissing (not run, or the run raised — e.g. D4): {missing}")
    print(f"\nCSV -> {write_rows_csv(OUT_CSV, rows)}")
    print(f"CSV -> {write_rows_csv(CHORD_CSV, chord)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
