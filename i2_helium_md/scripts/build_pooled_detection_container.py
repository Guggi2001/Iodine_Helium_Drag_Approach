"""Pooled detection figures-container builder (pair-preserving block layout).

Builds an N = 5000 pooled `detection.npz` figures container from five
battery members, preserving the `(i, i + num_molecules)` fragment-pairing
convention the covariance recipes require (`pair_correlation.py` block
layout): per-ion arrays are concatenated as
``[all fragment-1 blocks | all fragment-2 blocks]`` and the ragged event
CSR is permuted per ion accordingly. Sequential concatenation silently
factorizes the angular pair covariance — the 2026-07-22 user-caught bug
recorded in `drag_migration_log_tier2.md` ("Pooled detection_summary
figures generated"); this committed builder replaces the session-local
one so the layout lesson is enforced by code, not memory.

Every invocation FIRST rebuilds the committed cubic container
(`bigc1v725pooled`) in memory from its five members and requires exact
array-for-array equality with the on-disk container (the builder oracle)
— only then does it build the target container. A built-in pairing
spot-check additionally verifies sampled molecules' fragment rows and
event slices against their source members.

The output is a FIGURES CONTAINER, not an MD run: `cfg.json` is member
s1's cfg with ``num_molecules`` set to the pooled count (the seed field
necessarily shows s1's seed; all member seeds are recorded in
`README_POOLED.txt`).

Invocation (target = the §6.7 lq battery pool):

    python scripts/build_pooled_detection_container.py
"""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    DetectionResult,
    load_detection_result,
    save_detection_result,
)
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.gen_tier2atlas_lqbattery import (  # noqa: E402
    BATTERY_MATRIX,
    atlas_run_dir_name,
)

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

RUNS_ROOT = PROJECT_ROOT / "data" / "runs"

# Builder oracle: these five members must rebuild this committed container
# exactly (array-for-array) before any target is built.
ORACLE_MEMBERS = [
    f"9A_drag_shared_pure_cubic_N1000_tier2probe_conf270_bigc1v725s{k}"
    for k in (1, 2, 3, 4, 5)
]
ORACLE_CONTAINER = (
    "9A_drag_shared_pure_cubic_N5000_tier2probe_conf270_bigc1v725pooled"
)

# Target: the §6.7 lq battery pool (atlas namespace: no `_tier2_`, no
# `tier2probe` substring).
TARGET_MEMBERS = [atlas_run_dir_name(spec.label) for spec in BATTERY_MATRIX]
TARGET_SEEDS = [spec.seed for spec in BATTERY_MATRIX]
TARGET_CONTAINER = "9A_drag_shared_lq_N5000_tier2atlas_conf270_qccbigpooled"

PAIR_SPOTCHECK_MOLECULES = 60      # sampled molecules for the pairing oracle
PAIR_SPOTCHECK_SEED = 20260730

SCALAR_FIELDS = ("num_molecules", "t_handover_ps", "detection_time_ps")


def _field_groups() -> tuple[list[str], list[str]]:
    """(per-ion field names, flat event field names) from the dataclass."""
    per_ion, events = [], []
    for f in dataclasses.fields(DetectionResult):
        if f.name in SCALAR_FIELDS or f.name == "event_offsets":
            continue
        (events if f.name.startswith("event_") else per_ion).append(f.name)
    return per_ion, events


def pool_members(members: list[DetectionResult]) -> DetectionResult:
    """Pair-preserving pool of equal-N members (block layout, CSR rebuild)."""
    n = int(members[0].num_molecules)
    if any(int(m.num_molecules) != n for m in members):
        raise ValueError("members disagree on num_molecules")
    for name in ("t_handover_ps", "detection_time_ps"):
        vals = {float(getattr(m, name)) for m in members}
        if len(vals) != 1:
            raise ValueError(f"members disagree on scalar {name}: {vals}")

    per_ion_fields, event_fields = _field_groups()
    pooled: dict[str, object] = {
        "num_molecules": n * len(members),
        "t_handover_ps": float(members[0].t_handover_ps),
        "detection_time_ps": float(members[0].detection_time_ps),
    }
    for name in per_ion_fields:
        blocks_1 = [np.asarray(getattr(m, name))[:n] for m in members]
        blocks_2 = [np.asarray(getattr(m, name))[n:] for m in members]
        pooled[name] = np.concatenate(blocks_1 + blocks_2)

    # Event CSR: gather each ion's slice in the pooled ion order.
    order = [(mi, i) for mi in range(len(members)) for i in range(n)]
    order += [(mi, i + n) for mi in range(len(members)) for i in range(n)]
    counts = np.empty(len(order), dtype=np.int64)
    slices: list[tuple[int, int, int]] = []
    for row, (mi, i) in enumerate(order):
        off = members[mi].event_offsets
        s, e = int(off[i]), int(off[i + 1])
        counts[row] = e - s
        slices.append((mi, s, e))
    offsets = np.zeros(len(order) + 1,
                       dtype=np.asarray(members[0].event_offsets).dtype)
    offsets[1:] = np.cumsum(counts)
    pooled["event_offsets"] = offsets
    for name in event_fields:
        flats = [np.asarray(getattr(members[mi], name))[s:e]
                 for mi, s, e in slices]
        pooled[name] = np.concatenate(flats)

    result = DetectionResult(**pooled)  # type: ignore[arg-type]

    # Pairing oracle: sampled molecules' rows + event slices must equal the
    # source member's rows for BOTH fragments at the block-layout indices.
    rng = np.random.default_rng(PAIR_SPOTCHECK_SEED)
    total = n * len(members)
    for mol in rng.choice(total, size=min(PAIR_SPOTCHECK_MOLECULES, total),
                          replace=False):
        mi, j = divmod(int(mol), n)
        for frag, (pool_i, member_i) in enumerate(
                (((int(mol)), j), (int(mol) + total, j + n))):
            m = members[mi]
            for name in per_ion_fields:
                if not np.array_equal(
                        np.asarray(getattr(result, name))[pool_i],
                        np.asarray(getattr(m, name))[member_i]):
                    raise AssertionError(
                        f"pairing oracle: {name} mismatch at pooled ion "
                        f"{pool_i} (member {mi} ion {member_i})")
            ps, pe = int(result.event_offsets[pool_i]), int(
                result.event_offsets[pool_i + 1])
            ms, me = int(m.event_offsets[member_i]), int(
                m.event_offsets[member_i + 1])
            for name in event_fields:
                if not np.array_equal(
                        np.asarray(getattr(result, name))[ps:pe],
                        np.asarray(getattr(m, name))[ms:me]):
                    raise AssertionError(
                        f"pairing oracle: event {name} mismatch at pooled "
                        f"ion {pool_i}")
    return result


def _load_members(names: list[str]) -> list[DetectionResult] | None:
    missing = [nm for nm in names
               if not (RUNS_ROOT / nm / "detection.npz").exists()]
    if missing:
        for nm in missing:
            print(f"    missing {nm}")
        return None
    return [load_detection_result(RUNS_ROOT / nm / "detection.npz")
            for nm in names]


def run_builder_oracle() -> bool:
    """Rebuild the committed cubic container; require exact equality."""
    members = _load_members(ORACLE_MEMBERS)
    if members is None:
        print("*** BUILDER ORACLE cannot run (cubic members missing) ***")
        return False
    rebuilt = pool_members(members)
    committed = load_detection_result(
        RUNS_ROOT / ORACLE_CONTAINER / "detection.npz")
    bad = []
    for f in dataclasses.fields(DetectionResult):
        got = np.asarray(getattr(rebuilt, f.name))
        want = np.asarray(getattr(committed, f.name))
        if got.shape != want.shape or not np.array_equal(got, want):
            bad.append(f.name)
    if bad:
        print("*** BUILDER ORACLE FAILED — rebuilt cubic container differs "
              f"from the committed one on: {bad} ***")
        return False
    print("builder oracle OK (committed bigc1v725pooled rebuilt "
          "array-for-array from its five members).")
    return True


def main() -> None:
    if not run_builder_oracle():
        return

    print("target members:")
    members = _load_members(TARGET_MEMBERS)
    if members is None:
        print("target battery incomplete — container not built.")
        return
    pooled = pool_members(members)

    target_dir = RUNS_ROOT / TARGET_CONTAINER
    if (target_dir / "detection.npz").exists():
        print(f"target container already exists: {target_dir} — not "
              "overwriting (delete it to rebuild).")
        return

    run = RunDirectory(target_dir)
    cfg = RunDirectory(RUNS_ROOT / TARGET_MEMBERS[0]).load_cfg()
    cfg = dataclasses.replace(cfg, num_molecules=int(pooled.num_molecules))
    run.save_cfg(cfg)
    save_detection_result(pooled, target_dir / "detection.npz")
    (target_dir / "README_POOLED.txt").write_text(
        "FIGURES CONTAINER, not an MD run.\n\n"
        f"Pooled from the 5 x N=1000 members (seeds {TARGET_SEEDS}):\n"
        + "".join(f"  {nm}\n" for nm in TARGET_MEMBERS)
        + "\ncfg.json is member 1's cfg with num_molecules set to the "
        "pooled count\n(the seed field shows only member 1's seed).\n\n"
        "LAYOUT REQUIREMENT (do not rebuild sequentially): per-ion arrays\n"
        "are [all fragment-1 blocks | all fragment-2 blocks] so the\n"
        "(i, i + num_molecules) fragment-pairing convention of the\n"
        "covariance recipes holds; the event CSR is permuted per ion.\n"
        "Sequential concatenation factorizes the angular pair covariance\n"
        "(2026-07-22 lesson, drag_migration_log_tier2.md). Built by the\n"
        "committed scripts/build_pooled_detection_container.py (builder\n"
        "oracle: exact rebuild of bigc1v725pooled).\n",
        encoding="utf-8",
    )
    print(f"pooled container written -> {target_dir} "
          f"({int(pooled.num_molecules)} molecules, "
          f"{pooled.n_detected.size} fragments, "
          f"{pooled.event_time_ps.size} events)")


if __name__ == "__main__":
    main()
