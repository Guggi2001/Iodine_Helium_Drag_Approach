"""Run batteries of drag-model members and pool them.

One member = neutral stage -> ion stage to the fixed handover ``T_HANDOVER_PS``
-> detection stage, in its own run directory. A run directory is run when it
holds no artifact, skipped when it is complete and holds the same model, and
refused otherwise; nothing here overwrites a file.
"""

from __future__ import annotations

import dataclasses
import math
import zipfile
from collections.abc import Iterator, Sequence
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from concurrent.futures.process import BrokenProcessPool
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from ..config import SimConfig
from ..presets import same_model
from .detection_stage import (
    DetectionResult,
    load_detection_result,
    run_detection_stage,
    save_detection_result,
)
from .ion import (
    _NUM_2N_T_ARRAYS_ION,
    _decide_stride_ion,
    _internal_step_count_ion,
    run_ion_propagation,
)
from .neutral import run_neutral_propagation
from .run_directory import RunDirectory
from .stage1_events import STAGE1_EVENTS_FILENAME

#: End of the ion stage, where the detection stage takes over [ps].
T_HANDOVER_PS = 500.0
#: Thesis storage stride; it divides the last step id 50000, so the last stored
#: column is the final state at ``T_HANDOVER_PS``.
STORAGE_STRIDE = 16
#: Ion-checkpoint budget [bytes]: 3126 stored columns at N = 1000 give that stride.
MAX_BYTES_ION = _NUM_2N_T_ARRAYS_ION * 2 * 1000 * 8 * 3126

#: Files of a complete pool directory.
POOL_ARTIFACTS = ("cfg.json", "detection.npz", "README_POOLED.txt")

_SCALAR_FIELDS = ("num_molecules", "t_handover_ps", "detection_time_ps")
_PAIR_SPOTCHECK_MOLECULES = 60
_PAIR_SPOTCHECK_SEED = 20260730

#: Files of a complete member directory, in the order they are written.
ARTIFACTS = ("cfg.json", "neutral.npz", "ion.npz", STAGE1_EVENTS_FILENAME, "detection.npz")


#: Present while a process writes the directory; a crash leaves it behind.
LOCK_NAME = ".running"


class RunDirectoryMismatch(RuntimeError):
    """A run directory is partial, damaged, in use or holds a different model."""


def _claimed_by_another(run_dir: Path) -> RunDirectoryMismatch:
    return RunDirectoryMismatch(
        f"{run_dir}: {LOCK_NAME} is present. Another process is writing it, or a run "
        "crashed; if none is running, move the directory away.")


def _load_cfg(run_dir: Path) -> SimConfig:
    try:
        return RunDirectory(run_dir).load_cfg()
    except (ValueError, TypeError, KeyError) as exc:
        raise RunDirectoryMismatch(f"{run_dir}: cannot read cfg.json: {exc}") from exc


def is_finished(run_dir: Path, cfg: SimConfig, artifacts: Sequence[str] = ARTIFACTS) -> bool:
    """Whether ``run_dir`` already holds the finished run of ``cfg``.

    Parameters
    ----------
    run_dir : Path
        Directory to inspect; entries not named in ``artifacts`` are ignored.
    cfg : SimConfig
        Model the directory has to hold.
    artifacts : sequence of str
        File names of a complete run; must contain ``cfg.json``.

    Returns
    -------
    bool
        False if none of the artifacts exists, True if all do and the model is
        the same (:func:`i2_helium_md.presets.same_model`).

    Raises
    ------
    RunDirectoryMismatch
        If another process is writing the directory, some artifacts are missing
        or unreadable (named), or the model differs (fields named).
    """
    if (run_dir / LOCK_NAME).exists():
        raise _claimed_by_another(run_dir)
    return _holds_finished_run(run_dir, cfg, artifacts)


def _holds_finished_run(run_dir: Path, cfg: SimConfig, artifacts: Sequence[str]) -> bool:
    missing = [name for name in artifacts if not (run_dir / name).is_file()]
    if len(missing) == len(artifacts):
        return False
    if missing:
        raise RunDirectoryMismatch(
            f"{run_dir}: incomplete run, missing {missing}. Move the directory away "
            "to run it again.")
    unreadable = [name for name in artifacts
                  if name.endswith(".npz") and not zipfile.is_zipfile(run_dir / name)]
    if unreadable:
        raise RunDirectoryMismatch(
            f"{run_dir}: damaged run, truncated {unreadable}. Move the directory away "
            "to run it again.")
    stored = _load_cfg(run_dir)
    diff = same_model(stored, cfg)
    if diff:
        raise RunDirectoryMismatch(
            f"{run_dir}: holds a different model; stored vs requested: " + "; ".join(
                f"{name} {getattr(stored, name)!r} vs {getattr(cfg, name)!r}" for name in diff))
    return True


@contextmanager
def claimed(run_dir: Path, cfg: SimConfig, artifacts: Sequence[str]) -> Iterator[bool]:
    """Hold ``run_dir`` for writing, so that no second process starts in it.

    Parameters
    ----------
    run_dir, cfg, artifacts
        As for :func:`is_finished`.

    Yields
    ------
    bool
        False if another process finished this run before the claim: do not write.

    Raises
    ------
    RunDirectoryMismatch
        As for :func:`is_finished`.
    """
    run_dir.mkdir(parents=True, exist_ok=True)
    lock = run_dir / LOCK_NAME
    try:
        lock.touch(exist_ok=False)
    except FileExistsError:
        raise _claimed_by_another(run_dir) from None
    try:
        yield not _holds_finished_run(run_dir, cfg, artifacts)
    finally:
        lock.unlink()


def verify_handover_alignment(cfg: SimConfig) -> tuple[int, int]:
    """Check that the last stored ion column is the final state at ``T_HANDOVER_PS``.

    Parameters
    ----------
    cfg : SimConfig
        Reads the ion time grid and ``num_molecules``.

    Returns
    -------
    tuple of int
        Internal step count and storage stride.

    Raises
    ------
    ValueError
        If the last step is not at the handover, or the stride is not
        ``STORAGE_STRIDE`` or does not divide the last step id.
    """
    steps = _internal_step_count_ion(cfg)
    last_id = steps - 1
    if not math.isclose(last_id * cfg.dt_ion, T_HANDOVER_PS, abs_tol=1e-9):
        raise ValueError(
            f"last ion step {last_id} is at {last_id * cfg.dt_ion} ps, not at the "
            f"handover {T_HANDOVER_PS} ps")
    stride, _ = _decide_stride_ion(cfg.num_molecules, steps, MAX_BYTES_ION)
    if stride != STORAGE_STRIDE or last_id % stride:
        raise ValueError(
            f"storage stride {stride} (expected {STORAGE_STRIDE}) must divide the last "
            f"step {last_id}: the detection stage starts from the last stored column")
    return steps, stride


def run_member(cfg: SimConfig, run_dir: Path, *, smoke: bool = False) -> Path:
    """Run one member into ``run_dir``, or skip it if it is already there.

    Parameters
    ----------
    cfg : SimConfig
        Member config, e.g. a thesis preset.
    run_dir : Path
        Target directory; see :func:`is_finished` for skip and refusal.
    smoke : bool
        Skip the 500 ps handover checks (short test runs).

    Returns
    -------
    Path
        ``run_dir``, holding the files of ``ARTIFACTS``.

    Raises
    ------
    RunDirectoryMismatch
        If the directory is partial or holds a different model.
    ValueError
        If the handover is not aligned (not with ``smoke``).
    """
    run_dir = Path(run_dir)
    if is_finished(run_dir, cfg):
        print(f"[{run_dir.name}] skip (complete)", flush=True)
        return run_dir
    cfg.validate()
    if not smoke:
        verify_handover_alignment(cfg)
    with claimed(run_dir, cfg, ARTIFACTS) as start:
        if start:
            print(f"[{run_dir.name}] running, seed {cfg.seed}", flush=True)
            _run_stages(cfg, run_dir, smoke=smoke)
    return run_dir


def _run_stages(cfg: SimConfig, run_dir: Path, *, smoke: bool) -> None:
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION, record_events=True)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not smoke and not math.isclose(t_last, T_HANDOVER_PS, abs_tol=1e-9):
        raise RuntimeError(f"{run_dir}: last stored ion column at {t_last} ps, not at the "
                           f"handover {T_HANDOVER_PS} ps")
    detection = run_detection_stage(ion, cfg, save_path=run_dir / "detection.npz")
    detected = detection.detected_mask
    print(f"[{run_dir.name}] done: scored {int(detected.sum())}/{detected.size}", flush=True)


def run_battery(
    cfgs: dict[Path, SimConfig], *, concurrency: int = 2, smoke: bool = False,
) -> list[Path]:
    """Run several members, ``concurrency`` at a time.

    After a failure the started members finish and no new one starts. If a worker
    process is killed, the members running at that moment are left incomplete.

    Parameters
    ----------
    cfgs : dict of Path to SimConfig
        Member directory -> config.
    concurrency : int
        Processes; 1 runs in this process.
    smoke : bool
        See :func:`run_member`.

    Returns
    -------
    list of Path
        The member directories, in the order given.

    Raises
    ------
    ValueError
        If ``concurrency`` is below 1.
    RuntimeError
        Naming every failed member; the first failure is the ``__cause__``.
    """
    if concurrency < 1:
        raise ValueError(f"concurrency must be at least 1, got {concurrency}")
    failures: list[tuple[Path, BaseException]] = []
    queue = list(cfgs.items())
    if concurrency == 1:
        for run_dir, cfg in queue:
            try:
                run_member(cfg, run_dir, smoke=smoke)
            except Exception as exc:  # noqa: BLE001 - reported below with its directory
                failures.append((run_dir, exc))
                break
    else:
        # Submitted one by one: a pool prefetches queued jobs, which would then
        # start after a failure.
        with ProcessPoolExecutor(max_workers=concurrency) as pool:
            running: dict[Future, Path] = {}
            while running or (queue and not failures):
                while queue and not failures and len(running) < concurrency:
                    run_dir, cfg = queue.pop(0)
                    try:
                        running[pool.submit(run_member, cfg, run_dir, smoke=smoke)] = run_dir
                    except BrokenProcessPool as exc:
                        failures.append((run_dir, exc))
                done, _ = wait(running, return_when=FIRST_COMPLETED)
                for future in done:
                    run_dir = running.pop(future)
                    if future.exception() is not None:
                        failures.append((run_dir, future.exception()))
    if failures:
        raise RuntimeError("failed: " + "; ".join(
            f"{run_dir.name}: " + (
                "a worker process died (killed or out of memory); move the incomplete "
                "directories away" if isinstance(exc, BrokenProcessPool) else repr(exc))
            for run_dir, exc in failures)) from failures[0][1]
    return list(cfgs)


def _field_groups() -> tuple[list[str], list[str]]:
    """(per-ion, flat event) field names of ``DetectionResult``."""
    per_ion, events = [], []
    for f in dataclasses.fields(DetectionResult):
        if f.name in _SCALAR_FIELDS or f.name == "event_offsets":
            continue
        (events if f.name.startswith("event_") else per_ion).append(f.name)
    return per_ion, events


def pool_members(members: list[DetectionResult]) -> DetectionResult:
    """Pool equal-N members, keeping fragment pairs together.

    Per-ion arrays become ``[all fragment-1 blocks | all fragment-2 blocks]``,
    so ``(i, i + num_molecules)`` stays a pair; the ragged event record is
    permuted per ion. Plain concatenation would break the angular pair
    covariance. A sample of molecules is checked against its source members.

    Parameters
    ----------
    members : list of DetectionResult
        At least one; same ``num_molecules`` and times.

    Returns
    -------
    DetectionResult
        ``num_molecules`` summed over the members.

    Raises
    ------
    ValueError
        If there are no members or they disagree on size or times.
    """
    if not members:
        raise ValueError("no members to pool")
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
    rng = np.random.default_rng(_PAIR_SPOTCHECK_SEED)
    total = n * len(members)
    for mol in rng.choice(total, size=min(_PAIR_SPOTCHECK_MOLECULES, total),
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


def pool_battery(run_dirs: list[Path], out_dir: Path) -> Path:
    """Pool finished members into ``out_dir``, or leave an equal pool as it is.

    Parameters
    ----------
    run_dirs : list of Path
        Finished member directories of one model, distinct seeds. Order matters:
        it is the molecule order of the pool.
    out_dir : Path
        Pool directory; gets the files of ``POOL_ARTIFACTS``. ``cfg.json`` is the
        first member's with ``num_molecules`` summed.

    Returns
    -------
    Path
        ``out_dir``.

    Raises
    ------
    FileNotFoundError
        If a member has no ``detection.npz``.
    ValueError
        If there are no members, or two share a seed, or one is unseeded.
    RunDirectoryMismatch
        If members differ in more than the seed, a ``detection.npz`` is truncated,
        of another size than its config or the same as another member's, or
        ``out_dir`` holds another pool.
    """
    run_dirs, out_dir = [Path(d) for d in run_dirs], Path(out_dir)
    if not run_dirs:
        raise ValueError("no members to pool")
    for run_dir in run_dirs:
        if not (run_dir / "detection.npz").is_file():
            raise FileNotFoundError(f"{run_dir} has no detection.npz")
        if not zipfile.is_zipfile(run_dir / "detection.npz"):
            raise RunDirectoryMismatch(f"{run_dir}: detection.npz is truncated")
    cfgs = [_load_cfg(d) for d in run_dirs]
    for run_dir, cfg in zip(run_dirs, cfgs):
        diff = [name for name in same_model(cfgs[0], cfg) if name != "seed"]
        if diff:
            raise RunDirectoryMismatch(
                f"{run_dir} and {run_dirs[0].name} differ in more than the seed: {diff}")
    seeds = [cfg.seed for cfg in cfgs]
    bad_seeds = sorted({str(s) for s in seeds if s is None or seeds.count(s) > 1})
    if bad_seeds:
        raise ValueError(f"members to pool need distinct seeds; repeated or unset: {bad_seeds}")

    members = [load_detection_result(d / "detection.npz") for d in run_dirs]
    for k, (run_dir, cfg, member) in enumerate(zip(run_dirs, cfgs, members)):
        if member.num_molecules != cfg.num_molecules:
            raise RunDirectoryMismatch(
                f"{run_dir}: detection.npz has num_molecules {member.num_molecules}, "
                f"cfg.json {cfg.num_molecules}")
        for other_dir, other in zip(run_dirs[k + 1:], members[k + 1:]):
            if np.array_equal(member.vx_detected, other.vx_detected):
                raise RunDirectoryMismatch(
                    f"{run_dir} and {other_dir} hold the same ions under two seeds")
    pooled = pool_members(members)
    cfg = dataclasses.replace(cfgs[0], num_molecules=pooled.num_molecules)
    if is_finished(out_dir, cfg, POOL_ARTIFACTS):
        stored = load_detection_result(out_dir / "detection.npz")
        diff = [f.name for f in dataclasses.fields(DetectionResult)
                if not np.array_equal(getattr(pooled, f.name), getattr(stored, f.name))]
        if diff:
            raise RunDirectoryMismatch(f"{out_dir}: holds a different pool; differs in {diff}")
        print(f"[{out_dir.name}] skip (complete)", flush=True)
        return out_dir
    with claimed(out_dir, cfg, POOL_ARTIFACTS) as start:
        if not start:
            raise RunDirectoryMismatch(f"{out_dir}: another process pooled it meanwhile")
        RunDirectory(out_dir).save_cfg(cfg)
        save_detection_result(pooled, out_dir / "detection.npz")
        (out_dir / "README_POOLED.txt").write_text(
            "Pooled detection results, not an MD run.\n\nMembers:\n"
            + "".join(f"  {d.name} (seed {seed})\n" for d, seed in zip(run_dirs, seeds))
            + "\ncfg.json is the first member's with num_molecules summed. Per-ion arrays\n"
            "are [all fragment-1 blocks | all fragment-2 blocks], so (i, i + num_molecules)\n"
            "is a fragment pair.\n",
            encoding="utf-8",
        )
    print(f"[{out_dir.name}] pooled {len(run_dirs)} members", flush=True)
    return out_dir
