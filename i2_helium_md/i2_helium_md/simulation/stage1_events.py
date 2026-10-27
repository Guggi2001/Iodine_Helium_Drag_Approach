"""Stage-I (ion MD) event log -- ``stage1_events.npz`` (THESIS_DISCUSSION §16).

A dedicated sidecar artifact next to ``ion.npz``, **not** an IonCheckpoint
field (the checkpoint schema is unchanged). It records, per biphasic Stage-I
run, every discrete occupancy event and every crossing of the nominal droplet
surface at the internal step ``dt_ion``, independent of the storage stride:

* **Events** (one row per shed or pickup): ion index, internal step, time,
  kind, occupancy before / after, and the pre-step position the channels read
  (the gate / density depth of :func:`biphasic_step_fired`). A same-step
  shed + pickup on one ion gives two rows, shed first (the frozen draw order),
  so ``n_after`` of the shed is ``n_before`` of the pickup.
* **Crossings** (one row per sign change of the depth ``r - R``): ion index,
  internal step, the linearly interpolated crossing time, and the direction
  (+1 outward, -1 inward). Outside is ``depth > 0`` (the §12 / G13 rule).

Step convention: internal step ``k`` (1-based, the driver's ``internal_id``)
advances ``t_{k-1} -> t_k``. Events are applied at its pre-step seam, at time
``t_{k-1}`` and positions ``x(t_{k-1})``; the stored column of step ``k`` holds
the post-event, post-BAOAB state at ``t_k``. A crossing in step ``k`` lies in
``(t_{k-1}, t_k]``.

Recording reads only flags and positions the step already computed: it draws
no random numbers, so a run with the log is byte-identical to one without.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np

#: File name inside a run directory.
STAGE1_EVENTS_FILENAME: str = "stage1_events.npz"
#: Own schema version (independent of the IonCheckpoint schema).
STAGE1_EVENTS_SCHEMA_VERSION: int = 1
#: ``event_kind`` codes.
KIND_SHED: int = 0
KIND_PICKUP: int = 1


@dataclass(frozen=True)
class Stage1EventLog:
    """Complete Stage-I event and surface-crossing record of one run.

    Event arrays share one length E, crossing arrays one length C. Times in ps,
    positions and depths in Angstrom (depth = ``r - R``, negative inside).

    Attributes
    ----------
    num_atoms : int
        2N, the number of ions (valid ``*_ion`` range ``[0, num_atoms)``).
    dt_ps : float
        Internal step ``dt_ion`` [ps].
    event_ion, event_step : np.ndarray of int64, (E,)
    event_time_ps : np.ndarray of float64, (E,)
        Pre-step time ``t_{k-1}`` at which the event is applied.
    event_kind : np.ndarray of int8, (E,)
        :data:`KIND_SHED` or :data:`KIND_PICKUP`.
    event_n_before, event_n_after : np.ndarray of int16, (E,)
        Occupancy before / after this event (differ by exactly 1).
    event_x_angstrom, event_y_angstrom, event_z_angstrom, event_depth_angstrom :
        np.ndarray of float64, (E,) -- the pre-step position.
    crossing_ion, crossing_step : np.ndarray of int64, (C,)
    crossing_time_ps : np.ndarray of float64, (C,)
        Interpolated time of ``depth = 0`` inside the step.
    crossing_direction : np.ndarray of int8, (C,)
        +1 outward (inside -> outside), -1 inward.
    """

    num_atoms: int
    dt_ps: float
    event_ion: np.ndarray
    event_step: np.ndarray
    event_time_ps: np.ndarray
    event_kind: np.ndarray
    event_n_before: np.ndarray
    event_n_after: np.ndarray
    event_x_angstrom: np.ndarray
    event_y_angstrom: np.ndarray
    event_z_angstrom: np.ndarray
    event_depth_angstrom: np.ndarray
    crossing_ion: np.ndarray
    crossing_step: np.ndarray
    crossing_time_ps: np.ndarray
    crossing_direction: np.ndarray

    def __post_init__(self) -> None:
        ev = [f.name for f in fields(self) if f.name.startswith("event_")]
        cr = [f.name for f in fields(self) if f.name.startswith("crossing_")]
        for group in (ev, cr):
            lengths = {np.shape(getattr(self, k)) for k in group}
            if len(lengths) != 1 or len(next(iter(lengths))) != 1:
                raise ValueError(f"Stage1EventLog: {group} must be equal-length 1-D "
                                 f"arrays, got shapes {sorted(lengths)}")
        for name in ("event_ion", "crossing_ion"):
            a = np.asarray(getattr(self, name))
            if a.size and (a.min() < 0 or a.max() >= self.num_atoms):
                raise ValueError(f"Stage1EventLog: {name} outside [0, {self.num_atoms})")
        step = np.asarray(self.event_n_after, int) - np.asarray(self.event_n_before, int)
        want = np.where(np.asarray(self.event_kind) == KIND_SHED, -1, 1)
        if not np.array_equal(step, want):
            raise ValueError("Stage1EventLog: n_after - n_before must be -1 for a shed "
                             "and +1 for a pickup")
        if not np.isin(self.crossing_direction, (-1, 1)).all():
            raise ValueError("Stage1EventLog: crossing_direction must be +1 or -1")


class Stage1EventRecorder:
    """Accumulates :class:`Stage1EventLog` rows step by step (driver-side).

    Parameters
    ----------
    num_atoms : int
        2N.
    dt_ps : float
        Internal step [ps].
    """

    def __init__(self, num_atoms: int, dt_ps: float) -> None:
        self._num_atoms = int(num_atoms)
        self._dt_ps = float(dt_ps)
        self._ev: list[tuple] = []
        self._cr: list[tuple] = []

    def record_step(self, step: int, time_ps: float, n_before: np.ndarray,
                    fired_shed: np.ndarray, fired_pickup: np.ndarray,
                    x: np.ndarray, y: np.ndarray, z: np.ndarray,
                    depth: np.ndarray) -> None:
        """Book this step's sheds then pickups.

        ``n_before`` is the pre-step occupancy (2N,); ``x, y, z, depth`` the
        pre-step positions the channels read; ``fired_*`` boolean (2N,).
        """
        n0 = np.asarray(n_before).astype(np.int64)
        shed = np.asarray(fired_shed, bool)
        pick = np.asarray(fired_pickup, bool)
        for kind, mask, nb in ((KIND_SHED, shed, n0),
                               (KIND_PICKUP, pick, n0 - shed)):
            idx = np.flatnonzero(mask)
            if idx.size == 0:
                continue
            delta = -1 if kind == KIND_SHED else 1
            self._ev.append((idx, np.full(idx.size, step), np.full(idx.size, time_ps),
                             np.full(idx.size, kind), nb[idx], nb[idx] + delta,
                             x[idx], y[idx], z[idx], depth[idx]))

    def record_crossings(self, step: int, time_ps: float,
                         depth_before: np.ndarray, depth_after: np.ndarray) -> None:
        """Book sign changes of the depth over step ``(time_ps, time_ps + dt]``."""
        out_b, out_a = depth_before > 0.0, depth_after > 0.0
        idx = np.flatnonzero(out_b != out_a)
        if idx.size == 0:
            return
        d0, d1 = depth_before[idx], depth_after[idx]
        frac = d0 / (d0 - d1)   # d0 and d1 straddle 0 and differ: no division by 0
        self._cr.append((idx, np.full(idx.size, step), time_ps + frac * self._dt_ps,
                         np.where(out_a[idx], 1, -1)))

    def finish(self) -> Stage1EventLog:
        """Concatenate the booked rows into a validated :class:`Stage1EventLog`."""
        ev_types = (np.int64, np.int64, np.float64, np.int8, np.int16, np.int16,
                    np.float64, np.float64, np.float64, np.float64)
        cr_types = (np.int64, np.int64, np.float64, np.int8)
        ev = [np.concatenate([r[i] for r in self._ev]).astype(t) if self._ev
              else np.zeros(0, t) for i, t in enumerate(ev_types)]
        cr = [np.concatenate([r[i] for r in self._cr]).astype(t) if self._cr
              else np.zeros(0, t) for i, t in enumerate(cr_types)]
        return Stage1EventLog(self._num_atoms, self._dt_ps, *ev, *cr)


def replay_n_shell(log: Stage1EventLog, n_initial: np.ndarray,
                   steps: np.ndarray) -> np.ndarray:
    """Occupancy after each internal step in ``steps``, replayed from the log.

    Parameters
    ----------
    log : Stage1EventLog
    n_initial : np.ndarray, shape (2N,)
        Occupancy at ``t_0`` (the checkpoint's column 0).
    steps : np.ndarray of int, shape (S,)
        Internal step indices (0 = the initial state).

    Returns
    -------
    np.ndarray of int64, shape (2N, S)
        ``n_initial`` plus the net events with ``event_step <= steps[j]``.
    """
    n0 = np.asarray(n_initial).astype(np.int64)
    if n0.shape != (log.num_atoms,):
        raise ValueError(f"n_initial shape {n0.shape} != ({log.num_atoms},)")
    steps = np.asarray(steps, np.int64)
    delta = np.where(log.event_kind == KIND_SHED, -1, 1)
    out = np.empty((log.num_atoms, steps.size), np.int64)
    for j, k in enumerate(steps):
        m = log.event_step <= k
        out[:, j] = n0 + np.bincount(log.event_ion[m], delta[m], log.num_atoms).astype(np.int64)
    return out


def save_stage1_events(log: Stage1EventLog, path: str | Path) -> Path:
    """Write ``log`` to ``.npz`` (compressed, pickle-free). Returns the path."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = {f.name: np.asarray(getattr(log, f.name)) for f in fields(log)}
    payload["stage1_events_schema_version"] = np.asarray(STAGE1_EVENTS_SCHEMA_VERSION)
    np.savez_compressed(p, **payload)
    return p


def load_stage1_events(path: str | Path) -> Stage1EventLog:
    """Load a :class:`Stage1EventLog` (fail-loud on schema / fields / shapes)."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"stage-1 event log not found: {p}")
    with np.load(p, allow_pickle=False) as npz:
        raw = {k: npz[k] for k in npz.files}
    version = raw.pop("stage1_events_schema_version", None)
    if version is None or int(version) != STAGE1_EVENTS_SCHEMA_VERSION:
        raise ValueError(f"{p}: stage1_events_schema_version={version}, this code "
                         f"expects {STAGE1_EVENTS_SCHEMA_VERSION}")
    expected = {f.name for f in fields(Stage1EventLog)}
    if set(raw) != expected:
        raise ValueError(f"{p}: fields differ from Stage1EventLog (missing "
                         f"{sorted(expected - set(raw))}, unknown {sorted(set(raw) - expected)})")
    raw["num_atoms"] = int(raw["num_atoms"])
    raw["dt_ps"] = float(raw["dt_ps"])
    return Stage1EventLog(**raw)
