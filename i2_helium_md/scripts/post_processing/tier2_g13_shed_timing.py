"""When and where the ion sheds He: shed timing and exit vs detection (G13).

Reads the thesis production members (``{law}_N1000_s1`` ... ``_s6``,
``self_bound_gate_release = decoupled``, THESIS_FINAL_GOALS G13; ion.npz
kept) of B2 and, as contrast, B3 and the re-landed h405. Runs no MD.
THESIS_DISCUSSION §16.

Every He shed of every scored ion over the whole flight is timed, from two
complete event records:

* MD window (0-500 ps): ``stage1_events.npz`` (``simulation/stage1_events``),
  one row per shed / pickup at the internal step dt_ion = 0.01 ps, with the
  pre-step time and position at which the channels fire. Anchored: it
  replays the stored ``n_shell`` at every column. The frame reconstruction
  (one event per unit change between stored columns, 0.16 ps; a pickup and
  a shed in one interval cancel) is kept only as the measured undercount
  (``frac_stage1_*_hidden_by_frames``);
* detection stage (500 ps - 8.53 us): one shed per recorded event.

Droplet exit follows the §12.2 definition: depth = |r| - R (per-atom sharp
radius ``droplet_radii_angstrom``, the half-density point of the 3.5 A erf
profile), exit = the **last** outward crossing depth <= 0 -> depth > 0 in
[0, t_h], from the log's crossing rows at the internal step, time linearly
interpolated inside the step. A Stage I shed is **before exit** when its
step is <= the exit step (it fires at the pre-step seam, before the motion
across the surface); n_exit = n right after the exit step, so n_exit - n_h
= sheds after exit - pickups after exit exactly. The geometric booking
(depth <= 0 where the shed is applied) is written as its own column.
Anchored: the log's exit falls in the stored interval of the stored-depth
exit column for every scored ion.

Classes by detected n: bare (n = 0), n = 1, n >= 2 (every row over the
scored ions, ``detected_mask``). Where the committed ``never`` members of
the same seeds are on disk (B2, h405), the ions those runs flagged
``suppressed`` are also read as their own row ("formerly suppressed").

Kinetic energy at four points per ion (all laws run ``co_moving`` sheds:
a shed changes m, not v -- required, checked from cfg):

* **exit**: the stored Stage I state at the first stored column after the
  exit (0-0.16 ps late; its n is ``n_exit_col``);
* **handover** (t_h, = just before the Coulomb closure, which acts before
  the Stage II event loop at handover masses and positions): the last
  stored column;
* **post-closure**: ½ m_h |v_det|² -- the closure changes v at fixed m and
  the Stage II sheds keep v, so v_det is the post-closure velocity;
* **detection**: ``E_kin_detected_eV``.

The exact kinematic split ln(KE_det / KE_exit) = ln(m_det / m_exit)
[mass loss] + 2 ln(v_det / v_exit) [speed change, all forces lumped:
Coulomb, well, drag tail, pickup recoil] is written per class. No
per-force split.

Anchors (fail loudly): per member, the event-log replay and exit-column
checks above; n(t_h) - #detection events == n_det for
every scored ion; every scored ion has an exit; stored KE == ½mv² at exit,
handover and detection; closure KE_pc - KE_h == -(E_pot_det - E_pot_h -
Σ fold) (the per-ion closure amount, ``detection_stage`` docstring); Stage
II KE_det - KE_pc == -Σ ``event_dE_mass_transfer_eV`` (all to
``KE_TOL_EV``, scored ions). Pooled: scored count, n̄ and KE₁ mean
reproduce the law's row of ``g13_decoupled_table.csv`` and the bare share
its row of ``thesis_observables_extra.csv`` (4 dp, THESIS_DISCUSSION §15).

Outputs (commit with ``git add -f``):

* ``data/runs/h2b_forward_model/g13_shed_timing.csv`` -- law x class:
  shed-time quantiles, late-shed shares (per shed and per ion), shares
  before exit / beyond the gate width / in the detection stage;
* ``g13_exit_vs_detection.csv`` -- law x class: n at birth, exit and
  detection, sheds and pickups after exit, exit-time quantiles;
* ``g13_exit_detection_hist.csv`` -- law x n: scored ions at exit / at
  detection;
* ``g13_shed_chord.csv`` -- law x emission-axis chord bin;
* ``g13_ke_readouts.csv`` -- law x class: mean n and KE at exit /
  handover / post-closure / detection, closure and Stage II KE changes,
  the log split into mass and speed terms, speed ratios;
* figures ``data/figures/thesis_tier2/g13_shed_timing.png`` and
  ``g13_exit_vs_detection.png`` (``FIG_LAW`` only; not committed; style
  of ``tier2_b2_thesis_figure``: text width, B2 red). No experimental
  curve: these figures compare readout times within the model; the
  experiment comparison of the detected sizes is §15 / B2-1. The spatial
  coordinate is written r - R (positive outside), not d, because the
  thesis uses d = R - r.

Reference: Albrechtsen et al., Nature 623, 319 (2023), Methods + Extended
Data Fig. 4 (classical MD of Na+He_N -- a simulation, drawn grey):
dissociation "almost exclusively during the first tens
of picoseconds"; "the vast majority of dissociation events are finished
after 400 ps".

Usage::

    python scripts/post_processing/tier2_g13_shed_timing.py
"""

from __future__ import annotations

import csv
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# =============================================================================
# USER SETTINGS
# =============================================================================
LAWS: tuple[str, ...] = ("B2", "B3", "h405")
FIG_LAW: str = "B2"                       # the law drawn in the figures
ANCHOR_CHECK: bool = True
LATE_PS: tuple[float, ...] = (50.0, 100.0, 400.0, 500.0)    # 400 = the paper's bound
LONG_PS: tuple[float, ...] = (1e3, 1e4, 1e5, 1e6)           # per-ion last-shed tail
GATE_WIDTH_A: float = 3.5                 # erf density-gate width (B2 cfg drag_gate_steepness)
DEPTH_MARKS_A: tuple[float, ...] = (GATE_WIDTH_A, 3 * GATE_WIDTH_A)
KE_TOL_EV: float = 1e-12                  # identity anchors are exact (measured ≤ 3e-16 eV)
CHORD_EDGES_A: tuple[float, ...] =(0.0, 30.0, 40.0, 50.0, 60.0, 70.0, np.inf)
LIT_TENS_OF_PS: tuple[float, float] = (10.0, 100.0)  # "first tens of ps" (grey span)
LIT_BOUND_PS: float = 400.0                           # "vast majority finished"
T_XLIM_PS: tuple[float, float] = (1.0, 1.5e7)         # first shed > 1 ps; detector 8.53e6
DEPTH_XMAX_A: float = 150.0               # panel (c) cut; the share beyond is annotated

OUT_DIR: Path = PROJECT_ROOT / "data" / "figures" / "thesis_tier2"
FORMATS: tuple[str, ...] = ("png",)
DPI: int = 200
SHOW_FIGURE: bool = False
from scripts.post_processing.thesis_style import (  # noqa: E402  (shared thesis style)
    FIG_RC, FONT, LIT_COLOR, SIM_COLOR, TEXTWIDTH_IN, THESIS_NAMES, VIRIDIS, panel_letter,
    style_axis,
)

STYLE: dict = {                          # data colours: thesis_style.VIRIDIS samples
    "sim": SIM_COLOR, "sim_exit": VIRIDIS[0.70], "lit": LIT_COLOR,
    "before_exit": VIRIDIS[0.10], "after_exit_stage1": VIRIDIS[0.50],
    "stage2": VIRIDIS[0.85], "cmap": "viridis",
    "class": {"all": (VIRIDIS[0.00], "-"), "bare": (VIRIDIS[0.33], "--"),
              "n1": (VIRIDIS[0.60], "-"), "n2plus": (VIRIDIS[0.85], "-.")},
    **FONT, "fig_height_in": 5.4,
}

# =============================================================================
# Imports (after the settings so SHOW_FIGURE picks the backend)
# =============================================================================
import matplotlib  # noqa: E402

if not SHOW_FIGURE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

from i2_helium_md.simulation.detection_stage import load_detection_result  # noqa: E402
from i2_helium_md.simulation.ion_propagation_step import _E_kin_eV  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from i2_helium_md.simulation.stage1_events import (  # noqa: E402
    KIND_SHED, STAGE1_EVENTS_FILENAME, load_stage1_events, replay_n_shell,
)
from scripts.gen_tier2_detfix_battery import MEMBER_SEEDS  # noqa: E402
from scripts.gen_tier2_detfix_postpickup import (  # noqa: E402
    B2_C, B2_E0_EV, B2_TAU_PS, lq_run_dir, thesis_run_dir,
)
from scripts.gen_tier2_h405_pp35 import SCAN_POOL_CELL, run_dir_for  # noqa: E402
from scripts.post_processing.tier2_confirmation_score import (  # noqa: E402
    format_table, write_rows_csv,
)
from scripts.post_processing.tier2_detfix_postpickup_table import emission_chords  # noqa: E402
from scripts.post_processing.tier2atlas_geometry_table import RUNS_ROOT  # noqa: E402

OUT_ROOT = RUNS_ROOT / "h2b_forward_model"
TIMING_CSV = OUT_ROOT / "g13_shed_timing.csv"
EXIT_CSV = OUT_ROOT / "g13_exit_vs_detection.csv"
HIST_CSV = OUT_ROOT / "g13_exit_detection_hist.csv"
CHORD_CSV = OUT_ROOT / "g13_shed_chord.csv"
KE_CSV = OUT_ROOT / "g13_ke_readouts.csv"
SHED_CONVENTION = "co_moving"            # the post-closure read needs v unchanged by sheds
G13_TABLE_CSV = OUT_ROOT / "g13_decoupled_table.csv"
EXTRA_CSV = OUT_ROOT / "thesis_observables_extra.csv"
SUPPRESSED = "suppressed"
CLASS_NAMES = {"all": "All scored ions", "bare": r"$n_\mathrm{det} = 0$",
               "n1": r"$n_\mathrm{det} = 1$", "n2plus": r"$n_\mathrm{det} \geq 2$",
               "supp_never": "formerly suppressed"}


# =============================================================================
# Pure helpers (unit-tested in tests/test_tier2_g13_shed_timing.py)
# =============================================================================
def exit_columns(depth_angstrom: np.ndarray) -> np.ndarray:
    """Stored column of the droplet exit state, per ion (§12.2).

    Parameters
    ----------
    depth_angstrom : (M, T) float
        Depth |r| - R [Å] at every stored column (negative inside).

    Returns
    -------
    (M,) int
        First column after the **last** outward crossing (d <= 0 -> d > 0);
        0 for an ion outside at column 0 that never re-enters; -1 for an
        ion inside at the last column (no exit inside the window).
    """
    d = np.asarray(depth_angstrom, dtype=float)
    if d.ndim != 2 or d.shape[1] < 2:
        raise ValueError(f"depth must be (M, T >= 2), got {d.shape}")
    out = d > 0
    cross = (~out[:, :-1]) & out[:, 1:]
    has = cross.any(axis=1)
    last = cross.shape[1] - 1 - np.argmax(cross[:, ::-1], axis=1)
    j = np.where(has, last + 1, np.where(out[:, 0], 0, -1))
    return np.where(out[:, -1], j, -1)


def window_events(n_shell: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Shed and pickup events from a stored shell-size history.

    One event per unit change between stored columns (a pickup and a shed
    inside the same interval cancel and are not seen).

    Parameters
    ----------
    n_shell : (M, T) float
        Shell size per ion and stored column (integer-valued).

    Returns
    -------
    shed_ion, shed_col, pickup_ion, pickup_col : (K,) int each
        Owning ion and the stored column **after** the event (the event lies
        in (t[col - 1], t[col]]).
    """
    n = np.asarray(n_shell, dtype=float)
    if n.ndim != 2:
        raise ValueError(f"n_shell must be (M, T), got {n.shape}")
    dn = np.rint(np.diff(n, axis=1)).astype(int)

    def _expand(k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        i, c = np.nonzero(k)
        reps = k[i, c]
        return np.repeat(i, reps), np.repeat(c + 1, reps)

    si, sc = _expand(np.clip(-dn, 0, None))
    pi, pc = _expand(np.clip(dn, 0, None))
    return si, sc, pi, pc


def exit_from_crossings(crossing_ion: np.ndarray, crossing_step: np.ndarray,
                        crossing_time_ps: np.ndarray, crossing_direction: np.ndarray,
                        outside_at_start: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Droplet exit per ion from the exact surface-crossing record (§12.2).

    Parameters
    ----------
    crossing_ion, crossing_step, crossing_time_ps, crossing_direction : (C,)
        The ``stage1_events.npz`` crossing rows (direction +1 outward).
    outside_at_start : (M,) bool
        depth > 0 at t = 0.

    Returns
    -------
    exit_step : (M,) int
        Internal step of the **last** outward crossing; 0 for an ion born
        outside that never re-enters; -1 for an ion inside at the end.
    exit_time_ps : (M,) float
        Its interpolated crossing time; 0.0 / NaN for those two cases.
    """
    start = np.asarray(outside_at_start, dtype=bool)
    ion = np.asarray(crossing_ion, dtype=int)
    direc = np.asarray(crossing_direction, dtype=int)
    m = start.size
    outside_end = start.astype(int) + np.bincount(ion, direc, minlength=m).astype(int)
    if not np.isin(outside_end, (0, 1)).all():
        raise ValueError("crossings do not alternate in direction per ion")
    step = np.where(start, 0, -1).astype(int)
    time = np.where(start, 0.0, np.nan)
    out = np.flatnonzero(direc == 1)
    order = out[np.argsort(np.asarray(crossing_step)[out], kind="stable")]
    step[ion[order]] = np.asarray(crossing_step)[order]        # later rows overwrite
    time[ion[order]] = np.asarray(crossing_time_ps)[order]
    inside_end = outside_end == 0
    return np.where(inside_end, -1, step), np.where(inside_end, np.nan, time)


def ke_log_split(m0_kg: np.ndarray, v0: np.ndarray, m1_kg: np.ndarray,
                 v1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact split of ln(KE1 / KE0) for KE = ½ m v², per ion.

    Parameters
    ----------
    m0_kg, m1_kg : (I,) float
        Mass [kg] at the two points (> 0).
    v0, v1 : (I,) float
        Speed at the two points (> 0; any common unit).

    Returns
    -------
    mass_term, speed_term : (I,) float
        ln(m1 / m0) and 2 ln(v1 / v0); their sum is ln(KE1 / KE0), with no
        dependence on the order of the mass and speed changes.
    """
    m0, m1, a, b = (np.asarray(x, dtype=float) for x in (m0_kg, m1_kg, v0, v1))
    if np.any(m0 <= 0) or np.any(m1 <= 0) or np.any(a <= 0) or np.any(b <= 0):
        raise ValueError("masses and speeds must be positive")
    return np.log(m1 / m0), 2.0 * np.log(b / a)


# =============================================================================
# Data layer
# =============================================================================
@dataclass
class LawRead:
    """Pooled per-ion and per-shed arrays of one law (all 2N ions per member)."""

    scored: np.ndarray          # (I,) bool
    n_birth: np.ndarray         # (I,) shell size at t = 0
    n_exit: np.ndarray          # (I,) right after the exit step (exact, from the log)
    n_exit_col: np.ndarray      # (I,) at the KE exit column (first stored column after exit)
    n_det: np.ndarray           # (I,)
    t_exit_ps: np.ndarray       # (I,) interpolated last outward crossing
    t_last_ps: np.ndarray       # (I,) last shed time (0 if none)
    pickups_after_exit: np.ndarray   # (I,)
    stage1_sheds: np.ndarray    # (I,) Stage I sheds, log
    stage1_pickups: np.ndarray  # (I,) Stage I pickups, log
    frame_sheds: np.ndarray     # (I,) Stage I sheds seen between stored columns
    frame_pickups: np.ndarray   # (I,) Stage I pickups seen between stored columns
    chord_A: np.ndarray         # (I,)
    supp_never: np.ndarray | None    # (I,) bool, None without never members
    shed_ion: np.ndarray        # (K,) pooled ion index
    shed_t_ps: np.ndarray       # (K,)
    shed_before_exit: np.ndarray     # (K,) bool: Stage I step <= exit step
    shed_depth_A: np.ndarray    # (K,) r - R where applied; +inf in detection stage
    shed_in_detection: np.ndarray    # (K,) bool
    n_h: np.ndarray             # (I,) shell size at handover (last stored column)
    ke_exit_eV: np.ndarray      # (I,) KE at the exit column
    ke_h_eV: np.ndarray         # (I,) KE at handover = before the closure
    ke_pc_eV: np.ndarray        # (I,) KE after the closure, before Stage II sheds
    ke_det_eV: np.ndarray       # (I,) KE at detection
    m_exit_kg: np.ndarray       # (I,)
    m_h_kg: np.ndarray          # (I,)
    m_det_kg: np.ndarray        # (I,)
    v_exit: np.ndarray          # (I,) speed [Å/ps]
    v_h: np.ndarray             # (I,)
    v_det: np.ndarray           # (I,)


def ke_readouts(rd: Path, ion, det, je: np.ndarray, scored: np.ndarray) -> dict[str, np.ndarray]:
    """KE, mass and speed at exit / handover / post-closure / detection.

    ``je`` (I,) is the exit column (handover column for ions without exit).
    The identity anchors of the module docstring are checked over the
    scored ions; any violation raises.
    """
    conv = RunDirectory(rd).load_cfg().evaporation_shed_convention
    if conv != SHED_CONVENTION:
        raise RuntimeError(f"{rd.name}: shed convention {conv!r}; post-closure read needs "
                           f"{SHED_CONVENTION!r}")
    ii = np.arange(je.size)
    m = np.asarray(ion.mass_history_kg, dtype=float)

    def vel(c):
        return (ion.velocities_x[ii, c], ion.velocities_y[ii, c], ion.velocities_z[ii, c])

    v_ex, v_h = vel(je), vel(-1)
    v_d = (np.asarray(det.vx_detected), np.asarray(det.vy_detected), np.asarray(det.vz_detected))
    m_ex, m_h, m_det = m[ii, je], m[:, -1], np.asarray(det.mass_detected_kg, dtype=float)
    ke_ex = np.asarray(ion.E_kin_eV, dtype=float)[ii, je]
    ke_h = np.asarray(ion.E_kin_eV, dtype=float)[:, -1]
    ke_det = np.asarray(det.E_kin_detected_eV, dtype=float)
    ke_pc = _E_kin_eV(m_h, *v_d)

    ev_ion = np.repeat(ii, np.diff(np.asarray(det.event_offsets)))
    fold = np.bincount(ev_ion, np.asarray(det.event_dE_bind_fold_eV), minlength=ii.size)
    mass_transfer = np.bincount(ev_ion, np.asarray(det.event_dE_mass_transfer_eV), minlength=ii.size)
    closure_pot = np.asarray(det.E_pot_detected_eV) - np.asarray(ion.E_pot_eV)[:, -1] - fold
    checks = {
        "KE_exit == ½mv²": ke_ex - _E_kin_eV(m_ex, *v_ex),
        "KE_h == ½mv²": ke_h - _E_kin_eV(m_h, *v_h),
        "KE_det == ½mv²": ke_det - _E_kin_eV(m_det, *v_d),
        "closure dKE == -dE_pot(closure)": (ke_pc - ke_h) + closure_pot,
        "Stage II dKE == -Σ event_dE_mass_transfer": (ke_det - ke_pc) + mass_transfer,
    }
    for name, diff in checks.items():
        worst = float(np.max(np.abs(diff[scored]))) if scored.any() else 0.0
        if worst > KE_TOL_EV:
            raise RuntimeError(f"{rd.name}: anchor {name} off by {worst:.3e} eV")
    speed = lambda v: np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)  # noqa: E731
    return {"n_h": np.asarray(ion.n_shell, dtype=float)[:, -1],
            "ke_exit_eV": ke_ex, "ke_h_eV": ke_h, "ke_pc_eV": ke_pc, "ke_det_eV": ke_det,
            "m_exit_kg": m_ex, "m_h_kg": m_h, "m_det_kg": m_det,
            "v_exit": speed(v_ex), "v_h": speed(v_h), "v_det": speed(v_d)}


def never_dirs(law: str) -> list[Path] | None:
    """The committed ``never`` members of the same seeds, where on disk."""
    if law == "B2":
        return [lq_run_dir(m, B2_E0_EV, B2_TAU_PS, B2_C) for m in MEMBER_SEEDS]
    if law == "h405":
        v_c, tau, e0 = SCAN_POOL_CELL
        return [run_dir_for(m, e0, v_c, tau) for m in MEMBER_SEEDS]
    return None


def read_member(rd: Path, never_rd: Path | None) -> dict[str, np.ndarray]:
    """Per-ion and per-shed arrays of one member (anchors checked here)."""
    ion = RunDirectory(rd).load_ion()
    t = np.asarray(ion.time_ps, dtype=float)
    n = np.asarray(ion.n_shell, dtype=float)
    radii = np.asarray(ion.droplet_radii_angstrom, dtype=float)[:, None]
    depth = np.sqrt(ion.positions_x ** 2 + ion.positions_y ** 2
                    + ion.positions_z ** 2) - radii
    det = load_detection_result(rd / "detection.npz")
    scored = np.asarray(det.detected_mask, dtype=bool)
    n_det = np.asarray(det.n_detected, dtype=float)
    n_ions = n.shape[0]

    log = load_stage1_events(rd / STAGE1_EVENTS_FILENAME)
    steps = np.rint(t / log.dt_ps).astype(int)          # internal step of each stored column
    if not np.array_equal(replay_n_shell(log, n[:, 0], steps), n.astype(int)):
        raise RuntimeError(f"{rd.name}: event log does not replay the stored n_shell")

    exit_step, t_exit = exit_from_crossings(log.crossing_ion, log.crossing_step,
                                            log.crossing_time_ps, log.crossing_direction,
                                            depth[:, 0] > 0)
    if np.any(scored & (exit_step < 0)):
        raise RuntimeError(f"{rd.name}: {np.sum(scored & (exit_step < 0))} scored ions without exit")
    # KE exit column: first stored column at or after the exit step; it must be
    # the frame-based §12.2 exit column (the crossing lies in its interval).
    je = np.where(exit_step >= 0, np.searchsorted(steps, exit_step), t.size - 1)
    j_frames = exit_columns(depth)
    if np.any(scored & (je != j_frames)):
        raise RuntimeError(f"{rd.name}: log exit column != stored-depth exit column for "
                           f"{np.sum(scored & (je != j_frames))} scored ions")

    off = np.asarray(det.event_offsets)
    n_ev = np.diff(off)
    bad = scored & (np.abs(n[:, -1] - n_ev - n_det) > 1e-9)
    if bad.any():
        raise RuntimeError(f"{rd.name}: n(t_h) - events != n_det for {bad.sum()} scored ions")
    ev_ion = np.repeat(np.arange(n_ions), n_ev)
    ev_t = np.asarray(det.event_time_ps, dtype=float)

    is_shed = log.event_kind == KIND_SHED
    si, s_step = log.event_ion[is_shed], log.event_step[is_shed]
    pi, p_step = log.event_ion[~is_shed], log.event_step[~is_shed]
    shed_ion = np.concatenate([si, ev_ion])
    shed_t = np.concatenate([log.event_time_ps[is_shed], ev_t])
    t_last = np.zeros(n_ions)
    np.maximum.at(t_last, shed_ion, shed_t)
    up_to_exit = log.event_step <= exit_step[log.event_ion]
    net_to_exit = np.bincount(log.event_ion[up_to_exit],
                              np.where(is_shed, -1, 1)[up_to_exit], minlength=n_ions)
    fsi, _, fpi, _ = window_events(n)

    supp = None
    if never_rd is not None:
        nv = load_detection_result(never_rd / "detection.npz")
        supp = np.asarray(nv.state_reason).astype(str) == SUPPRESSED
    return {
        "scored": scored, "n_birth": n[:, 0], "n_exit": n[:, 0] + net_to_exit,
        "n_exit_col": n[np.arange(n_ions), je],
        "n_det": n_det, "t_exit_ps": t_exit, "t_last_ps": t_last,
        "pickups_after_exit": np.bincount(pi[p_step > exit_step[pi]], minlength=n_ions),
        "stage1_sheds": np.bincount(si, minlength=n_ions),
        "stage1_pickups": np.bincount(pi, minlength=n_ions),
        "frame_sheds": np.bincount(fsi, minlength=n_ions),
        "frame_pickups": np.bincount(fpi, minlength=n_ions),
        "chord_A": emission_chords(rd), "supp_never": supp,
        "shed_ion": shed_ion, "shed_t_ps": shed_t,
        "shed_before_exit": np.concatenate([s_step <= exit_step[si], np.zeros(ev_t.size, bool)]),
        "shed_depth_A": np.concatenate([log.event_depth_angstrom[is_shed],
                                        np.full(ev_t.size, np.inf)]),
        "shed_in_detection": np.concatenate([np.zeros(si.size, bool), np.ones(ev_t.size, bool)]),
        **ke_readouts(rd, ion, det, je, scored),
    }


def read_law(law: str) -> LawRead:
    """Pool the six members of ``law``; sheds of unscored ions dropped."""
    nv = never_dirs(law)
    parts, offset = [], 0
    for k, m in enumerate(MEMBER_SEEDS):
        p = read_member(thesis_run_dir(law, m), nv[k] if nv else None)
        p["shed_ion"] = p["shed_ion"] + offset
        offset += p["scored"].size
        parts.append(p)
    cat = {key: np.concatenate([p[key] for p in parts])
           for key in parts[0] if key != "supp_never"}
    cat["supp_never"] = (np.concatenate([p["supp_never"] for p in parts]) if nv else None)
    keep = cat["scored"][cat["shed_ion"]]
    for key in ("shed_ion", "shed_t_ps", "shed_before_exit", "shed_depth_A", "shed_in_detection"):
        cat[key] = cat[key][keep]
    return LawRead(**cat)


def check_anchor(law: str, r: LawRead) -> None:
    """Pooled scored count, n̄, KE₁ mean and bare share vs the committed §15
    rows (4 dp)."""
    with open(G13_TABLE_CSV, encoding="utf-8") as fh:
        row = next(x for x in csv.DictReader(fh) if x["label"] == law)
    with open(EXTRA_CSV, encoding="utf-8") as fh:
        extra = next(x for x in csv.DictReader(fh) if x["label"] == law)
    s = r.scored
    fresh = {"num_scored": float(s.sum()), "nbar_det": r.n_det[s].mean(),
             "bare_share": np.mean(r.n_det[s] == 0),
             "KE1_mean": r.ke_det_eV[s & (r.n_det == 1)].mean()}
    ref = {"num_scored": float(row["num_scored"]), "nbar_det": float(row["nbar_det"]),
           "bare_share": float(extra["bare_share"]), "KE1_mean": float(row["KE1_mean"])}
    for k in fresh:
        if round(fresh[k], 4) != round(ref[k], 4):
            raise RuntimeError(f"anchor {law} {k}: {fresh[k]} vs committed {ref[k]}")


def class_masks(r: LawRead) -> dict[str, np.ndarray]:
    """Per-ion masks (over the scored ions) of the final-state classes."""
    s = r.scored
    out = {"all": s.copy(), "bare": s & (r.n_det == 0), "n1": s & (r.n_det == 1),
           "n2plus": s & (r.n_det >= 2)}
    if r.supp_never is not None:
        out["supp_never"] = s & r.supp_never
    return out


def _q(x: np.ndarray, p: float) -> float:
    return float(np.percentile(x, p)) if x.size else float("nan")


def timing_row(law: str, cls: str, r: LawRead, m: np.ndarray) -> dict[str, Any]:
    """One ``g13_shed_timing.csv`` row: sheds of the ions in mask ``m``."""
    sel = m[r.shed_ion]
    t, ion = r.shed_t_ps[sel], r.shed_ion[sel]
    depth, n_ions = r.shed_depth_A[sel], int(m.sum())
    row: dict[str, Any] = {
        "law": law, "class": cls, "ions": n_ions,
        "share_of_scored": n_ions / r.scored.sum(),
        "sheds_total": int(t.size), "sheds_per_ion": t.size / max(n_ions, 1),
        "t_shed_q10_ps": _q(t, 10), "t_shed_median_ps": _q(t, 50),
        "t_shed_q90_ps": _q(t, 90), "t_shed_q99_ps": _q(t, 99),
    }
    for thr in LATE_PS:
        row[f"frac_sheds_after_{thr:g}ps"] = float(np.mean(t > thr)) if t.size else np.nan
    for thr in LATE_PS:
        row[f"frac_ions_shed_after_{thr:g}ps"] = np.unique(ion[t > thr]).size / max(n_ions, 1)
    for thr in LONG_PS:
        row[f"frac_ions_shed_after_{thr:g}ps"] = float(np.mean(r.t_last_ps[m] > thr))
    row["t_last_shed_median_ps"] = _q(r.t_last_ps[m], 50)
    row["frac_sheds_before_exit"] = float(np.mean(r.shed_before_exit[sel])) if t.size else np.nan
    row["frac_sheds_depth_le_0"] = float(np.mean(depth <= 0)) if t.size else np.nan
    for g in DEPTH_MARKS_A:
        row[f"frac_sheds_depth_gt_{g:g}A"] = float(np.mean(depth > g)) if t.size else np.nan
    row["frac_sheds_detection_stage"] = float(np.mean(r.shed_in_detection[sel])) if t.size else np.nan
    n1 = r.stage1_sheds[m].sum()
    row["frac_stage1_sheds_hidden_by_frames"] = (1.0 - r.frame_sheds[m].sum() / n1) if n1 else np.nan
    after = ~r.shed_before_exit[sel]
    row["t_after_exit_median_ps"] = _q(t[after] - r.t_exit_ps[ion[after]], 50)
    return row


def exit_row(law: str, cls: str, r: LawRead, m: np.ndarray) -> dict[str, Any]:
    """One ``g13_exit_vs_detection.csv`` row (ions in mask ``m``)."""
    sel = m[r.shed_ion]
    n_ions = max(int(m.sum()), 1)
    return {
        "law": law, "class": cls, "ions": int(m.sum()),
        "n_birth_mean": r.n_birth[m].mean(), "n_exit_mean": r.n_exit[m].mean(),
        "n_det_mean": r.n_det[m].mean(),
        "n_exit_minus_det_mean": (r.n_exit - r.n_det)[m].mean(),
        "sheds_before_exit_per_ion": np.sum(r.shed_before_exit[sel]) / n_ions,
        "sheds_after_exit_per_ion": np.sum(~r.shed_before_exit[sel]) / n_ions,
        "pickups_after_exit_per_ion": r.pickups_after_exit[m].sum() / n_ions,
        "stage1_pickups_per_ion": r.stage1_pickups[m].sum() / n_ions,
        "frac_stage1_pickups_hidden_by_frames": (
            1.0 - r.frame_pickups[m].sum() / r.stage1_pickups[m].sum()
            if r.stage1_pickups[m].sum() else np.nan),
        "frac_ions_n_changed_after_exit": float(np.mean(r.n_exit[m] != r.n_det[m])),
        "t_exit_q10_ps": _q(r.t_exit_ps[m], 10), "t_exit_median_ps": _q(r.t_exit_ps[m], 50),
        "t_exit_q90_ps": _q(r.t_exit_ps[m], 90),
        "frac_formerly_suppressed": (float(np.mean(r.supp_never[m]))
                                     if r.supp_never is not None else np.nan),
    }


def hist_rows(law: str, r: LawRead) -> list[dict[str, Any]]:
    s = r.scored
    n_max = int(max(r.n_exit[s].max(), r.n_det[s].max(), r.n_birth[s].max()))
    ce = np.bincount(r.n_exit[s].astype(int), minlength=n_max + 1)
    cd = np.bincount(r.n_det[s].astype(int), minlength=n_max + 1)
    return [{"law": law, "n": k, "ions_at_exit": int(ce[k]), "ions_at_detection": int(cd[k]),
             "frac_at_exit": ce[k] / s.sum(), "frac_at_detection": cd[k] / s.sum()}
            for k in range(n_max + 1)]


def chord_rows(law: str, r: LawRead) -> list[dict[str, Any]]:
    rows = []
    for lo, hi in zip(CHORD_EDGES_A[:-1], CHORD_EDGES_A[1:]):
        m = r.scored & (r.chord_A >= lo) & (r.chord_A < hi)
        if not m.any():
            continue
        sel = m[r.shed_ion]
        t = r.shed_t_ps[sel]
        rows.append({
            "law": law, "chord_lo_A": lo, "chord_hi_A": hi, "ions": int(m.sum()),
            "mean_n_det": r.n_det[m].mean(), "bare_share": float(np.mean(r.n_det[m] == 0)),
            "n1_share": float(np.mean(r.n_det[m] == 1)), "t_shed_median_ps": _q(t, 50),
            "frac_sheds_after_50ps": float(np.mean(t > 50)),
            "frac_sheds_after_400ps": float(np.mean(t > 400)),
            "frac_sheds_before_exit": float(np.mean(r.shed_before_exit[sel])),
            "t_exit_median_ps": _q(r.t_exit_ps[m], 50),
        })
    return rows


def ke_row(law: str, cls: str, r: LawRead, m: np.ndarray) -> dict[str, Any]:
    """One ``g13_ke_readouts.csv`` row (ions in mask ``m``); KE in eV."""
    mass_term, speed_term = ke_log_split(r.m_exit_kg[m], r.v_exit[m], r.m_det_kg[m], r.v_det[m])
    return {
        "law": law, "class": cls, "ions": int(m.sum()),
        "n_exit_col_mean": r.n_exit_col[m].mean(), "n_h_mean": r.n_h[m].mean(),
        "n_det_mean": r.n_det[m].mean(),
        "KE_exit_mean_eV": r.ke_exit_eV[m].mean(), "KE_h_mean_eV": r.ke_h_eV[m].mean(),
        "KE_postclosure_mean_eV": r.ke_pc_eV[m].mean(), "KE_det_mean_eV": r.ke_det_eV[m].mean(),
        "dKE_exit_to_h_mean_eV": (r.ke_h_eV - r.ke_exit_eV)[m].mean(),
        "dKE_closure_mean_eV": (r.ke_pc_eV - r.ke_h_eV)[m].mean(),
        "dKE_stage2_sheds_mean_eV": (r.ke_det_eV - r.ke_pc_eV)[m].mean(),
        "lnKE_ratio_det_exit_mean": (mass_term + speed_term).mean(),
        "lnKE_mass_term_mean": mass_term.mean(), "lnKE_speed_term_mean": speed_term.mean(),
        "v_det_over_v_exit_mean": (r.v_det / r.v_exit)[m].mean(),
        "v_det_over_v_h_mean": (r.v_det / r.v_h)[m].mean(),
    }


# =============================================================================
# Figures
# =============================================================================
TIME_MARK_LABELS: tuple[str, ...] = ("Handover (500 ps)", "Detector (8.53 µs)",
                                     "Literature bound (400 ps)")


def _time_marks(ax, t_h: float, t_det: float, *, label: bool, band: bool,
                ymax: float = 1.0) -> None:
    """Handover, detector and Albrechtsen 2023 400 ps lines (legend labels in
    ``TIME_MARK_LABELS`` order when ``label``); with ``band`` also the grey
    tens-of-ps band. ``ymax`` (axes fraction) stops the lines below a legend."""
    if band:
        lo, hi = LIT_TENS_OF_PS
        ax.axvspan(lo, hi, color=STYLE["lit"], alpha=0.15, lw=0)
    for x, color, lw, ls, lab in ((t_h, "k", 0.8, "--", TIME_MARK_LABELS[0]),
                                  (t_det, "k", 0.8, "-.", TIME_MARK_LABELS[1]),
                                  (LIT_BOUND_PS, STYLE["lit"], 1.0, ":", TIME_MARK_LABELS[2])):
        ax.axvline(x, ymax=ymax, color=color, lw=lw, ls=ls, label=lab if label else None)


def _ecdf(ax, x: np.ndarray, color: str, ls: str, label: str) -> None:
    xs = np.sort(x)
    ax.step(xs, np.arange(1, xs.size + 1) / xs.size, where="post", color=color, ls=ls,
            lw=1.3, label=label)


def plot_timing(r: LawRead, t_h: float, t_det: float) -> plt.Figure:
    """(a) shed-time density per decade, stacked by readout (Stage I before
    exit / Stage I after exit / Stage II); (b) shed-time ECDF per detected
    class; (c) per-ion ECDF of the last recorded evaporation event."""
    fig, ax = plt.subplot_mosaic([["a", "a"], ["b", "c"]],
                                 figsize=(TEXTWIDTH_IN, STYLE["fig_height_in"]),
                                 constrained_layout=True)
    edges =np.logspace(-1, np.log10(t_det * 1.05), 60)
    width = np.diff(np.log10(edges))
    n_sh = r.shed_t_ps.size
    groups = [
        (r.shed_before_exit, STYLE["before_exit"], "Stage I, before droplet exit"),
        (~r.shed_before_exit & ~r.shed_in_detection, STYLE["after_exit_stage1"],
         "Stage I, after droplet exit"),
        (r.shed_in_detection, STYLE["stage2"], "Stage II"),
    ]
    bottom = np.zeros(edges.size - 1)
    for g, color, lab in groups:
        h, _ = np.histogram(r.shed_t_ps[g], bins=edges)
        h = h / n_sh / width                    # dF/dlog10(t): fraction of all sheds / bin width in decades
        ax["a"].bar(edges[:-1], h, width=np.diff(edges), bottom=bottom, align="edge",
                    color=color, lw=0, label=f"{lab} ({np.mean(g) * 100:.1f} %)")
        bottom += h
    _time_marks(ax["a"], t_h, t_det, label=True, band=False, ymax=0.68)  # lines end below legend
    ax["a"].set_xscale("log")
    ax["a"].set_yscale("log")
    ax["a"].set_ylim(1e-5, 1e3)        # headroom for the legend
    ax["a"].set_xlim(*T_XLIM_PS)
    ax["a"].set_xlabel(r"Shed time $t$ / ps")
    ax["a"].set_ylabel("Fraction of sheds per decade")
    handles, labels = ax["a"].get_legend_handles_labels()
    order = [labels.index(x) for x in labels if x not in TIME_MARK_LABELS]
    order += [labels.index(x) for x in TIME_MARK_LABELS]     # second legend column
    ax["a"].legend([handles[i] for i in order], [labels[i] for i in order],
                   fontsize=STYLE["legend_fontsize"], ncol=2, loc="upper right",
                   frameon=True, facecolor="white", edgecolor="none", framealpha=1.0)
    ax["a"].set_yticks([1e-4, 1e-2, 1])
    masks = class_masks(r)
    for key in ("all", "bare", "n1", "n2plus"):
        color, ls = STYLE["class"][key]
        m = masks[key]
        _ecdf(ax["b"], r.shed_t_ps[m[r.shed_ion]], color, ls,
              f"{CLASS_NAMES[key]} ({m.sum()})")
        _ecdf(ax["c"], r.t_last_ps[m], color, ls, CLASS_NAMES[key])
    for k in ("b", "c"):
        _time_marks(ax[k], t_h, t_det, label=False, band=True)
        ax[k].set_xscale("log")
        ax[k].set_xlim(*T_XLIM_PS)
        ax[k].set_ylim(0, 1.02)
    ax["b"].text(0.93, 0.50, f"Median: {np.median(r.shed_t_ps):.1f} ps\n$t > 400$ ps: "
                 f"{np.mean(r.shed_t_ps > 400) * 100:.1f} % of sheds",
                 transform=ax["b"].transAxes, ha="right", va="center",
                 fontsize=STYLE["annot_fontsize"])
    ax["b"].set_xlabel(r"Shed time $t$ / ps")
    ax["b"].set_ylabel("Cumulative fraction of sheds")
    s = r.scored
    ax["c"].text(0.93, 0.30,
                 rf"$t_\mathrm{{last}} > 400$ ps: {np.mean(r.t_last_ps[s] > 400) * 100:.1f} % of ions"
                 "\n"
                 rf"$t_\mathrm{{last}} > 1$ µs: {np.mean(r.t_last_ps[s] > 1e6) * 100:.1f} % of ions",
                 transform=ax["c"].transAxes, ha="right", va="center",
                 fontsize=STYLE["annot_fontsize"])
    ax["c"].set_xlabel(r"Last recorded evaporation event $t_\mathrm{last}$ / ps")
    ax["c"].set_ylabel("Cumulative fraction of ions")
    ax["b"].legend(fontsize=STYLE["legend_fontsize"], frameon=False, loc="lower right")
    for k, a in ax.items():
        style_axis(a)
        panel_letter(a, k)
    return fig


def plot_exit(r: LawRead, law: str) -> plt.Figure:
    """(a) n at droplet exit vs n at detection (scored ions, both bars on
    the same normalisation) with the mean-n readout at exit / handover /
    detection; (b) per-ion n_exit vs n_det; (c) r - R at shed."""
    fig, ax = plt.subplot_mosaic([["a", "a"], ["b", "c"]],
                                 figsize=(TEXTWIDTH_IN, STYLE["fig_height_in"]),
                                 constrained_layout=True)
    s = r.scored
    rows = hist_rows(law, r)
    n = np.array([x["n"] for x in rows])
    fe = np.array([x["frac_at_exit"] for x in rows]) * 100      # percent of scored ions
    fd = np.array([x["frac_at_detection"] for x in rows]) * 100
    name = THESIS_NAMES[law]                                     # printed name (B2 -> LQ)
    ax["a"].bar(n - 0.2, fd, width=0.4, color=STYLE["sim"], label=f"{name}, at detection")
    ax["a"].bar(n + 0.2, fe, width=0.4, color=STYLE["sim_exit"],
                label=f"{name}, at droplet exit")
    ax["a"].set_xlabel(r"Shell size $n$")
    ax["a"].set_ylabel("Fraction of scored ions / %")
    ax["a"].set_xlim(-0.7, n.max() + 0.7)
    ax["a"].set_ylim(0, 1.3 * max(fe.max(), fd.max()))      # headroom for legend + readout
    ax["a"].xaxis.set_major_locator(MaxNLocator(integer=True))
    readout = "\n".join((
        rf"$\bar{{n}}_\mathrm{{exit}} = {r.n_exit[s].mean():.2f}$"
        f"  (median exit {np.median(r.t_exit_ps[s]):.1f} ps)",
        rf"$\bar{{n}}_\mathrm{{handover}} = {r.n_h[s].mean():.2f}$  (500 ps)",
        rf"$\bar{{n}}_\mathrm{{det}} = {r.n_det[s].mean():.2f}$  (8.53 µs)",
    ))
    ax["a"].text(0.36, 0.95, readout, transform=ax["a"].transAxes, va="top",
                 linespacing=1.5, fontsize=STYLE["annot_fontsize"])
    ax["a"].legend(fontsize=STYLE["legend_fontsize"], frameon=False, loc="upper right")

    hi = int(max(r.n_exit[s].max(), r.n_det[s].max())) + 1
    h, xe, ye = np.histogram2d(r.n_exit[s], r.n_det[s], bins=[np.arange(-0.5, hi + 0.5)] * 2)
    h = np.where(h > 0, h, np.nan)
    im = ax["b"].pcolormesh(xe, ye, h.T, cmap=STYLE["cmap"],norm=matplotlib.colors.LogNorm())
    ax["b"].plot([0, hi], [0, hi], color="0.5", lw=0.8, ls=":")
    cb = fig.colorbar(im, cax=ax["b"].inset_axes([1.03, 0.0, 0.05, 1.0]))
    cb.set_label("Ions", fontsize=STYLE["label_fontsize"])
    cb.ax.tick_params(labelsize=STYLE["tick_fontsize"])
    ax["b"].set_xlabel(r"Shell size at droplet exit $n_\mathrm{exit}$")
    ax["b"].set_ylabel(r"Shell size at detection $n_\mathrm{det}$")

    d = r.shed_depth_A[np.isfinite(r.shed_depth_A)]
    bins = np.arange(np.floor(d.min() / 2) * 2, min(d.max(), DEPTH_XMAX_A) + 2, 2.0)
    ax["c"].hist(d, bins=bins, color=STYLE["sim"], lw=0,
                 weights=np.full(d.size, 1.0 / r.shed_depth_A.size))
    ax["c"].axvline(0, color="k", lw=0.8)
    for g in DEPTH_MARKS_A:
        ax["c"].axvline(g, color="k", lw=0.8, ls=":")
    ax["c"].text(0.97, 0.95,
                 f"Before exit: {np.mean(r.shed_before_exit) * 100:.1f} %\n"
                 + "\n".join(f"$r-R > {g:g}$ Å: {np.mean(r.shed_depth_A > g) * 100:.1f} %"
                             for g in DEPTH_MARKS_A)
                 + f"\n$r-R > {DEPTH_XMAX_A:g}$ Å: "
                 f"{np.mean(r.shed_depth_A > DEPTH_XMAX_A) * 100:.1f} %\n(not shown)",
                 transform=ax["c"].transAxes, ha="right", va="top",
                 fontsize=STYLE["annot_fontsize"])
    ax["c"].set_xlabel(r"$r-R$ at shed / Å")
    ax["c"].set_ylabel("Fraction of all sheds per bin")
    for k, a in ax.items():
        style_axis(a)
        panel_letter(a, k)
    return fig


def main() -> int:
    warnings.simplefilter("ignore", RuntimeWarning)   # documented pairing warning
    timing, exits, hist, chords, kes = [], [], [], [], []
    reads: dict[str, LawRead] = {}
    for law in LAWS:
        r = read_law(law)
        if ANCHOR_CHECK:
            check_anchor(law, r)
        reads[law] = r
        for cls, m in class_masks(r).items():
            timing.append(timing_row(law, cls, r, m))
            exits.append(exit_row(law, cls, r, m))
            kes.append(ke_row(law, cls, r, m))
        hist += hist_rows(law, r)
        chords += chord_rows(law, r)
    show = ("law", "class", "ions", "sheds_per_ion", "t_shed_median_ps",
            "frac_sheds_after_50ps", "frac_sheds_after_100ps", "frac_sheds_after_400ps",
            "frac_ions_shed_after_400ps", "frac_sheds_before_exit")
    print("\n=== shed timing (pooled s1-s6, decoupled; scored ions) ===")
    print(format_table([{c: x[c] for c in show} for x in timing]))
    show = ("law", "class", "ions", "n_exit_mean", "n_det_mean", "t_exit_median_ps",
            "sheds_after_exit_per_ion", "pickups_after_exit_per_ion")
    print("\n=== exit vs detection ===")
    print(format_table([{c: x[c] for c in show} for x in exits]))
    show = ("law", "class", "n_h_mean", "KE_exit_mean_eV", "KE_h_mean_eV",
            "KE_postclosure_mean_eV", "KE_det_mean_eV", "lnKE_mass_term_mean",
            "lnKE_speed_term_mean", "v_det_over_v_exit_mean")
    print("\n=== KE readouts (exit / handover / post-closure / detection) ===")
    print(format_table([{c: x[c] for c in show} for x in kes]))
    for path, rows in ((TIMING_CSV, timing), (EXIT_CSV, exits), (HIST_CSV, hist),
                       (CHORD_CSV, chords), (KE_CSV, kes)):
        print(f"CSV -> {write_rows_csv(path, rows)}")

    cfg = RunDirectory(thesis_run_dir(FIG_LAW, next(iter(MEMBER_SEEDS)))).load_cfg()
    t_h, t_det = float(cfg.ion_simulation_time), float(cfg.detection_time_ps)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(FIG_RC):
        for stem, fig in (("g13_shed_timing", plot_timing(reads[FIG_LAW], t_h, t_det)),
                          ("g13_exit_vs_detection", plot_exit(reads[FIG_LAW], FIG_LAW))):
            for ext in FORMATS:
                path = OUT_DIR / f"{stem}.{ext}"
                fig.savefig(path, dpi=DPI)
                print(f"figure -> {path}")
    if SHOW_FIGURE:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
