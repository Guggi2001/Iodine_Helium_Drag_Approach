"""G8 bare-start channel screen — fragments born with no He shell (P1, Case 1).

THESIS_FINAL_GOALS.md G8 (frozen 2026-10-06). Two MD runs, production
pipeline unchanged (Stage I to 500 ps, E2 off, partner-aware closure,
N = 1000, seed s1). Each run is its parent's cfg with exactly one key moved,
``initial_shell_model: "density_tied" -> "bare"`` (n_0 = 0 for every
fragment, both fragments of every pair: P1). The parents already run the T6
``sigma_proportional`` onset law (inert for them: n_0 = 21 at production
depths), which gives E_int(0) = onset * Sigma(0)/Sigma(n*) = 0 for a bare
ion (no shell, no reservoir); the oracle checks the parent's law.

Parents:

* ``lq``   — the G4.14 A5 cell (c 2.0, tau 6.4, E0 0.3155, 3.5 A gate), s1;
* ``h405`` — the production detfix s1 member (14.2 A gate).

Oracles before any MD (``--dry-run`` runs them only):

1. the rebuilt parent cfg equals the parent run's stored ``cfg.json``;
2. the parent runs ``sigma_proportional``, and the bare cfg diffs from it in
   exactly ``initial_shell_model``;
3. handover alignment (the production generator's check, reused).

After MD: the neutral stage reproduces the parent's ``neutral.npz``
bit-exactly (same births), and ``bare_summary.npz`` (per ion: n at handover,
max n, pickup exposure int rho_hat dt, exit time, max E_int, speed at
handover) is written before ``ion.npz`` is dropped (disk; the scorer reads
cfg.json + detection.npz + the summary).

Usage::

    python scripts/gen_tier2_bare_screen.py --dry-run
    python scripts/gen_tier2_bare_screen.py            # both laws, 2 at a time
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
from pathlib import Path
import sys
from typing import Optional

for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.gen_tier2_detfix_battery import (  # noqa: E402
    MAX_BYTES_ION, build_member, member_run_dir_name, verify_handover_alignment,
)
from scripts.gen_tier2_g413_ridge_battery import LQ_BUNDLE  # noqa: E402
from scripts.gen_tier2_g414_tauline import (  # noqa: E402
    ION_DROPPED_MARKER, build as g414_build, confirm_cells,
    run_dir_for as g414_run_dir_for,
)
from scripts.gen_tier2atlas_g4finals import CASE, N, VARIANT  # noqa: E402
from scripts.gen_tier2atlas_linclone import cfg_field_diff  # noqa: E402
from scripts.tier0_common import run_dir_name  # noqa: E402
from scripts.tier2_common import cfg_diff_vs_reference  # noqa: E402
from i2_helium_md.config import SimConfig, drag_gate_steepness  # noqa: E402
from i2_helium_md.physics.helium_density import rho_he_ratio  # noqa: E402

# =============================================================================
# USER SETTINGS (pre-registered G8; FIXED)
# =============================================================================
LAWS = ("lq", "h405")
PARENT_MEMBER = "s1"
BARE_KEYS = {"initial_shell_model": "bare",
             "internal_energy_partition_law": "sigma_proportional"}
EXIT_RHO_HAT = 0.5          # exit time = last stored time with rho_hat >= 0.5
MAX_CONCURRENCY = 2
SUMMARY_NAME = "bare_summary.npz"
_REQUIRED_ARTIFACTS = ("cfg.json", "neutral.npz", "detection.npz", SUMMARY_NAME)


def _lq_parent_cell():
    return next(c for c in confirm_cells() if c.label.endswith(PARENT_MEMBER))


def parent_cfg(law: str) -> SimConfig:
    """The parent's cfg, rebuilt by its own generator."""
    if law == "lq":
        return g414_build(_lq_parent_cell())
    if law == "h405":
        return build_member(PARENT_MEMBER)
    raise KeyError(f"unknown law {law!r}; expected one of {LAWS}")


def parent_run_dir(law: str) -> Path:
    if law == "lq":
        return g414_run_dir_for(_lq_parent_cell())
    return PROJECT_ROOT / "data" / "runs" / member_run_dir_name(PARENT_MEMBER)


def build(law: str) -> SimConfig:
    cfg = dataclasses.replace(parent_cfg(law), **BARE_KEYS)
    cfg.validate()
    return cfg


def run_dir_for(law: str) -> Path:
    bundle = LQ_BUNDLE if law == "lq" else VARIANT
    name = run_dir_name(CASE, bundle, N, run_tag=f"detfix_conf270_g8bare{law}{PARENT_MEMBER}_th500")
    if "_tier2_" in name or "tier2probe" in name:
        raise AssertionError(f"namespace violated by {name!r}")
    return PROJECT_ROOT / "data" / "runs" / name


def verify(law: str) -> list[str]:
    """Oracles 1 and 2. Returns the bare-vs-parent diff keys."""
    parent = parent_cfg(law)
    stored = cfg_diff_vs_reference(parent, parent_run_dir(law) / "cfg.json",
                                   context=f"G8 parent {law}")
    if stored:
        raise AssertionError(f"[{law}] rebuilt parent differs from its cfg.json in {stored}")
    if parent.internal_energy_partition_law != BARE_KEYS["internal_energy_partition_law"]:
        raise AssertionError(f"[{law}] parent onset law is "
                             f"{parent.internal_energy_partition_law!r}, not sigma_proportional")
    diff = cfg_field_diff(parent, build(law))
    if diff != ["initial_shell_model"]:
        raise AssertionError(f"[{law}] bare cfg diff {diff}; pre-registered ['initial_shell_model']")
    return diff


def ion_summary(ion, cfg: SimConfig) -> dict[str, np.ndarray]:
    """Per-ion bare-class mechanics from the in-memory ion checkpoint.

    Returns arrays of shape ``(2N,)``: ``n_handover`` / ``n_max`` (shell
    count at the last stored column / its maximum), ``exposure_ps``
    (trapezoid int rho_hat dt over the stored columns, rho_hat at the cfg's
    gate steepness), ``t_exit_ps`` (last stored time with rho_hat >=
    :data:`EXIT_RHO_HAT`; NaN if never), ``E_int_max_eV`` and
    ``speed_handover_A_per_ps``.
    """
    t = np.asarray(ion.time_ps, dtype=float)
    r = np.sqrt(np.asarray(ion.positions_x) ** 2 + np.asarray(ion.positions_y) ** 2
                + np.asarray(ion.positions_z) ** 2)
    radii = np.asarray(ion.droplet_radii_angstrom, dtype=float).reshape(-1, 1)
    rho = rho_he_ratio(r - radii, steepness=drag_gate_steepness(cfg))
    exposure = np.trapezoid(rho, t, axis=1)
    inside = rho >= EXIT_RHO_HAT
    last_in = np.where(inside.any(axis=1),
                       inside.shape[1] - 1 - np.argmax(inside[:, ::-1], axis=1), -1)
    t_exit = np.where(last_in >= 0, t[np.clip(last_in, 0, None)], np.nan)
    n_shell = np.asarray(ion.n_shell, dtype=float)
    speed = np.sqrt(np.asarray(ion.velocities_x)[:, -1] ** 2
                    + np.asarray(ion.velocities_y)[:, -1] ** 2
                    + np.asarray(ion.velocities_z)[:, -1] ** 2)
    return {"n_handover": n_shell[:, -1], "n_max": n_shell.max(axis=1),
            "exposure_ps": exposure, "t_exit_ps": t_exit,
            "E_int_max_eV": np.asarray(ion.E_int_eV, dtype=float).max(axis=1),
            "speed_handover_A_per_ps": speed}


def _run_is_complete(run_dir: Path) -> bool:
    return all((run_dir / n).exists() for n in _REQUIRED_ARTIFACTS)


def run_law(law: str) -> str:
    from i2_helium_md.simulation.detection_stage import run_detection_stage
    from i2_helium_md.simulation.ion import run_ion_propagation
    from i2_helium_md.simulation.neutral import run_neutral_propagation
    from i2_helium_md.simulation.run_directory import RunDirectory

    verify(law)
    cfg = build(law)
    steps, stride = verify_handover_alignment(cfg)
    run_dir = run_dir_for(law)
    if _run_is_complete(run_dir):
        print(f"[{law}] skip (complete) -> {run_dir.name}", flush=True)
        return law
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    print(f"[{law}] bare start, seed {cfg.seed}, {steps} steps (stride {stride}) "
          f"-> {run_dir.name}", flush=True)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    ref = RunDirectory(parent_run_dir(law)).load_neutral()
    for f in ("positions_x", "positions_y", "positions_z", "droplet_radii"):
        if not np.array_equal(getattr(neutral, f), getattr(ref, f)):
            raise AssertionError(f"[{law}] neutral {f} differs from the parent's")
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False,
                              max_bytes=MAX_BYTES_ION)
    t_last = float(np.asarray(ion.time_ps)[-1])
    if not math.isclose(t_last, cfg.ion_simulation_time - 0.005, abs_tol=1e-9):
        raise AssertionError(f"[{law}] last stored column at {t_last} ps")
    if not np.all(np.asarray(ion.n_shell)[:, 0] == 0.0):
        raise AssertionError(f"[{law}] t0 shell is not bare")
    det = run_detection_stage(ion, cfg, save_path=run.root / "detection.npz")
    np.savez(run.root / SUMMARY_NAME, **ion_summary(ion, cfg))
    del ion, neutral
    n_scored = int(np.asarray(det.detected_mask).sum())
    if n_scored <= 0:
        raise AssertionError(f"[{law}] no scored ions; ion.npz kept")
    ion_path = run.root / "ion.npz"
    (run.root / ION_DROPPED_MARKER).write_text(json.dumps(
        {"ion_npz_bytes": ion_path.stat().st_size, "detection_scored": n_scored,
         "reason": "G8: disk (T: full); scorer reads cfg.json + detection.npz + "
                   f"{SUMMARY_NAME}"}, indent=2), encoding="utf-8")
    ion_path.unlink()
    print(f"[{law}] done: scored {n_scored}/{det.detected_mask.size}", flush=True)
    return law


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("laws", nargs="*", help=f"subset of {LAWS}")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    jobs = list(args.laws) or list(LAWS)
    unknown = sorted(set(jobs) - set(LAWS))
    if unknown:
        raise SystemExit(f"unknown law(s) {unknown}")
    for law in jobs:
        diff = verify(law)
        steps, stride = verify_handover_alignment(build(law))
        print(f"[{law}] guards PASSED: parent = stored cfg.json; diff {diff}; "
              f"{steps} steps, stride {stride} -> {run_dir_for(law).name}", flush=True)
    if args.dry_run:
        return 0
    from concurrent.futures import ProcessPoolExecutor, as_completed
    failures = 0
    with ProcessPoolExecutor(max_workers=min(MAX_CONCURRENCY, len(jobs))) as poolx:
        futures = {poolx.submit(run_law, j): j for j in jobs}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as exc:  # noqa: BLE001
                failures += 1
                print(f"[{futures[fut]}] FAILED: {exc!r}", flush=True)
    print(f"{len(jobs) - failures}/{len(jobs)} runs complete.", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
