"""Generate the Tier-2 staircase capability probe run matrix (pre-F5).

The probe (``docs/drag_port/Tier2/TIER2_STAIRCASE_PROBE_PLAN.md``) answers, at
small N and before any F5 production spend, whether *any* in-band
(kappa, picture, tau) combination reproduces the anchored 9 A shell evolution
21->19->14 -- and, more broadly, what range of post-relaxation terminal shells
the biphasic mechanism can express. The Phase-D bridge showed the single pinned
point sheds ~0.7 He vs 7 (kinetic bottleneck, the s-1 = 59 RRK exponent;
``TIER2_PHASE_D_BRIDGE_FINDINGS.md`` §2); this extends that existence test
across the full lever space named there.

**Existence probe, reported not auto-adjudicated.** TDDFT is not ground truth
(experiment arbitrates at Tier 2 via the size distribution): a miss everywhere
in-band is the RRK-dof mechanism-level OQ to surface, never a silent nu/s
retune; a landing region points the N=500 Stage-1/2 campaign, which is what
formally resolves the F5 gate.

The **active** USER SETTINGS are the **Addendum-A s_eff mini-probe** (12 runs):
the full 45-run kappa x picture x tau probe already executed (outcome (b): no
in-band point lands, freeze at n ~ 20 -- the miss is *kinetic*, the s-1 = 59 RRK
exponent), so kappa and picture are pinned to the bridge point and the RRK
effective-dof ``s`` is swept instead (``S_EFF_GRID``), the falsification lever
the mechanism-level OQ named. Restoring the "FULL PROBE" USER-SETTINGS block
(with ``S_EFF_GRID = [None]``) reproduces the delivered 45-run probe
byte-identically.

Design choices (plan §2 + Addendum A):

* **s_eff sweep** ``{None, 30, 20, 12, 8, 5}`` x ``tau {6.55, 16.5}`` at the
  pinned bridge point. ``None`` = the per-n ``s = 3n-3`` control (reproduces the
  bridge freeze); a constant ``s`` raises the RRK rate ``k = nu*(1-x)^(s-1)`` by
  cutting the exponent. The tau arm tests s<->tau separability (first-shed
  *timing* stays gate-open-governed ~ t_x; magnitude is s-governed).
* **f_int = 0.5 (the bridge pin, NOT the F2 Stage-1 floor).** f_int is
  timing-only (it moves t_x, not the cascade budget Sigma(21)); at 0.5,
  t_x ~ 5 ps -- aligned with the anchored first shed. The floor would open the
  gate at t ~ 0 and corrupt every staircase-*timing* comparison.
* **N = 50, bridge seed/window** -- the (kappa=1, statistical_mixture,
  tau=6.55, s_eff=None) grid point IS the Phase-D bridge configuration, so its
  ion stage must reproduce the bridge numbers exactly (a free wiring oracle).
* **Probe namespace** ``_tier2probe_``: disjoint from the F3 campaign glob
  ``*_tier2_*`` in both directions (namespace lock test). A set ``s_eff``
  appends ``_sNN.NN``; the per-n ``None`` appends nothing (tag byte-identity).
* **0.80 eV only** -- the staircase anchor exists only at the 9 A validation
  condition, so the generator refuses any other budget (including the
  sanctioned production 2.70).

The matching scorer is
``scripts/post_processing/tier2_staircase_probe_report.py``.

Usage::

    python scripts/gen_tier2_staircase_probe.py
"""

from __future__ import annotations

from pathlib import Path
import sys


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"                     # the staircase anchor is 9 A-only
VARIANT = "shared_pure_cubic"
N = 50                          # bridge scale: mean n(t) is stable here

# Run parameters inherit the Phase-D bridge exactly (wiring-oracle point).
ION_TIME_PS = 30.0
DT_ION_PS = 0.01
SEED = 20260604                 # single fixed seed across all grid points

# --- probe grid ---------------------------------------------------------------
# ACTIVE: the Addendum-A s_eff mini-probe (12 runs). The 45-run kappa x picture x
# tau probe already executed (outcome (b): no in-band point lands, freeze at
# n ~ 20; drag_migration_log_tier2.md "Staircase probe EXECUTED"), so kappa and
# picture are pinned to the bridge point and the RRK effective-dof s is swept
# instead (the falsification lever the OQ named). To restore the full 45-run
# probe, swap in the "FULL PROBE" values below and set S_EFF_GRID = [None].
KAPPA_GRID = [1.0]                             # bridge pin (kappa proved inverted/capped)
PICTURE_LIST = ["statistical_mixture"]         # bridge pin (picture near-degenerate)
TAU_GRID_PS = [6.55, 16.5]                     # geometric mid + slow-cooling band edge
S_EFF_GRID = [None, 30, 20, 12, 8, 5]          # None = per-n s=3n-3 control; then constant s

# FULL PROBE (delivered 45-run kappa x picture x tau grid; s_eff = per-n):
#     KAPPA_GRID   = [0.5, 1.0, 2.0, 4.0, 8.0]
#     PICTURE_LIST = ["statistical_mixture", "x2_only", "cooling_relaxed"]
#     TAU_GRID_PS  = [2.6, 6.55, 16.5]
#     S_EFF_GRID   = [None]

# --- pinned knobs (bridge point; f_int deliberately NOT the Stage-1 floor) --
BUDGET_EV = 0.80                # validation only -- see the guard below
F_INT = 0.5                     # bridge pin: t_x ~ 5 ps, timing-faithful
F_RET = 0.1                     # prior (not identifiable at 9 A only)
LAMBDA0_PER_PS = 0.9            # pickup live (bridge central value)

# --- E2 relaxation stage (the "after relaxation" leg of the probe) ----------
RELAXATION_TIME_PS = None       # None -> the experimental 8530 ns cap
RELAXATION_FORCES = "coulomb"   # the two I+ fragments still repel

# --- safety ------------------------------------------------------------------
SKIP_COMPLETED_RUNS = True      # resume: skip a run dir that already has all artifacts
OVERWRITE_EXISTING_RUN = False  # else refuse to clobber a partial/existing run dir


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np  # noqa: E402

from scripts.tier2_common import (  # noqa: E402
    EXPERIMENTAL_RELAXATION_TIME_PS,
    VALIDATION_BUDGET_EV,
    build_biphasic_cfg,
    tier2_probe_run_dir_name,
)
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


# A probe run dir is "complete" (safe to skip on resume) only with all four.
_REQUIRED_ARTIFACTS: tuple[str, ...] = ("cfg.json", "neutral.npz", "ion.npz", "relaxation.npz")


def _run_is_complete(run_dir: Path) -> bool:
    """True when the run dir holds every probe artifact (resume guard)."""
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def probe_grid_points() -> list[tuple[str, float, float, float | None]]:
    """Enumerate ``(picture, kappa, tau_ps, evap_rrk_dof)`` for the probe grid.

    ``evap_rrk_dof`` (the Addendum-A RRK-dof lever) is carried as the fourth
    element; ``None`` = the per-n ``s = 3n-3`` convention (the in-grid control).
    With ``S_EFF_GRID = [None]`` the grid reduces to the delivered
    ``(picture, kappa, tau)`` product (each tuple's fourth element ``None``).
    """
    return [
        (picture, float(kappa), float(tau_ps), s_eff)
        for picture in PICTURE_LIST
        for kappa in KAPPA_GRID
        for tau_ps in TAU_GRID_PS
        for s_eff in S_EFF_GRID
    ]


def build_probe(project_root: Path) -> list[tuple[str, SimConfig, Path]]:
    """Build ``(label, cfg, run_dir)`` for every probe grid point (no propagation).

    Refuses any budget other than the 0.80 eV validation stamp: the anchored
    21->19->14 staircase (the probe's comparator) is read from the 9 A TDDFT
    loss curve at that condition only, so a probe at another budget would be
    scored against an anchor that does not apply.
    """
    if not np.isclose(BUDGET_EV, VALIDATION_BUDGET_EV):
        raise ValueError(
            f"BUDGET_EV={BUDGET_EV} -- the staircase probe runs only at the "
            f"0.80 eV validation budget (the 9 A anchored staircase does not "
            "apply elsewhere; the production sweep is F5's job, not the probe's)."
        )
    relaxation_time_ps = (
        EXPERIMENTAL_RELAXATION_TIME_PS
        if RELAXATION_TIME_PS is None
        else RELAXATION_TIME_PS
    )
    scheduled: list[tuple[str, SimConfig, Path]] = []
    for picture, kappa, tau_ps, s_eff in probe_grid_points():
        cfg = build_biphasic_cfg(
            CASE,
            VARIANT,
            num_molecules=N,
            ion_time_ps=ION_TIME_PS,
            dt_ion_ps=DT_ION_PS,
            seed=SEED,
            lambda0_per_ps=LAMBDA0_PER_PS,
            f_int=F_INT,
            f_ret=F_RET,
            picture=picture,
            kappa=kappa,
            tau_ps=tau_ps,
            evap_rrk_dof=s_eff,
            coulomb_available_eV=BUDGET_EV,
            relaxation_time_ps=relaxation_time_ps,
            relaxation_forces=RELAXATION_FORCES,
        )
        run_dir = project_root / "data" / "runs" / tier2_probe_run_dir_name(
            CASE,
            VARIANT,
            N,
            picture=picture,
            kappa=kappa,
            lambda0_per_ps=LAMBDA0_PER_PS,
            f_int=F_INT,
            f_ret=F_RET,
            tau_ps=tau_ps,
            budget_eV=BUDGET_EV,
            evap_rrk_dof=s_eff,
        )
        s_tag = "per-n" if s_eff is None else f"s={s_eff:.2f}"
        label = (
            f"{CASE} {VARIANT} N={N} probe {picture} "
            f"k={kappa:.2f} tau={tau_ps:.2f} {s_tag}"
        )
        scheduled.append((label, cfg, run_dir))
    return scheduled


def _run_one(label: str, cfg, run_dir: Path) -> None:
    """Write one probe run directory: neutral -> ion -> relaxation (E2)."""
    if run_dir.exists():
        # Resume: a fully-written run dir is skipped so a re-run recovers only
        # the points that failed. A *partial* dir (crash mid-run) is not
        # complete, so it still trips the overwrite guard.
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{label}] skip (complete) -> {run_dir}")
            return
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(
                f"{run_dir} already exists and is not a complete run; set "
                "OVERWRITE_EXISTING_RUN=True to regenerate (probe runs are "
                "scored by the generating code version -- regeneration is the "
                "recovery path, not migration)."
            )
    cfg.validate()
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)

    print(
        f"[{label}] budget={cfg.coulomb_available_eV:.2f} eV "
        f"picture={cfg.ladder_electronic_picture} kappa={cfg.ladder_steepness:.2f} "
        f"tau={cfg.internal_energy_cooling_tau_ps:.2f} ps "
        f"f_int={cfg.internal_energy_partition_fraction:.2f} "
        f"f_ret={cfg.internal_energy_retained_fraction:.2f} "
        f"lambda0={cfg.pickup_rate_coefficient:.2f}/ps "
        f"s_eff={cfg.evap_rrk_dof if cfg.evap_rrk_dof is not None else 'per-n'}"
    )
    print(f"[{label}] neutral propagation ...")
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)

    print(f"[{label}] ion propagation ...")
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)

    print(f"[{label}] relaxation stage (E2) ...")
    relax = run_relaxation_stage(ion, cfg, save_path=run.root / "relaxation.npz")

    n_mean = float(relax.terminal_n.mean())
    print(
        f"[{label}] done -> {run_dir} "
        f"(ion {ion.time_ps.size} steps; relaxed {relax.time_relaxed_ps:.2f} ps, "
        f"terminal n_mean={n_mean:.2f})"
    )


def main() -> int:
    """Generate every grid point of the staircase probe (USER SETTINGS)."""
    scheduled = build_probe(PROJECT_ROOT)
    print(
        f"Tier-2 staircase probe: {len(scheduled)} grid point(s), "
        f"budget={BUDGET_EV} eV, N={N} (existence probe -- reported, "
        "not auto-adjudicated)."
    )
    for label, cfg, run_dir in scheduled:
        _run_one(label, cfg, run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
