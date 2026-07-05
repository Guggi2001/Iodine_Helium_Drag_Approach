"""Generate the Tier-2 staged calibration campaign run matrix (Slice F2).

Produces one self-describing run directory per grid point under ``data/runs``,
each carrying ``cfg.json`` / ``neutral.npz`` / ``ion.npz`` / ``relaxation.npz``
(the E2 post-ejection relaxation checkpoint at the experimental timescale). The
matching scoreboard is ``scripts/post_processing/tier2_size_distribution_table.py``
(Slice F3).

**Campaign scope (user, 2026-07-03): the 9 A case only.** There is no reference
shell evolution for the 18 A droplet, so the plan's 9/18 A density-contrast route
to ``f_ret`` is unavailable and ``f_ret`` is held at a prior (not swept). Recorded
in ``docs/drag_port/Tier2/drag_migration_log_tier2.md`` (Phase F entry).

**Staged operation (MASS §6.11 / plan §2).** This is a *parameterized* generator,
not an auto-adjudicator: the operator runs Stage 1 (the ``KAPPA_GRID`` x
``PICTURE_LIST`` co-fit at ``F_INT=None`` -> the per-point self-unbound floor),
scores it with F3/F4, hand-pins the winning ``kappa``/``picture`` into the USER
SETTINGS, and re-runs for the Stage-2 ``tau`` sweep (set ``KAPPA_GRID``/
``PICTURE_LIST`` to the singleton winner and enumerate ``TAU_PS``). The 0.80 eV
validation budget runs first; F5 flips ``BUDGET_EV`` to 2.70 only after it lands.

Bridge finding (``TIER2_PHASE_D_BRIDGE_FINDINGS.md`` §2): the pinned point
under-sheds ~10x and the miss is *kinetic* (the s-1 = 59 RRK exponent), so the
``KAPPA_GRID`` deliberately spans the sharp-cliff end. If no kappa/picture/tau in
the bands lands the 21->14 staircase, that is the mechanism-level RRK-dof open
question the findings name -- surface it, do not silently retune nu or s.

Usage::

    python scripts/gen_tier2_runs.py
"""

from __future__ import annotations

from pathlib import Path
import sys


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"                     # 9 A only (no 18 A reference shell evolution)
VARIANT = "shared_pure_cubic"
N = 500                         # larger-N single-seed for a stable size histogram

# Run parameters inherit Tier-1a/the bridge exactly.
ION_TIME_PS = 30.0
DT_ION_PS = 0.01
SEED = 20260604                 # single fixed seed across all grid points

# --- Stage-1 co-fit grid (the two genuinely-free knobs) ----------------------
KAPPA_GRID = [0.5, 1.0, 2.0, 4.0, 8.0]         # gradual -> sharp cliff
PICTURE_LIST = ["statistical_mixture", "x2_only", "cooling_relaxed"]

# --- pinned / bounded knobs --------------------------------------------------
BUDGET_EV = 0.80                # validation first; F5 flips to 2.70 for production
F_INT = None                    # None -> per-point self-unbound floor (Stage-1 timing pin)
F_RET = 0.1                     # prior (not identifiable at 9 A only)
TAU_PS = 6.55                   # Stage-1 fixed mid; enumerate for the Stage-2 sweep
LAMBDA0_PER_PS = 0.9            # pickup live (bridge central value)

# --- E2 relaxation stage -----------------------------------------------------
RELAXATION_TIME_PS = None       # None -> the experimental 8530 ns cap (set below)
RELAXATION_FORCES = "coulomb"   # the two I+ fragments still repel post-dissociation

# --- regime-axis + safety ----------------------------------------------------
TOTAL_STRIP = False             # F5 secondary regime-axis variant (tag only here)
SKIP_COMPLETED_RUNS = True      # resume: skip a run dir that already has all artifacts
OVERWRITE_EXISTING_RUN = False  # else refuse to clobber a partial/existing run dir


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tier2_common import (  # noqa: E402
    EXPERIMENTAL_RELAXATION_TIME_PS,
    PRODUCTION_BUDGET_EV,
    VALIDATION_BUDGET_EV,
    build_biphasic_cfg,
    tier2_run_dir_name,
)
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.physics.internal_energy_budget import f_int_floor  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.relaxation_stage import run_relaxation_stage  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


# Sanctioned per-ion Coulomb budgets: 0.80 eV validation / 2.70 eV production.
# The config field carries "NO hard refuse" by design (plan §8, a provenance
# stamp) -- so the generator guards the operator setting instead, catching a
# typo (e.g. 8.0 for 0.80) before it silently poisons every run's RRK gate.
_SANCTIONED_BUDGETS_EV: tuple[float, ...] = (VALIDATION_BUDGET_EV, PRODUCTION_BUDGET_EV)

# A campaign run dir is "complete" (safe to skip on resume) only with all four.
_REQUIRED_ARTIFACTS: tuple[str, ...] = ("cfg.json", "neutral.npz", "ion.npz", "relaxation.npz")


def _run_is_complete(run_dir: Path) -> bool:
    """True when the run dir holds every campaign artifact (resume guard)."""
    return all((run_dir / name).exists() for name in _REQUIRED_ARTIFACTS)


def resolve_f_int(picture: str, kappa: float) -> float:
    """Return the grid point's ``f_int``: the USER ``F_INT`` or the floor.

    ``F_INT=None`` pins ``f_int`` at the picture/kappa-keyed self-unbound floor
    ``Sigma(n*)/E_avail`` (the Stage-1 timing pin: gate open from onset, so the
    co-fit isolates kappa x picture from the shed *timing*; bridge findings §2.3).
    A literal ``F_INT`` overrides the floor for a conscious timing comparison.
    """
    if F_INT is not None:
        return float(F_INT)
    return f_int_floor(e_avail_eV=BUDGET_EV, picture=picture, kappa=kappa)


def campaign_grid_points() -> list[tuple[str, float, float]]:
    """Enumerate ``(picture, kappa, f_int)`` for the staged grid (USER SETTINGS)."""
    return [
        (picture, float(kappa), resolve_f_int(picture, float(kappa)))
        for picture in PICTURE_LIST
        for kappa in KAPPA_GRID
    ]


def build_campaign(project_root: Path) -> list[tuple[str, SimConfig, Path]]:
    """Build ``(label, cfg, run_dir)`` for every grid point (no propagation)."""
    if BUDGET_EV not in _SANCTIONED_BUDGETS_EV:
        raise ValueError(
            f"BUDGET_EV={BUDGET_EV} is not a sanctioned budget "
            f"{_SANCTIONED_BUDGETS_EV} (0.80 validation / 2.70 production). "
            "Edit this guard consciously for an exploratory budget."
        )
    relaxation_time_ps = (
        EXPERIMENTAL_RELAXATION_TIME_PS
        if RELAXATION_TIME_PS is None
        else RELAXATION_TIME_PS
    )
    scheduled: list[tuple[str, object, Path]] = []
    for picture, kappa, f_int in campaign_grid_points():
        cfg = build_biphasic_cfg(
            CASE,
            VARIANT,
            num_molecules=N,
            ion_time_ps=ION_TIME_PS,
            dt_ion_ps=DT_ION_PS,
            seed=SEED,
            lambda0_per_ps=LAMBDA0_PER_PS,
            f_int=f_int,
            f_ret=F_RET,
            picture=picture,
            kappa=kappa,
            tau_ps=TAU_PS,
            coulomb_available_eV=BUDGET_EV,
            relaxation_time_ps=relaxation_time_ps,
            relaxation_forces=RELAXATION_FORCES,
        )
        run_dir = project_root / "data" / "runs" / tier2_run_dir_name(
            CASE,
            VARIANT,
            N,
            picture=picture,
            kappa=kappa,
            lambda0_per_ps=LAMBDA0_PER_PS,
            f_int=f_int,
            f_ret=F_RET,
            tau_ps=TAU_PS,
            budget_eV=BUDGET_EV,
            total_strip=TOTAL_STRIP,
        )
        label = f"{CASE} {VARIANT} N={N} biphasic {picture} k={kappa:.2f}"
        scheduled.append((label, cfg, run_dir))
    return scheduled


def _run_one(label: str, cfg, run_dir: Path) -> None:
    """Write one campaign run directory: neutral -> ion -> relaxation (E2)."""
    if run_dir.exists():
        # Resume: a fully-written run dir is skipped so a re-run recovers only the
        # points that failed. A *partial* dir (crash mid-run) is not complete, so
        # it still trips the overwrite guard rather than being silently kept.
        if SKIP_COMPLETED_RUNS and _run_is_complete(run_dir):
            print(f"[{label}] skip (complete) -> {run_dir}")
            return
        if not OVERWRITE_EXISTING_RUN:
            raise FileExistsError(
                f"{run_dir} already exists and is not a complete run; set "
                "OVERWRITE_EXISTING_RUN=True to regenerate (campaign runs are "
                "scored by the generating code version -- regeneration is the "
                "recovery path, not migration)."
            )
    cfg.validate()
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)

    print(
        f"[{label}] form={cfg.drag_form} "
        f"E_bind={cfg.binding_energy_I_ion_eV:.4f} eV "
        f"budget={cfg.coulomb_available_eV:.2f} eV "
        f"picture={cfg.ladder_electronic_picture} kappa={cfg.ladder_steepness:.2f} "
        f"f_int={cfg.internal_energy_partition_fraction:.3f} "
        f"f_ret={cfg.internal_energy_retained_fraction:.2f} "
        f"tau={cfg.internal_energy_cooling_tau_ps:.2f} ps "
        f"lambda0={cfg.pickup_rate_coefficient:.2f}/ps"
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
    """Generate every grid point of the staged campaign (USER SETTINGS)."""
    scheduled = build_campaign(PROJECT_ROOT)
    print(f"Tier-2 campaign: {len(scheduled)} grid point(s), budget={BUDGET_EV} eV.")
    for label, cfg, run_dir in scheduled:
        _run_one(label, cfg, run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
