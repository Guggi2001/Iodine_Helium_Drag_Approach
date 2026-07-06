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

**Staged operation (MASS §6.11 / plan §2; Stage 1 re-scoped 2026-07-06).** This
is a *parameterized* generator, not an auto-adjudicator. The original Stage-1
kappa x picture co-fit was **superseded by the s_eff promotion** (Phase F plan,
F2 re-scope NB; ``drag_migration_log_tier2.md`` "s_eff PROMOTED"): the pre-F5
staircase probe proved both knobs near-flat on the staircase (kappa inverted +
normalisation-capped for 21->14; picture magnitude-degenerate) and landed the
anchored staircase at a constant ``s_eff ~ 8`` (Bounded, band ~[5, 20]). Stage 1
is therefore the **s_eff x tau co-fit** (``S_EFF_GRID`` x ``TAU_GRID_PS``, with
the per-n ``None`` arm as the in-grid classical-limit control), kappa/picture
pinned at the bridge point, and ``f_int`` pinned at the 0.5 landing prior
(timing-degenerate with the picture via Sigma(21) -> t_x; reported, not fit).
The operator runs the grid, scores it with F3/F4, and hand-pins winners for any
follow-up sweep. The 0.80 eV validation budget runs first; F5 flips
``BUDGET_EV`` to 2.70 only after it lands. tau-sensitivity (the old Stage 2) is
read from the in-grid tau dimension by F4 -- no separate stage needed.

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

# --- Stage-1 co-fit grid (re-scoped 2026-07-06: s_eff x tau) -----------------
# kappa/picture are pinned at the bridge point -- both proven near-flat on the
# staircase (kappa inverted/normalisation-capped, picture magnitude-degenerate;
# probe + cross-check, drag_migration_log_tier2.md). The promoted Bounded s_eff
# (band ~[5, 20], staircase landing [8, 12]) x the tau band is the live co-fit;
# the per-n None arm is the in-grid classical-limit control.
KAPPA_GRID = [1.0]                             # bridge pin
PICTURE_LIST = ["statistical_mixture"]         # bridge pin
TAU_GRID_PS = [2.6, 6.55, 16.5]                # Bounded band [2.6, 16.5] ps
S_EFF_GRID = [None, 5.0, 8.0, 12.0, 16.0, 20.0]  # None = per-n s=3n-3 control

# ORIGINAL STAGE 1 (superseded 2026-07-06; the kappa x picture co-fit at fixed
# tau -- restore by swapping these in and setting S_EFF_GRID = [None]):
#     KAPPA_GRID   = [0.5, 1.0, 2.0, 4.0, 8.0]
#     PICTURE_LIST = ["statistical_mixture", "x2_only", "cooling_relaxed"]
#     TAU_GRID_PS  = [6.55]

# --- pinned / bounded knobs --------------------------------------------------
BUDGET_EV = 0.80                # validation first; F5 flips to 2.70 for production
F_INT = 0.5                     # landing-prior pin (t_x ~ 5 ps; the s_eff ~ 8 landing
                                # was established here); None -> per-point floor
F_RET = 0.1                     # prior (not identifiable at 9 A only)
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


def campaign_grid_points() -> list[tuple[str, float, float, float | None, float]]:
    """Enumerate ``(picture, kappa, tau_ps, s_eff, f_int)`` for the staged grid.

    ``s_eff`` is the promoted constant RRK effective-dof knob (2026-07-06);
    ``None`` = the per-n ``s = 3n-3`` classical-limit control. ``f_int`` is
    resolved per point (literal ``F_INT`` or the picture/kappa floor).
    """
    return [
        (
            picture,
            float(kappa),
            float(tau_ps),
            None if s_eff is None else float(s_eff),
            resolve_f_int(picture, float(kappa)),
        )
        for picture in PICTURE_LIST
        for kappa in KAPPA_GRID
        for tau_ps in TAU_GRID_PS
        for s_eff in S_EFF_GRID
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
    for picture, kappa, tau_ps, s_eff, f_int in campaign_grid_points():
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
            tau_ps=tau_ps,
            evap_rrk_dof=s_eff,
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
            tau_ps=tau_ps,
            budget_eV=BUDGET_EV,
            evap_rrk_dof=s_eff,
            total_strip=TOTAL_STRIP,
        )
        s_tag = "per-n" if s_eff is None else f"s={s_eff:.2f}"
        label = (
            f"{CASE} {VARIANT} N={N} biphasic {picture} "
            f"k={kappa:.2f} tau={tau_ps:.2f} {s_tag}"
        )
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
    """Generate every grid point of the staged campaign (USER SETTINGS)."""
    scheduled = build_campaign(PROJECT_ROOT)
    print(f"Tier-2 campaign: {len(scheduled)} grid point(s), budget={BUDGET_EV} eV.")
    for label, cfg, run_dir in scheduled:
        _run_one(label, cfg, run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
