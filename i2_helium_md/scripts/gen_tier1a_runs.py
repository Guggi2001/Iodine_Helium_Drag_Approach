"""Generate the Tier-1a fixed/null + anchored-discrete t* sweep.

Produces four self-describing run directories under ``data/runs``:

* fixed null;
* ``anchored_discrete`` with ``t*`` = 0.5 ps;
* ``anchored_discrete`` with ``t*`` = 5.0 ps;
* ``anchored_discrete`` with ``t*`` = 9.0 ps.

This is the run-production half of the Tier-1a section-5/9 deliverable. The matching
scorer is ``scripts/post_processing/tier1a_rmse_table.py``.

Usage::

    python scripts/gen_tier1a_runs.py
"""

from __future__ import annotations

from pathlib import Path
import sys


# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"
VARIANT = "shared_pure_cubic"
N = 50

ION_TIME_PS = 30.0
DT_ION_PS = 0.01
SEED = 20260604
T_STAR_VALUES_PS = [0.5, 5.0, 9.0]

COEFF_OVERRIDES: dict[str, float] = {}
E_BIND_OVERRIDE = None


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tier0_common import build_drag_cfg  # noqa: E402
from scripts.tier1a_common import (  # noqa: E402
    build_anchored_cfg,
    tier1a_run_dir_name,
)
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _run_one(label: str, cfg, run_dir: Path) -> None:
    """Write one Tier-1a run directory using the standard pipeline."""
    cfg.validate()
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)

    print(
        f"[{label}] form={cfg.drag_form} "
        f"coeffs={dict(cfg.drag_coefficients.coefficients)} "
        f"E_bind={cfg.binding_energy_I_ion_eV:.4f} eV "
        f"scenario={cfg.mass_scenario} t*={cfg.t_star_ps:.3f} ps"
    )
    print(f"[{label}] neutral propagation ...")
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)

    print(f"[{label}] ion propagation ...")
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)

    print(
        f"[{label}] done -> {run_dir} "
        f"(neutral {neutral.time_ps.size} steps, ion {ion.time_ps.size} steps)"
    )


def main() -> int:
    """Generate the fixed null and anchored-discrete t* sweep."""
    fixed_cfg = build_drag_cfg(
        CASE,
        VARIANT,
        num_molecules=N,
        ion_time_ps=ION_TIME_PS,
        dt_ion_ps=DT_ION_PS,
        seed=SEED,
        coeff_overrides=COEFF_OVERRIDES,
        e_bind_override=E_BIND_OVERRIDE,
    )
    fixed_dir = PROJECT_ROOT / "data" / "runs" / tier1a_run_dir_name(
        CASE, VARIANT, N, "fixed", None
    )
    _run_one(f"{CASE} {VARIANT} N={N} fixed", fixed_cfg, fixed_dir)

    for t_star_ps in T_STAR_VALUES_PS:
        cfg = build_anchored_cfg(
            CASE,
            VARIANT,
            t_star_ps=t_star_ps,
            num_molecules=N,
            ion_time_ps=ION_TIME_PS,
            dt_ion_ps=DT_ION_PS,
            seed=SEED,
            coeff_overrides=COEFF_OVERRIDES,
            e_bind_override=E_BIND_OVERRIDE,
        )
        run_dir = PROJECT_ROOT / "data" / "runs" / tier1a_run_dir_name(
            CASE, VARIANT, N, "anchored_discrete", t_star_ps
        )
        _run_one(
            f"{CASE} {VARIANT} N={N} anchored t*={t_star_ps:.1f}",
            cfg,
            run_dir,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
