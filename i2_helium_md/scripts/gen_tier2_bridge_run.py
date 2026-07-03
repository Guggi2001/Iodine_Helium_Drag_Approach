"""Generate the Phase-D bridge run (Slice Z): 9 Å / 0.80 eV ``biphasic``.

Produces one self-describing run directory under ``data/runs`` — the
generative counterpart of the Tier-1a anchored runs (same conditions, same
Tier-0 drag bundle, same N/time/dt/seed; only ``mass_scenario`` and the
pinned priored knob values differ; plan §0/§3). The matching report script is
``scripts/post_processing/tier2_bridge_report.py``.

Usage::

    python scripts/gen_tier2_bridge_run.py
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

# Run parameters inherit Tier-1a exactly (gen_tier1a_runs.py).
ION_TIME_PS = 30.0
DT_ION_PS = 0.01
SEED = 20260604

COEFF_OVERRIDES: dict[str, float] = {}
E_BIND_OVERRIDE = None


# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.tier2_common import (  # noqa: E402
    build_biphasic_cfg,
    tier2_bridge_run_dir_name,
)
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


def _run_one(label: str, cfg, run_dir: Path) -> None:
    """Write one bridge run directory using the standard pipeline."""
    cfg.validate()
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)

    print(
        f"[{label}] form={cfg.drag_form} "
        f"coeffs={dict(cfg.drag_coefficients.coefficients)} "
        f"E_bind={cfg.binding_energy_I_ion_eV:.4f} eV "
        f"scenario={cfg.mass_scenario} "
        f"lambda0={cfg.pickup_rate_coefficient:.2f}/ps "
        f"f_int={cfg.internal_energy_partition_fraction:.2f} "
        f"f_ret={cfg.internal_energy_retained_fraction:.2f} "
        f"tau={cfg.internal_energy_cooling_tau_ps:.2f} ps "
        f"kappa={cfg.ladder_steepness:.1f} "
        f"picture={cfg.ladder_electronic_picture}"
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
    """Generate the single pinned-point bridge run."""
    cfg = build_biphasic_cfg(
        CASE,
        VARIANT,
        num_molecules=N,
        ion_time_ps=ION_TIME_PS,
        dt_ion_ps=DT_ION_PS,
        seed=SEED,
        coeff_overrides=COEFF_OVERRIDES,
        e_bind_override=E_BIND_OVERRIDE,
    )
    run_dir = PROJECT_ROOT / "data" / "runs" / tier2_bridge_run_dir_name(
        CASE, VARIANT, N
    )
    _run_one(f"{CASE} {VARIANT} N={N} bridge biphasic", cfg, run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
