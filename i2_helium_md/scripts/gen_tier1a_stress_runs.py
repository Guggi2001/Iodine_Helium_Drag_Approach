"""Generate Tier-1a onset-violent stripping stress runs.

Produces one fixed null and four diagnostic onset-strip stress run directories
under ``data/runs``. These are sensitivity tests, not TDDFT-anchored physical
Tier-1a runs.

Usage::

    python scripts/gen_tier1a_stress_runs.py
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
OVERWRITE_EXISTING_RUN = False

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
    TIER1A_STRESS_N_FINAL_VALUES,
    TIER1A_STRESS_T_STRIP_PS,
    build_onset_strip_cfg,
    tier1a_run_dir_name,
    tier1a_stress_run_dir_name,
)
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


_PROTECTED_OUTPUTS = ("cfg.json", "neutral.npz", "ion.npz")


def _raise_if_existing_outputs(run_dir: Path) -> None:
    """Protect production run outputs from accidental overwrite."""
    if OVERWRITE_EXISTING_RUN:
        return

    existing = [name for name in _PROTECTED_OUTPUTS if (run_dir / name).exists()]
    if existing:
        existing_text = ", ".join(existing)
        raise FileExistsError(
            f"Refusing to overwrite existing run directory {run_dir}: "
            f"found {existing_text}. Set OVERWRITE_EXISTING_RUN = True to rerun."
        )


def _run_one(label: str, cfg, run_dir: Path) -> None:
    """Write one run directory using the standard pipeline."""
    _raise_if_existing_outputs(run_dir)
    cfg.validate()
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)

    scenario_text = f"scenario={cfg.mass_scenario}"
    if cfg.mass_scenario != "fixed" and cfg.anchor_mode == "onset_strip":
        scenario_text = (
            f"{scenario_text} anchor_mode={cfg.anchor_mode} "
            f"t_strip={cfg.t_star_ps:.3f} ps"
        )
        if cfg.anchor_n_final is not None:
            scenario_text = f"{scenario_text} n_final={cfg.anchor_n_final}"

    print(
        f"[{label}] form={cfg.drag_form} "
        f"coeffs={dict(cfg.drag_coefficients.coefficients)} "
        f"E_bind={cfg.binding_energy_I_ion_eV:.4f} eV "
        f"{scenario_text}"
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

    for n_final in TIER1A_STRESS_N_FINAL_VALUES:
        cfg = build_onset_strip_cfg(
            CASE,
            VARIANT,
            n_final=n_final,
            t_strip_ps=TIER1A_STRESS_T_STRIP_PS,
            num_molecules=N,
            ion_time_ps=ION_TIME_PS,
            dt_ion_ps=DT_ION_PS,
            seed=SEED,
            coeff_overrides=COEFF_OVERRIDES,
            e_bind_override=E_BIND_OVERRIDE,
        )
        run_dir = PROJECT_ROOT / "data" / "runs" / tier1a_stress_run_dir_name(
            CASE,
            VARIANT,
            N,
            n_final=n_final,
            t_strip_ps=TIER1A_STRESS_T_STRIP_PS,
        )
        _run_one(
            f"{CASE} {VARIANT} N={N} onset-strip n_final={n_final}",
            cfg,
            run_dir,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
