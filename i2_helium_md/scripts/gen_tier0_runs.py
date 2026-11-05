"""Generate one Tier-0 drag run of a thesis law.

A full-duration fixed-mass run for one droplet geometry (9 A or 18 A), written to
``data/runs`` (gitignored). The law is ``LQ`` (B2), ``LQ_eff`` (B3) or ``CC``
(h405), or with ``RIDGE_C`` an m = 21 calibration-ridge point. Score it with
``scripts/post_processing/tier0_drag_comparison.py``.

Usage::

    # edit USER SETTINGS, then:
    python scripts/gen_tier0_runs.py
"""

from __future__ import annotations

from pathlib import Path
import sys

# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"      # droplet geometry: "9A" or "18A"
LAW = "LQ"       # "LQ" | "LQ_eff" | "CC"
RIDGE_C = None   # c [amu/A] of a calibration_ridge point (LAW = "LQ"), or None

N = 50           # ensemble size
ION_TIME_PS = 30.0
DT_ION_PS = 0.01
SEED = 20260604

E_BIND_OVERRIDE = None                     # eV, or None for the law's own binding
M_EFF_AMU = 126.90 + 21 * 4.0026           # fixed ion mass I+ + 21 He (G11.2) [amu]

# Marker appended to the run-dir name; set the same value in
# tier0_drag_comparison.py to score the run.
RUN_TAG = "B2m21"

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from scripts.tier0_common import build_tier0_cfg, run_dir_name  # noqa: E402


def run_tier0(cfg: SimConfig, run_dir: Path) -> None:
    """Run the neutral and ion stages of ``cfg`` into ``run_dir``.

    Writes ``cfg.json``, ``neutral.npz`` and ``ion.npz``, replacing any earlier
    run in that directory.
    """
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)


def main() -> int:
    cfg = build_tier0_cfg(
        CASE, LAW, num_molecules=N, ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        seed=SEED, m_eff_amu=M_EFF_AMU, e_bind_override=E_BIND_OVERRIDE,
        ridge_c=RIDGE_C,
    )
    run_dir = PROJECT_ROOT / "data" / "runs" / run_dir_name(CASE, LAW, N, RUN_TAG)
    print(f"[{CASE} {LAW} N={N}] form={cfg.drag_form} "
          f"coeffs={dict(cfg.drag_coefficients.coefficients)} "
          f"E_bind={cfg.binding_energy_I_ion_eV:.4f} eV m={cfg.m_eff_amu:.4f} amu")
    run_tier0(cfg, run_dir)
    print(f"[{CASE} {LAW} N={N}] done -> {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
