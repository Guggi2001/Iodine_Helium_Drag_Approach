"""Generate Tier-0 drag runs: the five thesis runs, or one run of your own.

Full-duration fixed-mass runs, written to ``data/runs`` (gitignored). A finished run
with the same model is skipped; a partial one or one with another model is refused.

* ``RUN_SET = "thesis"``: the runs behind the thesis figures
  (``scripts/tier0_runs.py``). The settings below it are not used.
* ``RUN_SET = "custom"``: one run for one droplet geometry (9 A or 18 A). The law is
  ``LQ`` (B2), ``LQ_eff`` (B3) or ``CC`` (h405), or with ``RIDGE_C`` an m = 21
  calibration-ridge point. Score it with
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

RUN_SET = "thesis"   # "thesis" | "custom"

# --- custom run ---
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

from i2_helium_md.simulation.battery import run_neutral_and_ion  # noqa: E402
from i2_helium_md.simulation.thesis_runs import RUNS_ROOT  # noqa: E402
from scripts.tier0_common import build_tier0_cfg, run_dir_name  # noqa: E402
from scripts.tier0_runs import tier0_runs  # noqa: E402

OUTPUT_ROOT = RUNS_ROOT


def main() -> int:
    """Make the selected runs."""
    if RUN_SET not in ("thesis", "custom"):
        raise SystemExit(f"RUN_SET must be 'thesis' or 'custom', got {RUN_SET!r}")
    if RUN_SET == "thesis":
        for name, cfg in tier0_runs().items():
            run_neutral_and_ion(cfg, OUTPUT_ROOT / name)
        return 0
    cfg = build_tier0_cfg(
        CASE, LAW, num_molecules=N, ion_time_ps=ION_TIME_PS, dt_ion_ps=DT_ION_PS,
        seed=SEED, m_eff_amu=M_EFF_AMU, e_bind_override=E_BIND_OVERRIDE,
        ridge_c=RIDGE_C,
    )
    run_dir = OUTPUT_ROOT / run_dir_name(CASE, LAW, N, RUN_TAG)
    print(f"[{CASE} {LAW} N={N}] form={cfg.drag_form} "
          f"coeffs={dict(cfg.drag_coefficients.coefficients)} "
          f"E_bind={cfg.binding_energy_I_ion_eV:.4f} eV m={cfg.m_eff_amu:.4f} amu")
    run_neutral_and_ion(cfg, run_dir)
    print(f"[{CASE} {LAW} N={N}] -> {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
