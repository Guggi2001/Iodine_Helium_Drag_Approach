"""Generate the Tier-0 drag runs (ephemeral; data/runs is gitignored).

Produces full-duration (20 ps, dt=0.01) drag runs at a chosen N for the 9A and
18A cases, used to set + confirm the Tier-0 acceptance thresholds. Not a
committed artifact -- the small mean-series CSV exported by
``tier0_drag_comparison.py`` is the committed regression reference.

Usage::

    python scripts/gen_tier0_runs.py <case> <N>
    # case in {9A, 18A}; N e.g. 50 or 2000
"""

from __future__ import annotations

from pathlib import Path
import sys

# =============================================================================
# USER SETTINGS
# =============================================================================

CASE = "9A"      # choose: "9A" or "18A"
N = 50          # choose e.g. 50 or 2000

ION_TIME_PS = 30.0
DT_ION_PS = 0.01
SEED = 20260604

# =============================================================================
# PROJECT IMPORT SETUP
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md import (  # noqa: E402
    single_pulse_N2000_drag,
    single_pulse_N2000_18Angst_drag,
)
from i2_helium_md.simulation.ion import run_ion_propagation  # noqa: E402
from i2_helium_md.simulation.neutral import run_neutral_propagation  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402


_BUILDERS = {
    "9A": single_pulse_N2000_drag,
    "18A": single_pulse_N2000_18Angst_drag,
}

# =============================================================================
# MAIN SCRIPT
# =============================================================================

if CASE not in _BUILDERS:
    raise ValueError(f"CASE must be one of {sorted(_BUILDERS)}, got {CASE!r}")

cfg = _BUILDERS[CASE](
    num_molecules=N,
    ion_simulation_time=ION_TIME_PS,
    dt_ion=DT_ION_PS,
    seed=SEED,
)

cfg.validate()

run_dir = PROJECT_ROOT / "data" / "runs" / f"{CASE}_drag_tier0_N{N}"
run = RunDirectory(run_dir)
run.save_cfg(cfg)

print(f"[{CASE} N={N}] neutral propagation ...")
neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)

print(f"[{CASE} N={N}] ion propagation ...")
ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)

print(
    f"[{CASE} N={N}] done -> {run_dir}  "
    f"(neutral {neutral.time_ps.size} steps, ion {ion.time_ps.size} steps)"
)