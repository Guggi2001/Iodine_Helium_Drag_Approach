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

# Deterministic Tier-0 envelope: full duration, full dt, fixed seed, noise off.
ION_TIME_PS = 20.0
DT_ION_PS = 0.01
SEED = 20260604


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    case = sys.argv[1]
    n = int(sys.argv[2])
    if case not in _BUILDERS:
        raise SystemExit(f"case must be one of {sorted(_BUILDERS)}, got {case!r}")

    cfg = _BUILDERS[case](
        num_molecules=n,
        ion_simulation_time=ION_TIME_PS,
        dt_ion=DT_ION_PS,
        seed=SEED,
    )
    cfg.validate()

    run_dir = PROJECT_ROOT / "data" / "runs" / f"{case}_drag_tier0_N{n}"
    run = RunDirectory(run_dir)
    run.save_cfg(cfg)

    print(f"[{case} N={n}] neutral propagation ...")
    neutral = run_neutral_propagation(cfg, run_dir=run, verbose=False)
    print(f"[{case} N={n}] ion propagation ...")
    ion = run_ion_propagation(cfg, neutral, run_dir=run, verbose=False)
    print(
        f"[{case} N={n}] done -> {run_dir}  "
        f"(neutral {neutral.time_ps.size} steps, ion {ion.time_ps.size} steps)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
