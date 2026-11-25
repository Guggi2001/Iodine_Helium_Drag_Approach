"""Tier-0 drag runs behind the thesis figures: configs and reading.

Runs are made by ``scripts/gen_tier0_runs.py``; the figure scripts only
read them (:func:`finished_tier0_run`). Drag laws: force a*v + c*v^2 with a in
amu/ps, c in amu/Angstrom; bindings in eV.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import NamedTuple

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md import load_drag_bundle  # noqa: E402
from i2_helium_md.config import SimConfig  # noqa: E402
from i2_helium_md.simulation.battery import TWO_STAGE_ARTIFACTS, is_finished  # noqa: E402
from i2_helium_md.simulation.run_directory import RunDirectory  # noqa: E402
from i2_helium_md.simulation.thesis_runs import RUNS_ROOT, require_run  # noqa: E402
from scripts.extraction.method_b_g412_ridge import (  # noqa: E402
    BOUND_DT_ION_PS,
    BOUND_ION_TIME_PS,
    BOUND_N,
    G11_M_EFF_AMU,
    SEED,
    check_point,
    point_cfg,
)
from scripts.tier0_common import M_EFF_AMU, build_tier0_cfg  # noqa: E402

DRIVER = "gen_tier0_runs.py"

# Calibration runs (LQ law, N = 50).
CALIBRATION_N = 50
CALIBRATION_SEED = 20260604
CALIBRATION_DT_ION_PS = 0.01
CALIBRATION_M_EFF_AMU = 126.90 + 21 * 4.0026   # I+ + 21 He (G11.2) [amu]


class CalibrationRun(NamedTuple):
    """One Tier-0 calibration run."""

    name: str                         # run-directory name
    case: str                         # droplet geometry, "9A" or "18A"
    ion_time_ps: float                # ion-stage duration [ps]
    e_bind_eV: float | None = None    # binding [eV]; None = the LQ bundle's

    def cfg(self) -> SimConfig:
        """The config of this run."""
        return build_tier0_cfg(
            self.case, "LQ", num_molecules=CALIBRATION_N, ion_time_ps=self.ion_time_ps,
            dt_ion_ps=CALIBRATION_DT_ION_PS, seed=CALIBRATION_SEED,
            m_eff_amu=CALIBRATION_M_EFF_AMU, e_bind_override=self.e_bind_eV)


STUCK_18A = CalibrationRun("tier0_18A_B2_stuck", "18A", 50.0, e_bind_eV=0.308)
B2_18A = CalibrationRun("tier0_18A_B2", "18A", 30.0)
B2_9A = CalibrationRun("tier0_9A_B2", "9A", 30.0)
CALIBRATION_RUNS = (STUCK_18A, B2_18A, B2_9A)

_B2 = check_point("B2")
_B3 = load_drag_bundle("LQ_eff", expected_m_eff_amu=M_EFF_AMU)
#: lq laws at m = 21 / 3.5 A gate: label -> (a [amu/ps], c [amu/A], E_bind [eV]).
LQ_LAWS = {
    "B2": (_B2["a_amu_per_ps"], _B2["c_amu_per_A"], _B2["e_bind_eV"]),
    "B3": (_B3.coefficients["a"], _B3.coefficients["c"],
           _B3.effective_binding_energy_I_ion_eV),
}


def h405_9a_cfg() -> SimConfig:
    """Radial 9 A Tier-0 config of the re-landed h405 (G4.12 bound settings, m = 21, 3.5 A)."""
    return build_tier0_cfg("9A", "CC", num_molecules=BOUND_N,
                           ion_time_ps=BOUND_ION_TIME_PS, dt_ion_ps=BOUND_DT_ION_PS,
                           seed=SEED, m_eff_amu=G11_M_EFF_AMU)


def bound_run_name(law: str) -> str:
    """Directory name of the 9 A exit-bound run of ``law`` ("B2", "B3" or "h405")."""
    if law != "h405" and law not in LQ_LAWS:
        raise KeyError(f"unknown law {law!r}; expected one of {[*LQ_LAWS, 'h405']}")
    return f"tier0_9A_{law}"


def tier0_runs() -> dict[str, SimConfig]:
    """Tier-0 runs the thesis figures read: directory name -> config.

    The 9 A calibration run of B2 is also its exit-bound run.
    """
    runs = {run.name: run.cfg() for run in CALIBRATION_RUNS}
    a, c, e_bind = LQ_LAWS["B3"]
    runs[bound_run_name("B3")] = point_cfg(c, a, e_bind, mass_amu=G11_M_EFF_AMU)
    runs[bound_run_name("h405")] = h405_9a_cfg()
    return runs


def finished_tier0_run(name: str, *, root: Path = RUNS_ROOT) -> RunDirectory:
    """The finished run ``name``, a key of :func:`tier0_runs`; runs nothing.

    Raises
    ------
    KeyError
        Unknown name; the message lists the choices.
    FileNotFoundError
        The run does not exist; the message names the driver that makes it.
    RunDirectoryMismatch
        The directory is partial or damaged, is being written, or holds another model.
    """
    runs = tier0_runs()
    if name not in runs:
        raise KeyError(f"unknown Tier-0 run {name!r}; expected one of {list(runs)}")
    path = require_run(root / name, DRIVER)
    if not is_finished(path, runs[name], TWO_STAGE_ARTIFACTS):
        raise FileNotFoundError(f"{path} holds no run; make it with scripts/{DRIVER}")
    return RunDirectory(path)
