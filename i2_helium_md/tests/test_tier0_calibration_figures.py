"""The run table of the five Tier-0 calibration figures (config only, no propagation)."""

from __future__ import annotations

import struct

import pytest

from scripts.post_processing.tier0_calibration_figures import OUT_DIR, RUNS, build_cfg
from scripts.tier0_common import run_dir_name


def _png_size(path):
    with open(path, "rb") as fh:
        return struct.unpack(">II", fh.read(24)[16:24])


def test_the_five_figures():
    names = {r.velocity_png for r in RUNS} | {r.force_png for r in RUNS if r.force_png}
    assert names == {"18A_Stuck_Velocity_R.png", "18A_Stuck_Force_Energy.png",
                     "18A_B2.png", "18A_B2_Force_Energy.png", "9A_B2.png"}


def test_committed_pngs_keep_the_paired_size():
    """Each pair sits side by side in the thesis, so every PNG must be 600 x 450 px."""
    for run in RUNS:
        for name in filter(None, (run.velocity_png, run.force_png)):
            assert _png_size(OUT_DIR / name) == (600, 450)


def test_stuck_and_b2_runs_get_distinct_directories():
    names = [run_dir_name(r.case, "LQ", 50, r.run_tag) for r in RUNS]
    assert len(set(names)) == len(names)


@pytest.mark.parametrize("run", RUNS, ids=lambda r: f"{r.case}_{r.run_tag}")
def test_run_configs(run):
    cfg = build_cfg(run)
    assert cfg.drag_coefficients.coefficients == {"a": 18.007189923460135, "c": 1.6}
    assert cfg.m_eff_amu == 126.90 + 21 * 4.0026
    assert (cfg.num_molecules, cfg.seed) == (50, 20260604)
    assert cfg.ion_simulation_time == run.ion_time_ps
    assert cfg.binding_energy_I_ion_eV == (run.e_bind_eV or 0.11727328616561973)
