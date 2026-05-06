"""Tests for scripts/run_single_pulse.py."""

from __future__ import annotations

import pytest

from scripts.run_single_pulse import main
from i2_helium_md.simulation.run_directory import RunDirectory


def test_tiny_single_pulse_run_writes_run_directory(tmp_path):
    run_path = tmp_path / "single_pulse"

    exit_code = main(
        [
            "--run-dir",
            str(run_path),
            "--num-molecules",
            "2",
            "--seed",
            "123",
            "--ion-simulation-time",
            "0.02",
        ]
    )

    assert exit_code == 0

    run = RunDirectory(run_path)
    assert run.has_cfg()
    assert run.has_neutral()
    assert run.has_ion()

    cfg = run.load_cfg()
    assert cfg.num_molecules == 2
    assert cfg.seed == 123
    assert cfg.ion_simulation_time == pytest.approx(0.02)

    neutral = run.load_neutral()
    ion = run.load_ion()
    assert neutral.num_molecules == 2
    assert ion.num_molecules == 2
    assert neutral.time_ps.size == 2
    assert ion.time_ps.size == 2


def test_existing_outputs_require_force(tmp_path):
    run_path = tmp_path / "existing"
    assert main(
        [
            "--run-dir",
            str(run_path),
            "--num-molecules",
            "2",
            "--seed",
            "1",
            "--ion-simulation-time",
            "0.02",
        ]
    ) == 0

    with pytest.raises(SystemExit, match="--force"):
        main(
            [
                "--run-dir",
                str(run_path),
                "--num-molecules",
                "2",
                "--seed",
                "2",
                "--ion-simulation-time",
                "0.02",
            ]
        )


def test_force_allows_existing_outputs(tmp_path):
    run_path = tmp_path / "force"
    base_args = [
        "--run-dir",
        str(run_path),
        "--num-molecules",
        "2",
        "--ion-simulation-time",
        "0.02",
    ]

    assert main([*base_args, "--seed", "1"]) == 0
    assert main([*base_args, "--seed", "2", "--force"]) == 0

    cfg = RunDirectory(run_path).load_cfg()
    assert cfg.seed == 2
