"""Coverage for the Tier-1a RMSE-table orchestration layer.

The production generator writes under ``data/runs`` and can be expensive. These
tests stay on small temporary run directories and exercise the reusable helpers
that the scripts call.
"""

from __future__ import annotations

import csv
import math

import numpy as np
import pytest

from i2_helium_md.physics.constants import U
from i2_helium_md.physics.shell_schedule import complex_mass_amu
from i2_helium_md.postprocess import HedftTrajectory, SmoothedSpeedReference
from i2_helium_md.simulation.checkpoint import IonCheckpoint
from i2_helium_md.simulation.run_directory import RunDirectory


def _tiny_ion_checkpoint() -> IonCheckpoint:
    """A deterministic one-molecule v6 checkpoint with simple RMSE oracles."""
    time = np.array([0.0, 1.0, 2.0, 3.0])
    zeros = np.zeros((2, time.size))

    positions_x = zeros.copy()
    positions_y = zeros.copy()
    positions_z = zeros.copy()
    positions_x[0, :] = 0.0
    positions_x[1, :] = np.array([1.0, 2.0, 3.0, 4.0])

    velocities_x = zeros.copy()
    velocities_y = zeros.copy()
    velocities_z = zeros.copy()
    velocities_x[1, :] = np.array([10.0, 20.0, 30.0, 40.0])

    # Four-term ledger residual: system energy is [0, 1, 0, -2] eV per molecule.
    e_kin = zeros.copy()
    e_kin[0, :] = np.array([0.0, 1.0, 0.0, -2.0])

    n_shell = np.array(
        [
            [21, 20, 19, 14],
            [21, 20, 19, 14],
        ],
        dtype=float,
    )

    mass_kg = np.full(2, complex_mass_amu(14) * U)
    mass_history_kg = (126.90 + n_shell * 4.0026) * U

    return IonCheckpoint(
        num_molecules=1,
        time_ps=time,
        positions_x=positions_x,
        positions_y=positions_y,
        positions_z=positions_z,
        velocities_x=velocities_x,
        velocities_y=velocities_y,
        velocities_z=velocities_z,
        positions_final_x=positions_x[:, -1],
        positions_final_y=positions_y[:, -1],
        positions_final_z=positions_z[:, -1],
        velocities_final_x=velocities_x[:, -1],
        velocities_final_y=velocities_y[:, -1],
        velocities_final_z=velocities_z[:, -1],
        mass_kg=mass_kg,
        mass_final_kg=mass_kg,
        mass_history_kg=mass_history_kg,
        droplet_radii_angstrom=np.full(2, 30.0),
        E_kin_eV=e_kin,
        E_pot_eV=zeros.copy(),
        E_dissip_eV=zeros.copy(),
        E_mass_transfer_eV=zeros.copy(),
        n_shell=n_shell,
        b_ion_outside=np.zeros(1, dtype=bool),
        relative_loss_per_ps=zeros.copy(),
        number_of_collisions=np.zeros((2, time.size), dtype=int),
        temperature_diagnostic=np.zeros((time.size, 3)),
        mass_scenario="anchored_discrete",
    )


def _hedft_reference() -> HedftTrajectory:
    time = np.array([0.0, 1.0, 2.0, 3.0])
    distance = np.array([1.0, 2.0, 3.0, 4.0])
    zeros = np.zeros_like(time)
    return HedftTrajectory(
        time_ps=time,
        v1_magnitude_Aps=zeros,
        v2_magnitude_Aps=np.array([10.0, 20.0, 30.0, 40.0]),
        v1_z_Aps=zeros,
        v2_z_Aps=zeros,
        v1_x_Aps=zeros,
        v2_x_Aps=zeros,
        v1_y_Aps=zeros,
        v2_y_Aps=zeros,
        distance_A=distance,
        x1_A=zeros,
        x2_A=distance,
        y1_A=zeros,
        y2_A=zeros,
        z1_A=zeros,
        z2_A=zeros,
        droplet_radius_A=9.0,
        source_path=__file__,
    )


def _smoothed_reference() -> SmoothedSpeedReference:
    return SmoothedSpeedReference(
        time_ps=np.array([1.0, 2.0, 3.0]),
        speed_Aps=np.array([20.0, 30.0, 50.0]),
        source_path=__file__,
    )


def test_build_anchored_cfg_validates_and_sets_tier1a_fields():
    from scripts.tier1a_common import build_anchored_cfg

    cfg = build_anchored_cfg(
        "9A",
        "shared_pure_cubic",
        t_star_ps=9.0,
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=123,
    )

    cfg.validate()
    assert cfg.mass_scenario == "anchored_discrete"
    assert cfg.t_star_ps == 9.0
    assert cfg.anchor_mode == "time"
    assert cfg.coulomb_available_eV == 0.80
    assert cfg.allow_inconsistent_mass_pairing is True
    assert cfg.mass_initial_amu == pytest.approx(complex_mass_amu(21))


def test_build_onset_strip_cfg_sets_stress_fields():
    from scripts.tier1a_common import build_onset_strip_cfg

    cfg = build_onset_strip_cfg(
        "9A",
        "shared_pure_cubic",
        n_final=2,
        t_strip_ps=0.5,
        num_molecules=2,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=123,
    )

    assert cfg.mass_scenario == "anchored_discrete"
    assert cfg.anchor_mode == "onset_strip"
    assert cfg.anchor_n_final == 2
    assert cfg.t_star_ps == 0.5
    assert cfg.mass_initial_amu == pytest.approx(complex_mass_amu(21))
    assert cfg.allow_inconsistent_mass_pairing is True


def test_tier1a_run_names_share_tier0_convention():
    from scripts.tier1a_common import tier1a_run_dir_name, tier1a_run_tag

    assert tier1a_run_tag("fixed", None) == "tier1a_fixed"
    assert tier1a_run_tag("anchored_discrete", 5.0) == "tier1a_anchored_continuous_t5.0"
    assert (
        tier1a_run_dir_name("9A", "shared_pure_cubic", 50, "anchored_discrete", 5.0)
        == "9A_drag_shared_pure_cubic_N50_tier1a_anchored_continuous_t5.0"
    )


def test_tier1a_stress_run_names_are_distinct():
    from scripts.tier1a_common import tier1a_stress_run_dir_name, tier1a_stress_run_tag

    assert (
        tier1a_stress_run_tag(n_final=0, t_strip_ps=0.5)
        == "tier1a_stress_onset_strip_n0_t0.5"
    )
    assert (
        tier1a_stress_run_dir_name("9A", "shared_pure_cubic", 50, n_final=0, t_strip_ps=0.5)
        == "9A_drag_shared_pure_cubic_N50_tier1a_stress_onset_strip_n0_t0.5"
    )


def test_score_tier1a_run_emits_rich_table_row(tmp_path):
    from scripts.post_processing.tier1a_rmse_table import (
        TIER1A_TABLE_COLUMNS,
        score_tier1a_run,
        write_rows_csv,
    )

    run = RunDirectory(tmp_path / "run")
    run.save_ion(_tiny_ion_checkpoint())

    row = score_tier1a_run(
        run,
        case="9A",
        variant="shared_pure_cubic",
        n=1,
        run_tag="tier1a_anchored_continuous_t5.0",
        scenario="anchored_discrete",
        t_star_ps=5.0,
        hedft=_hedft_reference(),
        smoothed=_smoothed_reference(),
        window_start_ps=1.0,
    )

    assert list(row) == TIER1A_TABLE_COLUMNS
    assert row["case"] == "9A"
    assert row["variant"] == "shared_pure_cubic"
    assert row["N"] == 1
    assert row["scenario"] == "anchored_discrete"
    assert row["t_star_ps"] == 5.0
    assert row["window_start_ps"] == 1.0
    assert row["window_end_ps"] == 3.0
    assert row["n_scored_R"] == 3
    assert row["n_scored_v2_smoothed"] == 3
    assert row["RMSE_R_raw_A"] == pytest.approx(0.0)
    assert row["RMSE_v2_smoothed_Aps"] == pytest.approx(math.sqrt(100.0 / 3.0))
    assert row["v2_smoothed_mean_ratio"] == pytest.approx(
        np.mean([20.0 / 20.0, 30.0 / 30.0, 40.0 / 50.0])
    )
    assert row["ledger_max_resid_eV"] == pytest.approx(2.0)
    assert row["n_sheds"] == 7
    assert row["n_shell_start"] == 21
    assert row["n_shell_end"] == 14

    csv_path = tmp_path / "tier1a.csv"
    write_rows_csv(csv_path, [row])
    with csv_path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert rows[0]["RMSE_v2_smoothed_Aps"] == str(row["RMSE_v2_smoothed_Aps"])


def test_build_tier1a_trajectory_figure_plots_only_selected_v2_traces():
    import matplotlib

    matplotlib.use("Agg")

    from scripts.post_processing.tier1a_rmse_table import (
        Tier1aRunRecord,
        build_trajectory_figure,
    )

    ion = _tiny_ion_checkpoint()
    records = [
        Tier1aRunRecord(
            label="fixed",
            run_tag="tier1a_fixed",
            scenario="fixed",
            t_star_ps=None,
            ion=ion,
            cfg=None,
            row={},
        ),
        Tier1aRunRecord(
            label="anchored t*=0.5",
            run_tag="tier1a_anchored_continuous_t0.5",
            scenario="anchored_discrete",
            t_star_ps=0.5,
            ion=ion,
            cfg=None,
            row={},
        ),
        Tier1aRunRecord(
            label="anchored t*=5.0",
            run_tag="tier1a_anchored_continuous_t5.0",
            scenario="anchored_discrete",
            t_star_ps=5.0,
            ion=ion,
            cfg=None,
            row={},
        ),
    ]

    fig = build_trajectory_figure(
        records,
        _hedft_reference(),
        _smoothed_reference(),
        window=(1.0, 3.0),
        selected_t_star_ps=5.0,
        title="Tier-1a test",
    )

    assert len(fig.axes) == 1
    velocity_labels = {line.get_label() for line in fig.axes[0].lines}
    assert "CEEMDAN+SG |v2|" in velocity_labels
    assert "HeDFT |v2|" in velocity_labels
    assert "fixed MD mean |v2|" in velocity_labels
    assert "anchored t*=5.0 MD mean |v2|" in velocity_labels
    assert "anchored t*=0.5 MD mean |v2|" not in velocity_labels
    assert all("|v1|" not in label and " R" not in label for label in velocity_labels)


def test_export_tier1a_mean_series_writes_one_csv_per_record(tmp_path):
    from scripts.post_processing.tier1a_rmse_table import (
        Tier1aRunRecord,
        export_tier1a_mean_series,
    )

    records = [
        Tier1aRunRecord(
            label="fixed",
            run_tag="tier1a_fixed",
            scenario="fixed",
            t_star_ps=None,
            ion=_tiny_ion_checkpoint(),
            cfg=None,
            row={},
        ),
        Tier1aRunRecord(
            label="anchored t*=5.0",
            run_tag="tier1a_anchored_continuous_t5.0",
            scenario="anchored_discrete",
            t_star_ps=5.0,
            ion=_tiny_ion_checkpoint(),
            cfg=None,
            row={},
        ),
    ]

    paths = export_tier1a_mean_series(tmp_path, records)

    assert [p.name for p in paths] == [
        "tier1a_fixed_mean_trajectory.csv",
        "tier1a_anchored_continuous_t5.0_mean_trajectory.csv",
    ]
    text = paths[1].read_text(encoding="utf-8")
    assert "Tier-1a MD ensemble-mean trajectory" in text
    assert "run_tag=tier1a_anchored_continuous_t5.0" in text
    assert "time_ps,mean_distance_A,mean_speed_I1_Aps,mean_speed_I2_Aps" in text
