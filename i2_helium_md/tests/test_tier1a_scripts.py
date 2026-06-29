"""Coverage for the Tier-1a RMSE-table orchestration layer.

The production generator writes under ``data/runs`` and can be expensive. These
tests stay on small temporary run directories and exercise the reusable helpers
that the scripts call.
"""

from __future__ import annotations

import csv
from dataclasses import replace
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


def _tiny_ion_checkpoint_with_shell_end(
    n_final: int,
    *,
    mass_scenario: str = "anchored_discrete",
) -> IonCheckpoint:
    ion = _tiny_ion_checkpoint()
    n_shell = ion.n_shell.copy()
    n_shell[:, -1] = n_final
    return replace(ion, n_shell=n_shell, mass_scenario=mass_scenario)


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
    assert cfg.coulomb_available_eV == 0.80
    assert cfg.mass_initial_amu == pytest.approx(complex_mass_amu(21))
    assert cfg.allow_inconsistent_mass_pairing is True


@pytest.mark.parametrize("n_final", [True, "2", 2.9, -1, 21])
def test_build_onset_strip_cfg_rejects_invalid_n_final(n_final):
    from scripts.tier1a_common import build_onset_strip_cfg

    with pytest.raises(ValueError, match="n_final"):
        build_onset_strip_cfg(
            "9A",
            "shared_pure_cubic",
            n_final=n_final,
            t_strip_ps=0.5,
            num_molecules=2,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
        )


@pytest.mark.parametrize("t_strip_ps", [-1, math.nan, math.inf, 0.51])
def test_build_onset_strip_cfg_rejects_invalid_t_strip_ps(t_strip_ps):
    from scripts.tier1a_common import build_onset_strip_cfg

    with pytest.raises(ValueError, match="t_strip_ps"):
        build_onset_strip_cfg(
            "9A",
            "shared_pure_cubic",
            n_final=2,
            t_strip_ps=t_strip_ps,
            num_molecules=2,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
        )


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


@pytest.mark.parametrize("n_final", [True, "2", 2.9, -1, 21])
def test_tier1a_stress_run_tag_rejects_invalid_n_final(n_final):
    from scripts.tier1a_common import tier1a_stress_run_tag

    with pytest.raises(ValueError, match="n_final"):
        tier1a_stress_run_tag(n_final=n_final, t_strip_ps=0.5)


@pytest.mark.parametrize("t_strip_ps", [-1, math.nan, math.inf, 0.51])
def test_tier1a_stress_run_tag_rejects_invalid_t_strip_ps(t_strip_ps):
    from scripts.tier1a_common import tier1a_stress_run_tag

    with pytest.raises(ValueError, match="t_strip_ps"):
        tier1a_stress_run_tag(n_final=2, t_strip_ps=t_strip_ps)


def test_tier1a_stress_generator_refuses_existing_outputs(tmp_path, monkeypatch):
    from scripts import gen_tier1a_stress_runs as script

    run_dir = tmp_path / "existing_run"
    run_dir.mkdir()
    (run_dir / "cfg.json").write_text("{}", encoding="utf-8")
    monkeypatch.setattr(script, "OVERWRITE_EXISTING_RUN", False)

    with pytest.raises(FileExistsError, match="OVERWRITE_EXISTING_RUN"):
        script._run_one("existing", object(), run_dir)


def test_tier1a_stress_generator_main_schedules_expected_runs(tmp_path, monkeypatch):
    from scripts import gen_tier1a_stress_runs as script
    from scripts.tier1a_common import (
        TIER1A_STRESS_N_FINAL_VALUES,
        TIER1A_STRESS_T_STRIP_PS,
        tier1a_run_dir_name,
        tier1a_stress_run_dir_name,
    )

    scheduled = []

    def capture_run(label, cfg, run_dir):
        scheduled.append((label, cfg, run_dir))

    monkeypatch.setattr(script, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(script, "_run_one", capture_run)

    assert script.main() == 0

    expected_dirs = [
        tmp_path
        / "data"
        / "runs"
        / tier1a_run_dir_name(script.CASE, script.VARIANT, script.N, "fixed", None)
    ]
    expected_dirs.extend(
        tmp_path
        / "data"
        / "runs"
        / tier1a_stress_run_dir_name(
            script.CASE,
            script.VARIANT,
            script.N,
            n_final=n_final,
            t_strip_ps=TIER1A_STRESS_T_STRIP_PS,
        )
        for n_final in TIER1A_STRESS_N_FINAL_VALUES
    )

    assert [run_dir for _, _, run_dir in scheduled] == expected_dirs
    assert scheduled[0][1].mass_scenario == "fixed"
    assert [cfg.anchor_n_final for _, cfg, _ in scheduled[1:]] == list(
        TIER1A_STRESS_N_FINAL_VALUES
    )
    assert all(cfg.anchor_mode == "onset_strip" for _, cfg, _ in scheduled[1:])


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


def test_score_tier1a_stress_run_emits_stress_columns(tmp_path):
    from scripts.post_processing.tier1a_stress_table import (
        TIER1A_STRESS_TABLE_COLUMNS,
        score_tier1a_stress_run,
    )
    from scripts.tier1a_common import build_onset_strip_cfg

    run = RunDirectory(tmp_path / "stress")
    run.save_ion(_tiny_ion_checkpoint())
    run.save_cfg(
        build_onset_strip_cfg(
            "9A",
            "shared_pure_cubic",
            n_final=14,
            t_strip_ps=0.5,
            num_molecules=1,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
        )
    )

    row = score_tier1a_stress_run(
        run,
        case="9A",
        variant="shared_pure_cubic",
        n=1,
        run_tag="tier1a_stress_onset_strip_n14_t0.5",
        n_final=14,
        t_strip_ps=0.5,
        hedft=_hedft_reference(),
        smoothed=_smoothed_reference(),
        window_start_ps=1.0,
    )

    assert list(row) == TIER1A_STRESS_TABLE_COLUMNS
    assert row["stress_family"] == "onset_strip"
    assert row["n_final_requested"] == 14
    assert row["t_strip_ps"] == 0.5
    assert row["n_removed"] == 7


def test_score_tier1a_stress_run_requires_cfg_json(tmp_path):
    from scripts.post_processing.tier1a_stress_table import score_tier1a_stress_run

    run = RunDirectory(tmp_path / "stress")
    run.save_ion(_tiny_ion_checkpoint())

    with pytest.raises(ValueError, match="missing cfg.json"):
        score_tier1a_stress_run(
            run,
            case="9A",
            variant="shared_pure_cubic",
            n=1,
            run_tag="tier1a_stress_onset_strip_n14_t0.5",
            n_final=14,
            t_strip_ps=0.5,
            hedft=_hedft_reference(),
            smoothed=_smoothed_reference(),
            window_start_ps=1.0,
        )


def _write_tier1a_stress_collect_fixture(
    tmp_path,
    *,
    stress_ion: IonCheckpoint | None = None,
    stress_cfg=None,
    fixed_ion: IonCheckpoint | None = None,
    fixed_cfg=None,
    write_fixed_cfg: bool = True,
    write_stress_cfg: bool = True,
):
    from scripts.tier1a_common import (
        build_onset_strip_cfg,
        tier1a_run_dir_name,
        tier1a_stress_run_dir_name,
    )
    from scripts.tier0_common import build_drag_cfg

    case = "9A"
    variant = "shared_pure_cubic"
    n = 1
    run_root = tmp_path / "data" / "runs"

    if fixed_ion is None:
        fixed_ion = replace(_tiny_ion_checkpoint(), mass_scenario="fixed")
    if fixed_cfg is None:
        fixed_cfg = build_drag_cfg(
            case,
            variant,
            num_molecules=1,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
        )
    fixed = RunDirectory(
        run_root / tier1a_run_dir_name(case, variant, n, "fixed", None)
    )
    fixed.save_ion(fixed_ion)
    if write_fixed_cfg:
        fixed.save_cfg(fixed_cfg)

    if stress_ion is None:
        stress_ion = _tiny_ion_checkpoint_with_shell_end(0)
    if stress_cfg is None:
        stress_cfg = build_onset_strip_cfg(
            case,
            variant,
            n_final=0,
            t_strip_ps=0.5,
            num_molecules=1,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
        )
    stress = RunDirectory(
        run_root
        / tier1a_stress_run_dir_name(case, variant, n, n_final=0, t_strip_ps=0.5)
    )
    stress.save_ion(stress_ion)
    if write_stress_cfg:
        stress.save_cfg(stress_cfg)

    return case, variant, n


def _collect_stress_fixture(script, tmp_path, case, variant, n, **kwargs):
    return script.collect_tier1a_stress_records(
        project_root=tmp_path,
        case=case,
        variant=variant,
        n=n,
        n_final_values=kwargs.pop("n_final_values", (0,)),
        t_strip_ps=kwargs.pop("t_strip_ps", 0.5),
        window_start_ps=1.0,
        **kwargs,
    )


@pytest.mark.parametrize(
    "missing_cfg_kwargs, match",
    [
        ({"write_fixed_cfg": False}, "missing cfg.json.*tier1a_fixed"),
        ({"write_stress_cfg": False}, "missing cfg.json.*tier1a_stress_onset_strip_n0_t0.5"),
    ],
)
def test_collect_tier1a_stress_records_requires_cfg_json(
    tmp_path,
    monkeypatch,
    missing_cfg_kwargs,
    match,
):
    from scripts.post_processing import tier1a_stress_table as script

    case, variant, n = _write_tier1a_stress_collect_fixture(
        tmp_path,
        **missing_cfg_kwargs,
    )
    monkeypatch.setattr(
        script,
        "_load_references",
        lambda project_root, case: (_hedft_reference(), _smoothed_reference()),
    )

    with pytest.raises(ValueError, match=match):
        _collect_stress_fixture(script, tmp_path, case, variant, n)


@pytest.mark.parametrize(
    "ion_kwargs, match",
    [
        (
            {"stress_ion": _tiny_ion_checkpoint_with_shell_end(0, mass_scenario="fixed")},
            "tier1a_stress_onset_strip_n0_t0.5.*ion.mass_scenario",
        ),
        (
            {"fixed_ion": replace(_tiny_ion_checkpoint(), mass_scenario="anchored_discrete")},
            "tier1a_fixed.*ion.mass_scenario",
        ),
    ],
)
def test_collect_tier1a_stress_records_rejects_ion_mass_scenario_mismatch(
    tmp_path,
    monkeypatch,
    ion_kwargs,
    match,
):
    from scripts.post_processing import tier1a_stress_table as script

    case, variant, n = _write_tier1a_stress_collect_fixture(
        tmp_path,
        **ion_kwargs,
    )
    monkeypatch.setattr(
        script,
        "_load_references",
        lambda project_root, case: (_hedft_reference(), _smoothed_reference()),
    )

    with pytest.raises(ValueError, match=match):
        _collect_stress_fixture(script, tmp_path, case, variant, n)


def test_collect_tier1a_stress_records_rejects_shell_endpoint_mismatch(
    tmp_path,
    monkeypatch,
):
    from scripts.post_processing import tier1a_stress_table as script

    ion = _tiny_ion_checkpoint_with_shell_end(0)
    n_shell = ion.n_shell.copy()
    n_shell[:, -1] = 2
    case, variant, n = _write_tier1a_stress_collect_fixture(
        tmp_path,
        stress_ion=replace(ion, n_shell=n_shell),
    )
    monkeypatch.setattr(
        script,
        "_load_references",
        lambda project_root, case: (_hedft_reference(), _smoothed_reference()),
    )

    with pytest.raises(
        ValueError,
        match="tier1a_stress_onset_strip_n0_t0.5.*n_shell_end",
    ):
        _collect_stress_fixture(script, tmp_path, case, variant, n)


@pytest.mark.parametrize(
    "cfg_patch, match",
    [
        ({"anchor_mode": "time"}, "anchor_mode"),
        ({"t_star_ps": 5.0}, "t_star_ps"),
    ],
)
def test_collect_tier1a_stress_records_rejects_anchor_cfg_mismatch(
    tmp_path,
    monkeypatch,
    cfg_patch,
    match,
):
    from scripts.post_processing import tier1a_stress_table as script
    from scripts.tier1a_common import build_onset_strip_cfg

    cfg = build_onset_strip_cfg(
        "9A",
        "shared_pure_cubic",
        n_final=0,
        t_strip_ps=0.5,
        num_molecules=1,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=123,
    )
    case, variant, n = _write_tier1a_stress_collect_fixture(
        tmp_path,
        stress_cfg=replace(cfg, **cfg_patch),
    )
    monkeypatch.setattr(
        script,
        "_load_references",
        lambda project_root, case: (_hedft_reference(), _smoothed_reference()),
    )

    with pytest.raises(ValueError, match=match):
        _collect_stress_fixture(script, tmp_path, case, variant, n)


@pytest.mark.parametrize(
    "collect_kwargs, match",
    [
        ({"n_final_values": (True,)}, "n_final"),
        ({"n_final_values": ("0",)}, "n_final"),
        ({"t_strip_ps": math.nan}, "t_strip_ps"),
    ],
)
def test_collect_tier1a_stress_records_preserves_strict_input_validation(
    tmp_path,
    monkeypatch,
    collect_kwargs,
    match,
):
    from scripts.post_processing import tier1a_stress_table as script

    case, variant, n = _write_tier1a_stress_collect_fixture(tmp_path)
    monkeypatch.setattr(
        script,
        "_load_references",
        lambda project_root, case: (_hedft_reference(), _smoothed_reference()),
    )

    with pytest.raises(ValueError, match=match):
        _collect_stress_fixture(script, tmp_path, case, variant, n, **collect_kwargs)


def test_score_tier1a_stress_run_rejects_nonuniform_shell_endpoint(tmp_path):
    from scripts.post_processing.tier1a_stress_table import score_tier1a_stress_run
    from scripts.tier1a_common import build_onset_strip_cfg

    ion = _tiny_ion_checkpoint_with_shell_end(0)
    n_shell = ion.n_shell.copy()
    n_shell[1, -1] = 1
    run = RunDirectory(tmp_path / "stress")
    run.save_ion(replace(ion, n_shell=n_shell))
    run.save_cfg(
        build_onset_strip_cfg(
            "9A",
            "shared_pure_cubic",
            n_final=0,
            t_strip_ps=0.5,
            num_molecules=1,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
        )
    )

    with pytest.raises(ValueError, match="n_shell end values are not uniform"):
        score_tier1a_stress_run(
            run,
            case="9A",
            variant="shared_pure_cubic",
            n=1,
            run_tag="tier1a_stress_onset_strip_n0_t0.5",
            n_final=0,
            t_strip_ps=0.5,
            hedft=_hedft_reference(),
            smoothed=_smoothed_reference(),
            window_start_ps=1.0,
        )


def test_score_tier1a_stress_run_rejects_fractional_shell_endpoint(tmp_path):
    from scripts.post_processing.tier1a_stress_table import score_tier1a_stress_run
    from scripts.tier1a_common import build_onset_strip_cfg

    ion = _tiny_ion_checkpoint_with_shell_end(0)
    n_shell = ion.n_shell.copy()
    n_shell[:, -1] = 0.25
    run = RunDirectory(tmp_path / "stress")
    run.save_ion(replace(ion, n_shell=n_shell))
    run.save_cfg(
        build_onset_strip_cfg(
            "9A",
            "shared_pure_cubic",
            n_final=0,
            t_strip_ps=0.5,
            num_molecules=1,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
        )
    )

    with pytest.raises(ValueError, match="n_shell end values must be integer-valued"):
        score_tier1a_stress_run(
            run,
            case="9A",
            variant="shared_pure_cubic",
            n=1,
            run_tag="tier1a_stress_onset_strip_n0_t0.5",
            n_final=0,
            t_strip_ps=0.5,
            hedft=_hedft_reference(),
            smoothed=_smoothed_reference(),
            window_start_ps=1.0,
        )


def test_collect_tier1a_stress_records_rejects_cfg_metadata_mismatch(tmp_path, monkeypatch):
    from scripts.post_processing import tier1a_stress_table as script
    from scripts.tier1a_common import (
        build_onset_strip_cfg,
        tier1a_run_dir_name,
        tier1a_stress_run_dir_name,
    )
    from scripts.tier0_common import build_drag_cfg

    case = "9A"
    variant = "shared_pure_cubic"
    n = 1
    run_root = tmp_path / "data" / "runs"
    fixed_ion = _tiny_ion_checkpoint_with_shell_end(14, mass_scenario="fixed")
    stress_ion = _tiny_ion_checkpoint_with_shell_end(0)

    fixed = RunDirectory(
        run_root / tier1a_run_dir_name(case, variant, n, "fixed", None)
    )
    fixed.save_ion(fixed_ion)
    fixed.save_cfg(
        build_drag_cfg(
            case,
            variant,
            num_molecules=1,
            ion_time_ps=0.02,
            dt_ion_ps=0.01,
            seed=123,
        )
    )

    mismatched = RunDirectory(
        run_root
        / tier1a_stress_run_dir_name(case, variant, n, n_final=0, t_strip_ps=0.5)
    )
    mismatched.save_ion(stress_ion)
    cfg = build_onset_strip_cfg(
        case,
        variant,
        n_final=2,
        t_strip_ps=0.5,
        num_molecules=1,
        ion_time_ps=0.02,
        dt_ion_ps=0.01,
        seed=123,
    )
    mismatched.save_cfg(cfg)

    monkeypatch.setattr(
        script,
        "_load_references",
        lambda project_root, case: (_hedft_reference(), _smoothed_reference()),
    )

    with pytest.raises(
        ValueError,
        match="tier1a_stress_onset_strip_n0_t0.5.*anchor_n_final",
    ):
        script.collect_tier1a_stress_records(
            project_root=tmp_path,
            case=case,
            variant=variant,
            n=n,
            n_final_values=(0,),
            t_strip_ps=0.5,
            window_start_ps=1.0,
        )


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


def test_build_tier1a_stress_figure_plots_fixed_plus_selected_case():
    import matplotlib

    matplotlib.use("Agg")

    from scripts.post_processing.tier1a_stress_table import (
        Tier1aStressRunRecord,
        build_stress_trajectory_figure,
    )

    ion = _tiny_ion_checkpoint()
    records = [
        Tier1aStressRunRecord(
            label="fixed",
            run_tag="tier1a_fixed",
            n_final=None,
            t_strip_ps=None,
            ion=ion,
            cfg=None,
            row={},
        ),
        Tier1aStressRunRecord(
            label="onset strip n=0",
            run_tag="tier1a_stress_onset_strip_n0_t0.5",
            n_final=0,
            t_strip_ps=0.5,
            ion=ion,
            cfg=None,
            row={},
        ),
        Tier1aStressRunRecord(
            label="onset strip n=2",
            run_tag="tier1a_stress_onset_strip_n2_t0.5",
            n_final=2,
            t_strip_ps=0.5,
            ion=ion,
            cfg=None,
            row={},
        ),
    ]

    fig = build_stress_trajectory_figure(
        records,
        _hedft_reference(),
        _smoothed_reference(),
        window=(1.0, 3.0),
        selected_n_final=0,
        title="Tier-1a stress test",
    )

    labels = {line.get_label() for line in fig.axes[0].lines}
    assert "fixed MD mean |v2|" in labels
    assert "onset strip n=0 MD mean |v2|" in labels
    assert "onset strip n=2 MD mean |v2|" not in labels


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
