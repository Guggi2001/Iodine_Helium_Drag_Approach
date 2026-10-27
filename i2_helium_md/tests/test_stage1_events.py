"""Stage-I event log (``simulation/stage1_events.py``, THESIS_DISCUSSION §16).

* Recorder unit behaviour: a same-step shed + pickup gives two rows in the
  frozen shed-then-pickup order with a chained occupancy; crossings carry the
  linearly interpolated time and direction.
* Driver: recording is RNG-free (the checkpoint is bit-identical with and
  without it), the log replays every stored ``n_shell`` column exactly, and
  the crossings match the stored depth signs.
* Guards and the npz round trip.
"""

from __future__ import annotations

from dataclasses import fields, replace

import numpy as np
import pytest

from i2_helium_md.physics.constants import U
from i2_helium_md.simulation.ion import run_ion_propagation
from i2_helium_md.simulation.run_directory import RunDirectory
from i2_helium_md.simulation.stage1_events import (
    KIND_PICKUP,
    KIND_SHED,
    STAGE1_EVENTS_FILENAME,
    Stage1EventRecorder,
    load_stage1_events,
    replay_n_shell,
    save_stage1_events,
)
from tests.test_biphasic_step import _biphasic_driver_cfg


def _surface_neutral():
    """Two ion pairs in a 3 A droplet: pair A straddles the centre (2.6 A apart,
    inside, exits during the Coulomb explosion), pair B starts outside."""
    from i2_helium_md.simulation.checkpoint import NeutralCheckpoint

    N, two_N = 2, 4
    x = np.array([1.3, -1.3, 1.3, -1.3])
    z = np.array([0.0, 0.0, 6.0, 6.0])
    col = lambda a: a.reshape(two_N, 1)  # noqa: E731
    zeros = np.zeros((two_N, 1))
    return NeutralCheckpoint(
        num_molecules=N, time_ps=np.zeros(1),
        positions_x=col(x), positions_y=zeros.copy(), positions_z=col(z),
        velocities_x=zeros.copy(), velocities_y=zeros.copy(), velocities_z=zeros.copy(),
        mass_kg=np.full(two_N, 202.95 * U), droplet_radii=np.full(two_N, 3.0),
        r0=np.zeros(N), E_kin_eV=zeros.copy(), E_pot_eV=zeros.copy(),
        E_initial_eV=np.zeros(N), E_dissip_eV=zeros.copy(),
        L_droplet_eV_ps=zeros.copy(),
    )


def _cfg():
    # 2 ps, decoupled release so the outside ions shed; a high pickup rate so
    # pickups (and same-step pairs) occur inside a 2 ps toy run.
    return replace(_biphasic_driver_cfg(ion_simulation_time=2.0,
                                        self_bound_gate_release="decoupled"),
                   pickup_rate_coefficient=50.0)


@pytest.fixture(scope="module")
def logged_run(tmp_path_factory):
    cfg = _cfg()
    run = RunDirectory(tmp_path_factory.mktemp("logged"))
    ck = run_ion_propagation(cfg, _surface_neutral(), run_dir=run,
                             record_events=True, max_bytes=60_000)
    return cfg, ck, load_stage1_events(run.root / STAGE1_EVENTS_FILENAME)


def test_recorder_same_step_shed_then_pickup_chains_occupancy():
    rec = Stage1EventRecorder(num_atoms=3, dt_ps=0.01)
    pos = np.array([1.0, 2.0, 3.0])
    rec.record_step(5, 0.04, np.array([21.0, 20.0, 19.0]),
                    fired_shed=np.array([True, False, True]),
                    fired_pickup=np.array([True, True, False]),
                    x=pos, y=pos, z=pos, depth=pos - 2.0)
    log = rec.finish()
    np.testing.assert_array_equal(log.event_ion, [0, 2, 0, 1])
    np.testing.assert_array_equal(log.event_kind,
                                  [KIND_SHED, KIND_SHED, KIND_PICKUP, KIND_PICKUP])
    np.testing.assert_array_equal(log.event_n_before, [21, 19, 20, 20])
    np.testing.assert_array_equal(log.event_n_after, [20, 18, 21, 21])
    np.testing.assert_array_equal(log.event_depth_angstrom, [-1.0, 1.0, -1.0, 0.0])


def test_recorder_crossing_time_is_interpolated():
    rec = Stage1EventRecorder(num_atoms=3, dt_ps=0.01)
    rec.record_crossings(7, 0.06, depth_before=np.array([-1.0, 2.0, 0.0]),
                         depth_after=np.array([3.0, -2.0, -1.0]))
    log = rec.finish()
    np.testing.assert_array_equal(log.crossing_ion, [0, 1])   # 0 -> -1 stays inside
    np.testing.assert_array_equal(log.crossing_direction, [1, -1])
    np.testing.assert_allclose(log.crossing_time_ps, [0.0625, 0.065], rtol=1e-12)


def test_recording_is_rng_free_checkpoint_bit_identical(logged_run, tmp_path):
    cfg, ck_logged, _ = logged_run
    ck_plain = run_ion_propagation(cfg, _surface_neutral(),
                                   run_dir=RunDirectory(tmp_path), max_bytes=60_000)
    assert not (tmp_path / STAGE1_EVENTS_FILENAME).exists()
    for f in fields(ck_plain):
        a, b = getattr(ck_plain, f.name), getattr(ck_logged, f.name)
        if isinstance(a, np.ndarray):
            assert np.array_equal(a, b, equal_nan=True), f.name
        else:
            assert a == b, f.name


def test_log_replays_every_stored_n_shell_column(logged_run):
    cfg, ck, log = logged_run
    stride = round((ck.time_ps[1] - ck.time_ps[0]) / cfg.dt_ion)
    assert stride > 1                                  # self-check: strided storage
    assert np.any(log.event_kind == KIND_SHED)         # self-check: events happened
    assert np.any(log.event_kind == KIND_PICKUP)
    steps = np.rint(ck.time_ps / cfg.dt_ion).astype(int)
    replay = replay_n_shell(log, ck.n_shell[:, 0], steps)
    np.testing.assert_array_equal(replay, ck.n_shell.astype(int))


def test_crossings_match_stored_depth_signs(logged_run):
    cfg, ck, log = logged_run
    depth = (np.sqrt(ck.positions_x ** 2 + ck.positions_y ** 2 + ck.positions_z ** 2)
             - ck.droplet_radii_angstrom[:, None])
    assert np.any(log.crossing_direction == 1)         # self-check: pair A exits
    steps = np.rint(ck.time_ps / cfg.dt_ion).astype(int)
    for j, k in enumerate(steps):
        m = log.crossing_step <= k
        net = np.bincount(log.crossing_ion[m], log.crossing_direction[m], log.num_atoms)
        outside = (depth[:, 0] > 0) + net               # 0 inside / 1 outside
        np.testing.assert_array_equal(outside.astype(bool), depth[:, j] > 0)
    t0 = (log.crossing_step - 1) * cfg.dt_ion
    assert np.all((log.crossing_time_ps >= t0 - 1e-12)
                  & (log.crossing_time_ps <= t0 + cfg.dt_ion + 1e-12))


def test_event_rows_carry_pre_step_time(logged_run):
    cfg, _, log = logged_run
    np.testing.assert_allclose(log.event_time_ps, (log.event_step - 1) * cfg.dt_ion,
                               atol=1e-9)


def test_round_trip(logged_run, tmp_path):
    _, _, log = logged_run
    back = load_stage1_events(save_stage1_events(log, tmp_path / "e.npz"))
    for f in fields(log):
        assert np.array_equal(getattr(back, f.name), getattr(log, f.name)), f.name


def test_record_events_needs_run_dir():
    with pytest.raises(ValueError, match="record_events"):
        run_ion_propagation(_cfg(), _surface_neutral(), record_events=True)


def test_record_events_refuses_fixed_mass(tmp_path):
    cfg = replace(_cfg(), mass_scenario="fixed")
    with pytest.raises(ValueError, match="record_events"):
        run_ion_propagation(cfg, _surface_neutral(), run_dir=RunDirectory(tmp_path),
                            record_events=True)
