"""Tests for compare_neutral_to_hedft in compare_trajectories.py."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.compare_trajectories import (
    NeutralComparison,
    compare_neutral_to_hedft,
)
from i2_helium_md.postprocess.hedft_loader import HedftTrajectory
from i2_helium_md.simulation.checkpoint import NeutralCheckpoint


def _make_neutral(*, t_md: np.ndarray, vz1: np.ndarray, vz2: np.ndarray,
                  z01: float = 0.0, z02: float = 0.0,
                  num_molecules: int = 2) -> NeutralCheckpoint:
    n = num_molecules
    t_steps = t_md.size
    # Replicate the same v(t) for every molecule of each atom.
    vz_full = np.zeros((2 * n, t_steps))
    vz_full[:n] = vz1[None, :]
    vz_full[n:] = vz2[None, :]
    # Position from cumulative trapezoid of velocity for self-consistency.
    z1 = z01 + np.concatenate(
        ([0.0], np.cumsum(0.5 * (vz1[1:] + vz1[:-1]) * np.diff(t_md)))
    )
    z2 = z02 + np.concatenate(
        ([0.0], np.cumsum(0.5 * (vz2[1:] + vz2[:-1]) * np.diff(t_md)))
    )
    pos_z = np.zeros((2 * n, t_steps))
    pos_z[:n] = z1[None, :]
    pos_z[n:] = z2[None, :]
    diag_zero = np.zeros((2 * n, t_steps))
    masses_kg = np.full(2 * n, 127.0) * U_KG
    return NeutralCheckpoint(
        num_molecules=n,
        time_ps=t_md,
        positions_x=np.zeros((2 * n, t_steps)),
        positions_y=np.zeros((2 * n, t_steps)),
        positions_z=pos_z,
        velocities_x=np.zeros((2 * n, t_steps)),
        velocities_y=np.zeros((2 * n, t_steps)),
        velocities_z=vz_full,
        mass_kg=masses_kg,
        droplet_radii=np.full(2 * n, 30.0),
        r0=np.zeros(n),
        E_kin_eV=diag_zero,
        E_pot_eV=diag_zero,
        E_initial_eV=np.zeros(n),
        E_dissip_eV=diag_zero,
        L_droplet_eV_ps=diag_zero,
    )


def _make_hedft(t_ref: np.ndarray, v1z: np.ndarray, v2z: np.ndarray,
                tmp_path: Path) -> HedftTrajectory:
    return HedftTrajectory(
        time_ps=t_ref,
        v1_magnitude_Aps=np.abs(v1z),
        v2_magnitude_Aps=np.abs(v2z),
        v1_z_Aps=v1z,
        v2_z_Aps=v2z,
        v1_x_Aps=np.zeros_like(v1z),
        v2_x_Aps=np.zeros_like(v2z),
        distance_A=np.linspace(2.0, 5.0, t_ref.size),
        droplet_radius_A=9.0,
        source_path=tmp_path / "synth.csv",
    )


class TestCompareNeutralToHedft:
    def test_constant_velocity_zero_rmse(self, tmp_path):
        t = np.linspace(0.0, 4.0, 21)
        v1 = np.full(t.size, 0.5)
        v2 = np.full(t.size, -0.5)
        neutral = _make_neutral(t_md=t, vz1=v1, vz2=v2)
        hedft = _make_hedft(t, v1, v2, tmp_path)

        cmp1 = compare_neutral_to_hedft(neutral, hedft, atom="I1")
        assert isinstance(cmp1, NeutralComparison)
        assert cmp1.atom == "I1"
        assert cmp1.rmse_velocity_Aps == pytest.approx(0.0, abs=1e-12)
        assert cmp1.rmse_position_A == pytest.approx(0.0, abs=1e-12)

    def test_offset_velocity_yields_known_rmse(self, tmp_path):
        t = np.linspace(0.0, 4.0, 41)
        v_md = np.full(t.size, 1.0)
        v_ref = np.full(t.size, 1.5)
        neutral = _make_neutral(t_md=t, vz1=v_md, vz2=v_md)
        hedft = _make_hedft(t, v_ref, v_ref, tmp_path)

        cmp1 = compare_neutral_to_hedft(neutral, hedft, atom="I1")
        assert cmp1.rmse_velocity_Aps == pytest.approx(0.5, rel=1e-6)
        # Position drift over t_max-t_min = 4 ps with a 0.5 A/ps gap.
        # Last sample of each is z_md_start + integral of v -> drift accumulates.
        assert cmp1.rmse_position_A > 0.0

    def test_no_overlap_raises(self, tmp_path):
        t_md = np.linspace(0.0, 1.0, 11)
        t_ref = np.linspace(5.0, 6.0, 11)
        v_md = np.zeros(t_md.size)
        v_ref = np.zeros(t_ref.size)
        neutral = _make_neutral(t_md=t_md, vz1=v_md, vz2=v_md)
        hedft = _make_hedft(t_ref, v_ref, v_ref, tmp_path)
        with pytest.raises(ValueError, match="overlap"):
            compare_neutral_to_hedft(neutral, hedft, atom="I1")

    def test_invalid_atom(self, tmp_path):
        t = np.linspace(0.0, 1.0, 5)
        v = np.zeros(t.size)
        neutral = _make_neutral(t_md=t, vz1=v, vz2=v)
        hedft = _make_hedft(t, v, v, tmp_path)
        with pytest.raises(ValueError, match="atom"):
            compare_neutral_to_hedft(neutral, hedft, atom="I3")  # type: ignore[arg-type]
