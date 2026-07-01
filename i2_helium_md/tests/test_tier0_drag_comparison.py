"""Tier-0 drag-comparison regression gate (TIER0_COMPARISON_spec §7).

Pins the Tier-0 acceptance thresholds so later validation tiers cannot silently
degrade the in-window form match. Anchored on the **18 A** case, which passes
Tier-0 cleanly; the **9 A** case carries a residual diagnosed as windowing +
bubble-mode (TIER0_FINDINGS.md), still under investigation after the
same-smoothed comparison, and is asserted finite-only, not gated.

What is committed
-----------------
Tiny ensemble-mean trajectory CSVs (``data/reference/drag/<case>/tier0/
md_mean_trajectory_N50.csv``) from the from-onset reduced-N (N=50) drag runs.
These four series (time, mean distance, mean |v| of I1/I2) fully reproduce every
Tier-0 scored quantity, since ``compare_*`` reduce to means over molecules. The
N=50 numbers were confirmed to transfer to N=2000 (18 A: dist RMSE 2.512 vs
2.512 A, mean|v| 0.163 vs 0.163 A/ps), so the committed reduced-N reference is a
faithful stand-in (data/runs is gitignored; no production checkpoint is stored).

What this test does / does not catch
-------------------------------------
It pins the committed reference + the windowed-comparison machinery + the named
thresholds (and is the auditable home of the threshold numbers). Re-validating
the drag *physics* after a change is done by regenerating the run
(``scripts/gen_tier0_runs.py`` + the post_processing harness), not in-test.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.physics.constants import U
from i2_helium_md.simulation.checkpoint import IonCheckpoint
from i2_helium_md.postprocess.compare_trajectories import (
    compare_distance,
    compare_velocity_magnitude,
)
from i2_helium_md.postprocess.hedft_loader import load_hedft_trajectory


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REF_ROOT = PROJECT_ROOT / "data" / "reference"
DRAG_ROOT = REF_ROOT / "drag"


# ===========================================================================
# Tier-0 acceptance thresholds (manufactured 2026-06-04, branch drag_implementation)
# ---------------------------------------------------------------------------
# Provenance: set from the 18 A from-onset N=50 run (in-window distance RMSE
# 2.512 A, mean(I1,I2) |v| RMSE 0.163 A/ps), confirmed to transfer to N=2000
# (2.512 A / 0.163 A/ps -- essentially identical). Margins ~20% (distance) and
# ~50% (velocity). See TIER0_FINDINGS.md for the full verdict and the 9 A
# different-regime finding.
# ===========================================================================
TIER0_18A_DISTANCE_RMSE_MAX_A = 3.0
TIER0_18A_MEAN_VELOCITY_RMSE_MAX_APS = 0.25
# Loose symmetry sanity bound (diagnostic, not a physics gate). 18 A split is
# ~0.005 A/ps; this only guards against a gross per-atom asymmetry.
TIER0_18A_VELOCITY_SPLIT_MAX_APS = 0.1


def _load_mean_series(case: str):
    """Load the committed ensemble-mean trajectory CSV for ``case``."""
    path = DRAG_ROOT / case / "tier0" / "md_mean_trajectory_N50.csv"
    data = np.loadtxt(path, delimiter=",", comments="#")
    t, dist, v1, v2 = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    return t, dist, v1, v2


def _read_window(case: str) -> tuple[float, float]:
    params = json.loads(
        (DRAG_ROOT / case / "linear_and_cubic" / "fit_parameters.json").read_text()
    )
    return float(params["t_start"]), float(params["t_end"])


def _ion_from_mean_series(t, dist, v1, v2) -> IonCheckpoint:
    """Rebuild a 1-molecule IonCheckpoint whose molecule-means reproduce the series.

    I1 at +dist/2 and I2 at -dist/2 along x (so the mean I-I distance equals
    ``dist``); per-atom speed carried in vx (so |v| of I1/I2 equals v1/v2). Both
    reductions are taken over the single molecule, so ``compare_*`` recover the
    stored series exactly.
    """
    t = np.asarray(t, float)
    T = t.size
    pos_x = np.empty((2, T))
    pos_x[0] = +dist / 2.0
    pos_x[1] = -dist / 2.0
    zeros = np.zeros((2, T))
    vel_x = np.empty((2, T))
    vel_x[0] = v1
    vel_x[1] = v2
    return IonCheckpoint(
        num_molecules=1,
        time_ps=t,
        positions_x=pos_x, positions_y=zeros.copy(), positions_z=zeros.copy(),
        velocities_x=vel_x, velocities_y=zeros.copy(), velocities_z=zeros.copy(),
        positions_final_x=np.zeros(2), positions_final_y=np.zeros(2),
        positions_final_z=np.zeros(2),
        velocities_final_x=np.zeros(2), velocities_final_y=np.zeros(2),
        velocities_final_z=np.zeros(2),
        mass_kg=np.full(2, 202.953908 * U),
        mass_final_kg=np.full(2, 202.953908 * U),
        mass_history_kg=np.full((2, T), 202.953908 * U),
        droplet_radii_angstrom=np.full(2, 27.936),
        E_kin_eV=zeros.copy(), E_pot_eV=zeros.copy(),
        E_dissip_eV=zeros.copy(), E_mass_transfer_eV=zeros.copy(),
        E_int_eV=zeros.copy(),
        n_shell=zeros.copy(),
        b_ion_outside=np.zeros(1, dtype=bool),
        relative_loss_per_ps=zeros.copy(),
        number_of_collisions=np.zeros((2, T), dtype=int),
        temperature_diagnostic=np.full((T, 3), np.nan, dtype=float),
    )


def _score(case: str):
    t, dist, v1, v2 = _load_mean_series(case)
    ion = _ion_from_mean_series(t, dist, v1, v2)
    hedft = load_hedft_trajectory(REF_ROOT / f"{case}_All_Data.csv")
    window = _read_window(case)
    d = compare_distance(ion, hedft, window=window)
    a = compare_velocity_magnitude(ion, hedft, atom="I1", window=window)
    b = compare_velocity_magnitude(ion, hedft, atom="I2", window=window)
    return d, a, b


class TestTier018A:
    """18 A is the Tier-0 anchor: linear_cubic reproduces TDDFT in-window."""

    def test_distance_rmse_within_threshold(self):
        d, _, _ = _score("18A")
        assert d.rmse <= TIER0_18A_DISTANCE_RMSE_MAX_A, (
            f"18 A in-window distance RMSE {d.rmse:.4f} A exceeds Tier-0 gate "
            f"{TIER0_18A_DISTANCE_RMSE_MAX_A} A"
        )

    def test_mean_velocity_rmse_within_threshold(self):
        _, a, b = _score("18A")
        mean_v_rmse = 0.5 * (a.rmse + b.rmse)
        assert mean_v_rmse <= TIER0_18A_MEAN_VELOCITY_RMSE_MAX_APS, (
            f"18 A in-window mean(I1,I2) |v| RMSE {mean_v_rmse:.4f} A/ps exceeds "
            f"Tier-0 gate {TIER0_18A_MEAN_VELOCITY_RMSE_MAX_APS} A/ps"
        )

    def test_velocity_split_symmetry_sane(self):
        _, a, b = _score("18A")
        split = abs(a.rmse - b.rmse)
        assert split <= TIER0_18A_VELOCITY_SPLIT_MAX_APS, (
            f"18 A I1-I2 velocity split {split:.4f} A/ps is unexpectedly large "
            "(symmetry sanity diagnostic)"
        )


class TestTier09A:
    """9 A carries a diagnosed residual: recorded finite, not gated.

    Recorded numbers (2026-06-05, current clean-window gamma a=24.876/b=2.085,
    from-onset N=50, window [2.67, 6.0]): in-window distance RMSE ~9.17 A,
    mean(I1,I2) |v| RMSE ~1.29 A/ps. The from-onset distance is inflated by the
    uncalibrated pre-t* transient; the t*-seeded clean-form |v2| residual is
    0.88 A/ps raw and 0.40 A/ps against the same-smoothed reference -- bubble-mode
    confirmed as the dominant raw contributor, with a ~0.40 A/ps residual still
    under investigation (TIER0_FINDINGS.md). Not a drag-form failure, and the
    earlier extraction-frame reading is withdrawn.
    """

    def test_produces_finite_windowed_rmse(self):
        d, a, b = _score("9A")
        assert np.isfinite(d.rmse)
        assert np.isfinite(a.rmse) and np.isfinite(b.rmse)
        # Loose sanity envelope only (documents the known divergence; not a gate).
        assert 0.0 < d.rmse < 20.0
        assert 0.0 < 0.5 * (a.rmse + b.rmse) < 3.0
