"""Tests for i2_helium_md/postprocess/ihe_ked.py.

The loader tests run against the real frozen reference under
``data/reference/ihe_ked/`` (small CSVs, committed) plus synthetic
corrupted files in tmp_path for the failure paths.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.ihe_ked import (
    IHeKedReference,
    load_ihe_ked_reference,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IHE_KED_DIR = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
REFERENCE_CSV = IHE_KED_DIR / "IHe_KED_reference.csv"


# ===========================================================================
# Helpers
# ===========================================================================
def _make_ion(
    *,
    num_molecules: int,
    final_speeds_per_atom: np.ndarray,
    masses_amu_per_atom: np.ndarray,
    b_outside: np.ndarray | None = None,
    num_steps: int = 4,
) -> IonCheckpoint:
    """Build a tiny IonCheckpoint with prescribed final |v| and final masses.

    Copied verbatim from ``tests/test_velocity_distribution.py`` (private
    test helpers are not shared across test modules). Atom index ``j`` and
    ``j + num_molecules`` belong to molecule ``j`` (concatenation layout,
    NOT adjacent pairing) -- see ``select_final_mass_gate``'s
    ``np.concatenate([b_ion_outside, b_ion_outside])``.

    The full velocity time series is filled with zeros (not used by the
    histogram code, which only reads ``velocities_final_*``).
    """
    n = num_molecules
    if final_speeds_per_atom.shape != (2 * n,):
        raise AssertionError("final_speeds_per_atom must have shape (2N,)")
    if masses_amu_per_atom.shape != (2 * n,):
        raise AssertionError("masses_amu_per_atom must have shape (2N,)")

    if b_outside is None:
        b_outside = np.ones(n, dtype=bool)

    # Distribute the speed in a simple way so each atom's speed magnitude
    # equals the prescribed value: put it all on velocities_final_x.
    vfx = final_speeds_per_atom.astype(float)
    vfy = np.zeros(2 * n)
    vfz = np.zeros(2 * n)

    masses_kg = masses_amu_per_atom.astype(float) * U_KG

    diag_zero = np.zeros((2 * n, num_steps))
    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, num_steps),
        positions_x=np.zeros((2 * n, num_steps)),
        positions_y=np.zeros((2 * n, num_steps)),
        positions_z=np.zeros((2 * n, num_steps)),
        velocities_x=np.zeros((2 * n, num_steps)),
        velocities_y=np.zeros((2 * n, num_steps)),
        velocities_z=np.zeros((2 * n, num_steps)),
        positions_final_x=np.zeros(2 * n),
        positions_final_y=np.zeros(2 * n),
        positions_final_z=np.zeros(2 * n),
        velocities_final_x=vfx,
        velocities_final_y=vfy,
        velocities_final_z=vfz,
        mass_kg=masses_kg.copy(),
        mass_final_kg=masses_kg,
        mass_history_kg=np.broadcast_to(
            masses_kg[:, None], (2 * n, num_steps)
        ).copy(),
        droplet_radii_angstrom=np.full(2 * n, 30.0),
        E_kin_eV=diag_zero,
        E_pot_eV=diag_zero,
        E_dissip_eV=diag_zero,
        E_mass_transfer_eV=diag_zero,
        E_int_eV=diag_zero,
        n_shell=diag_zero,
        b_ion_outside=np.asarray(b_outside, dtype=bool),
        relative_loss_per_ps=diag_zero,
        number_of_collisions=np.zeros((2 * n, num_steps), dtype=int),
        temperature_diagnostic=np.full((num_steps, 3), np.nan, dtype=float),
    )


class TestLoadIHeKedReference:
    def test_real_file_contract(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        assert isinstance(ref, IHeKedReference)
        np.testing.assert_array_equal(ref.n, np.arange(18))
        # Frozen reference value (README: I+ scale anchor 3.706 eV).
        assert ref.mean_KE_eV[0] == pytest.approx(3.70569, abs=1e-5)
        # Mass model column matches m(n) = 126.90 + 4.0026 n.
        np.testing.assert_allclose(
            ref.mass_center_u, 126.90 + 4.0026 * ref.n, atol=1e-4,
        )
        # noiseLimited is 0 for all 18 fragments (README).
        assert not ref.noise_limited.any()
        assert ref.source_path == REFERENCE_CSV.resolve()

    def test_gold_mask_is_calib_limited_set(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        expected_gold = {0, 3, 4, 5, 6, 7, 8, 9, 10, 12}
        assert set(ref.n[ref.gold_mask].tolist()) == expected_gold

    def test_point_err_is_stat_sys_quadrature(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        np.testing.assert_allclose(
            ref.point_err_eV,
            np.sqrt(
                ref.stat_err_mean_KE_eV ** 2 + ref.sys_err_mean_KE_eV ** 2
            ),
            rtol=1e-12,  # tight: pure arithmetic identity
        )

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_ihe_ked_reference(tmp_path / "nope.csv")

    def test_missing_column_raises(self, tmp_path):
        p = tmp_path / "bad.csv"
        p.write_text("n,label\n0,I+He_0\n", encoding="ascii")
        with pytest.raises(ValueError):
            load_ihe_ked_reference(p)

    def test_non_contiguous_n_raises(self, tmp_path):
        header = (
            "n,label,massCenter_u_per_e,N_counts,N_eff,meanKE_eV,"
            "modeKE_eV,medianKE_eV,sigmaKE_eV,statErr_meanKE_eV,"
            "sysErr_meanKE_eV,calibSyst_frac,conditionSyst_frac,"
            "bgOffShift_eV,dominantError,noiseLimited\n"
        )
        row = "5,I+He_5,146.9,100,90,1.0,0.9,0.95,0.5,0.01,0.02,0.04,0.06,0.0,calib,0\n"
        p = tmp_path / "bad.csv"
        p.write_text(header + row, encoding="ascii")
        with pytest.raises(ValueError, match="contiguous"):
            load_ihe_ked_reference(p)

    def test_unknown_dominant_error_raises(self, tmp_path):
        header = (
            "n,label,massCenter_u_per_e,N_counts,N_eff,meanKE_eV,"
            "modeKE_eV,medianKE_eV,sigmaKE_eV,statErr_meanKE_eV,"
            "sysErr_meanKE_eV,calibSyst_frac,conditionSyst_frac,"
            "bgOffShift_eV,dominantError,noiseLimited\n"
        )
        row = "0,I+He_0,126.9,100,90,1.0,0.9,0.95,0.5,0.01,0.02,0.04,0.06,0.0,vibes,0\n"
        p = tmp_path / "bad.csv"
        p.write_text(header + row, encoding="ascii")
        with pytest.raises(ValueError, match="dominantError"):
            load_ihe_ked_reference(p)


from i2_helium_md.postprocess.ihe_ked import (  # noqa: E402 (grouped here)
    IHeKedCurve,
    load_ihe_ked_curve,
)


class TestLoadIHeKedCurve:
    def test_real_n0_contract(self):
        curve = load_ihe_ked_curve(IHE_KED_DIR, 0)
        assert isinstance(curve, IHeKedCurve)
        assert curve.n == 0
        # Axes are finite and strictly ascending.
        assert np.all(np.isfinite(curve.E_eV))
        assert np.all(np.diff(curve.E_eV) > 0)
        assert np.all(np.diff(curve.v_mps) > 0)
        # n=0 Abel-center spike cut: 3-D columns are NaN below 0.4 eV,
        # while the 2-D columns cover the full detector range.
        low = curve.E_eV < 0.4
        assert np.all(np.isnan(curve.signal_3d_Pv[low]))
        assert np.all(np.isfinite(curve.signal_2d_Pv))
        # Above the cut the 3-D reconstruction is real data.
        assert np.isfinite(curve.signal_3d_Pv[~low]).any()

    def test_all_five_fragments_load(self):
        for n in range(5):
            curve = load_ihe_ked_curve(IHE_KED_DIR, n)
            assert curve.n == n

    def test_out_of_range_n_raises(self):
        with pytest.raises(ValueError, match="n must be"):
            load_ihe_ked_curve(IHE_KED_DIR, 5)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_ihe_ked_curve(tmp_path, 0)

    def test_missing_column_raises(self, tmp_path):
        p = tmp_path / "IHe_KED_curves_n0.csv"
        p.write_text("E_eV,v_mps\n0.1,100\n0.2,140\n", encoding="ascii")
        with pytest.raises(ValueError):
            load_ihe_ked_curve(tmp_path, 0)


from i2_helium_md.physics.constants import EV  # noqa: E402 (grouped here)
from i2_helium_md.physics.shell_schedule import complex_mass_amu  # noqa: E402
from i2_helium_md.postprocess.ihe_ked import (  # noqa: E402
    FragmentMeanKE,
    fragment_gate_counts,
    fragment_mean_kinetic_energy,
    speed_mps_of_energy_eV,
)


def _energy_eV_of_speed_Aps(speed_Aps: float, mass_amu: float) -> float:
    """Independent re-derivation: E = 1/2 m v^2, v in m/s (1 A/ps = 100 m/s)."""
    from i2_helium_md.physics.constants import U as _U
    v_mps = speed_Aps * 100.0
    return 0.5 * mass_amu * _U * v_mps * v_mps / EV


class TestFragmentMeanKineticEnergy:
    def test_analytic_mean_and_stat_err(self):
        # Two n=1 atoms (131 amu, gated by m(1)=130.9026) at 2 and 4 A/ps;
        # one n=0 atom that must not contaminate the gate.
        ion = _make_ion(
            num_molecules=2,
            final_speeds_per_atom=np.array([2.0, 4.0, 3.0, 0.0]),
            masses_amu_per_atom=np.array([131.0, 131.0, 127.0, 127.0]),
        )
        result = fragment_mean_kinetic_energy(ion, 1)
        assert isinstance(result, FragmentMeanKE)
        m1 = complex_mass_amu(1)
        e1 = _energy_eV_of_speed_Aps(2.0, m1)
        e2 = _energy_eV_of_speed_Aps(4.0, m1)
        expected_mean = 0.5 * (e1 + e2)
        # Analytical port: tight tolerance.
        assert result.mean_KE_eV == pytest.approx(expected_mean, rel=1e-12)
        expected_err = np.std([e1, e2], ddof=1) / np.sqrt(2.0)
        assert result.stat_err_mean_KE_eV == pytest.approx(
            expected_err, rel=1e-12,
        )
        assert result.num_atoms_used == 2
        assert result.mass_amu == pytest.approx(m1)

    def test_v_of_mean_E_round_trip(self):
        ion = _make_ion(
            num_molecules=1,
            final_speeds_per_atom=np.array([3.0, 0.0]),
            masses_amu_per_atom=np.array([131.0, 127.0]),
        )
        result = fragment_mean_kinetic_energy(ion, 1)
        # A single atom at 3 A/ps: v(<E>) must be exactly 300 m/s.
        assert result.v_of_mean_E_mps == pytest.approx(300.0, rel=1e-12)
        # And the standalone converter agrees.
        assert speed_mps_of_energy_eV(
            result.mean_KE_eV, result.mass_amu
        ) == pytest.approx(300.0, rel=1e-12)

    def test_single_atom_has_zero_stat_err(self):
        ion = _make_ion(
            num_molecules=1,
            final_speeds_per_atom=np.array([3.0, 0.0]),
            masses_amu_per_atom=np.array([131.0, 127.0]),
        )
        result = fragment_mean_kinetic_energy(ion, 1)
        assert result.stat_err_mean_KE_eV == 0.0

    def test_empty_gate_raises(self):
        ion = _make_ion(
            num_molecules=1,
            final_speeds_per_atom=np.array([3.0, 0.0]),
            masses_amu_per_atom=np.array([127.0, 127.0]),
        )
        with pytest.raises(ValueError, match="No atoms"):
            fragment_mean_kinetic_energy(ion, 1)


class TestFragmentGateCounts:
    def test_counts_per_gate(self):
        ion = _make_ion(
            num_molecules=2,
            final_speeds_per_atom=np.array([1.0, 2.0, 3.0, 4.0]),
            masses_amu_per_atom=np.array([127.0, 131.0, 131.0, 135.0]),
        )
        counts = fragment_gate_counts(ion, np.arange(4))
        np.testing.assert_array_equal(counts, [1, 2, 1, 0])

    def test_outside_filter_applies(self):
        ion = _make_ion(
            num_molecules=2,
            final_speeds_per_atom=np.array([1.0, 2.0, 3.0, 4.0]),
            masses_amu_per_atom=np.array([131.0, 131.0, 131.0, 131.0]),
            b_outside=np.array([True, False]),
        )
        counts = fragment_gate_counts(ion, np.arange(3))
        np.testing.assert_array_equal(counts, [0, 2, 0])


import importlib.util
import sys

import matplotlib


def _load_run_summary_module():
    matplotlib.use("Agg", force=True)
    script = (
        PROJECT_ROOT / "scripts" / "post_processing" / "plot_run_summary.py"
    )
    spec = importlib.util.spec_from_file_location(
        "plot_run_summary_under_test", script
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _smoke_ion() -> IonCheckpoint:
    """A tiny ensemble populating the n = 0..2 gates with plausible speeds."""
    speeds = np.array([24.0, 20.0, 10.0, 8.0, 6.0, 5.0])  # A/ps
    masses = np.array([127.0, 127.0, 131.0, 131.0, 135.0, 135.0])
    return _make_ion(
        num_molecules=3,
        final_speeds_per_atom=speeds,
        masses_amu_per_atom=masses,
    )


class TestRunSummaryIHeKedSections:
    def test_mean_energy_section_builds(self):
        import matplotlib.pyplot as plt
        mod = _load_run_summary_module()
        ked_ref = load_ihe_ked_reference(REFERENCE_CSV)
        fig = mod._section_ihe_ked_mean_energy(_smoke_ion(), ked_ref)
        assert fig is not None
        plt.close("all")

    def test_curves_sections_build_both_representations(self):
        import matplotlib.pyplot as plt
        mod = _load_run_summary_module()
        ked_ref = load_ihe_ked_reference(REFERENCE_CSV)
        for representation in ("2d", "3d"):
            fig = mod._section_ihe_ked_curves(
                _smoke_ion(), IHE_KED_DIR, ked_ref, representation
            )
            assert fig is not None
            # 2x2 grid, n = 4 panel dropped -- no dead sixth axis.
            assert len(fig.axes) == 4
        plt.close("all")

    def test_mass_spectrum_with_abundance_builds(self):
        import matplotlib.pyplot as plt
        from i2_helium_md.postprocess import load_he_abundance_reference
        mod = _load_run_summary_module()
        abundance = load_he_abundance_reference(
            PROJECT_ROOT / "data" / "reference"
            / "integrated_i_he_abundance.csv"
        )
        fig = mod._section_mass_spectrum(_smoke_ion(), abundance)
        assert fig is not None
        plt.close("all")

    def test_old_vmi_section_is_gone(self):
        mod = _load_run_summary_module()
        assert not hasattr(mod, "_section_radial_velocity")
        assert not hasattr(mod, "VMI_REF_HE_PATH")


class TestNanAwareMovingMean:
    def test_all_finite_matches_convolve_ratio(self):
        mod = _load_run_summary_module()
        rng = np.random.default_rng(0)
        v = rng.normal(size=37)
        window = 5
        expected = np.convolve(
            v, np.ones(window), "same"
        ) / np.convolve(np.ones_like(v), np.ones(window), "same")
        out = mod._nan_aware_moving_mean(v, window)
        assert out == pytest.approx(expected, rel=1e-12)

    def test_all_finite_hand_computed_five_point(self):
        mod = _load_run_summary_module()
        v = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        out = mod._nan_aware_moving_mean(v, 3)
        # centered window=3, 'same'/shrink edges:
        # idx0: mean(1,2) = 1.5; idx1: mean(1,2,3)=2; idx2: mean(2,3,4)=3;
        # idx3: mean(3,4,5)=4; idx4: mean(4,5)=4.5
        expected = np.array([1.5, 2.0, 3.0, 4.0, 4.5])
        assert out == pytest.approx(expected, rel=1e-12)

    def test_even_window_hand_computed_centering(self):
        # Production window is 10 (even) -- pin the even-window centering
        # direction: numpy's 'same' convolution puts window/2 samples BEFORE
        # and window/2 - 1 AFTER the current index, matching MATLAB
        # movmean's documented even-k convention ("centered about the
        # current and previous elements").
        mod = _load_run_summary_module()
        v = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        out = mod._nan_aware_moving_mean(v, 4)
        # idx0: mean(1,2)=1.5; idx1: mean(1,2,3)=2; idx2: mean(1,2,3,4)=2.5;
        # idx3: mean(2,3,4,5)=3.5; idx4: mean(3,4,5)=4
        expected = np.array([1.5, 2.0, 2.5, 3.5, 4.0])
        assert out == pytest.approx(expected, rel=1e-12)

    def test_interior_nan_run_preserves_nan_and_neighbor_means(self):
        mod = _load_run_summary_module()
        v = np.array([1.0, 2.0, np.nan, np.nan, 5.0, 6.0, 7.0])
        out = mod._nan_aware_moving_mean(v, 3)
        nan_mask = np.isnan(v)
        assert np.array_equal(np.isnan(out), nan_mask)
        # idx1 (value 2.0): window covers idx0,1,2 -> finite {1.0, 2.0} only.
        assert out[1] == pytest.approx(1.5, rel=1e-12)
        # idx4 (value 5.0): window covers idx3,4,5 -> finite {5.0, 6.0} only.
        assert out[4] == pytest.approx(5.5, rel=1e-12)

    def test_window_one_returns_input_unchanged(self):
        mod = _load_run_summary_module()
        v = np.array([1.0, np.nan, 3.0, np.nan, 5.0])
        out = mod._nan_aware_moving_mean(v, 1)
        assert np.isnan(out[1]) and np.isnan(out[3])
        finite_mask = np.isfinite(v)
        assert out[finite_mask] == pytest.approx(v[finite_mask], rel=1e-12)

    def test_window_zero_raises(self):
        mod = _load_run_summary_module()
        with pytest.raises(ValueError):
            mod._nan_aware_moving_mean(np.array([1.0, 2.0, 3.0]), 0)
