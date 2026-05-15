"""Tests for the paper-v2 I+He comparison post-processing helpers."""

from __future__ import annotations

import json
import importlib.util
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import pytest
from scipy.io import savemat

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.paper_v2 import (
    PAPER_V2_IMAGE_BINS_APS,
    PaperV2RadialReference,
    PaperV2VelocityCurve,
    PaperV2VelocityMap,
    PaperV2VMIImageReference,
    PaperV2VMIPolarImageReference,
    load_paper_v2_phi_reference,
    load_paper_v2_radial_references,
    load_paper_v2_vmi_image_reference,
    load_paper_v2_vmi_polar_image_reference,
    paper_v2_velocity_map,
)
from i2_helium_md.postprocess.paper_v2_plotting import build_radial_figure
from i2_helium_md.simulation.checkpoint import IonCheckpoint

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PAPER_V2_SCRIPT = PROJECT_ROOT / "scripts" / "post_processing" / "plot_paper_v2.py"


def _import_paper_v2_script():
    spec = importlib.util.spec_from_file_location("plot_paper_v2_orientation_test", PAPER_V2_SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _make_ion(
    *,
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    masses_amu: np.ndarray,
    b_outside: np.ndarray | None = None,
) -> IonCheckpoint:
    n2 = vx.size
    if n2 % 2:
        raise AssertionError("Need 2N atoms")
    n = n2 // 2
    if b_outside is None:
        b_outside = np.ones(n, dtype=bool)
    masses_kg = masses_amu.astype(float) * U_KG
    steps = 3
    zeros = np.zeros((n2, steps))
    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, steps),
        positions_x=zeros.copy(),
        positions_y=zeros.copy(),
        positions_z=zeros.copy(),
        velocities_x=zeros.copy(),
        velocities_y=zeros.copy(),
        velocities_z=zeros.copy(),
        positions_final_x=np.zeros(n2),
        positions_final_y=np.zeros(n2),
        positions_final_z=np.zeros(n2),
        velocities_final_x=vx.astype(float),
        velocities_final_y=vy.astype(float),
        velocities_final_z=vz.astype(float),
        mass_kg=masses_kg.copy(),
        mass_final_kg=masses_kg,
        mass_history_kg=np.broadcast_to(masses_kg[:, None], (n2, steps)).copy(),
        droplet_radii_angstrom=np.full(n2, 30.0),
        E_kin_eV=zeros.copy(),
        E_pot_eV=zeros.copy(),
        E_dissip_eV=zeros.copy(),
        E_mass_attach_defect_eV=zeros.copy(),
        b_ion_outside=np.asarray(b_outside, dtype=bool),
        relative_loss_per_ps=zeros.copy(),
        number_of_collisions=np.zeros((n2, steps), dtype=int),
        temperature_diagnostic=np.full((steps, 3), np.nan),
    )


def test_vmi_image_reference_loader_validates_npz_fields_and_shapes(tmp_path):
    path = tmp_path / "iplus_he_600mw_43569_vmi_image.npz"
    np.savez(
        path,
        vx_Aps=np.array([-1.0, 0.0, 1.0]),
        vy_Aps=np.array([-2.0, 2.0]),
        intensity=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
    )
    (tmp_path / "iplus_he_600mw_43569_vmi_image.json").write_text(
        json.dumps({"measurement_id": 43569, "power_mw": 600}),
        encoding="ascii",
    )

    ref = load_paper_v2_vmi_image_reference(path)

    assert isinstance(ref, PaperV2VMIImageReference)
    np.testing.assert_allclose(ref.vx_Aps, [-1.0, 0.0, 1.0])
    np.testing.assert_allclose(ref.vy_Aps, [-2.0, 2.0])
    np.testing.assert_allclose(ref.intensity, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    assert ref.metadata["measurement_id"] == 43569

    missing = tmp_path / "missing_field.npz"
    np.savez(missing, vx_Aps=np.array([0.0]), intensity=np.array([[1.0]]))
    with pytest.raises(ValueError, match="vy_Aps"):
        load_paper_v2_vmi_image_reference(missing)

    bad_shape = tmp_path / "bad_shape.npz"
    np.savez(
        bad_shape,
        vx_Aps=np.array([0.0, 1.0]),
        vy_Aps=np.array([0.0]),
        intensity=np.array([[1.0], [2.0]]),
    )
    with pytest.raises(ValueError, match="intensity shape"):
        load_paper_v2_vmi_image_reference(bad_shape)

    constant_axis = tmp_path / "constant_axis.npz"
    np.savez(
        constant_axis,
        vx_Aps=np.zeros(3),
        vy_Aps=np.array([-1.0, 1.0]),
        intensity=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
    )
    with pytest.raises(ValueError, match="zero range"):
        load_paper_v2_vmi_image_reference(constant_axis)


def test_vmi_image_reference_loader_accepts_matlab_mat_exports(tmp_path):
    path = tmp_path / "iplus_he_high_snr_vmi_image.mat"
    vx_grid, vy_grid = np.meshgrid(
        np.array([-1.0, 0.0, 1.0]),
        np.array([-2.0, 2.0]),
    )
    savemat(
        path,
        {
            "vx_Aps": vx_grid,
            "vy_Aps": vy_grid,
            "intensity": np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]),
        },
    )

    ref = load_paper_v2_vmi_image_reference(path)

    np.testing.assert_allclose(ref.vx_Aps, vx_grid)
    np.testing.assert_allclose(ref.vy_Aps, vy_grid)
    np.testing.assert_allclose(ref.intensity, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])


def test_experimental_image_plot_passes_matplotlib_ready_grids_unchanged(monkeypatch):
    module = _import_paper_v2_script()
    vx_grid, vy_grid = np.meshgrid(
        np.array([-1.0, 0.0, 1.0]),
        np.array([-2.0, 2.0]),
    )
    intensity = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    ref = PaperV2VMIImageReference(
        vx_Aps=vx_grid,
        vy_Aps=vy_grid,
        vx_mps=vx_grid * 100.0,
        vy_mps=vy_grid * 100.0,
        intensity=intensity,
        metadata={},
        source_path=Path("synthetic.mat"),
    )
    fig, ax = plt.subplots()
    captured = {}
    original_pcolormesh = ax.pcolormesh

    def capture_pcolormesh(x, y, c, *args, **kwargs):
        captured["x"] = x
        captured["y"] = y
        captured["c"] = c
        return original_pcolormesh(x, y, c, *args, **kwargs)

    monkeypatch.setattr(ax, "pcolormesh", capture_pcolormesh)
    monkeypatch.setattr(plt, "colorbar", lambda *args, **kwargs: None)

    module._draw_experimental_image(ax, ref)

    np.testing.assert_allclose(captured["x"], vx_grid * 100.0)
    np.testing.assert_allclose(captured["y"], vy_grid * 100.0)
    np.testing.assert_allclose(captured["c"], intensity)
    assert ax.get_xlabel() == "v_x / m/s"
    assert ax.get_ylabel() == "v_y / m/s"
    plt.close(fig)


def test_radial_reference_loader_reads_directory_and_labels(tmp_path):
    ref_dir = tmp_path / "paper_v2"
    ref_dir.mkdir()
    (ref_dir / "iplus_he_160mw_43556_radial.csv").write_text(
        "v_Aps,signal_arb\n0.0,9.0\n1.0,9.0\n",
        encoding="ascii",
    )
    (ref_dir / "iplus_he_600mw_43569_radial.csv").write_text(
        "v_Aps,signal_arb\n0.0,8.0\n1.0,8.0\n",
        encoding="ascii",
    )
    (ref_dir / "iplus_gas_300mw_43562_radial.csv").write_text(
        "v_Aps,signal_arb\n0.0,1.0\n1.0,0.5\n",
        encoding="ascii",
    )
    (ref_dir / "iplus_he_high_snr_radial.csv").write_text(
        "v_Aps,signal_arb\n0.0,2.0\n1.0,1.0\n",
        encoding="ascii",
    )

    refs = load_paper_v2_radial_references(ref_dir)

    assert [ref.label for ref in refs] == [
        "I+ gas 300 mW (43562)",
        "I+He high-SNR",
        "I+He 160 mW (43556)",
        "I+He 600 mW (43569)",
    ]
    np.testing.assert_allclose(refs[1].velocity_Aps, [0.0, 1.0])
    np.testing.assert_allclose(refs[1].signal_arb, [2.0, 1.0])


def test_radial_figure_title_identifies_raw_2d_vmi_comparison():
    radial_ref = PaperV2RadialReference(
        velocity_Aps=np.array([0.0, 1.0]),
        velocity_mps=np.array([0.0, 100.0]),
        signal_arb=np.array([0.0, 1.0]),
        label="raw VMI reference",
        source_path=Path("synthetic_radial.csv"),
    )
    velocity_curve = PaperV2VelocityCurve(
        mass_amu=131.0,
        bin_centers_Aps=np.array([0.0, 1.0]),
        bin_centers_mps=np.array([0.0, 100.0]),
        bin_edges_Aps=np.array([0.0, 0.5, 1.5]),
        bin_edges_mps=np.array([0.0, 50.0, 150.0]),
        counts=np.array([0.0, 1.0]),
        smoothed=np.array([0.0, 1.0]),
        normalised=np.array([0.0, 1.0]),
        num_atoms_used=2,
        smoothing_window=15,
    )

    fig = build_radial_figure([radial_ref], velocity_curve)
    try:
        assert fig.axes[0].get_title() == (
            "2-D detector-plane speed vs raw VMI radial profile"
        )
    finally:
        plt.close(fig)


def test_phi_reference_loader_enforces_columns(tmp_path):
    path = tmp_path / "iplus_he_high_snr_phi.csv"
    path.write_text("phi_rad,signal_arb\n0.0,1.0\n3.14,0.5\n", encoding="ascii")

    ref = load_paper_v2_phi_reference(path)

    np.testing.assert_allclose(ref.phi_rad, [0.0, 3.14])
    np.testing.assert_allclose(ref.signal_arb, [1.0, 0.5])

    bad = tmp_path / "bad_phi.csv"
    bad.write_text("angle,signal\n0.0,1.0\n", encoding="ascii")
    with pytest.raises(ValueError, match="phi_rad.*signal_arb"):
        load_paper_v2_phi_reference(bad)


def test_velocity_map_matches_matlab_nearest_bin_selection():
    ion = _make_ion(
        vx=np.array([-35.0, 0.09, 3.0, 5.0]),
        vy=np.array([-35.0, 0.09, 4.0, 12.0]),
        vz=np.array([999.0, 999.0, 999.0, 999.0]),
        masses_amu=np.array([131.0, 131.0, 127.0, 127.0]),
        b_outside=np.array([True, True]),
    )

    image = paper_v2_velocity_map(ion, mass_amu=131.0)

    np.testing.assert_allclose(image.velocity_bins_Aps, PAPER_V2_IMAGE_BINS_APS)
    assert image.counts.shape == (PAPER_V2_IMAGE_BINS_APS.size, PAPER_V2_IMAGE_BINS_APS.size)
    assert image.num_atoms_used == 2
    assert image.counts.sum() == 2
    idx_neg = int(np.where(np.isclose(PAPER_V2_IMAGE_BINS_APS, -35.0))[0][0])
    idx_zero = int(np.where(np.isclose(PAPER_V2_IMAGE_BINS_APS, 0.0))[0][0])
    assert image.counts[idx_neg, idx_neg] == 1
    assert image.counts[idx_zero, idx_zero] == 1


def test_color_norm_floor_clips_below_fraction_of_max():
    module = _import_paper_v2_script()

    legacy = module._color_norm(np.array([1.0, 10.0, 100.0]))
    assert legacy.vmin == pytest.approx(0.0)
    assert legacy.vmax == pytest.approx(80.0)

    floored = module._color_norm(
        np.array([1.0, 10.0, 100.0]), noise_floor_fraction=0.2
    )
    assert floored.vmin == pytest.approx(20.0)
    assert floored.vmax == pytest.approx(100.0)

    assert module._color_norm(np.array([]), noise_floor_fraction=0.2) is None
    assert module._color_norm(np.array([0.0, 0.0])) is None


def test_polar_image_reference_loader_validates_npz_fields_and_shapes(tmp_path):
    path = tmp_path / "iplus_he_high_snr_vmi_polar_image.npz"
    phi = np.linspace(0.0, 2.0 * np.pi, 6, endpoint=False)
    v_radius_mps = np.array([0.0, 100.0, 200.0, 300.0])
    intensity = np.arange(phi.size * v_radius_mps.size, dtype=float).reshape(
        (phi.size, v_radius_mps.size)
    )
    np.savez(
        path,
        phi_rad=phi,
        v_radius_mps=v_radius_mps,
        intensity_polar=intensity,
    )
    (tmp_path / "iplus_he_high_snr_vmi_polar_image.json").write_text(
        json.dumps({"channel": "I+He high-SNR processed VMI (polar)"}),
        encoding="ascii",
    )

    ref = load_paper_v2_vmi_polar_image_reference(path)

    assert isinstance(ref, PaperV2VMIPolarImageReference)
    np.testing.assert_allclose(ref.phi_rad, phi)
    np.testing.assert_allclose(ref.v_radius_mps, v_radius_mps)
    np.testing.assert_allclose(ref.v_radius_Aps, v_radius_mps / 100.0)
    np.testing.assert_allclose(ref.intensity, intensity)
    assert ref.metadata["channel"] == "I+He high-SNR processed VMI (polar)"

    missing_phi = tmp_path / "missing_phi.npz"
    np.savez(missing_phi, v_radius_mps=v_radius_mps, intensity_polar=intensity)
    with pytest.raises(ValueError, match="phi_rad"):
        load_paper_v2_vmi_polar_image_reference(missing_phi)

    missing_v = tmp_path / "missing_v.npz"
    np.savez(missing_v, phi_rad=phi, intensity_polar=intensity)
    with pytest.raises(ValueError, match="v_radius_mps"):
        load_paper_v2_vmi_polar_image_reference(missing_v)

    bad_shape = tmp_path / "bad_shape.npz"
    np.savez(
        bad_shape,
        phi_rad=phi,
        v_radius_mps=v_radius_mps,
        intensity_polar=intensity.T,
    )
    with pytest.raises(ValueError, match="intensity_polar shape"):
        load_paper_v2_vmi_polar_image_reference(bad_shape)

    constant_v = tmp_path / "constant_v.npz"
    np.savez(
        constant_v,
        phi_rad=phi,
        v_radius_mps=np.zeros_like(v_radius_mps),
        intensity_polar=intensity,
    )
    with pytest.raises(ValueError, match="zero range"):
        load_paper_v2_vmi_polar_image_reference(constant_v)


def test_polar_image_reference_loader_accepts_matlab_mat_exports(tmp_path):
    path = tmp_path / "iplus_he_high_snr_vmi_polar_image.mat"
    phi = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
    v_radius_mps = np.linspace(0.0, 350.0, 5)
    intensity = np.outer(np.cos(phi) ** 2 + 1.0, v_radius_mps + 1.0)
    savemat(
        path,
        {
            "phi_rad": phi,
            "v_radius_mps": v_radius_mps,
            "intensity_polar": intensity,
        },
    )

    ref = load_paper_v2_vmi_polar_image_reference(path)

    np.testing.assert_allclose(ref.phi_rad, phi)
    np.testing.assert_allclose(ref.v_radius_Aps, v_radius_mps / 100.0)
    np.testing.assert_allclose(ref.intensity, intensity)


def test_polar_histogram_helper_matches_reference_axes():
    module = _import_paper_v2_script()
    rng = np.random.default_rng(20260515)
    n_atoms = 200
    speed = rng.uniform(0.0, 2.0, size=n_atoms)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n_atoms)
    vx = speed * np.cos(phi)
    vy = speed * np.sin(phi)
    vz = np.zeros(n_atoms)
    masses_amu = np.full(n_atoms, 131.0)
    ion = _make_ion(
        vx=np.concatenate([vx, vx]),
        vy=np.concatenate([vy, vy]),
        vz=np.concatenate([vz, vz]),
        masses_amu=np.concatenate([masses_amu, masses_amu]),
        b_outside=np.ones(n_atoms, dtype=bool),
    )
    phi_rad = np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
    v_radius_Aps = np.linspace(0.0, 1.8, 9)
    polar_ref = PaperV2VMIPolarImageReference(
        phi_rad=phi_rad,
        v_radius_Aps=v_radius_Aps,
        v_radius_mps=v_radius_Aps * 100.0,
        intensity=np.zeros((phi_rad.size, v_radius_Aps.size)),
        metadata={},
        source_path=Path("synthetic_polar.npz"),
    )

    polar_hist = module._polar_histogram_matched_to_reference(ion, polar_ref)

    assert polar_hist.phi_centers_rad.size == phi_rad.size
    assert polar_hist.v_centers_Aps.size == v_radius_Aps.size
    expected_v_max = float(v_radius_Aps[-1] + (v_radius_Aps[1] - v_radius_Aps[0]))
    assert polar_hist.v_edges_Aps[-1] == pytest.approx(expected_v_max)
    assert polar_hist.counts.shape == (v_radius_Aps.size, phi_rad.size)


def test_polar_image_figure_renders_two_panels_with_aligned_axes():
    module = _import_paper_v2_script()
    phi = np.linspace(0.0, 2.0 * np.pi, 6, endpoint=False)
    v_radius_Aps = np.linspace(0.0, 1.0, 5)
    intensity = np.outer(np.cos(phi) ** 2, v_radius_Aps + 0.5)
    polar_ref = PaperV2VMIPolarImageReference(
        phi_rad=phi,
        v_radius_Aps=v_radius_Aps,
        v_radius_mps=v_radius_Aps * 100.0,
        intensity=intensity,
        metadata={},
        source_path=Path("synthetic_polar.npz"),
    )

    rng = np.random.default_rng(0)
    n_atoms = 80
    speed = rng.uniform(0.0, 0.9, size=n_atoms)
    angle = rng.uniform(0.0, 2.0 * np.pi, size=n_atoms)
    ion = _make_ion(
        vx=np.concatenate([speed * np.cos(angle), speed * np.cos(angle)]),
        vy=np.concatenate([speed * np.sin(angle), speed * np.sin(angle)]),
        vz=np.zeros(2 * n_atoms),
        masses_amu=np.full(2 * n_atoms, 131.0),
        b_outside=np.ones(n_atoms, dtype=bool),
    )
    polar_hist = module._polar_histogram_matched_to_reference(ion, polar_ref)

    fig = module._build_polar_image_figure(polar_ref, polar_hist)
    try:
        ax_exp, ax_sim = fig.axes[:2]
        assert ax_exp.get_xlim() == pytest.approx((0.0, 2.0 * np.pi))
        assert ax_sim.get_xlim() == pytest.approx((0.0, 2.0 * np.pi))
        assert ax_exp.get_ylim() == pytest.approx((0.0, polar_ref.v_radius_mps[-1]))
        assert ax_sim.get_ylim() == pytest.approx((0.0, polar_ref.v_radius_mps[-1]))
        assert "experimental polar VMI" in ax_exp.get_title()
        assert "simulated polar histogram" in ax_sim.get_title()
        assert ax_exp.collections, "experimental panel should plot a pcolormesh"
        assert ax_sim.collections, "simulated panel should plot a pcolormesh"
    finally:
        plt.close(fig)


def test_simulated_map_plot_transposes_internal_storage_for_physical_axes():
    module = _import_paper_v2_script()
    bins = np.array([-1.0, 0.0, 1.0])
    counts = np.array(
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]
    )
    velocity_map = PaperV2VelocityMap(
        velocity_bins_Aps=bins,
        counts=counts,
        mass_amu=131.0,
        num_atoms_used=9,
    )
    fig, ax = plt.subplots()

    module._draw_simulated_map(ax, velocity_map)

    plotted = ax.collections[0].get_array()
    np.testing.assert_allclose(np.asarray(plotted).reshape(counts.shape), counts.T)
    plt.close(fig)
