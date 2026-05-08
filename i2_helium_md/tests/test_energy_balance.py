"""Tests for the energy_balance recipe helpers and ion checkpoint v5."""

from __future__ import annotations

import numpy as np
import pytest

from i2_helium_md.physics.collisions import (
    CollisionDiagnostics,
    apply_collision,
    temperature_diagnostic_from_collision,
)
from i2_helium_md.physics.constants import U
from i2_helium_md.postprocess import (
    ion_energy_totals,
    mass_spectrum,
    neutral_energy_totals,
    phi_histogram,
)
from i2_helium_md.simulation.checkpoint import (
    IonCheckpoint,
    NeutralCheckpoint,
    load_ion_checkpoint,
    save_ion_checkpoint,
)


# ===========================================================================
# Helper builders
# ===========================================================================
def _neutral(n: int = 2, t: int = 3) -> NeutralCheckpoint:
    rng = np.random.default_rng(0)
    return NeutralCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, t),
        positions_x=np.zeros((2 * n, t)),
        positions_y=np.zeros((2 * n, t)),
        positions_z=np.zeros((2 * n, t)),
        velocities_x=np.zeros((2 * n, t)),
        velocities_y=np.zeros((2 * n, t)),
        velocities_z=np.zeros((2 * n, t)),
        mass_kg=np.full(2 * n, 127 * U),
        droplet_radii=np.full(2 * n, 30.0),
        r0=np.full(n, 5.0),
        E_kin_eV=rng.standard_normal((2 * n, t)),
        E_pot_eV=rng.standard_normal((2 * n, t)),
        E_initial_eV=rng.standard_normal(n),
        E_dissip_eV=rng.standard_normal((2 * n, t)),
        L_droplet_eV_ps=np.zeros((2 * n, t)),
    )


def _ion(
    n: int = 2,
    t: int = 3,
    *,
    velocities_final_x: np.ndarray | None = None,
    velocities_final_y: np.ndarray | None = None,
    mass_final_amu: np.ndarray | None = None,
    temperature_diagnostic: np.ndarray | None = None,
) -> IonCheckpoint:
    rng = np.random.default_rng(1)
    if velocities_final_x is None:
        velocities_final_x = np.zeros(2 * n)
    if velocities_final_y is None:
        velocities_final_y = np.zeros(2 * n)
    if mass_final_amu is None:
        mass_final_amu = np.full(2 * n, 127.0)
    if temperature_diagnostic is None:
        temperature_diagnostic = np.full((t, 3), np.nan, dtype=float)

    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, t),
        positions_x=np.zeros((2 * n, t)),
        positions_y=np.zeros((2 * n, t)),
        positions_z=np.zeros((2 * n, t)),
        velocities_x=np.zeros((2 * n, t)),
        velocities_y=np.zeros((2 * n, t)),
        velocities_z=np.zeros((2 * n, t)),
        positions_final_x=np.zeros(2 * n),
        positions_final_y=np.zeros(2 * n),
        positions_final_z=np.zeros(2 * n),
        velocities_final_x=velocities_final_x,
        velocities_final_y=velocities_final_y,
        velocities_final_z=np.zeros(2 * n),
        mass_kg=np.full(2 * n, 127 * U),
        mass_final_kg=mass_final_amu * U,
        mass_history_kg=np.full((2 * n, t), 127 * U),
        droplet_radii_angstrom=np.full(2 * n, 30.0),
        E_kin_eV=rng.standard_normal((2 * n, t)),
        E_pot_eV=rng.standard_normal((2 * n, t)),
        E_dissip_eV=rng.standard_normal((2 * n, t)),
        E_mass_attach_defect_eV=rng.standard_normal((2 * n, t)),
        b_ion_outside=np.zeros(n, dtype=bool),
        relative_loss_per_ps=np.zeros((2 * n, t)),
        number_of_collisions=np.zeros((2 * n, t), dtype=int),
        temperature_diagnostic=temperature_diagnostic,
    )


# ===========================================================================
# Energy totals
# ===========================================================================
class TestEnergyTotals:
    def test_neutral_E_system_is_sum_of_components(self):
        ckpt = _neutral(n=3, t=4)
        totals = neutral_energy_totals(ckpt)
        np.testing.assert_allclose(
            totals.E_system_eV,
            totals.E_kin_eV + totals.E_pot_eV + totals.E_dissip_eV,
        )
        assert totals.E_kin_eV.shape == (4,)
        assert totals.E_mass_attach_defect_eV is None

    def test_ion_totals_per_molecule_division(self):
        ckpt = _ion(n=2, t=4)
        totals = ion_energy_totals(ckpt)

        np.testing.assert_allclose(
            totals.E_kin_eV, np.sum(ckpt.E_kin_eV, axis=0) / 2.0,
        )
        np.testing.assert_allclose(
            totals.E_system_eV,
            totals.E_kin_eV + totals.E_pot_eV
            + totals.E_dissip_eV + totals.E_mass_attach_defect_eV,
        )


# ===========================================================================
# Phi histogram
# ===========================================================================
class TestPhiHistogram:
    def test_uniform_angles_close_to_uniform_density(self):
        n = 5
        # Cover full circle: phi_sim = atan2(vy,vx) + pi.
        thetas = np.linspace(-np.pi, np.pi, 2 * n, endpoint=False)
        ckpt = _ion(
            n=n, t=2,
            velocities_final_x=np.cos(thetas),
            velocities_final_y=np.sin(thetas),
        )
        h = phi_histogram(ckpt, bin_width_rad=2.0 * np.pi / 10.0)
        assert h.bin_centers_rad.size == 10
        assert h.counts.sum() == 2 * n
        # Density integrates (sum * width) to ~1.
        np.testing.assert_allclose(
            h.density.sum() * (2 * np.pi / 10.0), 1.0, atol=1e-12,
        )

    def test_mass_selection_filters_atoms(self):
        n = 4
        masses = np.array([127.0, 131.0, 131.0, 127.0,
                           127.0, 131.0, 131.0, 131.0])
        ckpt = _ion(
            n=n, t=2,
            velocities_final_x=np.ones(2 * n),
            velocities_final_y=np.zeros(2 * n),
            mass_final_amu=masses,
        )
        h_all = phi_histogram(ckpt, mass_amu=None)
        h_131 = phi_histogram(ckpt, mass_amu=131.0)
        assert h_all.num_atoms_used == 2 * n
        assert h_131.num_atoms_used == int(np.sum(masses == 131.0))


# ===========================================================================
# Mass spectrum
# ===========================================================================
class TestMassSpectrum:
    def test_known_masses_bin_correctly(self):
        masses = np.array([127.0, 127.1, 130.9, 131.0, 134.5, 135.0])
        ckpt = _ion(n=3, t=2, mass_final_amu=masses)
        spec = mass_spectrum(ckpt, bin_width_amu=1.0)

        # Peaks at 127 (2 atoms), 131 (2 atoms), 135 (2 atoms).
        for peak, expected in [(127.0, 2), (131.0, 2), (135.0, 2)]:
            idx = int(np.argmin(np.abs(spec.bin_centers_amu - peak)))
            assert spec.counts[idx] == expected, (
                f"peak {peak}: expected {expected}, got {spec.counts[idx]}"
            )
        assert int(spec.counts.sum()) == masses.size


# ===========================================================================
# Temperature diagnostic recipe
# ===========================================================================
class TestTemperatureDiagnosticRecipe:
    def test_no_collision_returns_nan(self):
        n = 4
        diag = CollisionDiagnostics(
            b_collision=np.zeros(n, dtype=bool),
            COSTHETA=np.ones(n),
            COStheta_lab=np.ones(n),
            rho=np.full(n, 0.25),
            E0_eV=np.ones(n),
            E1_eV=np.ones(n),
        )
        td = temperature_diagnostic_from_collision(diag)
        assert td.shape == (3,)
        assert np.all(np.isnan(td))

    def test_three_components_match_scalar_formula(self):
        # Two atoms collide, two do not. The diagnostic uses lab-frame
        # cos(theta), not COM-frame -- so we feed COStheta_lab values
        # picked to produce easy-to-verify means.
        b = np.array([True, False, True, False])
        E0 = np.array([1.0, 1.0, 2.0, 1.0])
        E1 = np.array([0.7, 1.0, 1.0, 1.0])
        rho = np.array([0.5, 1.0, 0.25, 1.0])
        cos_theta_com = np.array([0.0, 1.0, -0.5, 1.0])
        cos_theta_lab = np.array([0.5, 1.0, 0.99, 1.0])

        diag = CollisionDiagnostics(
            b_collision=b,
            COSTHETA=cos_theta_com,
            COStheta_lab=cos_theta_lab,
            rho=rho,
            E0_eV=E0, E1_eV=E1,
        )
        td = temperature_diagnostic_from_collision(diag)

        expected_T_actual = np.mean([0.7 / 1.0, 1.0 / 2.0])
        expected_T_mass = np.mean(
            [(1 + 0.5 ** 2) / (1 + 0.5) ** 2,
             (1 + 0.25 ** 2) / (1 + 0.25) ** 2]
        )
        expected_theta = np.mean([np.arccos(0.5), np.arccos(0.99)])

        np.testing.assert_allclose(td[0], expected_T_actual)
        np.testing.assert_allclose(td[1], expected_T_mass)
        np.testing.assert_allclose(td[2], expected_theta)

    def test_apply_collision_diagnostic_path_returns_extra_tuple_element(self):
        rng = np.random.default_rng(0)
        n = 6
        v = rng.standard_normal(n) * 5.0
        masses = np.full(n, 127.0)
        b = np.array([True, True, False, False, True, False])

        result4 = apply_collision(
            vx=v, vy=v, vz=v,
            masses_amu=masses,
            b_collision=b,
            scatter_mass_amu=4.0,
            rng=np.random.default_rng(0),
        )
        result5 = apply_collision(
            vx=v, vy=v, vz=v,
            masses_amu=masses,
            b_collision=b,
            scatter_mass_amu=4.0,
            rng=np.random.default_rng(0),
            return_diagnostics=True,
        )
        assert len(result4) == 4
        assert len(result5) == 5
        coll_diag = result5[4]
        assert isinstance(coll_diag, CollisionDiagnostics)
        np.testing.assert_array_equal(coll_diag.b_collision, b)


# ===========================================================================
# Schema v5 round-trip and v4 reject
# ===========================================================================
class TestSchemaV5:
    def test_v5_round_trip_preserves_temperature_diagnostic(self, tmp_path):
        rng = np.random.default_rng(0)
        td = rng.standard_normal((3, 3))
        td[1, :] = np.nan  # row with no collision
        ckpt = _ion(n=2, t=3, temperature_diagnostic=td)

        path = save_ion_checkpoint(ckpt, tmp_path / "i.npz")
        loaded = load_ion_checkpoint(path)
        # NaN-aware compare
        a = ckpt.temperature_diagnostic
        b = loaded.temperature_diagnostic
        assert a.shape == b.shape
        eq = (np.isnan(a) & np.isnan(b)) | (a == b)
        assert eq.all()

    def test_temperature_diagnostic_wrong_shape_rejected(self, tmp_path):
        from i2_helium_md import single_pulse_N2000
        cfg = single_pulse_N2000(seed=0)
        # Build a valid v5 ion checkpoint matching cfg.num_molecules but
        # with a bogus temperature_diagnostic shape.
        n = cfg.num_molecules
        t = 4
        bad = _ion(n=n, t=t)
        # Re-build with malformed td shape.
        bad_dict = {
            "num_molecules": bad.num_molecules,
            "time_ps": bad.time_ps,
            "positions_x": bad.positions_x,
            "positions_y": bad.positions_y,
            "positions_z": bad.positions_z,
            "velocities_x": bad.velocities_x,
            "velocities_y": bad.velocities_y,
            "velocities_z": bad.velocities_z,
            "positions_final_x": bad.positions_final_x,
            "positions_final_y": bad.positions_final_y,
            "positions_final_z": bad.positions_final_z,
            "velocities_final_x": bad.velocities_final_x,
            "velocities_final_y": bad.velocities_final_y,
            "velocities_final_z": bad.velocities_final_z,
            "mass_kg": bad.mass_kg,
            "mass_final_kg": bad.mass_final_kg,
            "mass_history_kg": bad.mass_history_kg,
            "droplet_radii_angstrom": bad.droplet_radii_angstrom,
            "E_kin_eV": bad.E_kin_eV,
            "E_pot_eV": bad.E_pot_eV,
            "E_dissip_eV": bad.E_dissip_eV,
            "E_mass_attach_defect_eV": bad.E_mass_attach_defect_eV,
            "b_ion_outside": bad.b_ion_outside,
            "relative_loss_per_ps": bad.relative_loss_per_ps,
            "number_of_collisions": bad.number_of_collisions,
            "temperature_diagnostic": np.zeros((t, 2)),  # wrong: must be (t, 3)
        }
        ckpt_bad = IonCheckpoint(**bad_dict)
        path = save_ion_checkpoint(ckpt_bad, tmp_path / "bad.npz")
        with pytest.raises(ValueError, match="temperature_diagnostic"):
            load_ion_checkpoint(path, cfg=cfg)

    def test_v4_file_missing_temperature_diagnostic_rejected(self, tmp_path):
        ckpt = _ion(n=2, t=3)
        # Save normally, then strip the temperature_diagnostic field
        # from the .npz to mimic a v4 file. We also bump schema_version
        # back to 4 so the loader's version check fires first.
        path = save_ion_checkpoint(ckpt, tmp_path / "v4.npz")

        with np.load(path, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files if k != "temperature_diagnostic"}
        data["schema_version"] = np.asarray(4)
        np.savez_compressed(path, **data)

        with pytest.raises(ValueError, match="schema_version"):
            load_ion_checkpoint(path)
