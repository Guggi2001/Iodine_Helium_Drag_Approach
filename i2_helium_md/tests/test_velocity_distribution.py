"""Tests for i2_helium_md/postprocess/velocity_distribution.py."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.physics.constants import U as U_KG
from i2_helium_md.postprocess.velocity_distribution import (
    FinalVelocityHistogram,
    VmiReference,
    compute_final_velocity_histogram,
    load_vmi_reference,
)
from i2_helium_md.simulation.checkpoint import IonCheckpoint


PROJECT_ROOT = Path(__file__).resolve().parents[1]
VMI_HE = PROJECT_ROOT / "data" / "reference" / "vmi_iplus_he.csv"
VMI_GAS = PROJECT_ROOT / "data" / "reference" / "vmi_iplus_gas.csv"


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
        E_mass_attach_defect_eV=diag_zero,
        b_ion_outside=np.asarray(b_outside, dtype=bool),
        relative_loss_per_ps=diag_zero,
        number_of_collisions=np.zeros((2 * n, num_steps), dtype=int),
    )


# ===========================================================================
# VMI loader
# ===========================================================================
class TestLoadVmiReference:
    def test_round_trip_synthetic(self, tmp_path):
        path = tmp_path / "vmi.csv"
        path.write_text(
            "v_Aps,signal_arb\n0.5,1.2\n1.0,2.4\n1.5,3.6\n",
            encoding="ascii",
        )
        ref = load_vmi_reference(path)

        assert isinstance(ref, VmiReference)
        np.testing.assert_array_equal(ref.velocity_Aps, [0.5, 1.0, 1.5])
        np.testing.assert_array_equal(ref.signal_arb, [1.2, 2.4, 3.6])
        assert ref.source_path == path.resolve()

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="VMI reference"):
            load_vmi_reference(tmp_path / "nope.csv")

    def test_bad_header_raises(self, tmp_path):
        path = tmp_path / "vmi.csv"
        path.write_text("velocity,intensity\n0.5,1.0\n", encoding="ascii")
        with pytest.raises(ValueError, match="expected.*v_Aps.*signal_arb"):
            load_vmi_reference(path)


@pytest.mark.skipif(
    not (VMI_HE.exists() and VMI_GAS.exists()),
    reason="real VMI references not present",
)
class TestRealVmi:
    def test_real_vmi_he_loads(self):
        ref = load_vmi_reference(VMI_HE)
        assert ref.velocity_Aps.size > 100
        assert ref.signal_arb.size == ref.velocity_Aps.size
        assert np.all(np.isfinite(ref.velocity_Aps))

    def test_real_vmi_gas_loads(self):
        ref = load_vmi_reference(VMI_GAS)
        assert ref.velocity_Aps.size > 100
        assert ref.signal_arb.size == ref.velocity_Aps.size


# ===========================================================================
# Histogram
# ===========================================================================
class TestComputeFinalVelocityHistogram:
    def test_filters_by_mass(self):
        """Half the atoms at 131 amu, half at 135 amu; each call selects one half."""
        n = 6
        speeds = np.full(2 * n, 5.0)  # all the same speed for simplicity
        masses = np.array(
            [131, 131, 131, 135, 135, 135,    # I1 atoms
             131, 131, 131, 135, 135, 135],   # I2 atoms (same molecule layout)
            dtype=float,
        )
        ion = _make_ion(
            num_molecules=n,
            final_speeds_per_atom=speeds,
            masses_amu_per_atom=masses,
        )

        h131 = compute_final_velocity_histogram(
            ion, mass_amu=131.0, num_bins=14, v_max_Aps=14.0
        )
        h135 = compute_final_velocity_histogram(
            ion, mass_amu=135.0, num_bins=14, v_max_Aps=14.0
        )

        assert h131.num_atoms_used == 6
        assert h135.num_atoms_used == 6
        assert h131.counts.sum() == 6
        assert h135.counts.sum() == 6
        # All speeds == 5.0, bin width 1.0 -> the bin covering [5,6) holds them.
        assert h131.counts[5] == 6
        assert h135.counts[5] == 6

    def test_respects_outside_flag(self):
        """b_ion_outside = [True, False] -> only atoms of molecule 0 contribute."""
        n = 2
        speeds = np.array([3.0, 4.0,    # I1 atoms
                           5.0, 6.0])   # I2 atoms
        masses = np.full(2 * n, 131.0)
        ion = _make_ion(
            num_molecules=n,
            final_speeds_per_atom=speeds,
            masses_amu_per_atom=masses,
            b_outside=np.array([True, False]),
        )

        with_outside = compute_final_velocity_histogram(
            ion, mass_amu=131.0, num_bins=10, v_max_Aps=10.0,
            require_outside=True,
        )
        without_outside = compute_final_velocity_histogram(
            ion, mass_amu=131.0, num_bins=10, v_max_Aps=10.0,
            require_outside=False,
        )

        # With outside: only molecule 0 -> atoms 0 (I1, |v|=3.0) and 2 (I2, |v|=5.0).
        assert with_outside.num_atoms_used == 2
        assert with_outside.counts.sum() == 2

        # Without: all 4 atoms.
        assert without_outside.num_atoms_used == 4
        assert without_outside.counts.sum() == 4

    def test_density_is_counts_over_bin_width(self):
        n = 4
        speeds = np.full(2 * n, 5.0)
        masses = np.full(2 * n, 131.0)
        ion = _make_ion(
            num_molecules=n,
            final_speeds_per_atom=speeds,
            masses_amu_per_atom=masses,
        )

        h = compute_final_velocity_histogram(
            ion, mass_amu=131.0, num_bins=10, v_max_Aps=10.0
        )

        bin_width = h.bin_edges_Aps[1] - h.bin_edges_Aps[0]
        np.testing.assert_allclose(h.density, h.counts / bin_width)

    def test_no_atoms_match_mass_raises(self):
        n = 2
        ion = _make_ion(
            num_molecules=n,
            final_speeds_per_atom=np.full(2 * n, 5.0),
            masses_amu_per_atom=np.full(2 * n, 131.0),
        )
        with pytest.raises(ValueError, match="No atoms passed"):
            compute_final_velocity_histogram(ion, mass_amu=200.0)

    def test_bin_centers_and_edges_consistent(self):
        n = 1
        ion = _make_ion(
            num_molecules=n,
            final_speeds_per_atom=np.array([5.0, 5.0]),
            masses_amu_per_atom=np.full(2, 131.0),
        )
        h = compute_final_velocity_histogram(
            ion, mass_amu=131.0, num_bins=8, v_max_Aps=16.0
        )
        np.testing.assert_allclose(
            h.bin_centers_Aps,
            0.5 * (h.bin_edges_Aps[:-1] + h.bin_edges_Aps[1:]),
        )

    def test_invalid_args_raise(self):
        n = 1
        ion = _make_ion(
            num_molecules=n,
            final_speeds_per_atom=np.array([5.0, 5.0]),
            masses_amu_per_atom=np.full(2, 131.0),
        )
        with pytest.raises(ValueError, match="num_bins"):
            compute_final_velocity_histogram(ion, mass_amu=131.0, num_bins=0)
        with pytest.raises(ValueError, match="v_max_Aps"):
            compute_final_velocity_histogram(
                ion, mass_amu=131.0, v_max_Aps=0.0
            )

    def test_dataclass_type(self):
        n = 1
        ion = _make_ion(
            num_molecules=n,
            final_speeds_per_atom=np.array([5.0, 5.0]),
            masses_amu_per_atom=np.full(2, 131.0),
        )
        h = compute_final_velocity_histogram(ion, mass_amu=131.0)
        assert isinstance(h, FinalVelocityHistogram)
        assert h.mass_amu == 131.0
