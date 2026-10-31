"""``legacy_ked_shapes_figure.py`` -- B2-3 on a legacy ion checkpoint, n = 0-3.

The script reuses the B2 data and draw layer (recipe parity is locked in
``test_tier2_b2_thesis_figure``); here: the legacy overrides (panels, ranges,
colour, legend name), that ``pv_data`` reads a plain ``IonCheckpoint`` (no
detection stage), and the points rows. Synthetic checkpoint, committed ihe_ked
reference, Agg backend, figures closed, nothing written to disk.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from i2_helium_md.physics.constants import U as U_KG  # noqa: E402
from i2_helium_md.postprocess import fragment_mean_kinetic_energy, load_ihe_ked_reference  # noqa: E402
from i2_helium_md.simulation.checkpoint import IonCheckpoint  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "post_processing" / "legacy_ked_shapes_figure.py"
IHE_KED_DIR = PROJECT_ROOT / "data" / "reference" / "ihe_ked"


def _make_ion(*, num_molecules: int, final_speeds_per_atom: np.ndarray,
              masses_amu_per_atom: np.ndarray, num_steps: int = 4) -> IonCheckpoint:
    """Tiny IonCheckpoint with prescribed final |v| [Å/ps] and final masses.

    Copied from ``tests/test_ihe_ked.py`` (private test helpers are not shared
    across test modules); every molecule outside the droplet.
    """
    n = num_molecules
    masses_kg = masses_amu_per_atom.astype(float) * U_KG
    zero = np.zeros((2 * n, num_steps))
    return IonCheckpoint(
        num_molecules=n,
        time_ps=np.linspace(0.0, 1.0, num_steps),
        positions_x=zero, positions_y=zero, positions_z=zero,
        velocities_x=zero, velocities_y=zero, velocities_z=zero,
        positions_final_x=np.zeros(2 * n), positions_final_y=np.zeros(2 * n),
        positions_final_z=np.zeros(2 * n),
        velocities_final_x=final_speeds_per_atom.astype(float),
        velocities_final_y=np.zeros(2 * n), velocities_final_z=np.zeros(2 * n),
        mass_kg=masses_kg.copy(), mass_final_kg=masses_kg,
        mass_history_kg=np.broadcast_to(masses_kg[:, None], (2 * n, num_steps)).copy(),
        droplet_radii_angstrom=np.full(2 * n, 30.0),
        E_kin_eV=zero, E_pot_eV=zero, E_dissip_eV=zero, E_mass_transfer_eV=zero,
        E_int_eV=zero, n_shell=zero,
        b_ion_outside=np.ones(n, dtype=bool),
        relative_loss_per_ps=zero,
        number_of_collisions=np.zeros((2 * n, num_steps), dtype=int),
        temperature_diagnostic=np.full((num_steps, 3), np.nan),
    )


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("legacy_ked_shapes_under_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def ion():
    # 8 molecules / 16 ions: four ions per n = 0..3 (legacy integer masses)
    rng = np.random.default_rng(7)
    masses = np.repeat([127.0, 131.0, 135.0, 139.0], 4)
    return _make_ion(num_molecules=8, final_speeds_per_atom=rng.uniform(4.0, 10.0, 16),
                     masses_amu_per_atom=masses)


@pytest.fixture(scope="module")
def ked_ref():
    return load_ihe_ked_reference(IHE_KED_DIR / "IHe_KED_reference.csv")


class TestSettings:
    def test_legacy_overrides_on_b2_3(self, mod):
        st = mod.legacy_style()
        assert st["n_values"] == (0, 1, 2, 3)
        assert st["sim_label"] == "Hard-sphere model"
        assert st["sim"]["color"] == mod.LEGACY_COLOR
        assert set(st["e_xlim_by_n"]) == set(st["n_values"])
        # the B2-3 recipe keys are inherited, not redefined
        for key in ("representation", "x_axis", "e_min_eV", "ylim"):
            assert st[key] == mod.FIG_B2_3[key]

    def test_legacy_colour_is_not_the_b2_colour(self, mod):
        from scripts.post_processing.thesis_style import SIM_COLOR
        assert mod.LEGACY_COLOR != SIM_COLOR


class TestFigure:
    def test_reads_a_plain_ion_checkpoint(self, mod, ion, ked_ref):
        panels = mod.legacy_panels(ion, ked_ref, (0, 1, 2, 3))
        assert [d.n for d in panels] == [0, 1, 2, 3]
        assert [d.num_sim for d in panels] == [4, 4, 4, 4]

    def test_points_rows_are_the_panel_numbers(self, mod, ion, ked_ref):
        panels = mod.legacy_panels(ion, ked_ref, (0, 1, 2, 3))
        rows = mod.points_rows(panels)
        assert [r["n"] for r in rows] == [0, 1, 2, 3]
        for r in rows:
            sim = fragment_mean_kinetic_energy(ion, r["n"]).mean_KE_eV
            assert r["mean_E_sim_eV"] == pytest.approx(sim, rel=1e-12)
            idx = int(np.flatnonzero(np.asarray(ked_ref.n) == r["n"])[0])
            assert r["mean_E_exp_eV"] == pytest.approx(ked_ref.mean_KE_eV[idx])

    def test_four_panels_with_legacy_ranges_and_label(self, mod, ion, ked_ref):
        st = mod.legacy_style()
        fig = mod.build_b2_3(mod.legacy_panels(ion, ked_ref, st["n_values"]), st)
        try:
            axes = fig.axes
            assert len(axes) == 4
            assert [ax.get_xlim() for ax in axes] == [st["e_xlim_by_n"][n]
                                                      for n in st["n_values"]]
            labels = [t.get_text() for leg in fig.legends for t in leg.get_texts()]
            assert "Simulation (Hard-sphere model)" in labels
            sim = [ln for ln in axes[0].lines
                   if ln.get_label() == "Simulation (Hard-sphere model)"]
            assert len(sim) == 1 and sim[0].get_color() == mod.LEGACY_COLOR
        finally:
            plt.close(fig)
