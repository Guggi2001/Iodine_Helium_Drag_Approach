"""``legacy_ke_size_figure.py`` -- mean KE and size bars of a legacy run, n = 0 included.

The per-n numbers come from the shared mass-gate helpers
(``fragment_mean_kinetic_energy`` / ``fragment_gate_counts``) and the KE panel
is the B2-2 ``draw_ke_mean``; here: the data products on a synthetic ion
checkpoint (n = 0 present, a fragment beyond the window), the window
normalisation, and that the figure starts both panels at n = 0. Committed
ihe_ked / abundance references; Agg backend, figures closed, nothing written.
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

from i2_helium_md.postprocess import (  # noqa: E402
    fragment_mean_kinetic_energy,
    load_he_abundance_reference,
    load_ihe_ked_reference,
)
from tests.test_legacy_ked_shapes_figure import _make_ion  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "post_processing" / "legacy_ke_size_figure.py"
IHE_KED_DIR = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
ABUNDANCE = PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
N_MAX = 20


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("legacy_ke_size_under_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def ion():
    # 6 molecules / 12 ions: n = 0 x4, 1 x3, 2 x2, 5 x1 (single-ion bin), and
    # two n = 30 ions beyond the n <= 20 window (legacy integer masses 127 + 4n)
    n = np.array([0, 0, 0, 0, 1, 1, 1, 2, 2, 5, 30, 30])
    rng = np.random.default_rng(3)
    return _make_ion(num_molecules=6, final_speeds_per_atom=rng.uniform(3.0, 9.0, 12),
                     masses_amu_per_atom=127.0 + 4.0 * n)


@pytest.fixture(scope="module")
def refs():
    return (load_ihe_ked_reference(IHE_KED_DIR / "IHe_KED_reference.csv"),
            load_he_abundance_reference(ABUNDANCE))


class TestData:
    def test_ke_stops_where_the_experiment_does(self, mod, refs):
        # ions at n = 17 / 18 / 20: KE keeps n <= 17 (the reference's last row),
        # the size bars keep the whole n <= 20 window
        n = np.array([0, 17, 18, 20])
        ion = _make_ion(num_molecules=2, final_speeds_per_atom=np.full(4, 5.0),
                        masses_amu_per_atom=127.0 + 4.0 * n)
        assert mod.ke_n_max(refs[0]) == 17
        assert mod.legacy_ke_data(ion, refs[0]).n.tolist() == [0, 17]
        size = mod.legacy_size_data(ion, refs[1], N_MAX)
        assert size.count[[18, 20]].tolist() == [1, 1]

    def test_ke_starts_at_n0_and_is_the_gate_mean(self, mod, ion, refs):
        ke = mod.legacy_ke_data(ion, refs[0])
        assert ke.n.tolist() == [0, 1, 2, 5]           # empty gates skipped, n = 30 out
        for k, mean, se in zip(ke.n, ke.mean_eV, ke.se_eV):
            ref = fragment_mean_kinetic_energy(ion, int(k))
            assert mean == pytest.approx(ref.mean_KE_eV, rel=1e-12)
            assert se == pytest.approx(ref.stat_err_mean_KE_eV, rel=1e-12)
        # single-ion bin is a display-only thin bin (B2-2 convention)
        assert ke.scored_mask.tolist() == [True, True, True, False]

    def test_size_fractions_over_the_window(self, mod, ion, refs):
        size = mod.legacy_size_data(ion, refs[1], N_MAX)
        assert size.n.tolist() == list(range(N_MAX + 1))
        assert size.frac_sim.sum() == pytest.approx(1.0)
        assert size.frac_exp.sum() == pytest.approx(1.0)
        assert size.frac_sim[0] == pytest.approx(4 / 10)        # 10 ions in the window
        assert size.frac_sim_all_ions[0] == pytest.approx(4 / 12)
        assert 1.0 - size.frac_sim_all_ions.sum() == pytest.approx(2 / 12)
        # the experimental n = 0 share is the abundance reference's, untouched
        assert size.frac_exp[0] == pytest.approx(0.435234834342411)

    def test_points_rows(self, mod, ion, refs):
        ke = mod.legacy_ke_data(ion, refs[0])
        rows = mod.points_rows(ke, mod.legacy_size_data(ion, refs[1], N_MAX))
        assert len(rows) == N_MAX + 1 and rows[0]["n"] == 0
        assert rows[0]["mean_E_exp_eV"] == pytest.approx(3.59357)
        assert rows[3]["mean_E_sim_eV"] == ""                    # empty gate n = 3
        assert rows[N_MAX]["mean_E_exp_eV"] == ""                # exp. KE ends at n = 17


class TestFigure:
    def test_both_panels_start_at_n0(self, mod, ion, refs):
        ke = mod.legacy_ke_data(ion, refs[0])
        size = mod.legacy_size_data(ion, refs[1], N_MAX)
        fig = mod.build_legacy_ke_size(ke, size)
        try:
            ax_e, ax_s = fig.axes
            assert ax_e.get_shared_x_axes().joined(ax_e, ax_s)
            assert ax_e.get_yscale() == "log"
            # experiment errorbar includes the bare-ion row (n_min = 0)
            exp_x = [c.lines[0].get_xdata() for c in ax_e.containers
                     if c.get_label() == r"Experiment $\bar{E}$"]
            assert len(exp_x) == 1 and exp_x[0][0] == 0
            heights = [p.get_height() for p in ax_s.patches]
            np.testing.assert_allclose(heights, 100.0 * np.concatenate(
                [size.frac_exp, size.frac_sim]))
            assert ax_s.get_xlim()[0] < 0
        finally:
            plt.close(fig)

    def test_split_layout(self, mod, ion, refs):
        ke = mod.legacy_ke_data(ion, refs[0])
        size = mod.legacy_size_data(ion, refs[1], N_MAX)
        fig = mod.build_legacy_ke_split_size(ke, size)
        try:
            a, b, c = fig.axes
            split = mod.FIG_LEGACY_KE_SIZE["split_n"]
            assert a.get_yscale() == b.get_yscale() == "linear"
            assert a.get_title() == f"$n$ = 0–{split - 1}"
            assert b.get_title() == f"$n$ = {split}–17"

            def sim_n(ax):
                return [x for cont in ax.containers
                        if cont.get_label().startswith("Simulation (")
                        for x in cont.lines[0].get_xdata()]
            # n = 0, 1, 2 scored in (a); the single-ion n = 5 bin is a thin
            # (open) marker in (b), drawn as a plain line
            assert sim_n(a) == [0, 1, 2] and sim_n(b) == []
            assert any(list(ln.get_xdata()) == [5] for ln in b.lines)
            assert len(c.patches) == 2 * (N_MAX + 1)          # the size bars
            assert len(fig.legends) == 1                       # shared (a) + (b) legend
        finally:
            plt.close(fig)
