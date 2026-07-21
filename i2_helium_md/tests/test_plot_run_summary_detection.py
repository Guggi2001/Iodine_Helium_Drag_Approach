"""Detection-stage sections of ``plot_run_summary.py`` (adjudication D4,
option (a) — log entry "PRODUCTION POINT ADJUDICATED", 2026-07-21).

Tier-2 confirmation runs carry ``relaxation.npz``/``detection.npz``; the
run summary gains three sections rendering the **detected** read
(retained-excluded, solvated renormalized, n = 1 median anchor) with every
physics convention imported from ``i2_helium_md.postprocess
.tier2_confirmation`` (rule 1 — no duplicated conventions):

* ``detected_size_distribution`` — scored histogram (incl. the n = 0
  suppressed/bare bin) + solvated branch vs the abundance reference;
* ``detected_ihe_ked_mean_energy`` — detected per-n mean KE vs the
  ihe_ked reference under the committed error model, χ²_med and the
  mean-legacy χ² annotated;
* ``detected_ke_anatomy`` — per-bin profiled residuals z (the §4z read).

Gating contract (load-bearing): a run dir **without** ``detection.npz``
renders exactly as before — the loader helper returns ``None`` and no
detected section enters the section list. Legacy runs are byte-identical.

Synthetic inputs reuse the ``test_tier2_confirmation`` builders — tiny
in-memory ``DetectionResult`` / reference objects, no production
artifacts, no figures written (Agg backend, figures closed).
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

from tests import test_tier2_confirmation as t2c  # noqa: E402  (builders)

from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
)
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    save_detection_result,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "post_processing" / "plot_run_summary.py"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location(
        "plot_run_summary_detection_tests", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["plot_run_summary_detection_tests"] = module
    spec.loader.exec_module(module)
    return module


def _synthetic_detection():
    """12 ions: 1 retained (excluded), 2 suppressed (bin 0), 9 solvated —
    n = 1..3 carry >= 2 ions each so the KE scorer (min_count = 2) keeps
    them; KE values are unit-scale eV."""
    return t2c.make_detection(
        n_detected=[21.0, 21.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0,
                    3.0, 9.0],
        state_reason=["suppressed", "suppressed"] + ["frozen"] * 9
        + ["droplet_retained"],
        ke_eV=[3.0, 3.1, 1.00, 1.05, 0.95, 1.10, 0.70, 0.72, 0.68, 0.50,
               0.52, 42.0],
    )


def _synthetic_read():
    return read_confirmation_detection(_synthetic_detection(), label="syn")


def _synthetic_refs():
    abundance = t2c.make_abundance_reference(
        n=[0, 1, 2, 3, 4],
        fraction=[0.30, 0.30, 0.20, 0.15, 0.05],
    )
    ked = t2c.make_ked_reference(
        n=[1, 2, 3],
        mean_KE_eV=[1.10, 0.70, 0.50],
        median_KE_eV=[0.98, 0.70, 0.50],
    )
    return abundance, ked


class TestDetectionLoaderGating:
    def test_returns_none_without_detection_artifact(self, mod, tmp_path):
        assert mod._load_detection_read(tmp_path) is None

    def test_reads_detection_artifact_retained_excluded(self, mod, tmp_path):
        save_detection_result(_synthetic_detection(), tmp_path / "detection.npz")
        read = mod._load_detection_read(tmp_path)
        assert read is not None
        assert read.num_ions == 12
        assert read.num_scored == 11  # droplet_retained excluded
        assert read.fraction[0] == pytest.approx(2.0 / 11.0)  # suppressed


class TestDetectedSectionBuilders:
    def test_size_distribution_renders_with_reference(self, mod):
        abundance, _ = _synthetic_refs()
        fig = mod._section_detected_size_distribution(
            _synthetic_read(), abundance
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_size_distribution_renders_without_reference(self, mod):
        # mass-spectrum precedent: sim-only rendering, never a skip
        fig = mod._section_detected_size_distribution(_synthetic_read(), None)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_ked_mean_energy_renders_and_annotates_chi2(self, mod):
        _, ked = _synthetic_refs()
        fig = mod._section_detected_ked_mean_energy(_synthetic_read(), ked)
        assert isinstance(fig, plt.Figure)
        text = " ".join(
            t.get_text() for ax in fig.axes for t in ax.texts
        ) + " " + " ".join(t.get_text() for t in fig.texts)
        assert "med" in text  # both chi^2 columns annotated
        plt.close(fig)

    def test_ked_mean_energy_skips_without_reference(self, mod):
        with pytest.raises(mod._SectionSkipped):
            mod._section_detected_ked_mean_energy(_synthetic_read(), None)

    def test_ke_anatomy_renders(self, mod):
        _, ked = _synthetic_refs()
        fig = mod._section_detected_ke_anatomy(_synthetic_read(), ked)
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_ke_anatomy_skips_without_reference(self, mod):
        with pytest.raises(mod._SectionSkipped):
            mod._section_detected_ke_anatomy(_synthetic_read(), None)


class TestMedianAnchorConvention:
    def test_ked_section_uses_median_anchor_at_n1(self, mod):
        """The annotated chi2_med must equal the committed scorer's
        median-anchor value (not the mean-legacy one) — the I77 operative
        convention, locked against the module the section imports from."""
        from i2_helium_md.postprocess.tier2_confirmation import (
            score_ke_curve_vs_reference,
        )
        _, ked = _synthetic_refs()
        read = _synthetic_read()
        want_med = score_ke_curve_vs_reference(
            read, ked, n_min=1, min_count=2, include_sim_se=True,
            n1_anchor="median",
        ).chi2_profiled
        fig = mod._section_detected_ked_mean_energy(read, ked)
        text = " ".join(
            t.get_text() for ax in fig.axes for t in ax.texts
        ) + " " + " ".join(t.get_text() for t in fig.texts)
        assert f"{want_med:.1f}" in text
        plt.close(fig)
