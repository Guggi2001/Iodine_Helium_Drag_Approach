"""``plot_detection_summary.py`` — the detection-stage figure set
(`TIER2_DETECTION_SUMMARY_SPEC.md` §2.2 / §4, delivered 2026-07-22).

Locks: the fail-loud no-legacy-mode contract, the §4 roster order, the
metadata cover, and the newly-populated n = 0 (suppressed/bare, RQ3)
speed panel rendered from the ``DetectedEnsembleView``. Synthetic
detections only; the committed ihe_ked reference directory is used for
the curve files (same idiom as ``test_ihe_ked``); Agg backend, figures
closed, nothing written to disk by the section builders.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from tests.test_detected_view import make_rich_detection  # noqa: E402

from i2_helium_md.postprocess import (  # noqa: E402
    detected_ensemble_view,
    load_ihe_ked_reference,
)
from i2_helium_md.postprocess.tier2_confirmation import (  # noqa: E402
    read_confirmation_detection,
)
from i2_helium_md.simulation.detection_stage import (  # noqa: E402
    save_detection_result,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    PROJECT_ROOT / "scripts" / "post_processing" / "plot_detection_summary.py"
)
IHE_KED_DIR = PROJECT_ROOT / "data" / "reference" / "ihe_ked"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location(
        "plot_detection_summary_under_test", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _detection():
    # 4 molecules / 8 ions: solvated n=1/2 (frozen), one suppressed pair
    # member (rides at handover n=6 on disk, scores bare), one retained,
    # one time_exhausted at n=1; speeds unit-scale A/ps.
    return make_rich_detection(
        [1.0, 2.0, 1.0, 2.0, 6.0, 9.0, 1.0, 1.0],
        ["frozen", "frozen", "frozen", "frozen",
         "suppressed", "droplet_retained", "time_exhausted", "frozen"],
        vx=[3.0, 2.0, 3.5, 2.1, 8.0, 99.0, 4.0, 3.2],
        vy=[4.0, 1.0, 3.0, 1.5, 6.0, 99.0, 0.0, 3.8],
        vz=[0.0, 0.5, 0.2, 0.3, 0.0, 99.0, 1.0, 0.1],
    )


ROSTER = (
    "detection_metadata",
    "detected_size_distribution",
    "detected_ihe_ked_mean_energy",
    "detected_ke_anatomy",
    "detected_ihe_ked_curves_3d",
    "detected_ihe_ked_curves_2d",
    "detected_mass_resolved_velocities",
    "detected_paper_v2_vmi",
    "detected_paper_v2_polar",
    "detected_paper_cov_radial_distribution",
    "detected_paper_cov_phi_distribution",
    "detected_paper_cov_angular_pair_cov",
    "detected_paper_cov_radial_pair_cov",
    "detected_paper_cov_pair_cov_traces",
)


class TestFailLoud:
    def test_fails_loudly_without_detection_npz(self, mod, tmp_path):
        with pytest.raises(FileNotFoundError, match="detection.npz"):
            mod.main(tmp_path)


class TestRoster:
    def test_section_order_matches_spec(self, mod, tmp_path):
        det = _detection()
        read = read_confirmation_detection(det, label="syn")
        view = detected_ensemble_view(det)
        sections = mod._sections(
            run_dir=tmp_path, cfg=None, det=det, read=read, view=view,
            ked_dir=None, ked_ref=None, abundance=None,
            paper_v2_ref_dir=None, paper_cov_ref_dir=None,
        )
        assert tuple(label for label, _ in sections) == ROSTER

    def test_reference_less_sections_skip(self, mod, tmp_path):
        det = _detection()
        read = read_confirmation_detection(det, label="syn")
        view = detected_ensemble_view(det)
        sections = dict(mod._sections(
            run_dir=tmp_path, cfg=None, det=det, read=read, view=view,
            ked_dir=None, ked_ref=None, abundance=None,
            paper_v2_ref_dir=None, paper_cov_ref_dir=None,
        ))
        for label in ("detected_paper_v2_vmi",
                      "detected_paper_cov_angular_pair_cov",
                      "detected_ihe_ked_curves_3d"):
            with pytest.raises(mod._SectionSkipped):
                sections[label]()


class TestMetadataCover:
    def test_metadata_section_renders_stage_identity(self, mod, tmp_path):
        det = _detection()
        read = read_confirmation_detection(det, label="syn")
        cfg = {"drag_form": "capped_cubic",
               "internal_energy_cooling_tau_ps": 3.2,
               "dissociation_ladder": "rq4graded",
               "drag_coefficients": {"coefficients": {"v_c": 7.25,
                                                      "p_tail": -1.0}}}
        fig = mod._section_detection_metadata(
            tmp_path, cfg, det, read
        )
        assert isinstance(fig, plt.Figure)
        text = " ".join(t.get_text() for ax in fig.axes for t in ax.texts)
        assert "t_detect" in text
        assert "suppressed" in text
        assert "droplet_retained" in text
        assert "capped_cubic" in text
        plt.close(fig)


class TestDetectedSpeedPanels:
    def test_n0_panel_carries_the_suppressed_channel(self, mod):
        # the RQ3 gain: the n = 0 speed panel is populated by the
        # suppressed/bare channel, which the legacy summary never could
        if not IHE_KED_DIR.exists():
            pytest.skip("committed ihe_ked reference not present")
        view = detected_ensemble_view(_detection())
        ked_ref = load_ihe_ked_reference(
            IHE_KED_DIR / "IHe_KED_reference.csv"
        )
        fig = mod._section_ihe_ked_curves(
            view, IHE_KED_DIR, ked_ref, "3d",
            stage_note=mod._DETECTED_NOTE,
        )
        n0_ax = fig.axes[0]
        labels = [t.get_text() for t in n0_ax.get_legend().get_texts()]
        assert any("simulation (N=1)" in lab for lab in labels), labels
        note = " ".join(t.get_text() for t in fig.texts)
        assert "suppressed" in note
        plt.close(fig)

    def test_mass_resolved_renders_from_view(self, mod):
        view = detected_ensemble_view(_detection())
        fig = mod._section_mass_resolved(
            view, stage_note=mod._DETECTED_NOTE
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)


class TestArtifactRoundTrip:
    def test_load_path_via_saved_detection(self, mod, tmp_path):
        save_detection_result(_detection(), tmp_path / "detection.npz")
        (tmp_path / "cfg.json").write_text("{}", encoding="utf-8")
        det, read, view, cfg = mod._load_inputs(tmp_path)
        assert view.num_molecules == 4
        assert read.num_ions == 8
        assert np.isnan(view.mass_final_kg[5])  # retained forced out
        assert cfg == {}
