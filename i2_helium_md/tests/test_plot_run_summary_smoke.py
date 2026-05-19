"""Smoke test for scripts/post_processing/plot_run_summary.py.

Imports the script as a module, monkey-patches ``plt.show`` and
``Figure.savefig``, overrides the module-level ``USER SETTINGS``
constants for the run under test, and calls ``module.main()``. Skipped
when no usable run is present.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

import pytest  # noqa: E402
import numpy as np  # noqa: E402
from scipy.io import savemat  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "post_processing" / "plot_run_summary.py"

EXPERIMENTAL_RUN = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
HEDFT_RUN_9A = PROJECT_ROOT / "data" / "runs" / "9A_hedft_comparison"
HEDFT_RUN_18A = PROJECT_ROOT / "data" / "runs" / "18A_hedft_comparison"

REF_9A = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
REF_18A = PROJECT_ROOT / "data" / "reference" / "18A_All_Data.csv"
VMI_HE = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_he.csv"
VMI_GAS = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_gas.csv"


def _import_script():
    spec = importlib.util.spec_from_file_location(
        "plot_run_summary_under_test", SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["plot_run_summary_under_test"] = module
    spec.loader.exec_module(module)
    return module


def _patch_io(monkeypatch):
    saved_figures = []

    def capture_pdf_figure(self, fig, *args, **kwargs):
        saved_figures.append([
            {"title": ax.get_title(), "xlim": ax.get_xlim()}
            for ax in fig.axes
        ])

    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    monkeypatch.setattr(plt.Figure, "savefig", lambda *a, **k: None)
    monkeypatch.setattr(PdfPages, "savefig", capture_pdf_figure)
    return saved_figures


def _set_user_settings(
    monkeypatch,
    module,
    *,
    run_dir,
    out_dir,
    hedft_ref=None,
    vmi_ref_he=None,
    vmi_ref_gas=None,
    vmi_ref_he_high_snr=None,
    paper_v2_ref_dir=None,
    paper_cov_ref_dir=None,
):
    monkeypatch.setattr(module, "RUN_DIR", run_dir)
    monkeypatch.setattr(module, "OUT_DIR", out_dir)
    monkeypatch.setattr(module, "HEDFT_REF_PATH", hedft_ref)
    monkeypatch.setattr(module, "VMI_REF_HE_PATH", vmi_ref_he)
    monkeypatch.setattr(module, "VMI_REF_GAS_PATH", vmi_ref_gas)
    monkeypatch.setattr(module, "VMI_REF_HE_HIGH_SNR_PATH", vmi_ref_he_high_snr)
    monkeypatch.setattr(module, "PAPER_V2_REFERENCE_DIR", paper_v2_ref_dir)
    monkeypatch.setattr(module, "PAPER_COV_REFERENCE_DIR", paper_cov_ref_dir)
    monkeypatch.setattr(module, "SHOW_FIGURES", False)


def _write_synthetic_paper_cov_reference(reference_dir: Path) -> Path:
    reference_dir.mkdir(parents=True, exist_ok=True)
    n_theta = 16
    n_v = 20
    cov_angular = np.zeros((n_theta, n_theta), dtype=float)
    cov_radial = np.zeros((n_v, n_v), dtype=float)
    cov_angular[3, 11] = cov_angular[11, 3] = 1.0
    cov_radial[5, 7] = cov_radial[7, 5] = 1.0
    path = reference_dir / "iplus_he_covariance.mat"
    savemat(
        path,
        {
            "cov_angular": cov_angular,
            "cov_radial": cov_radial,
            "theta_centers_rad": np.linspace(-np.pi, np.pi, n_theta, endpoint=False),
            "velocity_centers_mps": np.linspace(0.0, 2500.0, n_v),
        },
    )
    return path


@pytest.mark.skipif(
    not (HEDFT_RUN_9A / "ion.npz").exists() or not REF_9A.exists(),
    reason="9A HeDFT run or reference CSV not present",
)
def test_run_summary_on_hedft_9a(monkeypatch, tmp_path):
    plt.close("all")
    module = _import_script()
    _patch_io(monkeypatch)
    _set_user_settings(
        monkeypatch, module,
        run_dir=HEDFT_RUN_9A,
        out_dir=tmp_path,
        hedft_ref=REF_9A,
        paper_v2_ref_dir=tmp_path,
    )
    rc = module.main()
    assert rc == 0
    plt.close("all")


@pytest.mark.skipif(
    not (EXPERIMENTAL_RUN / "ion.npz").exists()
    or not VMI_HE.exists() or not VMI_GAS.exists(),
    reason="experimental run or VMI reference CSVs not present",
)
def test_run_summary_on_experimental(monkeypatch, tmp_path):
    plt.close("all")
    module = _import_script()
    saved_figures = _patch_io(monkeypatch)
    paper_cov_ref_dir = tmp_path / "paper_cov"
    _write_synthetic_paper_cov_reference(paper_cov_ref_dir)
    _set_user_settings(
        monkeypatch, module,
        run_dir=EXPERIMENTAL_RUN,
        out_dir=tmp_path,
        vmi_ref_he=VMI_HE,
        vmi_ref_gas=VMI_GAS,
        paper_v2_ref_dir=tmp_path,
        paper_cov_ref_dir=paper_cov_ref_dir,
    )
    rc = module.main()
    assert rc == 0
    titles = [ax["title"] for figure in saved_figures for ax in figure]
    assert any(
        title == "3-D speed vs Abel-inverted VMI radial distribution"
        for title in titles
    )
    assert "(c) velocity distribution comparison" in titles
    assert "(f) phi angular distribution" in titles
    assert "(d) experimental angular pair covariance" in titles
    assert "(e) experimental radial pair covariance" in titles
    assert "(g) angular pair-cov trace" in titles
    assert "(h) radial pair-cov trace" in titles
    assert "2-D detector-plane speed vs raw VMI radial profile" not in titles
    assert "paper v2 azimuthal distribution" not in titles
    assert not any(title.startswith("Angular pair covariance (") for title in titles)
    assert any(
        title == "Final ion mass spectrum"
        for title in titles
    )
    plt.close("all")
