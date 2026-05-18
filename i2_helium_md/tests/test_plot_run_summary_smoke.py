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
):
    monkeypatch.setattr(module, "RUN_DIR", run_dir)
    monkeypatch.setattr(module, "OUT_DIR", out_dir)
    monkeypatch.setattr(module, "HEDFT_REF_PATH", hedft_ref)
    monkeypatch.setattr(module, "VMI_REF_HE_PATH", vmi_ref_he)
    monkeypatch.setattr(module, "VMI_REF_GAS_PATH", vmi_ref_gas)
    monkeypatch.setattr(module, "VMI_REF_HE_HIGH_SNR_PATH", vmi_ref_he_high_snr)
    monkeypatch.setattr(module, "PAPER_V2_REFERENCE_DIR", paper_v2_ref_dir)
    monkeypatch.setattr(module, "SHOW_FIGURES", False)


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
    _set_user_settings(
        monkeypatch, module,
        run_dir=EXPERIMENTAL_RUN,
        out_dir=tmp_path,
        vmi_ref_he=VMI_HE,
        vmi_ref_gas=VMI_GAS,
        paper_v2_ref_dir=tmp_path,
    )
    rc = module.main()
    assert rc == 0
    assert any(
        any(
            ax["title"] == "3-D speed vs Abel-inverted VMI radial distribution"
            for ax in figure
        )
        for figure in saved_figures
    )
    assert any(
        any(
            ax["title"] == "Final ion mass spectrum"
            for ax in figure
        )
        for figure in saved_figures
    )
    plt.close("all")
