"""Smoke test for scripts/post_processing/plot_run_summary.py.

Imports the script as a module, monkey-patches ``plt.show`` and
``Figure.savefig``, and runs ``main(argv=[...])`` on each available
real run directory. Skipped when no usable run is present.
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
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    monkeypatch.setattr(plt.Figure, "savefig", lambda *a, **k: None)
    monkeypatch.setattr(PdfPages, "savefig", lambda self, *a, **k: None)


@pytest.mark.skipif(
    not (HEDFT_RUN_9A / "ion.npz").exists() or not REF_9A.exists(),
    reason="9A HeDFT run or reference CSV not present",
)
def test_run_summary_on_hedft_9a(monkeypatch, tmp_path):
    plt.close("all")
    module = _import_script()
    _patch_io(monkeypatch)
    rc = module.main(argv=[
        str(HEDFT_RUN_9A),
        "--hedft-ref", str(REF_9A),
        "--out-dir", str(tmp_path),
        "--no-show",
    ])
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
    _patch_io(monkeypatch)
    rc = module.main(argv=[
        str(EXPERIMENTAL_RUN),
        "--vmi-ref-he", str(VMI_HE),
        "--vmi-ref-gas", str(VMI_GAS),
        "--out-dir", str(tmp_path),
        "--no-show",
    ])
    assert rc == 0
    plt.close("all")
