"""Smoke tests for the post-processing plotting scripts.

Loads each script as a module via importlib, monkeypatches ``plt.show``,
runs ``main()``, and asserts that the expected figure count is created.

Gated on the real run + reference data being present so CI without
those artifacts skips cleanly.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import pytest  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = PROJECT_ROOT / "scripts" / "post_processing"
HEDFT_SCRIPT = SCRIPT_DIR / "plot_hedft_comparison.py"
EXPERIMENTAL_SCRIPT = SCRIPT_DIR / "plot_experimental_comparison.py"

HEDFT_ION_NPZ = (
    PROJECT_ROOT
    / "data"
    / "runs"
    / "9A_hedft_comparison"
    / "ion.npz"
)
EXPERIMENTAL_ION_NPZ = (
    PROJECT_ROOT
    / "data"
    / "runs"
    / "single_pulse_droplet_log_droplet"
    / "ion.npz"
)
HEDFT = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
VMI_HE = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_he.csv"
VMI_GAS = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_gas.csv"


pytestmark = pytest.mark.skipif(
    not (
        HEDFT_ION_NPZ.exists()
        and EXPERIMENTAL_ION_NPZ.exists()
        and HEDFT.exists()
        and VMI_HE.exists()
        and VMI_GAS.exists()
    ),
    reason="real run or reference data not present",
)


def _import_script(path: Path, module_name: str) -> object:
    """Load a script as a module without running ``main()``."""
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def test_hedft_comparison_main_produces_two_figures(monkeypatch):
    plt.close("all")
    module = _import_script(
        HEDFT_SCRIPT,
        "plot_hedft_comparison_under_test",
    )

    monkeypatch.setattr(plt, "show", lambda *a, **k: None)

    rc = module.main()
    assert rc == 0
    assert len(plt.get_fignums()) == 2
    plt.close("all")


def test_experimental_comparison_main_produces_one_figure(monkeypatch):
    plt.close("all")
    module = _import_script(
        EXPERIMENTAL_SCRIPT,
        "plot_experimental_comparison_under_test",
    )

    monkeypatch.setattr(plt, "show", lambda *a, **k: None)

    rc = module.main()
    assert rc == 0
    assert len(plt.get_fignums()) == 1
    plt.close("all")
