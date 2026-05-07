"""Smoke test for scripts/plot_hedft_comparison.py.

Loads the script as a module via importlib (so the ``USER SETTINGS`` at
the top get evaluated normally), monkeypatches ``plt.show``, runs
``main()`` once with the VMI tile and once without, and asserts that
two figures get created in each case.

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
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "plot_hedft_comparison.py"
ION_NPZ = PROJECT_ROOT / "data" / "runs" / "single_pulse_N_2000" / "ion.npz"
HEDFT = PROJECT_ROOT / "data" / "reference" / "9A_All_Data.csv"
VMI_HE = PROJECT_ROOT / "data" / "reference" / "vmi_iplus_he.csv"
VMI_GAS = PROJECT_ROOT / "data" / "reference" / "vmi_iplus_gas.csv"


pytestmark = pytest.mark.skipif(
    not (
        ION_NPZ.exists()
        and HEDFT.exists()
        and VMI_HE.exists()
        and VMI_GAS.exists()
    ),
    reason="real run or reference data not present",
)


def _import_script() -> object:
    """Load the script as a module without running ``main()``."""
    spec = importlib.util.spec_from_file_location(
        "plot_hedft_comparison_under_test", SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["plot_hedft_comparison_under_test"] = module
    spec.loader.exec_module(module)
    return module


def test_main_with_vmi_tile_produces_two_figures(monkeypatch):
    plt.close("all")
    module = _import_script()

    monkeypatch.setattr(module, "SHOW_VMI_TILE", True)
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)

    rc = module.main()
    assert rc == 0
    assert len(plt.get_fignums()) == 2
    plt.close("all")


def test_main_without_vmi_tile_produces_two_figures(monkeypatch):
    plt.close("all")
    module = _import_script()

    monkeypatch.setattr(module, "SHOW_VMI_TILE", False)
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)

    rc = module.main()
    assert rc == 0
    assert len(plt.get_fignums()) == 2
    plt.close("all")
