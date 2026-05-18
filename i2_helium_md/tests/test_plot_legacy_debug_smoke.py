"""Smoke tests for the legacy-debug post-processing plot scripts.

Each test imports the script as a module, monkeypatches plt.show and
plt.savefig (so figures are not written during test runs), invokes
``main()``, and asserts the expected figure count.

Gated on the realistic experimental-condition run being present so CI
without that artifact skips cleanly. The temperature-diagnostic test
additionally requires a v5 ion checkpoint with at least one
non-NaN row.
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

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = PROJECT_ROOT / "scripts" / "post_processing"

NEUTRAL_SCRIPT = SCRIPT_DIR / "plot_neutral_energy_balance.py"
ION_ENERGY_SCRIPT = SCRIPT_DIR / "plot_ion_energy_balance.py"
ION_TEMP_SCRIPT = SCRIPT_DIR / "plot_ion_temperature_diagnostic.py"
PAPER_SCRIPT = SCRIPT_DIR / "plot_paper_v3.py"
PAPER_V2_SCRIPT = SCRIPT_DIR / "plot_paper_v2.py"
PAPER_V4_SCRIPT = SCRIPT_DIR / "plot_paper_v4.py"

EXPERIMENTAL_RUN = (
    PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
)
EXPERIMENTAL_NEUTRAL = EXPERIMENTAL_RUN / "neutral.npz"
EXPERIMENTAL_ION = EXPERIMENTAL_RUN / "ion.npz"
VMI_HE = PROJECT_ROOT / "data" / "reference" /  "vmi_summary" / "vmi_iplus_he.csv"
VMI_GAS = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_gas.csv"


pytestmark = pytest.mark.skipif(
    not (
        EXPERIMENTAL_NEUTRAL.exists()
        and EXPERIMENTAL_ION.exists()
        and VMI_HE.exists()
        and VMI_GAS.exists()
    ),
    reason="experimental-condition run or VMI references not present",
)


def _import_script(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _patch_io(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    monkeypatch.setattr(plt.Figure, "savefig", lambda *a, **k: None)


def _override_run_dir(module, attr: str = "RUN_DIR"):
    """Repoint a script's RUN_DIR at the experimental run directory."""
    setattr(module, attr, EXPERIMENTAL_RUN)


def test_neutral_energy_balance_smoke(monkeypatch):
    plt.close("all")
    module = _import_script(NEUTRAL_SCRIPT, "plot_neutral_energy_balance_under_test")
    _patch_io(monkeypatch)
    _override_run_dir(module)

    rc = module.main()
    assert rc == 0
    assert len(plt.get_fignums()) == 1
    plt.close("all")


def test_ion_energy_balance_smoke(monkeypatch):
    plt.close("all")
    module = _import_script(ION_ENERGY_SCRIPT, "plot_ion_energy_balance_under_test")
    _patch_io(monkeypatch)
    _override_run_dir(module)

    rc = module.main()
    assert rc == 0
    assert len(plt.get_fignums()) == 1
    plt.close("all")


def _ion_has_diagnostic_rows(path: Path) -> bool:
    try:
        with np.load(path, allow_pickle=False) as z:
            if "temperature_diagnostic" not in z.files:
                return False
            td = z["temperature_diagnostic"]
            return bool(np.isfinite(td[:, 0]).any())
    except Exception:
        return False


def test_ion_temperature_diagnostic_smoke(monkeypatch):
    if not _ion_has_diagnostic_rows(EXPERIMENTAL_ION):
        pytest.skip(
            "ion.npz lacks usable temperature_diagnostic rows "
            "(likely a pre-v5 file or a run with no collisions)"
        )
    plt.close("all")
    module = _import_script(ION_TEMP_SCRIPT, "plot_ion_temperature_under_test")
    _patch_io(monkeypatch)
    _override_run_dir(module)

    rc = module.main()
    assert rc == 0
    assert len(plt.get_fignums()) == 1
    plt.close("all")


def test_paper_figure_smoke(monkeypatch):
    plt.close("all")
    module = _import_script(PAPER_SCRIPT, "plot_paper_v3_under_test")
    _patch_io(monkeypatch)
    _override_run_dir(module)

    rc = module.main(["--no-show"])
    assert rc == 0
    assert len(plt.get_fignums()) == 2
    legend_labels = [
        text.get_text()
        for fig_num in plt.get_fignums()
        for ax in plt.figure(fig_num).axes
        if ax.get_legend() is not None
        for text in ax.get_legend().get_texts()
    ]
    assert any("(296:297)" in label for label in legend_labels)
    plt.close("all")


def test_paper_v2_figure_smoke(monkeypatch, tmp_path):
    plt.close("all")
    module = _import_script(PAPER_V2_SCRIPT, "plot_paper_v2_under_test")
    _patch_io(monkeypatch)
    _override_run_dir(module)

    rc = module.main(["--no-show", "--reference-dir", str(tmp_path)])
    assert rc == 0
    assert len(plt.get_fignums()) == 5
    plt.close("all")


def test_paper_v4_figure_smoke(monkeypatch, tmp_path):
    plt.close("all")
    module = _import_script(PAPER_V4_SCRIPT, "plot_paper_v4_under_test")
    _patch_io(monkeypatch)
    _override_run_dir(module)

    rc = module.main(["--no-show", "--reference-dir", str(tmp_path)])
    assert rc == 0
    assert len(plt.get_fignums()) == 3
    plt.close("all")
