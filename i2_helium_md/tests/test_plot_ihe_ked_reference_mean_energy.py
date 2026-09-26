"""Non-interactive smoke tests for the experiment-only ihe_ked mean-KE figure.

The script is the *experimental reference chapter* figure: the
detection-summary mean-KE panel with the simulation overlay removed. These
tests check that it renders from the frozen reference alone (no run dir),
honors the n floor, and never draws a simulation series.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT = (
    PROJECT_ROOT
    / "scripts"
    / "post_processing"
    / "plot_ihe_ked_reference_mean_energy.py"
)
_REFERENCE_CSV = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked" / "IHe_KED_reference.csv"
)

from i2_helium_md.postprocess import load_ihe_ked_reference


def _load_script_module():
    sys.path.insert(0, str(_SCRIPT.parent))  # its `summary_sections` import
    try:
        spec = importlib.util.spec_from_file_location(
            "plot_ihe_ked_reference_mean_energy", _SCRIPT
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(_SCRIPT.parent))
    return module


@pytest.fixture(scope="module")
def script_module():
    return _load_script_module()


@pytest.fixture(scope="module")
def ked_ref():
    if not _REFERENCE_CSV.exists():
        pytest.skip(f"reference CSV missing: {_REFERENCE_CSV}")
    return load_ihe_ked_reference(_REFERENCE_CSV)


def test_renders_two_panels_without_a_run_directory(script_module, ked_ref):
    fig = script_module.build_reference_mean_energy_figure(ked_ref)
    try:
        assert len(fig.axes) == 2
        # Left panel carries the experiment series and both correlated
        # bands -- and no simulation series at all.
        labels = fig.axes[0].get_legend_handles_labels()[1]
        assert any("experiment" in lbl for lbl in labels)
        assert sum("band (correlated)" in lbl for lbl in labels) == 2
        assert not any("simulation" in lbl for lbl in labels)
    finally:
        plt.close(fig)


def test_n_min_zero_adds_the_bare_fragment(script_module, ked_ref):
    """n = 0 (bare I+) is measured data; the default floor hides it."""
    x_default = _scatter_x(script_module, ked_ref, n_min=1)
    x_with_bare = _scatter_x(script_module, ked_ref, n_min=0)
    assert x_default.min() == pytest.approx(1.0)
    assert x_with_bare.min() == pytest.approx(0.0)
    assert x_with_bare.size == x_default.size + 1


def test_rejects_an_empty_fragment_selection(script_module, ked_ref):
    with pytest.raises(ValueError, match="no reference rows"):
        script_module.build_reference_mean_energy_figure(
            ked_ref, n_min=int(ked_ref.n.max()) + 1
        )


def _scatter_x(script_module, ked_ref, *, n_min: int) -> np.ndarray:
    """x values of every errorbar-plotted point across both panels."""
    fig = script_module.build_reference_mean_energy_figure(
        ked_ref, n_min=n_min
    )
    try:
        xs = [
            np.asarray(container.lines[0].get_xdata())
            for ax in fig.axes
            for container in ax.containers
        ]
        return np.concatenate(xs) if xs else np.empty(0)
    finally:
        plt.close(fig)
