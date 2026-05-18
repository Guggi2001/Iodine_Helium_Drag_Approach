"""Smoke tests for plot_paper_cov.py.

Two cases:

1. With a synthetic experimental covariance reference present: expect
   six figures (VMI, radial distribution, phi distribution, angular cov,
   radial cov, pair-cov traces).
2. With no experimental reference: expect three figures (VMI, radial
   distribution, and phi distribution; the three cov-derived figures
   are skipped).

The phi distribution figure always renders -- if the optional 2-D VMI
image reference is missing it draws only the simulated curve.

Gated on the realistic experimental-condition run being present so CI
without that artifact skips cleanly.
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
from scipy.io import savemat  # noqa: E402


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = PROJECT_ROOT / "scripts" / "post_processing"
PAPER_COV_SCRIPT = SCRIPT_DIR / "plot_paper_cov.py"

EXPERIMENTAL_RUN = PROJECT_ROOT / "data" / "runs" / "single_pulse_droplet"
EXPERIMENTAL_ION = EXPERIMENTAL_RUN / "ion.npz"


pytestmark = pytest.mark.skipif(
    not EXPERIMENTAL_ION.exists(),
    reason="experimental-condition run not present",
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
    setattr(module, attr, EXPERIMENTAL_RUN)


def _write_synthetic_cov_reference(reference_dir: Path) -> Path:
    """Tiny synthetic covariance reference that satisfies the loader contract."""
    reference_dir.mkdir(parents=True, exist_ok=True)
    n_theta = 16
    n_v = 20
    cov_angular = np.zeros((n_theta, n_theta), dtype=float)
    cov_radial = np.zeros((n_v, n_v), dtype=float)
    cov_angular[3, 11] = cov_angular[11, 3] = 1.0
    cov_radial[5, 7] = cov_radial[7, 5] = 1.0
    theta = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)
    velocity_mps = np.linspace(0.0, 2500.0, n_v)
    out_path = reference_dir / "iplus_he_covariance.mat"
    savemat(
        out_path,
        {
            "cov_angular": cov_angular,
            "cov_radial": cov_radial,
            "theta_centers_rad": theta,
            "velocity_centers_mps": velocity_mps,
        },
    )
    return out_path


def test_paper_cov_with_reference_produces_six_figures(monkeypatch, tmp_path):
    plt.close("all")
    module = _import_script(PAPER_COV_SCRIPT, "plot_paper_cov_with_ref_under_test")
    _patch_io(monkeypatch)
    _override_run_dir(module)

    cov_dir = tmp_path / "paper_cov"
    _write_synthetic_cov_reference(cov_dir)
    paper_v2_dir = tmp_path / "paper_v2_empty"
    paper_v2_dir.mkdir()

    rc = module.main([
        "--no-show",
        "--reference-dir", str(cov_dir),
        "--paper-v2-reference-dir", str(paper_v2_dir),
    ])
    assert rc == 0
    # VMI + radial + phi + angular cov + radial cov + cov traces.
    assert len(plt.get_fignums()) == 6
    plt.close("all")


def test_paper_cov_without_reference_skips_cov_figures(monkeypatch, tmp_path):
    plt.close("all")
    module = _import_script(PAPER_COV_SCRIPT, "plot_paper_cov_no_ref_under_test")
    _patch_io(monkeypatch)
    _override_run_dir(module)

    empty_cov_dir = tmp_path / "paper_cov_empty"
    empty_cov_dir.mkdir()
    paper_v2_dir = tmp_path / "paper_v2_empty"
    paper_v2_dir.mkdir()

    rc = module.main([
        "--no-show",
        "--reference-dir", str(empty_cov_dir),
        "--paper-v2-reference-dir", str(paper_v2_dir),
    ])
    assert rc == 0
    # VMI comparison + radial distribution + phi distribution (sim-only).
    # The two heatmap cov figures and the cov-trace figure are skipped.
    assert len(plt.get_fignums()) == 3
    plt.close("all")
