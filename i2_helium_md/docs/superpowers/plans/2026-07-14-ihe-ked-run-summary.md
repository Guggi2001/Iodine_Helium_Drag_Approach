# IHe KED Run-Summary Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire the frozen `data/reference/ihe_ked/` experimental reference (per-fragment ⟨E⟩ table n = 0…17 + trusted 2-D/3-D curves n = 0…4) into `scripts/post_processing/plot_run_summary.py`, add a mass-spectrum vs. experimental-abundance side-by-side, and retire the old `vmi_summary` overlay section.

**Architecture:** New package module `i2_helium_md/postprocess/ihe_ked.py` (loaders + sim-side ⟨E⟩/gate helpers) built on existing primitives (`csv_contract.validate_columns`, `physics.shell_schedule.complex_mass_amu`, a mass-gate helper extracted from `velocity_distribution.py`). The plotting script stays a thin wrapper. Spec: `docs/superpowers/specs/2026-07-14-ihe-ked-run-summary-design.md`.

**Tech Stack:** Python 3.14, numpy, pandas, matplotlib, pytest.

## Global Constraints

- Test command (Python is NOT on PATH):
  `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q <file>`
- Work on branch `drag_implementation`. The working tree has unrelated modified files — `git add` ONLY the files named in each task's commit step, never `git add -A`.
- Frozen dataclasses for all return types; loud `FileNotFoundError`/`ValueError` validation (CLAUDE.md rule 4).
- Units in names: `_eV`, `_mps`, `_Aps`, `_amu` (CLAUDE.md rule 3).
- No duplicate physics: mass model comes from `i2_helium_md.physics.shell_schedule.complex_mass_amu` (m(n) = 126.90 + 4.0026·n); constants `EV`, `U` from `i2_helium_md.physics.constants`.
- Tests must not generate figures on disk or production-sized checkpoints; figure smoke tests use the Agg backend and `plt.close`.
- Do NOT touch: physics modules, checkpoint schema, `velocity_distribution.py` beyond the exact Task-1 change, `vmi_summary` data, `load_vmi_reference`.
- Reference data facts used by tests (verified against the real files): `IHe_KED_reference.csv` has 18 rows, n = 0…17, `meanKE_eV` at n=0 is `3.70569`; curve CSVs are `IHe_KED_curves_n{0..4}.csv` with header `E_eV,v_mps,signal_2d_Pv,signal_2d_PE,signal_3d_Pv,signal_3d_PE`; NaN appears in 3-D columns only (below 0.4 eV for n=0, below 0.01 eV otherwise). Gold points (dominantError == `calib`) are n = 0, 3, 4, 5, 6, 7, 8, 9, 10, 12.

---

### Task 1: Mass-gate helper + projected-speed option in `velocity_distribution.py`

Extract the mass-gate selection (currently inline in `compute_final_velocity_histogram`) into a public helper so `ihe_ked.py` can reuse it (DRY), and add a `projected` flag so the same histogram function serves the 2-D detector-projection overlay (in-plane speed √(vx²+vy²), the `paper_v2_velocity_curve` convention).

**Files:**
- Modify: `i2_helium_md/postprocess/velocity_distribution.py` (function `compute_final_velocity_histogram`, lines ~145–250)
- Modify: `i2_helium_md/postprocess/__init__.py` (export `select_final_mass_gate`)
- Test: `tests/test_velocity_distribution.py`

**Interfaces:**
- Consumes: `IonCheckpoint` fields `mass_final_kg (2N,)`, `b_ion_outside (N,)`, `velocities_final_x/y/z (2N,)`.
- Produces: `select_final_mass_gate(ion, *, mass_amu: float, mass_tolerance_amu: float = 0.5, require_outside: bool = True) -> np.ndarray` (bool, shape `(2N,)`); `compute_final_velocity_histogram(..., projected: bool = False)`. Tasks 4–5 rely on both.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_velocity_distribution.py` (the `_make_ion` helper already exists there):

```python
from i2_helium_md.postprocess.velocity_distribution import select_final_mass_gate


class TestSelectFinalMassGate:
    def test_mass_and_outside_mask(self):
        ion = _make_ion(
            num_molecules=2,
            final_speeds_per_atom=np.array([1.0, 2.0, 3.0, 4.0]),
            masses_amu_per_atom=np.array([131.0, 131.0, 127.0, 135.0]),
            b_outside=np.array([True, False]),
        )
        mask = select_final_mass_gate(ion, mass_amu=131.0)
        # atoms 0,1 belong to molecule 0 (outside); atoms 2,3 to molecule 1.
        np.testing.assert_array_equal(mask, [True, True, False, False])

    def test_require_outside_false_keeps_inside_atoms(self):
        ion = _make_ion(
            num_molecules=2,
            final_speeds_per_atom=np.array([1.0, 2.0, 3.0, 4.0]),
            masses_amu_per_atom=np.array([131.0, 127.0, 127.0, 131.0]),
            b_outside=np.array([True, False]),
        )
        mask = select_final_mass_gate(
            ion, mass_amu=131.0, require_outside=False,
        )
        np.testing.assert_array_equal(mask, [True, False, False, True])

    def test_reference_mass_model_inside_default_tolerance(self):
        # The ihe_ked mass model m(1) = 130.9026 amu must gate sim atoms
        # whose rounded final mass is 131 amu (|131 - 130.9026| < 0.5).
        ion = _make_ion(
            num_molecules=1,
            final_speeds_per_atom=np.array([5.0, 6.0]),
            masses_amu_per_atom=np.array([131.0, 127.0]),
        )
        mask = select_final_mass_gate(ion, mass_amu=130.9026)
        np.testing.assert_array_equal(mask, [True, False])


class TestProjectedHistogram:
    def test_projected_uses_only_in_plane_components(self):
        # One atom with vx=3, vy=4, vz=12: |v| = 13, in-plane = 5.
        n = 1
        ion = _make_ion(
            num_molecules=n,
            final_speeds_per_atom=np.array([0.0, 0.0]),
            masses_amu_per_atom=np.array([131.0, 127.0]),
        )
        object.__setattr__(ion, "velocities_final_x", np.array([3.0, 0.0]))
        object.__setattr__(ion, "velocities_final_y", np.array([4.0, 0.0]))
        object.__setattr__(ion, "velocities_final_z", np.array([12.0, 0.0]))

        h3 = compute_final_velocity_histogram(
            ion, mass_amu=131.0, num_bins=28, v_max_Aps=28.0,
        )
        h2 = compute_final_velocity_histogram(
            ion, mass_amu=131.0, num_bins=28, v_max_Aps=28.0,
            projected=True,
        )
        assert h3.counts[13] == 1  # bin [13, 14)
        assert h2.counts[5] == 1   # bin [5, 6)
        assert h2.num_atoms_used == 1
```

Note: if `IonCheckpoint` is not a frozen dataclass, replace the three `object.__setattr__` lines with plain attribute assignment — check `i2_helium_md/simulation/checkpoint.py` first and use whichever the class requires.

- [ ] **Step 2: Run tests to verify they fail**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_velocity_distribution.py -k "SelectFinalMassGate or Projected"`
Expected: FAIL / ERROR with `ImportError: cannot import name 'select_final_mass_gate'`.

- [ ] **Step 3: Implement**

In `i2_helium_md/postprocess/velocity_distribution.py`, add above `compute_final_velocity_histogram`:

```python
def select_final_mass_gate(
    ion: IonCheckpoint,
    *,
    mass_amu: float,
    mass_tolerance_amu: float = 0.5,
    require_outside: bool = True,
) -> np.ndarray:
    """Per-atom boolean mask selecting the final-mass gate.

    The single mass-gate convention shared by every mass-gated diagnostic
    (final-velocity histograms, ihe_ked fragment moments, abundance
    comparison): ``round(mass_final_kg / U)`` within ``mass_tolerance_amu``
    of ``mass_amu``, optionally AND-ed with the per-molecule
    ``b_ion_outside`` flag broadcast to both atoms.

    Parameters
    ----------
    ion
        Ion-stage checkpoint; reads ``mass_final_kg`` (shape ``(2N,)``)
        and ``b_ion_outside`` (shape ``(N,)``).
    mass_amu
        Target atomic mass in amu (may be non-integer, e.g. the ihe_ked
        mass model m(n) = 126.90 + 4.0026 n).
    mass_tolerance_amu
        Half-width of the acceptance window in amu.
    require_outside
        If ``True`` (default), both atoms of a molecule pass only when its
        ``b_ion_outside`` flag is set.

    Returns
    -------
    np.ndarray
        Boolean mask, shape ``(2N,)``.
    """
    mass_amu_per_atom = np.round(np.asarray(ion.mass_final_kg) / U_KG)
    mask = (
        np.abs(mass_amu_per_atom - float(mass_amu))
        <= float(mass_tolerance_amu)
    )
    if require_outside:
        outside_per_atom = np.concatenate(
            [ion.b_ion_outside, ion.b_ion_outside]
        ).astype(bool)
        mask = mask & outside_per_atom
    return mask
```

In `compute_final_velocity_histogram`:
1. Add parameter `projected: bool = False` (after `require_outside`).
2. Replace the inline selection block (`mass_amu_per_atom = ...` through `select = mass_mask`) with:

```python
    select = select_final_mass_gate(
        ion,
        mass_amu=mass_amu,
        mass_tolerance_amu=mass_tolerance_amu,
        require_outside=require_outside,
    )
```

3. Replace the `speed_final = ...` expression with:

```python
    vx = np.asarray(ion.velocities_final_x)[select]
    vy = np.asarray(ion.velocities_final_y)[select]
    vz = np.asarray(ion.velocities_final_z)[select]
    if projected:
        # In-plane detector projection sqrt(vx^2 + vy^2) -- the same
        # plane convention as paper_v2_velocity_curve.
        speed_final = np.sqrt(vx * vx + vy * vy)
    else:
        speed_final = np.sqrt(vx * vx + vy * vy + vz * vz)
```

4. Document `projected` in the docstring: "If ``True``, bin the in-plane projected speed √(vx²+vy²) instead of |v| — the 2-D detector-projection observable (`signal_2d_Pv` in the ihe_ked curves)."

In `i2_helium_md/postprocess/__init__.py`: add `select_final_mass_gate` to the `velocity_distribution` import block and to `__all__` (next to `compute_final_velocity_histogram`).

- [ ] **Step 4: Run the full module test to verify pass + no regression**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_velocity_distribution.py`
Expected: all PASS (old behavior unchanged: default `projected=False`).

- [ ] **Step 5: Commit**

```bash
git add i2_helium_md/postprocess/velocity_distribution.py i2_helium_md/postprocess/__init__.py tests/test_velocity_distribution.py
git commit -m "feat(postprocess): extract select_final_mass_gate, add projected-speed histograms"
```

---

### Task 2: `ihe_ked.py` — reference-table loader

**Files:**
- Create: `i2_helium_md/postprocess/ihe_ked.py`
- Modify: `i2_helium_md/postprocess/__init__.py`
- Test: `tests/test_ihe_ked.py` (create)

**Interfaces:**
- Consumes: `validate_columns` from `.csv_contract`.
- Produces: `IHeKedReference` (frozen dataclass of column arrays + `point_err_eV` / `gold_mask` properties) and `load_ihe_ked_reference(path: str | Path) -> IHeKedReference`. Task 5 relies on attribute names exactly as written here.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_ihe_ked.py`:

```python
"""Tests for i2_helium_md/postprocess/ihe_ked.py.

The loader tests run against the real frozen reference under
``data/reference/ihe_ked/`` (small CSVs, committed) plus synthetic
corrupted files in tmp_path for the failure paths.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from i2_helium_md.postprocess.ihe_ked import (
    IHeKedReference,
    load_ihe_ked_reference,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IHE_KED_DIR = PROJECT_ROOT / "data" / "reference" / "ihe_ked"
REFERENCE_CSV = IHE_KED_DIR / "IHe_KED_reference.csv"


class TestLoadIHeKedReference:
    def test_real_file_contract(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        assert isinstance(ref, IHeKedReference)
        np.testing.assert_array_equal(ref.n, np.arange(18))
        # Frozen reference value (README: I+ scale anchor 3.706 eV).
        assert ref.mean_KE_eV[0] == pytest.approx(3.70569, abs=1e-5)
        # Mass model column matches m(n) = 126.90 + 4.0026 n.
        np.testing.assert_allclose(
            ref.mass_center_u, 126.90 + 4.0026 * ref.n, atol=1e-4,
        )
        # noiseLimited is 0 for all 18 fragments (README).
        assert not ref.noise_limited.any()
        assert ref.source_path == REFERENCE_CSV.resolve()

    def test_gold_mask_is_calib_limited_set(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        expected_gold = {0, 3, 4, 5, 6, 7, 8, 9, 10, 12}
        assert set(ref.n[ref.gold_mask].tolist()) == expected_gold

    def test_point_err_is_stat_sys_quadrature(self):
        ref = load_ihe_ked_reference(REFERENCE_CSV)
        np.testing.assert_allclose(
            ref.point_err_eV,
            np.sqrt(
                ref.stat_err_mean_KE_eV ** 2 + ref.sys_err_mean_KE_eV ** 2
            ),
            rtol=1e-12,  # tight: pure arithmetic identity
        )

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_ihe_ked_reference(tmp_path / "nope.csv")

    def test_missing_column_raises(self, tmp_path):
        p = tmp_path / "bad.csv"
        p.write_text("n,label\n0,I+He_0\n", encoding="ascii")
        with pytest.raises(ValueError):
            load_ihe_ked_reference(p)

    def test_non_contiguous_n_raises(self, tmp_path):
        header = (
            "n,label,massCenter_u_per_e,N_counts,N_eff,meanKE_eV,"
            "modeKE_eV,medianKE_eV,sigmaKE_eV,statErr_meanKE_eV,"
            "sysErr_meanKE_eV,calibSyst_frac,conditionSyst_frac,"
            "bgOffShift_eV,dominantError,noiseLimited\n"
        )
        row = "5,I+He_5,146.9,100,90,1.0,0.9,0.95,0.5,0.01,0.02,0.04,0.06,0.0,calib,0\n"
        p = tmp_path / "bad.csv"
        p.write_text(header + row, encoding="ascii")
        with pytest.raises(ValueError, match="contiguous"):
            load_ihe_ked_reference(p)

    def test_unknown_dominant_error_raises(self, tmp_path):
        header = (
            "n,label,massCenter_u_per_e,N_counts,N_eff,meanKE_eV,"
            "modeKE_eV,medianKE_eV,sigmaKE_eV,statErr_meanKE_eV,"
            "sysErr_meanKE_eV,calibSyst_frac,conditionSyst_frac,"
            "bgOffShift_eV,dominantError,noiseLimited\n"
        )
        row = "0,I+He_0,126.9,100,90,1.0,0.9,0.95,0.5,0.01,0.02,0.04,0.06,0.0,vibes,0\n"
        p = tmp_path / "bad.csv"
        p.write_text(header + row, encoding="ascii")
        with pytest.raises(ValueError, match="dominantError"):
            load_ihe_ked_reference(p)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_ihe_ked.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'i2_helium_md.postprocess.ihe_ked'`.

- [ ] **Step 3: Implement the module + loader**

Create `i2_helium_md/postprocess/ihe_ked.py`:

```python
"""Loaders + sim-side helpers for the I+He_n kinetic-energy reference.

Reads the frozen experimental export under ``data/reference/ihe_ked/``
(exporter: MATLAB ``export_IHe_KED_reference.m``, measurement 17.10.24,
FLIR camera; full provenance in that folder's README/provenance JSON):

* ``IHe_KED_reference.csv`` -- per-fragment first moments <E> of the 3-D
  kinetic-energy distribution for n = 0..17, with the full error model
  (per-point stat + analysis-systematic errors, two correlated fractional
  scale bands, trust flags).
* ``IHe_KED_curves_n{0..4}.csv`` -- trusted curves for n = 0..4: 2-D
  detector projection and 3-D reconstruction, each per unit speed and per
  unit energy, peak-normalized, on a dual ``E_eV``/``v_mps`` axis pair.

Comparison conventions (reference README, binding for consumers):

* compare **mean-to-mean**, never mean-to-peak -- the distributions are
  strongly skewed;
* the MD 3-D ``|v|`` histogram maps onto ``signal_3d_Pv``; the in-plane
  projected speed maps onto ``signal_2d_Pv``;
* mass model ``m(n) = 126.90 + 4.0026 n`` u
  (:func:`i2_helium_md.physics.shell_schedule.complex_mass_amu` -- matches
  the reference's ``massCenter_u_per_e`` column);
* both sides are droplet-rest-frame;
* curves are independently peak-normalized: shapes only, never amplitudes.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..physics.constants import EV, U as U_KG
from ..physics.shell_schedule import complex_mass_amu
from ..simulation.checkpoint import IonCheckpoint
from .csv_contract import validate_columns
from .velocity_distribution import select_final_mass_gate


_REFERENCE_COLUMNS: tuple[str, ...] = (
    "n",
    "label",
    "massCenter_u_per_e",
    "N_counts",
    "N_eff",
    "meanKE_eV",
    "modeKE_eV",
    "medianKE_eV",
    "sigmaKE_eV",
    "statErr_meanKE_eV",
    "sysErr_meanKE_eV",
    "calibSyst_frac",
    "conditionSyst_frac",
    "bgOffShift_eV",
    "dominantError",
    "noiseLimited",
)

# Legal dominantError tags (reference COLUMNS.md).
_LEGAL_DOMINANT_ERROR: frozenset[str] = frozenset(
    {"stat", "analysis", "calib", "bg-structural"}
)

# Fragments with trusted full curves (IHe_KED_curves_n{0..4}.csv).
CURVE_N_MAX: int = 4


@dataclass(frozen=True)
class IHeKedReference:
    """Per-fragment <E> reference table (one array entry per fragment n).

    All energies in eV. ``calib_syst_frac`` / ``condition_syst_frac`` are
    fractional **correlated** scale bands (they shift the whole curve, not
    individual points; never fold them into per-point error bars).

    Attributes
    ----------
    n : np.ndarray, shape (Nn,), int
        He count, contiguous ascending from 0 (Nn = 18 in the frozen file).
    label : np.ndarray, shape (Nn,), str
        Fragment name ``I+He_n``.
    mass_center_u : np.ndarray, shape (Nn,), float
        Fragment mass in u: 126.90 + 4.0026 n.
    N_counts, N_eff : np.ndarray, shape (Nn,), float
        Raw gate counts / background-subtracted in-mask counts.
    mean_KE_eV : np.ndarray, shape (Nn,), float
        THE reference values -- first moments of the unsmoothed 3-D P(E).
    mode_KE_eV, median_KE_eV, sigma_KE_eV : np.ndarray, shape (Nn,), float
        Shape descriptors (sigma is a physical width, NOT an error bar).
    stat_err_mean_KE_eV, sys_err_mean_KE_eV : np.ndarray, shape (Nn,), float
        Per-point uncertainties of the mean; combine in quadrature.
    calib_syst_frac, condition_syst_frac : np.ndarray, shape (Nn,), float
        Correlated fractional scale bands (multiply by ``mean_KE_eV``).
    bg_off_shift_eV : np.ndarray, shape (Nn,), float
        Diagnostic-only background-off sensitivity; not an uncertainty.
    dominant_error : np.ndarray, shape (Nn,), str
        Largest error source per row (one of stat/analysis/calib/
        bg-structural).
    noise_limited : np.ndarray, shape (Nn,), bool
        True would flag <E> as an upper bound (0 for all frozen rows).
    source_path : Path
        Resolved path of the loaded CSV.
    """

    n: np.ndarray
    label: np.ndarray
    mass_center_u: np.ndarray
    N_counts: np.ndarray
    N_eff: np.ndarray
    mean_KE_eV: np.ndarray
    mode_KE_eV: np.ndarray
    median_KE_eV: np.ndarray
    sigma_KE_eV: np.ndarray
    stat_err_mean_KE_eV: np.ndarray
    sys_err_mean_KE_eV: np.ndarray
    calib_syst_frac: np.ndarray
    condition_syst_frac: np.ndarray
    bg_off_shift_eV: np.ndarray
    dominant_error: np.ndarray
    noise_limited: np.ndarray
    source_path: Path

    @property
    def point_err_eV(self) -> np.ndarray:
        """Per-point error of the mean: sqrt(stat^2 + sys^2) [eV]."""
        return np.sqrt(
            self.stat_err_mean_KE_eV ** 2 + self.sys_err_mean_KE_eV ** 2
        )

    @property
    def gold_mask(self) -> np.ndarray:
        """Calib-limited ("gold") fragments: a disagreement there is real
        physics beyond the two correlated bands (reference README)."""
        return self.dominant_error == "calib"


def load_ihe_ked_reference(path: str | Path) -> IHeKedReference:
    """Load ``IHe_KED_reference.csv`` into an :class:`IHeKedReference`.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        Missing/extra columns; non-finite numeric entries; ``n`` not
        contiguous ascending from 0; a non-positive ``meanKE_eV``; a
        fractional band outside [0, 1]; an unknown ``dominantError`` tag;
        or ``noiseLimited`` not in {0, 1}.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"IHe KED reference file not found: {p.resolve()}"
        )

    frame = pd.read_csv(p)
    validate_columns(
        list(frame.columns),
        _REFERENCE_COLUMNS,
        file_label=f"IHe KED reference {p.name}",
    )

    def _finite(name: str) -> np.ndarray:
        values = np.asarray(frame[name].to_numpy(), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"IHe KED reference {p.name} column {name!r} has "
                f"non-finite entries: {values}."
            )
        return values

    n_float = _finite("n")
    n_int = np.rint(n_float).astype(int)
    if n_float.size == 0 or np.any(n_float != n_int) or not np.array_equal(
        n_int, np.arange(n_int.size)
    ):
        raise ValueError(
            f"IHe KED reference {p.name} n column must be contiguous "
            f"integers ascending from 0; got {n_float.tolist()}."
        )

    mean_KE_eV = _finite("meanKE_eV")
    if np.any(mean_KE_eV <= 0.0):
        raise ValueError(
            f"IHe KED reference {p.name} has non-positive meanKE_eV: "
            f"{mean_KE_eV[mean_KE_eV <= 0.0]}."
        )

    for frac_name in ("calibSyst_frac", "conditionSyst_frac"):
        frac = _finite(frac_name)
        if np.any((frac < 0.0) | (frac > 1.0)):
            raise ValueError(
                f"IHe KED reference {p.name} {frac_name} outside [0, 1]: "
                f"{frac}."
            )

    dominant = np.asarray(frame["dominantError"].to_numpy(), dtype=str)
    unknown = set(dominant.tolist()) - _LEGAL_DOMINANT_ERROR
    if unknown:
        raise ValueError(
            f"IHe KED reference {p.name} has unknown dominantError tags "
            f"{sorted(unknown)}; legal: {sorted(_LEGAL_DOMINANT_ERROR)}."
        )

    noise = _finite("noiseLimited")
    if not np.all(np.isin(noise, (0.0, 1.0))):
        raise ValueError(
            f"IHe KED reference {p.name} noiseLimited must be 0 or 1; "
            f"got {noise}."
        )

    return IHeKedReference(
        n=n_int,
        label=np.asarray(frame["label"].to_numpy(), dtype=str),
        mass_center_u=_finite("massCenter_u_per_e"),
        N_counts=_finite("N_counts"),
        N_eff=_finite("N_eff"),
        mean_KE_eV=mean_KE_eV,
        mode_KE_eV=_finite("modeKE_eV"),
        median_KE_eV=_finite("medianKE_eV"),
        sigma_KE_eV=_finite("sigmaKE_eV"),
        stat_err_mean_KE_eV=_finite("statErr_meanKE_eV"),
        sys_err_mean_KE_eV=_finite("sysErr_meanKE_eV"),
        calib_syst_frac=_finite("calibSyst_frac"),
        condition_syst_frac=_finite("conditionSyst_frac"),
        bg_off_shift_eV=_finite("bgOffShift_eV"),
        dominant_error=dominant,
        noise_limited=noise.astype(bool),
        source_path=p.resolve(),
    )
```

(The imports `EV`, `U_KG`, `complex_mass_amu`, `IonCheckpoint`, `select_final_mass_gate` are used by Tasks 3–4; leave them in place now so the module header is final. If flake8 blocks unused imports at this stage, add them in Task 4 instead.)

In `i2_helium_md/postprocess/__init__.py` add:

```python
from .ihe_ked import (
    IHeKedReference,
    load_ihe_ked_reference,
)
```

and append `"IHeKedReference", "load_ihe_ked_reference"` to `__all__`.

- [ ] **Step 4: Run tests to verify pass**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_ihe_ked.py`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add i2_helium_md/postprocess/ihe_ked.py i2_helium_md/postprocess/__init__.py tests/test_ihe_ked.py
git commit -m "feat(postprocess): ihe_ked reference-table loader with full error model"
```

---

### Task 3: `ihe_ked.py` — trusted-curve loader

**Files:**
- Modify: `i2_helium_md/postprocess/ihe_ked.py`
- Modify: `i2_helium_md/postprocess/__init__.py`
- Test: `tests/test_ihe_ked.py`

**Interfaces:**
- Produces: `IHeKedCurve` and `load_ihe_ked_curve(directory: str | Path, n: int) -> IHeKedCurve`. Task 5 uses attributes `v_mps`, `signal_2d_Pv`, `signal_3d_Pv` and the `IHe_KED_curves_n{n}.csv` naming convention.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ihe_ked.py`:

```python
from i2_helium_md.postprocess.ihe_ked import (  # noqa: E402 (grouped here)
    IHeKedCurve,
    load_ihe_ked_curve,
)


class TestLoadIHeKedCurve:
    def test_real_n0_contract(self):
        curve = load_ihe_ked_curve(IHE_KED_DIR, 0)
        assert isinstance(curve, IHeKedCurve)
        assert curve.n == 0
        # Axes are finite and strictly ascending.
        assert np.all(np.isfinite(curve.E_eV))
        assert np.all(np.diff(curve.E_eV) > 0)
        assert np.all(np.diff(curve.v_mps) > 0)
        # n=0 Abel-center spike cut: 3-D columns are NaN below 0.4 eV,
        # while the 2-D columns cover the full detector range.
        low = curve.E_eV < 0.4
        assert np.all(np.isnan(curve.signal_3d_Pv[low]))
        assert np.all(np.isfinite(curve.signal_2d_Pv))
        # Above the cut the 3-D reconstruction is real data.
        assert np.isfinite(curve.signal_3d_Pv[~low]).any()

    def test_all_five_fragments_load(self):
        for n in range(5):
            curve = load_ihe_ked_curve(IHE_KED_DIR, n)
            assert curve.n == n

    def test_out_of_range_n_raises(self):
        with pytest.raises(ValueError, match="n must be"):
            load_ihe_ked_curve(IHE_KED_DIR, 5)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_ihe_ked_curve(tmp_path, 0)

    def test_missing_column_raises(self, tmp_path):
        p = tmp_path / "IHe_KED_curves_n0.csv"
        p.write_text("E_eV,v_mps\n0.1,100\n0.2,140\n", encoding="ascii")
        with pytest.raises(ValueError):
            load_ihe_ked_curve(tmp_path, 0)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_ihe_ked.py -k Curve`
Expected: FAIL with `ImportError: cannot import name 'IHeKedCurve'`.

- [ ] **Step 3: Implement**

Append to `i2_helium_md/postprocess/ihe_ked.py`:

```python
_CURVE_COLUMNS: tuple[str, ...] = (
    "E_eV",
    "v_mps",
    "signal_2d_Pv",
    "signal_2d_PE",
    "signal_3d_Pv",
    "signal_3d_PE",
)


@dataclass(frozen=True)
class IHeKedCurve:
    """One trusted per-fragment curve file (n = 0..4).

    Both representations on one shared axis pair. Each signal column is
    independently peak-normalized (smoothed envelope = 1): compare shapes,
    never amplitudes. NaN marks excluded regions (3-D columns only:
    Abel-center spike / low-E cut); matplotlib renders them as gaps.

    Attributes
    ----------
    n : int
        He count of the fragment (0..4).
    E_eV : np.ndarray, shape (M,)
        Kinetic energy axis (mass-independent detector mapping).
    v_mps : np.ndarray, shape (M,)
        Speed axis of the mass-m(n) complex; differs per fragment at
        equal E.
    signal_2d_Pv, signal_2d_PE : np.ndarray, shape (M,)
        2-D detector projection per unit speed / per unit energy. The
        coordinate is the projected in-plane speed (sim counterpart:
        ``compute_final_velocity_histogram(..., projected=True)``).
    signal_3d_Pv, signal_3d_PE : np.ndarray, shape (M,)
        Reconstructed 3-D speed / energy distribution (sim counterpart:
        the 3-D ``|v|`` histogram).
    source_path : Path
        Resolved path of the loaded CSV.
    """

    n: int
    E_eV: np.ndarray
    v_mps: np.ndarray
    signal_2d_Pv: np.ndarray
    signal_2d_PE: np.ndarray
    signal_3d_Pv: np.ndarray
    signal_3d_PE: np.ndarray
    source_path: Path


def load_ihe_ked_curve(directory: str | Path, n: int) -> IHeKedCurve:
    """Load ``IHe_KED_curves_n{n}.csv`` for one fragment n in 0..4.

    Raises
    ------
    ValueError
        ``n`` outside 0..``CURVE_N_MAX``; missing/extra columns; a
        non-finite axis value; a non-ascending axis; or a signal column
        with no finite entries at all.
    FileNotFoundError
        If the curve file does not exist in ``directory``.
    """
    if not (0 <= int(n) <= CURVE_N_MAX):
        raise ValueError(
            f"n must be in 0..{CURVE_N_MAX} (trusted curves); got {n}."
        )
    p = Path(directory) / f"IHe_KED_curves_n{int(n)}.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"IHe KED curve file not found: {p.resolve()}"
        )

    frame = pd.read_csv(p)
    validate_columns(
        list(frame.columns),
        _CURVE_COLUMNS,
        file_label=f"IHe KED curve {p.name}",
    )

    def _axis(name: str) -> np.ndarray:
        values = np.asarray(frame[name].to_numpy(), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"IHe KED curve {p.name} axis {name!r} has non-finite "
                f"entries."
            )
        if np.any(np.diff(values) <= 0.0):
            raise ValueError(
                f"IHe KED curve {p.name} axis {name!r} must be strictly "
                f"ascending."
            )
        return values

    def _signal(name: str) -> np.ndarray:
        values = np.asarray(frame[name].to_numpy(), dtype=float)
        if not np.isfinite(values).any():
            raise ValueError(
                f"IHe KED curve {p.name} signal {name!r} has no finite "
                f"entries."
            )
        if np.any(np.isinf(values)):
            raise ValueError(
                f"IHe KED curve {p.name} signal {name!r} contains inf."
            )
        return values

    return IHeKedCurve(
        n=int(n),
        E_eV=_axis("E_eV"),
        v_mps=_axis("v_mps"),
        signal_2d_Pv=_signal("signal_2d_Pv"),
        signal_2d_PE=_signal("signal_2d_PE"),
        signal_3d_Pv=_signal("signal_3d_Pv"),
        signal_3d_PE=_signal("signal_3d_PE"),
        source_path=p.resolve(),
    )
```

In `i2_helium_md/postprocess/__init__.py`: extend the `.ihe_ked` import block and `__all__` with `IHeKedCurve`, `load_ihe_ked_curve`.

- [ ] **Step 4: Run tests to verify pass**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_ihe_ked.py`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add i2_helium_md/postprocess/ihe_ked.py i2_helium_md/postprocess/__init__.py tests/test_ihe_ked.py
git commit -m "feat(postprocess): ihe_ked trusted-curve loader (n=0..4, 2D/3D, dual axes)"
```

---

### Task 4: `ihe_ked.py` — sim-side ⟨E⟩ and gate-count helpers

**Files:**
- Modify: `i2_helium_md/postprocess/ihe_ked.py`
- Modify: `i2_helium_md/postprocess/__init__.py`
- Test: `tests/test_ihe_ked.py`

**Interfaces:**
- Consumes: `select_final_mass_gate` (Task 1), `complex_mass_amu`, constants `EV`, `U_KG`.
- Produces (Task 5 relies on these exact names):
  - `speed_mps_of_energy_eV(energy_eV, mass_amu) -> float | np.ndarray`
  - `FragmentMeanKE(n, mass_amu, mean_KE_eV, stat_err_mean_KE_eV, v_of_mean_E_mps, num_atoms_used)`
  - `fragment_mean_kinetic_energy(ion, n, *, mass_tolerance_amu=0.5, require_outside=True) -> FragmentMeanKE` (raises `ValueError` on an empty gate)
  - `fragment_gate_counts(ion, n_values, *, mass_tolerance_amu=0.5, require_outside=True) -> np.ndarray` (int, same shape as `n_values`)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ihe_ked.py`. Copy the `_make_ion` helper from `tests/test_velocity_distribution.py` verbatim (lines 28–90 there; it builds a tiny synthetic `IonCheckpoint` from `final_speeds_per_atom`, `masses_amu_per_atom`, `b_outside`) — tests must not share private helpers across test modules, so the copy is intentional. Add the imports it needs (`from i2_helium_md.physics.constants import U as U_KG`, `from i2_helium_md.simulation.checkpoint import IonCheckpoint`).

```python
from i2_helium_md.physics.constants import EV
from i2_helium_md.physics.shell_schedule import complex_mass_amu
from i2_helium_md.postprocess.ihe_ked import (  # noqa: E402
    FragmentMeanKE,
    fragment_gate_counts,
    fragment_mean_kinetic_energy,
    speed_mps_of_energy_eV,
)


def _energy_eV_of_speed_Aps(speed_Aps: float, mass_amu: float) -> float:
    """Independent re-derivation: E = 1/2 m v^2, v in m/s (1 A/ps = 100 m/s)."""
    from i2_helium_md.physics.constants import U as _U
    v_mps = speed_Aps * 100.0
    return 0.5 * mass_amu * _U * v_mps * v_mps / EV


class TestFragmentMeanKineticEnergy:
    def test_analytic_mean_and_stat_err(self):
        # Two n=1 atoms (131 amu, gated by m(1)=130.9026) at 2 and 4 A/ps;
        # one n=0 atom that must not contaminate the gate.
        ion = _make_ion(
            num_molecules=2,
            final_speeds_per_atom=np.array([2.0, 4.0, 3.0, 0.0]),
            masses_amu_per_atom=np.array([131.0, 131.0, 127.0, 127.0]),
        )
        result = fragment_mean_kinetic_energy(ion, 1)
        assert isinstance(result, FragmentMeanKE)
        m1 = complex_mass_amu(1)
        e1 = _energy_eV_of_speed_Aps(2.0, m1)
        e2 = _energy_eV_of_speed_Aps(4.0, m1)
        expected_mean = 0.5 * (e1 + e2)
        # Analytical port: tight tolerance.
        assert result.mean_KE_eV == pytest.approx(expected_mean, rel=1e-12)
        expected_err = np.std([e1, e2], ddof=1) / np.sqrt(2.0)
        assert result.stat_err_mean_KE_eV == pytest.approx(
            expected_err, rel=1e-12,
        )
        assert result.num_atoms_used == 2
        assert result.mass_amu == pytest.approx(m1)

    def test_v_of_mean_E_round_trip(self):
        ion = _make_ion(
            num_molecules=1,
            final_speeds_per_atom=np.array([3.0, 0.0]),
            masses_amu_per_atom=np.array([131.0, 127.0]),
        )
        result = fragment_mean_kinetic_energy(ion, 1)
        # A single atom at 3 A/ps: v(<E>) must be exactly 300 m/s.
        assert result.v_of_mean_E_mps == pytest.approx(300.0, rel=1e-12)
        # And the standalone converter agrees.
        assert speed_mps_of_energy_eV(
            result.mean_KE_eV, result.mass_amu
        ) == pytest.approx(300.0, rel=1e-12)

    def test_single_atom_has_zero_stat_err(self):
        ion = _make_ion(
            num_molecules=1,
            final_speeds_per_atom=np.array([3.0, 0.0]),
            masses_amu_per_atom=np.array([131.0, 127.0]),
        )
        result = fragment_mean_kinetic_energy(ion, 1)
        assert result.stat_err_mean_KE_eV == 0.0

    def test_empty_gate_raises(self):
        ion = _make_ion(
            num_molecules=1,
            final_speeds_per_atom=np.array([3.0, 0.0]),
            masses_amu_per_atom=np.array([127.0, 127.0]),
        )
        with pytest.raises(ValueError, match="No atoms"):
            fragment_mean_kinetic_energy(ion, 1)


class TestFragmentGateCounts:
    def test_counts_per_gate(self):
        ion = _make_ion(
            num_molecules=2,
            final_speeds_per_atom=np.array([1.0, 2.0, 3.0, 4.0]),
            masses_amu_per_atom=np.array([127.0, 131.0, 131.0, 135.0]),
        )
        counts = fragment_gate_counts(ion, np.arange(4))
        np.testing.assert_array_equal(counts, [1, 2, 1, 0])

    def test_outside_filter_applies(self):
        ion = _make_ion(
            num_molecules=2,
            final_speeds_per_atom=np.array([1.0, 2.0, 3.0, 4.0]),
            masses_amu_per_atom=np.array([131.0, 131.0, 131.0, 131.0]),
            b_outside=np.array([True, False]),
        )
        counts = fragment_gate_counts(ion, np.arange(3))
        np.testing.assert_array_equal(counts, [0, 2, 0])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_ihe_ked.py -k "Fragment"`
Expected: FAIL with `ImportError: cannot import name 'FragmentMeanKE'`.

- [ ] **Step 3: Implement**

Append to `i2_helium_md/postprocess/ihe_ked.py`:

```python
def speed_mps_of_energy_eV(energy_eV, mass_amu):
    """Speed [m/s] of a mass-``mass_amu`` complex with kinetic energy
    ``energy_eV`` [eV]: ``v = sqrt(2 E / m)``. Scalar or array in,
    same shape out."""
    return np.sqrt(2.0 * np.asarray(energy_eV, dtype=float) * EV
                   / (float(mass_amu) * U_KG))


@dataclass(frozen=True)
class FragmentMeanKE:
    """Simulated per-fragment mean kinetic energy (droplet rest frame).

    The MD counterpart of one ``mean_KE_eV`` row of the reference table:
    the first moment of the mass-gated final-speed ensemble, computed with
    the reference's mass model m(n) so the comparison is mean-to-mean on
    identical conventions.

    Attributes
    ----------
    n : int
        He count of the fragment gate.
    mass_amu : float
        m(n) = 126.90 + 4.0026 n used for both the gate and the energy.
    mean_KE_eV : float
        Ensemble mean of E = 1/2 m(n) |v_final|^2 [eV].
    stat_err_mean_KE_eV : float
        Standard error of the mean, sample-std(ddof=1)/sqrt(N); 0.0 for
        N = 1 (a single atom has no scatter estimate).
    v_of_mean_E_mps : float
        sqrt(2 <E> / m(n)) [m/s] -- the marker position on the curve
        overlays. NOT the mean speed <|v|>.
    num_atoms_used : int
        Atoms that passed the mass + outside gate.
    """

    n: int
    mass_amu: float
    mean_KE_eV: float
    stat_err_mean_KE_eV: float
    v_of_mean_E_mps: float
    num_atoms_used: int


def fragment_mean_kinetic_energy(
    ion: IonCheckpoint,
    n: int,
    *,
    mass_tolerance_amu: float = 0.5,
    require_outside: bool = True,
) -> FragmentMeanKE:
    """Simulated <E> of the I+He_n fragment gate (mean-to-mean observable).

    Raises
    ------
    ValueError
        If no atoms pass the mass + outside gate (callers typically catch
        this and skip the fragment), or ``n < 0``.
    """
    if int(n) < 0:
        raise ValueError(f"n must be >= 0, got {n}.")
    mass_amu = float(complex_mass_amu(int(n)))
    select = select_final_mass_gate(
        ion,
        mass_amu=mass_amu,
        mass_tolerance_amu=mass_tolerance_amu,
        require_outside=require_outside,
    )
    num_used = int(np.count_nonzero(select))
    if num_used == 0:
        raise ValueError(
            f"No atoms in the I+He_{int(n)} gate "
            f"(mass {mass_amu:.4f} +/- {mass_tolerance_amu} amu"
            f"{', outside only' if require_outside else ''})."
        )

    speed_mps = 100.0 * np.sqrt(
        np.asarray(ion.velocities_final_x)[select] ** 2
        + np.asarray(ion.velocities_final_y)[select] ** 2
        + np.asarray(ion.velocities_final_z)[select] ** 2
    )
    energy_eV = 0.5 * mass_amu * U_KG * speed_mps ** 2 / EV
    mean_eV = float(energy_eV.mean())
    stat_err_eV = (
        float(energy_eV.std(ddof=1) / np.sqrt(num_used))
        if num_used > 1 else 0.0
    )
    return FragmentMeanKE(
        n=int(n),
        mass_amu=mass_amu,
        mean_KE_eV=mean_eV,
        stat_err_mean_KE_eV=stat_err_eV,
        v_of_mean_E_mps=float(speed_mps_of_energy_eV(mean_eV, mass_amu)),
        num_atoms_used=num_used,
    )


def fragment_gate_counts(
    ion: IonCheckpoint,
    n_values,
    *,
    mass_tolerance_amu: float = 0.5,
    require_outside: bool = True,
) -> np.ndarray:
    """Atom count per I+He_n mass gate, for the abundance comparison.

    Uses the same gate convention as every other mass-gated diagnostic
    (``select_final_mass_gate`` at m(n)); gates are disjoint because the
    He spacing (4.0026 amu) exceeds twice the default tolerance. Note the
    Tier-2 scoring path (:func:`~i2_helium_md.postprocess.size_distribution.
    compute_terminal_shell_distribution`) reads the v7 ``n_shell`` field
    instead; for v7 biphasic runs the two agree because mass jumps track
    ``n_shell``, while this mass-gate variant also works for pre-v7 and
    hard-sphere runs.

    Returns
    -------
    np.ndarray
        Integer counts, same shape as ``n_values``.
    """
    n_arr = np.asarray(n_values, dtype=int)
    counts = np.zeros(n_arr.shape, dtype=int)
    for i, n in enumerate(n_arr.ravel()):
        select = select_final_mass_gate(
            ion,
            mass_amu=float(complex_mass_amu(int(n))),
            mass_tolerance_amu=mass_tolerance_amu,
            require_outside=require_outside,
        )
        counts.ravel()[i] = int(np.count_nonzero(select))
    return counts
```

In `i2_helium_md/postprocess/__init__.py`: extend the `.ihe_ked` block and `__all__` with `FragmentMeanKE`, `fragment_mean_kinetic_energy`, `fragment_gate_counts`, `speed_mps_of_energy_eV`.

- [ ] **Step 4: Run tests to verify pass**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_ihe_ked.py`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add i2_helium_md/postprocess/ihe_ked.py i2_helium_md/postprocess/__init__.py tests/test_ihe_ked.py
git commit -m "feat(postprocess): sim-side fragment <E> and gate-count helpers for ihe_ked"
```

---

### Task 5: Rewire `plot_run_summary.py` (new sections, retire vmi_summary overlay)

**Files:**
- Modify: `scripts/post_processing/plot_run_summary.py`
- Test: `tests/test_ihe_ked.py` (smoke tests at the end)

**Interfaces:**
- Consumes everything Tasks 1–4 produced, plus existing `load_he_abundance_reference`, `mass_spectrum`, `moving_mean`, `normalise_trace`, `compute_final_velocity_histogram`.
- Produces: section builders `_section_ihe_ked_mean_energy(ion, ked_ref)`, `_section_ihe_ked_curves(ion, ked_dir, ked_ref, representation)`, `_section_mass_spectrum(ion, abundance)` — the smoke tests call these directly.

- [ ] **Step 1: Write the failing smoke tests**

Append to `tests/test_ihe_ked.py`:

```python
import importlib.util
import sys

import matplotlib


def _load_run_summary_module():
    matplotlib.use("Agg", force=True)
    script = (
        PROJECT_ROOT / "scripts" / "post_processing" / "plot_run_summary.py"
    )
    spec = importlib.util.spec_from_file_location(
        "plot_run_summary_under_test", script
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _smoke_ion() -> IonCheckpoint:
    """A tiny ensemble populating the n = 0..2 gates with plausible speeds."""
    speeds = np.array([24.0, 20.0, 10.0, 8.0, 6.0, 5.0])  # A/ps
    masses = np.array([127.0, 127.0, 131.0, 131.0, 135.0, 135.0])
    return _make_ion(
        num_molecules=3,
        final_speeds_per_atom=speeds,
        masses_amu_per_atom=masses,
    )


class TestRunSummaryIHeKedSections:
    def test_mean_energy_section_builds(self):
        import matplotlib.pyplot as plt
        mod = _load_run_summary_module()
        ked_ref = load_ihe_ked_reference(REFERENCE_CSV)
        fig = mod._section_ihe_ked_mean_energy(_smoke_ion(), ked_ref)
        assert fig is not None
        plt.close("all")

    def test_curves_sections_build_both_representations(self):
        import matplotlib.pyplot as plt
        mod = _load_run_summary_module()
        ked_ref = load_ihe_ked_reference(REFERENCE_CSV)
        for representation in ("2d", "3d"):
            fig = mod._section_ihe_ked_curves(
                _smoke_ion(), IHE_KED_DIR, ked_ref, representation
            )
            assert fig is not None
        plt.close("all")

    def test_mass_spectrum_with_abundance_builds(self):
        import matplotlib.pyplot as plt
        from i2_helium_md.postprocess import load_he_abundance_reference
        mod = _load_run_summary_module()
        abundance = load_he_abundance_reference(
            PROJECT_ROOT / "data" / "reference"
            / "integrated_i_he_abundance.csv"
        )
        fig = mod._section_mass_spectrum(_smoke_ion(), abundance)
        assert fig is not None
        plt.close("all")

    def test_old_vmi_section_is_gone(self):
        mod = _load_run_summary_module()
        assert not hasattr(mod, "_section_radial_velocity")
        assert not hasattr(mod, "VMI_REF_HE_PATH")
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_ihe_ked.py -k RunSummary`
Expected: FAIL — `_section_ihe_ked_mean_energy` does not exist; `_section_mass_spectrum` has the wrong signature; `_section_radial_velocity` still exists.

- [ ] **Step 3: Edit `scripts/post_processing/plot_run_summary.py`**

**(a) Module docstring:** in the consolidated-scripts list, replace the `simulation_image_only_trajectories.m` line's companion entry for the VMI overlay — change the line

```
* ``post_process_single_pulse_paper_v3.m``     -> 1D and 2D polar VMI panels
```

to stay as is, and ADD after the list:

```
The legacy ``simulation_image.m`` velocity overlay against the
``vmi_summary`` CSVs was retired 2026-07-14 in favor of the frozen
``data/reference/ihe_ked/`` per-fragment reference (mean-KE table +
trusted 2-D/3-D curves for n = 0..4); see
``docs/superpowers/specs/2026-07-14-ihe-ked-run-summary-design.md``.
```

**(b) Imports:** in the `from i2_helium_md.postprocess import (...)` block: remove `load_vmi_reference`; add `fragment_gate_counts`, `fragment_mean_kinetic_energy`, `load_he_abundance_reference`, `load_ihe_ked_curve`, `load_ihe_ked_reference`, `speed_mps_of_energy_eV` (keep alphabetical order). Add after that block:

```python
from i2_helium_md.physics.shell_schedule import complex_mass_amu  # noqa: E402
```

**(c) Plot-tuning constants:** remove `MASS_I_HE_AMU`, `MASS_I_HE2_AMU` if now unused elsewhere in the file — check first: `PAPER_V2_MASS_AMU` (131.0) and `_section_mass_resolved` use them; `_section_mass_resolved` keeps `MASS_I_HE_AMU`/`MASS_I_HE2_AMU`, so KEEP both. Add:

```python
# ihe_ked overlay panels: cover the full detector range (crop edge
# 354 px = 3052 m/s at m = 127; reference README).
IHE_KED_HIST_V_MAX_APS = 31.0
IHE_KED_HIST_NUM_BINS = int(round(IHE_KED_HIST_V_MAX_APS / HIST_BIN_WIDTH_APS))
IHE_KED_PLOT_V_MAX_MPS = 3100.0
```

**(d) USER SETTINGS:** delete the `VMI_REF_HE_PATH`, `VMI_REF_GAS_PATH`, `VMI_REF_HE_HIGH_SNR_PATH` settings (and their comment lines, including the commented `VMI_REF_*` lines in the 9 A HeDFT block — replace those two commented lines with `# IHE_KED_REFERENCE_DIR = None` and `# ABUNDANCE_REF_PATH = None`). Add:

```python
# Directory holding the frozen I+He_n kinetic-energy reference
# (IHe_KED_reference.csv + IHe_KED_curves_n{0..4}.csv). ``None`` skips
# the ihe_ked mean-energy and curve-overlay sections.
IHE_KED_REFERENCE_DIR: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "ihe_ked"
)

# Experimental I+He_n abundance CSV for the mass-spectrum side-by-side.
# ``None`` keeps the plain simulated mass spectrum.
ABUNDANCE_REF_PATH: Path | None = (
    PROJECT_ROOT / "data" / "reference" / "integrated_i_he_abundance.csv"
)
```

**(e) `main()`:** remove the `vmi_ref_*` path handling, `load_vmi_reference` calls, and the `vmi_ref_*` fields of `args`. Add:

```python
    ihe_ked_dir = Path(IHE_KED_REFERENCE_DIR) if IHE_KED_REFERENCE_DIR else None
    abundance_ref_path = Path(ABUNDANCE_REF_PATH) if ABUNDANCE_REF_PATH else None
```

add `ihe_ked_dir=ihe_ked_dir, abundance_ref_path=abundance_ref_path,` to the `SimpleNamespace`, and after the `hedft = ...` load:

```python
    ked_ref = (
        load_ihe_ked_reference(ihe_ked_dir / "IHe_KED_reference.csv")
        if ihe_ked_dir else None
    )
    abundance = (
        load_he_abundance_reference(abundance_ref_path)
        if abundance_ref_path else None
    )
```

In the `sections` list (ion branch): replace

```python
                ("mass_spectrum",
                 lambda: _section_mass_spectrum(ion)),
                ("radial_velocity_with_vmi",
                 lambda: _section_radial_velocity(
                     ion, vmi_he, vmi_gas, vmi_he_high_snr)),
```

with

```python
                ("mass_spectrum",
                 lambda: _section_mass_spectrum(ion, abundance)),
                ("ihe_ked_mean_energy",
                 lambda: _section_ihe_ked_mean_energy(ion, ked_ref)),
                ("ihe_ked_curves_3d",
                 lambda: _section_ihe_ked_curves(
                     ion, ihe_ked_dir, ked_ref, "3d")),
                ("ihe_ked_curves_2d",
                 lambda: _section_ihe_ked_curves(
                     ion, ihe_ked_dir, ked_ref, "2d")),
```

**(f) `_section_metadata`:** replace the three `args.vmi_ref_*` blocks with:

```python
    if args.ihe_ked_dir:
        refs.append(f"IHe KED (mean-KE + curves): {args.ihe_ked_dir}")
    if args.abundance_ref_path:
        refs.append(f"I+He_n abundance: {args.abundance_ref_path}")
```

**(g) Delete `_section_radial_velocity` entirely** (lines ~423–513).

**(h) Replace `_section_mass_spectrum`:**

```python
def _section_mass_spectrum(ion, abundance) -> plt.Figure:
    """Final ion mass spectrum; side-by-side with the experimental
    I+He_n abundance when the reference is configured.

    Sim fractions use the same mass gates as the velocity diagnostics
    (m(n) +/- 0.5 amu, outside ions only), so a detected-fraction bar is
    directly the mass-gated ensemble the other ihe_ked sections draw from.
    """
    if abundance is None:
        spec = mass_spectrum(ion, bin_width_amu=1.0)
        fig, ax = plt.subplots(figsize=(8.0, 4.0), constrained_layout=True)
        ax.bar(spec.bin_centers_amu, spec.counts, width=0.9,
               edgecolor="black", linewidth=0.5)
        ax.set(title="Final ion mass spectrum",
               xlabel="m / u", ylabel="count")
        ax.set_xlim(left=MASS_I - 1, right=MASS_SPECTRUM_MAX_AMU + 1)
        ax.set_xticks(np.arange(MASS_I, MASS_SPECTRUM_MAX_AMU + 1, 4))
        return fig

    counts = fragment_gate_counts(ion, abundance.n)
    total = counts.sum()
    if total == 0:
        raise _SectionSkipped("no outside ions in any I+He_n mass gate")
    sim_percent = 100.0 * counts / total
    exp_percent = 100.0 * abundance.ion_fraction

    fig, ax = plt.subplots(figsize=(9.5, 4.5), constrained_layout=True)
    width = 0.4
    ax.bar(abundance.n - width / 2, exp_percent, width=width,
           color="tab:blue", label="experiment (ionPercent)")
    ax.bar(abundance.n + width / 2, sim_percent, width=width,
           color="tab:red", label=f"simulation (N={int(total)} atoms)")
    ax.set(title="Final I$^+$He$_n$ size distribution vs experimental abundance",
           xlabel="n (attached He atoms)", ylabel="fraction / %")
    ax.set_xticks(abundance.n)
    ax.legend(frameon=False)
    return fig
```

**(i) Add the two new section builders** (place after `_section_mass_spectrum`):

```python
def _section_ihe_ked_mean_energy(ion, ked_ref) -> plt.Figure:
    """Mean kinetic energy per fragment: sim vs experiment, mean-to-mean.

    Error model per the reference README: per-point error is
    sqrt(stat^2 + sys^2); the calibration and condition bands are
    correlated (they shift the whole experimental curve coherently) and
    are drawn as shaded envelopes, never folded into point errors. Gold
    points are calib-limited: a disagreement there is real physics beyond
    the two bands.
    """
    if ked_ref is None:
        raise _SectionSkipped("IHE_KED_REFERENCE_DIR is None")

    sim_points = []
    for n in ked_ref.n:
        try:
            sim_points.append(fragment_mean_kinetic_energy(ion, int(n)))
        except ValueError:
            continue
    fig, ax = plt.subplots(figsize=(9.5, 5.0), constrained_layout=True)

    mean = ked_ref.mean_KE_eV
    for frac, label, alpha in (
        (ked_ref.calib_syst_frac, "calibration band (correlated)", 0.20),
        (ked_ref.condition_syst_frac, "condition band (correlated)", 0.12),
    ):
        ax.fill_between(ked_ref.n, mean * (1.0 - frac), mean * (1.0 + frac),
                        color="tab:blue", alpha=alpha, linewidth=0,
                        label=label)
    ax.errorbar(ked_ref.n, mean, yerr=ked_ref.point_err_eV, fmt="o",
                color="tab:blue", markersize=4, capsize=2,
                label=r"experiment $\langle E\rangle$ (stat $\oplus$ sys)")
    gold = ked_ref.gold_mask
    ax.plot(ked_ref.n[gold], mean[gold], "o", markersize=10,
            markerfacecolor="none", markeredgecolor="goldenrod",
            markeredgewidth=1.5, label="gold points (calib-limited)")

    if sim_points:
        ax.errorbar(
            [p.n for p in sim_points],
            [p.mean_KE_eV for p in sim_points],
            yerr=[p.stat_err_mean_KE_eV for p in sim_points],
            fmt="s", linestyle="--", color="tab:red", markersize=5,
            capsize=2, label=r"simulation $\langle E\rangle$",
        )
        for p in sim_points:
            ax.annotate(f"N={p.num_atoms_used}",
                        (p.n, p.mean_KE_eV),
                        textcoords="offset points", xytext=(0, 7),
                        fontsize=7, color="tab:red", ha="center")

    ax.set_yscale("log")
    ax.set(title=r"I$^+$He$_n$ mean kinetic energy (mean-to-mean)",
           xlabel="n (attached He atoms)",
           ylabel=r"$\langle E\rangle$ / eV")
    ax.set_xticks(ked_ref.n)
    ax.legend(frameon=False, fontsize=8)
    return fig


def _section_ihe_ked_curves(ion, ked_dir, ked_ref, representation) -> plt.Figure:
    """Per-fragment speed-distribution overlays for n = 0..4.

    representation: "3d" overlays the sim 3-D |v| histogram on the
    reconstructed ``signal_3d_Pv``; "2d" overlays the sim in-plane
    projected speed sqrt(vx^2+vy^2) on the detector projection
    ``signal_2d_Pv``. Reference curves are peak-normalized (smoothed
    envelope = 1): compare shapes, never amplitudes. Vertical markers sit
    at v(<E>) -- not <v> -- on both sides; the reference NaN cuts render
    as gaps.
    """
    if ked_dir is None or ked_ref is None:
        raise _SectionSkipped("IHE_KED_REFERENCE_DIR is None")
    if representation not in ("2d", "3d"):
        raise ValueError(f"representation must be '2d' or '3d', "
                         f"got {representation!r}")

    titles = {
        "3d": ("3-D speed distributions P(v) vs I$^+$He$_n$ reference "
               "(peak-normalized, shapes only)"),
        "2d": ("2-D detector projections (in-plane speed) vs I$^+$He$_n$ "
               "reference (peak-normalized, shapes only)"),
    }
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.0),
                             constrained_layout=True)
    axes = axes.ravel()
    axes[5].set_axis_off()
    fig.suptitle(titles[representation])

    for n in range(5):
        ax = axes[n]
        curve = load_ihe_ked_curve(ked_dir, n)
        ref_signal = getattr(curve, f"signal_{representation}_Pv")
        ax.plot(curve.v_mps, ref_signal, color="tab:blue", linewidth=1.2,
                label="experiment")

        mass_amu = float(complex_mass_amu(n))
        sim_note = None
        try:
            hist = compute_final_velocity_histogram(
                ion, mass_amu=mass_amu,
                num_bins=IHE_KED_HIST_NUM_BINS,
                v_max_Aps=IHE_KED_HIST_V_MAX_APS,
                projected=(representation == "2d"),
            )
        except ValueError:
            sim_note = "sim: no atoms in gate"
        else:
            sim_density = normalise_trace(
                moving_mean(hist.density, HIST_SMOOTHING_WINDOW)
            )
            ax.plot(hist.bin_centers_mps, sim_density, "--",
                    color="tab:red", linewidth=1.4,
                    label=f"simulation (N={hist.num_atoms_used})")

        # v(<E>) markers (labelled v(<E>), not <v>): experiment from the
        # reference table, simulation from the mass-gated ensemble.
        v_exp = speed_mps_of_energy_eV(
            float(ked_ref.mean_KE_eV[n]), mass_amu,
        )
        ax.axvline(v_exp, color="tab:blue", linestyle=":", linewidth=1.0,
                   label=r"exp $v(\langle E\rangle)$")
        try:
            sim_mean = fragment_mean_kinetic_energy(ion, n)
        except ValueError:
            pass
        else:
            ax.axvline(sim_mean.v_of_mean_E_mps, color="tab:red",
                       linestyle=":", linewidth=1.0,
                       label=r"sim $v(\langle E\rangle)$")

        if sim_note:
            ax.annotate(sim_note, (0.97, 0.9), xycoords="axes fraction",
                        ha="right", fontsize=8, color="tab:red")
        ax.set(title=f"n = {n}", xlim=(0.0, IHE_KED_PLOT_V_MAX_MPS),
               ylim=(0.0, 1.25))
        ax.set_xlabel("v / m/s")
        ax.set_ylabel("signal / arb. units")
        if n == 0:
            ax.legend(frameon=False, fontsize=7)
    return fig
```

- [ ] **Step 4: Run the smoke tests, then the full test file**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q tests/test_ihe_ked.py`
Expected: all PASS.

- [ ] **Step 5: Manual end-to-end check (only if a finished run directory exists)**

If `data/runs/single_pulse_droplet` (or another finished run) is present, run the script once and confirm the new PNGs appear and the old `radial_velocity_with_vmi.png` is no longer produced:

`& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' scripts/post_processing/plot_run_summary.py`

Expected console lines include `wrote mass_spectrum`, `wrote ihe_ked_mean_energy`, `wrote ihe_ked_curves_3d`, `wrote ihe_ked_curves_2d`. If no run directory exists, skip — the smoke tests already exercise the builders.

- [ ] **Step 6: Commit**

```bash
git add scripts/post_processing/plot_run_summary.py tests/test_ihe_ked.py
git commit -m "feat(run-summary): ihe_ked mean-energy + curve-overlay sections, abundance side-by-side; retire vmi_summary overlay"
```

---

### Task 6: Full test run + delivery log entry

**Files:**
- Modify: `docs/drag_port/Tier2/drag_migration_log_tier2.md` (append delivery note)

- [ ] **Step 1: Run the full test suite**

Run: `& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' -m pytest -q`
Expected: all PASS (no regressions outside the touched modules). Report any pre-existing failures verbatim without fixing unrelated code.

- [ ] **Step 2: Append the delivery note to the Tier-2 log**

Read the tail of `docs/drag_port/Tier2/drag_migration_log_tier2.md` first and match its entry format (dated heading + short prose). Content to record:

> **2026-07-14 — run-summary comparison layer upgraded to the ihe_ked reference.** New `postprocess/ihe_ked.py` (reference-table loader with the full error model, trusted-curve loader n = 0…4, sim-side `fragment_mean_kinetic_energy` / `fragment_gate_counts` on the shared `select_final_mass_gate` gate). `plot_run_summary.py`: new `ihe_ked_mean_energy` (mean-to-mean, correlated bands, gold points), `ihe_ked_curves_3d` / `ihe_ked_curves_2d` (5-panel overlays with `v(⟨E⟩)` markers; sim 2-D side is the projected in-plane speed), mass spectrum side-by-side with `integrated_i_he_abundance.csv`. Retired the `vmi_summary` overlay section (different-campaign data, ⟨E⟩ ~16 % lower; CSVs and loader kept for `plot_experimental_comparison.py`). Spec: `docs/superpowers/specs/2026-07-14-ihe-ked-run-summary-design.md`; tests: `tests/test_ihe_ked.py`.

- [ ] **Step 3: Commit**

```bash
git add docs/drag_port/Tier2/drag_migration_log_tier2.md
git commit -m "docs(tier2): log ihe_ked run-summary comparison-layer delivery"
```
