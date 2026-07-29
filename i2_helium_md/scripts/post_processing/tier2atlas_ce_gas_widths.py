"""Atlas (C) pre-step P2 — Abel gas channel-width read (zero MD, read-only).

Reads the per-channel KER means and widths (sigma_Q2, sigma_Q3, E_single)
from the committed Abel-inverted gas export
``data/reference/vmi_summary/vmi_iplus_gas.csv`` (measurement 43632,
pyabel ``rIbeta()`` 3-D speed distribution I(v); README pipeline) — the
design-doc §3.5 priors (Hatherly 200 fs per-I⁺ sigma 0.55 / 1.23 eV,
E_single Bounded [0.3, 0.8] eV) are sharpened by the experiment's own
calibration frame before probe registration
(`TIER2_CE_CHANNEL_EXIT_STRIP_DESIGN.md` §9 P2).

Conventions:

- Energy per I⁺ fragment: E = 1/2 * m_I * v^2 with m_I = MASS_I_ION_AMU
  (126.90 u — the committed shell-schedule iodine-ion mass; the §3.5l
  read-(iii) convention).
- KED: P(E) = I(v) * dv/dE = I(v) / (m_I * v) on the export's native grid
  (the Jacobian-correct 1-D energy distribution; normalization free).
- Width fit: bounded least squares (scipy) of a three-Gaussian model
  {slow/single, Q2, Q3} on P(E) over the fit window. The slow component is
  ONE effective Gaussian: the finer 0.27/0.49/0.70 sub-structure quoted in
  the §3.5l record is below this 72-point export's resolving power and is
  deliberately not multi-fit here (over-parameterized on this grid).

Oracle (gate before any new number): the local maxima of the exported I(v)
must reproduce the committed session peak readings 2.26 / 4.11 eV
(RQ8 NB / §3.5l rim-rider withdrawal record) within one grid step.

Stated limitations (carried into the findings):

- The measured sigmas contain the instrument response, the Abel smoothing
  (movmean 3x3 pre-inversion), and the 4 % ``calibSyst_frac`` — they are
  UPPER bounds on the intrinsic channel widths.
- The gas frame is field-free vacuum CE; in-droplet channel widths inherit
  these priors by design assumption (Hatherly: fraction channel-independent).

Pure reader — runs nothing, mutates nothing.

Invocation:

    python scripts/post_processing/tier2atlas_ce_gas_widths.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass

import numpy as np
from scipy.optimize import curve_fit

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from i2_helium_md.physics.constants import EV, MASS_I_ION_AMU, U  # noqa: E402

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

GAS_CSV = PROJECT_ROOT / "data" / "reference" / "vmi_summary" / "vmi_iplus_gas.csv"

# Committed session peak readings (RQ8 NB 2026-07-29; the §3.5l rim-rider
# withdrawal re-confirmed them as the repo-frame channel positions).
COMMITTED_PEAKS_EV = (2.26, 4.11)

FIT_E_MIN_EV = 0.15       # below: Abel-center residue + the v = 0 bin
FIT_E_MAX_EV = 6.5        # above: detector edge tail

# Three-Gaussian bounds: (weight, mu [eV], sigma [eV]) per component.
COMP_BOUNDS = {
    "single": {"mu": (0.20, 0.95), "sigma": (0.05, 0.60)},
    "Q2":     {"mu": (1.60, 2.80), "sigma": (0.15, 1.20)},
    "Q3":     {"mu": (3.40, 5.00), "sigma": (0.30, 2.00)},
}

# Design-doc §3.5 priors printed next to the fit (context, not gates).
HATHERLY_SIGMA_EV = {"Q2": 0.55, "Q3": 1.23}
DESIGN_E_SINGLE_BAND = (0.3, 0.8)
CALIB_SYST_FRAC = 0.04

# ---------------------------------------------------------------------------

_MI_KG = MASS_I_ION_AMU * U


def _energy_eV(v_mps: np.ndarray) -> np.ndarray:
    return 0.5 * _MI_KG * v_mps ** 2 / EV


def _three_gauss(E, w1, mu1, s1, w2, mu2, s2, w3, mu3, s3):
    out = np.zeros_like(E)
    for w, mu, s in ((w1, mu1, s1), (w2, mu2, s2), (w3, mu3, s3)):
        out = out + w * np.exp(-0.5 * ((E - mu) / s) ** 2)
    return out


def main() -> None:
    with open(GAS_CSV, newline="", encoding="utf-8") as fh:
        rows = [(float(r["v_mps"]), float(r["signal_arb"]))
                for r in csv.DictReader(fh)]
    v = np.asarray([r[0] for r in rows])
    iv = np.asarray([r[1] for r in rows])
    E = _energy_eV(v)
    print(f"{GAS_CSV.name}: {v.size} points, v 0-{v.max():.0f} m/s "
          f"(E 0-{E.max():.2f} eV per I+ at m_I {MASS_I_ION_AMU} u); "
          f"grid step {np.diff(v).mean():.1f} m/s.")

    # ---- Oracle: I(v) local maxima reproduce the committed 2.26 / 4.11 eV.
    interior = np.arange(1, v.size - 1)
    is_max = (iv[interior] > iv[interior - 1]) & (iv[interior] >= iv[interior + 1])
    peaks_idx = interior[is_max]
    # Ignore sub-threshold ripples: keep maxima above 10 % of the global max.
    peaks_idx = peaks_idx[iv[peaks_idx] >= 0.10 * iv.max()]
    peaks_E = E[peaks_idx]
    fast = peaks_E[peaks_E > 1.2]
    bad = []
    for target in COMMITTED_PEAKS_EV:
        if fast.size == 0:
            bad.append(f"no I(v) maximum found near {target} eV")
            continue
        j = int(np.argmin(np.abs(fast - target)))
        # one grid step in E at the peak position
        k = peaks_idx[np.argmin(np.abs(peaks_E - fast[j]))]
        dE = abs(E[min(k + 1, E.size - 1)] - E[max(k - 1, 0)]) / 2.0
        if abs(fast[j] - target) > dE:
            bad.append(f"nearest I(v) maximum to {target} eV sits at "
                       f"{fast[j]:.3f} (one-step tol {dE:.3f})")
    if bad:
        print("*** P2 PEAK ORACLE FAILED — do not read the widths ***")
        for line in bad:
            print(f"    {line}")
        return
    print(f"Peak oracle OK: I(v) maxima at "
          f"{', '.join(f'{x:.2f}' for x in sorted(fast))} eV reproduce the "
          f"committed {COMMITTED_PEAKS_EV} within one grid step.")

    # ---- The Jacobian-correct KED and the three-Gaussian width fit.
    sel = (E >= FIT_E_MIN_EV) & (E <= FIT_E_MAX_EV) & (v > 0)
    Ef = E[sel]
    ked = iv[sel] / v[sel]          # P(E) ∝ I(v)/(m v); constant m dropped
    ked = ked / ked.max()

    names = ("single", "Q2", "Q3")
    lo, hi, p0 = [], [], []
    for name, mu0, s0 in (("single", 0.5, 0.2), ("Q2", 2.2, 0.5),
                          ("Q3", 4.1, 1.0)):
        b = COMP_BOUNDS[name]
        lo += [0.0, b["mu"][0], b["sigma"][0]]
        hi += [np.inf, b["mu"][1], b["sigma"][1]]
        p0 += [0.5, mu0, s0]
    popt, pcov = curve_fit(_three_gauss, Ef, ked, p0=p0,
                           bounds=(lo, hi), maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    resid = ked - _three_gauss(Ef, *popt)
    rms = float(np.sqrt(np.mean(resid ** 2)))

    print(f"\n=== P2 three-Gaussian KED fit ({FIT_E_MIN_EV}-{FIT_E_MAX_EV} eV, "
          f"{Ef.size} points, RMS resid {rms:.4f} of peak) ===")
    weights = popt[0::3] * popt[2::3]        # area ∝ w * sigma
    weights = weights / weights.sum()
    for i, name in enumerate(names):
        w, mu, s = popt[3 * i], popt[3 * i + 1], popt[3 * i + 2]
        dmu, ds = perr[3 * i + 1], perr[3 * i + 2]
        line = (f"{name:>7}: mu {mu:.3f} +/- {dmu:.3f} eV "
                f"(calibSyst +/- {CALIB_SYST_FRAC * mu:.3f}), "
                f"sigma {s:.3f} +/- {ds:.3f} eV, area share {weights[i]:.3f}")
        if name in HATHERLY_SIGMA_EV:
            line += f"   [Hatherly prior sigma {HATHERLY_SIGMA_EV[name]}]"
        else:
            line += f"   [design E_single band {DESIGN_E_SINGLE_BAND}]"
        print(line)

    print("\nNotes: sigmas are UPPER bounds on intrinsic channel widths "
          "(instrument + Abel smoothing + calib folded in); the slow "
          "component is one effective Gaussian (sub-structure below this "
          "export's resolving power); area shares are NOT channel weights "
          "(2-D-projection-free but still detection-filtered — the OQ-C "
          "weights stay anchored on the covariance + power series).")


if __name__ == "__main__":
    main()
