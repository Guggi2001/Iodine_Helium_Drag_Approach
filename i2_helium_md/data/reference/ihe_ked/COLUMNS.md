# Column dictionary

Data dictionary for the CSV files in this folder. Conventions, error model,
and provenance details: see `README.md` and `IHe_KED_reference.provenance.json`.

---

## `IHe_KED_reference.csv` — one row per fragment n = 0…17

### Identity

| column | meaning |
|---|---|
| `n` | number of attached He atoms (I⁺Heₙ) |
| `label` | human-readable fragment name, `I+He_n` |
| `massCenter_u_per_e` | fragment mass in u (per elementary charge): 126.90 + 4.0026·n |

### Statistics

| column | meaning |
|---|---|
| `N_counts` | raw centroided counts in the TOF mass gate (whole image) |
| `N_eff` | background-subtracted counts inside the signal mask — the effective sample size behind `statErr` |

### Values (all in eV; all computed within the signal mask)

| column | meaning |
|---|---|
| `meanKE_eV` | **the reference value**: first moment ⟨E⟩ of the 3-D kinetic-energy distribution P(E), computed from the *unsmoothed* distribution. This is what an MD forward model should compare its own ⟨E⟩ against. |
| `modeKE_eV` | location of the **peak** of P(E) (mask-smoothed curve, same smoothing the signal mask uses). ⚠️ The distributions are strongly skewed — for clusters the mean lies well *above* the peak (n=1: mean ≈ 1.4× mode), for I⁺ slightly *below* it. Use this, not the mean, for "the spectrum peaks at …" statements. |
| `medianKE_eV` | median of the unsmoothed P(E) — robust central value, typically between mode and mean (the mode is grid-quantized from the smoothed curve, so at high n the median can fall slightly below it — n = 8, 11, 13, 17) |
| `sigmaKE_eV` | physical width (standard deviation) of the unsmoothed P(E). **Not an error bar** — it describes how broad the distribution is. |

### Uncertainties (of the *mean*, not of the shape)

| column | meaning | how to use |
|---|---|---|
| `statErr_meanKE_eV` | counting statistics, σ_E/√N_eff | add in quadrature per point |
| `sysErr_meanKE_eV` | analysis systematic: half-range of ⟨E⟩ over a 22-member ensemble varying every analysis choice | add in quadrature per point |
| `calibSyst_frac` | velocity-calibration scale band (fractional, E ∝ vf²) | multiply by `meanKE_eV`; **correlated** — shifts all fragments together |
| `conditionSyst_frac` | measurement-condition reproducibility (fractional), from the same-day 300 mW cross-check | multiply by `meanKE_eV`; **correlated** — shifts all fragments together |
| `bgOffShift_eV` | **diagnostic only**: how much ⟨E⟩ would move if the global-background subtraction were switched off (a rejected method). Not an uncertainty; flags where the bg choice matters (essentially only n=1). | do not fold into error bars |

**Recipe:** per-point error = √(`statErr`² + `sysErr`²); then draw
`calibSyst_frac`·⟨E⟩ and `conditionSyst_frac`·⟨E⟩ as two separate bands that
move the whole curve coherently (they do not enter point-to-point scatter).

### Flags

| column | meaning |
|---|---|
| `dominantError` | largest of {stat, analysis, calib, bg-structural} in absolute eV. `bg-structural` = the value hinges on the background choice. (`conditionSyst` is deliberately excluded — at 6% it would tag nearly every point.) |
| `noiseLimited` | 1 = treat ⟨E⟩ as an upper bound (N_eff < 100 or errors > 50% of the mean). **0 for all 18 fragments** — every value is a genuine measurement. |

---

## `IHe_KED_spectra_n0_n1.csv` — full spectra for n = 0 and n = 1

| column | meaning |
|---|---|
| `E_eV` | kinetic energy, native (uniform-in-r) grid shared by both fragments |
| `signal_IHe0` | **unsmoothed** P(E) of I⁺ (n=0), normalized so the smoothed envelope peaks at 1; `NaN` = excluded region (E < 0.01 eV, and E < 0.4 eV Abel-center spike) |
| `signal_IHe1` | same for I⁺He (n=1); `NaN` below 0.01 eV |

Note P(E)·dE is the conserved measure; on this grid dE ∝ √E, so for moments
weight each bin by P(E)·√E (equivalently: work with P(v) on the uniform-v grid).

---

## `IHe_KED_curves_n0.csv` … `IHe_KED_curves_n4.csv` — trusted curves, n = 0…4

Full reference curves for the five verified low-n fragments, one file per
fragment: both representations (2-D detector projection / 3-D reconstruction)
in both measures (per unit speed / per unit energy) on a common grid.

| column | meaning | plot against |
|---|---|---|
| `E_eV` | kinetic energy (mass-independent detector mapping, ½·127·(r·vf)²) | — |
| `v_mps` | speed of the mass-m(n) complex, √(2E/m(n)) — this axis differs per fragment even at equal E | — |
| `signal_2d_Pv` | 2-D detector projection, density **per unit speed**: angle-integrated radial profile of the background-subtracted, recentered detector image, *not* Abel-inverted (interpolated from its native 1-px grid) | `v_mps` |
| `signal_2d_PE` | the same projection as a density **per unit energy** (× dv/dE ∝ 1/v) | `E_eV` |
| `signal_3d_Pv` | the reconstructed **3-D speed distribution** P(v) (pyabel rIbeta) — the measure all reference moments use | `v_mps` |
| `signal_3d_PE` | the 3-D **kinetic-energy distribution** P(E) = P(v)·dv/dE — same convention as the spectra file | `E_eV` |

Both signal representations share this one axis pair because Abel inversion
redistributes intensity across radii without remapping the radius axis (a
shell at speed v₀ sits at the same radius in both curves: rim of the
projection = peak of the reconstruction). The *interpretation* differs:
against the `_3d_` columns the coordinate is the **true** speed |v| / energy
of the ion; against the `_2d_` columns it is the **projected in-plane**
component v·sinθ / the corresponding apparent energy (an energetic ion flying
along the TOF axis contributes to the 2-D curve near zero — the center fill).

Each signal column is normalized **independently** so its smoothed envelope
peaks at 1 in the signal region (n=0: E > 1 eV; else E > 0.01 eV) — raw bins
may exceed 1; compare shapes, never amplitudes across columns.

**Ranges.** `NaN` marks excluded regions in the **3-D columns only**:
E < 0.01 eV (n ≥ 1) or E < 0.4 eV (n = 0 — the Abel-center spike, an
inversion artifact). The **2-D columns cover the full detector range**: the
raw projection is clean down to r = 0 (verified for n=0; occasional small
negatives are background-subtraction noise). The upper end of every curve is
the detector/crop edge, 354 px ≈ 6.13 eV ≈ 3052 m/s at m = 127. The
on-the-day lab-book plots show this same edge at ≈3576 m/s because they used
the legacy velocity factor 10.0995 — same pixels, superseded calibration.

⚠️ A `_Pv` column plotted against `E_eV` is only an axis relabel, *not* the
energy density — peak positions and apparent widths differ between the `_Pv`
and `_PE` versions of the same data (Jacobian).

⚠️ `signal_2d_PE` is elevated at low E: per unit energy the projection profile
is ∝ the ring-averaged image *brightness*, and the center of a VMI image is
filled by fast ions projected inward (every shell covers the whole disk). It
tends to a **finite center-fill plateau** at E → 0 — do not read the low-E
elevation as slow-ion signal. The 3-D columns are immune (P(v) ∝ v² near 0).

---

## `IHe_KED_crosscheck_300mW.csv` — independent 300 mW series vs. the reference

| column | meaning |
|---|---|
| `n`, `fn_300mW` | fragment and its 300 mW filenumber (43614–43631) |
| `N_eff_300mW` | effective counts in the 300 mW measurement |
| `meanKE_300mW_eV`, `statErr_300mW_eV` | ⟨E⟩ of the 300 mW series, same pipeline (central config only) |
| `bgOffShift_300mW_eV` | bg-off diagnostic for the 300 mW series |
| `meanKE_ref600_eV`, `refErr_stat_sys_eV` | the reference value and its per-point error (stat ⊕ sys) |
| `diff_eV`, `diff_frac` | 300 mW minus reference, absolute and fractional |
| `z_indicative` | diff / combined error. *Indicative only*: analysis systematics are correlated between the series, and the calibration band cancels (same camera). |
