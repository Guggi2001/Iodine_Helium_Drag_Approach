# Drag_function

Extracts a closed-form **drag law** `F_drag(v)` for an iodine ion moving through
a helium nanodroplet, using TDDFT reference trajectories of I₂⁺ Coulomb
explosions as the source of truth. The resulting law is intended to replace the
hard-sphere collision model in the sibling MD code
(`../i2_helium_md/physics/collisions.py`).

## Physical idea

For each iodine atom the equation of motion in the droplet is

```
m_eff · a(t) = F_Coulomb(R(t)) + F_drag(v(t))
```

so the drag force is obtained by **force balance**:

```
F_drag(t) = m_eff · a(t) − F_Coulomb(R(t))
```

- `R(t)` is the instantaneous internuclear distance (TDDFT input).
- `F_Coulomb = q₁q₂ / R²` (in `amu·Å/ps²` after unit conversion).
- `m_eff = M_EFF_AMU ≈ 179.912 amu` is the effective mass of an iodine atom
  dressed by its co-moving helium shell (set in `drag_calculation.py`).

Because both iodines in the dimer feel the *same* Coulomb repulsion, the
Coulomb term is known analytically from `R(t)`; whatever remains in
`m_eff·a − F_C` is, by construction, helium-induced drag.

### Why a spline for `a(t)`?

Numerically differentiating a noisy `v(t)` amplifies noise by `~1/Δt` per
derivative — `a(t) = dv/dt` would be unusable. Instead `a(t)` is taken as the
**analytic second derivative of a cubic smoothing spline** (`spline_k=3`) fit
to a denoised `v(t)`. The spline's smoothing parameter `s` controls how closely
it follows the data (see "Parameter sensitivity" below).

A **trusted-interior mask** drops the first/last `truncate_points = 500`
samples (~0.5 ps) before fitting. Spline boundary derivatives are unreliable,
and the early window is dominated by the bubble-exit transient — both ranges
would otherwise bias the drag fit.

## Smoothing the reference velocity

The raw TDDFT `v(t)` carries two noise sources: a high-frequency numerical
component, and a 1.2 ps quasi-periodic oscillation from the helium-bubble
breathing mode. Three smoothing tiers are available; they can be used
individually or stacked.

### 1. Savitzky–Golay (`fourier_analysis.sg_smooth_v`, l.22)

```python
v_s, resid, dt = sg_smooth_v(t, v, window_length, polyorder=3)
```

`scipy.signal.savgol_filter` with `mode="interp"`. Default `polyorder=3` in the
signature; **production runs use `polyorder=1`** with these window lengths
(samples, dt ≈ 1 fs):

| case | window_length | t-window where SG is applied |
| ---  | ---           | ---                          |
| 9 Å  | 3901          | [2.67, 8.5] ps |
| 18 Å | 3401          | [4.543, 8.0] ps |

Returned `resid = v − v_s` is what feeds the FFT diagnostic below.

### 2. CEEMDAN (`smoothing_ceemdam.CeemdanConfig`, l.86)

```python
@dataclass
class CeemdanConfig:
    trials: int = 50
    noise_width: float = 0.15
    random_seed: int = 0
    max_imfs: Optional[int] = None
```

Production overrides: `trials=200`, `noise_width=0.2`. After decomposition,
IMFs are filtered by their dominant period:

- `target_period_ps = 1.2` (the bubble mode) with `rel_tol = 0.25`
  → drop any IMF with mean period in **[0.9, 1.5] ps**.
- `allow_second_imf = True`, `second_imf_amp_ratio = 0.35` → also drop the
  *next* oscillatory IMF if it carries ≥35 % of the matched IMF's amplitude
  (handles split modes).

### 3. Hybrid CEEMDAN + SG (`smoothing_pipeline.py`)

`CeemdanParams(trials=200, noise_width=0.2)` peels the bubble IMF off first,
then `SgParams(window_length=1401, polyorder=1)` smooths the residual.
`ModeSelectionParams` (l.120-141) carries the period-matching knobs above.

### FFT diagnostic that justifies the cut

`fourier_analysis.dominant_periods_fft` (l.118-177) reports the top
`n_peaks = 3` spectral peaks above `fmin = 0.05` with
`peak_min_prominence_ratio = 0.05 · max(power)`. The 1.2 ps target period was
read off this plot for the current dataset — **regenerate this diagnostic
before re-using the IMF-drop rule on a different droplet size or charge
state.** The `if __name__ == "__main__":` block in `fourier_analysis.py`
(l.412-423) sweeps `wls = [1001, 2501, 3401]` and overlays the residual
spectra so the SG window can be picked on the same plot.

## Choosing the drag-extraction window `[t*, t_end]`

The bubble-exit transient (≲ few ps after ionisation) violates the
single-iodine force-balance assumption, so the fit window must start after the
ion has stabilised in the surrounding helium. `transient_exclusion.find_t_star_stationary_residual`
(l.35-152) finds that crossover automatically:

1. SG-detrend `v(t)` with `trend_win_ps = 0.8` (polyorder 2). Subtract the
   trend → residual `r(t)`.
2. Compute rolling-RMS of `r(t)` over `rms_win_ps = 0.5`.
3. In the late-time `late_window` (where the system is known to be stationary),
   build a robust baseline band:
   `median(rRMS) ± k · MAD`, `k = 2.5`, `MAD = 1.4826 · median(|x − median|)`.
4. Walk forward from `t = 0`. **`t*` is the earliest time at which the
   rolling RMS stays inside the band continuously for `sustain_ps = 0.8` ps**,
   capped at `t_out`.

Per-case wrappers in `transient_exclusion.py:157-182`:

| case | `t_out` (ps) | `late_window` (ps) | extracted window |
| ---  | ---          | ---                | ---              |
| 9 Å  | 9.0          | (5.0, 8.5)         | ~[2.67, 8.5] ps |
| 18 Å | 8.5          | (5.0, 8.0)         | ~[4.54, 8.0] ps |

The function also returns a diagnostics dict (`late_med`, `late_mad`, `band_lo`,
`band_hi`, `trend_win`, `rms_win`) for plotting the band over the residual.

## Parameter sensitivity / sweeps

Production fits are tested for stability under both smoothing knobs *and*
window choice. The formal sweep is `RobustnessSweep`
(`smoothing_pipeline.py:299-310`):

| knob          | sweep values                  | what it tests                          |
| ---           | ---                           | ---                                    |
| `noise_widths`| `(0.15, 0.25)`                | CEEMDAN noise-width robustness         |
| `sg_windows`  | `(1401,)`                     | placeholder — expand for new datasets  |
| `trials_list` | `(100, 200, 400)`             | CEEMDAN ensemble convergence           |
| `seeds`       | `range(10)`                   | 10 independent realisations → 1σ band  |
| `t_min_list`  | `(2.5, 2.7, 3.0, 3.5, 4.0)`   | sensitivity to `t*`                    |
| `t_max_list`  | `(8.0, 8.5, 9.0)`             | sensitivity to `t_end`                 |

The 10-seed sweep is what produces the uncertainty band on the fitted drag
coefficients. The `t_min/t_max` lists are the way to overrule the automatic
`t*` and confirm the result is not an artifact of one specific window. Each
module has an `if __name__ == "__main__":` block that runs these sweeps and
plots them side-by-side — that block is the intended re-tuning entry point.

### Spline-`s` heuristic

In `drag_calculation.py:253-271`, if `DragExtractionSettings.spline_s is None`:

```
σ = 1.4826 · MAD(Δv)        # robust noise scale per sample, in Å/ps
s = N · σ²                  # total smoothing budget for UnivariateSpline
```

with a fallback `s = 1e-12` if `σ ≈ 0` (already smooth). Override
`DragExtractionSettings(spline_s=...)` to fix `s` when the heuristic is too
tight (wiggly `a(t)`) or too loose (flat `F_drag`).

## Fit variants and extracted results

Both fits are selected via `DragExtractionSettings.fit_variant`. Defaults
(`drag_calculation.py:194-214`):

```python
@dataclass(frozen=True)
class DragExtractionSettings:
    q1: int = 1
    q2: int = 1
    meff_amu: float = M_EFF_AMU         # 179.912 amu
    spline_k: int = 3
    spline_s: Optional[float] = None    # auto if None (see above)
    truncate_points: int = 500
    truncate_time_ps: Optional[float] = None
    show_plots: bool = True
    fit_power_law: bool = True
    fit_variant: int = 1                # 1 = power law (default), 2 = linear+cubic
    case: int = None                    # 9 or 18, for plot titles
    plot_only_drag_with_fit: bool = False
```

### Variant 1 — power law `|F_drag| = γ · v^n` (l.485-529)

- Solved as **log-linear least squares**: `log|F| = log γ + n · log v` via
  `np.linalg.lstsq` over the positive `(v, |F|)` pairs on the trusted interior.
  Signs are restored only when plotting.
- No explicit `p0` / bounds — log-linear regression is closed-form.
- Returns `fit_result = {"gamma": γ, "n": n}`.
- R² is printed on the scatter overlay but **not** stored in `fit_result`.

This is the default (`fit_variant = 1`) and is the form intended to drop into
`../i2_helium_md/physics/collisions.py`.

### Variant 2 — linear + cubic `F_drag = a · v + b · v³` (l.536-597)

- Non-linear fit:
  `popt, pcov = curve_fit(drag_model, v_fit, F_fit, p0=p0, maxfev=20000)`.
- Initial guesses:
  - `a₀` = slope of the linear fit on the slowest ~10 % of sorted `v`.
  - `b₀ = median((F − a₀·v) / v³)` on the fastest 10 % of `v`.
  - No explicit bounds — relies on the LM optimiser.
- 1σ uncertainties from `sqrt(diag(pcov))`.
- R² is computed *and returned* for this variant:
  `R² = 1 − Σ(F − F_pred)² / Σ(F − mean(F))²`.
- Returns `fit_result = {"a": a, "b": b, "a_err": a_err, "b_err": b_err}`.
  Units: `a` in amu / ps, `b` in amu · ps / Å².

## Data inputs

Local TDDFT data lives outside the repo (paths in
`config_utils_local/config.py`):

```
…/Masterarbeit/Drag_Calculation/Data_DFT/{9A,18A}/
    data_vabs.csv      # |v| of iodine 1
    data_vabs2.csv     # |v| of iodine 2
    R1-R2.csv          # internuclear distance
    9A_All_Data.csv    # pre-processed bundle used by io.load_data
    18A_All_Data.csv
```

Pre-processed CSV columns (units: Å, ps, Å/ps):
`Time_ps, V1_mag, V2_mag, V1_z, V2_z, V1_x, V2_x, R_distance`.

The two cases (`9A`, `18A`) correspond to two helium-bubble sizes around the
I₂⁺ at the moment of ionisation.

## Module map

```
drag_function/
├── reference_data.py     # load raw TDDFT CSVs, build cleaned dataset, sanity plots
├── io.py                 # load pre-processed CSVs (load_data), 2×2 reference plots
├── fourier_analysis.py   # Savitzky–Golay smoothing + FFT diagnostics
├── smoothing_ceemdam.py  # CEEMDAN IMF / period utilities (PyEMD)
├── smoothing_pipeline.py # CEEMDAN+SG hybrid denoise w/ sensitivity sweeps
├── transient_exclusion.py# find_t_star_stationary_residual() — skips bubble-exit transient
└── drag_calculation.py   # main_calculation(): spline a(t), F_C, force balance, fit
config_utils_local/
└── config.py             # PATH9A / PATH18A — local data locations
cleaned_data.csv          # exported pre-processed dataset
```

## Pipeline

```
raw TDDFT CSV
    │  reference_data.py / io.load_data
    ▼
v(t), R(t) per atom
    │  fourier_analysis.sg_smooth_v  /  smoothing_pipeline (CEEMDAN+SG)
    ▼
denoised v(t)
    │  transient_exclusion.find_t_star_stationary_residual
    ▼
window [t*, t_end]
    │  drag_calculation.main_calculation:
    │    cubic smoothing spline → a(t)
    │    F_C = q1q2/R²
    │    F_drag = m_eff·a − F_C
    │    curve_fit on |F_drag| vs v  (variant 1 or 2)
    ▼
fit_result dict + diagnostic plots
```

## Output

`main_calculation(t_ps, R_A, v_Aps, settings)` returns a dict with:

- `t_ps`, `R_A`, `v_data_Aps`, `v_spline_Aps`, `a_spline_Aps2`
- `F_C_amuAps2`, `F_inert_amuAps2`, `F_drag_amuAps2`
- `trusted_mask` — boolean interior selector
- `fit_result` — `{γ, n}` (variant 1) or `{a, b, a_err, b_err}` (variant 2)
- `spline_s` — smoothing parameter actually used

## How to run

```python
from drag_function.io import load_data
from drag_function.fourier_analysis import sg_smooth_v
from drag_function.drag_calculation import main_calculation, DragExtractionSettings
from config_utils_local.config import PATH9A

d = load_data(PATH9A)
t, v, R = d["t"], d["v1"], d["R"]
v_sg, _, _ = sg_smooth_v(t, v, window_length=3901, polyorder=1)

out = main_calculation(
    t, R, v_sg,
    DragExtractionSettings(case=9, truncate_points=500, fit_variant=1),
)
print(out["fit_result"])
```

Each module also contains an `if test:` / `if __name__ == "__main__":`
block that runs interactive diagnostic plots — useful for re-tuning the
smoother on a new dataset.

## Status

- No `pyproject.toml`, `requirements.txt`, or pytest suite in this folder.
- Runtime deps: `numpy`, `scipy`, `pandas`, `matplotlib`, `EMD-signal`.
- Data paths are hard-coded for the author's machines (`home = True/False`
  switch in `config_utils_local/config.py`).
- Tests for the eventual MD integration live in `../i2_helium_md/tests/`.

## Related

- `../i2_helium_md/PHYSICS_BASELINE.md` — frozen baseline of the MD code
  that will consume the drag law produced here.
