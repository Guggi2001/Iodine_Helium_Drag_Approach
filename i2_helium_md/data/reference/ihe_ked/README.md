# I⁺Heₙ kinetic-energy reference (experiment → MD validation)

Per-fragment kinetic-energy distributions of I⁺Heₙ (n = 0…17), reduced to a
mean energy with a full error budget, for comparison against the molecular-
dynamics simulation (H.2b forward model).

**Source of truth:** the MATLAB analysis path in
`matfile_data_scripts/A_state_paper_figures_single_pulse/fragment scan/171024/`.
Exporter: `export_IHe_KED_reference.m`. Measurement date 17.10.24, FLIR camera.

---

## Files in this folder

| file | contents |
|---|---|
| `IHe_KED_reference.csv` | the per-fragment table (one row per n) |
| `IHe_KED_spectra_n0_n1.csv` | UNSMOOTHED normalized P(E) for n=0 and n=1 (common native E axis; NaN marks excluded regions: E < 0.01 eV for both, E < 0.4 eV for n=0 only) |
| `IHe_KED_curves_n0.csv` … `n4.csv` | trusted curves for n = 0…4: 2-D detector projection AND 3-D reconstruction, each in BOTH measures (per unit v and per unit E), peak-normalized, dual E/v axes (see `COLUMNS.md`) |
| `IHe_KED_reference.provenance.json` | machine-readable record of every convention below |
| `COLUMNS.md` | data dictionary — what every column in every CSV means and how to use it |
| `IHe_KED_crosscheck_300mW.csv` | per-fragment ⟨E⟩ of the independent 300 mW series vs. the reference (same pipeline; see cross-check section) |
| `IHe_KED_crosscheck_300mW.png` | the comparison figure for the row above |
| `Experimental_details/New Section 1.one` | lab-book pages (OneNote) for 17.10.24: measurement tables, TOF gates, on-the-day plots |
| `README.md` | this file |

To deploy: copy the CSVs, the provenance JSON, `COLUMNS.md`, and this README
into the MD repo's `data/reference/ihe_ked/` (the crosscheck PNG and
`Experimental_details/` stay here).

---

## What the numbers mean (conventions)

- **Measurement series (lab book, 17.10.24):** the exported fragments are the
  day's final production series at **600 mW probe power**: 43634 (I⁺),
  43636 (IHe⁺), 43640…43654 (IHe₃⁺…IHe₁₇⁺). The first IHe₂⁺ take (43637) was
  contaminated by an accidental pump-probe measurement and **remeasured as
  43656** — that is why the filenumber list looks irregular. 43655 is a noise
  measurement (MCP front 0 kV) taken because 600 mW shows a "substantial
  background for all masses"; the lab book prescribes subtracting it, which
  this export does (scaled per fragment, see Background below). An independent
  **300 mW series** (43614–43631, lower intrinsic high-E background, ~half
  signal) exists as a cross-check — see below.
- **Quantity:** lab-frame translational KE of the *detected* complex,
  ½·m(n)·v², per detected ion (NOT kinetic-energy release, NOT pair-CM energy).
- **Mass model:** `m(n) = 126.90 + 4.0026·n` u (matches the MD mass gates).
- **Density:** pyabel's (rbasex, reg=pos) speed distribution `I(r)` from
  `distr.rIbeta()` — the properly v²- and sin θ-weighted **3-D KED**, not the
  2-D slice density and not the projection.
- **Moments:** the grid is uniform in r, so dv is constant and
  `⟨E⟩ = Σ E·P(v) / Σ P(v)` **is** the proper first moment
  `∫E·P(E)dE / ∫P(E)dE`. This is what an MD forward model should compute from
  its own P(E) — no pipeline quirks to replicate.
- **Signal mask:** the *contiguous* supra-threshold region of the
  movmean-smoothed ('omitnan'), peak-normalized P(E) containing the peak
  (threshold 0.05 central; low-E cut 0.01 eV; n=0 additionally E > 0.4 eV to
  exclude the Abel center spike). Isolated supra-threshold islands at high E
  are residual-background artifacts with large leverage under proper dE
  weighting and are excluded. Smoothing is used ONLY for the mask —
  **moments and sigmaKE come from the unsmoothed distribution.**
- **Mass-independent KE:** for singly-charged VMI the detector radius maps to
  energy independently of mass, `E = ½·127·(r·vf)²`. ⟨E⟩ is computed on each
  fragment's **native** radial grid (no cross-fragment resampling; the n=0/n=1
  grid identity for the spectra file is asserted, not assumed).
- **Frame:** droplet rest frame. The VMI image is centered on the beam spot
  (`VM_center` + per-fragment recentering), i.e. the droplet-beam velocity is
  removed, so ⟨E⟩ agrees with the MD droplet-rest-frame speed to first order.
- **Background:** the global background (measurement 43655) is scaled per-cluster
  over the signal-free 4–8 eV window and subtracted; the bare I⁺ (n=0) uses its
  own local background (fn+1) instead. (The window is signal-free for every
  fragment — the n=1 tail ends ≳6σ below 4 eV.)
- **Calibration:** `velocity_factor = 8.6178` — the **FLIR + 7.95 kV** value.
  ⚠️ NOT 10.0995 (that is the IDS-camera value; using it inflates every energy
  by (10.0995/8.6178)² ≈ 1.37). Rule of thumb: filenumbers > ~2100 are FLIR.

### ⚠️ Not comparable to the published fragment figure

These values are **not** on the same footing as the paper figure's ⟨E⟩ axis,
for two independent reasons: (1) the paper pipeline computed moments of the
2-D slice density without the dE Jacobian (a legacy convention that biases
⟨E⟩ low by ~25–30%), and (2) the published axis appears to use the IDS
velocity factor (+37% if applied to this FLIR data). For scale: ⟨E⟩(I⁺) is
**3.71 eV** here, ~2.81 eV under the legacy measure at the FLIR factor, and
~2.5 eV in the published figure. For MD validation use THIS file; reconcile
the paper separately.

---

## Error model

Five contributions, deliberately kept separate because they behave differently:

| column | meaning | size |
|---|---|---|
| `statErr_meanKE_eV` | counting error, σ_E/√N_eff | ~0.1–1.1% (statistics never limit) |
| `sysErr_meanKE_eV` | **analysis systematic** — half-range of ⟨E⟩ over a 22-member ensemble that varies *every* analysis choice (bg window, threshold incl. the paper's κ-tuned 0.01, low-E cut, I⁺ Abel-spike cut, smoothing) | 0.9–15.5%, trend grows with n (not monotone) |
| `calibSyst_frac` | fractional **coherent scale** band (E ∝ v²) | 0.04 (see caveat) |
| `conditionSyst_frac` | fractional **measurement-condition reproducibility** band, from the 300 mW cross-check (coherent +5–7% cluster shift, −6.3% for n=0 between the two same-day series) | 0.06 |
| `bgOffShift_eV` | **diagnostic only** — shift if the global-bg subtraction is turned off | large only for n=1 (~29%); elsewhere \|shift\| ≤ 8% (−8.0% at n=17, < 5% for all other n) |

**Shape ≠ error.** The distributions are strongly skewed, so the mean is *not*
where P(E) peaks: for the clusters the mean lies well above the peak (n=1:
mean ≈ 1.4× mode, skewness ≈ +1.8), for I⁺ slightly below it. This is a
property of the distribution, not an uncertainty — the CSV carries it as data
in the `modeKE_eV` and `medianKE_eV` columns. All error columns refer to the
**mean**. Compare MD moment-to-moment (mean to mean), never mean-to-peak.

**Combine as:** per-point error = `√(statErr² + sysErr²)`, then apply
`calibSyst_frac · meanKE` and `conditionSyst_frac · meanKE` as **separate
correlated bands** (each shifts all fragments together, so neither enters the
point-to-point scatter). `conditionSyst_frac` is not folded into
`dominantError` — at 6% it would tag nearly every point and hide the internal
error structure the tags are meant to show.

`bgOffShift` is a *rejected-method* sensitivity (bg-off leaves cluster tails
in, and under the proper v²-weighted measure those tails carry even more
leverage), so it is reported but **not** folded into `sysErr`. It flags where
the background choice matters — essentially only n=1.

`dominantError` names the largest of {stat, analysis, calib, bg-structural} per
row. `noiseLimited` = 1 would mean treat ⟨E⟩ as an upper bound — it is **0 for
all 18 fragments** (N_eff 3275–194000), so every value is a genuine measurement,
including the sub-0.1-eV high-n points.

---

## Trust guidance for the comparison

- **Gold points (calib-limited): n = 0, 3, 4, 5, 6, 7, 8, 9, 10, 12.** Internal
  analysis systematic is below the calibration floor → a disagreement here is
  real physics — *beyond* the two correlated bands (4% calib, 6% condition),
  which shift the whole curve, not individual points.
- **I⁺ scale anchor (n = 0):** 3.706 eV, good to ~1% ⊕ 4% ⊕ 6%
  (analysis ⊕ calib ⊕ condition; the last from the −6.3% shift at 300 mW).
- **n = 1:** background-structural (bg-off shifts it +29%); report with the bg
  caveat, not as a clean point.
- **n = 2 and n = 16:** nominally bg-structural, but all error sources are
  ~4% there — treat as ordinary points with a mild bg caveat.
- **Analysis-limited (n = 11, 13, 14, 15, 17):** 5–15.5%; weight accordingly.
- Statistics limit nothing here.

---

## Cross-check: independent 300 mW series

The same day's **300 mW** series (43614–43631, ~half the count rate, much lower
intrinsic high-E background) was run through the *identical* pipeline
(`crosscheck_300mW_series.m` → `IHe_KED_crosscheck_300mW.csv/.png`).
The calibration band cancels (same camera/vf); the analysis systematics are
correlated between the series, so the z column is indicative.

- **Decay shape confirmed at the few-% level for all n.** No fragment deviates
  by more than ~8% (max +8.0% at n=14); nothing resembling the 2–3× failures
  the old pipeline had.
- **n = 1 validated.** The reference's weakest point (bg-structural, +29%
  bg-off shift) reproduces to **+0.02%** at 300 mW, where the bg-off
  sensitivity drops to +8.7%. The n=1 value does not, in fact, hinge on the
  background treatment as much as the 600 mW diagnostic alone suggests.
- **A coherent inter-series offset of ~+5–7% for n ≥ 2** (300 mW above the
  reference; z up to ~4 for the tightest gold points). This is a genuine
  measurement-condition systematic (probe power and/or the 600 mW background
  regime), *not* captured by the analysis-choice ensemble. It is carried in
  the reference CSV as the `conditionSyst_frac = 0.06` correlated band (see
  Error model above).
- **n = 0 is the outlier: −6.3%** (300 mW *below* the reference, formally
  z ≈ −6 against its ~1% analysis error). I⁺ ran at ~81 counts/frame at
  600 mW vs ~11 at 300 mW, so saturation/space-charge or true power dependence
  of the bare-I⁺ channel are plausible; the sign is opposite to the cluster
  offset. The I⁺ scale anchor is therefore good to ~1% ⊕ 4% (calib) *within*
  the 600 mW condition, but carries a ~6% condition sensitivity.

---

## Regenerating

```
matlab -batch "run('export_IHe_KED_reference.m')"   % ~10 min (5 Abel passes)
```
Runs headless; writes all files in this folder. The script pins its own
`velocity_factor`, `blocksize_default`, `abel_inv_post`, and `abel_inv_method`,
validates that cached `VMIdata_*.mat` match the pinned blocksize AND share one
grid (blocksize can match while the crop/VM_center vintage differs — that
mixed-cache state crashes `subtract_processed_data`), reprocessing any stale
ones, and asserts the n=0/n=1 grid identity before writing the spectra file —
so it is independent of MATLAB session state and of caches left by other
figure scripts. The cross-check regenerates the same way:

```
matlab -batch "run('crosscheck_300mW_series.m')"    % ~4 min (2 Abel passes)
```

---

## Open items / caveats

1. **Calibration uncertainty is a placeholder.** `calibSyst_frac = 0.04` assumes
   2% on `velocity_factor`. Replace with the real FLIR CEI-fit precision (+ the
   2031.8 m/s reference-velocity model error) by editing `u_vf_frac` and re-running.
   When it tightens, the "calib"-tagged points tighten together and some may flip
   to "analysis".
2. **Mass gates: centers recovered, widths still unknown.** The TOF gate delay
   (D3) centers are absent from the data info files but WERE found in the lab
   book (`Experimental_details/`) and are now in the provenance JSON
   (`mass_gate.D3_gate_center_ns`): 8530 (I⁺), 8650, 8730, 8800, 8880, 8953,
   9030, 9105, 9175, 9248, 9317, 9386, 9458, 9527, 9595, 9662, 9729, 9790 ns
   (I⁺He₁₇). They sit on the observed I⁺Heₙ TOF peaks with √m-consistent
   spacing → unit-mass selection confirmed. Gate *widths* were not logged;
   `massWindow*` columns therefore remain omitted.
3. **Paper reconciliation (independent of this export).** See the "not
   comparable" box above: the published fragment figure differs from this
   reference by the legacy moment convention AND (apparently) the IDS velocity
   factor. Worth reconciling on the paper side.
4. **κ of ⟨E⟩(n) not re-fit.** The earlier consistency check (κ ≈ 2.5 vs paper
   2.27, ±0.4) was done under the legacy measure; if κ is needed it must be
   re-fit to these corrected means.

---

## Changelog — 2026-07-14 deployment audit (docs only, values unchanged)

A dedicated review of the MD-repo deployment (`i2_helium_md/data/reference/
ihe_ked/`): byte-parity, consumer-side code review, and numeric re-validation
of every documented claim against the shipped data. **Data verified clean** —
all deployed files byte-identical to this folder, reference means reproduce
from the shipped spectra to 0.00% (n=0: 3.7057, n=1: 1.3017 eV),
`dominantError` matches the recomputed argmax for all 18 rows, curve axes
exactly E = ½·m(n)·v². Fixes applied:

- The vmi_summary byte-identical note (promised in the changelog below) was
  actually missing from the MD repo — now written (`vmi_summary/README.md` +
  the data/reference README known-issues list).
- This README and `COLUMNS.md` are now part of the deployment (the deployed
  `COLUMNS.md` referenced a README that wasn't shipped, leaving the trust
  tiers unreachable from the MD repo).
- Provenance `measurement_ids.fragments` reordered to ascending n = 0…17 with
  an explicit `fragments_order` field (was descending — an index-zip trap
  against the ascending `D3_gate_center_ns` array); exporter updated to match.
- Doc-precision corrections found by re-measuring every claim: bgOffShift
  bound is "≤8% elsewhere (−8.0% at n=17)", not "<5%"; median lies between
  mode and mean only *typically* (grid-quantized smoothed mode; n=8,11,13,17
  violate); crop edge corrected 352 px/6.05 eV/3030 m/s → 354 px/6.13 eV/
  3052 m/s everywhere incl. the exporter; crosscheck "≤8%" annotated with its
  max (+8.0% at n=14); sysErr "grows with n" is a trend, not monotone.
- Overview plot script: P(v)-panel dashed marker now labeled v(⟨E⟩) (it is
  not ⟨v⟩); display smoothing no longer bleeds into the NaN-cut regions; PNG
  regenerated and gitignored in the MD repo (generated-figure policy).

---

## Changelog — 2026-07-14 lab-book audit, validation & trusted-curve release

A full verification-and-hardening pass. **Reference ⟨E⟩ values unchanged
throughout** — everything below is validation, error modeling, shape
annotation, and new exported artifacts.

### Lab-book audit (data selection verified)

Recovered the 17.10.24 lab book (`Experimental_details/`, binary OneNote) and
cross-checked every filenumber: the export uses the day's final 600 mW
production series, correctly swapping the contaminated IHe₂⁺ take 43637
(accidental pump-probe; 43638/43639 were the pump-/probe-only checks) for its
remeasure 43656, and subtracting the 43655 noise measurement (MCP front 0 kV)
exactly as the lab book prescribes. Added to provenance: probe power (600 mW),
series chronology, effusive references, D3 mass-gate centers (closing open
item 2). Also recovered: the same-day 300 mW series 43614–43631 as an
independent cross-check candidate.

### 300 mW cross-check (values validated)

Ran the identical pipeline on the independent 300 mW series
(`crosscheck_300mW_series.m` → CSV/PNG, section above). Decay shape confirmed
≤8% everywhere (max +8.0% at n=14); n=1 reproduces to +0.02% (its bg-off sensitivity drops
29% → 8.7% at 300 mW, softening the bg-structural concern); a coherent
+5–7% inter-series offset for n ≥ 2 (n=0: −6.3%) is now carried as the
`conditionSyst_frac = 0.06` correlated-band column.

### Abel inversion validated (I⁺ shape is real)

The broad 2-D I⁺ projection → rim-peaked 3-D KED (peak ≈ 4.9 eV) was verified
to be correct, not an inversion artifact: (1) a manual shell integral
r²⟨F|sinθ|⟩ of the inverted slice peaks at the same 2740 m/s as pyabel's
rIbeta; (2) synthetic round-trips through the identical method
(rbasex reg=pos) recover a broad *mid*-peaked P(v) with no rim artifact
(1.7% rms) and an I⁺-like rim-peaked P(v) exactly (1.5%), whose forward
projection reproduces the observed broad plateau.

### Shape annotation (mean ≠ peak)

The distributions are strongly skewed (n=1: mean ≈ 1.46× mode; I⁺ skewed the
other way) → new `modeKE_eV`/`medianKE_eV` columns, `shape_note` in the
provenance, and the "Shape ≠ error" paragraph in the Error model. Errors
refer to the mean; compare MD mean-to-mean, never mean-to-peak.

### Trusted curves for n = 0…4

New `IHe_KED_curves_n0..n4.csv`: both representations (2-D detector
projection / 3-D reconstruction) in both measures (per unit v / per unit E)
on one shared axis pair (`E_eV`, `v_mps`; v is per-fragment mass-corrected).
Axis semantics documented (2-D: projected in-plane speed / apparent energy;
3-D: true speed/energy). Low-E NaN cuts apply to the **3-D columns only**
(Abel-center spike = inversion artifact); the 2-D columns ship the full
detector range after verifying the raw projection is clean to r = 0.
`signal_2d_PE` tends to a finite center-fill *plateau* at low E (projected
fast ions — not slow-ion signal; an interim note wrongly called this a 1/√E
divergence and was corrected). The upper end of every curve is the crop edge
354 px ≈ 6.13 eV ≈ 3052 m/s; the on-the-day lab-book plots show that same
edge at ≈ 3576 m/s because they used the legacy velocity factor 10.0995.

### Infrastructure & deployment

Second cache guard in the exporter (shared-grid check — blocksize alone
missed mixed crop/center cache vintages, which crash
`subtract_processed_data`). New `COLUMNS.md` data dictionary. Everything
deployed to the MD repo at `i2_helium_md/data/reference/ihe_ked/` (indexed in
its README, with an energy-space units-exception note); overview plot
generator `data/reference/scripts/experimental_reference_fragments.py` →
`ihe_ked/experimental_reference_fragments.png`. Also noted there: the two
`vmi_summary` I⁺He CSVs are byte-identical (the "high-SNR" MAT is the same
4-run average, ~23k counts total vs 63k in this reference) — so the
`vmi_summary` I⁺He curve is a single different-campaign (45xxx) measurement,
whose core agrees with this reference while its ⟨E⟩ sits ~16% lower
(condition-dependent tail weight).

---

## Changelog — 2026-07 correctness overhaul

A multi-agent code review of the export (2026-07-13) found systematic defects;
all were fixed and the reference regenerated (2026-07-14). Summary of what
changed and why, so future readers can interpret diffs against older copies:

### The central bug: the moments were not moments of P(E)

The legacy pipeline (inherited from the paper scripts) computed
`⟨E⟩ = Σ E·y / Σ y` with `y = radial_distribution / v` on the r-uniform grid.
Two stacked errors:

1. `radial_distribution` (the angle-summed polar image, one r Jacobian) is the
   **2-D slice radial density** ∝ v·⟨F⟩(v), NOT the speed distribution
   ρ(v) ∝ v²·⟨F⟩(v). The textbook relation ρ(E) = ρ(v)·dv/dE ∝ ρ(v)/√E is
   correct, but it was applied to the wrong function — one factor of v short
   (plus missing sin θ weighting for anisotropic channels).
2. The moment sum omitted the bin-width Jacobian dE ∝ √E of the r-uniform grid.

Verified empirically against the data itself: the proper shell integral of the
Abel-inverted slice (r²·Σ F·|sin θ|) matches pyabel's `distr.rIbeta()` speed
distribution to 7–9% rms, while `radial_distribution` is off by 36–57% with
log-slope ≈ −1 (the missing factor of v). **Fix:** the export now uses pyabel's
`I(r)` directly (propagated through `abel_invert_processed_VMI` as
`speed_r`/`speed_distribution`); on the uniform-r grid dv is constant, so
`⟨E⟩ = Σ E·P(v) / Σ P(v)` is the exact first moment of P(E). Net effect: all
means shifted coherently UP ~25–32% (n=0: 2.81 → 3.71 eV); κ/decay shape barely
affected (which is why the paper never noticed).

### Secondary bugs fixed

- `movmean` without `'omitnan'` silently widened the low-E cut by ~5 bins
  (effective 0.0196 eV instead of the documented 0.01 eV) and clipped real
  high-n signal at 0.011–0.019 eV.
- The n=0 Abel-spike cut zeroed bins BEFORE smoothing, leaking spike signal
  back below 0.4 eV and into the mean. Cuts are now masks applied to an
  unsmoothed copy.
- `sigmaKE` was computed from the smoothed distribution (boxcar-inflated
  ~6% at high n, propagating into statErr). Now from the unsmoothed one;
  smoothing is used ONLY to build the signal mask.
- The bg-scaling window [2 6] eV overlapped real n=1 signal → common-mode
  over-subtraction invisible to sysErr. Now [4 8] eV (signal-free, matching
  the paper's κ-sweep choice) with signal-free ensemble variants.
- The systematic ensemble excluded the paper's κ-tuned threshold 0.01; the
  ensemble (now 22 members) includes it, so sysErr brackets the published
  analysis.
- **Contiguous signal mask:** under the corrected v²-weighted measure,
  isolated supra-threshold noise islands at high E acquire huge leverage
  (diagnosed: thr=0.01 admitted 2–4 islands up to ~4 eV and inflated high-n
  means 2–3×). The mask is now the contiguous supra-threshold region
  containing the peak — which is what "signal region" always meant.

### Reproducibility guards added

Pinned `abel_inv_method='pyabel'`; cached VMIdata validated against the pinned
blocksize (stale caches reprocessed); n=0/n=1 grid identity asserted before
writing the spectra CSV; clear errors when output CSVs are locked (e.g. open
in Excel). Also added: `plot_velocities` flag in the exporter (P(v) figure for
n = 0…4 with dual velocity/energy axes and ⟨E⟩ markers →
`IHe_velocity_distributions.png`).

### Related fix in the MD repo

The same slice-vs-ρ(v) bug affected the `vmi_summary` references in
`Iodine_Helium_Drag_Approach/i2_helium_md/data/reference/` (Abel-inverted
radial_distribution overlaid on the simulation's 3-D |v| histogram). Fixed
2026-07-14: those CSVs now also export pyabel `I(v)`. The paper_v2/v3/v4/cov
reference families were audited and are NOT affected (projection-space exports
vs projected-speed histograms — self-consistent; do not "fix" them).
