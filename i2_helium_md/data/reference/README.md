# Reference data layout

This directory holds normalized, Python-consumable copies of the
experimental VMI and HeDFT trajectory references used by the
post-processing scripts.

## Subfolders

```text
paper_v2/       <- post_process_single_pulse_paper_IplusHe_comparison.m
                   (export_paper_v2_reference_data.m)
paper_v3/       <- post_process_single_pulse_paper_v3.m
                   (export_paper_v3_reference_data.m)
paper_v4/       <- post_process_single_pulse_paper_v4.m
                   (export_paper_v4_reference_data.m)
vmi_summary/    <- export_vmi_reference_data.m (legacy Abel-inverted summary
                   consumed by plot_experimental_comparison.py and the
                   consolidated plot_run_summary.py)
ihe_ked/        <- export_IHe_KED_reference.m in the VMI repo (fragment scan
                   171024): per-fragment I+He_n (n=0..17) kinetic-energy
                   reference with full error budget, for the H.2b forward
                   model. See ihe_ked/COLUMNS.md for the data dictionary.
scripts/        <- the MATLAB exporters that generate the above CSV/MAT files,
                   plus optional Python verification helpers.
```

Top-level files:

```text
9A_All_Data.csv                 <- HeDFT trajectory reference (9 angstrom droplet)
18A_All_Data.csv                <- HeDFT trajectory reference (normalized 18 angstrom droplet)
integrated_i_he_abundance.csv   <- experimental I+He_n size distribution (Tier-2 arbiter)
```

## `integrated_i_he_abundance.csv` — I+He_n size distribution (Tier-2 arbiter)

The measured I+He_n integer-n abundance the Tier-2 generative mass mechanism is
falsified against (loaded by `postprocess/abundance_loader.py ::
load_he_abundance_reference`; scored by the E4 integer-support Wasserstein
comparison). Data contract verified against the file 2026-07-03:

```text
column                    unit     meaning
n                         -        He-shell count; integer, contiguous 0..20
label                     -        species label: I^+, I^+He_1, ..., I^+He_20
massCenter_u_per_e        u/e      mass-window centre (integer I mass: 127 at n=0,
                                   step 4.0026 = He; NOT the sim-side 126.90 base)
massWindowLower_u_per_e   u/e      mass-window lower edge
massWindowUpper_u_per_e   u/e      mass-window upper edge
ionCounts                 -        per-species ion intensity (non-negative;
                                   fractional values, column sum ~ 0.19 -- a
                                   normalized intensity, NOT raw detector
                                   counts; the absolute normalization is part
                                   of the open provenance item below)
ionPercent                %        ionCounts as a percentage of the column
                                   total; sums to 100
```

21 rows (n = 0..20). The loader normalizes `ionPercent` by its own sum into
`ion_fraction` (sums to exactly 1.0). Note the reference support stops at n = 20
while the simulation's legal support runs to n = 21 (Langmuir cap n* = 21); the
E4 comparison zero-fills the union support, it does not clip.

Consumer (in-repo): `data/reference/scripts/plotting_histogram.py`.

**Provenance: to be completed.** This file ships without a documented producer
(no MATLAB/Python exporter under `scripts/` generates it, and no measurement-ID
record exists in-repo). The verified column/unit/normalization contract above is
authoritative for loading; the *measurement provenance* (source measurement IDs,
the producing script, calibration/scaling steps) is an open item — to be
supplied by the data owner.

## `ihe_ked/` — per-fragment kinetic-energy reference (H.2b arbiter)

Deployed 2026-07-14 from the VMI repo
(`matfile_data_scripts/A_state_paper_figures_single_pulse/fragment scan/171024/
reference_export/`, exporter `export_IHe_KED_reference.m`; regenerate there,
then re-copy). Contents:

```text
IHe_KED_reference.csv             <- one row per n = 0..17: meanKE_eV (the
                                     reference value), modeKE_eV/medianKE_eV
                                     (shape; distributions are strongly skewed),
                                     sigmaKE_eV, four error terms, flags
IHe_KED_spectra_n0_n1.csv         <- unsmoothed P(E) for n=0 and n=1
IHe_KED_curves_n0..n4.csv         <- trusted curves for n=0..4: raw 2-D detector
                                     projection AND 3-D speed distribution,
                                     peak-normalized, dual axes (E_eV + v_mps;
                                     v_mps is per-fragment mass-corrected)
IHe_KED_crosscheck_300mW.csv      <- independent same-day 300 mW series vs the
                                     reference (origin of conditionSyst_frac)
IHe_KED_reference.provenance.json <- machine-readable conventions/provenance
COLUMNS.md                        <- data dictionary for every CSV above,
                                     incl. the error-combination recipe
README.md                         <- full conventions, error model, and
                                     per-fragment trust guidance (gold points)
```

Overview plot: `scripts/experimental_reference_fragments.py` renders all five
trusted fragment curves in the four representation/measure combinations
(2-D/3-D × per-v/per-E) → `ihe_ked/experimental_reference_fragments.png`.

Key contract points: the reference value is the **mean** of the 3-D KED
(pyabel rIbeta speed distribution, proper dE weighting) in the droplet rest
frame; compare MD mean-to-mean, never mean-to-peak (see `shape_note` in the
provenance). Errors: per-point = sqrt(stat^2 + sys^2), plus TWO correlated
bands (`calibSyst_frac` = 0.04, `conditionSyst_frac` = 0.06) that shift all
fragments together. ⚠️ Units exception: this family is ENERGY-space (eV
columns, unit-suffixed), not m/s — the convention below applies to the
velocity-space families.

## On-disk units convention

All reference CSVs use **m/s** as the canonical velocity unit. CSV
columns are named with explicit unit suffixes (`v_mps`, `vx_mps`,
`vy_mps`).

Python loaders convert m/s → Å/ps (divide by 100) at load time so
downstream plotting and binning code can stay in the documented Å/ps
conventions (CLAUDE.md "Known Plotting Conventions"). Until the MATLAB
exporters have been re-run, some on-disk files may still carry the
legacy `v_Aps` column; the loaders detect this and convert without
warning.

## Known issues

- **43556 image center disagreement.** `paper_v2/iplus_he_160mw_43556_radial.csv`
  uses VMI center `[524.5297, 380.8430]`; `paper_v4/iplus_he_160mw_43556_radial.csv`
  uses `[509.3664, 387.6409]`. Same raw measurement, two different radial
  profiles. We do not currently have a reliable VMI-centering procedure
  to choose between the two. Treat overlays of `paper_v2/43556` and
  `paper_v4/43556` as approximate until centering is verified.
- **`paper_v3` I+He curves are the high-SNR co-add, not 43563.** Files
  there are named `iplus_he_high_snr_{radial,phi}.csv` to reflect this.
  If you need the standalone 43563 radial, use
  `paper_v4/iplus_he_300mw_43563_radial.csv` (which is also not the same
  measurement as the high-SNR co-add).
- **`vmi_summary/` is a different reduction pipeline.** See
  `vmi_summary/README.md`. It runs Abel inversion, image smoothing, and a
  mass correction not applied in the paper-era exports, and it averages
  raw measurements that no paper-era export uses. Its two I⁺He CSVs are
  **byte-identical** (one 45xxx-campaign reference, not two — see the
  warning in that README).

## Forbidden

Per CLAUDE.md, do not delete reference data or commit large checkpoints
into this directory. Keep reference CSVs small and inspectable.
