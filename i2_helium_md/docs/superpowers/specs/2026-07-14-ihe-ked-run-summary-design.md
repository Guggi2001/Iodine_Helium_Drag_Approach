# Design: I⁺Heₙ KED reference integration into the run summary

**Date:** 2026-07-14
**Status:** approved (user, this date)
**Scope:** post-processing comparison layer only — no physics, no checkpoint
schema, no propagation changes.

## Goal

Replace the retired `vmi_summary` experimental overlay in
`scripts/post_processing/plot_run_summary.py` with comparisons against the new
frozen experimental reference `data/reference/ihe_ked/` (per-fragment mean
kinetic energies with a full error model for n = 0…17, plus trusted 2-D/3-D
curves for n = 0…4), and add a mass-spectrum vs. experimental-abundance
side-by-side.

## Decisions (user-confirmed)

1. Per-fragment curve overlays for n = 0…4, **both** representations:
   2-D detector projection and 3-D reconstruction.
2. Mass-spectrum section becomes a side-by-side against the Tier-2 abundance
   reference `data/reference/integrated_i_he_abundance.csv`.
3. The mean-energy table gets its **own section** (⟨E⟩ vs n, full error model)
   **and** supplies `v(⟨E⟩)` marker lines on the five curve overlays.
4. The old `vmi_summary` overlay section (`radial_velocity_with_vmi`) is
   retired; **no gas-phase I⁺ trace** — only the respective per-fragment
   curves are overlaid. The `vmi_summary` CSVs and `load_vmi_reference` stay
   (other scripts use them; deleting reference data is forbidden).
5. Overlays use the **per-unit-v measure on an m/s speed axis** (sim
   histograms are naturally per-unit-v; the `_PE` columns are the same data
   under a Jacobian). Mass-spectrum comparison uses **grouped bars per n**.

## Architecture

New package module **`i2_helium_md/postprocess/ihe_ked.py`** (chosen over
script-local loaders and over extending `velocity_distribution.py`, whose
contract is 2-column VMI CSVs). `plot_run_summary.py` only builds figures from
the module API. Loaders are testable and reusable by future Tier-2 report
scripts.

### Module API

- `load_ihe_ked_reference(path) -> IHeKedReference`
  Frozen dataclasses; one row per fragment with
  `n, massCenter_u_per_e, N_counts, N_eff, meanKE_eV, modeKE_eV, medianKE_eV,
  sigmaKE_eV, statErr_meanKE_eV, sysErr_meanKE_eV, calibSyst_frac,
  conditionSyst_frac, bgOffShift_eV, dominantError, noiseLimited`.
  Validates 18 rows and required columns; fails loudly.
- `load_ihe_ked_curve(dir, n) -> IHeKedCurve`
  Columns `E_eV, v_mps, signal_2d_Pv, signal_2d_PE, signal_3d_Pv,
  signal_3d_PE`; NaN cuts preserved as gaps; n restricted to 0…4.
- `fragment_mean_kinetic_energy(ion, n) -> FragmentMeanKE`
  Sim-side ⟨E⟩ in eV: mass-gate selection reused from the existing
  final-velocity helper; `E = ½·m(n)·|v_final|²` with
  `m(n) = 126.90 + 4.0026·n` u (the reference's mass model); returns mean,
  statistical error `σ_E/√N`, and gate count N. Both sim and reference are
  droplet-rest-frame → mean-to-mean comparison is direct (never
  mean-to-peak, per the reference README).

### Run-summary sections

- **`ihe_ked_mean_energy`** — ⟨E⟩ vs n. Experimental points with per-point
  error `√(statErr² + sysErr²)`; the calibration (4%) and condition (6%)
  fractional bands drawn as two separate shaded correlated envelopes (they
  shift the whole curve, not points); gold points (n = 0, 3–10, 12) visually
  distinguished. Sim ⟨E⟩ for every populated mass gate with statistical error
  and per-gate N annotation. Log-y (values span 3.7 → <0.1 eV).
- **`ihe_ked_curves_3d`** — one figure, 5 panels (n = 0…4): sim mass-gated 3-D
  |v| histogram (15-bin moving mean, peak-normalized — mirroring the
  reference's smoothed-envelope-peaks-at-1 convention) over `signal_3d_Pv` on
  the `v_mps` axis. Vertical markers at experimental and sim `v(⟨E⟩)`
  (labelled `v(⟨E⟩)`, not ⟨v⟩).
- **`ihe_ked_curves_2d`** — same layout: sim projected in-plane speed
  `√(vx² + vy²)` (the `paper_v2_velocity_curve` convention) over
  `signal_2d_Pv`, same markers.
- **`mass_spectrum`** — upgraded: grouped bars per n, sim fraction (%) next to
  `ionPercent`, via the existing `load_he_abundance_reference`.
- **Retired:** `radial_velocity_with_vmi` section and the three `VMI_REF_*`
  user settings. New settings: `IHE_KED_REFERENCE_DIR`
  (default `data/reference/ihe_ked`; `None` skips) and `ABUNDANCE_REF_PATH`
  (default `data/reference/integrated_i_he_abundance.csv`; `None` keeps the
  plain mass spectrum). Metadata section lists the new references.

## Comparison semantics (owned trade-offs)

- Reference curves are **independently peak-normalized** → overlays compare
  shapes only, never amplitudes across fragments or columns. Cross-fragment
  amplitude information enters exclusively through the mass-spectrum /
  abundance section.
- Low-statistics mass gates make sim ⟨E⟩ noisy; sections annotate N per gate
  rather than hide it. Gates with zero ions are skipped, not drawn at 0.
- Unit conversion: sim velocities are Å/ps; reference `v_mps` is m/s
  (1 Å/ps = 100 m/s). Reuse the existing conversion path in the
  velocity-histogram helpers.

## Testing

New `tests/test_ihe_ked.py`:

- reference-table loader against the real CSV (18 rows; n=0 mean = 3.7057 eV);
- curve loader: column contract, NaN handling in 3-D columns, n out of range;
- missing-file / missing-column loud failures;
- `fragment_mean_kinetic_energy` on a tiny synthetic checkpoint with an
  analytically known ⟨E⟩ (tight tolerance — analytical port);
- non-interactive figure-builder smoke test.

No figures or production-sized checkpoints generated in tests.

## Out of scope

Abel inversion, gas-phase comparison, pump-probe, re-fitting κ of ⟨E⟩(n),
2-D VMI image sections beyond the existing paper-v2 surface, any change to
`velocity_distribution.py`, `vmi_summary` data, or physics modules.
