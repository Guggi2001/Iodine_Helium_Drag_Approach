# Post-processing strategy: connecting the MD simulation to experimental VMI

This document explains how the molecular-dynamics simulation output is turned
into a signal that can be compared against experimentally measured VMI
(velocity-map imaging) images, what physical assumptions are required, and
why each in-scope post-processing script applies its own histogram recipe.

It consolidates two findings:

1. **The conceptual gameplan** (Strategy A vs Strategy B) — how the same set
   of MD trajectories is mapped onto two complementary observables, what
   role Abel inversion plays, and which assumptions are baked in.
2. **The script-by-script recipe differences** between
   `plot_paper_v2.py`, `plot_paper_v3.py`, `plot_paper_v4.py`, and
   `plot_experimental_comparison.py`. The radial velocity curves these
   scripts produce can look very different *for the same run* — the
   differences are intentional MATLAB-faithful conventions, not porting
   bugs.

The companion documents in this directory cover the implementation of
individual modules (`velocity_distribution_module.md`,
`compare_trajectories_module.md`, `energy_balance_module.md`, etc.) and the
overall porting plan (`post_processing_port_plan.md`).

## 1 — Physics being compared

In the lab, a UV pump pulse photodissociates I₂ inside a helium nanodroplet.
A femtosecond probe pulse ionises the iodine fragments, which are then
accelerated by a VMI spectrometer onto a position-sensitive detector. The
key facts:

- The ion's position on the detector encodes its velocity components in
  the **detector plane** (here `vx`, `vy`).
- The velocity component along the spectrometer axis (`vz`) is integrated
  out by the projection — every `vz` slice of the 3-D velocity sphere
  contributes to the same detector pixel.
- A time-of-flight (TOF) electronic gate downstream selects only ions of
  a specific mass-to-charge ratio.

The MD simulation, by contrast, knows the **full 3-D velocity vector**
`(vx, vy, vz)` per atom at the end of propagation, the per-atom mass, and
a per-molecule boolean `b_ion_outside` flagging whether the ion has
actually escaped the droplet. There is no instrument projection or
detector blur in the simulation.

The comparison problem reduces to a single question:

> *Do we project the simulation down to match the raw VMI 2-D image, or
> do we Abel-invert the experimental image up to match the simulation's
> 3-D distribution?*

The codebase implements **both** strategies in parallel.

## 2 — Coordinate conventions

| Quantity | Convention |
| --- | --- |
| Laser polarisation axis | **z** (cos² angular distribution lives along z) |
| VMI detector plane | **(x, y)** — `vz` is the integrated-out axis |
| Simulation final velocities | `vx_final`, `vy_final`, `vz_final` per atom |
| Simulation mass | `mass_final_kg` per atom |
| Simulation escape flag | `b_ion_outside` per molecule |
| MD propagation length | `ion_simulation_time = 20 ps` (`i2_helium_md/config.py`) |

Final MD velocities are taken to be **asymptotic detector-arrival
velocities** — after 20 ps of ion propagation the ions have escaped the
droplet and no longer interact with the helium density, so their
velocities at the MD endpoint equal the velocities they would have at
the detector under free flight.

## 3 — Strategy A: project the simulation, compare against raw VMI

Used by `plot_paper_v2.py`, `plot_paper_v3.py`, `plot_paper_v4.py`.

```text
simulation final 3-D velocities
    │
    ├── filter: b_ion_outside == True              (only ions that escaped droplet)
    ├── filter: round(mass_amu) in {131} or {127,131} (TOF mass gating)
    │
    ▼
project to detector plane:  speed_2D = sqrt(vx² + vy²)    (vz integrated out)
    │
    ▼
1-D histogram with paper-specific bins / movmean / baseline rule
    │
    ▼
plotted against `data/reference/paper_v{2,3,4}/*_radial.csv`,
which are 1-D radial profiles extracted directly from the raw 2-D
VMI detector image (NOT Abel-inverted).
```

The comparison is therefore **2-D detector signal vs 2-D detector
signal**. Strategy A:

- never assumes cylindrical symmetry,
- never propagates inversion noise,
- keeps both sides of the comparison in instrument-native coordinates,
- includes whatever instrument/geometry projection bias is present in
  the experimental projection, since the simulation is projected the
  same way (so those biases cancel at first order).

This is the path used by the **publishable paper figures**: it is the
most defensible quantity because no inversion assumption is required.

## 4 — Strategy B: Abel-invert the experimental image, compare against full 3-D simulation

Used by `plot_experimental_comparison.py` and
`postprocess/velocity_distribution.py`.

```text
simulation final 3-D velocities
    │
    ├── filter: b_ion_outside == True
    ├── filter: TOF mass gating
    │
    ▼
full 3-D speed:  speed_3D = sqrt(vx² + vy² + vz²)
    │
    ▼
1-D histogram (`0:0.04:26` A/ps, movmean 15, baseline-subtracted)
    │
    ▼
plotted against `data/reference/vmi_summary/vmi_iplus_he.csv` (and gas),
which were produced by `data/reference/scripts/export_vmi_reference_data.m`.
That MATLAB script:
    • calls `abel_invert_processed_VMI(...)` on the measured image,
    • applies a velocity calibration factor `vf_single ≈ 8.6178`
      (m/s per pixel-radius unit; encodes 7.95 kV repeller voltage),
    • applies a sqrt(127/131) mass-velocity correction for the IHe⁺
      branch (KE-to-speed scaling at fixed kinetic energy).
```

The comparison is therefore **3-D radial distribution vs 3-D radial
distribution**.

### Why Abel inversion is required here

A VMI image is the Abel projection of the 3-D velocity distribution
onto the detector plane:

```text
I_2D(x, y) = ∫ I_3D(x, y, z) dz
```

under the assumption that the 3-D distribution is **cylindrically
symmetric** about the polarisation axis (`z`). Given that symmetry, the
inverse Abel transform recovers `I_3D(r, θ)` — equivalently, the
central slice through the 3-D velocity sphere. The radial profile of
that central slice is what is directly comparable to the simulation's
3-D `sqrt(vx² + vy² + vz²)` distribution.

Without inversion, the raw 2-D radial profile contains contributions
from every `vz` slice and is **not** the same physical quantity as a
3-D radial cut. That is why Strategy A pairs a 2-D simulation speed
with a raw 2-D reference, and Strategy B pairs a 3-D simulation speed
with an Abel-inverted reference. **Mixing strategies is a unit
error.**

## 5 — Recipe differences across the four scripts

Even within Strategy A, the three paper scripts use different
smoothing windows and one of them (v4) baseline-subtracts. Strategy B
uses its own bins and window. All four recipes are MATLAB-faithful.

| Step | `plot_paper_v2.py` | `plot_paper_v3.py` | `plot_paper_v4.py` | `plot_experimental_comparison.py` |
|---|---|---|---|---|
| Strategy | A (2-D projected) | A (2-D projected) | A (2-D projected) | B (3-D vs Abel-inverted) |
| Module called | `postprocess/paper_v2.py` | `postprocess/paper_v3.py` | `postprocess/paper_v4.py` | `postprocess/velocity_distribution.py` |
| Source checkpoint | `ion.npz` of selected run | `ion.npz` of selected run | `ion.npz` of selected run | `ion.npz` of `single_pulse_droplet` |
| Outside-droplet filter | `b_ion_outside == 1` | `b_ion_outside == 1` | `b_ion_outside == 1` | `require_outside=True` (default) |
| Mass filter | `round(mass) == 131` (IHe⁺) | `round(mass) == 131` (IHe⁺) | `round(mass) ∈ {127, 131}` (I⁺ + IHe⁺) | per `velocity_distribution.py` defaults |
| Speed formula | `sqrt(vx² + vy²)` (2-D) | `sqrt(vx² + vy²)` (2-D) | `sqrt(vx² + vy²)` (2-D) | **`sqrt(vx² + vy² + vz²)` (3-D)** |
| Bin edges (A/ps) | `np.arange(0, 35.05, 0.05)` → 700 bins | `np.arange(0, 35.05, 0.05)` → 700 bins | `np.arange(0, 35.05, 0.05)` → 700 bins | **`0:0.04:26`** → 650 bins |
| `movmean` window | **10** bins | **20** bins | **40** bins | **15** bins |
| Baseline subtraction | none | none | **`smoothed − min(smoothed)`** | **`smoothed − min(smoothed)`** (via `normalise_trace`) |
| Final normalisation | `y / max(y)` | `y / max(y)` | `y / max(y)` after baseline shift | `y / max(y)` after baseline shift |
| x-axis range plotted | 0–3500 m/s | 0–3500 m/s | 0–3500 m/s | **0–2800 m/s** |
| Experimental overlay CSVs | `data/reference/paper_v2/*_radial.csv` (4 traces) | `data/reference/paper_v3/*_radial.csv` (2 traces incl. timescan) | `data/reference/paper_v4/*_radial.csv` (1 trace) | `data/reference/vmi_summary/vmi_iplus_he.csv`, `vmi_iplus_gas.csv` |
| MATLAB provenance | `post_process_single_pulse_paper_IplusHe_comparison.m` (`movmean` 10) | `post_process_single_pulse_paper_v3.m` (`movmean` 20) | `post_process_single_pulse_paper_v4.m` (`movmean` 40, baseline shift) | `simulation_image.m` velocity overlay, see `CLAUDE.md` §"Known Plotting Conventions" |

### Ranked impact on the simulation curve shape

For a *given* finished run, what makes the simulation radial curve
look different across these four scripts, in order of visual impact:

1. **Strategy choice (2-D vs 3-D speed).** Strategy A peaks at lower
   speeds than Strategy B because the 2-D projection systematically
   removes `vz` contributions. This is the single largest difference
   between a paper panel and the experimental-comparison panel.
2. **`movmean` window (10 → 15 → 20 → 40 bins).** With 0.05 A/ps bins
   the windows are 0.5 / 0.75 / 1.0 / 2.0 A/ps respectively; v4 is
   ~4× more smoothed than v2 and loses most sub-peak structure.
3. **Baseline subtraction.** v4 and `plot_experimental_comparison`
   subtract `min(smoothed)` before max-normalising, which lifts the
   floor and stretches the peak relative to v2 and v3.
4. **Mass filter.** v4 admits both `127` (I⁺) and `131` (IHe⁺); v2
   and v3 admit only `131`. The two ion populations have different
   speed distributions — this is real physics, not smoothing.
5. **Bin width and range.** 0.04 vs 0.05 A/ps, range 26 vs 35 A/ps.
   Small effect on resolution and high-velocity tail.
6. **Reference CSVs and plotted x-range.** Affect only the
   *experimental* overlay, not the simulation curve.

Within a single strategy, comparing across scripts is therefore a
comparison of smoothing/baseline conventions, not of physics. Across
strategies, the projection-vs-inversion choice dominates.

## 6 — Assumptions baked into the comparison

1. **Cylindrical symmetry about z (Strategy B only).** Required for
   the inverse Abel transform of the experimental reference. The MD
   ensemble has this symmetry naturally because the initial molecular
   orientation is sampled isotropically and the only symmetry-breaking
   direction is the laser polarisation along z.

2. **Asymptotic free-flight velocities.** After 20 ps of ion
   propagation, MD-final velocities are taken to equal detector
   arrival velocities. The `b_ion_outside` flag excludes molecules
   still trapped in the droplet at the end of MD; their final
   velocities are not yet asymptotic.

3. **TOF mass gating ↔ simulation mass filter.** The experimental
   detector sees only ions whose flight time matches the gate window,
   i.e. only ions of a specific m/z. The simulation matches this by
   selecting `round(mass) == 131` (IHe⁺) or `127` (I⁺). Mass is
   per-atom in `IonCheckpoint`; both atoms of a molecule must pass
   the filter.

4. **Detector response is already corrected in the experimental
   data.** The simulation histogram is unweighted; the experimental
   processing pipeline (in `imported_vmi_reference.py` /
   `export_vmi_reference_data.m`) is responsible for any
   gain/geometry corrections.

5. **No further fragmentation after MD ends.** The simulation assumes
   ion identity (I⁺ vs IHe⁺) is frozen at the MD endpoint. Any
   post-droplet evaporation of the He partner — converting a mass-131
   ion into a mass-127 ion in flight to the detector — is not
   modelled.

6. **Velocity calibration of the experimental reference is correct.**
   The reference relies on `vf_single = 8.6178` m/s per pixel-radius
   unit (set by repeller voltage 7.95 kV and VMI focal geometry).
   Calibration drift shifts the experimental peak horizontally
   relative to the simulation and is not corrected for in the Python
   pipeline.

7. **Image centering of the experimental reference is correct.**
   Off-centre VMI images broaden the radial profile artificially.
   The MATLAB export uses hard-coded image centres
   `[524.5, 380.8]` (IHe⁺) and `[482.9, 392.5]` (I⁺ gas). Commit
   `306aaeb` flagged this as an open question — there is no
   automated verification that these centres remain optimal for the
   present image set.

8. **Anisotropy (β) is recovered from φ, not from the full 3-D
   vector.** `polar_velocity.py` fits `f(φ) = a + b·cos(φ - φ₀)²`
   and extracts `β = 2b / (2a + b)` (range `[-1, +2]`). Because the
   simulation has full 3-D velocities, no Abel inversion is needed
   on the simulation side for this analysis.

## 7 — Practical implications for review

When eyeballing run-summary figures across scripts:

- **Same script, different runs** — peak shape differences are run
  physics. Compare directly.
- **Same run, paper script vs paper script** — peak position is the
  same, smoothness and baseline differ. This is convention, not
  physics.
- **Paper script vs `plot_experimental_comparison`** — peak position
  *will* differ. The paper script peaks at the 2-D projected speed,
  the experimental-comparison script peaks at the 3-D speed. Do not
  read this as disagreement; it is the projection-vs-inversion
  choice.
- **Comparison against MATLAB** — for each script, compare against
  its own MATLAB ancestor (e.g. `plot_paper_v2.py` against
  `post_process_single_pulse_paper_IplusHe_comparison.m`). Cross-
  script MATLAB comparison would re-introduce the strategy mismatch.

## 8 — Critical files

Simulation side:

- `i2_helium_md/simulation/ion.py` — final velocities and `b_ion_outside` flag.
- `i2_helium_md/simulation/checkpoint.py` — `IonCheckpoint` schema (mass, velocities, outside flag, temperature diagnostic).
- `i2_helium_md/config.py` — `ion_simulation_time = 20 ps` propagation length.

Strategy A (2-D vs raw 2-D VMI):

- `i2_helium_md/postprocess/paper_v2.py`, `paper_v3.py`, `paper_v4.py`
- `scripts/post_processing/plot_paper_v2.py`, `plot_paper_v3.py`, `plot_paper_v4.py`
- `data/reference/paper_v{2,3,4}/*_radial.csv`
- `data/reference/scripts/export_paper_v{2,3,4}_reference_data.m`

Strategy B (3-D vs Abel-inverted):

- `i2_helium_md/postprocess/velocity_distribution.py`
- `scripts/post_processing/plot_experimental_comparison.py`
- `data/reference/vmi_summary/vmi_iplus_he.csv`, `vmi_iplus_gas.csv`
- `data/reference/scripts/export_vmi_reference_data.m` (calls `abel_invert_processed_VMI`)
- `data/reference/scripts/imported_vmi_reference.py` (velocity factor, centering)

Anisotropy:

- `i2_helium_md/postprocess/polar_velocity.py`

Shared smoothing / normalisation helpers:

- `i2_helium_md/postprocess/_smoothing.py` (`movmean`, `normalise_trace`)

## 9 — Open questions

- The export-script provenance of
  `data/reference/paper_v{2,3,4}/*_radial.csv` should be re-verified
  to confirm those CSVs are raw (un-inverted) 2-D radial profiles.
  This conclusion follows from the absence of an
  `abel_invert_processed_VMI` call in the v2/v3/v4 export scripts,
  but a direct read of each `.m` file is the authoritative check.
- The image-centering coordinates in the VMI summary export are
  hard-coded and were flagged in commit `306aaeb`. No automated
  check verifies they remain optimal for the present image set.
- Velocity-factor calibration (`vf_single = 8.6178`) is a single
  fixed constant; any drift in repeller voltage or detector
  geometry between data sets manifests as a global horizontal
  shift of the experimental curve relative to the simulation.
- Numerical RMSE of the Python simulation curves against
  `data/reference/paper_v{2,3,4}/iplus_he_*_radial.csv` is not yet
  measured. If a port bug is suspected this would be the
  authoritative check.
