# The `velocity_distribution.py` module — a walkthrough

## What problem does this file solve?

The bottom tile of the legacy
`legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
simulation_image.m` figure (lines 159-256) overlays four 1-D
velocity-space curves at the detector:

1. experimental I⁺ from gas phase,
2. experimental I⁺He from helium droplets,
3. simulated I⁺He fragments (mass = 131 amu = 127 I + 4 He),
4. simulated I⁺He₂ fragments (mass = 135 amu = 127 I + 8 He).

The two simulated curves come from histogramming `|v_final|` for the
subset of ion-stage atoms that match a target final mass and that
exited the droplet (`b_ion_outside`).

This module factors that out of the plotting layer:

- `load_vmi_reference(path)` reads either
  `data/reference/vmi_iplus_he.csv` or `vmi_iplus_gas.csv`
  (already Abel-inverted 1-D spectra exported from the legacy MATLAB
  pipeline);
- `compute_final_velocity_histogram(ion, mass_amu=...)` ports the
  MATLAB selection + binning step, replacing the legacy `bayes_hist`
  call with `np.histogram` (we drop the Bayesian uncertainty bars
  because the plot doesn't show them).

## Public API

```python
from i2_helium_md.postprocess import (
    VmiReference, load_vmi_reference,
    FinalVelocityHistogram, compute_final_velocity_histogram,
)

# Experimental references
he  = load_vmi_reference("data/reference/vmi_iplus_he.csv")
gas = load_vmi_reference("data/reference/vmi_iplus_gas.csv")
# he.velocity_Aps  -> shape (M,)  Å/ps
# he.signal_arb    -> shape (M,)  raw signal

# Simulated histograms
ion = RunDirectory("data/runs/single_pulse_N_2000").load_ion()
sim_he  = compute_final_velocity_histogram(ion, mass_amu=131.0)
sim_he2 = compute_final_velocity_histogram(ion, mass_amu=135.0)
# sim_he.bin_centers_Aps, sim_he.density   -> arrays for plotting
# sim_he.num_atoms_used                    -> scalar diagnostic
```

## The mass + outside filter

```python
mass_amu_per_atom = round(ion.mass_final_kg / U)            # 2N atoms
mass_mask         = abs(mass_amu_per_atom - mass_amu) <= mass_tolerance_amu
outside_per_atom  = concat([b_ion_outside, b_ion_outside])  # broadcast N -> 2N
select            = mass_mask & outside_per_atom            # if require_outside
```

`b_ion_outside` is per-molecule; both atoms of an "outside" molecule
contribute. The default `mass_tolerance_amu=0.5` is a physical 1-amu
window centred on the requested mass — large enough to absorb any
floating-point noise in `mass_final_kg / U`, small enough to keep
131 amu and 135 amu populations disjoint.

`require_outside=False` is available for diagnostic runs where the
simulation hasn't been propagated long enough for ions to leave the
droplet, so every molecule's flag is still `False`.

`ValueError` if zero atoms pass — the caller must respond rather than
silently producing an all-zero distribution.

## Histogram contract

| Field | Meaning | Units |
|---|---|---|
| `bin_centers_Aps` | bin centres | angstrom/ps |
| `bin_edges_Aps`   | bin edges, length B+1 | angstrom/ps |
| `counts`          | raw atom counts per bin | dimensionless |
| `density`         | `counts / bin_width` | atoms·ps/angstrom |
| `mass_amu`        | filter target | amu |
| `num_atoms_used`  | filter cardinality | dimensionless |

Density is what the plotting layer divides by `density.max()` to
overlay against the unit-normalised experimental curves.

## Differences vs the legacy MATLAB

- No Bayesian uncertainty bars (`bayes_hist`). The plotting layer only
  renders the central line, so we use `np.histogram`.
- No `movmean(..., 15)` smoothing pass on the simulated curves. The
  raw histogram is what `compute_final_velocity_histogram` returns;
  smoothing, if desired, belongs to a downstream call.
- No 3D-vs-2D projection switch. We bin the full 3D speed
  `sqrt(vx² + vy² + vz²)`; the legacy `abel_inv_post=true` branch did
  the same, so the default behaviour matches the published figures.

## Out of scope

Experimental references already arrive Abel-inverted (see
`post_processing_strategy.md` §4 Strategy B). Plotting lives in
`scripts/post_processing/plot_hedft_comparison.py` and
`plot_experimental_comparison.py`. Project-wide scope rules live in
`CLAUDE.md`.
