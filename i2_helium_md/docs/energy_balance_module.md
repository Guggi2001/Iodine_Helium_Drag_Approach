# The `postprocess/energy_balance.py` module

## What problem does this file solve?

The legacy MATLAB single-pulse pipeline produces three live-debug
figures (drawn at the end of `vmi_sim_3d_neutral_propa_HeDFT_mimic.m`
and `vmi_sim_3d_ion_propa.m`) plus the `post_process_single_pulse_paper_v3.m`
paper figure. We reproduce them post-hoc from the saved checkpoints
rather than inside the propagation drivers, so the simulation never
opens a window and the recipes stay testable in isolation.

This module owns the **pure functions** that the four legacy-debug
plot scripts call. No plotting, no I/O.

## Position in the dependency chain

```
NeutralCheckpoint    IonCheckpoint
        │                │
        ▼                ▼
postprocess/energy_balance.py   ← THIS MODULE
        │
        ▼
scripts/post_processing/
    plot_neutral_energy_balance.py
    plot_ion_energy_balance.py
    plot_ion_temperature_diagnostic.py
    plot_paper_figure.py
```

## Public API

```python
from i2_helium_md.postprocess import (
    EnergyTotals, neutral_energy_totals, ion_energy_totals,
    PhiHistogram, phi_histogram,
    MassSpectrum, mass_spectrum,
)
```

### `neutral_energy_totals(ckpt)` → `EnergyTotals`

Sum-over-atoms traces. Mirrors `vmi_sim_3d_neutral_propa_HeDFT_mimic.m`
line 965: each component is `sum(E_*, 1)` and `E_system = E_kin +
E_pot + E_dissip`. The `E_mass_attach_defect_eV` field is `None` for
the neutral stage.

### `ion_energy_totals(ckpt)` → `EnergyTotals`

Per-molecule traces (`sum / num_molecules`). Mirrors
`vmi_sim_3d_ion_propa.m` line 898. `E_system` includes
`E_mass_attach_defect_eV` so the running total is conserved up to
Verlet symplectic drift.

### `phi_histogram(ckpt, *, bin_width_rad=0.05, mass_amu=None, mass_tolerance_amu=0.5)` → `PhiHistogram`

Histograms `np.mod(arctan2(vy_final, vx_final) + pi, 2*pi)`. Mirrors
`post_process_single_pulse_paper_v3.m` line 314. Optional mass
selection via `mass_amu ± mass_tolerance_amu` (in amu) on
`mass_final_kg`.

### `mass_spectrum(ckpt, *, bin_width_amu=1.0)` → `MassSpectrum`

Histograms `mass_final_kg / U` with bin edges at half-integer amu so
that 127.0 lands in the bin centred on 127, not on a bin edge.
Mirrors `post_process_single_pulse_paper_v3.m` line 397
(`histogram(data_ion.mass_i(:,end)/u)`).

## Companion module: `_smoothing.py`

`postprocess/_smoothing.py` exposes two shared helpers:

- `moving_mean(values, window)` -- MATLAB-style `movmean` with
  shortened endpoint windows. Used to apply the `15`-bin smoothing
  the legacy paper figure applies to simulation curves.
- `normalise_trace(values)` -- baseline-subtract and scale to unit
  maximum.

Both `scripts/post_processing/plot_experimental_comparison.py` and
`scripts/post_processing/plot_paper_figure.py` import from here so
the smoothing convention has a single source of truth.

## Why post-hoc rather than live mid-run

CLAUDE.md restricts changes to neutral and ion propagation physics.
Because every input array (`E_kin_eV`, `E_pot_eV`, `E_dissip_eV`,
`L_droplet_eV_ps` for neutral; the same plus
`E_mass_attach_defect_eV` for ion) is already saved in the v2 / v4
checkpoints, the energy-balance figures are reproducible without
touching the propagation modules.

The temperature-diagnostic figure is the one exception: the legacy
MATLAB `diagnostic_array` row is built inside the collision block
from quantities (lab-frame `COStheta_lab`, `rho`, pre/post-collision
energies) that were not previously persisted. To reproduce the figure
faithfully we added a `(num_steps, 3)` `temperature_diagnostic` field
to `IonCheckpoint` (schema v4 -> v5) populated during ion propagation.
The third column is the lab-frame angle (mean `arccos(COStheta_lab)`),
matching MATLAB `vmi_sim_3d_ion_propa.m:561` where `theta` is built
from the lab-frame post-smearing cosine. See `docs/checkpoint_module.md`
and `docs/ion_propagation_step_module.md` for the capture path.

## Testing

`tests/test_energy_balance.py` covers:

- `E_system == E_kin + E_pot + E_dissip[ + E_mass_attach_defect]`
  on synthetic checkpoints,
- per-molecule division for the ion totals,
- uniform-angle phi histogram (bin count, density integrates to 1),
- mass selection filter passes the right atoms,
- 1-amu mass spectrum bin alignment (127 / 131 / 135 amu peaks
  land in the integer-centred bin, not on an edge),
- temperature-diagnostic recipe matches the scalar formula on
  hand-picked inputs and returns NaN when no atom collided,
- `apply_collision(return_diagnostics=True)` returns the 5-tuple,
- v5 round-trip preserves NaN sentinels,
- malformed `(T, 2)` shapes are rejected at load time,
- a v4 file (synthetic, missing `temperature_diagnostic`) raises the
  standard "schema_version=4, this code expects 5" error.
