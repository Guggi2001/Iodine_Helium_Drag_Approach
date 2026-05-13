# Port decision: `post_process_single_pulse.m`

Source:
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse.m`

## Decision

**Not ported.** The script is the pre-version-tagged ancestor of the
paper_v2 / v3 / v4 post-processing scripts. Its in-scope
simulation-side panels are fully covered by the existing Python ports
and contribute no new comparison value; its remaining panels are out
of scope per `i2_helium_md/CLAUDE.md`.

This note exists so a future agent does not redo the same analysis or
re-port a superseded ancestor.

A short-lived port lived under
`i2_helium_md/postprocess/single_pulse_base.py` and
`i2_helium_md/scripts/post_processing/plot_single_pulse_base.py`. It
worked but added no measurement-vs-simulation insight beyond what
`plot_paper_v3.py` and `plot_paper_v2.py` already provide, so it was
removed in favour of this decision note.

## Already-ported Python equivalents

- `i2_helium_md/scripts/post_processing/plot_paper_v2.py`
  ← `post_process_single_pulse_paper_IplusHe_comparison.m`
- `i2_helium_md/scripts/post_processing/plot_paper_v3.py`
  ← `post_process_single_pulse_paper_v3.m`
- `i2_helium_md/scripts/post_processing/plot_paper_v4.py`
  ← `post_process_single_pulse_paper_v4.m`

Helper modules:

- `i2_helium_md/postprocess/paper_v2.py`,
  `i2_helium_md/postprocess/paper_v3.py`,
  `i2_helium_md/postprocess/paper_v4.py`
- `i2_helium_md/postprocess/velocity_2d.py` (`velocity_density_2d`,
  whose docstring already cites this exact MATLAB script and whose
  defaults match `velocity_bins = -22:0.5:22`).
- `i2_helium_md/postprocess/_smoothing.py` (`moving_mean`).

## Recipe-by-recipe classification

Script structure:

- Lines 4–56: effusive (gas-phase) branch.
- Lines 57–128: non-effusive I+He droplet branch (external MATLAB
  references + experimental quadrant fit + 3-D polar surface).
- Lines 130–204: shared simulation-side tail.

| Script recipe | MATLAB lines | Python equivalent | Status |
|---|---|---|---|
| Projected speed `sqrt(vx² + vy²)` | 139, 247 | `paper_v3_velocity_curve` (uses `velocities_final_x/y`) | DUPLICATE — same physics, narrower v-range |
| Velocity bin edges `0:0.05:22 A/ps` | 143 | `paper_v3` uses `0:0.05:35`; same recipe, just a different upper bound | DUPLICATE — already covered by the wider v3 grid |
| `movmean(h, 20)` velocity smoothing | 150 | `PAPER_V3_VELOCITY_SMOOTHING_WINDOW = 20` | DUPLICATE |
| `vd_ion = vd_ion - min(vd_ion)` baseline subtract before max-divide | 151 | not in v3; baseline shift only, not a new recipe | NO COMPARISON VALUE |
| `phi = atan2(vy, vx) + π`, `0:0.1:2π` edges | 158, 160 | `paper_v3_phi_curve` uses bin width 0.05 (finer); same phi recipe otherwise | DUPLICATE — finer v3 bins subsume the 0.1 grid |
| 2-D `(vx, vy)` velocity map, `[-22:0.5:22]` | 178–195 | `velocity_density_2d` defaults match exactly | DUPLICATE |
| 3-D `scatter3(vx, vy, vz)` of all ions | 136 | not separately ported; see note below | NO COMPARISON VALUE |
| Experimental radial overlay (timescan probe-only + I+He) | 91–96 | reuse `load_paper_v3_radial_reference` with `data/reference/paper_v3/` CSVs in `plot_paper_v3.py` | DUPLICATE — same CSVs, same overlay |
| 3-D `surf(rr, phiphi, image_polar_select)` (experimental polar image) | 51–55, 120–126 | not ported | OUT-OF-SCOPE — `res.image_polar` is an experimental MATLAB construct delegated to `data/reference/` exports per `CLAUDE.md`; visualization-only legacy utility |
| `cos(phi)^n + offset` fit on experimental quadrant slice | 33–44, 103–115 | not ported | OUT-OF-SCOPE — operates on experimental polar image, delegated to reference exports |
| Effusive (gas-phase) branch | 4–56 | not ported | OUT-OF-SCOPE — `CLAUDE.md` forbids effusive dynamics without explicit request |
| External `plot_processed_VMI`, `mean_timescan_2d_VMI([296:297], …)` | 7, 71, 116 | delegated to `data/reference/` CSV / MAT exports | OUT-OF-SCOPE — `CLAUDE.md` data contract |

## Notes on the candidates that looked novel

- **Baseline subtract before max-normalise** (line 151,
  `vd_ion = vd_ion - min(vd_ion)`). Shifts the baseline of the
  simulated trace before normalising. Not a new recipe — at most a
  cosmetic choice. If a future visual comparison shows the v3
  simulated curve sitting visibly above the experimental baseline
  near `v ≈ 0`, apply this min-subtract step as a focused v3
  refinement instead of porting this whole script.
- **3-D scatter of `(vx, vy, vz)` at the final state** (line 136).
  Pure visualisation of the same data already shown by the 2-D
  `(vx, vy)` map (`velocity_density_2d`). Adds no quantitative
  comparison surface, no new physics, and no MATLAB / experimental
  reference. If wanted ad hoc, three lines in a notebook on the
  existing `IonCheckpoint.velocities_final_x/y/z` arrays produce it.

## What this note is not

- It is not a port. No Python module, plot script, test, reference
  CSV, run directory, or checkpoint is added or modified by this
  decision. Any earlier `single_pulse_base.py` /
  `plot_single_pulse_base.py` files are removed.
- It is not authorisation to port the effusive branch, the
  experimental polar surface, or the experimental quadrant
  `cos(phi)^n + offset` fit.

If a real Python / MATLAB figure mismatch later surfaces on the
existing v2 / v3 / v4 ports, treat that as a focused refinement of
the affected paper script and reference the relevant row above
rather than re-porting `post_process_single_pulse.m`.
