# Port decisions: superseded MATLAB ancestor scripts

This note records why two pre-version-tagged ancestor scripts under
`legacy_matlab_repository/single_pulse_simulation/` were **not ported**.
It exists so a future agent does not redo the analysis or re-port a
superseded ancestor.

Both ancestors were the precursors of the active paper_v2 / v3 / v4
scripts; their in-scope simulation-side panels are fully covered by the
existing Python ports, and their remaining content is either out of
scope per `CLAUDE.md` or already delegated to CSV / MAT reference
exports under `data/reference/`.

## Already-ported Python equivalents (shared by both ancestors)

| Python | MATLAB descendant |
|---|---|
| `scripts/post_processing/plot_paper_v2.py` | `post_process_single_pulse_paper_IplusHe_comparison.m` |
| `scripts/post_processing/plot_paper_v3.py` | `post_process_single_pulse_paper_v3.m` |
| `scripts/post_processing/plot_paper_v4.py` | `post_process_single_pulse_paper_v4.m` |

Helper modules: `i2_helium_md/postprocess/paper_v2.py`, `paper_v3.py`,
`paper_v4.py`, `velocity_2d.py`, `polar_velocity.py`,
`energy_balance.py` (`mass_spectrum`), `_smoothing.py`.

## 1. `post_process_single_pulse.m`

Source:
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse.m`

Script structure:

- Lines 4–56: effusive (gas-phase) branch.
- Lines 57–128: non-effusive I+He droplet branch (external MATLAB
  references + experimental quadrant fit + 3-D polar surface).
- Lines 130–204: shared simulation-side tail.

A short-lived port lived under `postprocess/single_pulse_base.py` and
`scripts/post_processing/plot_single_pulse_base.py`; it worked but added
no measurement-vs-simulation insight beyond what `plot_paper_v3.py` and
`plot_paper_v2.py` already provide, so it was removed in favour of this
decision note.

| Script recipe | MATLAB lines | Python equivalent | Status |
|---|---|---|---|
| Projected speed `sqrt(vx² + vy²)` | 139, 247 | `paper_v3_velocity_curve` (uses `velocities_final_x/y`) | DUPLICATE — same physics, narrower v-range |
| Velocity bin edges `0:0.05:22 A/ps` | 143 | `paper_v3` uses `0:0.05:35`; same recipe, just a different upper bound | DUPLICATE — already covered by the wider v3 grid |
| `movmean(h, 20)` velocity smoothing | 150 | unified `PAPER_V3_VELOCITY_SMOOTHING_WINDOW = 15` in the Python port (see `post_processing_strategy.md` §5) | NEAR-DUPLICATE — different window, equivalent physics |
| `vd_ion = vd_ion - min(vd_ion)` baseline shift | 151 | not in v3; baseline shift only, not a new recipe | NO COMPARISON VALUE |
| `phi = atan2(vy, vx) + π`, `0:0.1:2π` edges | 158, 160 | `paper_v3_phi_curve` uses bin width 0.05 (finer); same phi recipe otherwise | DUPLICATE — finer v3 bins subsume the 0.1 grid |
| 2-D `(vx, vy)` velocity map, `[-22:0.5:22]` | 178–195 | `velocity_density_2d` defaults match exactly | DUPLICATE |
| 3-D `scatter3(vx, vy, vz)` of all ions | 136 | not separately ported (visualisation-only) | NO COMPARISON VALUE |
| Experimental radial overlay (timescan + I+He) | 91–96 | `data/reference/paper_v3/*.csv` via `plot_paper_v3.py` | DUPLICATE |
| 3-D `surf(rr, phiphi, image_polar_select)` | 51–55, 120–126 | not ported | OUT-OF-SCOPE — experimental MATLAB construct, visualization-only |
| `cos(phi)^n + offset` fit on quadrant slice | 33–44, 103–115 | not ported | OUT-OF-SCOPE — operates on experimental polar image |
| Effusive (gas-phase) branch | 4–56 | not ported | OUT-OF-SCOPE per `CLAUDE.md` |
| External `plot_processed_VMI`, `mean_timescan_2d_VMI` | 7, 71, 116 | delegated to `data/reference/` exports | OUT-OF-SCOPE — data contract |

Candidates that looked novel:

- **Baseline subtract before max-normalise** (line 151). Shifts the
  simulated trace baseline before normalising. Apply as a focused v3
  refinement if a future visual comparison shows the curve sitting
  above the experimental baseline near `v ≈ 0`; do not re-port the
  whole script.
- **3-D `scatter3` of `(vx, vy, vz)`** (line 136). Pure visualisation of
  the same data already shown by the 2-D `(vx, vy)` map; three lines
  in a notebook on the existing checkpoint arrays produce it ad hoc.

## 2. `post_process_single_pulse_paper.m`

Source:
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper.m`

Script structure:

- Lines 36–103: effusive (gas-phase) branch.
- Lines 105–205: non-effusive I+He droplet branch (external MATLAB
  references).
- Lines 212–357: shared simulation-side tail.

| BASE recipe | MATLAB lines | Python equivalent | Status |
|---|---|---|---|
| Mass selection (mass=131 droplet) | 219–235 | `energy_balance.mass_spectrum`, `paper_v3._paper_v3_selection` | DUPLICATE |
| Projected speed `sqrt(vx² + vy²)` at t=0 | 247 | `paper_v3_velocity_curve` | DUPLICATE |
| Velocity bin edges `0:0.05:35 A/ps` | 251 | `PAPER_V3_VELOCITY_BIN_WIDTH_APS = 0.05`, `PAPER_V3_VELOCITY_MAX_APS = 35.0` | DUPLICATE |
| `movmean(h, 20)` velocity smoothing | 258 | unified `PAPER_V3_VELOCITY_SMOOTHING_WINDOW = 15` | NEAR-DUPLICATE — window unified to 15 across paper_v2/v3/v4 |
| `phi = atan2(vy, vx) + π`, `0:0.05:2π`, `movmean(h, 15)` | 283–295 | `paper_v3_phi_curve` (`PAPER_V3_PHI_BIN_WIDTH_RAD = 0.05`, `PAPER_V3_PHI_SMOOTHING_WINDOW = 15`) | DUPLICATE |
| 2-D vx/vy map, `[-35, 35]`, bin 0.5 A/ps | 301–331 | `paper_v2_velocity_map` | DUPLICATE |
| Side-by-side simulated vs. experimental 2-D VMI | 333–347 | `plot_paper_v2.py` main figure | DUPLICATE |
| `cos(phi)^n + offset` fit on radial-slice angular curve | 184–188 | not ported (experimental quadrant slice) | OUT-OF-SCOPE |
| Quadrant-masked radial slice (`b_phi`, `mean(image_polar(b_phi, :), 1)`) | 73–102, 175–201 | not ported | OUT-OF-SCOPE |
| 3-D `surf(rr, phiphi, image_polar_select)` | 102, 197–201 | not ported | OUT-OF-SCOPE — visualization-only |
| External `plot_processed_VMI`, `subtract_processed_data`, `mean_timescan_2d_VMI` | 116–134 | delegated to `data/reference/` exports | OUT-OF-SCOPE — data contract |
| Effusive (gas-phase) branch | 36–103, 214–217 | not ported | OUT-OF-SCOPE per `CLAUDE.md` |

### Intentional divergences from v3

The Python v3 port is faithful to `post_process_single_pulse_paper_v3.m`,
which itself differs from this BASE script in two subtle places.
Recorded here so a future agent recognises the divergence as
intentional, not as a porting regression.

1. **BASE subtracts the minimum before max-normalising the simulated
   projected-velocity curve** (lines 258–260):
   ```matlab
   vd_ion = movmean(h, 20);
   vd_ion = vd_ion - min(vd_ion);
   plot(centers_velocity * 100, vd_ion / max(vd_ion));
   ```
   `paper_v3.py` does max-only normalisation, matching the v3 MATLAB
   source. If a future visual comparison shows the Python v3 simulated
   curve sitting visibly above the experimental baseline near
   `v ≈ 0`, apply the min subtraction as a focused v3 refinement
   instead of porting BASE.
2. **BASE re-smooths the first three reference traces on `fig1` with a
   second `movmean(10)` after plotting** (lines 268–277). The Python
   v3 pipeline delivers reference traces already smoothed in the CSV
   export and does not double-smooth. Apply the second pass only if
   a future visual comparison shows the experimental references too
   noisy near small velocities.

## What this note is not

- Not a port. No Python module, plot script, test, reference CSV, run
  directory, or checkpoint is added or modified by these decisions.
- Not authorisation to port the effusive branch, the experimental
  polar surface, the quadrant-masked angular slice, or the
  `cos(phi)^n + offset` experimental fit.

If a real Python / MATLAB figure mismatch later surfaces on the
existing v2 / v3 / v4 ports, treat it as a focused refinement of the
affected paper script and reference the relevant row in section 1 or
the divergence note in section 2 rather than re-porting either
ancestor.
