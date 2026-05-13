# Port decision: `post_process_single_pulse_paper.m`

Source:
`legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper.m`

## Decision

**Not ported.** The script is the pre-version-tagged ancestor of the
v2 / v3 / v4 paper post-processing scripts and contains no in-scope
post-processing recipe that the existing Python ports do not already
cover. Its remaining content is either out of scope per
`i2_helium_md/CLAUDE.md` or already delegated to CSV / MAT reference
exports under `data/reference/`.

This note exists so a future agent does not redo the same analysis.

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
- `i2_helium_md/postprocess/energy_balance.py` (`mass_spectrum`)
- `i2_helium_md/postprocess/velocity_2d.py`
- `i2_helium_md/postprocess/polar_velocity.py`
- `i2_helium_md/postprocess/_smoothing.py` (`moving_mean`)

## Recipe-by-recipe classification

BASE script structure:

- Lines 36–103: effusive (gas-phase) branch.
- Lines 105–205: non-effusive I+He droplet branch (external MATLAB
  references).
- Lines 212–357: shared simulation-side tail.

| BASE recipe | MATLAB lines | Python equivalent | Status |
|---|---|---|---|
| Mass selection (mass=131 droplet) | 219–235 | `energy_balance.mass_spectrum`, `paper_v3._paper_v3_selection` | DUPLICATE |
| Projected speed `sqrt(vx² + vy²)` at t=0 | 247 | `paper_v3_velocity_curve` (uses `velocities_final_x/y`) | DUPLICATE |
| Velocity bin edges `0:0.05:35 A/ps` | 251 | `PAPER_V3_VELOCITY_BIN_WIDTH_APS = 0.05`, `PAPER_V3_VELOCITY_MAX_APS = 35.0` | DUPLICATE |
| `movmean(h, 20)` velocity smoothing | 258 | `PAPER_V3_VELOCITY_SMOOTHING_WINDOW = 20` | DUPLICATE |
| `phi = atan2(vy, vx) + π`, `0:0.05:2π` edges, `movmean(h, 15)` | 283–295 | `paper_v3_phi_curve` (`PAPER_V3_PHI_BIN_WIDTH_RAD = 0.05`, `PAPER_V3_PHI_SMOOTHING_WINDOW = 15`) | DUPLICATE |
| 2-D vx/vy map, `[-35, 35]`, bin 0.5 A/ps | 301–331 | `paper_v2_velocity_map` | DUPLICATE |
| Side-by-side simulated vs. experimental 2-D VMI | 333–347 | `plot_paper_v2.py` main figure | DUPLICATE |
| `cos(phi)^n + offset` fit on radial-slice angular curve | 184–188 | (operates on the experimental quadrant slice, not on simulation) | OUT-OF-SCOPE — experimental processing delegated to reference exports |
| Quadrant-masked radial slice (`b_phi`, `mean(image_polar(b_phi, :), 1)`) | 73–102, 175–201 | not ported | OUT-OF-SCOPE — experimental, intentionally skipped |
| 3-D `surf(rr, phiphi, image_polar_select)` | 102, 197–201 | not ported | OUT-OF-SCOPE — visualization-only legacy utility per `CLAUDE.md` |
| External `plot_processed_VMI`, `subtract_processed_data`, `mean_timescan_2d_VMI([296:297], …)` | 116–134 | delegated to `data/reference/` CSV / MAT exports | OUT-OF-SCOPE — `CLAUDE.md` data contract |
| Effusive (gas-phase) branch | 36–103, 214–217 | not ported | OUT-OF-SCOPE — effusive dynamics per `CLAUDE.md` |

## Conventions in BASE that intentionally diverge from v3

The Python v3 port is faithful to `post_process_single_pulse_paper_v3.m`,
which itself differs from the BASE script in two subtle places. These
are recorded here so a future agent can recognise the divergence as
intentional, not as a porting regression.

1. **BASE subtracts the minimum before max-normalising the simulated
   projected-velocity curve.**
   - BASE lines 258–260:
     ```matlab
     vd_ion = movmean(h, 20);
     vd_ion = vd_ion - min(vd_ion);
     plot(centers_velocity * 100, vd_ion / max(vd_ion));
     ```
   - `paper_v3.py` `matlab_max_normalise` performs max-only normalisation
     (`y / max(y)`), matching the v3 MATLAB source. No min subtraction.
   - This is a baseline shift on the simulated trace, not a new recipe.
   - If a future visual comparison shows the Python v3 simulated curve
     sitting visibly above the experimental baseline near `v ≈ 0`, apply
     the min subtraction as a focused v3 refinement, not as a BASE port.

2. **BASE re-smooths the first three reference traces on `fig1` with a
   second `movmean(10)` after plotting.**
   - BASE lines 268–277 iterate over the first three line handles on
     `fig1` and replace `YData` with `movmean(YData, 10)`.
   - The Python v3 pipeline delivers reference traces already smoothed
     in the CSV export and does not double-smooth.
   - This is a presentation detail on the experimental reference, not a
     new simulation recipe.
   - If a future visual comparison shows the Python v3 experimental
     reference traces too noisy near small velocities, apply the
     `movmean(10)` second pass as a focused v3 refinement.

## What this note is not

- It is not a port. No Python module, plot script, test, reference CSV,
  run directory, or checkpoint is added or modified by this decision.
- It is not authorisation to port the effusive branch, the
  quadrant-masked angular slice, or the `cos(phi)^n + offset`
  experimental fit.

If a real Python / MATLAB figure mismatch later surfaces, treat that as
a focused v3 (or v2 / v4) refinement and reference the relevant
"Conventions in BASE that intentionally diverge from v3" item above.
