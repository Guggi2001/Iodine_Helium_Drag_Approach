# B2 Coulomb-explosion timestep-refinement benchmark (deterministic, single pair)

2026-10-20. Preserves audit finding **F1** of the §4.4 numerical audit
(B2 production path). Script `scripts/bench_b2_coulomb_dt.py`; results
`data/runs/b2_coulomb_dt_benchmark/b2_coulomb_dt_benchmark.csv` (40 rows).

**What this is not.** It is a deterministic benchmark of one pair at fixed mass.
It is separate from the six-member stochastic B2 timestep-refinement study, and
nothing here is a B2 ensemble result. The −1.62 meV-per-pair 5-term ledger
residual quoted next to F1 in the audit was read from the proxy run
`…g413t64e310s1_th500` (lq c 2.0, `raw` sampler, E₀ 0.310 eV). It is **not a B2
result** and is not reproduced here.

## Configuration

| item | value |
|---|---|
| force | production `make_ion_accel_fn`: FD pair Coulomb (h = 10⁻⁴ Å, `interactions.py:183`) + erf droplet well (14.2 Å, E_bind 0.11727 eV) |
| integrator | production `make_ion_baoab_step` (B-A-O-A-B), T_eff = 0 (noise off) |
| drag | `off` (γ ≡ 0) or `B2`: `linear_quadratic`, a 18.007, c 1.6, `erf_independent` gate 3.5 Å (bundle `srho3p5A18A_linear_quadratic`) |
| cfg source | `gen_tier2_detfix_postpickup.build_lq_member("s1", B2_E0_EV, B2_TAU_PS, B2_C)`. Diff against the committed B2 s1 run `cfg.json` = ∅ (checked 2026-10-20) |
| pair | 2 I⁺, born at rest, R₀ = 2.666 Å, centre at depth −35.00 Å in R_drop = 53.47 Å |
| geometry | `radial`: axis along the droplet radius (one ion outward, one inward through the droplet). `tangential`: axis perpendicular (ions symmetric) |
| mass | fixed, `complex_mass_amu(21)` = 210.9546 amu per ion. **No** pickup / evaporation / cooling / mass jumps |
| time | 40 ps. By then both ions have left helium in every case (final depths 129–576 Å) |
| dt series | 0.02, 0.01, 0.005, 0.0025, 0.00125 ps (0.02 kept from the original series) |
| observable | final per-ion KE [eV] (production `_amu_ang2_ps2_to_eV`) |

## Method

- `dKE_vs_finest` = KE(dt) − KE(0.00125).
- `ratio_vs_finest` = dKE(dt) / dKE(dt/2). It is biased near the finest step,
  because the reference's own error enters both terms.
- **Reference-free** `increment` = KE(dt/2) − KE(dt), and `increment_ratio` =
  increment(dt) / increment(dt/2).
- **dt = 0.00125 ps is a numerical reference, not an exact solution.**

## Results (from the committed CSV; tangential ions are identical, shown once)

| geometry | drag | ion | dt [ps] | KE [eV] | KE − KE(0.00125) [eV] | rel. | ratio vs finest | increment KE(dt/2)−KE(dt) [eV] | increment ratio |
|---|---|---|---|---|---|---|---|---|---|
| radial | off | outward | 0.02 | 2.573532662 | -3.1385e-03 | -1.22e-03 | 4.065 | 2.3664e-03 | 4.021 |
| radial | off | outward | 0.01 | 2.575899058 | -7.7210e-04 | -3.00e-04 | 4.204 | 5.8846e-04 | 4.005 |
| radial | off | outward | 0.005 | 2.576487519 | -1.8364e-04 | -7.13e-05 | 5.001 | 1.4692e-04 | 4.001 |
| radial | off | outward | 0.0025 | 2.576634442 | -3.6719e-05 | -1.43e-05 | — | 3.6719e-05 | — |
| radial | off | outward | 0.00125 | 2.576671160 | 0.0000e+00 | 0.00e+00 | — | — | — |
| radial | off | inward | 0.02 | 2.575031292 | -3.1368e-03 | -1.22e-03 | 4.065 | 2.3651e-03 | 4.021 |
| radial | off | inward | 0.01 | 2.577396438 | -7.7170e-04 | -2.99e-04 | 4.204 | 5.8815e-04 | 4.005 |
| radial | off | inward | 0.005 | 2.577984589 | -1.8354e-04 | -7.12e-05 | 5.001 | 1.4685e-04 | 4.001 |
| radial | off | inward | 0.0025 | 2.578131434 | -3.6699e-05 | -1.42e-05 | — | 3.6699e-05 | — |
| radial | off | inward | 0.00125 | 2.578168133 | 0.0000e+00 | 0.00e+00 | — | — | — |
| radial | B2 | outward | 0.02 | 1.065784067 | -2.2823e-03 | -2.14e-03 | 3.064 | 1.5374e-03 | 3.114 |
| radial | B2 | outward | 0.01 | 1.067321501 | -7.4486e-04 | -6.97e-04 | 2.967 | 4.9379e-04 | 2.763 |
| radial | B2 | outward | 0.005 | 1.067815286 | -2.5107e-04 | -2.35e-04 | 3.470 | 1.7872e-04 | 2.470 |
| radial | B2 | outward | 0.0025 | 1.067994010 | -7.2350e-05 | -6.77e-05 | — | 7.2350e-05 | — |
| radial | B2 | outward | 0.00125 | 1.068066359 | 0.0000e+00 | 0.00e+00 | — | — | — |
| radial | B2 | inward | 0.02 | 0.183582078 | -1.1211e-03 | -6.07e-03 | 2.776 | 7.1718e-04 | 2.799 |
| radial | B2 | inward | 0.01 | 0.184299260 | -4.0391e-04 | -2.19e-03 | 2.735 | 2.5622e-04 | 2.495 |
| radial | B2 | inward | 0.005 | 0.184555479 | -1.4769e-04 | -8.00e-04 | 3.282 | 1.0269e-04 | 2.282 |
| radial | B2 | inward | 0.0025 | 0.184658172 | -4.4999e-05 | -2.44e-04 | — | 4.4999e-05 | — |
| radial | B2 | inward | 0.00125 | 0.184703171 | 0.0000e+00 | 0.00e+00 | — | — | — |
| tangential | off | either (symmetric) | 0.02 | 2.574278194 | -3.1377e-03 | -1.22e-03 | 4.065 | 2.3658e-03 | 4.021 |
| tangential | off | either (symmetric) | 0.01 | 2.576643978 | -7.7190e-04 | -2.99e-04 | 4.204 | 5.8831e-04 | 4.005 |
| tangential | off | either (symmetric) | 0.005 | 2.577232287 | -1.8359e-04 | -7.12e-05 | 5.001 | 1.4688e-04 | 4.001 |
| tangential | off | either (symmetric) | 0.0025 | 2.577379171 | -3.6709e-05 | -1.42e-05 | — | 3.6709e-05 | — |
| tangential | off | either (symmetric) | 0.00125 | 2.577415881 | 0.0000e+00 | 0.00e+00 | — | — | — |
| tangential | B2 | either (symmetric) | 0.02 | 0.603430271 | -1.8115e-03 | -2.99e-03 | 2.899 | 1.1867e-03 | 2.938 |
| tangential | B2 | either (symmetric) | 0.01 | 0.604616991 | -6.2477e-04 | -1.03e-03 | 2.830 | 4.0398e-04 | 2.606 |
| tangential | B2 | either (symmetric) | 0.005 | 0.605020972 | -2.2079e-04 | -3.65e-04 | 3.357 | 1.5501e-04 | 2.357 |
| tangential | B2 | either (symmetric) | 0.0025 | 0.605175982 | -6.5778e-05 | -1.09e-04 | — | 6.5778e-05 | — |
| tangential | B2 | either (symmetric) | 0.00125 | 0.605241760 | 0.0000e+00 | 0.00e+00 | — | — | — |

## Interpretation

- **At production dt = 0.01 ps, relative to dt = 0.00125 ps:**
  - Drag off: KE is low by **−0.772 meV per ion** (−3.0·10⁻⁴), with no dependence on geometry.
  - B2 drag: **−0.745 / −0.404 / −0.625 meV** (outward / inward / tangential), i.e. **−0.07 % / −0.22 % / −0.10 %** of the post-exit KE.
  - Every deviation is negative (one-signed).
- **Drag off.** The increment ratios are 4.02 → 4.01 → 4.00, so the conservative (Coulomb + well) error falls ≈ 4× per halving. Within this series dt = 0.00125 is converged to the 10⁻⁵ eV level.
- **B2 drag on.** The increment ratios fall with dt: 3.11 → 2.47 (outward), 2.80 → 2.28 (inward), 2.94 → 2.36 (tangential).
  - The drag O-step adds a component that converges more slowly than the conservative part, consistent with the frozen-γ O-step. Its order is not determined here.
  - Because the ratio is still falling, the dt = 0.00125 reference itself carries an unmeasured error. The listed dt = 0.01 deficits therefore probably understate the error relative to dt → 0.
  - A rough extrapolation, not a measurement: the missing part is of order the last increment, 4.5–7.2·10⁻⁵ eV, i.e. ~10–20 % of the dt = 0.01 value.
- **Dissipated energy** E_dissip (pair) changes by ≤ 2.8·10⁻³ eV over the series (CSV column).
- **Calibration context.** The drag law was extracted with this same integrator at dt = 0.01 (`extraction/trajectory_matching.py:168`), so part of this bias is absorbed into the calibrated (a, c).

## Scope limits

- **Fixed masses, no stochastic jumps.** The coupling between event timing and dt is not tested here.
- **One birth geometry** (median droplet and depth) and three axis orientations.
- **No convergence order** is established for the complete B2 model (variable mass, Bernoulli events, density-gated cooling and pickup). This note gives none.

## Provenance and recovery

- **Origin.** The script is the recovered 2026-10-20 session script (`dtconv.py`, scratchpad), cleaned up.
- **Recovery check.** `--session-original` restores the session's literals:
  - mass 210.95858 amu (atomic I − mₑ + 21 He at 4.002602),
  - birth depth −34.99 Å,
  - eV constant 9648.533212 amu·Å²/ps².

  In that mode it reproduces every value the session printed, to all printed digits.
- **Committed numbers** use the production mass (210.9546 amu), the production eV conversion and depth −35.00 Å. They differ from the session values by ≤ 1.1·10⁻⁷ eV in the dt = 0.01 deficits, so the audit's F1 figures are unchanged at their quoted precision.

## Reproduce / verify

```powershell
& 'C:\Users\user\AppData\Local\Programs\Python\Python314\python.exe' scripts/bench_b2_coulomb_dt.py --check
```

`--check` reruns the benchmark (~1 min) and compares it string-exactly against the
committed CSV. Verified on 2026-10-20: `CHECK PASSED: 40 rows`.
