# Python vs MATLAB neutral-propagation comparison

This directory contains a deterministic side-by-side comparison of the
Python neutral-propagation code against the legacy MATLAB
implementation.

## What's tested

A simple, hand-built initial condition with no RNG involvement:

- **1 I2 molecule** with the bond aligned along the x-axis
- **Bond length: 3.0 Å** (past the Morse equilibrium at 2.666 Å, far
  from the Xdip perturbation at 9 Å) — atoms see a strong attractive
  Morse force and oscillate around the equilibrium.
- **Both atoms at rest** at t=0
- **Single He droplet of fixed radius 27.97 Å** (= MATLAB's
  `2.22 * 2000^(1/3)`)
- **Collisions disabled** by setting the geometric scattering cross
  section to 0 (deterministic; no RNG involved)
- **Xdip turned off** in both implementations (cleaner Morse curve)
- **100 steps** of dt = 0.01 ps each → 1 ps total propagation

The simulation should produce a clean Morse oscillation. Both
implementations use velocity-Verlet integration, so they should agree
to numerical precision modulo small differences in physical constants.

## Files

- `run_python.py` — Python script that produces `python_trajectory.csv`
- `run_matlab.m`  — MATLAB script that produces `matlab_trajectory.csv`
- `compare.py`    — Loads both CSVs and reports deviations

## How to run

### 1. Run the Python side

From the project root:

```bash
python scripts/compare_with_matlab/run_python.py
```

This writes `python_trajectory.csv` next to the script.

### 2. Run the MATLAB side

You need the legacy MATLAB code on your machine. The `run_matlab.m`
script expects to find it at a path stored in the `MATLAB_LEGACY_DIR`
variable at the top of the file. Edit that to point to your copy of
`Iodine_Helium_Simulation/`.

Then from the same directory:

```matlab
run_matlab
```

(or `matlab -batch "run_matlab"` from the command line)

This writes `matlab_trajectory.csv` next to the script.

### 3. Compare

```bash
python scripts/compare_with_matlab/compare.py
```

The script prints:
- Side-by-side values at selected timesteps
- Max absolute / relative deviations across all 101 rows
- A pass/fail verdict on whether the agreement is within ~0.5%

## What "agreement" means here

Python and MATLAB use slightly different physical constants:

| Constant | Python (CODATA 2022 / SI 2019)       | MATLAB legacy        | Relative diff |
|----------|---------------------------------------|----------------------|---------------|
| EV [J]   | 1.602176634e-19 (exact)               | 1.602e-19 (rounded)  | ~110 ppm      |
| U  [kg]  | 1.66053906892e-27                     | 1.66053907e-27       | ~0.05 ppm     |

The leapfrog converts forces from eV/Å to Å/ps² using `EV`, so the
~110 ppm difference compounds across the 100 steps. We expect:
- **Positions**: scaled-relative deviation ~110 ppm (well within the constants budget)
- **Energies**: scaled-relative deviation ~110-300 ppm (similar)
- **Total energy drift over the run**: ~1–2% from leapfrog symplectic
  error in both implementations (similar in both)

If the deviation is **much larger** (more than ~0.5%), that indicates a
real bug in either implementation, not just constant differences.

## A note on the t=0 E_pot fix

The legacy MATLAB code in `vmi_sim_3d_neutral_propa_HeDFT_mimic.m`
contains a small bug: line 476 sets `E_pot(:,1)` (the t=0 column) to
the **droplet potential only**, while line 885 sets all subsequent
columns to **droplet + half partner Morse**. This causes a multi-eV
discontinuity in `E_pot` between t=0 and t=1 in the legacy output.

In our Python port, we identified this bug via an energy-conservation
test (`E_kin + E_pot + E_dissip` should be conserved across timesteps)
and **fixed it** in `build_initial_state` per project principle #10
("don't preserve legacy approximations").

**The MATLAB script in this folder uses the FIXED formula at t=0** —
it explicitly calls `add_partner_interaction` to get `E_pot_partner_0`
and includes it in the t=0 energy. This makes the comparison test
internally consistent (both implementations use the same correct
formula at t=0) and lets the energy-conservation drift be small.

If you want to reproduce the legacy bug for verification, change the
t=0 block in `run_matlab.m` to omit the partner Morse — you'll see a
multi-eV jump in `E_pot` between t=0 and t=1 that matches what the
legacy code produces.

### Why "scaled-relative" not pure-relative?

The standard relative-error metric `|diff| / |value|` blows up when
the underlying value passes through zero. In our test scenario the
atoms pass through the Morse equilibrium (~0.64 ps) where the
potential energy momentarily approaches zero. A 110 ppm constant
difference + tiny phase offset means Python reports `E_pot ≈ 2e-4`
and MATLAB reports `E_pot ≈ 1e-4` at that instant — physically
identical, but the relative error reads as 100% which is misleading.

The scaled-relative metric divides by the **typical magnitude** of
the column over the whole trajectory (`max |column|`), which gives
a meaningful "fraction of typical signal" number that stays well-
behaved across zero crossings. The compare script also flags any
column that passes near zero so you know when this matters.

## What's NOT tested

This comparison **deliberately** avoids:

- The droplet-size sampler (RNG-dependent)
- Initial-state assembly (RNG-dependent)
- Hard-sphere collisions (RNG-dependent)
- The Xdip Gaussian dip (we set it to off in both)
- Mass attachment (`attach_he`, ion-only)
- Effusive dynamics

These can't be cross-referenced deterministically because Python's
`numpy.random.default_rng` and MATLAB's `rng` produce fundamentally
different sequences. Statistical tests of those would be a separate
exercise.
