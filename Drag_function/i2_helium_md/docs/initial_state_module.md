# The `simulation/initial_state.py` module

## What problem does this file solve?

The neutral propagation needs a complete physical state at `t=0` and pre-allocated
trajectory arrays for the rest of the run. This module produces that initial
state by:

1. Sampling droplet sizes and radii
2. Sampling molecule positions inside the droplets
3. Sampling molecular orientations + bond lengths
4. Computing initial velocities from the laser parameters
5. Combining everything into a fully-allocated `NeutralCheckpoint` with
   column 0 populated and the rest zero

The rest of the simulation (driver in `simulation/neutral.py`, not yet
written as of this commit) writes columns 1..num_steps-1 in place by
time-stepping.

This module replaces the variable-initialization block of
`vmi_sim_3d_neutral_propa_HeDFT_mimic.m` (lines ~284-483, excluding the
optional DFT pre-fill which is left unimplemented for now).

## Position in the dependency chain

```
   physics.constants ─────── droplet_radius_bulk_angstrom
                              EV, HC, MASS_I_AMU, U
                                  │
   sampling.droplet_sizes ────────┤
   sampling.radial_positions ─────┤
   sampling.orientations ─────────┤
   physics.potentials ────────────┘
                                  │
                                  ▼
                       initial_state.build_initial_state
                                  │
                                  ▼
                        NeutralCheckpoint with col 0 filled
```

## Public API

```python
from i2_helium_md.simulation.initial_state import build_initial_state

ckpt = build_initial_state(cfg, num_steps=200, rng=np.random.default_rng(0))
# ckpt.positions_x  shape (2N, 200), col 0 populated, rest zero
# ckpt.velocities_x shape (2N, 200), col 0 populated, rest zero
# ckpt.E_kin_eV     shape (2N, 200), col 0 populated, rest zero
# ckpt.r0           shape (N,)      molecule centre radial distance
# ckpt.E_initial_eV shape (N,)      photon energy delivered per molecule
# ...
```

## What's inside

### Step 1: Droplet sizes
- If `cfg.use_single_droplet_size`: assign every molecule the same N
  (typical HeDFT comparison run).
- Otherwise: call `sample_droplet_sizes(cfg, mode="post_pickup")`.

### Step 2: Droplet radii
Use `physics.constants.droplet_radius_bulk_angstrom(N)`, which uses
the bulk helium density (NOT the 0.8x density used in pickup
sampling). See that function's docstring for the legacy
inconsistency.

### Step 3: Molecule centre positions (β, γ, r0)
Pure radial+angular: r0 from `sample_radial_positions`, β/γ from
`sample_orientations`'s position-angle outputs.

### Step 4: Axis orientation (α, δ) + bond length
From `sample_orientations`. In single-pulse mode the axis is
cos²-weighted (laser polarisation); otherwise isotropic. Bond length
= R0_GS + N(0, deltaR0).

### Step 5: Initial speed v0

Two complications worth flagging:

**(a) The variable name `E_initial` in MATLAB is misleading.** It's
``hc/lambda * eV`` and ends up in **joules**, not eV, despite the name.
We use the same formula but compute `E_initial_J` with explicit units;
the eV-valued `E_initial_eV` per-molecule field on the checkpoint
stores `hc/lambda` (in eV) for reference.

**(b) The mean velocity formula has a `partner_interaction` branch:**

* If `partner_interaction=True`: `mean_v = sqrt(E_initial / m)`.
  The full photon energy goes into kinetic energy.
* Otherwise: `mean_v = sqrt((E_initial - E_diss) / m)`.
  Subtract the dissociation energy first.

This is because with the partner interaction (Morse) included, the
dissociation energy is naturally subtracted as the atoms separate
along the Morse repulsive wall; without it, the kinetic-energy
budget already has E_diss removed.

**(c) FWHM velocity spread:**

```
fwhm_v = (1/2) * fwhm_E_eV * EV / sqrt(E_initial_J * m_kg) / 100
```

The `/100` converts m/s -> Å/ps. The `EV` factor mixes with `eV` and
sometimes-J variables; we follow the MATLAB literally and document
this oddness here rather than refactor it.

**(d) Single-pulse override:** in single-pulse mode (`cfg.single_pulse=True`)
all initial velocities are set to zero, regardless of laser parameters.
This matches MATLAB and makes physical sense: the single-pulse
simulation freezes time at the moment of photoexcitation, so atoms
haven't yet moved.

### Step 6: Atomic xyz from molecule centre + axis + bond

For each molecule i:

```
x_centre = r0[i] · cos(β[i]) · sin(γ[i])
y_centre = r0[i] · sin(β[i]) · sin(γ[i])
z_centre = r0[i] · cos(γ[i])

x_atom1[i] = x_centre + cos(α[i])·sin(δ[i]) · bond/2
y_atom1[i] = y_centre + sin(α[i])·sin(δ[i]) · bond/2
z_atom1[i] = z_centre + cos(δ[i])           · bond/2

x_atom2[i] = x_centre - cos(α[i])·sin(δ[i]) · bond/2
y_atom2[i] = y_centre - sin(α[i])·sin(δ[i]) · bond/2
z_atom2[i] = z_centre - cos(δ[i])           · bond/2
```

The MATLAB writes `sin(δ+π) = -sin(δ)` and `cos(δ+π) = -cos(δ)`
explicitly to make atom 2's offset opposite from atom 1; we use the
algebraic simplification.

Velocities are similarly constructed: along the molecular axis with
shared speed v0, atom 2 has the opposite sign.

### Step 7-9: Allocate trajectory arrays

All `(2N, num_steps)` arrays are allocated as zeros, with column 0
filled from steps 5-6 above. Static arrays (`mass_kg`, `droplet_radii`)
are sized `(2N,)` by tiling the per-molecule values.

`E_kin[:, 0]` and `E_pot[:, 0]` are computed from positions and
velocities at column 0. `E_dissip[:, 0]` and `L_droplet[:, 0]` start
at zero (cumulative quantities).

**`E_pot[:, 0]` includes the partner Morse term**, matching the
convention used by `propagation_step.py` for all subsequent timesteps:

```
E_pot_per_atom = droplet_potential_atom(r - R_droplet)
                  + Morse_pair / 2          # half-half split
```

The legacy MATLAB code OMITS the Morse term at t=0 (line 476 of
`vmi_sim_3d_neutral_propa_HeDFT_mimic.m`) but INCLUDES it from t=1
onward (line 885), causing a multi-eV discontinuity between t=0 and
t=1 that broke energy conservation tests. We treat that as a legacy
bug and include the Morse term at t=0 here.

`E_initial_eV` is per-molecule (size N) and stores the photon energy
in eV for diagnostic and bookkeeping purposes.

## Per-atom vs per-molecule arrays

| Array | Shape | Notes |
|---|---|---|
| `positions_*` | `(2N, T)` | Full trajectory |
| `velocities_*` | `(2N, T)` | Full trajectory |
| `mass_kg` | `(2N,)` | Per-atom mass (= MASS_I_AMU * U for both atoms) |
| `droplet_radii` | `(2N,)` | Tiled per-molecule radius |
| `E_kin_eV` | `(2N, T)` | Per-atom KE |
| `E_pot_eV` | `(2N, T)` | Per-atom PE (droplet + half partner) |
| `E_dissip_eV` | `(2N, T)` | Per-atom cumulative dissipation |
| `L_droplet_eV_ps` | `(2N, T)` | Per-atom cumulative pathlength inside droplet |
| `r0` | `(N,)` | Per-molecule initial radial distance |
| `E_initial_eV` | `(N,)` | Per-molecule photon energy (eV) |
| `time_ps` | `(T,)` | Time axis |

The 2N layout convention (atom 1 at indices 0..N-1, atom 2 at N..2N-1)
is preserved throughout the codebase. To recover per-molecule values:

```python
N = ckpt.num_molecules
E_kin_per_molecule = ckpt.E_kin_eV[:N] + ckpt.E_kin_eV[N:]
```

## Departures from MATLAB

1. **Single function returning a dataclass** instead of dropping
   `x_components`, `y_components`, ..., `mass`, `time`, `E_kin`,
   `E_pot`, `E_dissip` into the global workspace.

2. **Validates `num_steps >= 1`** with a clear ValueError. MATLAB would
   silently produce empty arrays.

3. **No `b_invalid_spawn` flag.** MATLAB tracks atoms that spawn
   outside their droplet but doesn't actually use that flag for
   anything in the propagation loop. We omit it; the test suite
   verifies atoms spawn near the droplet centre as expected from
   the radial sampler.

4. **No DFT pre-fill (yet).** The `custom_DFT_start` block of MATLAB
   is intentionally left out. When that path is needed (HeDFT
   comparison runs with TD-HeDFT initial conditions), we'll add a
   separate function `apply_dft_prefill(state, cfg)` that's invoked
   optionally by the driver. Until then, calling with
   `cfg.custom_DFT_start=True` will need the driver to raise
   NotImplementedError (not this function's responsibility).

5. **No debug plots.** MATLAB scatters `(x0, y0, z0)` in a figure;
   our smoke tests verify the geometry instead.
