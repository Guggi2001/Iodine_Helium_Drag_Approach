# The `radial_positions.py` module — a walkthrough

## What problem does this file solve?

Once we know **how big** each helium droplet is (Step 7), we still need to
decide **where inside the droplet** each I₂ molecule sits before the laser
fires. The molecule is in thermal equilibrium with the droplet at
temperature `T_particles_K` (≈0.4 K for a real droplet). So we draw its
radial position from the Boltzmann distribution in the droplet solvation
potential.

This module does that draw.

## Position in the dependency chain

```
config.py
   ↓
physics/constants.py
physics/potentials.py        (provides droplet_potential)
   ↓
sampling/radial_positions.py  ← THIS MODULE
   ↓
simulation/neutral.py        (uses radial samples to set initial positions)
```

It depends on `droplet_potential` from `physics/`. It does **not** depend
on droplet sizes — instead it accepts an arbitrary array of droplet radii
and samples one position per droplet.

## The physics

We sample `r` (distance from droplet center) from

$$p(r) \propto r^2 \cdot \exp\left(-\frac{U_{\text{drop}}(r - R)}{k_B T}\right)$$

Two factors:

### 1. The `r²` Jacobian
In 3D, the spherical volume element is `4π r² dr`. Even for a uniform
density (no potential), there's just **more space** at larger `r`. This
is why the peak of the radial distribution is **not at r=0** even when
the molecule "wants" to sit at the droplet center.

### 2. The Boltzmann factor `exp(-U/kT)`
The droplet potential `U(r-R)` is approximately zero deep inside the
droplet and rises sharply to the binding energy as `r` crosses the
droplet surface. At low temperatures, the Boltzmann factor strongly
penalizes positions near or outside the surface, so molecules
concentrate inside.

## What the temperature does

Three regimes:

| T | Sampler behaviour |
|---|---|
| **T → 0**  | Tight peak near `r ≈ 0` — molecule sits at the droplet center, only the `r²` Jacobian gives any spread |
| **T moderate (≈ 4 K)** | Peak shifts outward, distribution broadens — molecule has thermal energy to climb the soft well |
| **T → ∞** | `exp(-U/kT) ≈ 1`, so `p(r) ∝ r²` — uniform-density spherical sampling |

This behaviour is illustrated in `radial_positions_visualization.png`
(generated alongside the test suite). The histogram of sampled radii
matches the analytical density to within Monte Carlo noise.

## Public API

```python
from i2_helium_md.sampling.radial_positions import sample_radial_positions

r = sample_radial_positions(cfg, droplet_radii)
# r.shape == (cfg.num_molecules,)
```

Typical use:

```python
from i2_helium_md import single_pulse_N2000
from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes
from i2_helium_md.sampling.radial_positions import sample_radial_positions

cfg = single_pulse_N2000(num_molecules=500, seed=42, T_particles_K=0.4)

# Step 7: pick droplet sizes
N = sample_droplet_sizes(cfg, mode="post_pickup")
droplet_radii = (3 * N / (4 * np.pi * 0.0218 * 0.8)) ** (1/3)  # A

# Step 8: place each molecule inside its droplet
r = sample_radial_positions(cfg, droplet_radii)
```

## Internal walkthrough

### `sample_radial_positions(cfg, droplet_radii, ...)`

The user-facing function. Groups molecules by their **unique** droplet
radius and runs the rejection sampler once per group, then unscatters the
results back to per-molecule order.

The `np.unique` + boolean indexing approach is much faster than running
the sampler N times for N molecules when many droplets share the same
radius (e.g. the `single_droplet_size = 2000` case where every droplet
is identical — we run the sampler exactly once and split the result).

Reproducibility comes from `cfg.seed` if set; otherwise a fresh
`np.random.default_rng()` is used. The same `rng` is shared across all
unique radii, so the sequence is deterministic regardless of how many
distinct sizes appear.

### `_radial_probability_density(r, droplet_radius, T, steepness, binding)`

Evaluates the unnormalized Boltzmann density
`r² · exp(-U(r-R) / (k_B T))` on an array of `r` values.

Used both for finding the sampling envelope (the maximum value of `p` on
a fine grid, used as `y_max` in rejection sampling) and for evaluating
acceptance probabilities at proposal points.

The Boltzmann constant is converted to eV/K (`K_B / EV`) so that all
energy quantities live in eV — matching the convention used by
`droplet_potential`.

### `_sample_radial_for_one_droplet(...)`

Rejection sampler for a single droplet radius:

1. **Build envelope on a grid** (`r_step = 0.01 Å` by default).
2. **Estimate acceptance rate** as `mean(p) / max(p)` on the grid. This
   is more efficient than the MATLAB version's hardcoded 1000-proposals-
   per-batch loop, because we know roughly how many proposals to make
   upfront.
3. **Vectorized rejection sample** in batches sized by the estimate
   plus a 5% safety factor.
4. **Top up** if we under-shoot (rare with the safety factor).

A `RuntimeError` is raised if the density is identically zero (which
would mean the configuration is unphysical — for instance, T = 0 with
an unboundable potential).

## Departures from MATLAB

These are intentional and documented:

1. **No plotting code.** The legacy function had inline matplotlib-like
   plots for each droplet radius. Visualization belongs in the
   postprocessing layer.

2. **Energy unit consistency.** MATLAB went meV → J → eV → meV in three
   different places via the awkward `[0:0.001:E_max]/1000*eV` form. We
   work in eV throughout, with explicit conversion only at the
   `K_B / EV` step.

3. **Single normalization step.** MATLAB normalized `p(E)` once over an
   energy axis, then renormalized `p_radius(r)` over the radial axis.
   The first normalization was unused — we drop it.

4. **Vectorized batch sizing.** MATLAB used a fixed 1000 proposals per
   batch in a `while` loop. We estimate the acceptance rate from the
   grid-sampled density and propose the right number of samples in one
   batch.

5. **`np.unique` deduplication.** MATLAB did this too (`unique_radii`)
   but our version cleanly preserves per-molecule order through
   `np.unique(..., return_inverse=True)`. The MATLAB version
   appended results in unique-radius order, requiring an external
   permutation to map back.

## Testing

The test file `tests/test_radial_positions.py` covers:

- **The density function**: zero at the origin (because of `r²`),
  positive inside the droplet, suppressed far outside, and `r²`-shaped
  in the high-T limit.
- **The single-droplet sampler**: returns the correct count, samples
  stay in `[0, 2R)`, seeded reproducibility, and at `T = 0.4 K` more
  than 90% of samples land inside the droplet.
- **The public sampler**: returns one sample per molecule, handles
  heterogeneous droplet sizes correctly (each slot uses the right
  `R`), reproducible with a fixed seed, and accepts 2D droplet-radius
  inputs (e.g. column vectors from MATLAB-style code).

A visual verification image (`radial_positions_visualization.png`)
shows the sampled histograms tracking the analytical density at three
temperatures (0.4 K, 4 K, 50 K).

## Common pitfalls

1. **Don't expect the peak at `r = 0`** even at T = 0. The `r²` Jacobian
   shifts the peak outward by a few Å.

2. **Energy units**: `cfg.binding_energy_molecule_meV` is in **meV**
   (matching the MATLAB convention), but `droplet_potential` expects
   **eV**. The function divides by 1000 internally.

3. **Heterogeneous droplet radii**: passing 1000 random radii means we
   run the rejection sampler 1000 times. For best performance, share
   radii (e.g. via post-pickup sampling that yields integer-rounded
   sizes) or use `use_single_droplet_size = True`.

## Future improvements

1. **Inverse CDF sampling.** For one-shot sampling tasks, building the
   CDF on a grid and interpolating its inverse is faster than rejection
   sampling and has uniform accuracy. This becomes worthwhile if we
   start sampling > 10⁵ molecules per run.

2. **Direct sampling without the grid.** The Boltzmann density is
   smooth enough that we could compute `r² exp(-U/kT)` analytically
   without the envelope grid. This would remove the `r_step`
   parameter from the public API.

3. **Optional `phi`, `theta` sampling.** Currently we only sample radial
   positions — angular components are sampled separately by the main
   simulation loop. Combining them here would simplify the call site.

## References

- The "place a thermal point particle inside a smooth potential well"
  problem is standard textbook physics; we follow the formulation used
  in the original MATLAB code, which in turn cites the helium droplet
  literature for the potential shape. See `docs/physics_background.md`
  for the connection to the larger simulation pipeline.
