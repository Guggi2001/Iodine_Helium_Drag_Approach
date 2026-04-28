# The `orientations.py` module

## What problem does this file solve?

For each molecule in the ensemble, we need five random scalars at the
start of a simulation:

| Symbol | Role | Distribution |
|---|---|---|
| `beta` (β) | azimuth of molecule centre position | uniform on [0, 2π) |
| `gamma` (γ) | polar angle of molecule centre position | uniform on the sphere |
| `alpha` (α) | azimuth of molecular axis | uniform OR cos²-weighted |
| `delta` (δ) | polar angle of molecular axis | uniform OR cos²-weighted |
| bond length | I-I distance with zero-point fluctuation | `R0_GS + N(0, deltaR0)` |

Position angles are always uniform on the sphere — molecules nucleate
randomly inside the droplet. Axis angles depend on the laser geometry:
in single-pulse mode the linearly polarised pump excites molecules
preferentially when their axis is parallel to the polarisation vector,
yielding a `cos²φ` weighting along the polarisation axis.

This module replaces the orientation-sampling block of
`vmi_sim_3d_neutral_propa_HeDFT_mimic.m` (lines ~225-275).

## Public API

```python
from i2_helium_md.sampling.orientations import sample_orientations

out = sample_orientations(
    num_molecules=2000,
    R0_GS_angstrom=9.0,        # 2.666 for ground state, 9.0 for HeDFT comparison
    deltaR0_angstrom=0.0,      # zero-point width (0 for fixed bond)
    anisotropic=True,          # True for single-pulse laser polarisation
    rng=np.random.default_rng(42),
)
out.alpha                       # axis azimuth, shape (N,)
out.delta                       # axis polar, shape (N,)
out.beta, out.gamma             # position angles
out.bond_length_angstrom        # I-I distance
```

Returns a frozen `MolecularOrientations` dataclass.

## What's inside

```
sample_orientations(...)
├─ _sample_uniform_sphere_angles()    # for (beta, gamma) always
├─ _sample_uniform_sphere_angles()    # for (alpha, delta) if isotropic
└─ _sample_anisotropic_axis_angles()  # for (alpha, delta) if anisotropic
```

### Uniform-sphere sampling

Uses inverse-CDF sampling: azimuth uniform in [0, 2π) and `cos(polar)`
uniform in [-1, 1]. The `cos(polar)` step is what makes the
distribution uniform on the sphere — naively sampling `polar ~ Uniform(0, π)`
oversamples the poles. The MATLAB code makes the same fix; an earlier
buggy version that didn't is preserved as a comment in the source.

### Anisotropic axis sampling (cos²φ rejection)

By symmetry, define the polarisation axis along x. A candidate axis
unit vector ``u = (cos α sin δ, sin α sin δ, cos δ)`` is accepted with
probability ``|u_x|² = (cos α sin δ)²``. Implemented by batched
rejection: the natural acceptance rate is 1/3, so we propose ~3.5N
candidates per round and loop until we have N acceptances. A
`max_oversample_factor=20` safety guard prevents infinite loops if
something is wrong with the RNG.

Reference for the cos²φ form: *Molecular reorientation during
dissociative multiphoton ionisation*, PRA, 1993.

### Bond length

Just `R0_GS + deltaR0 * standard_normal(N)`. With `deltaR0 = 0`
(default in all our input files) the bond length is exact. The
parameter is exposed for completeness — future work may use a nonzero
zero-point width.

## Statistical signatures (used as regression tests)

**Uniform on sphere:** `<cos²(polar)>` integrated over the sphere is
1/3. This holds for `gamma` always, and for `delta` in isotropic mode.

**Cos²-weighted:** `<u_x²>` over the sphere with weight `u_x²` is

```
<u_x²> = ∫ x² · x² dΩ / ∫ x² dΩ = (1/5) / (1/3) = 3/5 = 0.6
```

So the test `<(cos α · sin δ)²> ≈ 0.6` in anisotropic mode (vs ≈ 1/3
in isotropic mode) is the diagnostic that catches whether the
rejection sampling is wired correctly. With n=50 000 samples we get
`0.600 ± 0.005` in practice.

## Why the conversion to atomic xyz is NOT in this module

The MATLAB code computes atomic xyz coordinates by combining:

1. Molecule centre position from `(r0, beta, gamma)`
2. Axis offset from `(alpha, delta)` and bond length

That conversion uses the project's 2N array layout convention (atom 1
at indices 0..N-1, atom 2 at indices N..2N-1) which is owned by
`simulation/neutral.py`. Keeping the orientations sampler purely
"return angles and lengths" makes it reusable, smaller, and easier to
test without bringing in 2N-array layout conventions.

## Departures from MATLAB

1. **Returns a dataclass instead of five separate arrays.** Makes the
   contract explicit and guarantees consistent length.

2. **Validates inputs.** `num_molecules <= 0` and `deltaR0 < 0` raise
   immediately rather than producing nonsensical output later.

3. **Bounded rejection sampling.** MATLAB does `while length(alpha) < N`
   with no safety guard. We cap at 20 batches; for cos² with ~33%
   acceptance, even one batch covers the typical case.

4. **Single-shot batches instead of one-at-a-time.** MATLAB also
   batches via `rand(size(r0))` proposals, so this matches.
