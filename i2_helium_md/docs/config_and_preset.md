# The `config.py` and `presets.py` modules

## What problem do these files solve?

The legacy MATLAB code scattered ~36 `global` variables across
`run_simulation.m`, `physical_constants.m`, and a family of
`inputfiles_*/*.m` preset scripts. Two modules replace that:

- **`config.py`** — one strongly-typed `SimConfig` dataclass holding *every*
  tunable parameter, plus the validation that fails fast on a nonsensical
  configuration. Rule of thumb: every physical parameter lives here; nothing
  else in the codebase should carry a tunable magic number.
- **`presets.py`** — small builder functions, one per legacy input file, that
  return a `SimConfig` pre-filled for a named scenario. You start from a preset
  and override the few fields you need.

```python
from i2_helium_md import SimConfig, single_pulse_N2000
cfg = single_pulse_N2000(num_molecules=500, seed=123)
cfg.validate()
```

This document describes both modules in general **and** the **Slice 3**
drag-port additions they now carry: the `SimConfig` drag surface, the
config-load guard, and the coefficient loader.

Authoritative specs:
- `SLICE3_GOALS_config_and_guard.md` — this slice's contract.
- `DRAG_PORT_DESIGN_DECISIONS.md` — all architectural choices (§3.3, §5.5,
  §6.5, §6.6 are the ones Slice 3 implements).
- `docs/physics/drag_module.md` — Slice 1, the `DragCoefficients` type these
  modules consume.

## Position in the dependency chain

```
physics/constants.py        physics/drag.py  (DragCoefficients, form tags)
        ↓                          ↓
config.py  ←──────────────────────┘     SimConfig + DragForm/... aliases + check_drag_config
        ↓
presets.py  ←── data/reference/drag/<case>/   load_drag_coefficients + drag presets
        ↓
simulation drivers (neutral / ion stages), scripts, tests
```

`config.py` imports the `DragCoefficients` *type* and the form-tag constants
from `physics/drag.py`. There is **no import cycle**: `drag.py` depends only on
numpy/scipy, never on config. `presets.py` is the **only** drag-port file that
touches disk — `physics/` stays I/O-free (CLAUDE.md rule 6).

---

# Part 1 — `config.py`

## `SimConfig` — the parameter container

A single `@dataclass` whose defaults reproduce
`inputfiles_dft_comparison/single_pulse_N2000.m` combined with the constants
in `run_simulation.m`. Fields are grouped by concern (reproducibility, time
grid, laser, ensemble, droplet, collisions, ion, drag, output) and every
physical field encodes its unit in the name or a comment (`t_max_neutral` in
ps, `m_eff_amu` in amu, …).

### Derived quantities

A handful of `@property` accessors convert or compute and are **not** set by
hand: `num_timesteps_neutral`, `v_limit_angstrom_per_ps`,
`binding_energy_I_atom_eV`, `binding_energy_molecule_meV`, `E_min_eV`. They
encode the MATLAB conversion recipes (e.g. `E_min = (127·u)·v_limit²/2/eV`).

### `validate()` — fail fast

The single validation entry point. Existing checks: the MATLAB "all neutrals
will escape" warning when `E_min > binding_energy`, positive `num_molecules`
and timesteps, a legal `hard_sphere_collision_mode`. **Slice 3 adds one line at
the end** — a call to the separable `check_drag_config(self)` guard (below).
Nothing else in `validate()` changed.

## The drag surface (Slice 3, all declared now)

Slice 3 adds ~18 fields to `SimConfig`, **all with inert defaults**: a config
left untouched behaves exactly as before (no drag, hard-sphere path runs). The
full set is declared now — even fields no Tier-0 code reads yet — because the
enums *are* the interchangeability design (`DRAG_PORT_DESIGN_DECISIONS.md`),
and declaring them is what makes the cross-check apparatus real. The
temporary declared-but-unread fields are tracked in CLAUDE.md's
"Slice 3 declared-field exception" table (field → activating slice), so a
future no-dead-code audit does not flag them.

### Enum aliases — named `Literal` unions

Each enum is a module-scope `Literal` alias, matching the existing
`CollisionMode = Literal[1, 2, 3]` house style. A named alias is the single
source of truth for an enum's *members* (greppable, all members declared in one
place) without introducing `enum.Enum`. Field defaults are the bare inert
string values.

| alias | members | default (inert) |
|---|---|---|
| `DragForm` | `linear_cubic`, `linear_quadratic`, `threshold`, `power_law` | `linear_cubic` |
| `DragSpatialGate` | `density_proportional`, `erf_tied`, `erf_independent`, `sharp` | `density_proportional` |
| `MassScenario` | `fixed`, `biphasic`, `anchored_discrete` | `fixed` |
| `NoiseForm` | `none`, `multiplicative_local_fdt`, `empirical_residual` | `none` |
| `NoiseCalibration` | `hard_sphere_variance`, `tddft_residual`, `strict_fdt_bath` | `hard_sphere_variance` |
| `NoiseGeometry` | `longitudinal`, `isotropic`, `anisotropic` | `longitudinal` |
| `NoiseLowVBehavior` | `vanish`, `blend_to_isotropic` | `vanish` |
| `PickupRateForm` | `density_only`, `sweeping`, `dwell_time` | `density_only` |
| `PickupOccupancyCap` | `langmuir`, `none` | `langmuir` |
| `HeCaptureVelocity` | `at_rest`, `thermal` | `at_rest` |
| `ValidationHistogramMetric` | `wasserstein`, `chi2`, `ks` | `wasserstein` |

> **Default policy:** every form-selector defaults to its **inert** member, not
> the design's "primary." The design's "primary" means "first hypothesis to
> run," not "default when unspecified." So `mass_scenario=fixed` (not
> `anchored_discrete`) and `noise_form=none` (not `multiplicative_local_fdt`)
> by default.

> **The `Literal`-over-`enum.Enum` tradeoff and its recovery.** A `Literal`
> alias is **not enforced at runtime** — Python will happily store
> `drag_form="linaer_cubic"` (a typo) without complaint, where an `enum.Enum`
> would have rejected it. That static typo-catching is recovered at runtime by
> the guard's **reject arms**: `check_drag_config` rejects any
> `drag_form ∉ _KNOWN_DRAG_FORMS` and any `mass_scenario` outside its members
> (see the guard section). So the only structural advantage `enum.Enum` had is
> bought back without the new pattern.

### The fields

**Tier-0-live** — read by the guard now and/or the Slice 4 stepper:

| field | type / default | role |
|---|---|---|
| `drag_form` | `DragForm` = `linear_cubic` | guard (dissipativity); Slice 4 stepper build |
| `drag_coefficients` | `Optional[DragCoefficients]` = `None` | guard; Slice 4 γ build |
| `drag_spatial_gate` | `DragSpatialGate` = `density_proportional` | Slice 4 gate build (collapses to erf-tied, §5.5) |
| `drag_gate_steepness` | `float` (Å) = `14.2` (= `potential_steepness`) | Slice 4 gate width |
| `mass_scenario` | `MassScenario` = `fixed` | guard |
| `m_eff_amu` | `float` (amu) = `202.953908` | guard; loader cross-check; Slice 4 |
| `mass_initial_amu` | `float` (amu) = `202.953908` (= `m_eff_amu`) | Slice 4 stepper mass |
| `allow_inconsistent_mass_pairing` | `bool` = `False` | guard refuse→warn downgrade |
| `drag_low_v_floor` | `float` (Å/ps) = `0.0` | inert for `linear_cubic`; only `power_law` n<0 |

> **`drag_low_v_floor = 0.0` is inert for now and carries no meaning.** It is the
> low-velocity floor that regularises a `power_law` drag with `n < 0` (where
> `γ → ∞` as `v → 0`). The realised `linear_cubic` form is finite at rest
> (`γ → g·a`) and never reads it, and the real `power_law` export has `n ≈ +2`
> (also regular at `v=0`), so the `0.0` default is a placeholder — it is neither
> consumed nor validated by anything at Tier 0. It activates only if a
> hypothetical `n < 0` `power_law` re-extraction is ever wired in.

**Deferred** — declared now, no Tier-0 reader; activated by a later slice:

`noise_form`, `noise_calibration`, `noise_geometry`, `noise_low_v_behavior`
(Tier 3); `pickup_rate_coefficient` (λ₀), `pickup_occupancy_exponent` (p)
(Tier-2 Phase-B Slice P; read by the Phase-C driver — the enum selectors
`pickup_rate_form` / `pickup_occupancy_cap` / `he_capture_velocity` and
`helium_density_profile` are **live** at their slice via the config-load guards);
`validation_histogram_metric` (Tier 2).

> `mass_initial_amu` defaults to `m_eff_amu` (the value correct for the only
> Tier-0-runnable scenario, `fixed`) rather than `None`, to avoid an
> untested-dead "resolve initial mass from scenario" branch. Presets override it
> per scenario once those scenarios become runnable.

## `check_drag_config(cfg)` — the config-load guard

A **separate module-level function** (not a method), called from
`validate()`. Keeping it separable makes it unit-testable against constructed
configs without running a simulation.

### 0. `drag_form` typo-recovery + consistency — runs first

Before anything else, and **unconditionally** (before the
`drag_coefficients is None` early return), the guard rejects an unrecognised
`drag_form`:

```
cfg.drag_form ∉ _KNOWN_DRAG_FORMS   →   raise ValueError("unknown drag_form ...")
```

`_KNOWN_DRAG_FORMS` reuses the `drag.py` form tags (no duplicate string
literals). This is the runtime recovery for the typo-catching that `Literal`
gives up (above); it mirrors the `else: raise` arm the `mass_scenario` dispatch
already had. It fires even on a non-drag config, so a typo never passes
silently; the default `linear_cubic` passes unchanged.

Then the guard **no-ops when `drag_coefficients is None`** (the inert default),
so existing/non-drag configs are otherwise unaffected. When coefficients **are**
present, it first cross-checks form agreement:

```
coeffs.form ≠ cfg.drag_form          →   raise ValueError("...does not match coefficient form...")
```

The Slice 4 stepper builds from `cfg.drag_form` while consuming `coeffs`, so the
two must name the same form. (Note `coeffs.form` is itself construction-validated
by `DragCoefficients`, so the two checks together mean the dissipativity
dispatch below is keyed on a value proven to be a known member — hence its own
`else: raise` is defensive `# pragma: no cover`.) It then folds the two
substantive checks.

### 1. Per-form dissipativity (§3.3) — the only live Tier-0 refusal

A drag force must oppose motion at all operating speeds. For `linear_cubic`:

```
require a > 0                              # dissipative at low v
turnover  v† = √(−a/b)  is real only if b < 0
  b > 0 ⇒ no real turnover ⇒ vacuously dissipative everywhere  (assert-and-skip)
```

This is the **only** branch that exercises a live, non-vacuous refusal on
Tier-0-reachable input — it validates the actual coefficients in hand. The
"max trajectory speed" a `b<0` turnover would need to be tested against is
**unsourced** (`DRAG_PORT_DESIGN_DECISIONS.md` §6.10 open item) and is
deliberately *not* invented; it stays dormant while both extracted cases have
`b>0`. The reserved forms' branches (`linear_quadratic`: `a>0,c≥0`;
`threshold`: `F_sat>0,v0>0`; `power_law`: `γ>0`) are written for completeness
but **unreachable** — `physics/drag.py` raises `NotImplementedError` for those
forms upstream before a config carrying them could run.

### 2. Mass ↔ coefficient consistency (§6.5)

The extraction force balance `F_drag = m(t)·a − F_C` constrains only the
*combination* `m·a`. A drag law fit at constant `m_eff` is self-consistent
**only** when re-applied at constant `m_eff`; an evolving-mass scenario needs
coefficients re-extracted under that scenario's `m(t)`. So `mass_scenario` and
`drag_coefficients` are a **coupled pair**, enforced here:

```
mass_scenario == fixed:
    require extraction_mass_model == "constant"
    require |extraction_mass_amu − m_eff_amu| ≤ _MASS_COEFFICIENT_CONSISTENCY_TOL_AMU
mass_scenario ∈ {biphasic, anchored_discrete}:
    require extraction_mass_model == "time_resolved"
inconsistent  →  raise ValueError
              →  (or warnings.warn, if allow_inconsistent_mass_pairing=True)
```

At Tier 0 only the `fixed` branch is reachable **through the production path**
(the evolving scenarios have no coefficients to load yet), but the whole branch
set is written and tested now — that is the slice's reason for existing.

> **The non-`fixed` branch is exercised by synthetic construction, not left
> written-but-untested-green.** Because no loader produces evolving-scenario
> coefficients yet, the only way to reach the non-`fixed` arm before Slice 4 is
> to hand-build a config: an evolving `mass_scenario` paired with a `constant`
> bundle (`_constant_coeffs()` in the tests). Both forks are asserted:
> - **refuse** — `mass_scenario="anchored_discrete"` + constant coeffs →
>   `ValueError` (matches `"time_resolved"`);
> - **override → warn** — `mass_scenario="biphasic"` + constant
>   coeffs + `allow_inconsistent_mass_pairing=True` → `RuntimeWarning`, no raise.
>
> The two cover both evolving members. (Tier-1a's `anchored_discrete` run takes
> the override→warn arm in production, on the §6.6 mid-window defence.) The dispatch's final `else: raise "unknown mass_scenario"` is the
> `mass_scenario` analogue of the `drag_form` typo arm; it is unreachable given
> the `Literal` members and is marked `# pragma: no cover`.

### `_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU = 8.0`

A **named module-level constant**, the 2-He edge of the §6.5 "~1–2 He" band
(2 × 4.0026). It is a *physical* statement — the drag curve is mass-insensitive
within this band (§6.6) — **not** a user knob. (1-He / 4.0 amu was the stricter
alternative; 8.0 is the looser, safer-against-false-refuse choice. Moot at
Tier 0, where it compares 202.954 to 202.954.)

> **This is NOT the loader's tolerance.** The guard's 8-amu band is a *physics*
> tolerance. The loader (Part 2) has a separate **exact-match** provenance
> check. They check different things and are kept as distinct named constants
> with distinct tests — see "Two tolerances" below.

---

# Part 2 — `presets.py`

## Preset builders

Each builder mirrors one legacy `inputfiles_*/*.m` and returns a `SimConfig`
via `SimConfig(...)` + `dataclasses.replace(cfg, **overrides)`. The existing
three are unchanged by Slice 3:

- `single_pulse_N2000` — canonical 9 Å He-DFT comparison, 2000 molecules.
- `single_pulse_N2000_18Angst` — 18 Å variant (larger I–I distance, weaker I⁺
  cross-section, weaker ion binding, lower attachment probability).
- `single_pulse_droplet_distribution` — source-condition droplet-size sampler,
  8000 molecules, ground-state I₂ distance.

## `REFERENCE_DRAG_ROOT` — the data anchor

```python
REFERENCE_DRAG_ROOT = Path(__file__).resolve().parents[1] / "data" / "reference" / "drag"
```

`parents[1]` is the repo root (this file is `i2_helium_md/i2_helium_md/`).
Encoded **once** here so the `data/reference/drag/` layout is not a scattered
literal. There is deliberately **no `DragCase` enum and no `drag_case_dir`
helper**: the case (9 Å vs 18 Å) is the droplet geometry the *preset* already
encodes, so each drag preset builds its own `REFERENCE_DRAG_ROOT/<case>/...`
path — the case→dir fact lives once, in which preset you call.

## `load_drag_coefficients(coeff_dir, *, expected_m_eff_amu)` → `DragCoefficients`

A thin, content-validating, provenance-enforcing loader. It lives here (not in
`physics/`) to keep the physics layer I/O-free. It is **case-agnostic** — it
takes a directory, not a "9 Å / 18 Å" key.

Steps:
1. Read `coeff_dir/fit_parameters.json` (`{a, b, a_err, b_err, meff_amu}`).
2. Refuse on any failure mode, loudly: missing file (`FileNotFoundError`
   naming the path), malformed JSON, missing required keys, or **provenance
   mismatch** (`ValueError`).
3. Construct and return a `DragCoefficients` with `form="linear_cubic"`,
   `coefficients={"a", "b"}`, `extraction_mass_model="constant"`, and
   `extraction_mass_amu` **stamped from the JSON** — the single source of truth,
   so the bundle's provenance can never silently diverge from the file. The
   bundle's own `__post_init__` re-validates arity.

### Provenance check — exact match

```
abs(json_meff − expected_m_eff_amu) > 1e-6   →   refuse
```

`expected_m_eff_amu` is what the *caller* (the preset) expects this case to
carry. This is a **plumbing/provenance identity** — the same number flowing two
ways should be exactly equal; any real difference is a wiring bug. It enforces
at load time the concern that the stamped extraction mass really is the mass
the fit ran under.

## Drag-enabled presets

Two new presets wire the real coefficients (the only place coefficients are
actually loaded in production):

- `single_pulse_N2000_drag` — builds on `single_pulse_N2000`, loads the **9 Å**
  `linear_and_cubic` coefficients, sets `drag_coefficients`, `m_eff_amu`,
  `mass_initial_amu`.
- `single_pulse_N2000_18Angst_drag` — the **18 Å** analog.

Both pin `mass_scenario="fixed"`, so `validate()` → `check_drag_config` passes
(constant coeffs at the matching `m_eff`). **No behavioral change in Slice 3:**
the hard-sphere collision fields are inherited unchanged from the base presets,
so the collision path still runs — the loaded coefficients are validated but
**not yet consumed by any stepper** (that is Slice 4). Existing presets are not
mutated in place; the drag presets are *additional*.

## Two tolerances — distinct, do not collapse

| | where | what it compares | tolerance |
|---|---|---|---|
| **Loader provenance** | `presets.py` | JSON `meff_amu` vs. preset `expected_m_eff_amu` | exact (1e-6) |
| **Guard mass band** | `config.py` `_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU` | `extraction_mass_amu` vs. `m_eff_amu` under `fixed` | ~8 amu (physics) |

At Tier 0 both compare 202.954 to 202.954 and pass trivially, but they answer
different questions (is the file the one I asked for? vs. is the law's mass
close enough to the run's mass to be physically valid?). Folding them into one
would invite deleting the wrong one. Kept separate with separate tests.

---

## What's NOT here (Slice 3 scope fence)

- **No driver wiring / no behavioral change.** The collision path still runs;
  Slice 4 swaps it. Loaded coefficients are validated, not consumed.
- **No `physics/` I/O.** The loader lives here; `drag.py` / `baoab.py` stay
  pure.
- **No checkpoint changes / energy rename.** `IonCheckpoint` v6 comes with mass
  dynamics (Slice ≥4).
- **No active noise.** `noise_form="none"` default; the noise fields are inert.
- **No mutation of existing presets** or existing `SimConfig` fields, and no
  change to `validate()` beyond the one `check_drag_config` call.

## Regression-test signatures

Locked in by `tests/test_drag_config.py` (config-construction + loader only,
**no run, no trajectory**):

| Quantity | Expected | Tolerance |
|---|---|---|
| Inert defaults — default `SimConfig` has `mass_scenario=fixed`, `noise_form=none`, `drag_coefficients=None`, `drag_form=linear_cubic` | enforced | exact |
| `mass_initial_amu` defaults to `m_eff_amu`; `drag_gate_steepness` to `potential_steepness` | enforced | exact |
| Default config validates (guard no-ops when no coefficients) | passes | — |
| Existing presets carry no drag coefficients; core hard-sphere fields unchanged | enforced | exact |
| Guard pass — `fixed` + constant coeffs at / within band of `m_eff` | passes | exact |
| Guard refuse — `fixed` + constant coeffs > 8 amu off | raises `ValueError` | exact |
| Guard warn-not-refuse — same, `allow_inconsistent_mass_pairing=True` | warns `RuntimeWarning` | exact |
| Guard refuse — evolving scenario + constant coeffs (synthetic) | raises `ValueError` ("time_resolved") | exact |
| Guard override→warn — evolving scenario + constant coeffs + `allow_inconsistent_mass_pairing` | warns `RuntimeWarning`, no raise | exact |
| Guard dissipativity — `linear_cubic` with `a ≤ 0` | raises `ValueError` | exact |
| Guard turnover assert-and-skip — `b > 0` ⇒ no real `v†` | passes | exact |
| Guard typo-reject — unrecognised `drag_form` (default config, no coeffs) | raises `ValueError` ("unknown drag_form") | exact |
| Guard typo-reject — unrecognised `drag_form` with coefficients / via `validate()` | raises | exact |
| Guard form consistency — `drag_form` ≠ `coeffs.form` | raises ("does not match coefficient form") | exact |
| Guard form match — member `drag_form` + matching coeffs | passes | exact |
| Loader → valid bundle from real 9 Å / 18 Å dirs | constructs | — |
| Loader stamps `extraction_mass_amu` from JSON | == JSON value | exact |
| Loader provenance mismatch (JSON ≠ expected) | raises `ValueError` | exact |
| Loader missing file / malformed JSON / missing keys | raises (path-naming) | exact |
| Loader vs. guard tolerances distinct — 0.001 amu passes guard but fails loader | both behaviours | exact |
| Drag presets validate and carry `linear_cubic` coefficients; keep hard-sphere fields | enforced | exact |
| Enum completeness — all members present on each `Literal` alias | enforced | exact |
