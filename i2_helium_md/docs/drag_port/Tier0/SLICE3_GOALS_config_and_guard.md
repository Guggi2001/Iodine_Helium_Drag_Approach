# Slice 3 — Goals: `SimConfig` Drag Surface, Config-Load Guard, Coefficient Loader

**Status:** Specification. Config contract + validation logic + loader contract
only. No implementation code until `[PROCEED TO IMPLEMENTATION]`.

**Scope:** Slice 3 *only* — the declarative + validation layer of the drag
port. It adds the `SimConfig` drag fields (decisions doc §1–§6 surfaces), the
config-load `check_drag_config(cfg)` guard (§6.5 mass↔coefficient consistency
**and** §3.3 per-form dissipativity), and a thin JSON→`DragCoefficients`
loader in the presets/config layer. It consumes Slice 1's `DragCoefficients`
type and is consumed by Slice 4's driver.

**Slice 3 produces no behavioral change.** The collision path still runs (Slice
4 swaps it). The deliverable is: *a `SimConfig` that can express a valid Tier-0
drag setup and refuses an invalid one* — pure declaration + validation +
loading, verified entirely by config-construction and loader tests. The
hard-sphere path stays **bit-identical** for anyone not opting into a
drag-enabled preset.

**Why now (3 after 2).** Slice 2's stepper takes its parameters as plain
arguments; Slice 3 is what will let Slice 4's driver *source* those arguments
from `cfg`. It must land before Slice 4 (which reads the fields and calls the
loader) and after Slice 2 (which defined what parameters the stepper needs).

---

## 1. The one framing fact: this slice touches `config.py`, a frozen file

Unlike Slice 1 (new file) and Slice 2 (one surgical extraction), Slice 3 is
**inherently additive to `config.py`**, which the baseline lists under "avoid
changing." The slice is therefore disciplined to: *additive fields only, with
safe inert defaults, plus a separable validation function* — no edits to
existing fields, no change to existing `validate()` behaviour beyond invoking
the new guard, no change to existing presets.

---

## 2. The `SimConfig` drag fields (D1 — all declared now)

All ~18 fields are declared in this slice even though most are inert at Tier 0,
because the enums *are* the interchangeability design and declaring them is what
makes the cross-check apparatus real. The temporary dead surface is documented
in CLAUDE.md as a field→activating-slice table with an expiry (§8), **not**
blanket-exempted from the no-dead-code rule.

**Enum *types* are defined fully** (all members), even where some members are
unreachable at Tier 0 — the §6.5 guard's refusal logic is *about* the
non-`fixed` members, so they must exist as types to be referenced.

### 2.1 Tier-0-live fields (a Tier-0 path reads these)

| field | type / default | reader |
|---|---|---|
| `drag_form` | enum, default `linear_cubic` | guard (dissipativity); Slice 4 stepper build |
| `drag_coefficients` | `DragCoefficients` or `None`; default `None` | guard; Slice 4 `gamma_fn` build |
| `drag_spatial_gate` | enum, default `density_proportional` (collapses to erf-tied, §5.5) | Slice 4 gate build |
| `drag_gate_steepness` | float (Å), default = `potential_steepness` (14.2) | Slice 4 gate build |
| `mass_scenario` | enum, default `fixed` (**inert, not primary**) | guard |
| `m_eff_amu` | float (amu), default ≈ 202.954 | guard; loader cross-check; Slice 4 |
| `mass_initial_amu` | float (amu), default = `m_eff_amu` (A-option ii) | Slice 4 stepper `m` |
| `allow_inconsistent_mass_pairing` | bool, default `False` | guard (refuse→warn downgrade) |
| `drag_low_v_floor` | float (Å/ps); declared, **inert** for `linear_cubic` | (only `power_law`, deferred) |

### 2.2 Deferred fields (declared now, no Tier-0 reader; CLAUDE.md table records
the activating slice)

| field | type / inert default | activates |
|---|---|---|
| `noise_form` | enum, default `none` (**inert, not the primary `multiplicative_local_fdt`**) | Slice ≥4 / Tier 3 |
| `noise_calibration` | enum, default `hard_sphere_variance` (unread until noise on) | Tier 3 |
| `noise_geometry` | enum, default `longitudinal` (unread until noise on) | Tier 3 |
| `noise_low_v_behavior` | enum, default `vanish` (anisotropic-only) | Tier 3 |
| `mass_rate_form` | enum, default `density_only` (unread until scenario ≠ fixed) | Tier 1 |
| `mass_rate_coefficient` | float, default `0.0` | Tier 1 |
| `mass_relaxation_tau_ps` | float, default `0.0` (biphasic only) | Tier 1 |
| `helium_density_profile` | placeholder/`None` (promotes G4 beyond erf) | future |
| `validation_histogram_metric` | enum, default `wasserstein` | Tier 2 (no Tier-0 consumer — deferred, §C) |

**Default policy (A, confirmed):** every form-selector field defaults to its
**inert** member, not the design's "primary." The design's "primary" means
"first hypothesis to run," not "default when unspecified." A config left
untouched must do nothing surprising: `mass_scenario=fixed`, `noise_form=none`.

**`mass_initial_amu` default (A option ii, confirmed):** defaults to
`m_eff_amu`, the value correct for the only Tier-0-runnable scenario (`fixed`).
This avoids a `None`-resolution branch (resolve-from-scenario) that would itself
be untested-dead at Tier 0. Presets override per scenario when those scenarios
become runnable (127 under A, ~211 under B — §2.8 of the decisions doc).

---

## 3. `check_drag_config(cfg)` — the config-load guard (D2 + B)

A **separate function** (D2), called from `SimConfig.validate()`. Single
validation entry point preserved; guard logic separable and unit-testable
against constructed configs. Folds **two** checks (B):

### 3.1 Mass↔coefficient consistency (§6.5)

```
mass_scenario == fixed:
    require drag_coefficients.extraction_mass_model == "constant"
    require |extraction_mass_amu − m_eff_amu| ≤ TOL_amu        (§3.3 below)
    else → refuse (hard error) | warn if allow_inconsistent_mass_pairing
mass_scenario ∈ {scenario_A_accretion, scenario_B_stripping, biphasic}:
    require drag_coefficients.extraction_mass_model == "time_resolved"
    else → refuse | warn if allow_inconsistent_mass_pairing
```

At Tier 0 only the `fixed` branch is reachable (the evolving scenarios have no
coefficients to load yet), but the **whole branch set is written and tested
now** — that is the slice's reason for existing. Testing requires constructing
non-`fixed` configs with mismatched coefficients and asserting the refusal.
Not dead code: the guard *is* the deliverable that makes the §6.5 coupling
safe.

### 3.2 Per-form dissipativity (§3.3) — the live, non-vacuous branch at Tier 0

```
drag_form == linear_cubic:
    require a > 0
    turnover v† = √(−a/b): assert-and-skip (B, confirmed) —
        b > 0 ⇒ no real turnover ⇒ vacuously satisfied;
        assert no real turnover and RECORD that the max-trajectory-speed
        ceiling is unsourced, rather than inventing one (§7 open item).
drag_form == linear_quadratic:  require a > 0, c ≥ 0   (not reachable: NotImplemented upstream)
drag_form == threshold:         require F_sat > 0, v0 > 0   (not reachable)
drag_form == power_law:         require γ > 0   (not reachable)
```

**This is the only part of `check_drag_config` that exercises a live,
non-vacuous refusal on Tier-0-reachable input.** The mass-coupling branches
(§3.1) concern scenarios with no coefficients to load yet; the dissipativity
branch validates the *actual* `linear_cubic` coefficients in hand. The reserved
forms' dissipativity branches are written but unreachable (they raise
`NotImplementedError` in `drag.py` before a config carrying them could run);
they exist so the guard is complete when those forms are extracted.

### 3.3 The tolerance (D3 — hard-coded, named, documented)

`_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU = 8.0` — the 2-He end of the §6.5
"~1–2 He" band (≈ 2 × 4.0026). A *physical* statement (the drag curve is
mass-insensitive within this band, §6.6), **not** a user knob. Documented in
its definition with the §6.6 cross-reference.

> **Open sub-question (record, decide at implementation):** 1 He (~4 amu) vs.
> 2 He (~8 amu) as the band edge. 8.0 is the looser, safer-against-false-refuse
> choice; 4.0 is stricter. Trivially moot at Tier 0 (compares 202.954 to
> 202.954 → 0), but the test must construct a config differing by ~10 amu and
> assert refusal regardless of which edge, so the guard is not untested-green.

---

## 4. The coefficient loader (D4 + D5 + loader-case decision)

### 4.1 Signature and placement

A thin loader in the **presets/config layer** (D4) — `physics/` stays I/O-free
(rule 6). Case-agnostic, content-validating, provenance-enforcing:

```
load_drag_coefficients(coeff_dir: Path, *, expected_m_eff_amu: float)
    → DragCoefficients
```

- Reads `fit_parameters.json` from the given directory; constructs and returns
  a validated `DragCoefficients` (Slice 1 type, which self-validates arity).
- **The loader does not know or care about "9 Å vs 18 Å."** The case is the
  droplet geometry the *preset* already encodes; it is not a loader parameter,
  not an enum, not a user-facing key (loader-case decision, confirmed). The
  preset builds the path.

### 4.2 The case lives in the presets (D5 + loader-case)

- A named `REFERENCE_DRAG_ROOT` constant anchors the layout
  (`data/reference/drag/`) so the directory is not a scattered literal.
- Each drag-enabled preset owns its case→directory mapping because it *is* the
  geometry: the 9 Å drag preset builds
  `REFERENCE_DRAG_ROOT / "9A" / "linear_and_cubic"`, the 18 Å one builds the
  `18A` path. The case→dir fact is encoded **once** (in which preset you call),
  not duplicated in a parallel enum (rule 1).
- **No `DragCase` enum, no `drag_case_dir` helper.** Preset-free ad-hoc loading
  (e.g. a script comparing the two laws) constructs a path directly; the
  production path goes through presets.

### 4.3 Loader failure modes (D)

Refuse loudly, never silently proceed:

- `fit_parameters.json` missing → clear error naming the expected path.
- malformed / missing required keys (`a, b, a_err, b_err, meff_amu`) → error.
- **provenance disagreement:** JSON `meff_amu` ≠ `expected_m_eff_amu` → refuse.
  The loader **stamps `extraction_mass_amu` from the JSON** (single source of
  truth), so the bundle's provenance can never silently diverge from the file
  it came from. This makes the earlier-reviewed provenance concern a *load-time
  enforced* check rather than a trusted assumption.

### 4.4 Loader-tolerance vs. guard-tolerance — DISTINCT, do not collapse

Two mass comparisons exist and check *different* things; they need **separate
named constants and separate tests**:

- **Loader (§4.3):** JSON `meff_amu` vs. preset `expected_m_eff_amu`. A
  *provenance / plumbing identity* — the same number flowing two ways, should
  be **exactly** equal. Tolerance ~0 (exact match, or a float-eps guard);
  refuse on any real difference.
- **Guard (§3.3):** `extraction_mass_amu` vs. `m_eff_amu` under `fixed`. A
  *physics* tolerance — the ~8 amu mass-insensitivity band (§6.6).

At Tier 0 both compare 202.954 to 202.954 and both pass trivially, but folding
them into one check would look redundant and invite deletion of one. Keep
explicitly separate.

---

## 5. Presets (D5)

- `SimConfig` ships **inert defaults** (§2): `drag_form=linear_cubic`,
  `mass_scenario=fixed`, `noise_form=none`, `drag_coefficients=None`. Existing
  presets and the collision path are unaffected and remain bit-identical.
- A **new** drag-enabled preset (e.g. `single_pulse_N2000_drag`, 9 Å) is where
  coefficients are actually wired: it builds the `coeff_dir`, calls
  `load_drag_coefficients`, and sets `drag_coefficients` + `m_eff_amu`. The
  18 Å analog builds the `18A` path.
- Existing presets (`single_pulse_N2000`, etc.) are **not mutated in place** —
  the hard-sphere path for anyone not opting into drag stays untouched. (At end
  of Slice 3 the drag preset still runs the *collision* path too, since the
  driver swap is Slice 4; the drag preset's coefficients are loaded and
  validated but not yet consumed by a stepper.)

---

## 6. Acceptance criteria (test surface)

All tests are config-construction + loader tests; **no run, no trajectory.**

| Quantity | Expected | Tolerance |
|---|---|---|
| Inert defaults — a default `SimConfig` has `mass_scenario=fixed`, `noise_form=none`, `drag_coefficients=None`, `drag_form=linear_cubic` | enforced | exact |
| Bit-identical hard-sphere path — existing presets unchanged; a non-drag run is byte-identical to pre-Slice-3 | no behavioral change | exact / round-off |
| Guard pass — `fixed` + constant coeffs at matching `m_eff` | `check_drag_config` passes | exact |
| Guard refuse (mass band) — `fixed` + constant coeffs differing by ~10 amu (> TOL) | raises (hard error) | exact |
| Guard warn-not-refuse — same, with `allow_inconsistent_mass_pairing=True` | warns, does not raise | exact |
| Guard refuse (model mismatch) — non-`fixed` + constant coeffs | raises | exact |
| Guard dissipativity (live) — `linear_cubic` with `a ≤ 0` | raises | exact |
| Guard turnover assert-and-skip — `b > 0` ⇒ no real `v†`; passes with recorded note | passes | exact |
| Loader provenance — JSON `meff_amu` ≠ `expected_m_eff_amu` | raises | exact |
| Loader stamps from JSON — bundle `extraction_mass_amu` == JSON value, not preset value | enforced | exact |
| Loader missing/malformed JSON | raises clear path-naming error | exact |
| Loader → valid bundle — real 9 Å / 18 Å dirs produce arity-valid `DragCoefficients` | constructs | exact |
| Loader/guard tolerances separate — distinct named constants, distinct tests | enforced | exact |
| Enum completeness — all `mass_scenario` / `noise_*` / `drag_form` members defined | enforced | exact |

---

## 7. Open items surfaced (recorded, not blockers)

- **Mass-band tolerance edge: 1 He (4 amu) vs. 2 He (8 amu)** (§3.3). Proposed
  8.0; decide at implementation. Moot at Tier 0; the refusal test is
  edge-agnostic.
- **Max-trajectory-speed ceiling** for the `linear_cubic` turnover guard
  remains unsourced (§3.2 assert-and-skip). Dormant while `b > 0`; reactivated
  by any future `b < 0` re-extraction. Recorded, not sourced now.
- **Per-step closure rebuild reads these fields (Slice 4).** Slice 3 only makes
  the fields *available* and *valid*; Slice 4 sources them into the stepper
  build. Confirm the drag preset's loaded-but-unconsumed coefficients at end of
  Slice 3 are acceptable (they are — load-time validation is the point).

---

## 8. CLAUDE.md update required (D1 exception)

The declared-but-unread fields are a deliberate, time-limited exception to the
no-dead-code rule. CLAUDE.md gets a **field → activating-slice** table so the
exception has an expiry and a future rule-2 audit does not flag them:

- Tier-0-live (§2.1): read now or at Slice 4.
- Deferred (§2.2): `noise_*` → Tier 3; `mass_rate_*`, `mass_relaxation_tau_ps`
  → Tier 1; `helium_density_profile` → future G4; `validation_histogram_metric`
  → Tier 2.

The table names which slice removes each field from the exception.

---

## 9. Scope fence — what Slice 3 does NOT touch

- **No driver wiring / no behavioral change.** The collision path still runs;
  Slice 4 swaps it. Slice 3 fields are loaded and validated, not consumed by a
  stepper.
- **No `physics/` I/O.** The loader lives in the presets/config layer; `drag.py`
  and `baoab.py` stay pure.
- **No checkpoint changes / energy rename.** `IonCheckpoint` v6 comes with mass
  dynamics (Slice ≥4).
- **No active noise.** `noise_form=none` default; the fields are declared inert.
- **No mutation of existing presets.** New drag preset only.
- **No edits to existing `SimConfig` fields or `validate()` behaviour** beyond
  invoking `check_drag_config`.
- **No new physics.** Declaration + validation + loading only.

---

## 10. Definition of done

- All ~18 drag fields declared on `SimConfig` with inert defaults (§2); enum
  types complete.
- `check_drag_config(cfg)` implemented (mass↔coefficient consistency +
  per-form dissipativity), called from `validate()`; the turnover guard
  assert-and-skips with a recorded note.
- `_MASS_COEFFICIENT_CONSISTENCY_TOL_AMU` (≈ 8.0) named and documented;
  distinct from the loader's exact-match provenance check.
- `load_drag_coefficients(coeff_dir, *, expected_m_eff_amu)` in the
  presets/config layer: content-validating, stamps `extraction_mass_amu` from
  JSON, refuses on disagreement / missing / malformed.
- `REFERENCE_DRAG_ROOT` constant; case→dir mapping in the new drag preset(s);
  no enum, no helper.
- A new drag-enabled preset wires coefficients; existing presets and the
  hard-sphere path bit-identical.
- All §6 acceptance criteria green.
- CLAUDE.md field→activating-slice table added (§8).
- The §7 open items recorded in-repo.
