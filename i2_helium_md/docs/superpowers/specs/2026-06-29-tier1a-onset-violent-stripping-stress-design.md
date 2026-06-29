# Tier-1a Onset-Violent Stripping Stress Design

## Purpose

The continuous-velocity Tier-1a anchored shedding result appears small. This
stress slice asks a deliberately stronger diagnostic question:

Can mass loss alone produce a visibly different trajectory if the shell is
stripped violently at onset?

This is not a new physical TDDFT-anchored model. It is a sensitivity test that
keeps the corrected continuous-velocity mass update and replaces only the shell
schedule with a deliberately severe one-event stripping schedule.

## Scope

In scope:

- Add a separate stress schedule family named `onset_strip`.
- Use one violent onset time: `t_strip_ps = 0.5`.
- Run four final-shell endpoints:
  - `n_final = 14`: lose 7 He at once, matching the current Tier-1a endpoint.
  - `n_final = 2`: keep only 2 He.
  - `n_final = 1`: keep only 1 He.
  - `n_final = 0`: bare iodine ion, lose all He.
- Keep continuous velocity at the event:
  `v_plus = v_minus`.
- Book the kinetic energy carried away by the removed co-moving He:
  `dE_mass_transfer = +0.5 * (n_before - n_after) * m_He * |v|^2`.
- Keep the existing fixed null available for comparison.
- Report these runs as stress diagnostics, not as the main Tier-1a anchored RMSE
  table.

Out of scope:

- Reintroducing the cold-shed velocity kick.
- Adding `E_int`, RRK, ladder, or five-term energy-gated evaporation.
- Changing Tier-0 drag coefficients, checkpoints, or the current physical
  `anchored_discrete` continuous-velocity baseline.
- Treating stress endpoints as TDDFT-derived or physically adjudicated.

## Architecture

The stress family should be implemented as a schedule-family extension, not as a
physics rewrite.

The current production Tier-1a path has three separable pieces:

- `physics/shell_schedule.py`: builds event schedules.
- `physics/mass_jump.py`: applies continuous-velocity mass changes and cold-shed
  diagnostic bounds.
- `simulation/ion_propagation_step.py::shed_step`: applies scheduled mass events
  before the BAOAB step.

The stress implementation should add a one-event batch schedule that the existing
driver can consume after a small generalization from "remove one He" to "remove
`n_before - n_after` He".

## Data Model

The existing `ShedEvent` already carries `n_before`, `n_after`,
`mass_before_amu`, and `mass_after_amu`. The stress schedule can reuse this shape
if downstream code computes the removed count as:

`n_removed = n_before - n_after`

For the existing anchored schedule this remains `1`. For stress schedules it is
larger, e.g. `21` for `n_final = 0`.

No checkpoint schema change is required. The existing `n_shell` field records the
integer shell count, and `mass_history_kg` records the realized mass.

## Schedule Behavior

Add a builder such as:

`build_onset_strip_schedule(t_strip_ps=0.5, n_final=0)`

It returns a schedule-compatible object with:

- one event at `t_strip_ps`;
- `n_before = 21`;
- `n_after = n_final`;
- `mass_before_amu = complex_mass_amu(21)`;
- `mass_after_amu = complex_mass_amu(n_final)`;
- `n_of_t(t) = 21` before the event and `n_final` at/after the event.

Validation:

- `t_strip_ps >= 0`;
- `n_final` integer in `[0, 20]`;
- reject `n_final >= 21` because no stripping would occur.

The existing anchored schedule remains unchanged.

## Physics And Ledger

The stress path uses the same continuous-velocity interpretation as the corrected
Tier-1a production path:

`v_plus = v_minus`

For a batch event removing more than one He:

`dE_mass_transfer = +0.5 * n_removed * m_He * |v|^2`

The post-event mass is:

`m_plus = m_minus - n_removed * m_He`

This conserves total momentum and kinetic energy of:

`remaining complex + removed co-moving He`

while the tracked complex alone loses the carried-away kinetic energy.

## Run And Reporting Surface

Add a separate stress generator and scorer surface, rather than mixing stress rows
into the default Tier-1a RMSE table.

Suggested run tags:

- `tier1a_stress_onset_strip_n14_t0.5`
- `tier1a_stress_onset_strip_n2_t0.5`
- `tier1a_stress_onset_strip_n1_t0.5`
- `tier1a_stress_onset_strip_n0_t0.5`

Suggested scripts:

- `scripts/gen_tier1a_stress_runs.py`
- `scripts/post_processing/tier1a_stress_table.py`

The stress table can reuse the same metric columns as the Tier-1a RMSE table:

- traceability columns;
- raw distance RMSE;
- same-smoothed I2 speed RMSE;
- speed mean ratio;
- ledger residual;
- `n_sheds` or `n_removed`;
- `n_shell_start`;
- `n_shell_end`.

The main plot should mirror the current clean visualization:

- fixed null;
- one selected stress case;
- HeDFT `|v2|`;
- CEEMDAN+SG `|v2|`;
- only `|v2|` traces.

## Testing

Unit tests:

- Stress schedule emits exactly one event at `0.5 ps`.
- `n_of_t` is `21` before onset and `n_final` at/after onset.
- Batch continuous-velocity shed leaves velocities unchanged.
- Batch mass drop equals `(21 - n_final) * m_He`.
- Batch mass-transfer energy equals `+0.5 * (21 - n_final) * m_He * |v|^2`.
- Momentum and KE of remaining complex plus removed co-moving He are conserved.

Driver tests:

- Existing anchored schedule still removes one He per event.
- Stress event removes multiple He in one step.
- `n_shell` jumps directly from 21 to the requested endpoint.
- Ledger closes with positive mass-transfer bookkeeping.
- Fixed run remains independent of stress schedule parameters.

Script tests:

- Stress run tags are distinct from physical Tier-1a tags.
- Stress scorer uses temporary run dirs, not production `data/runs`.
- Selected stress plot shows fixed plus one chosen stress case.

## Documentation

Update:

- `TIER1A_IMPLEMENTATION_PLAN.md` with a new stress slice marked diagnostic.
- `drag_migration_log_tier1a.md` with the delivery record once implemented.

The documentation must state that `onset_strip` is a sensitivity stress test, not
a TDDFT-anchored physical conclusion.
