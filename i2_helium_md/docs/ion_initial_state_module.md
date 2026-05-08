# The `ion_initial_state.py` module

## What problem does this file solve?

The neutral stage runs from t=0 (atoms placed inside droplets) to
t=t_neutral_max (atoms have moved a bit, possibly equilibrated).  At
that point we model the laser pulse: every atom suddenly becomes
singly charged. The ion stage then propagates this charged system
under Coulomb repulsion (between the two atoms of each molecule),
the ion-droplet potential, hard-sphere collisions with droplet
helium, and possible mass attachment.

`build_initial_ion_state` is the bridge: it takes the neutral
checkpoint as input, picks one column as the start state,
ionizes everything (charge=+1), allocates the trajectory arrays for
the ion stage, and computes the t=0 energies. The driver
(`run_ion_propagation` -- Step 11d) then fills in the rest.

This is the analogue of `build_initial_state` for the neutral stage.

## Position in the dependency chain

```
neutral_ckpt ──┐
               ▼
   build_initial_ion_state(cfg, neutral_ckpt, ...) ──► IonCheckpoint (column 0 only)
                                                                │
                                                                ▼
                                                run_ion_propagation fills cols 1..T-1
```

## Public API

```python
from i2_helium_md.simulation.ion_initial_state import build_initial_ion_state

ion = build_initial_ion_state(
    cfg,
    neutral_ckpt,
    *,
    num_steps_ion=2000,    # = ion_simulation_time / dt_ion
    start_id=-1,           # default: last column of neutral_ckpt
    rng=None,              # reserved; not yet used
)
```

## What gets read from the neutral checkpoint

- `positions_{x,y,z}[:, start_id]` — atom positions at the chosen start
- `velocities_{x,y,z}[:, start_id]` — atom velocities at the chosen start
- `mass_kg` — initial atom mass (per atom; iodine = 127 amu in our scope)
- `droplet_radii` — per-atom droplet radius

We don't read the energy columns from neutral; ion E_kin/E_pot use
*ion* potentials and start fresh.

## What we explicitly compute at t=0

- **`E_kin_eV[:, 0]`** — `½ m v²` per atom (in eV, with v in m/s = Å/ps × 100).
- **`E_pot_eV[:, 0]`** — sum of two terms:
  1. **Ion-droplet potential** with `binding_energy_I_ion_eV` (typically
     0.3 eV vs the neutral atom's 318 K ≈ 0.027 eV — much deeper well).
  2. **Half-pair Coulomb** from `partner_interaction_ion`, which already
     splits the pair Coulomb energy as half-per-atom in 2N layout.
- **`mass_history_kg[:, 0]`** — initial mass (no attachment yet).

All other trajectory arrays start at zero: `E_dissip_eV`,
`E_mass_attach_defect_eV`, `number_of_collisions`, and
`relative_loss_per_ps`. The schema-v5 `temperature_diagnostic`
field is allocated as a `(num_steps, 3)` array of NaN; only rows
where the driver actually observes a collision are overwritten.

## Bug fixes vs. legacy MATLAB

The legacy `vmi_sim_3d_ion_propa.m` has two t=0 bookkeeping bugs:

| Line | Legacy formula | Bug |
|------|---------------|-----|
| 289  | `E_kin_ion(:,1) = mass_i.*(vx² + vy²)²/2/eV` | Missing vz; `(...)²` squares v² → gives v⁴ (off by huge factor) |
| 291  | `E_pot_ion(:,1) = droplet_potential(sqrt(x² + y²) - R)` | Missing z in radius; missing partner Coulomb entirely |

These bugs are silent in production because `E_pot_ion(:,1)` is only
used for diagnostic plotting at the end of the run. The Python port
**fixes both** per project principle #10, and `test_ion_initial_state.py`
contains regression tests that catch each bug if it returns.

## Out-of-scope features (raise `NotImplementedError`)

The two production input scripts (`single_pulse_N2000.m` and
`single_pulse_droplet_distribution.m`) leave four legacy switches at
their default-disabled values. We refuse to build an ion state with
any of these enabled to avoid silent wrong physics:

- `effusive_dynamics`
- `single_charge_ionization_allowed`
- `additional_droplet_charges > 0`
- `highly_charged_iodine`

The check fires at build time, before any expensive ion stepping.

## Pump-probe note (`start_id`)

The legacy MATLAB outer loop runs an entire ion propagation for
**every** neutral column (= every probe delay). For our scope (single
pulse, neutral checkpoint stores t=0 and t=t_max only), only the LAST
column is physically meaningful — that's the end-state of the neutral
evolution which then gets ionized.

`start_id` accepts other values for forward-compatibility, but the
driver currently runs only one ion propagation (since pump-probe is
out of scope).
