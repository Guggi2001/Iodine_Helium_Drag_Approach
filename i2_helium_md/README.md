# i2\_helium\_md

Modern Python port of Michael Stadlhofer's MATLAB molecular-dynamics code
for simulating iodine in helium nanodroplets (single-pulse, HeDFT-comparison
scope).

## Migration progress

|#|Python module|MATLAB source|Status|
|-|-|-|-|
|1|`physics/constants.py`|`physical\_constants.m`|✅ done|
|2|`config.py` (SimConfig)|\~36 MATLAB globals|✅ done|
|3|`presets.py`|`inputfiles\_dft\_comparison/single\_pulse\_N2000.m`|✅ done|
|4|`physics/potentials.py`|`droplet\_potential.m`, `get\_morse\_potential\_X.m`, `get\_morse\_potential\_I2plus.m`, `morse\_potential\_I2plus\_state\_select.m`|✅ done|
|5|`physics/interactions.py`|`atom\_interaction\_potential.m`, `ion\_interaction\_potential.m`, `add\_partner\_interaction.m`, `add\_partner\_interaction\_ion.m`|✅ done|
|6|`physics/leapfrog.py`|`frog\_step\_neutral.m`, `frog\_step\_ion.m`|✅ done|
|7|`sampling/droplet\_sizes.py`|`generate\_droplet\_sizes.m`, `get\_dropletsize.m`|✅ done|
|8|`sampling/radial\_positions.py`|`generate\_radial\_samples\_3d.m`|✅ done|
|9|`simulation/checkpoint.py`, `simulation/run\_directory.py`|`save('neutral\_propagation\_checkpoint', ...)`|✅ done|
|10|`simulation/neutral.py`, `sampling/orientations.py`, `physics/collisions.py`, `simulation/initial\_state.py`, `simulation/propagation\_step.py`|`vmi\_sim\_3d\_neutral\_propa\_HeDFT\_mimic.m`|✅ done|
|11|`simulation/ion.py`|`vmi\_sim\_3d\_ion\_propa.m`|🚧 11a + 11b + 11c done; 11d pending|
|12|`scripts/run\_single\_pulse.py`|`run\_simulation.m`|⏳|
|13|`postprocess/hedft\_loader.py` + `compare\_trajectories.py`|`simulation\_image\_only\_trajectories.m`|⏳|

## Documentation

* `docs/physics\_background.md` — physical model, potentials, and design rationale
* `docs/constants\_module.md` — walkthrough of the `constants.py` module
* `docs/interactions\_module.md` — walkthrough of the `interactions.py` module
* `docs/leapfrog\_module.md` — walkthrough of the `leapfrog.py` integrator
* `docs/droplet\_sizes\_module.md` — walkthrough of the `droplet\_sizes.py` sampler
* `docs/droplet\_sizes\_diagnostics\_module.md` — debugging plots for the pickup simulation
* `docs/radial\_positions\_module.md` — walkthrough of the `radial\_positions.py` sampler
* `docs/checkpoint\_module.md` — walkthrough of the `checkpoint.py` I/O module
* `docs/run\_directory\_module.md` — walkthrough of the `RunDirectory` convention layer
* `docs/orientations\_module.md` — walkthrough of the `orientations.py` angular sampler
* `docs/initial\_state\_module.md` — walkthrough of `build\_initial\_state` (Step 10c-i)
* `docs/propagation\_step\_module.md` — walkthrough of `neutral\_propagation\_step` (Step 10c-ii)
* `docs/neutral\_module.md` — walkthrough of `run\_neutral\_propagation` driver (Step 10c-iii)
* `docs/ion\_initial\_state\_module.md` — walkthrough of `build\_initial\_ion\_state` (Step 11b)
* `docs/ion\_propagation\_step\_module.md` — walkthrough of `ion\_propagation\_step` (Step 11c)
* `docs/collisions\_module.md` — walkthrough of the `collisions.py` hard-sphere physics
* `docs/migration\_log.md` — chronological record of decisions, deviations, and open questions

## Project layout

```
i2\_helium\_md\_py/
├── data/reference/              data files copied from legacy repo (see below)
├── docs/                        physics background + migration log
├── i2\_helium\_md/                Python package
│   ├── config.py                SimConfig dataclass
│   ├── presets.py               preset builders
│   ├── physics/                 constants, potentials, interactions, integrators
│   ├── sampling/                random samplers
│   ├── simulation/              neutral + ion propagation
│   └── postprocess/             HeDFT comparison plots
├── scripts/                     entry points
├── tests/                       smoke tests + pytest suite
└── pyproject.toml
```

## Data files needed in `data/reference/`

Copy these three files from the legacy MATLAB repo:

|Legacy path|→|New path|
|-|-|-|
|`HeDFT\_MD\_comparison\_neutral/custom\_start\_interpolating\_functions.mat`|→|`data/reference/hedft\_custom\_start.mat`|
|`single\_pulse\_simulation/HeDFT\_comparison/9Angstr\*ö\*m/data\_vabs2.csv`|→|`data/reference/hedft\_9A\_velocity.csv`|
|`single\_pulse\_simulation/HeDFT\_comparison/9Angstr\*ö\*m/R1-R2.csv`|→|`data/reference/hedft\_9A\_distance.csv`|

Once copied they're referenced via `SimConfig.data\_dir`, no hardcoded paths
anywhere in the code.

## Quickstart (after all steps done)

```python
from i2\_helium\_md import single\_pulse\_N2000
from i2\_helium\_md.simulation import run\_neutral, run\_ion

cfg = single\_pulse\_N2000(num\_molecules=500, seed=123)
neutral\_result = run\_neutral(cfg)
ion\_result = run\_ion(cfg, neutral\_result)
```

## Scope decisions (agreed with user)

* **In scope:** single-pulse neutral + ion dynamics; 9 Å HeDFT comparison
* **Out of scope:** pump-probe, effusive, VMI experimental comparison, 18 Å HeDFT
(no data in repo), Abel inversion, all image-processing utilities

## 

