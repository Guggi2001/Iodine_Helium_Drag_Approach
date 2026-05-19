# The `hedft_loader.py` module — a walkthrough

## What problem does this file solve?

The legacy MATLAB pipeline imported HeDFT/TDDFT reference trajectories
through three separate ad-hoc importers, each with its own delimiter and
column convention:

```
9Angstroem/importfile_v2.m         (commas, 4 cols, time + |v|)
9Angstroem/importfile_R1_R2.m      (spaces, 2 cols, time + I-I distance)
18Angstroem/importfile_marti.m     (multi-file, separate position/velocity)
```

The Python port consolidates the reference data into a single normalised
8-column CSV per droplet size:

```
data/reference/9A_All_Data.csv
data/reference/18A_All_Data.csv
```

Header (identical for both):

```
Time_ps, V1_mag, V2_mag, V1_z, V2_z, V1_x, V2_x, R_distance
```

`hedft_loader.load_hedft_trajectory(path)` reads either file into a
:class:`HedftTrajectory` dataclass. One loader handles both droplet
sizes; downstream comparison code is droplet-size agnostic.

## Public API

```python
from i2_helium_md.postprocess import (
    HedftTrajectory, load_hedft_trajectory,
)

traj = load_hedft_trajectory("data/reference/9A_All_Data.csv")
# traj.time_ps         -> shape (T,)  picoseconds
# traj.distance_A      -> shape (T,)  angstrom
# traj.v1_magnitude_Aps, traj.v2_magnitude_Aps  -> angstrom/ps
# traj.v1_x_Aps, traj.v2_x_Aps                  -> angstrom/ps
# traj.v1_z_Aps, traj.v2_z_Aps                  -> angstrom/ps
# traj.droplet_radius_A  -> 9.0 (inferred from the filename prefix)
# traj.source_path       -> resolved Path to the file
```

The droplet radius is inferred from the filename stem (`9A_...` -> 9.0,
`18A_...` -> 18.0). When the prefix is missing, pass it explicitly:

```python
traj = load_hedft_trajectory(custom_path, droplet_radius_A=12.0)
```

## Contract

| Concern | Behaviour |
|---|---|
| Missing file | `FileNotFoundError` with the resolved absolute path |
| Wrong header | `ValueError` listing missing/unexpected columns |
| Non-monotonic `Time_ps` | `ValueError` |
| < 2 time samples | `ValueError` |
| No filename prefix and no `droplet_radius_A` | `ValueError` |
| Override radius | Explicit `droplet_radius_A=` wins over the prefix |

## Units

All units match the rest of the pipeline (and the MATLAB reference):

- time: picoseconds
- distance: angstrom
- velocity: angstrom/ps

No conversion happens at load time. The downstream
`compare_trajectories` module assumes this contract.

## Out of scope

The loader does I/O only. It does not smooth, interpolate, plot, or
derive quantities beyond what is in the file. Comparison logic lives in
:mod:`i2_helium_md.postprocess.compare_trajectories`. Project-wide scope
rules live in `CLAUDE.md`.
