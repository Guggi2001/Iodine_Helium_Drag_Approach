# Extraction scripts

One driver remains here: `method_b_g412_ridge.py`, the 18 Å linear + quadratic
ridge and its 9 Å exit bound at m = I⁺ + 21 He. Its modes are `ridge_m21` and
`bound_m21`, and it writes to `data/reference/drag/calibration_ridge/`. The heavy
lifting lives in `i2_helium_md/extraction/trajectory_matching.py`
(`docs/extraction/trajectory_matching_module.md`).

The June Method-A / Method-B drivers, the G4.9 / G4.11 / G4.15 refits and the
m = 19 ridge were run to completion and their verdicts are recorded in
`docs/drag_port/Tier0/METHOD_B_trajectory_matching_extraction.md` and
`docs/drag_port/Tier2/THESIS_FINAL_GOALS.md`. They were removed in the drag-tree
cleanup with the data they wrote; both are in git tag `pre-drag-cleanup`:

```powershell
git checkout pre-drag-cleanup -- scripts/extraction data/reference/drag
```
