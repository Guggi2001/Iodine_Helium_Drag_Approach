# Thesis figures

Every figure in this folder and the script that writes it. Run scripts from the
`i2_helium_md` project root, for example:

```powershell
python scripts/post_processing/tier2_b2_thesis_figure.py
python data/reference/scripts/plotting_histogram.py
```

Shared style (fonts, colours, text width, printed law names):
[`thesis_style.py`](../../scripts/post_processing/thesis_style.py).

Inputs under `data/runs/` are not in git, so those figures can only be regenerated
where the runs exist. Figures marked *committed data* need only files in this repo.

## Next step

**`tier0_B2_18A.png` and `tier0_B2_9A.png` are not yet cleanly connected to a
script.** [`tier0_drag_comparison.py`](../../scripts/post_processing/tier0_drag_comparison.py)
draws this figure type (one `CASE = "18A"` or `"9A"` per run), but:

- it has no `savefig`, so both PNGs were saved by hand from the plot window;
- the code has changed since: the panel title and the smoothed-curve colour no
  longer match the saved PNGs, so a rerun does not reproduce them;
- both PNGs are still in the default matplotlib style at 720 × 540 px.

To connect them, the script needs to save into `thesis_microscopic_calibration/`
and draw them in `thesis_style`. `TDDFT_Reference.png` already uses their line
styles (see `thesis_rest/`), so it should stay in step with any change.

## `thesis_microscopic_calibration/`

| File | Script | Inputs |
|---|---|---|
| [`tier0_B2_18A.png`](thesis_microscopic_calibration/tier0_B2_18A.png) | **Not cleanly connected, see [Next step](#next-step)** | `data/runs/`, tag `B2m21` |
| [`tier0_B2_9A.png`](thesis_microscopic_calibration/tier0_B2_9A.png) | **Not cleanly connected, see [Next step](#next-step)** | `data/runs/`, tag `B2m21` |
| [`tier0_B2_B3_h405_series.csv`](thesis_microscopic_calibration/tier0_B2_B3_h405_series.csv) | [`tier0_thesis_figure.py`](../../scripts/post_processing/tier0_thesis_figure.py) | Committed data (runs the Tier-0 integrations) |
| [`tier0_B2_B3_h405_scores.csv`](thesis_microscopic_calibration/tier0_B2_B3_h405_scores.csv) | [`tier0_thesis_figure.py`](../../scripts/post_processing/tier0_thesis_figure.py) | Committed data (runs the Tier-0 integrations) |
| [`tier0_ridge_m21_laws.png`](thesis_microscopic_calibration/tier0_ridge_m21_laws.png) | [`tier0_law_figure.py`](../../scripts/post_processing/tier0_law_figure.py) | Committed data |

## `thesis_results_discussion/`

| File | Script | Inputs |
|---|---|---|
| [`b2_1_size_distribution.png`](thesis_results_discussion/b2_1_size_distribution.png) | [`tier2_b2_thesis_figure.py`](../../scripts/post_processing/tier2_b2_thesis_figure.py) | `data/runs/B2_N6000_pooled` |
| [`b2_2_mean_ke.png`](thesis_results_discussion/b2_2_mean_ke.png) | [`tier2_b2_thesis_figure.py`](../../scripts/post_processing/tier2_b2_thesis_figure.py) | `data/runs/B2_N6000_pooled` |
| [`b2_2_mean_ke_split.png`](thesis_results_discussion/b2_2_mean_ke_split.png) | [`tier2_b2_thesis_figure.py`](../../scripts/post_processing/tier2_b2_thesis_figure.py) | `data/runs/B2_N6000_pooled` |
| [`b2_3_ked_shapes.png`](thesis_results_discussion/b2_3_ked_shapes.png) | [`tier2_b2_thesis_figure.py`](../../scripts/post_processing/tier2_b2_thesis_figure.py) | `data/runs/B2_N6000_pooled` |
| [`b2_4_ked_shapes_2d.png`](thesis_results_discussion/b2_4_ked_shapes_2d.png) | [`tier2_b2_thesis_figure.py`](../../scripts/post_processing/tier2_b2_thesis_figure.py) | `data/runs/B2_N6000_pooled` |
| [`laws_1_size_distribution_b3.png`](thesis_results_discussion/laws_1_size_distribution_b3.png) | [`tier2_law_comparison_figure.py`](../../scripts/post_processing/tier2_law_comparison_figure.py) | `data/runs/{B2,B3}_N6000_pooled` |
| [`laws_1_size_distribution_h405.png`](thesis_results_discussion/laws_1_size_distribution_h405.png) | [`tier2_law_comparison_figure.py`](../../scripts/post_processing/tier2_law_comparison_figure.py) | `data/runs/{B2,h405}_N6000_pooled` |
| [`laws_2_mean_ke_b3.png`](thesis_results_discussion/laws_2_mean_ke_b3.png) | [`tier2_law_comparison_figure.py`](../../scripts/post_processing/tier2_law_comparison_figure.py) | `data/runs/{B2,B3}_N6000_pooled` |
| [`laws_2_mean_ke_h405.png`](thesis_results_discussion/laws_2_mean_ke_h405.png) | [`tier2_law_comparison_figure.py`](../../scripts/post_processing/tier2_law_comparison_figure.py) | `data/runs/{B2,h405}_N6000_pooled` |
| [`laws_3_ked_shapes_b3.png`](thesis_results_discussion/laws_3_ked_shapes_b3.png) | [`tier2_law_comparison_figure.py`](../../scripts/post_processing/tier2_law_comparison_figure.py) | `data/runs/{B2,B3}_N6000_pooled` |
| [`laws_3_ked_shapes_h405.png`](thesis_results_discussion/laws_3_ked_shapes_h405.png) | [`tier2_law_comparison_figure.py`](../../scripts/post_processing/tier2_law_comparison_figure.py) | `data/runs/{B2,h405}_N6000_pooled` |
| [`laws_4_size_with_bare_lq.png`](thesis_results_discussion/laws_4_size_with_bare_lq.png) | [`tier2_law_comparison_figure.py`](../../scripts/post_processing/tier2_law_comparison_figure.py) | `data/runs/B2_N6000_pooled` |
| [`laws_5_ked_n012_lq.png`](thesis_results_discussion/laws_5_ked_n012_lq.png) | [`tier2_law_comparison_figure.py`](../../scripts/post_processing/tier2_law_comparison_figure.py) | `data/runs/B2_N6000_pooled` |
| [`g13_shed_timing.png`](thesis_results_discussion/g13_shed_timing.png) | [`tier2_g13_shed_timing.py`](../../scripts/post_processing/tier2_g13_shed_timing.py) | `data/runs/{law}_N1000_s1…s6` |
| [`g13_exit_vs_detection.png`](thesis_results_discussion/g13_exit_vs_detection.png) | [`tier2_g13_shed_timing.py`](../../scripts/post_processing/tier2_g13_shed_timing.py) | `data/runs/{law}_N1000_s1…s6` |
| [`g6_chord_map_b2_only.jpg`](thesis_results_discussion/g6_chord_map_b2_only.jpg) | [`tier2_g6_chord_map.py`](../../scripts/post_processing/tier2_g6_chord_map.py) | `data/runs/B2_N1000_s1…s6` |
| [`g6_chord_map_with_h405.jpg`](thesis_results_discussion/g6_chord_map_with_h405.jpg) | [`tier2_g6_chord_map.py`](../../scripts/post_processing/tier2_g6_chord_map.py) | `data/runs/{B2,h405}_N1000_s1…s6` |
| [`tier2_B2_B3.png`](thesis_results_discussion/tier2_B2_B3.png) | [`tier0_law_figure.py`](../../scripts/post_processing/tier0_law_figure.py) | Committed data |
| [`tier2_B2_h405.png`](thesis_results_discussion/tier2_B2_h405.png) | [`tier0_law_figure.py`](../../scripts/post_processing/tier0_law_figure.py) | Committed data |

## `thesis_rest/`

| File | Script | Inputs |
|---|---|---|
| [`I+He_n_Abundance.png`](thesis_rest/I+He_n_Abundance.png) | [`plotting_histogram.py`](../reference/scripts/plotting_histogram.py) | Committed data |
| [`I+He_n_Abundance_with_bare.png`](thesis_rest/I+He_n_Abundance_with_bare.png) | [`plotting_histogram.py`](../reference/scripts/plotting_histogram.py) | Committed data |
| [`TDDFT_Reference.png`](thesis_rest/TDDFT_Reference.png) | [`visualize_trajectory_velocity.py`](../reference/scripts/visualize_trajectory_velocity.py) | Committed data |
| [`experimental_reference_fragments.png`](thesis_rest/experimental_reference_fragments.png) | [`experimental_reference_fragments.py`](../reference/scripts/experimental_reference_fragments.py) | Committed data |
| [`experimental_reference_fragments_velocity.png`](thesis_rest/experimental_reference_fragments_velocity.png) | [`experimental_reference_fragments.py`](../reference/scripts/experimental_reference_fragments.py) | Committed data |
| [`ihe_ked_reference_mean_energy.png`](thesis_rest/ihe_ked_reference_mean_energy.png) | [`plot_ihe_ked_reference_mean_energy.py`](../../scripts/post_processing/plot_ihe_ked_reference_mean_energy.py) | Committed data |
| [`ihe_ked_reference_mean_energy_with_bare.png`](thesis_rest/ihe_ked_reference_mean_energy_with_bare.png) | [`plot_ihe_ked_reference_mean_energy.py`](../../scripts/post_processing/plot_ihe_ked_reference_mean_energy.py) | Committed data |

## `thesis_legacy/`

| File | Script | Inputs |
|---|---|---|
| [`fig_continuum_drag.png`](thesis_legacy/fig_continuum_drag.png) | [`legacy_continuum_drag_figure.py`](../../scripts/post_processing/legacy_continuum_drag_figure.py) | `data/runs/{9A,18A}_hedft_comparison` |
| [`fig_continuum_drag.pdf`](thesis_legacy/fig_continuum_drag.pdf) | [`legacy_continuum_drag_figure.py`](../../scripts/post_processing/legacy_continuum_drag_figure.py) | `data/runs/{9A,18A}_hedft_comparison` |
| [`fig_continuum_drag_points.csv`](thesis_legacy/fig_continuum_drag_points.csv) | [`legacy_continuum_drag_figure.py`](../../scripts/post_processing/legacy_continuum_drag_figure.py) | `data/runs/{9A,18A}_hedft_comparison` |
| [`legacy_ke_and_size.png`](thesis_legacy/legacy_ke_and_size.png) | [`legacy_ke_size_figure.py`](../../scripts/post_processing/legacy_ke_size_figure.py) | `data/runs/single_pulse_droplet_long` |
| [`legacy_ke_split_and_size.png`](thesis_legacy/legacy_ke_split_and_size.png) | [`legacy_ke_size_figure.py`](../../scripts/post_processing/legacy_ke_size_figure.py) | `data/runs/single_pulse_droplet_long` |
| [`legacy_ke_and_size_points.csv`](thesis_legacy/legacy_ke_and_size_points.csv) | [`legacy_ke_size_figure.py`](../../scripts/post_processing/legacy_ke_size_figure.py) | `data/runs/single_pulse_droplet_long` |
| [`legacy_ked_shapes.png`](thesis_legacy/legacy_ked_shapes.png) | [`legacy_ked_shapes_figure.py`](../../scripts/post_processing/legacy_ked_shapes_figure.py) | `data/runs/single_pulse_droplet_long` |
| [`legacy_ked_shapes_points.csv`](thesis_legacy/legacy_ked_shapes_points.csv) | [`legacy_ked_shapes_figure.py`](../../scripts/post_processing/legacy_ked_shapes_figure.py) | `data/runs/single_pulse_droplet_long` |
