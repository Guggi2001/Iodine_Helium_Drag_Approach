from i2_helium_md.sampling.droplet_sizes_diagnostics import diagnose_pickup
from i2_helium_md import single_pulse_N2000

cfg = single_pulse_N2000(num_molecules=2000, p_source_mbar=40, T_source_K=14)
diagnostics, fig = diagnose_pickup(cfg, reduced_crosssection=False)
fig.show()

from i2_helium_md.sampling.droplet_sizes_diagnostics import plot_thesis_figure_3_2

fig = plot_thesis_figure_3_2()
fig.show()