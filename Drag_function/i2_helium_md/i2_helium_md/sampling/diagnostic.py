from i2_helium_md import single_pulse_N2000
from i2_helium_md.sampling.droplet_sizes import sample_droplet_sizes
import matplotlib.pyplot as plt

cfg = single_pulse_N2000(num_molecules=5000, seed=42, use_single_droplet_size=False)

N = sample_droplet_sizes(cfg, mode="post_pickup")    # shape (500,)
N_naive = sample_droplet_sizes(cfg, mode="raw")    # shape (500,)
# Create histogram
plt.figure(figsize=(10, 6))
plt.hist(N, bins=100, edgecolor='black', alpha=0.7)
plt.hist(N_naive, bins=100, edgecolor='black', alpha=0.7)
plt.xlabel('Droplet Size')
plt.ylabel('Frequency')
plt.title('Histogram of Droplet Sizes (Post-Pickup)')
plt.grid(axis='y', alpha=0.3)
plt.show()
