import pandas as pd
import matplotlib.pyplot as plt
import os
from config_utils_local import config as C
from drag_function import io
import numpy as np

dict9 = io.load_data(C.PATH9A)
dict18 = io.load_data(C.PATH18A)

vel_9 = pd.read_csv(C.PATH_9A_CLEANED)
vel_9_t = np.array(vel_9['time'])
vel_9_cleaned = np.array(vel_9['cleaned_SG'])
vel_18 = pd.read_csv(C.PATH_18A_CLEANED)
vel_18_t = np.array(vel_18['time'])
vel_18_cleaned = np.array(vel_18['cleaned_SG'])

bubble_exit_9 = 9.0
bubble_exit_18 = 8.5

scaling_factor_plot = 0.4
# Create figure with two subplots side by side
fig, (ax, ax18) = plt.subplots(1, 2, figsize=(24*scaling_factor_plot, 8*scaling_factor_plot))

# Plotting 9A
# Velocity Iodine 1
ax.plot(dict9['t'], dict9['v2'], color='#0072BD', lw=1.5, alpha=0.3)
ax.plot(vel_9_t, vel_9_cleaned, color='#A2142F', lw=1.5, label='Cleaned Velocity')
ax.set_title('Velocity: Top Iodine (9Å)')
ax.set_xlabel('t / ps')
ax.set_ylabel('v / Å/ps')
ax.axvline(x= bubble_exit_9, color='green', linestyle='--', label='approximate bubble exit')
ax.axvline(2.67, color='#D95319', lw=1.5, label='Violent onset threshold', ls='--')
ax.grid(True, alpha=0.3)
ax.legend()

# Plotting 18A
# Velocity Iodine 1
ax18.plot(dict18['t'], dict18['v2'], color='#0072BD', lw=1.5, alpha=0.3)
ax18.plot(vel_18_t, vel_18_cleaned, color='#A2142F', lw=1.5, label='Cleaned Velocity')
ax18.set_title('Velocity: Top Iodine (18Å)')
ax18.axvline(x=bubble_exit_18, color='green', linestyle='--', label='approximate bubble exit')
ax18.set_xlabel('t / ps')
ax18.set_ylabel('v / Å/ps')
ax18.axvline(4.54, color='#D95319', lw=1.5, label='Violent onset threshold', ls='--')
ax18.grid(True, alpha=0.3)
ax18.legend()

plt.tight_layout()
plt.show()
