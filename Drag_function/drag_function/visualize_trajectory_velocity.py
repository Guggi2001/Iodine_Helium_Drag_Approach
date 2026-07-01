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

scaling_factor_plot = 0.6
# Plotting 9A
fig9 = plt.figure(figsize=(12*scaling_factor_plot, 8*scaling_factor_plot))
fig9.suptitle('HeDFT Reference Results: 9Å', fontsize=16, fontweight='bold')
ax = fig9.add_subplot(111)

# Velocity Iodine 1
ax.plot(dict9['t'], dict9['v2'], color='#0072BD', lw=1.5, alpha=0.3)
ax.plot(vel_9_t, vel_9_cleaned, color='#A2142F', lw=1.5, label='Cleaned Velocity')
ax.set_title('Velocity: Top Iodine (9Å)')
ax.set_xlabel('t / ps')
ax.set_ylabel('v / Å/ps')
ax.axvline(2.67, color='#D95319', lw=1.5, label='Violent onset threshold', ls='--')
ax.grid(True, alpha=0.3)
ax.legend()

plt.show()

# Plotting 18A
fig18 = plt.figure(figsize=(12*scaling_factor_plot, 8*scaling_factor_plot))
fig18.suptitle('HeDFT Reference Results: 18Å', fontsize=16, fontweight='bold')
ax18 = fig18.add_subplot(111)

# Velocity Iodine 1
ax18.plot(dict18['t'], dict18['v2'], color='#0072BD', lw=1.5, alpha=0.3)
ax18.plot(vel_18_t, vel_18_cleaned, color='#A2142F', lw=1.5, label='Cleaned Velocity')
ax18.set_title('Velocity: Top Iodine (18Å)')
ax18.set_xlabel('t / ps')
ax18.set_ylabel('v / Å/ps')
ax18.axvline(4.54, color='#D95319', lw=1.5, label='Violent onset threshold', ls='--')
ax18.grid(True, alpha=0.3)
ax18.legend()

plt.show()
