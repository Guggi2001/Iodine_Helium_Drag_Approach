import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.integrate import cumulative_trapezoid
# Path Configuration
office = False
if office:
    base_path = r"T:\NextCloud_PaulGuggenbichler\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT\9A"
else:
    base_path = r"C:\Users\paulg\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT\9A"

# Load 9A Data (Comma Separated)
df_v1 = pd.read_csv(os.path.join(base_path, "data_vabs.csv"), header=None)
df_v2 = pd.read_csv(os.path.join(base_path, "data_vabs2.csv"), header=None)
df_R  = pd.read_csv(os.path.join(base_path, "R1-R2.csv"), header=None)

# Plotting 9A
fig9, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
fig9.suptitle('HeDFT Reference Results: 9Å', fontsize=16, fontweight='bold')

# Velocity Iodine 1
axs[0, 0].plot(df_v1[0], df_v1[1].abs(), color='#0072BD', lw=1.5)
axs[0, 0].set_title('Velocity: Iodine 1 (9Å)')
axs[0, 0].set_ylabel('v / Å/ps')
axs[0, 0].set_xlim(-2, 12)
axs[0, 0].grid(True, alpha=0.3)

# Velocity Iodine 2
axs[0, 1].plot(df_v2[0], df_v2[1].abs(), color='#D95319', lw=1.5)
axs[0, 1].set_title('Velocity: Iodine 2 (9Å)')
axs[0, 1].set_xlim(-2, 12)
axs[0, 1].grid(True, alpha=0.3)

# Internuclear Distance (spanning bottom)
gs = axs[1, 0].get_gridspec()
for ax in axs[1, :]: ax.remove()
ax_bottom = fig9.add_subplot(gs[1, :])
ax_bottom.plot(df_R[0], df_R[1], color='black', lw=2)
ax_bottom.set_title('Internuclear Distance R (9Å)')
ax_bottom.set_xlabel('t / ps')
ax_bottom.set_ylabel('R / Å')
ax_bottom.set_xlim(0, 5) # Matches your MATLAB screenshot zoom
ax_bottom.grid(True, alpha=0.3)

plt.show()
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
# Export Cleaned Data to CSV
# -------------------------------------------------------------------------------------
# --- Configuration ---
R0 = 9.0
output_file = os.path.join(base_path, "9A_All_Data_old.csv")

# --- 1. Load Data ---
df1 = pd.read_csv(os.path.join(base_path, "data_vabs.csv"), header=None)
df2 = pd.read_csv(os.path.join(base_path, "data_vabs2.csv"), header=None)
df_actual = pd.read_csv(os.path.join(base_path, "R1-R2.csv"), header=None)

# Raw data
t1, v1_raw = df1[0].values, df1[1].abs().values
t2, v2_raw = df2[0].values, df2[1].abs().values

# --- 2. Master Grid & Interpolation ---
# Create a high-res grid spanning the available data
t_master = np.linspace(0, max(t1.max(), t2.max()), 10000)

v1_interp = np.interp(t_master, t1, v1_raw)
v2_interp = np.interp(t_master, t2, v2_raw)

# --- 3. Integration ---
v_rel = v1_interp + v2_interp
R_recon = R0 + cumulative_trapezoid(v_rel, t_master, initial=0)

# --- 4. Clean CSV Export ---
export_df = pd.DataFrame({
    'Time_ps': t_master,
    'V1_mag': v1_interp,
    'V2_mag': v2_interp,
    'R_distance': R_recon
})
export_df.to_csv(output_file, index=False)
print(f"Success! Clean distance data exported to: {output_file}")

# --- 5. Plotting 9A-Verification & Interpolated Velocities ---

fig9, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
fig9.suptitle('HeDFT Reference Results: 9Å', fontsize=16, fontweight='bold')

# Velocity Iodine 1
axs[0, 0].plot(t_master, v1_interp, color='#0072BD', lw=1.5)
axs[0, 0].set_title('Interpolated Velocity: Iodine 1 (9Å)')
axs[0, 0].set_ylabel('v / Å/ps')
axs[0, 0].set_xlim(-2, 12)
axs[0, 0].grid(True, alpha=0.3)

# Velocity Iodine 2
axs[0, 1].plot(t_master, v2_interp, color='#D95319', lw=1.5)
axs[0, 1].set_title('Interpolated Velocity: Iodine 2 (9Å)')
axs[0, 1].set_xlim(-2, 12)
axs[0, 1].grid(True, alpha=0.3)

# Internuclear Distance (spanning bottom)
gs = axs[1, 0].get_gridspec()
for ax in axs[1, :]: ax.remove()
ax_bottom = fig9.add_subplot(gs[1, :])
ax_bottom.plot(t_master, R_recon, color='black', lw=2)
ax_bottom.plot(df_R[0], df_R[1], 'r-', lw=2, label='Actual Data (from files)')
ax_bottom.set_title('9Å Case: Distance Verification & Extension')
ax_bottom.set_xlabel('t / ps')
ax_bottom.set_ylabel('R / Å')
ax_bottom.set_xlim(0, 12) # Matches your MATLAB screenshot zoom
ax_bottom.grid(True, alpha=0.3)

plt.show()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
# 9 A case new: Load, Plotting & Exporting
# -------------------------------------------------------------------------------------

def reconstruct_R_from_v(v_w, t_w, R0):
    """
    Reconstruct R from smoothed velocity data using Savitzky-Golay filter.
    """
    return (R0 + cumulative_trapezoid(v_w, t_w, initial=0))


# --- Configuration ---
output_file = os.path.join(base_path, "9A_All_Data.csv")
step = 100  # Take every 100th point

def marti_loader(path):
    """Robust loader for large Marti files with scientific notation."""
    return np.genfromtxt(path, comments='#', invalid_raise=False)

# --- 1. Load Original Data ---
print("Loading large Marti files...")
v1_raw = marti_loader(os.path.join(base_path, "vimp.I_1"))
v2_raw = marti_loader(os.path.join(base_path, "vimp.I_2"))
r1_raw = marti_loader(os.path.join(base_path, "rimp.I_1"))
r2_raw = marti_loader(os.path.join(base_path, "rimp.I_2"))

# --- 2. Calculate Magnitudes and Distance ---
# Calculate R (Distance) from positions
min_len_r = min(len(r1_raw), len(r2_raw))
t_master = r1_raw[:min_len_r, 0]
R_actual = np.linalg.norm(r1_raw[:min_len_r, 1:4] - r2_raw[:min_len_r, 1:4], axis=1)

# Calculate V1 Magnitude
v1_mag = np.linalg.norm(v1_raw[:min_len_r, 1:4], axis=1)
v2_mag = np.linalg.norm(v2_raw[:min_len_r, 1:4], axis=1)
v1_z = v1_raw[:min_len_r, 3]
v2_z = v2_raw[:min_len_r, 3]
v1_x = v1_raw[:min_len_r, 1]
v2_x = v2_raw[:min_len_r, 1]


R_reconstructed = reconstruct_R_from_v(v1_mag+v2_mag, t_master[:min_len_r], 9)
R_reconstructed = R_reconstructed[::step]

# --- 3. Decimate (Sub-sampling) ---
t_final  = t_master[::step]
v1_final = v1_mag[::step]
v2_final = v2_mag[::step]
R_final  = R_actual[::step]
v1_z_final = v1_z[::step]
v2_z_final = v2_z[::step]
v1_x_final = v1_x[::step]
v2_x_final = v2_x[::step]

# --- 4. Export to Unified CSV ---
export_df = pd.DataFrame({
    'Time_ps': t_final,
    'V1_mag': v1_final,
    'V2_mag': v2_final,
    'V1_z': np.abs(v1_z_final),
    'V2_z': np.abs(v2_z_final),
    'V1_x': np.abs(v1_x_final),
    'V2_x': np.abs(v2_x_final),
    'R_distance': R_final
})

# Save to CSV
export_df.to_csv(output_file, index=False)
print(f"Export complete: {output_file}")
print(f"Rows reduced from {len(t_master)} to {len(t_final)}")

# --- 5. Verification Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax1.plot(t_final, v1_final, label='Atom 1 Speed')
ax1.plot(t_final, v2_final, label='Atom 2 Speed')
ax1.set_ylabel('v / Å/ps')
ax1.set_title('Decimated Velocity Data (Every 100th Point)')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.plot(t_final, R_final, color='black', label='Ground Truth Distance')
ax2.plot(t_final, R_reconstructed, color='red', label='Reconstructed Distance')
ax2.set_ylabel('R / Å')
ax2.set_xlabel('t / ps')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
# 18 A case: Load, Plotting & Exporting
# -------------------------------------------------------------------------------------


# --- Configuration ---
if office:
    base_path = r"T:\NextCloud_PaulGuggenbichler\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT\18A"
else:
    base_path = r"C:\Users\paulg\Dokumente\Studium\Masterarbeit\Drag_Calculation\Data_DFT\18A"
output_file = os.path.join(base_path, "18A_All_Data.csv")
step = 100  # Take every 100th point

def marti_loader(path):
    """Robust loader for large Marti files with scientific notation."""
    return np.genfromtxt(path, comments='#', invalid_raise=False)

# --- 1. Load Original Data ---
print("Loading large Marti files...")
v1_raw = marti_loader(os.path.join(base_path, "vimp.I_1"))
v2_raw = marti_loader(os.path.join(base_path, "vimp.I_2"))
r1_raw = marti_loader(os.path.join(base_path, "rimp.I_1"))
r2_raw = marti_loader(os.path.join(base_path, "rimp.I_2"))

# --- 2. Calculate Magnitudes and Distance ---
# Calculate R (Distance) from positions
min_len_r = min(len(r1_raw), len(r2_raw), len(v2_raw), len(v1_raw))  # Ensure we only use the overlapping portion of all datasets
t_master = r1_raw[:min_len_r, 0]
R_actual = np.linalg.norm(r1_raw[:min_len_r, 1:4] - r2_raw[:min_len_r, 1:4], axis=1)

# Calculate V1 Magnitude
v1_mag = np.linalg.norm(v1_raw[:min_len_r, 1:4], axis=1)
v1_z = v1_raw[:min_len_r, 3]
v1_x = v1_raw[:min_len_r, 1]
v2_mag = np.linalg.norm(v2_raw[:min_len_r, 1:4], axis=1)
v2_z = v2_raw[:min_len_r, 3]
v2_x = v2_raw[:min_len_r, 1]

R_reconstructed = reconstruct_R_from_v(v1_mag+v2_mag, t_master[:min_len_r], 18)
R_reconstructed = R_reconstructed[::step]

# --- 3. Decimate (Sub-sampling) ---
t_final  = t_master[::step]
v1_final = v1_mag[::step]
v2_final = v2_mag[::step]
v1_z_final = v1_z[::step]
v2_z_final = v2_z[::step]
v1_x_final = v1_x[::step]
v2_x_final = v2_x[::step]
R_final  = R_actual[::step]

# --- 4. Export to Unified CSV ---
export_df = pd.DataFrame({
    'Time_ps': t_final,
    'V1_mag': v1_final,
    'V2_mag': v2_final,
    'V1_z': np.abs(v1_z_final),
    'V2_z': np.abs(v2_z_final),
    'V1_x': np.abs(v1_x_final),
    'V2_x': np.abs(v2_x_final),
    'R_distance': R_final
})

# Save to CSV
export_df.to_csv(output_file, index=False)
print(f"Export complete: {output_file}")
print(f"Rows reduced from {len(t_master)} to {len(t_final)}")

# --- 5. Verification Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

ax1.plot(t_final, v1_final, label='Atom 1 Speed (Full)')
ax1.plot(t_final, v2_final, label='Atom 2 Speed (Shortened)')
ax1.set_ylabel('v / Å/ps')
ax1.set_title('Decimated Velocity Data (Every 100th Point)')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.plot(t_final, R_final, color='black', label='Ground Truth Distance')
ax2.plot(t_final, R_reconstructed, color='red', label='Reconstructed Distance')
ax2.set_ylabel('R / Å')
ax2.set_xlabel('t / ps')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
