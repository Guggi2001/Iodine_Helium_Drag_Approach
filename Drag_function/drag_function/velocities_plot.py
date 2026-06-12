from config_utils_local import config as C
from drag_function import io
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import pandas as pd
home = False
dict9 = io.load_data(C.PATH9A)
dict18 = io.load_data(C.PATH18A)


t18 = dict18["t"]
v18 = dict18["v1"]
R18 = dict18["R"]
mask = (t18 >= 4.54) & (t18 <= 14.76800)
t_18_w = t18[mask]
v_18_w = v18[mask]
R_18_w = R18[mask]
wl = 9401
polyorder = 1
v_18_SG = savgol_filter(v_18_w, window_length=wl, polyorder=polyorder, deriv=0, mode="interp")

# Create DataFrame with time, cleaned_SG, and IMF cleaned data
export_data = pd.DataFrame({
    'time': t_18_w,
    'cleaned_SG': v_18_SG,
})

# Export to CSV
export_data.to_csv('cleaned_data_18_long.csv', index=False)
print(f"Data exported to cleaned_data_long.csv")

plt.figure(figsize=(12,8))
plt.plot(t18, v18)
plt.plot(t_18_w,v_18_SG )
plt.show()




a = 3







v_9x = dict9['v2_x']
v_9y = dict9['v2_y']
v_9z = dict9['v2_z']
v_9 = dict9['v2']
t_9 = dict9['t']
a = 3
plt.figure(figsize=(12,8))
plt.plot(t_9, v_9x, label = r'$v_x$')
plt.plot(t_9, v_9y, label = r'$v_y$')
plt.plot(t_9, v_9z, label = r'$v_z$')
plt.plot(t_9, v_9, label = r'$v_magn$')
plt.title('9A Top Atom')
plt.legend()
plt.show()


v_9x_1 = dict9['v1_x']
v_9y_1 = dict9['v1_y']
v_9z_1 = dict9['v1_z']
v_9_1 = dict9['v1']

plt.figure(figsize=(12,8))
plt.plot(t_9, v_9x_1, label = r'$v_x$')
plt.plot(t_9, v_9y_1, label = r'$v_y$')
plt.plot(t_9, v_9z_1, label = r'$v_z$')
plt.plot(t_9, -v_9_1, label = r'$v_magn$')
plt.title('9A Bottom Atom')
plt.legend()
plt.show()


# --- 18A case: same visualizations as for 9A ---
v_18x = dict18['v2_x']
v_18y = dict18['v2_y']
v_18z = dict18['v2_z']
v_18 = dict18['v2']
t_18 = dict18['t']

plt.figure(figsize=(12,8))
plt.plot(t_18, v_18x, label = r'$v_x$')
plt.plot(t_18, v_18y, label = r'$v_y$')
plt.plot(t_18, v_18z, label = r'$v_z$')
plt.plot(t_18, v_18, label = r'$v_magn$')
plt.title('18A Top Atom')
plt.legend()
plt.show()


v_18x_1 = dict18['v1_x']
v_18y_1 = dict18['v1_y']
v_18z_1 = dict18['v1_z']
v_18_1 = dict18['v1']

plt.figure(figsize=(12,8))
plt.plot(t_18, v_18x_1, label = r'$v_x$')
plt.plot(t_18, v_18y_1, label = r'$v_y$')
plt.plot(t_18, v_18z_1, label = r'$v_z$')
plt.plot(t_18, -v_18_1, label = r'$v_magn$')
plt.title('18A Bottom Atom')
plt.legend()
plt.show()


R9 = dict9['R']
r1_9 = np.vstack((dict9['x1'], dict9['y1'], dict9['z1'])).T
r2_9 = np.vstack((dict9['x2'], dict9['y2'], dict9['z2'])).T
R9_test = np.linalg.norm(r1_9 - r2_9, axis=1)
plt.figure(figsize=(12,8))
plt.plot(t_9, R9, label = r'$R$')
plt.plot(t_9, R9_test, ls = '--', label = 'R from positions')
plt.title('9A Distance')
plt.legend()
plt.show()


R18 = dict18['R']
r1_18 = np.vstack((dict18['x1'], dict18['y1'], dict18['z1'])).T
r2_18 = np.vstack((dict18['x2'], dict18['y2'], dict18['z2'])).T
R18_test = np.linalg.norm(r1_18 - r2_18, axis=1)
plt.figure(figsize=(12,8))
plt.plot(t_18, R18, label = r'$R$')
plt.plot(t_18, R18_test, ls = '--', label = 'R from positions')
plt.title('18A Distance')
plt.legend()
plt.show()





from i2_helium_md.physics.potentials import droplet_potential

R = 56
r = np.linspace(0,100, 1000)
plt.figure(figsize=(12,8))
plt.plot(r-R, droplet_potential(r-R, steepness = 5, binding_energy = 0.3 ),label = r'$Steepness 5$')
plt.plot(r-R, droplet_potential(r-R, steepness = 9, binding_energy = 0.3 ),label = r'$Steepness 9$')
plt.plot(r-R, droplet_potential(r-R, steepness = 14.2, binding_energy = 0.3 ),label = r'$Steepness 14.2$')

plt.xlabel(r'$r-Rdroplet$')
plt.ylabel(r'$V$ / eV')
plt.title('Droplet Potential')
plt.legend()
plt.show()

if home:
    import numpy as np
    base_path_clean = r"C:\Users\paulg\Dokumente\GitHub\Iodine_Helium_Drag_Approach\i2_helium_md\data\reference\drag\9A\velocity_smoothed"
    cleaned = base_path_clean + r"\mean_velocity.csv"
    structured = np.genfromtxt(cleaned, delimiter=",", names=True, dtype=float)
    time_ps = np.atleast_1d(np.asarray(structured["time_ps"], dtype=float))
    speed_Aps = np.atleast_1d(np.asarray(structured["mean_velocity_Aps"], dtype=float))
    cleaned_SG = base_path_clean + r"\cleaned_data_long.csv"
    structured_SG = np.genfromtxt(cleaned_SG, delimiter=",", names=True, dtype=float)
    time_SG = np.atleast_1d(np.asarray(structured_SG["time"], dtype=float))
    velocity_SG = np.atleast_1d(np.asarray(structured_SG["cleaned_SG"], dtype=float))


    plt.figure(figsize=(12,8))
    plt.plot(time_ps, speed_Aps, label = r'$mean velocity MD$')
    plt.plot(time_SG, velocity_SG, label = r'$mean velocity SG$')
    plt.plot(t_9, v_9, label = r'$v_magn$')
    plt.title('9A Top Atom')
    plt.legend()
    plt.show()
    a = 3
