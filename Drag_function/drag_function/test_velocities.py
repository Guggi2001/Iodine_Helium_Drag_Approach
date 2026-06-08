from config_utils_local import config as C
from drag_function import io
import numpy as np
import matplotlib.pyplot as plt

dict9 = io.load_data(C.PATH9A)
dict18 = io.load_data(C.PATH18A)

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