from lammps import lammps
import lammps_logfile
from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np
import time

# Initialize LAMMPS instance
log_filename = "log.cnt"
lmp = lammps(cmdargs=["-log", log_filename, "-screen", "none"])
# Use mpi4py to get the rank, ensuring MPI is initialized correctly
if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"LAMMPS Version Date Code: {lmp.version()}")

lj_script = """
# 3D Lennard-Jones melt
# 1) Initialization
units real
atom_style molecular
boundary f f f
# 2) System definition
pair_style lj/cut 14.0
bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style harmonic
special_bonds lj 0.0 0.0 0.5
# 3) Settings
read_data unbreakable.data
include unbreakable.inc # Instead of copying content here showcasing: better to repeat

# In this tutorial, a deformation will be applied to the CNT by displacing the atoms located at its edges.
# To achieve this we will first isolate the atoms at the two edges and place them into groups named rtop and rbot
group carbon_atoms type 1 # includes all atoms of type 1 in group carbon_atoms
variable xmax equal bound(carbon_atoms,xmax)-0.5
variable xmin equal bound(carbon_atoms,xmin)+0.5
region rtop block ${xmax} INF INF INF INF INF   # x>x_max
region rbot block INF ${xmin} INF INF INF INF   # x<x_min
region rmid block ${xmin} ${xmax} INF INF INF INF # x_min < x < x_max

# Define 3 groups pf atoms for each 3 region
group cnt_top region rtop
group cnt_bot region rbot
group cnt_mid region rmid
set group cnt_top mol 1
set group cnt_bot mol 2
set group cnt_mid mol 3

# Delete some random atoms - to avoid to close to edges new region defined
variable xmax_del equal ${xmax}-2
variable xmin_del equal ${xmin}+2
region rdel block ${xmin_del} ${xmax_del} INF INF INF INF
group rdel region rdel
delete_atoms random fraction 0.02 no rdel NULL 2793 bond yes

# 4) Monitoring
variable n_carbon equal count(carbon_atoms)
variable n_top equal count(cnt_top)
variable n_bot equal count(cnt_bot)
variable n_mid equal count(cnt_mid)
write_dump all image carbon_nanotube.ppm element type size 1000 400 zoom 6 shiny &
    0.3 fsaa yes bond atom 0.8 view 0 90 box no 0.0 axes no 0.0 0.0 modify backcolor &
    white adiam 1 0.85 bdiam 1 1.0

# MOLECULAR DYNAMICS
# Re-setting the atom IDs is necessary before using the velocity command when atoms were deleted
reset_atoms id sort yes
velocity cnt_mid create 300 48455 mom yes rot yes

#Specify thermalization and dynamics system
fix mynve1 cnt_top nve
fix mynve2 cnt_bot nve
fix mynvt cnt_mid nvt temp 300 300 100
    
run 0 post no
"""

# Execute the commands line-by-line and time it
start_time = time.time()
lmp.commands_string(lj_script)
end_time = time.time()

elapsed_time = end_time - start_time
print(f"\n{'='*60}")
print(f"LAMMPS Calculation Time: {elapsed_time:.4f} seconds ({elapsed_time/60:.2f} minutes)")
print(f"{'='*60}\n")

# Extract variables (style 0 = global, type 0 = scalar)
n_carbon = int(lmp.extract_variable("n_carbon", 0, 0))
n_top = int(lmp.extract_variable("n_top", 0, 0))
n_bot = int(lmp.extract_variable("n_bot", 0, 0))
n_mid = int(lmp.extract_variable("n_mid", 0, 0))

print(f"carbon_atoms: {n_carbon}")
print(f"cnt_top: {n_top}")
print(f"cnt_bot: {n_bot}")
print(f"cnt_mid: {n_mid}")

# lmp.finalize() explicitly calls MPI_Finalize() before shutting down LAMMPS
lmp.finalize()

