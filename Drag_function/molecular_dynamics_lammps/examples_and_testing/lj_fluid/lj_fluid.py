from lammps import lammps
import lammps_logfile
from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np
import time

# Initialize LAMMPS instance
log_filename = "log.lj"
lmp_1 = lammps(cmdargs=["-log", log_filename, "-screen", "none"])
# Use mpi4py to get the rank, ensuring MPI is initialized correctly
if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"LAMMPS Version Date Code: {lmp_1.version()}")


lj_script = """
# 3D Lennard-Jones melt
# 1) Initialization
units lj
dimension 3
atom_style atomic
boundary p p p
# 2) System definition
region simbox block -20 20 -20 20 -20 20
create_box 2 simbox
create_atoms 1 random 1500 34134 simbox overlap 0.3
create_atoms 2 random 100 12756 simbox overlap 0.3
# 3) Settings
mass 1 1.0
mass 2 5.0
pair_style lj/cut 4.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 0.5 3.0
# 4) Monitoring
thermo 10
thermo_style custom step ke pe
# 5) Run

minimize 1.0e-6 1.0e-6 1000 10000

# PART B - MOLECULAR DYNAMICS
# 4) Monitoring
thermo 50
thermo_style custom step temp etotal pe ke press

# 5) Run
fix mynve all nve
fix mylgv all langevin 1.0 1.0 0.1 10917
timestep 0.005
dump mydmp all atom 100 dump.lj
run 15000
"""

# Execute the commands line-by-line and time it
start_time = time.time()
lmp_1.commands_string(lj_script)
end_time = time.time()

elapsed_time = end_time - start_time
print(f"\n{'='*60}")
print(f"LAMMPS Calculation Time: {elapsed_time:.4f} seconds ({elapsed_time/60:.2f} minutes)")
print(f"{'='*60}\n")




#--------------------------------------------------------------------------------------------------------------
#                                   Post processing
#--------------------------------------------------------------------------------------------------------------

# Reads the log file. If you have multiple 'run' or 'minimize' commands,
# it stores them as separate blocks.
log = lammps_logfile.File("log.lj")

# Get data from the first run/minimize block
step = log.get("Step", run_num=0)
pe = log.get("PotEng", run_num=0)
ke = log.get("KinEng", run_num=0)

step_md = log.get("Step", run_num=1)
pe_md = log.get("PotEng", run_num=1)
ke_md = log.get("KinEng", run_num=1)


# 4. Plot the Data in 2x2 layout
plt.figure(figsize=(12, 10))

# Top-left: Minimization Kinetic Energy
plt.subplot(2, 2, 3)
plt.plot(step, ke, color='red', linewidth=1.5)
plt.title('Minimization: Kinetic Energy (Always 0)')
plt.xlabel('Minimization Step')
plt.ylabel('KE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

# Top-right: Minimization Potential Energy
plt.subplot(2, 2, 1)
plt.ylim(np.min(pe) * 1.1, 1)
plt.plot(step, pe, color='blue', linewidth=1.5)
plt.title('Minimization: Potential Energy Descent')
plt.xlabel('Minimization Step')
plt.ylabel('PE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

# Bottom-left: MD Potential Energy
plt.subplot(2, 2, 2)
plt.plot(step_md, pe_md, color='blue', linewidth=1.5)
plt.title('MD: Potential Energy')
plt.xlabel('MD Step')
plt.ylabel('PE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

# Bottom-right: MD Kinetic Energy
plt.subplot(2, 2, 4)
plt.plot(step_md, ke_md, color='red', linewidth=1.5)
plt.title('MD: Kinetic Energy')
plt.xlabel('MD Step')
plt.ylabel('KE per atom')
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()

#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
#                                       Optimizing Simulation
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------


# Initialize LAMMPS instance (Want to overwrite as this is example)
log_filename = "log.lj_improved"
lmp_2 = lammps(cmdargs=["-log", log_filename, "-screen", "none"])
# Use mpi4py to get the rank, ensuring MPI is initialized correctly
if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"LAMMPS Version Date Code: {lmp_2.version()}")


lj_script = """
# 3D Lennard-Jones melt
# 1) Initialization
units lj
dimension 3
atom_style atomic
boundary p p p
# 2) System definition --> Want two separate regions
region simbox block -20 20 -20 20 -20 20
create_box 2 simbox
# for creating atoms
region cyl_in cylinder z 0 0 10 INF INF side in
region cyl_out cylinder z 0 0 10 INF INF side out
create_atoms 1 random 1000 34134 cyl_out
create_atoms 2 random 150 12756 cyl_in
# 3) Settings
mass 1 1.0
mass 2 5.0
pair_style lj/cut 4.0
pair_coeff 1 1 1.0 1.0
pair_coeff 2 2 0.5 3.0
# 4) Monitoring
thermo 10
thermo_style custom step ke pe
# 5) Run

minimize 1.0e-6 1.0e-6 1000 10000

# PART B - MOLECULAR DYNAMICS
# 4) Monitoring
thermo 50
thermo_style custom step temp etotal pe ke press

# 5) Run
fix mynve all nve
fix mylgv all langevin 1.0 1.0 0.1 10917
timestep 0.005
dump mydmp all atom 100 dump.lj_improved
run 15000

# 6) Save system
write_data improved.min.data

"""

# Execute the commands line-by-line and time it
start_time = time.time()
lmp_2.commands_string(lj_script)
end_time = time.time()

elapsed_time = end_time - start_time
print(f"\n{'='*60}")
print(f"LAMMPS Calculation Time: {elapsed_time:.4f} seconds ({elapsed_time/60:.2f} minutes)")
print(f"{'='*60}\n")

# lmp_1/2.finalize() explicitly calls MPI_Finalize() before shutting down LAMMPS
lmp_1.finalize()
lmp_2.finalize()