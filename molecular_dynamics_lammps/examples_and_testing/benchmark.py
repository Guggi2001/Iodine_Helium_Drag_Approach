from lammps import lammps
from mpi4py import MPI


# Initialize LAMMPS instance
lmp = lammps()
# Use mpi4py to get the rank, ensuring MPI is initialized correctly
if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"LAMMPS Version Date Code: {lmp.version()}")

# Define the LJ melt commands as a multi-line string
lj_script = """
# 3D Lennard-Jones melt
units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 1.44 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        0.3 bin
neigh_modify    delay 0 every 20 check no
fix             1 all nve
thermo          50
run             250
"""

# Execute the commands line-by-line
lmp.commands_string(lj_script)

# lmp.finalize() explicitly calls MPI_Finalize() before shutting down LAMMPS
lmp.finalize()

# For parallel runs, you need to type "mpiexec -np 4 python benchmark.py" in the terminal