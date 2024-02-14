import FDFD.FDFD2D as fdfd
import numpy as np
## Simulation domain parameters
Lx = 50e-3;
Ly = 50e-3;
ds = 0.15e-3;

## PML parameters
x_bnd = True; 
y_bnd = True;
pml_width = 15;

## Create simulation instance
my_sim = fdfd.FDFD_simDomain('my_sim', Lx, Ly, ds);
my_sim.snap_grid(); # Snap simualtion to grid

## Create perfectly matched layer
# x_bnd and y_bnd = True -> pml on x and y boundaries
# pml_width -> Width of the pml in number of nodes
my_sim.PML = fdfd.FDFD_upml(pml_width, x_bnd, y_bnd);


## Materials
# Initialize materials to vakuum and pml
materials = fdfd.FDFD_materials(my_sim); # Create material instance
materials.init_materials(my_sim); # Initialize materials to vacuum and pml

## Differential operators
diff_operators = fdfd.FDFD_operators(my_sim); # Initialize differential operators
diff_operators.gen_diff_operators(my_sim, float(10e9)); # Differential operators at 10 GHz
## Each frequency gives different differential operators !!!

## Assemble linear system. A[0] is for TM mode and A[1] is for TE mode
A = diff_operators.assemble_systems(materials);

# We have assambled the system Ae = 0. We need a source b.
# After that all that is left is to solve Ae = b and we have our 
# Electric fields in the simulation domina e = inv(A)*b;

source = fdfd.FDFD_source(my_sim);
source.plane_wave(my_sim, 1, np.pi/3, 10e9);
source.plot_source(my_sim);

tfsf = fdfd.FDFD_tfsf(my_sim, my_sim.PML.pml_width + 6)
tfsf.masking_fnc(my_sim);

b = tfsf.assemble_source(A, source.f_src);

"""
materials.init_materials(my_sim);

materials.visualize_materials();

"""

