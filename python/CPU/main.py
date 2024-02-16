import FDFD.FDFD2D as fdfd
import numpy as np
## Simulation domain parameters
Lx = 60e-3; # Width of simulation domain in meter
Ly = 60e-3; # Same sas Ly but height
ds = 0.2e-3; # Distance between adjecent nodes in the grid
f = float(70e9); # Frequency to simulated

## PML parameters
x_bnd = True; # Put PML on the x boundaries
y_bnd = True; # Put PML on the y boundaries
pml_width = 10; # Width of the pml in number of nodes

## Create simulation instance
my_sim = fdfd.FDFD_simDomain('my_sim', Lx, Ly, ds);
my_sim.snap_grid(); # Snap simualtion to grid

## Create perfectly matched layer
# x_bnd and y_bnd = True -> pml on x and y boundaries
# pml_width -> Width of the pml in number of nodes
my_sim.PML = fdfd.FDFD_upml(pml_width, x_bnd, y_bnd); # Create PML


## Materials
# Initialize materials to vakuum and pml
materials = fdfd.FDFD_materials(my_sim); # Create material instance
materials.init_materials(my_sim); # Initialize materials to vacuum and pml

a = 8e-3;
x0 = np.array([0.0,0.0]);
eps_r = np.array([2.25]);
mu_r = np.array([1.0]);

primitive_geometry = fdfd.FDFD_primitive_geometry(my_sim);
primitive_geometry.disc(materials,eps_r,mu_r,x0,a);

## Differential operators
diff_operators = fdfd.FDFD_operators(my_sim); # Initialize differential operators
diff_operators.gen_diff_operators(my_sim, f); # Differential operators at 10 GHz
## Each frequency gives different differential operators !!!


## Assemble linear system. A[0] is for TM mode and A[1] is for TE mode
A = diff_operators.assemble_systems(materials);


# We have assambled the system Ae = 0. We need a source b.


source = fdfd.FDFD_source(my_sim);
source.plane_wave(my_sim, 1, -np.pi/2, f); # Source is a plane wave angled 60 degrees
##source.plot_source(my_sim);



tfsf = fdfd.FDFD_tfsf(my_sim, my_sim.PML.pml_width + 6)
tfsf.masking_fnc(my_sim);
b = tfsf.assemble_source(A, source.f_src);

## Now we have the systm Ax = b;
# All that is left is to solve the system for x using a direct or iterative solver

solver = fdfd.FDFD_solver(A,b);
electric_fields_TM = solver.LU_solve_TM();

efield = fdfd.FDFD_fields(my_sim, electric_fields_TM);
efield.plot_fields(my_sim);





