
import FDFD.FDFD2D as fdfd;
Lx = 600e-3;
Ly = 6.2e-3;
ds = 0.1e-3;
pml_width = 10;
x_bnd = True;
y_bnd = True;

my_sim = fdfd.FDFD_simDomain('my_sim', Lx, Ly, ds)
my_sim.snap_grid();

diff_operators = fdfd.FDFD_operators(my_sim);

freq = 10e9;
diff_operators.FDFD_gen_operators(my_sim, freq);


