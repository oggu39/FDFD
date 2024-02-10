import cupy as cp
import cupyx.scipy as cps
import cupyx.scipy.sparse.linalg as cpsl


Nx = 10;
dx = 1/(Nx-1);
dx2 = dx**2;
k = 6;

# Sparse matrix
A = cps.sparse.csr_matrix((Nx,Nx));

# Source b
b = dx2*(-k)*cp.ones((Nx,));
b[0] = 0; b[-1] = 0;

# Fill sparse matrix
d1 = -2*cp.ones((Nx,));
d1[0] = 1; d1[-1] = 1;
d2 = cp.ones((Nx-1,));
d2[0] = 0;

d3 = cp.ones((Nx-1,));
d3[-1] = 0;

A.setdiag(d1,k=0);
A.setdiag(d2,k=1);
A.setdiag(d3,k=-1);

x = cpsl.spsolve(A,b);



