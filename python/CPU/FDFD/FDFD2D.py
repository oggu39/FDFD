import numpy as np
import scipy

class FDFD_constants:
        c0 = float(299792458); # m/s
        mu0 = 4*np.pi*1e-7; # H/m
        eps0 = 1/(c0*c0*mu0); # F/m
        Z0 = np.sqrt(mu0/eps0); # Free space impedence[ohm]


class FDFD_simDomain:
    def __init__(self, name, Lx, Ly, ds):
        self.name = name;
        self.Lx = Lx;
        self.Ly = Ly;
        self.ds = ds;
        self.Nx = False; # These values will be calculated once snap_grid() is run
        self.Ny = False;
        self.PML = False; # Should be assigned once pml is created
        
    def snap_grid(self):
        self.Nx = (np.ceil(self.Lx/self.ds) + 1).astype(int);
        self.Ny = (np.ceil(self.Ly/self.ds) + 1).astype(int);
        
        self.ds = self.Lx/(self.Nx-1); # Snap across x-boundary
        self.Ly = self.ds*(self.Ny-1); # Calculate new boundary of y
        



class FDFD_upml():
    def __init__(self, pml_width,x_bnd, y_bnd):
        self.pml_width = pml_width;
        self.x_bnd = y_bnd;
        self.y_bnd = y_bnd;   
        self.sgm_max = float(1);
        self.p = float(3);
        self.s_max=  float(2.5);
        
        
        
    def upml_sgm(self): # Ficticious conductivity in the pml
        w = np.linspace(0,1,self.pml_width); # Width along the pml
        return self.sgm_max*(np.sin(0.5*np.pi*w)**2);
    
    
    
    def upml_s0(self):
        if(self.pml_width > 1):
            w = np.linspace(0,1,self.pml_width); # Width along the pml
            return 1+self.s_max*(w**self.p); # w is length along the pml 0 - 1
        else:
            raise Exception("Pml must have a width more than one node");
            
            
            
    def upml_s(self): 
        if(self.pml_width > 1):
            sgm = self.upml_sgm();
            s0 = self.upml_s0();
            s = s0*(1-1j*FDFD_constants.Z0*sgm); # S-par of pml along the width of pml
        else:
            raise Exception("Pml must have a width more than one node")
            
        return s;
        
        
class FDFD_operators:
    
    def __init__(self, simDomain):
        # Initialise differential operator to empty sparse. Total number of nodes(unknowns) = Nx*Ny -> matrix size of Nx*Ny x Nx*Ny
        self.Dxe = scipy.sparse.lil_matrix((simDomain.Nx*simDomain.Ny, simDomain.Nx*simDomain.Ny)); 
        self.Dye = scipy.sparse.lil_matrix((simDomain.Nx*simDomain.Ny, simDomain.Nx*simDomain.Ny)); 
        self.Dxh = scipy.sparse.lil_matrix((simDomain.Nx*simDomain.Ny, simDomain.Nx*simDomain.Ny)); 
        self.Dyh = scipy.sparse.lil_matrix((simDomain.Nx*simDomain.Ny, simDomain.Nx*simDomain.Ny)); 
        
    def FDFD_gen_operators(self, simDomain,freq):
        Nx = simDomain.Nx; # Number of nodes in the x - direction should be integer
        Ny = simDomain.Ny; # Number of nodes in the y - direction should be integer
        ds = simDomain.ds; # Distance between nodes. dx = dy = ds(cubic grid is assumed)
        k0 = 2*np.pi*freq/FDFD_constants.c0; # Wave number
        
        
        self.Dxe.setdiag(-1,k=0);
        self.Dxe.setdiag(1,k=1);
        
        
        # Set dirichlet condition at the right boundary(x = Nx*ds)
        for q in range(Nx-1,Nx,Nx*Ny-1):
            self.Dxe[q,q+1] = 0;
        
        
        self.Dxe = self.Dxe.tocsr()/(k0*ds); # Compressed Sparse Row (CSR) format

        
        # Dirichlet condition for the y differential is implicitly fullfilled
        self.Dye.setdiag(-1,k=0);
        self.Dye.setdiag(1,k=Nx); 
        self.Dye = self.Dye.tocsr()/(k0*ds); # Compressed Sparse Row(CSR) format
        
        # Differential matricies for H is transpose of electric field 
        self.Dxh = -self.Dxe.T; # Should be CSC format(Compressed Sparse Column)
        self.Dyh = -self.Dye.T; # Should be CSC format
        
class FDFD_materials: 
    def __init__(self, simDomain):
        self.Nx = simDomain.Nx;
        self.Ny = simDomain.Ny;
        self.eps_xx = False;
        self.eps_yy = False;
        self.eps_zz = False;
        self.mu_xx = False;
        self.mu_yy = False;
        self.mu_zz = False;
        
        

        
        
        
        
        
        
        
        