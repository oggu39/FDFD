# Package needs numpy as scipy if run on cpu
import numpy as np
import scipy
import matplotlib.pyplot as plt


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
        
        
        
    def upml_sgm(self,w): # Ficticious conductivity in the pml at the width w ~(0,1)
        if(w >= 0 or w <= 1):
            return self.sgm_max*(np.sin(0.5*np.pi*w)**2);
        else:
            raise Exception("Outside of the width of the pml sgm")
    
    
    
    def upml_s0(self, w):
        if(w >= 0 and w <= 1):
            return 1+self.s_max*(w**self.p); # w is length along the pml 0 - 1
        else:
            print(w)
            raise Exception("Outside of the width of the pml s0");
            
            
            
    def upml_s(self,w): 
        if(w >= 0 or w <= 1):
            sgm = self.upml_sgm(w);
            s0 = self.upml_s0(w);
            s = s0*(1-1j*FDFD_constants.Z0*sgm); # S-par of pml along the width w of pml
            return s;
        else:
            raise Exception("Pml must have a width more than one node s")
            
        
        
        

class FDFD_materials: 
    def __init__(self, simDomain):
        self.Nx = simDomain.Nx;
        self.Ny = simDomain.Ny;
        self.eps_xx = scipy.sparse.lil_matrix((self.Nx*self.Ny,self.Nx*self.Ny)).astype(complex);
        self.eps_yy = scipy.sparse.lil_matrix((self.Nx*self.Ny,self.Nx*self.Ny)).astype(complex);
        self.eps_zz = scipy.sparse.lil_matrix((self.Nx*self.Ny,self.Nx*self.Ny)).astype(complex);
        self.mu_xx = scipy.sparse.lil_matrix((self.Nx*self.Ny,self.Nx*self.Ny)).astype(complex);
        self.mu_yy = scipy.sparse.lil_matrix((self.Nx*self.Ny,self.Nx*self.Ny)).astype(complex);
        self.mu_zz = scipy.sparse.lil_matrix((self.Nx*self.Ny,self.Nx*self.Ny)).astype(complex);
        
    # Initialize entire grid to vacuumm
    def init_materials(self, simDomain):
        
        Nx = self.Nx;
        Ny = self.Ny;
        pml_width = simDomain.PML.pml_width;
        eps_xx = np.ones((Nx*Ny,)).astype(complex);
        eps_yy = np.ones((Nx*Ny,)).astype(complex);
        eps_zz = np.ones((Nx*Ny,)).astype(complex);

        
        
        mu_xx = np.ones((Nx*Ny,)).astype(complex);
        mu_yy = np.ones((Nx*Ny,)).astype(complex);
        mu_zz = np.ones((Nx*Ny,)).astype(complex);
        
        sx = np.ones((Nx*Ny,)).astype(complex);
        sy = np.ones((Nx*Ny,)).astype(complex);
        # sz = np.ones((Nx*Ny,)); Isn't needed in 2D
        
        # Incorporate pmls on all sides
        # left side
        for i in range(0,pml_width):
            for j in range(0,Ny):
                w = (pml_width-i)/pml_width;
                q = i + j*Nx;
                sx[q] = simDomain.PML.upml_s(w);

        # right side
        for i in range(Nx-pml_width,Nx,1):
            for j in range(0,Ny):
                w = (i-(Nx-pml_width)+1)/pml_width;
                q = i + j*Nx;
                sx[q] = simDomain.PML.upml_s(w);
                
        
        # bottom side
        for i in range(0,Nx):
            for j in range(0,pml_width):
                w = (pml_width-j)/pml_width;
                q = i+j*Nx;
                sy[q] = simDomain.PML.upml_s(w);
        
        # top
        for i in range(0,Nx):
            for j in range(Ny-pml_width,Ny,1):
                w = (j-(Ny-pml_width)+1)/pml_width;
                q = i + j*Nx;
                sy[q] = simDomain.PML.upml_s(w);
                
        
        # Initialize materials in the entire grid to pml and vakum
        self.eps_xx.setdiag(eps_xx*sy/sx,k=0);
        self.eps_yy.setdiag(eps_yy*sx/sy,k=0);
        self.eps_zz.setdiag(eps_zz*sx*sy,k=0);
        
        self.mu_xx.setdiag(mu_xx*sy/sx,k=0);
        self.mu_yy.setdiag(mu_yy*sx/sy,k=0);
        self.mu_zz.setdiag(mu_zz*sx*sy,k=0);
   
                
    def get_materials(self):
        Nx = self.Nx;
        Ny = self.Ny;
        eps_xx = np.ones((Nx,Ny)).astype(complex);
        eps_yy = np.ones((Nx,Ny)).astype(complex);
        eps_zz = np.ones((Nx,Ny)).astype(complex);
        mu_xx = np.ones((Nx,Ny)).astype(complex);
        mu_yy = np.ones((Nx,Ny)).astype(complex);
        mu_zz = np.ones((Nx,Ny)).astype(complex);
        
        for i in range(0,Nx):
            for j in range(0,Ny):
                q = i+j*Nx;
                eps_xx[i,j] = self.eps_xx[q,q];
                eps_yy[i,j] = self.eps_yy[q,q];
                eps_zz[i,j] = self.eps_zz[q,q];
                
                mu_xx[i,j] = self.mu_xx[q,q];
                mu_yy[i,j] = self.mu_yy[q,q];
                mu_zz[i,j] = self.mu_zz[q,q];
                
        return [eps_xx, eps_yy, eps_zz, mu_xx, mu_yy, mu_zz];
    
    def visualize_materials(self):
        
        materials = self.get_materials();
        eps_z = np.abs(materials[2]);
        plt.pcolor(eps_z,vmin=0,vmax=2000);
        

        
        
        
class FDFD_operators:
    
    def __init__(self, simDomain):
        # Initialise differential operator to empty sparse. Total number of nodes(unknowns) = Nx*Ny -> matrix size of Nx*Ny x Nx*Ny
        self.Dxe = scipy.sparse.lil_matrix((simDomain.Nx*simDomain.Ny, simDomain.Nx*simDomain.Ny)).astype(complex); 
        self.Dye = scipy.sparse.lil_matrix((simDomain.Nx*simDomain.Ny, simDomain.Nx*simDomain.Ny)).astype(complex); 
        self.Dxh = scipy.sparse.lil_matrix((simDomain.Nx*simDomain.Ny, simDomain.Nx*simDomain.Ny)).astype(complex); 
        self.Dyh = scipy.sparse.lil_matrix((simDomain.Nx*simDomain.Ny, simDomain.Nx*simDomain.Ny)).astype(complex); 
        
    def gen_diff_operators(self, simDomain,freq):
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
        self.Dxh = -self.Dxe.T.tocsr(); # Should be CSR format(Compressed Sparse Row)
        self.Dyh = -self.Dye.T.tocsr(); # Should be CSR format
        
    def assemble_systems(self,materials):
            Ae = self.Dxh.multiply(materials.mu_yy.power(-1).multiply(self.Dxe)) + \
                 self.Dyh.multiply(materials.mu_xx.power(-1).multiply(self.Dye)) + \
                 materials.eps_zz;
            
            Ah = self.Dxe.multiply(materials.eps_yy.power(-1).multiply(self.Dxh)) + \
                 self.Dye.multiply(materials.eps_xx.power(-1).multiply(self.Dyh)) + \
                 materials.mu_zz;
            return [Ae,Ah];
        

class FDFD_source:
    def __init__(self,simDomain):
        self.f_src = np.ones((simDomain.Nx*simDomain.Ny,)).astype(complex);
        
    def plane_wave(self,simDomain,amp,theta,freq):
        k = 2*np.pi*freq/FDFD_constants.c0; # Magnitude of wave vector
        
        for i in range(0,simDomain.Nx):
            for j in range(0,simDomain.Ny):
                x = i*simDomain.ds;
                y = j*simDomain.ds;
                q = i+j*simDomain.Nx;
                self.f_src[q] = np.exp(-1j*k*(x*np.cos(theta) + y*np.sin(theta)));
                
                
                
        self.f_src = amp*self.f_src;
     
    def get_source(self,simDomain):
        Nx = simDomain.Nx;
        Ny = simDomain.Ny;
        f_src = np.ones((Nx,Ny)).astype(complex);
        for i in range(0,Nx):
            for j in range(0,Ny):
                q = i + j*Nx;
                f_src[i,j] = self.f_src[q];
                
        return f_src;
    
    def plot_source(self,simDomain):
        f_src = np.real(self.get_source(simDomain));
        plt.pcolor(f_src,vmin=-2,vmax=2);
    
    
class FDFD_tfsf:
    def __init__(self,simDomain,tfsf_width):
        Nx = simDomain.Nx; 
        Ny = simDomain.Ny;
        self.Q = scipy.sparse.lil_matrix((Nx*Ny,Nx*Ny));
        self.tfsf_width = tfsf_width; # TFSF width in number of nodes(must be integer). Should be at least 1-2 nodes larger than pml_width
        
    def masking_fnc(self,simDomain):
        Nx = simDomain.Nx; Ny = simDomain.Ny;
        
        # Left edge
        for i in range(0,self.tfsf_width):
            for j in range(0,Ny):
                q = i + j*Nx;
                self.Q[q,q] = 1;
        # Right edge
        for i in range(Nx-self.tfsf_width,Nx):
            for j in range(0,Ny):
                q = i+j*Nx;
                self.Q[q,q] = 1;
        # Bottom edge
        for i in range(0,Nx):
            for j in range(0,self.tfsf_width):
                q = i+j*Nx;
                self.Q[q,q] = 1;
        # Top edge
        for i in range(0,Nx):
            for j in range(Ny-self.tfsf_width,Ny):
                q = i+j*Nx;
                self.Q[q,q] = 1;
                
        self.Q = self.Q.tocsr(); # Convert to compressed sparse row format
    
    def assemble_source(self,A,f_src):
        
        return self.Q.multiply(A[0])-A[0].multiply(self.Q);
        
        
        