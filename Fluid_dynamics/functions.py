import numpy as np
import matplotlib.pyplot as plt
import sys
def equilibrium_function(size_flow_density,rho,w,v,e):
    """equilibrium_function(flow_density,rho,w,v,e):
    Description:
    ------------
    This function calculates the equilibrium distribution functions for the lattice boltzmann model for a d2q9 configuration
    
    Parameters:
    ----------
    size_flow_density: array of size (3,1)
        array with the size of equilibrium distribution functions for each velocity (size_flow_density = i.e (N,M,9) with N and M the 
        size of the fluid flow)
    rho: array of size (N,M)
        array containing the density of the fluid on each location in the flow region
    w: array of size (9,1)   
        array containing the weights for the equilibrium distribution
    v: array of size (N,M,2)
        array containing the x-velocities in  v(:,:,0) and y-velocities in v(:,:,1)
    e: array of size (9,2)
        array containing the velocity vectors with x-component on e(:,0) and y-component on e(:,1) 
        
    Results:
    --------    
    fi: array of size (N,M,9)
        array containing the equilibrium distribution functions ordered according to the e vectors
    """
    
    fi = np.zeros(size_flow_density)
    for i in range(9):
        fi[:,:,i] = w[i]*rho*(1+3*(v[:,:,0]*e[i,0]+v[:,:,1]*e[i,1])+
                             9/2*(v[:,:,0]*e[i,0]+v[:,:,1]*e[i,1])**2-
                             3/2*(v[:,:,0]**2+v[:,:,1]**2))
    
    
    return fi



def velocity_calc(flow_density,e,rho):
    """velocity_calc(flow_density,e,rho)
    Description:
    ------------
    This function calculates the velocity of the flow 
    
    Parameters:
    ----------    
    flow_density: array of size (N,M,9) 
        array containing the equilibrium distribution functions ordered according to the e vectors
    e: array of size (9,2)
        array containing the velocity vectors with x-component on e(:,0) and y-component on e(:,1) 
    rho: array of size (N,M)
        array containing the density of the fluid on each location in the flow region
    
    Results:
    --------
    u: array of size (N,M,2)
        array containing the x-velocities in  v(:,:,0) and y-velocities in v(:,:,1) 
    """
    a = np.shape(flow_density)
    u = np.zeros((a[0],a[1],2))
    
    for i in range(9):
        rho = rho + 1e-16*(rho<1e-16)
        u[:,:,0] += flow_density[:,:,i]*e[i,0]/rho
        u[:,:,1] += flow_density[:,:,i]*e[i,1]/rho
       
    return u

def compute_structure(N,M,structure,r):
    """compute_structure(N,M)
    Description:
    ------------
    This function sets the boundary conditions on the left and right side of a pipe flow
    
    Parameters:
    ----------
    N: integer
        number of lattice points in x-direction
    M: integer
        number of lattice points in y-direction
    structure: string
        String indicating the structure wanted in the pipe

       
    
    Results:
    --------         
    boundary: array of size (N,M)
        array with True on the boundary segments and False on all other sides
    """
    #Location of the square/cylinders if necessary
    xcyl = N/4  
    ycyl = M/2
    [X,Y] = np.meshgrid(np.linspace(0,M-1,M),np.linspace(0,N-1,N),)
    if structure == "cylinders":
        obst = (X-ycyl)**2 + (Y-xcyl)**2<r**2
        obst2 =(X-(ycyl-M/4))**2 + (Y-xcyl)**2<r**2
        obst3 =(X-(ycyl+M/4))**2 + (Y-xcyl)**2<r**2
        boundary = obst*np.ones((N,M)) + obst2*np.ones((N,M)) + obst3*np.ones((N,M))
        if r > M/8:
            plt.imshow(boundary)
            plt.show()
            sys.exit("The radius of the cylinders is too big, now they overlap, please reduce r")
    elif structure == "none":
        boundary = np.zeros((N,M))
            
    elif structure == "maze":
        boundary = np.zeros((N,M))
        boundary[int(N/4):int(N/4+N/20) ,0:int(3*M/4)] = 1
        boundary[int(3*N/4):int(3*N/4+N/20),0:int(3*M/4)] = 1
        boundary[int(N/2):int(N/2+N/20),int(M/4):] = 1
    elif structure == "square":
        obst = (np.abs(Y-xcyl) < r) & (np.abs(X-ycyl)<r)
        boundary = obst*np.ones((N,M))
    elif structure == "bifurcation":
        boundary = np.zeros((N,M))
        boundary[abs(X-M/2)<Y-N/3] = 1
        boundary[:,0:int(M/4)] = 0
        boundary[:,int(3*M/4):M] = 0
    else:
        sys.exit("Please enter a valid structure")


        
    boundary[:,0] = 1                        #define the outer boundaries of the pipe
    boundary[:,-1] = 1                       #define the outer boundaries of the pipe 
    
    return boundary == 1
    

def density_calc(flow_density):
    """"density_calc(flow_density)
    Description:
    ------------
    This function calculates the density of the fluid from the distribution functions of the Lattice boltzmann model 
    
    Parameters:
    ----------
    flow_density: array of size (N,M,9) 
        array containing the equilibrium distribution functions ordered according to the e vectors
    
    Results:
    --------
    rho: array of size (N,M)
        array containing the density of the fluid on each location in the flow region
  
    """
    rho = np.sum(flow_density, axis=2)
    return rho

def flow_fluid(flow_density,e,boundary):
    """flow_fluid(flow_density,e,boundary,bounce)
    Description:
    ------------
    This function flows each of the distribution functions to their new location according to a vector containg the flow direction
    
    Parameters:
    ----------
    flow_density: array of size (N,M,9) 
        array containing the equilibrium distribution functions ordered according to the e vectors
    e: array of size (9,2)
        array containing the velocity vectors with x-component on e(:,0) and y-component on e(:,1) 
    boundary: array of size (N,M)
        array with True on the boundary segments and False on all other sides
        
    Results:
    --------
    fi: array of size (N,M,9) 
        array containing the flowed equilibrium distribution functions ordered according to the e vectors after

    
    """
    fi = np.zeros(np.shape(flow_density))
    
    for i in range(9):
        fi[:, :, i] = np.roll(np.roll(flow_density[:, :, i], e[i, 0], axis=0), e[i, 1], axis=1) 
    #fi[-1,:, :] = 0   #outflow drain
    #fi[0, :, :] = laminarflow
    return fi
     


def bounce_back(flow_density_new,flow_density_old,boundary,bounce):
    """bounce_back(flow_density_new,flow_density_old,boundary,bounce)
    Description:
    -----------
    This function implements the bounce back boundary conditions at the positions given in boundary
        
    Paramters:
    ----------
    flow_density_new: array of size (N,M,9)
        array containing the distribution functions after the relaxation step
    flow_density_old: array of size (N,M,9)
        array containing the distribution functions before the relaxation step        
    boundary: array of size (N,M)
        array with True on the boundary segments and False on all other sides    
    bounce: array of size (2,9) 
        array containing the flow back direction of each distribution function where the direction in the column of row one (0 to 8)
        flows back to the direction in the column of row two. (i.e. [[0,1],[1,0]] means that distribution function 0 flows to 
        distribution function 1 and distribution function 1 flows to distribution function 0) 
    
    Results:
    --------
    flow_density_new: array of size (N,M,9) 
        array containing the distribution functions considering the boundary condition after the relaxation step and 
        """
    
    for i in range(9):
        flow_density_new[boundary, bounce[0,i]] = flow_density_old[boundary,bounce[1,i]]
    return flow_density_new

def relaxation(flow_density,f_eq,tau):
    """relaxation(flow_density,f_eq,tau)
    Description:
    ------------
    This function relax the densities of the distribution function (also known as collision step)
    
    Parameters:
    -----------
    flow_density: array of size (N,M,9) 
        array containing the distribution functions
    f_eq: array of size (N,M,9)
        array containing the equilibrium distribution functions
    tau: float
        relaxation time
    
    Results:
    --------
    f: array of size (N,M,9)
        array containing the relaxed distribution functions
    """
    
    f = flow_density-(flow_density-f_eq)/tau
    return f


def init_velocity(maxv,N,M):
    """init_velocity(maxv,N,M)
    Description:
    ------------
    This function computes the inflow velocity profile (laminar). It sets the velocity profile to
    a poiseuille flow profile in a tube.
    
    Parameters:
    -----------
    N: integer
        number of lattice points in x-direction
    M: integer
        number of lattice points in y-direction
    maxv: float
        number indicating the maximum of the velocity profile
    Results:
    --------
    v: array of size (N,M,2)
        array containing the inflow velocity profile
    """   
    
    
    [X,Y] = np.meshgrid(np.linspace(-1,1,M),np.linspace(-1,1,N),)
    v = (-X**2+1)*maxv
    return v

def rho_zou_he(flow_density,u,rho):
    """rho_zou_he(flow_density,u,rho):
    Description:
    ------------
    This function calculates the density at the inlet required for zou he boundary conditions. The densities are
    scaled to ensure mass conservation
    
    Parameters:
    -----------
    flow_density: array of size (N,M,9) 
        array containing the distribution functions
    u: array of size (M)   
        array containing the velocity profile at the inlet 
    """ 
        
    rho1 = np.sum(rho[0,:])
    rho_boundary = (1/(1-u))*(flow_density[0,:,0]+flow_density[0,:,3]+flow_density[0,:,7]+2*(flow_density[0,:,5]+flow_density[0,:,4]+flow_density[0,:,6]))
    rho[0,:] = rho_boundary/np.sum(rho_boundary)*rho1
    return rho
        
        
def zou_he_boundary(flow_density,f_eq):
    """zou_he_boundary(flow_density,f_eq):
    Description:
    ------------
    This function implements the zou he boundary conditions at the inlet of the structure
    
    Parameters:
    -----------
    flow_density: array of size (N,M,9) 
        array containing the distribution functions
    f_eq: array of size (N,M,9)
        array containing the equilibrium distribution functions
        
    Results:
    --------
    flow_density: array of size (N,M,9)
        array containing the distribution functions corrected for the zou he boundary conditions 
    """
    rho0 = np.sum(flow_density[0,:,[1,2,8]])
    rho1 = np.sum(flow_density[0,:,[5,6,4]] + f_eq[0,:,[1,2,8]] - f_eq[0,:,[5,6,4]])
    
    flow_density[0,:,[1,2,8]] = (flow_density[0,:,[5,6,4]] + f_eq[0,:,[1,2,8]] - f_eq[0,:,[5,6,4]])/rho1*rho0

    return flow_density
    
def outflow_boundary(flow_density):
    """outflow_boundary(flow_density):
    Description:
    ------------
    This function ensures that there is no backward flow in the system by stetting the backflow densities of the 
    of the outlet to the same values as the backward flow at the location just before the outlet. A scaling of the 
    density is performed to ensure that total mass of the system is conserved

    Parameters:
    -----------
    flow_density: array of size (N,M,9) 
        array containing the distribution functions
    
    Results:
    --------
    flow_density: array of size (N,M,9)
        array containing the distribution functions corrected for the outlet 
    """

    flow_density[-1,:,[4,5,6]] = flow_density[-2,:,[4,5,6]]
    

    return flow_density