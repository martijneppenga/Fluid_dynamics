import time
import numpy as np
import importlib
import functions
importlib.reload(functions)
from functions import *
import matplotlib.pyplot as plt
from matplotlib import cm
from IPython import display
import constants
importlib.reload(constants)
from constants import *

def start():
    
    #compute Reynolds number for the structures
    if structure == "none" or structure == "maze" or structure == "bifurcation":
        nu = M*maxv/Re
    elif structure == "cylinders" or structure == "square":
        nu = 2*r*maxv/Re
    else:
        sys.exit("Please choose a valid structure")

    #structural constants for D2Q9 configuration     
    tau      = (6*nu + 1)/2                                                                #relaxation parameter
    w        = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36])                         #weights for the different directions
    e        = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]])     #directions in which the lquid can flow [y,x]
    bounce = np.array([[0,1,2,3,4,5,6,7,8],[0,5,6,7,8,1,2,3,4]])                           #Array that pairs the opposite directions of e

    #initialize boundary positions, density, velocity profile, boundary velocity profile (u) and density for each flow direction (flow_density)
    boundary = compute_structure(N,M,structure,r)
    rho = np.ones((N,M))
    velocity = np.zeros((N,M,2))
    velocity[:,:,0] = init_velocity(maxv,N,M)
    velocity[boundary,:] = 0
    u = np.copy(velocity)
    flow_density = equilibrium_function((N,M,9),rho,w,velocity,e)
    
    #initialise variables for verification
    rho_tot = np.zeros(Numtimesteps)
    E = np.zeros(Numtimesteps)
    px = np.zeros(Numtimesteps)
    py = np.zeros(Numtimesteps)


    

    start = time.time()
    ## Main loop ##
    for t in range(Numtimesteps):
        
        #compute values of verfication variables
        rho_tot[t] = np.sum(rho)
        E[t]= np.sum((rho*(velocity[:,:,0] ** 2 + velocity[:,:,1] ** 2)))
        px[t]=np.sum(rho*velocity[:,:,0])
        py[t]=np.sum(rho*velocity[:,:,1])
        
        #compute flow at the outlet
        flow_density = outflow_boundary(flow_density)
        
        #compute the density and velocity of the flow
        rho = np.sum(flow_density, axis=2)
        velocity = velocity_calc(flow_density,e,rho)
        velocity[0,:,:] = u[0,:,:]

        #compute flow at the inlet and equilibrium distribution
        rho = rho_zou_he(flow_density,u[0,:,0],rho)
        f_eq = equilibrium_function((N,M,9),rho,w,velocity,e)
        flow_density = zou_he_boundary(flow_density,f_eq)


        #compute the relaxation of the flow densities (collision step) 
        flow_density_new = relaxation(flow_density,f_eq,tau)
        
        #compute the bounce back boundary conditions at the boundary
        flow_density_new = bounce_back(flow_density_new,flow_density,boundary,bounce)


        #flow each of the flow densities to the new positions
        flow_density = flow_fluid(flow_density_new,e,boundary)
        
        print("progress is: "+str(np.round(t/Numtimesteps*100,decimals=2))+"%", end = '\r')


    end = time.time()
    dt = end-start
    print("Time elapsed", dt)
    return E,px,py,rho_tot