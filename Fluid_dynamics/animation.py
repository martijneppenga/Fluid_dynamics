import numpy as np
import matplotlib.pyplot as plt

from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider
import sys
import matplotlib.patches as mpatches
from matplotlib import cm
from IPython import display
from functions import *
from constants import *

#compute Reynolds number for the structures
if structure == "none" or structure == "maze" or structure == "bifurcation":
    nu = M*maxv/Re
elif structure == "cylinders" or structure == "square":
    nu = 2*r*maxv/Re
else:
    sys.exit("Please choose a valid structure")

#structural constants for D2Q9 configuration 
tau = (6*nu + 1)/2                                                            #Relaxation parameter
w   = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36])                     #weights for the different directions
e   = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]]) #directions in which the liquid can flow [y,x] 
bounce = np.array([[0,1,2,3,4,5,6,7,8],[0,5,6,7,8,1,2,3,4]])                  #Array that pairs the opposite directions of e




                           

#initialize boundary positions, density, velocity profile, boundary velocity profile (u) and density for each flow direction (flow_density)
boundary = compute_structure(N,M,structure,r)
rho = np.ones((N,M))
velocity = np.zeros((N,M,2))
velocity[:,:,0] = init_velocity(maxv,N,M)
velocity[boundary,:] = 0
u = np.copy(velocity)
flow_density = equilibrium_function((N,M,9),rho,w,velocity,e)


#initialize figure for the animation
fig, ax = plt.subplots(1, figsize=(8,4), facecolor=(1,1,1))
fig.subplots_adjust(left=0, right=1, bottom=0)
xx, yy = np.meshgrid(np.linspace(-2,3,500), np.linspace(-1,2,500))

def make_frame(t):
    global flow_density
    for i in range(int(100*M/100)):
        
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
        
        #compute the velocity profile
        vabs =  np.sqrt(velocity[:,:,0]**2 + velocity[:,:,1]**2)
        vabs[boundary] = float('NaN')
    
    
    ax.clear()
    ax.axis('off')
    ax.set_title("Velocity Profile", fontsize=16)    

    #the varying weights make the points appear one after the other
    ax.imshow(vabs.transpose())

    

    return mplfig_to_npimage(fig)

animation = VideoClip(make_frame, duration = 10)
animation.write_gif(save_animation, fps=30)