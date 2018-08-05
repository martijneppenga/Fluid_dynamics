import numpy as np
import sys
# Flow constants
#---------------
Numtimesteps = 100000        # amount of cycles

N = 400                            # Lattice points in x-direction
M = 200                            # Lattice points in y-direction
maxv = 0.05                        #maximum velocity of the laminar inflow
Re = 150
#v_int = 0.01                      # maximum velocity of Poiseuille flow

#Structure inside the pipe
#---------
structure = "cylinders"                 #choose from "none", "maze", "cylinders","bifurcation" or "square"
save_animation = "tmp.gif"


r = M/10

