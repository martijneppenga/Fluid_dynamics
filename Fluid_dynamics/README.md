# Lattice Boltzmann Simulation for fluid flow d2q9 

This third project for for the computational physics course pertains the Lattice Boltzmann Method (LBM) for fluid dynamics. 
In this simulation, the d2q9 LBM is used to simulate the flow of a liquid in a straight pipe, with several options for structures in the pipe. This simulation uses bounce-back boundary conditions at the walls and Zou-He velocity boundary conditions for the inflow and outflow.


## Getting Started
In order to make an animation, simply run "animation.py"

In order to run the simulation to extract physical quantities and checks, copy all the .py files in the repository to your working directory. Change the input constants in constants.py, save the file, then start the simulation from the simulation module, and use the results to do the processing, as follows:

```
import importlib
import simulation
import data_processing

importlib.reload(simulation)
results = simulation.start()
importlib.reload(data_processing)
data_processing.start(*results)
```



### Options for constants

The following inputs can be adjusted in constants.py

Parameters of the system
* **N** the number of pixels in the direction of the flow
* **M** the number of pixels tranverse to the direction of the flow
* **maxv** the maximum speed of the Poiseuille Flow in the center of the pipe
* **structure** the type of structure to be present in the pipe. Options are "none", "square", "cylinders", "maze"or "bifurcation"
* **r** half of the diameter of the square or cylinders if "square" or "cylinders" was chosen
* **Re** The Reynolds number of the system, M * D / nu, where D is defined as the diameter of cylinders or square if "cylinders" or "square", or the diameter of the pipe if a different structure was chosen. Nu is the kinematic viscosity.
* **save_animation** The path for animation.py to save the .gif of the animation



### Obtained results
The output of the simulation contains the following, all meant to test conservation laws (energy, momentum, mass):
* **E**         The energy of the system
* **px, py**    The total momentum of the system in the x- and y-direction respectively
* **rho_tot**   The density of the system summed over all pixels


### Other folders (GIFs, Checks)
In the repository two other folders can be found, GIFs and Checks.

#### GIFs 
contains example animations for several structures, for different Reynolds numbers, resolutions, and structures. The naming in the files is done in the following way: structure_widthxheight_ReReynoldsnumber.
##### Square and Cylinders
For the "square" and "cylinders" structure, it is seen that for low Reynolds numbers, a stable laminar flow arises as expected. For higher Reynolds numbers, the flow gets more turbulent, also in the line of expectation. 
Furthermore, the animation was made resolution-invariant, which means that the time-scale and patterns of the flow are similar for every resolution for different Reynolds number. This is verified by seeing that for different resolutions, the gifs for the "cylinders" structures are similar for different resolutions.
Also check out the high resolution, 60FPS gif (permission for using it as screensaver is granted)
##### Maze 
For the maze structure, it is seen that for low Reynolds numbers a steady laminar flow arises, while for higher Reynolds numbers, the flow is more turbulent as expected. Unfortunately, simulating situations with even more turbulence is beyond the capabilities of the simulation.

##### Bifurcation
Within the possible parameters of our simulation, the flow in the bifurcation is always laminar. It is seen that the parabolic inflow results in two parabolic outflows with higher peaks due to mass conservation.

##### Equilibrium of the simulation
The simulation is initialised in such a way, that the velocity has a parabolic profile everywhere except inside the simulated structures (so even very close to the structures themselves). This is an non-physical situation, which creates shockwaves, clearly seen in "shockwave.gif".
After some time, these shock waves die out and the flow enters a form of equilibrium, which does not depend on the initial conditions (but only the boundary conditions of the in- and outflow).


#### Checks
Checks contains two folders:
* /conservation\_laws: plots are given for an example simulation of the energy, momentum, and mass of the system to show that they are conserved in the case of the "cylinders" structure. In the conservation of energy and momentum, it is clearly seen (as stated before) that the simulation does not start in an equilibrium. As seen in the Mass_conservation.png, mass is not perfectly conserved, probably due to model errors. However, as seen in Mass_conservation_rescaled, this is only a small deviation with respect to the total density
* /inflow\_type: a simulation is done in a simple pipe for a uniformly distributed inflow, aswell as a Poiseuille (parabolic) inflow. The comparison shows that at the end of the pipe, the profile is the same. This justifies the use of the Poisseuile inflow for all simulation to simulate a pipe that was extended at the start.

### Performance
Since all non-necessary for-loops were vectorised in this project and numpy was used for the matrix operations, it was seen that our simulation was able to compute high-resolution data (a grid of 800x400, with 120 000 timesteps takes aroundtakes 7.5h, meaning at least 1.4 million operations per second). Unfortunately there was no time for any checks on grid-size dependence of the computation time.

# Distribution of tasks
Once again, we mostly programmed with the three of us at the same time, so most of the simulation has been done by all of us, while brainstorming about how to tackle problems along the way. However, there was a different focus for all of us:
* Martijn focused on the boundary conditions and in and out flow
* Cyrus focused on some simulation structures and structuring the data
* Richard focused on making the animation and the simulation structures


## Authors
* Martijn Eppenga
* Richard Faasse
* Cyrus Tirband 

