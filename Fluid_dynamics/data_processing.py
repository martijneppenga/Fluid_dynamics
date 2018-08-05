import numpy as np
import matplotlib.pyplot as plt
import datetime
from constants import *
def start(E,px,py,rho_tot):
    

    plt.plot(E)
    plt.title('Energy versus Time')
    plt.xlabel('t (iterations)')
    plt.ylabel('E')
    plt.show()
    
    plt.plot(px)
    plt.title('X-component of momentum versus Time')
    plt.xlabel('t (iterations)')
    plt.ylabel('px')
    plt.show()
    
    plt.plot(py)
    plt.title('Y-component of momentum versus Time')
    plt.xlabel('t (iterations)')
    plt.ylabel('py')
    plt.show()

    plt.plot(rho_tot)
    plt.title('Density versus Time')
    plt.xlabel('t (iterations)')
    plt.ylabel('rho')
    plt.show()