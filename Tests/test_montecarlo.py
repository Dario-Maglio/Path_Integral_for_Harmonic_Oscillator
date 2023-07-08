import os

import numpy as np
import matplotlib.pyplot as plt

side = 480
beta = 4
DIM = 1000

#-------------------------------------------------------------------------------

def plot_metropolis():

    fig, axes = plt.subplots(1, 2, num="energy", figsize=(12, 12))

    axes[0].set_title(f"Renormalized energy | Beta: {beta} , Eta: {beta/side}")
    axes[0].set_ylabel(r'$U_{ren}$')
    axes[0].set_xlabel('$step$')

    axes[1].set_title("?")
    axes[1].set_ylabel(r'$ ? $')
    axes[1].set_xlabel(r'$step$')

    directory = f"../Data_simulations/Beta_{beta}"
    filename = f"side_{side}.dat" #_beta_{1:.6f}.dat".format(side, beta)
    file = os.path.join(directory, filename)
    print("Loading file " + file)

    x = [ i for i in range(DIM)]

    if os.path.isfile(file):
        ene = np.loadtxt(file, unpack='True')
        axes[0].plot(x, ene[:DIM])

    print("\nPlots of energy and magnetization: \n")
    plt.show()


#-------------------------------------------------------------------------------

if __name__ == '__main__':

    plot_metropolis()
