import os

import numpy as np
import matplotlib.pyplot as plt

side = 260
beta = 20
DIM = 1000

#-------------------------------------------------------------------------------

def plot_metropolis():

    fig, axes = plt.subplots(1, 2, num="energy", figsize=(12, 12))

    #--- Trajectory ------------------------------------------------------------

    axes[0].set_title(f"Trajectories | Beta: {beta}")
    axes[0].set_ylabel('$ x $')
    axes[0].set_xlabel('step')

    file = "test_traje_conf.dat"

    if os.path.isfile(file):
        print("Loading file " + file)
        traje = np.loadtxt(file, unpack='True')

    x = [ i for i in range(side)]
    t = [ 0, 100, 180, 199]
    for i in range(4):
        k = t[i]
        y = [traje[j+1][k] + 3*i for j in range(side)]
        axes[0].plot(x, y, label=f"Eta: {traje[0][k]}")

    axes[0].legend(loc="upper right")


    #--- Energy ----------------------------------------------------------------

    axes[1].set_title(f"Renormalized energy | Beta: {beta} , Eta: {beta/side:.4f}")
    axes[1].set_ylabel(r'$U_{ren}$')
    axes[1].set_xlabel('step')

    directory = f"../Data_simulations/Beta_{beta}"
    filename = f"side_{side}.dat" #_beta_{1:.6f}.dat".format(side, beta)
    file = os.path.join(directory, filename)

    if os.path.isfile(file):
        print("Loading file " + file)
        ene = np.loadtxt(file, unpack='True')

    x = [ i for i in range(DIM)]
    y = ene[:DIM]
    axes[1].plot(x, y)

    plt.show()


#-------------------------------------------------------------------------------

if __name__ == '__main__':

    plot_metropolis()
