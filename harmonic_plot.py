"""*****************************************************************************
*
* Plot program for the outcomes of the simulation
*
*****************************************************************************"""

import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats

#*******************************************************************************
# PARAMETERS OF THE SIMULATION
#
# SIDE_SEP = separation between the sides of different simulations.
#
#*******************************************************************************

CORREL_LENGTH = 40
# beta for ground state
BETA_GS = 50
# simulated sides
SIDE_SEP = 40
SIDE_MIN = 20
SIDE_MAX = 500
# simulated betas
betas = [1, 2, 3, 4, 7, 10, 15, 20]

plt.style.use('ggplot')

sides = np.arange(SIDE_MIN, SIDE_MAX+1, SIDE_SEP, dtype='int')

#--- Contents ------------------------------------------------------------------

def fit_fun(x, A, B):
    y = A + B * np.power(x, 2)
    return y

def fit_ene(x, A, B):
    y = A + B / (np.exp(x) - 1)
    return y

def fit_gap(x, A, B):
    y = B * np.exp(-np.multiply(A, x))
    return y


#--- Procedures ----------------------------------------------------------------

def plot_wavefun(beta, side):

    #---Load data
    filename = f"collect_side_{side}.dat"
    file_path = os.path.join("Data_simulations", f"Beta_{beta}")
    file_path = os.path.join(file_path, filename)
    print("Loading " + file_path)
    # load data from each side file
    if os.path.isfile(file_path):
        x = np.loadtxt(file_path, unpack='True')
        x = np.reshape(x, -1)
    else:
        raise FileNotFoundError("Error: File not found!")

    #---Plot
    title = f"GS | beta = {beta} , side = {side}"
    print("\nPlot " + title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.title(title)
    plt.ylabel('$ Pr(x) $')
    plt.xlabel('$ x $')
    # data points
    binwidth = 0.1
    n_bins = np.arange(-4, 4, binwidth)
    lab = 'MC positions'
    plt.hist(x, n_bins, density=True, ls='dashed', alpha=0.5, lw=3, label=lab)
    # expected distribution
    binwidth = 0.01
    n_bins = np.arange(-4, 4, binwidth)
    lab = 'GS wave func'
    y = np.random.normal(0, np.sqrt(0.5), 10000000)
    plt.hist(y, n_bins, density=True, ls='solid', alpha=0.5, lw=3, label=lab)
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()



def plot_energy(beta):

    #---Load data
    filename = "energy_data.dat"
    file_path = os.path.join("Data_simulations", f"Beta_{beta}")
    file_path = os.path.join(file_path, filename)
    print("Loading " + file_path)
    # load data from each side file
    if os.path.isfile(file_path):
        x, y, y_err = np.loadtxt(file_path, unpack='True')
    else:
        raise FileNotFoundError("Error: File not found!")

    #---Fit
    parameters, covariance = curve_fit(fit_fun, x, y, sigma=y_err)
    fit_a = parameters[0]
    fit_b = parameters[1]
    std_deviation = np.sqrt(np.diag(covariance))
    fit_da = std_deviation[0]
    fit_db = std_deviation[1]
    print("\nFit parameters:")
    print(f"{fit_a} ± {fit_da}\n{fit_b} ± {fit_db}")
    # reduced chi squared
    fit_y = fit_fun(x, *parameters)
    chisq = np.sum(np.power(((y - fit_y)/ y_err), 2))
    chisqrd = chisq / (len(x) - 3)
    print(f"Reduced chi squared: {chisqrd}")

    #---Plot
    title = f"Energy | beta = {beta}"
    print("\nPlot " + title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.title(title)
    plt.ylabel(r'$ U(\eta) $')
    plt.xlabel(r'$ \eta $')
    # points and function
    fit_x = np.linspace(min(x), max(x), 100)
    fit_y = fit_fun(fit_x, *parameters)
    fit_label = f'fit E_c = {fit_a:.4f} ± {fit_da:.4f}'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    sim_label = f'collected data'
    plt.errorbar(x, y, yerr=y_err, fmt='<',label=sim_label)
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

    return (fit_a , fit_da)

def plot_energy_beta(x, y, y_err):

    #---Fit
    parameters, covariance = curve_fit(fit_ene, x, y, sigma=y_err, bounds=((0, 0), (5, 5)))
    fit_a = parameters[0]
    fit_b = parameters[1]
    std_deviation = np.sqrt(np.diag(covariance))
    fit_da = std_deviation[0]
    fit_db = std_deviation[1]
    print("\nFit parameters:")
    print(f"{fit_a} ± {fit_da}\n{fit_b} ± {fit_db}")
    # reduced chi squared
    fit_y = fit_ene(x, *parameters)
    chisq = np.sum(np.power(((y - fit_y)/ y_err), 2))
    chisqrd = chisq / (len(x) - 3)
    print(f"Reduced chi squared: {chisqrd}")

    #---Plot
    title = f"Energy as a function of beta"
    print("\nPlot " + title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.title(title)
    plt.ylabel(r'$ U_{ren}(\beta) $')
    plt.xlabel(r'$ N \eta $')
    # points and function
    fit_x = np.linspace(min(x), max(x), 100)
    fit_y = fit_ene(fit_x, *parameters)
    fit_label = f'fit E_0 = {fit_a:.6f} ± {fit_da:.6f}'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    sim_label = f'computed energy'
    plt.errorbar(x, y, yerr=y_err, fmt='<',label=sim_label)
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_correlators(beta, label):

    eta_val = []
    gap_val = []
    gap_err = []

    #--- Load data
    file_path = os.path.join("Data_simulations", f"Beta_{beta}")
    file = os.path.join(file_path, f"correlations_{label}_data.dat")
    # load data from file
    print("Loading " + file)
    if os.path.isfile(file):
        data = np.loadtxt(file, unpack='True')
    else:
        raise FileNotFoundError("Error: File not found!")

    #--- Fit and plot correlator
    title = f"Correlator {label} | Beta = {beta}"
    print("\nPlot " + title + "\n")
    fig = plt.figure(title)
    plt.title(title)
    plt.ylabel(r'$ C(k) $')
    plt.xlabel(r'$ k $')

    for n in range(data.shape[1]):
        row = data[:, n]
        eta = row[0]
        row = row[1:]
        k = [(i + 1) for i in range(CORREL_LENGTH)]
        x = [val*eta for val in k]
        y_val = [row[2*i] for i in range(CORREL_LENGTH)]
        y_err = [row[2*i+1] for i in range(CORREL_LENGTH)]

        # fit data
        parameters, covariance = curve_fit(fit_gap, x, y_val, sigma=y_err)
        fit_a = parameters[0]
        fit_b = parameters[1]
        std_deviation = np.sqrt(np.diag(covariance))
        fit_da = std_deviation[0]
        fit_db = std_deviation[1]
        print("\nFit parameters:")
        print(f"{fit_a} ± {fit_da}\n{fit_b} ± {fit_db}")
        # reduced chi squared
        fit_y = fit_gap(x, *parameters)
        chisq = np.sum(np.power(((y_val - fit_y)/ y_err), 2))
        chisqrd = chisq / (len(x) - 3)
        print(f"Reduced chi squared: {chisqrd}")
        # save fit results
        eta_val.append(eta)
        gap_val.append(fit_a)
        gap_err.append(fit_da)
        # plot data
        sim_label = f'side {int(beta/eta + 0.5)}'
        plot = plt.errorbar(k, y_val, y_err, fmt='<', label=sim_label)
        color = plot[0].get_color()
        fit_x = np.linspace(min(x), max(x), 100)
        fit_y = fit_gap(fit_x, *parameters)
        plt.plot(fit_x/eta, fit_y, '-', color=color)

    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

    #--- Fit and plot energy gap
    parameters, covariance = curve_fit(fit_fun, eta_val, gap_val, sigma=gap_err)
    fit_a = parameters[0]
    fit_b = parameters[1]
    std_deviation = np.sqrt(np.diag(covariance))
    fit_da = std_deviation[0]
    fit_db = std_deviation[1]
    print("\n---Fit parameters Energy Gap:")
    print(f"{fit_a} ± {fit_da}\n{fit_b} ± {fit_db}")
    # reduced chi squared
    fit_y = fit_fun(eta_val, *parameters)
    chisq = np.sum(np.power(((gap_val - fit_y)/ gap_err), 2))
    chisqrd = chisq / (len(eta_val) - 3)
    print(f"Reduced chi squared: {chisqrd}")
    # plot details
    title = f"Energy gap {label} | Beta = {beta}"
    print("\nPlot " + title + "\n")
    fig = plt.figure(title)
    plt.title(title)
    plt.ylabel(r'$ \Delta E (\eta) $')
    plt.xlabel(r'$ \eta $')
    # plot data
    fit_x = np.linspace(min(eta_val), max(eta_val), 100)
    fit_y = fit_fun(fit_x, *parameters)
    plt.plot(fit_x, fit_y, '-', label=f'Fit gap = {fit_a:.4f} ± {fit_da:.4f}')
    plt.errorbar(eta_val, gap_val, yerr=gap_err, fmt='.')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    plot_wavefun(BETA_GS, int((SIDE_MIN + SIDE_MAX)/2))

    plot_correlators(BETA_GS, 1)

    plot_correlators(BETA_GS, 2)

    ene_val = []
    ene_err = []
    for beta in betas:
        y_val, y_err = plot_energy(beta)
        ene_val.append(y_val)
        ene_err.append(y_err)

    plot_energy_beta(betas, ene_val, ene_err)
