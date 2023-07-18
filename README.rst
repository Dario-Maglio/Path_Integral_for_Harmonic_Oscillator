=====================================
Path Integral for Harmonic Oscillator
=====================================

In the present repository, we describe the study of the quantum harmonic oscillator through numerical simulation. By adopting the Path Integral approach and implementing the Monte Carlo method, the internal energy of the system is studied as a function of temperature. The wave function of the ground state is determined, and the values of the energy gaps between the ground state and the first two excited levels are obtained.

Repository Structure
====================

The structure of the repository is as follows:

- ``class_lattice.h``:

  This header, located in the ``include`` directory, contains the lattice class from which the trajectory lattice is instantiated. The class includes a Metropolis update method and an incorporated Pseudo Random Number Generator.

- ``harmonic_main.cpp``:

  This program calls the simulation subroutine for several values of the temperature and sides, collecting internal energy, positions and correlators measures in the ``Data_simulation`` folder.

- ``harmonic_plot.py``:

  This program utilizes the data in the ``Data_simulation`` folder to plot the physical quantities of interest.

- ``Plots_and_fit``:

  All produced plots are stored in this folder.

- ``Tests``:

  This directory contains easy-to-use examples for testing the class methods and verifying that the Monte Carlo algorithm has thermalized.


Analysis Results
================

Here are some of the plots generated from the analysis:

- Monte Carlo trajectories and energies:

  .. image:: https://github.com/Dario-Maglio/Path_Integral_for_Harmonic_Oscillator/blob/54dc4be0294df21678a78ab28b849ae03f2e6852/Tests/test_montecarlo.png
     :align: center
     :width: 80%


- Internal energy as a function of the temperature:

  .. image:: https://github.com/Dario-Maglio/Path_Integral_for_Harmonic_Oscillator/blob/54dc4be0294df21678a78ab28b849ae03f2e6852/Plots_and_fit/Energy%20as%20a%20function%20of%20beta.png
     :align: center

- Ground state wavefunction:

  .. image:: https://github.com/Dario-Maglio/Path_Integral_for_Harmonic_Oscillator/blob/54dc4be0294df21678a78ab28b849ae03f2e6852/Plots_and_fit/GS%20%7C%20beta%20%3D%2050%20%2C%20side%20%3D%20260.png
     :align: center

- Two times correlator with t = nk :

  .. image:: https://github.com/Dario-Maglio/Path_Integral_for_Harmonic_Oscillator/blob/54dc4be0294df21678a78ab28b849ae03f2e6852/Plots_and_fit/Correlator%201%20%7C%20Beta%20%3D%2050.png
     :align: center


Feel free to explore the repository and use the provided programs for further analysis and investigation.

License
=======

This repository is licensed under the GNU General Public License v3.0 (GPL-3.0). 

See the LICENSE file for more information.
