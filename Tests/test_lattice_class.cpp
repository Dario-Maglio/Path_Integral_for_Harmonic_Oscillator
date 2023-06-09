/*******************************************************************************
*
* Test program for the Ising simulation
*
*******************************************************************************/

// g++ -I ../include test_lattice_class.cpp -o t_lattice.out

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>

#include <chrono>

// Import the Class lattice
#include "class_lattice.h"

using namespace std;

/*******************************************************************************
* PARAMETERS OF THE LATTICE
*
* SIDE = size of the lattice's side.
*
* G_FLAG = lattice's geometry flag.
*          0 and others not implemented yet.
*          1 for 1D periodic chain.
*          2 for 2D square with PBC.
*
* I_FLAG = lattice's initial configuration flag:
*          0 for cold initialization.
*          1 for hot (random) initialization.
*          2 for loading the previous configuration.
*
*******************************************************************************/

#define SIDE 260
#define G_FLAG 1
#define I_FLAG 0

/*******************************************************************************
* PARAMETERS OF THE SIMULATION
*
* BETA = adimensional reciprocal of the product of the temperature and k_B.
*
* EXTFIELD = adimensional intensity of the external magnetic field.
*
* I_DECORREL = MC-updates of the lattice between different measurements.
*
* LOOPS = timed MC-updates
*
*******************************************************************************/

#define ETA 1.
#define ETA_SEP 0.005
#define LOOPS 10000 // 100*I_DECORREL

//--- Main Test ----------------------------------------------------------------

int main(){
    /* Test the methods of the lattice Class. */

    int up_tot;
    string file_name;
    vector<double> lattice_old;
    ofstream file_traj;
    lattice trajectory(SIDE, G_FLAG, I_FLAG);

    // Initialise the lattice
    cout << "Running with eta = " << ETA << endl;
    file_name = "test_traje";
    if(I_FLAG == 2){
        // Loading configuration
        trajectory.load_configuration(file_name);
    } else {
        // Thermalization phase
        for(int i = 0; i < (1000); i++) trajectory.update(ETA);
    }
    cout << "Ready to go!" << endl;

    // Show initial configuration
    // trajectory.show_configuration();
    // Show the list of all nearest neighbors to each lattice site
    // trajectory.show_nearest_neighbors();

    // Test time MC updates
    auto start = chrono::steady_clock::now();
    for(int i = 0; i < LOOPS; i++) trajectory.update(ETA);
    auto end = chrono::steady_clock::now();
    // Print elapsed time
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time for " << LOOPS << " loops: ";
    cout << elapsed_seconds.count() << "s " << endl << endl;

    // Test efficiency MC updates and save trajectories
    file_traj.open(file_name + "_conf.dat");
    for(double eta_val = ETA; eta_val > 0.001; eta_val -= ETA_SEP) {
        up_tot = 0;
        for(int i = 0; i < LOOPS; i++){
            lattice_old = trajectory.latt_conf;
            trajectory.update(eta_val);
            for(int j = 0; j < SIDE; j++)
               if(trajectory.latt_conf[j] != lattice_old[j]) up_tot++;
        }
        cout << "Eta : " << eta_val << " | " << "Probability ratio: ";
        cout << double(up_tot) / (LOOPS * SIDE) << endl;

        file_traj << eta_val << " ";
        for(auto val : lattice_old) file_traj << val << " ";
        file_traj << endl;
    }
    file_traj.close();

    // Show final configuration
    // trajectory.show_configuration();
    // Save final configuration
    trajectory.save_configuration(file_name);
}
