/*******************************************************************************
*
* Main program for the computation of the ground state energy and wave function
*
*******************************************************************************/

// g++ -I ./ harmonic_main.cpp -o harmonic.out

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

// Import the Class lattice
#include "class_lattice.h"
// Import the bootstrap procedure
#include "data_analysis.h"

using namespace std;

/*******************************************************************************
* PARAMETERS OF THE SIMULATION
*
* BETA = 1 / (k_b * Temperature) := (SIDE * ETA) / (hbar * omega) .
*
* ETA = omega * lattice spacing .
*
* SIDE = (hbar * BETA) / lattice spacing := Time period / lattice spacing .
*
* MEASURES = desired number of measures for each value of SIDE.
*
* I_DECORREL = MC-updates of the lattice between different measurements.
*
* I_FLAG = lattice's initial configuration flag:
*          0 for cold initialization.
*          1 for hot (random) initialization.
*          2 for loading the previous configuration and append data.
*
* G_FLAG = lattice's geometry flag:
*          0 and others not implemented yet.
*          1 for 1D periodic chain.
*          2 for 2D square with PBC.
*
*******************************************************************************/

// sides to simulate
#define SIDE_SEP 40
#define SIDE_MIN 40
#define SIDE_MAX 480
// fix the Temperature
#define BETA_MAX 20
// number of measures to save
#define MEASURES 1000
// decorrelation between measures
#define I_DECORREL 100 // * V
// initialization flags
#define I_FLAG 2
#define G_FLAG 1

using namespace std;

const vector<int> betas = {1, 2, 3, 4, 7, 10, 15, BETA_MAX};

//--- Contents -----------------------------------------------------------------

void wave_function(int side, int beta){
    /* Collect x-coordinates for a given side and beta */

    ofstream file;
    string directory, name_file, name_file_data, name_file_state;
    double eta = double(beta) / side;
    lattice traj(side, G_FLAG, I_FLAG);

    // Define path data directory
    directory = "Data_simulations/Beta_" + to_string(beta) + "/";
    // Define name file last configuration of the lattice
    name_file_state = "collect_side_" + to_string(side);
    // Define name file simulation for a given side and beta
    name_file_data =  name_file_state + ".dat";
    cout << "Creation of " << name_file_data << " begins..." << endl;

    // Prepare the lattice and the data file for the simulation
    if(I_FLAG == 2){
        // Load last configuration of the lattice
        traj.load_configuration(directory + name_file_state);
        // Open the existing data file
        file.open(directory + name_file_data, ios_base::app);
    } else {
        // Thermalization phase
        for(int i = 0; i < (100*I_DECORREL); i++) traj.update(eta);
        // Open a new data file
        file.open(directory + name_file_data);
    }
    // Update ising and take measures
    for(int n = 0; n < MEASURES; n++){
        for(int i = 0; i < I_DECORREL; i++) traj.update(eta);
        for(auto val : traj.latt_conf) file << val << endl;
    }
    file.close();

    // We can stop the simulation when a file is completed
    traj.save_configuration(directory + name_file_state);
    cout << "Creation of " << name_file_data << " completed." << endl << endl;
}

void run_simulation(int side, int beta){
    /* MC-simulation for a given side and beta */

    ofstream file;
    string directory, name_file, name_file_data, name_file_state;
    double eta = double(beta) / side;
    lattice traj(side, G_FLAG, I_FLAG);

    // Define path data directory
    directory = "Data_simulations/Beta_" + to_string(beta) + "/";
    // Define name file last configuration of the lattice
    name_file_state = "side_" + to_string(side);
    // Define name file simulation for a given side and beta
    name_file_data =  name_file_state + ".dat";
    cout << "Creation of " << name_file_data << " begins..." << endl;

    // Prepare the lattice and the data file for the simulation
    if(I_FLAG == 2){
        // Load last configuration of the lattice
        traj.load_configuration(directory + name_file_state);
        // Open the existing data file
        file.open(directory + name_file_data, ios_base::app);
    } else {
        // Thermalization phase
        for(int i = 0; i < (100*I_DECORREL); i++) traj.update(eta);
        // Open a new data file
        file.open(directory + name_file_data);
    }
    // Update ising and take measures
    for(int n = 0; n < MEASURES; n++){
        for(int i = 0; i < I_DECORREL; i++) traj.update(eta);
        file << traj.energy(eta) << endl;
    }
    file.close();

    // We can stop the simulation when a file is completed
    traj.save_configuration(directory + name_file_state);
    cout << "Creation of " << name_file_data << " completed." << endl << endl;
}

void beta_simulation(int beta){
    /* Iterates the timed simulation with fixed beta over sides */

    // Initialize timing variables
    auto start = chrono::steady_clock::now(), end = start;
    chrono::duration<double> elapsed_sec = end - start;
    // Start simulation
    cout << "Beta = " << beta << endl;
    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){
        start = chrono::steady_clock::now();
        run_simulation(side, beta);
        end = chrono::steady_clock::now();
        elapsed_sec = end - start;
        cout << "Elapsed time: " << elapsed_sec.count() << "s " << endl << endl;
    }
}

void energy_analysis(int beta){
    /* Error analysis of the energy for a given beta */

    string file_path, file_name;
    ofstream file_data, file_analysis;
    vector<double> data;

    // create the data and logging files
    file_path = "Data_simulations/Beta_" + to_string(beta) + "/energy";
    file_name =  file_path + "_analysis.txt";
    cout << "Creating file: " << file_name << endl;
    file_analysis.open(file_name);
    file_name = file_path + "_data.dat";
    cout << endl << "Creating file: " << file_name << endl << endl;
    file_data.open(file_name);
    // begin loop over sides
    for(double side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){
        data = file_operations(side, beta, file_analysis);
        file_data << double(beta) / side << " ";
        for(auto val : data) file_data << val << " ";
        file_data << endl;
    }

    file_data.close();
    file_analysis.close();
}

//--- Main ---------------------------------------------------------------------

int main(){

    cout << "Collecting positions for beta = " << BETA_MAX  << " and side = ";
    cout << 200 << endl;
    wave_function(200, BETA_MAX);

    cout << "--- Starting simulations... " << endl;
    for(auto beta : betas)
        beta_simulation(beta);

    cout << "--- Starting error analysis... " << endl;
    for(auto beta : betas)
        energy_analysis(beta);

    cout << "The work is done." << endl << endl;
    return 0;
}
