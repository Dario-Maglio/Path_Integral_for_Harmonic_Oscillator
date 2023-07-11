/*******************************************************************************
*
* Lattice Class definition
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#ifndef LATTICE_CLASS_H
#define LATTICE_CLASS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <cmath>

using namespace std;

random_device device;

//--- Contents -----------------------------------------------------------------

class lattice {

private:
    int tot_lenght_;
    vector<vector<int>> nearest_neighbors_;

    // Define the PRNG
    mt19937_64 generator_;

public:
    vector<double> latt_conf;
    const int side_lenght, geometry_flag, initial_flag;

    lattice(const int &SIDE, const int &G_FLAG, const int &I_FLAG):
       /* Class lattice constructor **********************************

       Inputs:

       SIDE = number of sites per side.

       G_FLAG =  geometry_flag for the topology.
                 1 for 1D lattice with PBC.
                 2 for 2D square lattice with PBC.
                 else gives "Error: not implemented yet".

       I_FLAG = initial_flag for the initial condition.
                0 for initialising cold lattice.
                else initialise hot lattice.

        *************************************************************/
        side_lenght(SIDE),
        geometry_flag(G_FLAG),
        initial_flag(I_FLAG)
        {// CONSTRUCTION BEGIN

        mt19937_64 generator_(device());

        double random_number;
        vector<int> nearest_list;

        // Defining topology and nearest neighbors
        if (geometry_flag == 1) {
            tot_lenght_ = side_lenght;
            latt_conf.reserve(tot_lenght_);
            nearest_neighbors_.reserve(tot_lenght_);
            nearest_list.reserve(2);

            // Fill nearest_neighbors_ with nearest_list for each site
            for (int i = 0; i < tot_lenght_; i++){
                // Fill nearest_list for site i
                nearest_list.clear();
                nearest_list.push_back(i - 1);
                nearest_list.push_back(i + 1);

                // Implementing PBC
                if (i == 0){
                    // Correct nearest_list first element
                    nearest_list[0] = tot_lenght_-1;
                } else if (i == tot_lenght_-1){
                    // Correct nearest_list last element
                    nearest_list[1] = 0;
                }

                // Push nearest_list into the nearest_neighbors_ vector
                nearest_neighbors_.push_back(nearest_list);
            }
            nearest_list.clear();
        } else if (geometry_flag == 2) {
            tot_lenght_ = side_lenght * side_lenght;
            latt_conf.reserve(tot_lenght_);
            nearest_neighbors_.reserve(tot_lenght_);
            nearest_list.reserve(4);

            // Fill nearest_neighbors_ with nearest_list for each site
            for (int i = 0; i < tot_lenght_; i++){
                // Fill nearest_list for site i
                nearest_list.clear();
                nearest_list.push_back(i - 1);
                nearest_list.push_back(i + 1);
                nearest_list.push_back(i - side_lenght);
                nearest_list.push_back(i + side_lenght);

                // Implementing PBC
                if (i % side_lenght  == 0){
                    // Correct nearest_list first culomn elements
                    nearest_list[0] += side_lenght;
                } else if ((i + 1) % side_lenght  == 0){
                    // Correct nearest_list last culomn elements
                    nearest_list[1] -= side_lenght;
                }
                if (i < side_lenght){
                    // Correct nearest_list first row elements
                    nearest_list[2] += tot_lenght_;
                } else if(i >= tot_lenght_ - side_lenght){
                    // Correct nearest_list last row elements
                    nearest_list[3] -= tot_lenght_;
                }

                // Push nearest_list into the nearest_neighbors_ vector
                nearest_neighbors_.push_back(nearest_list);
            }
            nearest_list.clear();
        } else {
            cerr << "Error: geometry not implemented." << endl;
            exit(1);
        }

        // Initialising the lattice configuration
        if (initial_flag == 0) {
            for(int i = 0; i < tot_lenght_; i++)
                latt_conf.push_back(0.);
        } else {
            for (int i = 0; i < tot_lenght_; i++){
                if (rand_double() < 0.5){
                    latt_conf.push_back(-1.);
                } else {
                    latt_conf.push_back(1.);
                }
            }
        }
    }// CONSTRUCTION END

    int rand_int(){
        /* Generate a random index for the lattice */

        uniform_int_distribution<long int> random_int(0, tot_lenght_ - 1);
        return random_int(generator_);
    }

    double rand_double(){
        /* Generate a random double */

        uniform_real_distribution<double> random_num(0., 1.);
        return random_num(generator_);
    }

    void update(double eta){
        /* Update the state of the lattice with a MC step */

        int ind;
        double delta, force;

        // Metropolis-Hastings algorithm
        for(int it = 0; it < tot_lenght_; it++){
           ind = rand_int();
           delta = 1.5 * sqrt(eta) * (1. - 2 * rand_double());

           // MC step for the site with index ind
           force = - (delta + 2 * latt_conf[ind]) * (1. + pow(eta,2) / 2.);
           for(auto nn : nearest_neighbors_[ind]) force += latt_conf[nn];
           force = exp(force * (delta / eta));

           // accept or reject step
           if (rand_double() < force) latt_conf[ind] += delta;
        }
    }

    double energy(const double& eta){
        /* Compute the energy of the present configuration */

        double kin_ener = 0., pot_ener = 0.;

        for(int i = 0; i < tot_lenght_; i++) {
            if(i == tot_lenght_ - 1){
                kin_ener += pow(latt_conf[0] - latt_conf[i], 2);
            } else {
                kin_ener += pow(latt_conf[i + 1] - latt_conf[i], 2);
            }
            pot_ener += pow(latt_conf[i], 2);
        }
        kin_ener = - kin_ener / (tot_lenght_ * pow(eta, 2));
        pot_ener = pot_ener / tot_lenght_;

        return (kin_ener + pot_ener + 1./eta) / 2;
    }

    double correlator(int k){
        /* Compute the correlator function between i and i+k */

        int index = 0;
        double corr = 0.;

        for(int i = 0; i < tot_lenght_; i++) {
            index = i + k;
            if(index > tot_lenght_ - 1)
                index -= tot_lenght_ - 1;
            corr += latt_conf[index] * latt_conf[i];
        }
        corr = corr / tot_lenght_;

        return corr;

    }

    double correl_sqr(int k){
        /* Compute the correlator function between i and i+k */

        int index = 0;
        double corr = 0.;

        for(int i = 0; i < tot_lenght_; i++) {
            index = i + k;
            if(index > tot_lenght_ - 1)
                index -= tot_lenght_ - 1;
            corr += pow(latt_conf[index] * latt_conf[i], 2);
        }
        corr = corr / tot_lenght_;

        return corr;

    }

    //--- Config methods -------------------------------------------------------

    void save_configuration(string file_state = "test_traje"){
        /* Save the current configuration of lattice and generator */

        ofstream file;
        file.open(file_state + "_state.dat");
        for (int i = 0; i < tot_lenght_; i++){
            file << latt_conf[i] << endl;
        }
        file.close();

        ofstream file_seed;
        file_seed.open(file_state + "_seed.dat");
        file_seed << generator_;
        file_seed.close();

    }

    void load_configuration(string file_state = "test_traje"){
        /* Load configuration and seed from file */

        ifstream file_seed(file_state + "_seed.dat");
        cout << "Loading initial seed..." << endl;
        if (file_seed.is_open()) {
            file_seed >> generator_;
            file_seed.close();
        } else {
            cerr << "Error: unable to open the seed file." << endl;
            exit(1);
        }

        int i = 0;
        ifstream file(file_state + "_state.dat");
        cout << "Loading initial configuration..." << endl;
        if (file.is_open()) {
            while((file >> latt_conf[i]) && (i < tot_lenght_)) i++;
            file.close();
        } else {
            cerr << "Error: unable to open the file." << endl;
            exit(1);
        }

        if (i != tot_lenght_){
            cerr << "Error: different number of sites." << endl;
            exit(1);
        }

    }

    //--- Show methods ---------------------------------------------------------

    void show_configuration(){
        /* Print the lattice configuration */

        int value;
        cout << "Lattice configuration: " << endl;
        if (geometry_flag == 1){
            for (int i = 0; i < tot_lenght_; i++){
                cout << latt_conf[i] << " ";
            }
            cout << endl << endl;
        } else if (geometry_flag == 2){
            for (int i = 0; i < tot_lenght_; i++){
                cout << latt_conf[i] << " ";
                if ((i + 1) % side_lenght == 0){
                    cout << endl;
                }
            }
            cout << endl;
        }
    }

    void show_nearest_index(const int &index){
        /* Print the nearest neighbors of the site index */

        int size;
        size = nearest_neighbors_[index].size();
        for (int i = 0; i < size; i++){
            cout << nearest_neighbors_[index][i] << " ";
        }
    }

    void show_nearest_neighbors(){
        /* Print the list of all nearest neighbors to each lattice site */

        for (int i = 0; i < tot_lenght_; i++){
            cout << "Nearest neighbors to " << i << " are: ";
            show_nearest_index(i);
            cout << endl;
        }
        cout << endl;
    }

};

#endif
