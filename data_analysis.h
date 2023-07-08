/*******************************************************************************
*
* Data analysis program for the outcomes of the simulations
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#ifndef DATA_ANAL_H
#define DATA_ANAL_H

#include <iostream>
#include <fstream>
#include <chrono>

#include <cmath>
#include <vector>
#include <string>
#include <random>

using namespace std;

// Define the PRNG
#define SEED 42
mt19937_64 generator(SEED);

/*******************************************************************************
* PARAMETERS OF THE ANALYSIS
*
* MIN_CORR_LENGHT = min dim of the correlated blocks in the bootstrap algorithm.
*
* MAX_CORR_LENGHT = max dim of the correlated blocks in the bootstrap algorithm.
*
* NUM_FAKE_SAMP = number of fake samples in the bootstrap algorithm.
*
* DIM_FAKE_SAMP = dimension of the fake samples in the bootstrap algorithm.
*
*******************************************************************************/

#define MIN_CORR_LENGHT 1
#define MAX_CORR_LENGHT 10
#define NUM_FAKE_SAMP 30
#define DIM_FAKE_SAMP 1000

//--- Contents -----------------------------------------------------------------

double bootstrap_corr(vector<double>& x, ofstream &file){
    /* Compute sigma with the bootstrap algorithm */
    // Given a sample x, generate NUM_FAKE_SAMP samples via resampling
    // and compute average and variance of an estimator while varying the
    // correlation lenght from MIN to MAX_CORR_LENGHT and writing them
    // into a file. The function bootstrap returns the last std deviation.

    int rand_index;
    double est_val, est_ave, est_var;
    vector<double> fake_sample, measures;

    measures.reserve(NUM_FAKE_SAMP);
    fake_sample.reserve(DIM_FAKE_SAMP);
    uniform_int_distribution<long int> randomint(0, x.size());

    // iterate over different correlation lenghts
    for(int lenght = MIN_CORR_LENGHT; lenght <= MAX_CORR_LENGHT; lenght *= 2){

        // initialite the average of the
        // estimator over the fake samples
        est_ave = 0.;

        // evaluate the estimator over each fake sample
        // and save the result into the vector measures
        for(int i = 0; i < NUM_FAKE_SAMP; i++){

            // generate the i-fake_sample
            for (int j = 0; fake_sample.size() < DIM_FAKE_SAMP; j++){
                // extract the j-random_int and push the j-block
                rand_index = randomint(generator);
                for (int k = 0; k < lenght; k++){
                    if (rand_index + k == x.size()) break;
                    if (fake_sample.size() == DIM_FAKE_SAMP) break;
                    fake_sample.push_back(x[rand_index + k]);
                }
            }

            // evaluate the estimator and save the value
            est_val = 0;
            for(auto val: fake_sample) est_val += val;
            est_val = est_val / DIM_FAKE_SAMP;

            est_ave += est_val;
            measures.push_back(est_val);
            fake_sample.clear();
        }
        est_ave = est_ave / NUM_FAKE_SAMP;

        // compute the estimator variance over fake samples
        est_var = 0.;
        for (auto val : measures) est_var += pow(val - est_ave, 2);
        est_var = est_var / (NUM_FAKE_SAMP - 1);

        measures.clear();

        file << "correl lenght " << lenght << " -> ";
        file << est_ave << " ± " << sqrt(est_var) << endl;
    }

    return sqrt(est_var);
}

//--- File operations ----------------------------------------------------------

vector<double> file_operations(int side, int beta, ofstream &file_analysis){
    /* Operations over each beta file */
    // For a given side and beta, compute average energy and errors.

    int measures = 0;
    double ene, ene_ave = 0.;
    vector<double> energies, data;
    string file_name;

    file_name = "Data_simulations/Beta_" + to_string(beta) + "/";
    file_name += "side_" + to_string(side) + ".dat";
    ifstream file(file_name);

    if (file.is_open()) {
        // load data and compute averages
        while (file >> ene){
            measures++;
            ene_ave += ene;
            energies.push_back(ene);
        }
        file.close();
    } else {
        cout << "Error: unable to open the file." << endl;
        exit(1);
    }

    file_analysis << "-----------------------------" << endl;
    file_analysis << "side: " << side << " | beta: " << beta << endl;
    file_analysis << "-----------------------------" << endl;

    // compute and store the averages
    ene_ave = ene_ave / measures;
    data.push_back(ene_ave);

    // compute the errors with bootstrap algorithm
    file_analysis << " --- energy error: " << endl;
    data.push_back(bootstrap_corr(energies, file_analysis));

    file_analysis << "-----------------------------" << endl << endl;
    return data;
}

#endif
