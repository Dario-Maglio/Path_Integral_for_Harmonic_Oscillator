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
#include <sstream>
#include <string>
#include <chrono>

#include <cmath>
#include <vector>
#include <random>

using namespace std;

// Define the PRNG
#define SEED 42
mt19937_64 generator(SEED);

/*******************************************************************************
* PARAMETERS OF THE ANALYSIS
*
* CORREL_LENGTH = maximum value of k in correlators C_k.
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

// maximum correlation in a path
#define CORREL_LENGTH 40
// correlations in the MC-average
#define MIN_CORR_LENGHT 1
#define MAX_CORR_LENGHT 128
#define NUM_FAKE_SAMP 100
#define DIM_FAKE_SAMP 15000

//--- Contents -----------------------------------------------------------------

double bootstrap_corr(vector<double>& x, ofstream &file){
    /* Compute error with the bootstrap algorithm */
    // Given a sample x, generate NUM_FAKE_SAMP samples via resampling
    // and compute average and variance while varying the correlation
    // lenght from MIN to MAX_CORR_LENGHT and writing them into a file.
    // The function bootstrap returns the last std deviation.

    int rand_index;
    double est_val, est_ave, est_var;
    vector<double> fake_sample, measures;

    measures.reserve(NUM_FAKE_SAMP);
    fake_sample.reserve(DIM_FAKE_SAMP);
    uniform_int_distribution<int> randomint(0, x.size());

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
    /* Compute MC-average energy and error for a given beta and side */

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
    file_analysis << "result: " << data[0] << " ± " << data[1] << endl;

    file_analysis << "-----------------------------" << endl << endl;
    return data;
}

vector<double> gap_operations(int side, int beta, int label, ofstream &file_analysis){
    /* Compute MC-averages of C_k and errors for a given beta and side */

    int measures = 0;
    double value = 0., C_zero = 0.;
    vector<double> correl(2*CORREL_LENGTH, 0.);
    vector<double> values(CORREL_LENGTH, 0.);
    vector<vector<double>> data;
    string file_path, file_name, line;

    // define file path and name
    file_path = "Data_simulations/Beta_" + to_string(beta) + "/";
    file_name = "correlations_side_" + to_string(side) + "_";
    file_name += to_string(label) + ".dat";
    ifstream file(file_path + file_name);
    // load data and compute averages
    if (file.is_open()) {
        // get a line with k elements C_k
        while(getline(file, line)){
            measures++;
            istringstream sstream(line);
            sstream >> value;
            C_zero += value;
            // sum each C_k in correl[k]
            for(int k = 0; k < CORREL_LENGTH; k++){
                sstream >> value;
                correl[2*k] += value;
                values[k] = value;
            }
            data.push_back(values);
        }
        file.close();
    } else {
        // Stop execution: invalid file
        cerr << "Error: unable to open the file." << endl;
        exit(1);
    }

    // normalize each C_k and compute the error
    C_zero = C_zero / measures;
    values.resize(measures, 0.);
    for(int k = 0; k < CORREL_LENGTH; k++){

        // normalize the average
        correl[2*k] = (correl[2*k] / measures) - pow(C_zero, 2);

        // load measures for this k
        for(int i = 0; i < measures; i++)
            values[i] = data[i][k];

        // call the bootstrap algorithm
        file_analysis << "-----------------------------" << endl;
        file_analysis << "k_val: " << k + 1 << " | side: " << side << endl;
        file_analysis << "-----------------------------" << endl;
        correl[2*k+1] = bootstrap_corr(values, file_analysis);
        file_analysis << "result: " << correl[2*k] << " ± " << correl[2*k+1];
        file_analysis << endl << endl;
    }

    file_analysis << "-----------------------------" << endl << endl;
    return correl;
}

#endif
