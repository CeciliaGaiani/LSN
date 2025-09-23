/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "random.h"

using namespace std;

/****************************************************************
                      FUNCTION DECLARATIONS
*****************************************************************/

// Trial wave function (variational parameters = μ, 𝜎)
double psi_T(double x, double sigma, double mu) {
    double sigma2 = sigma*sigma;
    double exp_minus = exp(-((x-mu)*(x-mu))/(2.0*sigma2));
    double exp_plus = exp(-((x+mu)*(x+mu))/(2.0*sigma2));
    return exp_minus + exp_plus;
}

// Squared modulus of the trial wave function
double squared_psi_T(double x, double sigma, double mu) {
    double psi = psi_T(x, sigma, mu);
    return psi * psi;
}

// Second derivative of the trial wave function
double d2_psi_T(double x, double sigma, double mu) {
    double sigma2 = sigma*sigma;
    double sigma4 = sigma2*sigma2;

    double x_m = x-mu;
    double x_p = x+mu;

    double exp_m = exp(-x_m * x_m/(2.0 * sigma2));
    double exp_p = exp(-x_p * x_p/(2.0 * sigma2));

    double term_m = ((x_m*x_m)-sigma2)/sigma4*exp_m;
    double term_p = ((x_p*x_p)-sigma2)/sigma4*exp_p;

    return term_m + term_p;
}

// Hamiltonian expectation value (local energy)
double exp_H(double x, double sigma, double mu) {
    double psi = psi_T(x, sigma, mu);
    double d2psi = d2_psi_T(x, sigma, mu);
    double V = pow(x, 4) - 2.5 * x * x;

    return -0.5*(d2psi/psi) + V;
}

// Hamiltonian expectation value computed via Metropolis algorithm (single value)
double Metropolis_ave_H(double x_old, Random& rnd, int M, int N, double step, double sigma, double mu) {
    double psi2_old = squared_psi_T(x_old, sigma, mu);
    double H_old = exp_H(x_old, sigma, mu);

    int L = M / N; 
    double global_sum = 0.0;

    for (int i = 0; i < N; ++i) {
        double block_sum = 0.0;

        for (int j = 0; j < L; ++j) {
            double x_new = x_old + rnd.Rannyu(-step, step);
            double psi2_new = squared_psi_T(x_new, sigma, mu);
            double A = min(1.0, psi2_new / psi2_old);

            if (rnd.Rannyu() <= A) {
                x_old = x_new;
                psi2_old = psi2_new;
                H_old = exp_H(x_old, sigma, mu);
            }

            block_sum += H_old;
        }
        global_sum += block_sum / L;
    }
    return global_sum / N;
}

// Hamiltonian expectation value computed via Metropolis algorithm (fill a vector to use for data blocking)
vector<double> ave_H_prog(double x_old, Random& rnd, int M, int N, double step, double sigma, double mu, bool save_conf, ofstream & x_sampl) {
    
    int L = M / N;
    vector<double> H_mean(N, 0.0);

    double psi2_old = squared_psi_T(x_old, sigma, mu);
    double H_old = exp_H(x_old, sigma, mu);

    for (int i = 0; i < N; i++) {
        double sum = 0.0;

        for (int j = 0; j < L; j++) {
            double x_new = x_old + rnd.Rannyu(-step, step);
            double psi2_new = squared_psi_T(x_new, sigma, mu);
            double A = min(1.0, psi2_new / psi2_old);

            if (rnd.Rannyu() <= A) {
                x_old = x_new;
                psi2_old = psi2_new;
                H_old = exp_H(x_old, sigma, mu);
            }
            if (x_sampl && save_conf) x_sampl << x_old << endl;
            sum += H_old;
        }
        H_mean[i] = sum / L;
    }
    return H_mean;
}

// Compute progressive averages
vector<double> progressive_mean(vector<double>& H_block) {
    vector<double> mean(H_block.size(), 0.0);
    double sum = 0.0;
    for (size_t i = 0; i < H_block.size(); i++) {
        sum += H_block[i];
        mean[i] = sum / (i + 1);
    }
    return mean;
}

// Formula for the statistical error (standard deviation of the mean)
vector<double> progressive_error(vector<double>& mean, vector<double>& values) {
    vector<double> err(mean.size(), 0.0);
    double sum2 = 0.0;
    for (size_t i = 0; i < values.size(); i++) {
        sum2 += values[i] * values[i];
        double mean2 = sum2 / (i + 1);
        if (i == 0) err[i] = 0.0;
        else err[i] = sqrt((mean2 - mean[i] * mean[i]) / i);
    }
    return err;
}

/****************************************************************
                               MAIN
*****************************************************************/

int main() {
    
    // Random number generator initialization
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes.txt");
    if (Primes.is_open()) {     
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
        return 1;
    }
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        return 1;
    }

    // Output files initialization
    ofstream H_vs_T("H_vs_temperature.dat"), mu_sigma_acc("accepted_mu_sigma.dat"), best_E("E_optimal_values.dat"), acc_Metro_var("acceptance_Metro_var.dat"), x_sampl("sampled_configurations_x.dat");
    if (!H_vs_T || !mu_sigma_acc || !best_E || !acc_Metro_var || !x_sampl) {
        cerr << "PROBLEM: Unable to open output data files!" << endl;
        return 1;
    }

    int N_temp = 500;                           // Number of temperatures
    double T_old = 1, T_new, alpha = 0.99;      // Cooling
    double T_best;                              // Best values for T

    double mu_old = 1, mu_new;                // Starting value for μ
    double sigma_old = 1, sigma_new;          // Starting value for 𝜎
    mu_sigma_acc << mu_old << " " << sigma_old << endl;

    double mu_best = mu_old;                   
    double sigma_best = sigma_old;            // Best values for μ and 𝜎

    double H_min = 1e10;                      // GS energy  
   
    double delta_mu = 0.2;                    // μ variation
    double delta_sigma = 0.2;                 // 𝜎 variation
    
    double step = 0.5;                        // Metropolis step for <H>  
    double x_old = rnd.Rannyu(-step, step);   // Initial configuration                                           
    int N_samples = 15000;                    // Number of Metropolis attempts for position sampling
    int N_Metro = 200;                        // Number of Metropolis attempts for μ and 𝜎
    int N_blocks = 100;                       // Number of blocks for data blocking

    // Compute the local energy for the initial configuration
    double H_old = Metropolis_ave_H(x_old, rnd, N_samples, N_blocks, step, sigma_old, mu_old);
    for (int i_temp = 0; i_temp < N_temp; i_temp++) { // Cycle on temperatures
        T_new = alpha*T_old; // Exponential cooling
        acc_Metro_var << endl << "==================== Temperature n° " << i_temp+1 << ": " << T_new <<  " ====================" << endl;

        double accept_Metro = 0.0;
        double reject_Metro = 0.0;

        for (int i_Metro = 0; i_Metro<N_Metro; i_Metro++){ 

            // Variation of μ and 𝜎
            mu_new = mu_old + rnd.Rannyu(-delta_mu, delta_mu);
            sigma_new = sigma_old + rnd.Rannyu(-delta_sigma, delta_sigma);
            if(sigma_new < 0.05) sigma_new = 0.05;

            // Local energy for the new configuration
            double H_new = Metropolis_ave_H(x_old, rnd, N_samples, N_blocks, step, sigma_new, mu_new);
            
            // The proposal (μ, 𝜎) is accepted with Metropolis probability
            double dE = H_new - H_old;
            if (dE < 0 || rnd.Rannyu() < exp(-dE / T_new)) {
                mu_old = mu_new;
                sigma_old = sigma_new;
                H_old = H_new;
                // Update variables if the energy is lower
                if (H_new < H_min) {
                    H_min = H_new;
                    mu_best = mu_new;
                    sigma_best = sigma_new;
                    T_best = T_new;
                }
                accept_Metro++;
            } else {
                reject_Metro++;
            }
        }
        
        acc_Metro_var << "Acceptance for sigma/mu sampling: " << accept_Metro/(reject_Metro+accept_Metro)*100 << "%" << endl;
        acc_Metro_var << "Last accepted mu,sigma: " << mu_old << ", " << sigma_old << endl;
        
        // For each temperature, data blocking with the last μ and 𝜎
        vector<double> blocks = ave_H_prog(x_old, rnd, N_samples, N_blocks, step, sigma_old, mu_old, 0, x_sampl);
        vector<double> mean = progressive_mean(blocks);
        vector<double> err  = progressive_error(mean, blocks);
        
        // Write only last mean/err for H_vs_T
        H_vs_T << T_new << " " << mean.back() << " " << err.back() << "\n";

        mu_sigma_acc << mu_old << " " << sigma_old << endl;
               
        T_old = T_new; // Start with the new temperature
    }
    
    // Print the best values for μ and 𝜎
    cout << "\nBest found configuration:" << endl;
    cout << "mu = " << mu_best << ", sigma = " << sigma_best << ", E = " << H_min << ", T = " << T_best << endl;
    
    // Compute Hamiltonian expectation value for the best values of μ and 𝜎 (progressive averages)
    vector<double> blocks = ave_H_prog(x_old, rnd, 40000, 100, 0.7, sigma_best, mu_best, 1, x_sampl);
    vector<double> mean = progressive_mean(blocks);
    vector<double> err  = progressive_error(mean, blocks);
    for (size_t i = 0; i < mean.size(); ++i) {
        best_E << (i+1) << " " << mean[i] << " " << err[i] << "\n";
    }

    cout << endl << "Execution completed!" << endl;

    H_vs_T.close();
    mu_sigma_acc.close();
    best_E.close();
    acc_Metro_var.close();

    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
