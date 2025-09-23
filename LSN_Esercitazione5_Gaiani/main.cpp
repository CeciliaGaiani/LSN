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
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "random.h"

using namespace std;

/****************************************************************
                      FUNCTION DECLARATIONS
*****************************************************************/

// Struct defining a position
struct pos {
    double x, y, z;
};

// Function computing distance from the origin
double r(const pos& pos) {
    return sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);
}

//Function computing polar angle in spherical coordinates
double theta(const pos& pos) {
    return pos.x/r(pos);
}

// Square modulus for 1s wave function
double psi1s(const pos& pos) {
   return 1./M_PI * exp(-2 * r(pos)); 
}

// Square modulus for 2p wave function
double psi2p(const pos& pos) {
    double R = r(pos);
    return 1.0 / (32.0 * M_PI) * pos.z * pos.z * exp(-R);
}

// Metropolis algorithm for uniform transition probability (for a single block)
double MetropolisUniform(Random& rnd, pos& old, int L, double delta,
                       double& accept, double& reject, ofstream& out, string orb_type) {
    double r_sum = 0;                    
    for (int j = 0; j < L; j++) {
        pos prop = {
            old.x + rnd.Rannyu(-delta, delta),
            old.y + rnd.Rannyu(-delta, delta),
            old.z + rnd.Rannyu(-delta, delta)
        };

        double p_old, p_new;

        if(orb_type == "1s"){
            p_old = psi1s(old);
            p_new = psi1s(prop);
        } else {
            p_old = psi2p(old);
            p_new = psi2p(prop);
        }

        double A = min(1.0, p_new/p_old);
        if (rnd.Rannyu() <= A) {
            old = prop;
            accept++;
        } else {
            reject++;
        }

        out << old.x << " " << old.y << " " << old.z << endl;
        r_sum += r(old);
    }
    return r_sum;
}

// Metropolis algorithm for Gaussian transition probability (for a single block) 
double MetropolisGaussian(Random& rnd, pos& old, int L, double sigma,
                        double& accept, double& reject, ofstream& out, string orb_type) {
    double r_sum = 0;
    for (int j = 0; j < L; j++) {
        pos prop = {
            old.x + rnd.Gauss(0., sigma),
            old.y + rnd.Gauss(0., sigma),
            old.z + rnd.Gauss(0., sigma),
        };

        double p_old, p_new;

        if(orb_type == "1s"){
            p_old = psi1s(old);
            p_new = psi1s(prop);
        } else {
            p_old = psi2p(old);
            p_new = psi2p(prop);
        }

        double A = min(1.0, p_new/p_old);
        if (rnd.Rannyu() <= A) {
            old = prop;
            accept++;
        } else {
            reject++;
        }

        out << old.x << " " << old.y << " " << old.z << endl;
        r_sum += r(old);
    }
    return r_sum;
}

pos equilibration(Random& rnd, pos& old, int L, double par, string orb_type, string distr) {
    double accept=0, reject=0;
    for (int j = 0; j < L; j++) {
        pos prop;
        if(distr == "unif"){
            prop = {
                old.x + rnd.Rannyu(-par, par),
                old.y + rnd.Rannyu(-par, par),
                old.z + rnd.Rannyu(-par, par)
            };
        } else {
            prop = {
                old.x + rnd.Gauss(0., par),
                old.y + rnd.Gauss(0., par),
                old.z + rnd.Gauss(0., par)
            };
        }

        double p_old, p_new;
        if(orb_type == "1s"){
            p_old = psi1s(old);
            p_new = psi1s(prop);
        } else {
            p_old = psi2p(old);
            p_new = psi2p(prop);
        }

        double A = min(1.0, p_new/p_old);
        if (rnd.Rannyu() <= A) {
            old = prop;
            accept++;
        }else{
            reject++;
        }
    }
    cout << "Acceptance: " << accept/(accept+reject) << endl;
    return old;
}


// Formula for the statistical error (standard deviation of the mean) calculated using the data blocking method
double error(double ave, double ave2, int n) {
    if (n<2) {
        return 0;
    } else {
        return sqrt((ave2 - pow(ave, 2)) / n);
    }
  }   

/****************************************************************
                               MAIN
*****************************************************************/

int main() {

    // Random number generator initialization
    Random rnd;
    int seed[4];
    int p1, p2;

    ifstream Primes("Primes.txt");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
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
    }

    // Output files initialization (1s orbital)
    ofstream fout_unif("radius_unif.dat"), fout_gauss("radius_gauss.dat"), points_1s_unif("points_1s_unif.dat"), points_1s_gauss("points_1s_gauss.dat");
    if (!fout_unif) {
        cerr << "PROBLEM: Unable to open radius_unif.dat" << endl;
        return 1;
    }

    if (!fout_gauss) {
        cerr << "PROBLEM: Unable to open radius_gauss.dat" << endl;
        return 1;
    }

    if (!points_1s_unif) {
        cerr << "PROBLEM: Unable to open points_1s_unif.dat" << endl;
        return 1;
    }

    if (!points_1s_gauss) {
        cerr << "PROBLEM: Unable to open points_1s_Gauss.dat" << endl;
        return 1;
    }

    // Output files initialization (2p orbital)
    ofstream fout_unif_2p("radius_unif_2p.dat"), fout_gauss_2p("radius_gauss_2p.dat"), points_2p_unif("points_2p_unif.dat"), points_2p_gauss("points_2p_gauss.dat");
    if (!fout_unif_2p) {
        cerr << "PROBLEM: Unable to open radius_unif_2p.dat" << endl;
        return 1;
    }

    if (!fout_gauss_2p) {
        cerr << "PROBLEM: Unable to open radius_gauss_2p.dat" << endl;
        return 1;
    }

    if (!points_2p_unif) {
        cerr << "PROBLEM: Unable to open points_2p_unif.dat" << endl;
        return 1;
    }

    if (!points_2p_gauss) {
        cerr << "PROBLEM: Unable to open points_2p_unif.dat" << endl;
        return 1;
    }


    int M = 1000000;   // Number of throws
    int N = 100;       // Number of blocks
    int L = M / N;     // Number of throws per block

    // Vector (for block means) initialization
    vector<double> r_mean, r2_mean, r_mean_cum, r2_mean_cum, r_err;
    vector<double> r_mean_gauss, r2_mean_gauss, r_mean_cum_gauss, r2_mean_cum_gauss, r_err_gauss;
    vector<double> r_mean_2p, r2_mean_2p, r_mean_cum_2p, r2_mean_cum_2p, r_err_2p;
    vector<double> r_mean_2p_gauss, r2_mean_2p_gauss, r_mean_cum_2p_gauss, r2_mean_cum_2p_gauss, r_err_2p_gauss;

    double accept_1s_unif = 0, accept_1s_Gauss = 0, accept_2p_unif = 0, accept_2p_Gauss = 0;
    double reject_1s_unif = 0, reject_1s_Gauss = 0, reject_2p_unif = 0, reject_2p_Gauss = 0;


    // Starting positions
    pos start_1s = {0, 0, 0};
    pos start_2p = {0, 0, 1};

    // Equilibration
    pos start_1s_unif = equilibration(rnd, start_1s, 5000, 0.75, "1s", "unif");
    pos start_1s_gauss = equilibration(rnd, start_1s, 5000, 0.75, "1s", "gauss");
    pos start_2p_unif = equilibration(rnd, start_2p, 5000, 2.8, "2p", "unif");
    pos start_2p_gauss = equilibration(rnd, start_2p, 5000, 1.7, "2p", "gauss");

    for (int i = 0; i < N; i++) {

        double r_Metro_unif_1s = MetropolisUniform(rnd, start_1s_unif, L, 0.75, accept_1s_unif, reject_1s_unif, points_1s_unif, "1s");
        double r_Metro_Gauss_1s = MetropolisGaussian(rnd, start_1s_gauss, L, 0.75, accept_1s_Gauss, reject_1s_Gauss, points_1s_gauss, "1s");
        
        double r_Metro_unif_2p = MetropolisUniform(rnd, start_2p_unif, L, 2.8, accept_2p_unif, reject_2p_unif, points_2p_unif, "2p");
        double r_Metro_Gauss_2p = MetropolisGaussian(rnd, start_2p_gauss, L, 1.7, accept_2p_Gauss, reject_2p_Gauss, points_2p_gauss, "2p");;

        r_mean.push_back(r_Metro_unif_1s/L);
        r2_mean.push_back(pow(r_mean[i], 2));

        r_mean_gauss.push_back(r_Metro_Gauss_1s/L);
        r2_mean_gauss.push_back(pow(r_mean_gauss[i], 2));

        r_mean_2p.push_back(r_Metro_unif_2p/L);
        r2_mean_2p.push_back(pow(r_mean_2p[i], 2));

        r_mean_2p_gauss.push_back(r_Metro_Gauss_2p/L);
        r2_mean_2p_gauss.push_back(pow(r_mean_2p_gauss[i], 2));
    }

    cout << "Acceptance (1s, unif): " << (accept_1s_unif / (accept_1s_unif+reject_1s_unif)) * 100 << "%" << endl;
    cout << "Acceptance (1s, Gauss): " << (accept_1s_Gauss / (accept_1s_Gauss+reject_1s_Gauss)) * 100 << "%" << endl;
    cout << "Acceptance (2p, unif): " << (accept_2p_unif / (accept_2p_unif+reject_2p_unif )) * 100 << "%" << endl;
    cout << "Acceptance (2p, Gauss): " << (accept_2p_Gauss / (accept_2p_Gauss+reject_2p_Gauss)) * 100 << "%" << endl;

    for (int i=0; i<N; i++){
        double r_prog_1s_unif=0, r_prog_1s_Gauss=0, r_prog_2p_unif=0, r_prog_2p_Gauss=0;
        double r2_prog_1s_unif=0,  r2_prog_1s_Gauss=0, r2_prog_2p_unif=0, r2_prog_2p_Gauss=0;
     
        for (int j=0; j<i+1; j++){
           r_prog_1s_unif += r_mean[j];
           r2_prog_1s_unif += r2_mean[j];

           r_prog_1s_Gauss += r_mean_gauss[j];
           r2_prog_1s_Gauss += r2_mean_gauss[j];

           r_prog_2p_unif += r_mean_2p[j];
           r2_prog_2p_unif += r2_mean_2p[j];

           r_prog_2p_Gauss += r_mean_2p_gauss[j];
           r2_prog_2p_Gauss += r2_mean_2p_gauss[j];
        }

        r_mean_cum.push_back(r_prog_1s_unif/(i+1));
        r2_mean_cum.push_back(r2_prog_1s_unif/(i+1));
        r_err.push_back(error(r_mean_cum[i], r2_mean_cum[i], i));

        r_mean_cum_gauss.push_back(r_prog_1s_Gauss/(i+1));
        r2_mean_cum_gauss.push_back(r2_prog_1s_Gauss/(i+1));
        r_err_gauss.push_back(error(r_mean_cum_gauss[i], r2_mean_cum_gauss[i], i));

        r_mean_cum_2p.push_back(r_prog_2p_unif/(i+1));
        r2_mean_cum_2p.push_back(r2_prog_2p_unif/(i+1));
        r_err_2p.push_back(error(r_mean_cum_2p[i], r2_mean_cum_2p[i], i));

        r_mean_cum_2p_gauss.push_back(r_prog_2p_Gauss/(i+1));
        r2_mean_cum_2p_gauss.push_back(r2_prog_2p_Gauss/(i+1));
        r_err_2p_gauss.push_back(error(r_mean_cum_2p_gauss[i], r2_mean_cum_2p_gauss[i], i));
     
     }
        
    
    for (int i = 0; i < N; i++) {
        fout_unif << i + 1 << " " << r_mean_cum[i] << " " << r_err [i] << endl;
        fout_gauss << i + 1 << " " << r_mean_cum_gauss[i] << " " << r_err_gauss[i] << endl;
        fout_unif_2p << i + 1 << " " << r_mean_cum_2p[i] << " " << r_err_2p[i] << endl;
        fout_gauss_2p << i + 1 << " " << r_mean_cum_2p_gauss[i] << " " << r_err_2p_gauss[i] << endl;
     }

    cout << "Data files produced!" << endl;

    fout_unif.close();
    fout_gauss.close();
    fout_unif_2p.close();
    fout_gauss_2p.close();

    rnd.SaveSeed();
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