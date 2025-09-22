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
#include <cmath>
#include <algorithm>
#include "random.h"

using namespace std;

/****************************************************************
                      FUNCTION DECLARATIONS
*****************************************************************/

// Function computing the price value by sampling directly the final asset price S(T)
double S_dir(double S0, double T, double r, double sigma, double Z){
   return S0 * exp((r - 0.5 * sigma * sigma) * T + sigma * Z * sqrt(T));
}

// Function computing the discounted payoff of a call option at maturity given price S
double call(double r, double T, double S, double K){
    return exp(-r * T) * max(0.0, S - K);
}

// Function computing the discounted payoff of a put option at maturity given price S
double put(double r, double T, double S, double K){
    return exp(-r * T) * max(0.0, K - S);
}

// Formula for the statistical error (standard deviation of the mean) calculated using the data blocking method
double error(double ave, double ave2, int n){
   if (n<2) return 0;
   return sqrt((ave2 - ave * ave) / n);
}

// Function computing the price value by sampling the discretized GMB(r, 𝜎²) path of the asset price
double S_discr(double S0, double r, double sigma, double T, int steps, Random &rnd) {
    double dt = T / steps;
    double S = S0;

    for (int i = 0; i < steps; i++) {
        double Z = rnd.Gauss(0, 1);
        S *= exp((r - 0.5 * sigma * sigma) * dt + sigma * Z * sqrt(dt));
    }
    return S;
}

/****************************************************************
                               MAIN
*****************************************************************/

int main(int argc, char *argv[]){
   // Random number generator initialization
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes.txt");
   if (Primes.is_open()) {
      Primes >> p1 >> p2;
      Primes.close();
   } else {
      cerr << "PROBLEM: Unable to open Primes" << endl;
      return 1;
   }

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

   // Ouput files initialization
   ofstream foutC("Call.dat"), foutP("Put.dat");
   ofstream foutCi("Call_i.dat"), foutPi("Put_i.dat");

   if (!foutC || !foutP) {
      cerr << "Errore nell'apertura dei file di output!" << endl;
      return 1;
   }

   double S0 = 100;       // Price at time t=0
   double T = 1;          // Delivery time
   double K = 100;        // Strike price
   double r = 0.1;        // Risk free interest rate
   double sigma = 0.25;   // Volatility

   int M = 10000;         // Number of throws
   int N = 100;           // Number of blocks
   int L = M / N;         // Number of throws per block

   vector<double> C(N), C2(N), P(N), P2(N);
   vector<double> Ci(N), C2i(N), Pi(N), P2i(N);

   // Direct sampling of the asset price S(T)
   for (int i = 0; i < N; i++) {
      double c = 0, p = 0;
      for (int j = 0; j < L; j++) {
         double Z = rnd.Gauss(0, 1);
         double ST = S_dir(S0, T, r, sigma, Z);
         c += call(r, T, ST, K);
         p += put(r, T, ST, K);
      }
      C[i] = c/L;
      C2[i] = C[i]*C[i];
      P[i] = p/L;
      P2[i] = P[i]*P[i];
   }

   double sum_c = 0, sum_c2 = 0, sum_p = 0, sum_p2 = 0;
   for (int i = 0; i < N; i++) {
      sum_c += C[i]; 
      sum_c2 += C2[i];
      sum_p += P[i]; 
      sum_p2 += P2[i];
      foutC << i + 1 << " " << sum_c / (i + 1) << " " << error(sum_c / (i + 1), sum_c2 / (i + 1), i) << endl;
      foutP << i + 1 << " " << sum_p / (i + 1) << " " << error(sum_p / (i + 1), sum_p2 / (i + 1), i) << endl;
   }

   // Sampling the discretized GMB(r, 𝜎²) path of the asset price
   int steps = 100; // Divide [0, T] in 100 intervals
   for (int i = 0; i < N; i++) {
      double c = 0, p = 0;
      for (int j = 0; j < L; j++) {
         double ST = S_discr(S0, r, sigma, T, steps, rnd);
         c += call(r, T, ST, K);
         p += put(r, T, ST, K);
      }
      Ci[i] = c/L;
      C2i[i] = Ci[i]*Ci[i];
      Pi[i] = p/L;
      P2i[i] = Pi[i]*Pi[i];
   }

   double sum_ci = 0, sum_c2i = 0, sum_pi = 0, sum_p2i = 0;
   for (int i = 0; i < N; i++) {
      sum_ci += Ci[i]; 
      sum_c2i += C2i[i];
      sum_pi += Pi[i]; 
      sum_p2i += P2i[i];
      foutCi << i + 1 << " " << sum_ci / (i + 1) << " " << error(sum_ci / (i + 1), sum_c2i / (i + 1), i) << endl;
      foutPi << i + 1 << " " << sum_pi / (i + 1) << " " << error(sum_pi / (i + 1), sum_p2i / (i + 1), i) << endl;
   }

   foutC.close();
   foutP.close();
   foutCi.close();
   foutPi.close();

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