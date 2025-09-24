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
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){
   // Random number generator initialization
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes.txt");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   int M = 10000;
   double lambda = 1.0; // Exponential distribution scale factor
   double mean = 0.0, gamma = 1.0; // Cauchy-Lorentz distribution parameters

   vector<double> N = {1, 2, 10, 100}; // Different sample sizes for the sums (number of random variables to average)

   // Output file initialization (one for each n)
   for (int n : N) {
      ofstream out("N_" + std::to_string(n) + ".dat");
      if (!out) {
         cerr << "PROBLEM: unable to open file for N = " << n << endl;
         return 1;
       }
      // Sums and mean value calculation
       for (int i=0; i<M; i++){
         double sum_unif = 0;
         double sum_exp = 0;   
         double sum_lor = 0; 
      
         for(int j=0; j<n; j++){
            sum_unif += rnd.Rannyu()/n;
            sum_exp+=rnd.Exponent(lambda)/n;
            sum_lor+=rnd.CauchyLorentz(mean, gamma)/n; 
         }
            out << sum_unif << " " << sum_exp << " " << sum_lor << endl;
        }
      }
   
   cout << "Data files written succesfully!" << endl;

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
