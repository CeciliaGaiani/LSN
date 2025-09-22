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

/****************************************************************
                      FUNCTION DECLARATIONS
*****************************************************************/

// Formula for the statistical error (standard deviation of the mean) calculated using the data blocking method
double error(double ave, double ave2, int n){
   if (n == 0){
      return 0;
   } else{
      return sqrt((ave2-pow(ave,2))/n);
   }
 }

// Formula for the calculation of the value of 𝝌²: counts contains the calculated. Here counts[i] = observed value in each bin; exp = expected value for each bin (constant in this case)
 double chi_square(vector<double> counts, double exp){
   double sum = 0;
   for(int i=0; i<counts.size(); i++){
      sum+=pow(counts[i]-exp, 2);
   }
   return sum/exp;
 }

/****************************************************************
                               MAIN
*****************************************************************/

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

   // Output files initialization
   ofstream fout("measures.dat"), fout_chi("chi_square.dat");

   if (!fout || !fout_chi) {
      cerr << "Unable to open file measures.dat or chi_square.dat!" << endl;
      return 1;
  }

// MEAN VALUE AND VARIANCE (WITH CORRESPECTRIVE UNCERTAINTIES) CALCULATION

   int M = 10000; // Number of throws
   int N = 100;   // Number of blocks
   int L = M/N;   // Number of throws per block

   // Initialization of STL vectors to be used in the calculation of mean values and uncertanties
   vector<double> ave, ave2;
   vector<double> var, var2;
   vector<double> sum_prog,sum2_prog;
   vector<double> var_prog, var2_prog;
   vector<double> err_prog, err_var;
   
   for (int i=0; i<N; i++){    // Cycle on the number of blocks      
     double sum = 0; 
     double sum_var = 0;            
      for (int j=0; j<L; j++){ // Cycle on the number of throws per block
         double n = rnd.Rannyu(); 
         sum += n;
         sum_var += pow(n-0.5, 2);
      }
      double mean = sum / L;
      double var_block = sum_var / L;

      ave.push_back(mean);
      ave2.push_back(mean * mean); // We must store the square of the block mean, and then calculate the mean squared average over the blocks
      var.push_back(var_block);
      var2.push_back(var_block * var_block);
    }

   for (int i=0; i<N; i++){
      double sum = 0;
      double sigma = 0;
      double sum2 = 0;
      double sigma2 = 0;

      for (int j=0; j<i+1; j++){ // Perform progressive sums
        sum += ave[j];
        sigma += var[j];
        sum2 += ave2[j];
        sigma2 += var2[j];
       }
 
       sum_prog.push_back(sum/(i+1));
       sum2_prog.push_back(sum2/(i+1));
       var_prog.push_back((sigma)/(i+1));
       var2_prog.push_back((sigma2)/(i+1));
       err_prog.push_back(error(sum_prog[i], sum2_prog[i], i));
       err_var.push_back(error(var_prog[i], var2_prog[i], i));
    }

    for (int i=0; i<N; i++){
      fout << (i+1) << " " << sum_prog[i] << " " << " " << err_prog[i] << " " << " " << var_prog[i] << " " << err_var[i] << endl;
    }
   
   fout.close();

// 𝝌² CALCULATION (FOR A UNIFORM DISTRIBUTION)

   double N_throws = 10000;                  // Number of throws for 𝝌² calculation
   double N_bins = 100;                      // Number fo bins
   double N_att = 500;                       // Number of 𝝌² tests to be performed                              
   double exp_counts = N_throws / N_bins;    // Number of expected counts per bin
   double delta_int = 1.0 / N_bins;          // Bin amplitude

   for (int i = 0; i < N_att; i++) {
      vector<double> counts(N_bins, 0.);
      vector<double> extractions;

      for (int j = 0; j < N_throws; j++) {  // Generate N_throws random numbers uniformly distributed in [0,1) and store them in extractions
         extractions.push_back(rnd.Rannyu());
      }

      for (int i = 0; i < N_throws; i++) {
        for (int n_bin = 0; n_bin < N_bins; n_bin++) { // Assign each extraction to the correct bin by checking which sub-interval it falls into, then increment the count for that bin 
            if (extractions[i] < (n_bin + 1) * delta_int && extractions[i] >= n_bin * delta_int) {
                counts[n_bin]++;
                break;
            }
        }
    }

    double chi = chi_square(counts, exp_counts);
    fout_chi << chi << endl;
   }

   fout_chi.close();

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