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

// Function computing 𝜋
 double compute_pi(double N_hit, double N_tot, double L, double d){
   if (N_hit == 0) {
      cerr << "ERROR: N_hit is zero, division by zero avoided." << endl;
      return 0;
   }
   return (2*L*N_tot)/(d*N_hit);
}

// Funtion computing r in polar coordinates
 double r (double x, double y){
   return sqrt(x*x + y*y);
 }

// Funtion computing 𝜃 in polar coordinates
 double theta (double x, double y){
   return acos(x/r(x,y));
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

   // Output file initialization
   ofstream fout("pi.dat");

   if (!fout) {
      cerr << "Unable to open file pi.dat!" << endl;
      return 1;
  }

   int N_tot=500000;             // Number of throws
   int N=100;                    // Number of blocks
   int L=N_tot/N;                // Number of throws per block

   double d = 5;   // Distance between two adjacent lines
   double l = 4; // Needle length

   // Buffon's experiment requires d > L
    if(d<l){
      cout << "ERROR: d must be grater than L!" << endl;
      return -1;
    }

   // Initialization of STL vectors to be used in the calculation of mean values and uncertanties   vector<double> pi, pi2;
   vector<double> pi, pi2;
   vector<double> pi_prog, pi2_prog;
   vector<double> err_pi;

   for (int i = 0; i < N; i++) { // Cicle on the number of blocks
      double N_hit = 0;
      for (int j = 0; j < L; j++) { // Cycle on the number of throws per block
          double x, y;
          double pos = rnd.Rannyu(0, d/2.); // Extract the distance of the center of the needle from the closest line
  
          do {
              x = rnd.Rannyu(-1, 1);
              y = rnd.Rannyu(-1, 1);
          } while (r(x, y) > 1);
  
          double angle = acos(x/r(x, y)); // Extraction of the (random) angle that the needle forms with the x-axis
          double proj = (l/2) * sin(angle); // Projection of the half-legth of the needle along the y-axis
  
          if (proj>=pos) { // Check if the needle intersects the line
              N_hit++;
          }
      }
  
      double pi_estimate = compute_pi(N_hit, L, l, d);
      pi.push_back(pi_estimate);
      pi2.push_back(pow(pi_estimate, 2));
   }
   
   //Progressive averages
   for (int i=0; i<N; i++){
      double sum = 0;
      double sum2 = 0;

      for (int j=0; j<i+1; j++){
        sum +=pi[j];
        sum2 += pi2[j];
       }
 
       pi_prog.push_back(sum/(i+1));
       pi2_prog.push_back(sum2/(i+1));
       err_pi.push_back(error(pi_prog[i], pi2_prog[i], i));
   }

    for (int i=0; i<N; i++){
      fout << (i+1) << " " << pi_prog[i] << " " << " " << err_pi[i] << endl;
   }
   
   fout.close();

   cout << "Data file written succesfully!" << endl;
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