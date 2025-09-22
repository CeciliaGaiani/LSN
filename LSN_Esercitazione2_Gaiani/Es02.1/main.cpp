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
#include "random.h"

using namespace std;

/****************************************************************
                      FUNCTION DECLARATIONS
*****************************************************************/

// Cosine calculation
double coseno (double x){
   return (M_PI/2)*cos((M_PI/2)*x);
}

// Function d(x) chosen for importance sampling
double d_of_x(double x){
   return -2*x+2;
}

// Formula for the statistical error (standard deviation of the mean) calculated using the data blocking method
double error(double ave, double ave2, int n){
   if (n<2){
      return 0;
   } else{
      return sqrt((ave2-pow(ave,2))/n);
   }
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
   ofstream fout("integral.dat"), fout_sampl("integral_sampl.dat");

   if (!fout) {
      cerr << "PROBLEM: Unable to open integral.dat!" << endl;
      return 1;
  }

   int M=1000000; // Total number of throws
   int N=100;     // Number of blocks
   int L=M/N;     // Number of throws per block

   vector<double> ave, ave_sampl;
   vector<double> ave2, ave2_sampl;
   vector<double> sum_prog, sum_prog_sampl;
   vector<double> sum2_prog, sum2_prog_sampl;
   vector<double> err_prog, err_prog_sampl;

   for (int i=0; i<N; i++){ // Cycle on the number of blocks      
      double sum = 0;
      double sum_sampl = 0;         
       for (int j=0; j<L; j++){ // Cycle on the number of throws per block
         // Integral calculation sampling a uniform distribution in [0,1)
         double x = rnd.Rannyu();
         sum += coseno(x);

         // Integral calculation using importance sampling (sampling the linear distribution -2x+2)
         double x_sampl = rnd.Line();
         sum_sampl+=coseno(x_sampl)/d_of_x(x_sampl);
       }
       ave.push_back(sum/L);   
       ave_sampl.push_back(sum_sampl/L);    
       ave2.push_back(pow(ave[i],2));
       ave2_sampl.push_back(pow(ave_sampl[i],2));
     }

     // Progressive averages
     for (int i=0; i<N; i++){
      double sum = 0;
      double sum2 = 0;
      double sum_sampl = 0;
      double sum2_sampl = 0;

      for (int j=0; j<i+1; j++){
        sum += ave[j];
        sum2 += ave2[j];

        sum_sampl += ave_sampl[j];
        sum2_sampl += ave2_sampl[j];
       }
 
       sum_prog.push_back(sum/(i+1));
       sum2_prog.push_back(sum2/(i+1));
       err_prog.push_back(error(sum_prog[i], sum2_prog[i], i+1));

       sum_prog_sampl.push_back(sum_sampl/(i+1));
       sum2_prog_sampl.push_back(sum2_sampl/(i+1));
       err_prog_sampl.push_back(error(sum_prog_sampl[i], sum2_prog_sampl[i], i+1));
    }

    for (int i=0; i<N; i++){
      fout << (i+1) << " " << sum_prog[i] << " " << " " << err_prog[i] << endl;
      fout_sampl << (i+1) << " " << sum_prog_sampl[i] << " " << " " << err_prog_sampl[i] << endl;
    }
   
   fout.close();
 
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