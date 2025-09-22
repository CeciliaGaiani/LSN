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

// Function computing the squared distance between two points
double quad_dist(double x1, double x2, double y1, double y2, double z1, double z2) {
   return pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2);
}

// Formula for the statistical error (standard deviation of the mean) calculated using the data blocking method
double error(double ave, double ave2, int n) {
   if (n<2)
      return 0;
   else
      return sqrt((ave2 - pow(ave, 2)) / n);
}

// Function computing a discrete random walk in 3D
void RandomWalkDiscrete3D(int N_walks, int N_steps, int N_blocks, Random &rnd, vector<double> &r_prog, vector<double> &r_err_prog) {
   
   int L = N_walks / N_blocks; // The number of walks per block is defined as the total number of RW to be performed, divided by the number of blocks 

   vector<double> r2(N_steps, 0.0);   
   vector<double> r2_2(N_steps, 0.0);

   for (int i_block = 0; i_block < N_blocks; i_block++) { // Cycle on the number of blocks
      vector<double> r2_block(N_steps, 0.0);

      for (int i_walk = 0; i_walk < L; i_walk++) { // Cycle on the number of walks per block

         // RW is a vector of (N_steps + 1) positions in 3D space. Each position is a vector with 3 components (x, y, z), all initialized to zero.
         // RW[0] corresponds to the origin. At each step, a random direction among the 6 possible (±x, ±y, ±z) is selected,
         // and the walker moves by one lattice unit (a = 1) in that direction, updating the position relative to the previous one.

         vector<vector<double>> RW(N_steps + 1, vector<double>(3, 0.0));

         for (int i_step = 1; i_step <= N_steps; i_step++) { // Cycle on the number of steps per walk
            int val = int(rnd.Rannyu() * 6);

            if(val == 0){         //+x
               RW[i_step][0]=RW[i_step-1][0]+1;
               RW[i_step][1]=RW[i_step-1][1];
               RW[i_step][2]=RW[i_step-1][2];
            }else if(val == 1){   //-x
               RW[i_step][0]=RW[i_step-1][0]-1;
               RW[i_step][1]=RW[i_step-1][1];
               RW[i_step][2]=RW[i_step-1][2];
            }else if(val == 2){   //+y
               RW[i_step][0]=RW[i_step-1][0];
               RW[i_step][1]=RW[i_step-1][1]+1;
               RW[i_step][2]=RW[i_step-1][2];
            }else if(val == 3){   //-y
               RW[i_step][0]=RW[i_step-1][0];
               RW[i_step][1]=RW[i_step-1][1]-1;
               RW[i_step][2]=RW[i_step-1][2];
            }else if(val == 4){   //+z
               RW[i_step][0]=RW[i_step-1][0];
               RW[i_step][1]=RW[i_step-1][1];
               RW[i_step][2]=RW[i_step-1][2]+1;
            }else if(val == 5){   //-z
               RW[i_step][0]=RW[i_step-1][0];
               RW[i_step][1]=RW[i_step-1][1];
               RW[i_step][2]=RW[i_step-1][2]-1;
            }

            // Compute the squared distance from the origin after this step
            double d2 = quad_dist(RW[i_step][0], 0.0, RW[i_step][1], 0.0, RW[i_step][2], 0.0);

            // Accumulate d² in the block-average vector at index (i_step - 1)
            r2_block[i_step - 1] += d2;
         }
      }

      // Compute single block averages for each step
      for (int i_step = 0; i_step < N_steps; i_step++) {
         double r2_avg_block = r2_block[i_step]/L;
         r2[i_step] += r2_avg_block;
         r2_2[i_step] += r2_avg_block*r2_avg_block;
        }
    }

   // Compute global averages and statistical uncertainties for each step
   for (int i_step= 0; i_step < N_steps; i_step++) {
      double r2_avg = r2[i_step]/N_blocks; 
      double r2_2_avg = r2_2[i_step]/N_blocks;
      r_prog[i_step] = sqrt(r2_avg);
      r_err_prog[i_step] = 0.5 * error(r2_avg, r2_2_avg, N_blocks)/r_prog[i_step]; // Use of error propagation
    }
}

// Function computing a continuous random walk in 3D
void RandomWalkContinuum3D(int N_walks, int N_steps, int N_blocks, Random &rnd, vector<double> &r_prog, vector<double> &r_err_prog) {
   
   int L = N_walks / N_blocks;

   vector<double> r2(N_steps, 0.0);
   vector<double> r2_2(N_steps, 0.0);

   for (int i_block = 0; i_block < N_blocks; i_block++) {
      vector<double> r2_block(N_steps, 0.0);

      for (int i_walk = 0; i_walk < L; i_walk++) {
         vector<vector<double>> RW(N_steps + 1, vector<double>(3, 0.0));

         for (int i_step = 1; i_step <= N_steps; i_step++) {
         // θ is not sampled uniformly in [0, π] because the area on the surface of the sphere associated with an interval Δθ varies 
         // with sinθ. If θ were sampled uniformly, there would be a higher density of points near the poles (close to 0 and π)
            double cos_theta = rnd.Rannyu(-1, 1);
            double phi = rnd.Rannyu(0, 2 * M_PI);

            double sin_theta = sqrt(1 - cos_theta*cos_theta);
            RW[i_step][0] = RW[i_step - 1][0] + sin_theta * cos(phi);
            RW[i_step][1] = RW[i_step - 1][1] + sin_theta * sin(phi);
            RW[i_step][2] = RW[i_step - 1][2] + cos_theta;

            double d2 = quad_dist(RW[i_step][0], 0.0, RW[i_step][1], 0.0, RW[i_step][2], 0.0);
            r2_block[i_step - 1] += d2;
         }
      }
      // Compute single block averages for each step
      for (int i_step = 0; i_step < N_steps; i_step++) {
         double r2_avg_block = r2_block[i_step]/L;
         r2[i_step] += r2_avg_block;
         r2_2[i_step] += r2_avg_block*r2_avg_block;
        }
    }

   // Compute global averages and statistical uncertainties for each step
   for (int i_step= 0; i_step < N_steps; i_step++) {
      double r2_avg = r2[i_step]/N_blocks; 
      double r2_2_avg = r2_2[i_step]/N_blocks;
      r_prog[i_step] = sqrt(r2_avg);
      r_err_prog[i_step] = 0.5 * error(r2_avg, r2_2_avg, N_blocks)/r_prog[i_step]; // Use of error propagation
    }
}

/****************************************************************
                               MAIN
*****************************************************************/

int main(int argc, char *argv[]) {
   // Random number generator initialization
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes.txt");
   if (Primes.is_open()) {
      Primes >> p1 >> p2;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
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
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   int N_walks = 10000;
   int N_steps = 100;
   int N_blocks = 100;

   vector<double> r_prog_discr(N_steps, 0.0);
   vector<double> r_err_prog_discr(N_steps, 0.0);

   vector<double> r_prog_cont(N_steps, 0.0);
   vector<double> r_err_prog_cont(N_steps, 0.0);

   RandomWalkDiscrete3D(N_walks, N_steps, N_blocks, rnd, r_prog_discr, r_err_prog_discr);
   RandomWalkContinuum3D(N_walks, N_steps, N_blocks, rnd, r_prog_cont, r_err_prog_cont);

   ofstream r_mean_discr("RW_discrete_distance.dat"), r_mean_cont("RW_continuum_distance.dat");
   for (int i = 0; i < N_steps; i++) {
      r_mean_discr << i << " " << r_prog_discr[i] << " " << r_err_prog_discr[i] << endl;
      r_mean_cont << i << " " << r_prog_cont[i] << " " << r_err_prog_cont[i] << endl;
   }

   r_mean_discr.close();
   r_mean_cont.close();
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