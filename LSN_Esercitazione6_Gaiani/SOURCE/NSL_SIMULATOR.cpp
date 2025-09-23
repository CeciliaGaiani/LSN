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
#include "system.h"

using namespace std;

int main(int argc, char *argv[]) {

  System SYS;
  SYS.initialize();
  SYS.initialize_properties();

  // If _therma = 1, perform the equilibration test
  if (SYS.get_therma()) {
    SYS.block_reset(0); 
    for (int i = 0; i < SYS.get_nbl(); i++) {
      for (int j = 0; j < SYS.get_nsteps(); j++) { 
        SYS.step();
        SYS.measure();
      }
      SYS.averages(i + 1);
      SYS.block_reset(i + 1);
    }
    SYS.finalize();
  }

  // Otherwise, proceed with the simulation
  else {
    double step_T = 0.05;
    double T;

    // Temperature scaling
    for (T = 2.0; T >= 0.45; T -= step_T) {
      SYS.set_temp(T);
      SYS.reset_averages();

      // Thermalization
      for (int iblk = 0; iblk < 6000; iblk++) {
        SYS.step();
      }

      // Simulation and measures
      for (int i = 0; i < SYS.get_nbl(); i++) {
        SYS.block_reset(i+1);
        for (int j = 0; j < SYS.get_nsteps(); j++) {
          SYS.step();
          SYS.measure();
        }
        SYS.averages(i + 1);
        SYS.block_reset(i + 1);
      }

      SYS.finalize();
    }
  }

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
