/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __City__
#define __City__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h>
#include "random.h"

using namespace std;
using namespace arma;

// The City class allows to define elements of type “City” (genes of the genetic algorithm), characterized by their corresponding x and y coordinates, and methods for accessing and possibly modifying them

/***************************************************************************************************************
                                                DATA MEMBERS
***************************************************************************************************************/

class City {

private:
  double _x, _y;         // x and y coordinates of the city
  int _city_ind;         // Each city is characterized by a specific index, ranging from 1 to 34
 

/***************************************************************************************************************
                                                  METHODS
***************************************************************************************************************/

public: 

// Constructor (no arguments)
  City() {
    _x = 0.0;
    _y = 0.0;
    _city_ind = 0.0;
  };

// Constructor (arguments)
  City(double x, double y, int ind) {
    _x = x;
    _y = y;
    _city_ind = ind;
  };

// Destructor
  ~City(){;};

// Methods to access coordinates and index of a city
  double get_x () const {
    return _x;
  };

  double get_y () const {
    return _y;
  };

  int get_ind () const {
    return _city_ind;
  };

// Methods to assign coordinates or index to a city
  void set_x (double x) {
    _x = x;
  };

  void set_y (double y) {
    _y = y;
  };

  void set_ind (int ind) {
    _city_ind = ind;
  };

};

#endif // __City__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
