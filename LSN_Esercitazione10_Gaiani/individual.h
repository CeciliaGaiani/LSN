/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Individual__
#define __Individual__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h>
#include "city.h"
#include "random.h"

using namespace std;
using namespace arma;

// The individual class allows the construction of the succession of 34 cities to be visited (chromosomes of the genetic algorithm)
// An Individual is an Armadillo field, characterized by a number of cities _ncities (fixed to 34)
// Initially, the cities are randomly arranged in the Individual

class Individual {

/***************************************************************************************************************
                                                DATA MEMBERS
***************************************************************************************************************/

private:
  field <City> _individual;    
  int _ncities = 34;      // Number of cities to visit

/***************************************************************************************************************
                                                  METHODS
***************************************************************************************************************/

public:

// Constructor (no arguments)
  Individual(){
    _individual.set_size(_ncities);
    for(int i = 0; i < _ncities; i++){
      _individual(i) = City(0.0, 0.0, i+1);
    } 
  };

// Constructor (arguments)
  Individual(double x, double y){
    _individual.set_size(_ncities);
    for(int i = 0; i < _ncities; i++){
      _individual(i) = City(x, y, i+1);
    } 
  };

// Destructor  
   ~Individual(){;};

// Copy operator (the default one should be fine, however we implement it to make sure we have no problems running the algorithm)
  Individual(const Individual& other)
        : _individual(other._individual), _ncities(other._ncities) {};

// Assignment operator (the default one should be fine, however we implement it to make sure we have no problems running the algorithm)   
  Individual& operator=(const Individual& other) {
    if (this != &other) {
        _individual = other._individual;
        _ncities = other._ncities;
      }
    return *this;
  };


// Operator that initializes the individual with cities randomly distributed over a sphere of radius "radius"
  void initialize_individual_circumference (Random& rnd, double radius){
    _individual.set_size(_ncities);
    for (int i = 0; i < _ncities; ++i) {
      double theta = rnd.Rannyu(0, 2  * M_PI);
      double x = radius * cos(theta);
      double y = radius * sin(theta);
      _individual(i) = City(x, y, i + 1);
    }
  };

// Operator that initializes the individual with cities randomly distributed in a square of side "side"
  void initialize_individual_square (Random& rnd, double side){
    _individual.set_size(_ncities);
    for (int i = 0; i < _ncities; ++i) {
      double x = rnd.Rannyu(-side/2+0.1, side/2);
      double y = rnd.Rannyu(-side/2+0.1, side/2);
      _individual(i) = City(x, y, i + 1);
    }
  };

// Method that allows to initialize the individual reading data from a file
  void initialize_individual_file (const string& filename) {
    _individual.set_size(_ncities);
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Errore nell'apertura del file: " << filename << endl;
        return;
    }
    int i = 0;
    double x, y;
    while (infile >> x >> y) {
        _individual(i) = City(x, y, i + 1);
        i ++; 
    }
    infile.close();
  }; 

// Method for assigning the individual the city "city" to the i-th component
  void set_individual(const City& city, int i){
    if (i >= 0 && i < _ncities) {
        _individual(i) = city;
    } else {
        cerr << "[set_individual] Errore: indice fuori dal range!" << endl;
    }
  };

// Method for defining the individual's city number
  void set_ncities(int ncities){
    _ncities = ncities;
    _individual.set_size(ncities);
  }

// Method that allows to set _nind based on the number of data read from a file
  void set_ncities(const string& filename){
    ifstream file_in;
    file_in.open(filename);

    if (!file_in.is_open()) {
        cout << "Errore nell'apertura del file." << endl;
        return;
        } else {
          int conta = 0;
          while(!file_in.eof()) {
          string s;
		      conta++;
		      getline(file_in,s);
		    }
        _ncities = conta;
        _individual.set_size(conta);
        file_in.close();
     }
  };

// Method returning the i-th city of the individual
  const City& get_city(int i) const {
    return _individual(i);
  }; 

// Method returning the number of cities of the individual
  int get_ncities() const {
    return _ncities;
  };

// Method that returns the index of the city at the i-th position of the individual
  int get_ind_city(int i) const {
    return _individual(i).get_ind();
  };

// Method that allows the individual to insert the city “city” in the i-th position
  void set_city(int i, const City& city) {
    _individual(i) = city;
  };

// Compute the total distance of the path defined by the cities in the individual
  double get_distance() const {
    double d = 0;
    for(int i_city=0; i_city<_ncities-1; i_city++){
      d += pow(_individual(i_city).get_x() - _individual(i_city+1).get_x(), 2) + pow(_individual(i_city).get_y() - _individual(i_city+1).get_y(), 2);
    }
      d+= pow(_individual(_ncities-1).get_x() - _individual(0).get_x(), 2) + pow(_individual(_ncities-1).get_y() - _individual(0).get_y(), 2);
      return d;
  }; 

// Method that saves in a file index and coordinates of the individual cities 

void save_individual(const string filename){
  ofstream fout (filename);
    if (!fout) {
        cerr << "Errore: impossibile aprire i file di output!" << endl;
    }
    for (int i = 0; i < _ncities; ++i) {
        City c = get_city(i);
        fout << c.get_ind() << " " << c.get_x() << " " << c.get_y() << endl;
    }
  };

/***************************************************************************************************************
Definition of 4 possible mutations that can be realized on the individual
***************************************************************************************************************/
  
// Shift of two cities (the first city is fixed)
  void pair_permutation(Random& rnd){
    int i1 = static_cast<int>(rnd.Rannyu(1, _ncities));
    int i2 = static_cast<int>(rnd.Rannyu(1, _ncities));

    while (i1 == i2) {
        i2 = static_cast<int>(rnd.Rannyu(1, _ncities));
    }
    swap(_individual(i1), _individual(i2));
  };   

// Shift of a block of m elements of +n positions
  void block_shift_right(Random& rnd) {                              
    int i0 = static_cast<int>(rnd.Rannyu(1, _ncities - 2));             
    int m = static_cast<int>(rnd.Rannyu(1, _ncities - i0 - 1));        
    int n = static_cast<int>(rnd.Rannyu(1, _ncities - i0 - m));

    vector<City> block(m);
    for (int i = 0; i < m; i++) {
        block[i] = _individual(i0 + i);
    }

    for (int i = 0; i < n; i ++) {
        _individual(i0 + i) = _individual(i0 + m + i);
    }

    for (int i = 0; i < m; ++i) {
        _individual[i0 + n + i] = block[i];
    }
  };

// Swap of two blocks of m elements
  void block_swap(Random& rnd) {
    int i1 = static_cast<int>(rnd.Rannyu(1, _ncities - 2));
    int max_i2 = min(_ncities - 1, i1 + (_ncities - i1 - 1) / 2 + 1);

    int i2 = static_cast<int>(rnd.Rannyu(i1 + 1, max_i2 + 1));

    int dist = i2 - i1;
    int m_max = min(dist, _ncities - i2);

    int m = static_cast<int>(rnd.Rannyu(1, m_max + 1));

    for (int i = 0; i < m; ++i) {
       swap(_individual(i1 + i), _individual(i2 + i));
    }
  };


// Inversion of the order of appearance of the cities in the individual (the first and the last are set)
  void inversion (Random& rnd){
    int i1 = static_cast<int>(rnd.Rannyu(1, _ncities - 2));
    int i2 = static_cast<int>(rnd.Rannyu(i1 + 1, _ncities - 1));
    int dist = i2 - i1;
    int m = dist + 1;

    vector<City> block(m);
    for(int i = 0; i < m; i ++){
      block[i] = _individual(i2 - i);
    }
    for (int i = 0; i < m; ++i) {
        _individual[i1 + i] = block[i];
    }
  };

// Method performing a mutation process on an individual: according to the value of a random number uniformly distributed in [0,1), it realizes a mutation among the possible ones
  void mutate(Random& rnd) {
    double mutation = rnd.Rannyu();

    if (mutation < 0.45) block_shift_right(rnd);
    else if (mutation < 0.70) inversion(rnd);
    else if (mutation < 0.90) block_swap(rnd);
    else pair_permutation(rnd);
  };

/***************************************************************************************************************
================================================================================================================
***************************************************************************************************************/ 

/***************************************************************************************************************
Definition of functions to check if individuals in the population satisfy the constraints of the problem: 
each city must be visited only once, the first city (of index 1) must be untouched by the mutations
***************************************************************************************************************/

// Function that searches for the city of index 1 in the individual and possibly places it at the top
  void set_first_city() {
    if(_individual(0).get_ind()!=1){
      for(int i_city = 0; i_city < _ncities; i_city ++){
        if(_individual(i_city).get_ind() == 1){
          swap(_individual(i_city), _individual(0));
          break;
        }
      }
    }
  };

// Function that checks that in the vector _individual all cities have different index
  bool check_index_cities() const {
    for(int i = 0; i < _ncities - 1; i++) {
        for(int j = i + 1; j < _ncities; j++) {
            if(_individual(i).get_ind() == _individual(j).get_ind()) {
                return false; 
            }
        }
    }
    return true;
  };

// Function that checks whether the index city i is present in the individual
  bool contains_city(int ind) {
    for (int i = 0; i < _ncities; ++i) {
      if (_individual(i).get_ind() == ind) {
        return true;
      }
    }
    return false;
  };
  
// Method that corrects the individual if it presents repeated cities (I already assume that the check_index_cities function has found equal indexes)
  void fix_duplicates (Individual& ref_ind) {
    int ncities = get_ncities();
    vector<int> count(ncities+1, 0.0);
    vector<int> duplicate_positions;

  for (int i = 0; i < ncities; i++) {
    int ind = _individual(i).get_ind();    // Index of the i-th city
    count[ind]++;                          // The i-th component of the "count" vector is incremented by one each time if index i is found              
    if (count[ind] > 1) {
      duplicate_positions.push_back(i);    // Save the positions to fix
     }
   }

  vector<int> missing_indices;
  for (int i = 1; i <= ncities; i++) {
    if (count[i] == 0) { // If any of the count components is zero, it means that the city of corresponding index does not appear in the individual, so the position is saved
      missing_indices.push_back(i);
     }
   }

  // Cicle on the positions to fix
  size_t min_size = min(duplicate_positions.size(), missing_indices.size()); // Just in case duplicate_positions and missing_indices have different dimensions
  for(size_t i = 0; i < min_size; i++) {
    int pos = duplicate_positions[i]; // Position of the city to fix in the _individual vector
    int replace_index = missing_indices[i]; // Index of the missing city
    City replacement = ref_ind.get_city(replace_index - 1); // Search the reference individual for the city to replace
    set_individual(replacement, pos);
    }
  };

/***************************************************************************************************************
================================================================================================================
***************************************************************************************************************/   

// Method that extracts the sequence of city indexes and stores them in a vector
  vector<int> extract_city_indices(){
    vector <int> indexes (_ncities);
    for(int i = 0; i < _ncities; i ++){
      indexes[i] = _individual(i).get_ind();
    }
    return indexes;
  };

};

#endif // __Individual__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
