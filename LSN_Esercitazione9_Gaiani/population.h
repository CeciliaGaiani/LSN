/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Population__
#define __Population__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h>
#include "random.h"
#include "individual.h"

using namespace std;
using namespace arma;

// The Population class allows to define a population (Armadillo field) of Individuals 
// The generation contains a number of individuals equal to _nind and a reference individual _ref_ind 
// Each individual in the generation is obtained by applying 100 permutations to the reference individual
// The reference inidvidual also serves as a model with which to compare individuals in control functions

/***************************************************************************************************************
                                                DATA MEMBERS
***************************************************************************************************************/

class Population {

private:
  field <Individual> _population;
  Individual _ref_ind;
  int _nind;

/***************************************************************************************************************
                                                  METHODS
***************************************************************************************************************/

public:

// Constructor of a Population object-type
  Population(int nind, Individual& ind, Random& rnd) {
    _population.set_size(nind);
    _ref_ind = ind;
    _nind = nind;

    for (int i = 0; i < _nind; ++i) {
      Individual temp = ind; 
      for (int j = 0; j < 100; ++j) {
          temp.pair_permutation(rnd);
      }
      _population(i) = temp;
    }
  };

// Method returning the i-th element of the population
  Individual& get_individual(int i){
    return _population(i);
  };


// Function that assigns the population Individual "ind" to the i-th position
  void set_individual(const Individual& ind, int i) {
    _population(i) = ind;
  };

// Method returning the number of individuals of the population
  int get_nind() const { return _nind; }

// Algorithm for selecting the individuals:
// - p = 1: completely random selection, with no preference for the best or worst.
// - p < 1: tends to favor the poorest individuals, thus more exploration.
// - p > 1: tends to favor the best individuals, thus more exploitation. This will make best individuals reproduce more.
// The larger p is, the more likely individuals in the top positions (the best, if the population is ORIDNATED by increasing fitness) are chosen.  
  int selection(Random& rnd, double p){
    return int(_nind*pow(rnd.Rannyu(), p));
  };

// Method for sorting individuals in order of increasing distance (best individuals at the beginning)
  void sort_population() {
    vector<Individual> sorted_pop(_nind);

    // It is necessary to work with a vector in order to use the STL sort
    for (int i = 0; i < _nind; ++i) {
        sorted_pop[i] = _population(i);
    }

    sort(sorted_pop.begin(), sorted_pop.end(),
              [](const Individual& a, const Individual& b) {
                  return a.get_distance() < b.get_distance();
              });

    for (int i = 0; i < _nind; ++i) {
        _population(i) = sorted_pop[i];
    }
  };

// Method that checks if the individuals in the population satisfy the constraints of the problem, and eventually corrects problems
  void check_population(){  
    for(int i_ind = 0; i_ind < _nind; i_ind ++){
      if(_population(i_ind).get_ind_city(0) != 1){
        _population(i_ind).set_first_city();
      }
      if(!_population(i_ind).check_index_cities()){
        _population(i_ind).fix_duplicates(_ref_ind);
      }
    }
  };

// Method that saves to a file the path distance associated with the best individual
  void save_best_path(ofstream& fout){
    fout << _population(0).get_distance() << endl;
  };

// Method that computes the average distance of the first half of the population
  void compute_mean_distance(ofstream& fout){
    int n_mean = _nind/2;
    double mean_dist = 0;
    for(int i = 0; i < n_mean; i ++){
      mean_dist+=_population(i).get_distance();
    }
    fout << mean_dist/static_cast<double>(n_mean) << endl;
  };

// Implementation of the crossover algorithm
  void crossover(int i1, int i2, Random& rnd) {
    
    // Parent definition
    Individual& parent1 = _population(i1);
    Individual& parent2 = _population(i2);

    int ncities = parent1.get_ncities();

    // Definition of the interval that will remain unchanged
    int start = static_cast<int>(rnd.Rannyu(1, ncities - 1));  
    int end = static_cast<int>(rnd.Rannyu(start + 1, ncities));

    // Definition of offspring individuals, in which I copy the section [start, end]
    Individual child1 = parent1;
    Individual child2 = parent2;
    for (int i = start; i <= end; ++i) {
        child1.set_city(i, parent1.get_city(i));
        child2.set_city(i, parent2.get_city(i));
    }

    // In the completion of the offspring, it is necessary to complete the paths with the missing cities adding them in the order in which they appear in the consort

    // Complete child1 from parent2
    int idx = (end + 1) % ncities;
    for (int i = 0; i < ncities; ++i) {
        City c = parent2.get_city((end + 1 + i) % ncities);
        if (!child1.contains_city(c.get_ind())) {
            child1.set_city(idx, c);
            idx = (idx + 1) % ncities;
        }
    }

    // Complete child2 from parent1
    idx = (end + 1) % ncities;
    for (int i = 0; i < ncities; ++i) {
        City c = parent1.get_city((end + 1 + i) % ncities);
        if (!child2.contains_city(c.get_ind())) {
            child2.set_city(idx, c);
            idx = (idx + 1) % ncities;
        }
    }

    // Correcting any duplicates
    child1.fix_duplicates(_ref_ind);
    child2.fix_duplicates(_ref_ind);

    // Insert the offspring in the population, overwriting the parents
    _population(i1) = child1;
    _population(i2) = child2;
  };

// Method that creates a new generation, saving the 5 best individuals from the previous one, and applying, with a probability of 90% and 20% respectively, crossovers and mutations on the missing individuals 
void evolve_generation(Random& rnd, double p_sel) {
    field<Individual> new_population(_nind);
    
    sort_population();

    // Save the 5 best individuals of the previous generation
    int elitism_count = 5;
    for (int i = 0; i < elitism_count; i++) {
        new_population(i) = _population(i);
    }

    int idx = elitism_count;
    while (idx < _nind) { // Fill the ne population until the number _nind is reached
        int i1 = selection(rnd, p_sel);
        int i2 = selection(rnd, p_sel);
        while (i2 == i1) i2 = selection(rnd, p_sel); // Just in case i1 = i2 at the first extraction

        Individual parent1 = _population(i1);
        Individual parent2 = _population(i2);

        if (rnd.Rannyu() <= 0.9)
            crossover(i1, i2, rnd);

        if (rnd.Rannyu() < 0.1) parent1.mutate(rnd);
        if (rnd.Rannyu() < 0.1 && idx + 1 < _nind) parent2.mutate(rnd);

        new_population(idx++) = parent1;
        if (idx < _nind) {
            new_population(idx++) = parent2;
        }
    }

    _population = new_population;
    sort_population();
  };
 
};

#endif // __Population__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
