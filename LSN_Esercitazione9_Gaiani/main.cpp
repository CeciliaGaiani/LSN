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
#include "random.h"
#include "individual.h"
#include "population.h"

using namespace std;

int main() {
    Random rnd;
    string seed_file = "seed.in";
    string primes_file = "Primes.txt";
    rnd.Initialize(seed_file, primes_file); 

    ofstream fout_c_best("best_path_c.dat"), fout_s_best("best_path_s.dat");
    ofstream fout_c_mean("mean_path_c.dat"), fout_s_mean("mean_path_s.dat");

    const int ngen = 1000;
    const int npop = 500;
    const double p_sel = 2;

    Individual cities_circle;
    Individual cities_square;
    cities_circle.initialize_individual_circumference(rnd, 1);
    cities_square.initialize_individual_square(rnd, 1);

    Population population_circumference(npop, cities_circle, rnd);
    Population population_square(npop, cities_square, rnd);
    population_circumference.sort_population();
    population_square.sort_population();

    for (int gen = 0; gen < ngen; ++gen) {
        cout << "Generation n° " << gen + 1 << endl;

        population_circumference.evolve_generation(rnd, p_sel);
        population_square.evolve_generation(rnd, p_sel);

        population_circumference.save_best_path(fout_c_best);
        population_square.save_best_path(fout_s_best);

        population_circumference.compute_mean_distance(fout_c_mean);
        population_square.compute_mean_distance(fout_s_mean);
    }

    population_circumference.sort_population();
    population_square.sort_population();

    Individual best_circle = population_circumference.get_individual(0);
    Individual best_square = population_square.get_individual(0);

    best_circle.save_individual("best_path_circle.dat");
    best_square.save_individual("best_path_square.dat");

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
