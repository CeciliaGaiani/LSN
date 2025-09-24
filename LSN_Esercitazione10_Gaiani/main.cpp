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
#include <mpi.h>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "individual.h"
#include "population.h"

using namespace std;

int main(int argc, char* argv[]) {

    // MPI initialization
    int num_procs, rank; 
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Decide wether perform migration or not
    bool enable_migration = true; 

    // Random number generator initialization
    Random rnd;
    string seed_file = "seed.in";
    string primes_file = "Primes.txt";
    rnd.Initialize(seed_file, primes_file, rank); 

    const int ngen = 2000;   // Number of generations
    const int npop = 500;    // Number of individuals in each generation
    const int N_migr = 50;   // Each N_migr generations perform a migration
    const double p_sel = 2;  // Selection parameter

    string label = enable_migration ? "mig" : "nomig";

    // Output files
    ofstream best_distance_file("best_distance_rank_" + to_string(rank) + "_" + label + ".dat");
    ofstream mean_distance_file("mean_distance_rank_" + to_string(rank) + "_" + label + ".dat");

    // Read cities from file
    Individual cities;
    cities.set_ncities("cap_prov_ita.dat"); // Number of cities
    cities.initialize_individual_file("cap_prov_ita.dat"); // Read reference individual from file

    // Population
    Population population_prov(npop, cities, rnd);
    population_prov.sort_population();

    // Evolve generation
    for (int gen = 0; gen < ngen; ++gen) {
        population_prov.evolve_generation(rnd, p_sel);
        population_prov.sort_population();
        population_prov.check_population();

        population_prov.save_best_path(best_distance_file);
        population_prov.compute_mean_distance(mean_distance_file);
        Individual best_ind = population_prov.get_individual(0);

        if (enable_migration && gen % N_migr == 0 && gen > 0) {
            population_prov.migrate(best_ind, rank, num_procs);
        }
    }

    // For each process, save the best path after evolution
    population_prov.sort_population();
    Individual best_individual = population_prov.get_individual(0);
    double local_best_distance = best_individual.get_distance();

    // This vector will contain the best distances for each process
    vector<double> all_best_distances(num_procs);
    // All processes send their local_best_distance to rank 0, which collects them in all_best_distances
    MPI_Gather(&local_best_distance, 1, MPI_DOUBLE,
               all_best_distances.data(), 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // Search the rank corresponding to the best path
    int best_rank = 0;
    if (rank == 0) {
        double min_dist = all_best_distances[0];
        for (int i = 1; i < num_procs; i++) {
            if (all_best_distances[i] < min_dist) {
                min_dist = all_best_distances[i];
                best_rank = i;
            }
        }
    }

    // All processes must know which is the rank with the best path
    MPI_Bcast(&best_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // The best rank saves his results
    if (rank == best_rank) {
        string overall_path_file = "best_overall_path_" + label + ".dat";
        string overall_dist_file = "best_overall_distance_" + label + ".dat";

        population_prov.sort_population();
        Individual best_ind = population_prov.get_individual(0);
        best_ind.save_individual(overall_path_file);
        ofstream out(overall_dist_file);
        out << best_ind.get_distance() << endl;
        out.close();
    }

    MPI_Finalize();
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
