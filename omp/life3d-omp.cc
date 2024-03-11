#include <iostream>
using namespace std;

#include <cstring> // memcpy function
#include <iomanip> // setprecision function
#include <omp.h>

#include "grid.hh"

#define N_SPECIES 9

int max_cells[N_SPECIES + 1]; // Global array to store the maximum number of cells for each species...
int generation[N_SPECIES + 1]; //... and the generation in which it occurred

char next_state(int x, int y, int z, char ***grid, long long N) {
    // Routine to compute the next state of a cell based on its 26 neighbors and return the new state
    static int ax, ay, az;
    int sum = 0; // Variable to store the number of live neighbors
    int n_species[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // Array to store the count of neighbors for each species

    for(int i = -1; i <= 1; i++) { // Iterate over the 26 neighbors
        ax = (x + i) % N;
        for(int j = -1; j <= 1; j++) {
            ay = (y + j) % N;
            for(int k = -1; k <= 1; k++) {
                if (i == 0 && j == 0 && k == 0) continue;
                az = (z + k) % N;
                sum += ((int) grid[ax][ay][az] != 0); // Increment the number of live neighbors
                n_species[(int) grid[ax][ay][az]]++; // Increment the count of neighbors for the current species
            }
        }
    }

    x %= N; y %= N; z %= N; // Normalize the coordinates

    if ((int) grid[x][y][z]) { // If the cell is alive
        if (sum > 4 && sum < 14) return (char) grid[x][y][z];
        else return (char) 0;
    }
    else { // If the cell is dead
        if (sum > 6 && sum < 11) {
            int max = 1;
            for (int i = 2; i <= N_SPECIES; i++)
                if (n_species[i] > n_species[max]) max = i;
            return (char) max;
        }
        else return (char) 0;
    }
}

void print_result(char ***grid, long long N) {
    // Print the result to stdout
    for (int i = 1; i <= N_SPECIES; i++) {
        cout << i << " " << max_cells[i] << " " << generation[i] << endl;
    }
}

void init_max_cell() {

    // Initialize the max_cells and generation arrays (all elements are set to 0)
    for (int j = 1; j <= N_SPECIES; j++) {
        max_cells[j] = 0;
        generation[j] = 0;
    }
}

void count_cells(char ***grid, long long N) {

    // Parallelize the counting of cells in initial grid. The grid is divided in order to optimize the computation.
    // The reduction clause is used to update the max_cells array safely.
    #pragma omp parallel for reduction(+:max_cells)
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                max_cells[(int) grid[x][y][z]]++;
            }
        }
    }
}

void get_max(int *cells, int gen) {

    // Update the max_cells and generation arrays, based on the current generation cells count
    for (int i = 1; i <= N_SPECIES; i++) {
        if (cells[i] > max_cells[i]) {
            max_cells[i] = cells[i];
            generation[i] = gen;
        }
    }
}

void simulation(char ***grid, long long N, int generations) {
    char ***new_grid = alloc_grid(N), ***temp; // Allocate the new auxiliary grid

    // omp_set_num_threads(16); // Set the number of threads 
    init_max_cell(); // Initialize the max_cells and generation arrays
    count_cells(grid, N); // Count the initial cells (generation 0)

    for (int i = 0; i < generations; i++) {
        
        int cells[N_SPECIES + 1] = {0}; // Initialize the count of cells for each species in the current generation

        // Parallelize the simulation. The cube is divided in order to optimize the computation.
        // The reduction clause is used update the cells species count array safely.
        #pragma omp parallel for reduction(+:cells) 
        for (int x = 0; x < N; x++) {
            for (int y = 0; y < N; y++) {
                for (int z = 0; z < N; z++) {
                    new_grid[x][y][z] = next_state(x + N, y + N, z + N, grid, N);
                    cells[(int) new_grid[x][y][z]]++;
                }
            }
        }
        // There is an implicit barrier here

        get_max(cells, i + 1); // Update the max_cells and generation arrays
        
        // Swap the grid and the new_grid pointers
        temp = grid;
        grid = new_grid;
        new_grid = temp;
    }
}

int main(int argc, char *argv[]) {

    // Check if the input is valid
    if (argc != 5 || atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atof(argv[3]) < 0 || atof(argv[3]) > 1) {
        cerr << "Usage: " << argv[0] << " <generations (positive integer)> <N (positive integer)> <density (float between 0 and 1)> <seed (integer)>\n";
        return 1;
    }

    // Parse the input
    int generations = atoi(argv[1]);
    long long N = atoi(argv[2]);
    float density = atof(argv[3]);
    int seed = atoi(argv[4]);

    // Allocate and generate the initial grid
    char ***grid = gen_initial_grid(N, density, seed);

    // Run the simulation and measure the execution time
    double exec_time;
    exec_time = -omp_get_wtime();

	simulation(grid, N, generations);

	exec_time += omp_get_wtime();

    // Print the execution time to stderr
    cerr << fixed << setprecision(1) << exec_time << "s\n";

    // Print the result to stdout
    print_result(grid, N);
    return 0;
}
