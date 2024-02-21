#include <iostream>
using namespace std;

#include <omp.h>
#include <iomanip> // setprecision function
#include <cstring> // memcpy function

#include "grid.hh"

#define N_SPECIES 9

int max_cells[N_SPECIES + 1];
int generation[N_SPECIES + 1];

char next_state(int x, int y, int z, char ***grid, long long N) {
    static int ax, ay, az;
    int sum = 0;
    int n_species[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    for(int i = -1; i <= 1; i++) {
        ax = (x + i) % N;
        for(int j = -1; j <= 1; j++) {
            ay = (y + j) % N;
            for(int k = -1; k <= 1; k++) {
                if (i == 0 && j == 0 && k == 0) continue;
                az = (z + k) % N;
                sum += ((int) grid[ax][ay][az] != 0);
                n_species[(int) grid[ax][ay][az]]++;
            }
        }
    }

    x %= N; y %= N; z %= N;
    if ((int) grid[x][y][z]) {
        if (sum > 4 && sum < 14) return (char) grid[x][y][z];
        else return (char) 0;
    }
    else {
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
    for (int i = 1; i <= N_SPECIES; i++) {
        cout << i << " " << max_cells[i] << " " << generation[i] << endl;
    }
}

void init_max_cell() {
    for (int j = 1; j <= N_SPECIES; j++) {
        max_cells[j] = 0;
        generation[j] = 0;
    }
}

void count_cells(char ***grid, long long N) {
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                max_cells[(int) grid[x][y][z]]++;
            }
        }
    }
}

void get_max(int *cells, int gen) {
    for (int i = 1; i <= N_SPECIES; i++) {
        if (cells[i] > max_cells[i]) {
            max_cells[i] = cells[i];
            generation[i] = gen;
        }
    }
}

void simulation(char ***grid, long long N, int generations) {
    char ***new_grid = alloc_grid(N), ***temp;

    // omp_set_num_threads(16);
    
    init_max_cell();
    count_cells(grid, N);

    for (int i = 0; i < generations; i++) {
        
        int cells[N_SPECIES + 1];
        for (int j = 1; j <= N_SPECIES; j++) {
            cells[j] = 0;
        }

        for (int x = 0; x < N; x++) {
            for (int y = 0; y < N; y++) {
                for (int z = 0; z < N; z++) {
                    new_grid[x][y][z] = next_state(x + N, y + N, z + N, grid, N);
                    cells[(int) new_grid[x][y][z]]++;
                }
            }
        }

        get_max(cells, i + 1);
        
        temp = grid;
        grid = new_grid;
        new_grid = temp;

        // cout << "\033[2J\033[1;1H";
        // cout << "Generation " << i << endl;
        // print_result(new_grid, N);
    }
    
}

int main(int argc, char *argv[]) {
    if (argc != 5 || atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atof(argv[3]) < 0 || atof(argv[3]) > 1) {
        cerr << "Usage: " << argv[0] << " <generations (positive integer)> <N (positive integer)> <density (float between 0 and 1)> <seed (integer)>\n";
        return 1;
    }

    int generations = atoi(argv[1]);
    long long N = atoi(argv[2]);
    float density = atof(argv[3]);
    int seed = atoi(argv[4]);

    double exec_time;

    char ***grid = gen_initial_grid(N, density, seed);

    exec_time = -omp_get_wtime();

    simulation(grid, N, generations);

    exec_time += omp_get_wtime();

    cerr << fixed << setprecision(1) << exec_time << "s\n";

    print_result(grid, N);
    return 0;
}
