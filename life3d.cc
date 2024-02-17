#include <iostream>
using namespace std;

#include <omp.h>
#include <iomanip> // setprecision function
#include <cstring> // memcpy function

#include "grid.hh"

#define N_SPECIES 9

typedef struct {
    int n;
    int generation;
} cell;

cell max_cells[N_SPECIES + 1];

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
        cout << i << " " << max_cells[i].n << " " << max_cells[i].generation << endl;
    }
}

void init_cell(cell *cells, int i) {
    for (int j = 0; j <= N_SPECIES; j++) {
        cells[j].n = 0;
        cells[j].generation = i;
    }
}

void count_cells(char ***grid, long long N, cell *cells) {
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                cells[(int) grid[x][y][z]].n++;
            }
        }
    }
}

void get_max(cell *cells) {
    for (int i = 0; i <= N_SPECIES; i++) {
        if (cells[i].n > max_cells[i].n) {
            memcpy(&max_cells[i], &cells[i], sizeof(cell));
        }
    }
}

void simulation(char ***grid, long long N, int generations) {
    char ***new_grid = alloc_grid(N), ***temp;
    
    init_cell(max_cells, 0);
    count_cells(grid, N, max_cells);

    for (int i = 0; i < generations; i++) {
        cell cells[N_SPECIES + 1];
        init_cell(cells, i + 1);

        for (int x = 0; x < N; x++) {
            for (int y = 0; y < N; y++) {
                for (int z = 0; z < N; z++) {
                    new_grid[x][y][z] = next_state(x + N, y + N, z + N, grid, N);
                    cells[(int) new_grid[x][y][z]].n++;
                }
            }
        }

        get_max(cells);
        
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
