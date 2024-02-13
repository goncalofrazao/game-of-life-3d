#include <iostream>
using namespace std;

#include <omp.h>
#include <iomanip> // setprecision function
#include <unistd.h> // sleep function

#include "grid.hh"

#define N_SPECIES 9

char next_state(int x, int y, int z, char ***grid, long long N) {
    static int dx[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    static int dy[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1};
    static int dz[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1};

    int ax, ay, az;

    int sum = 0;

    // death cell
    if ((int) grid[x][y][z] == 0) {
        int n_species[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        
        // count neighbors and number of each species
        for (int i = 0; i < 26; i++) {
            if (sum > 10 || (i - sum) > 19) {
                break;
            }

            ax = (x + dx[i] + N) % N;
            ay = (y + dy[i] + N) % N;
            az = (z + dz[i] + N) % N;
            if ((int) grid[ax][ay][az] != 0) {
                sum++;
                n_species[(int) grid[ax][ay][az] - 1]++;
            }
        }

        // keep cell dead
        if (sum < 7 || sum > 10) {
            return (char) 0;
        }
        // birth cell
        else {
            int max = 0;
            for (int i = 1; i < N_SPECIES; i++) {
                if (n_species[i] > n_species[max]) {
                    max = i;
                }
            }
            return (char) (max + 1);
        }
    }
    // alive cell
    else {
        for (int i = 0; i < 26; i++) {
            if (sum > 13 || (i - sum) > 21) {
                break;
            }

            ax = (x + dx[i] + N) % N;
            ay = (y + dy[i] + N) % N;
            az = (z + dz[i] + N) % N;
            if ((int) grid[ax][ay][az] != 0) {
                sum++;
            }
        }

        // kill cell
        if (sum < 5 || sum > 13) {
            return (char) 0;
        }
        // keep cell alive
        else {
            return (char) grid[x][y][z];
        }
    }
}

void print_result(char ***grid, long long N) {
    // print board
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                cout << (((int) grid[x][y][z] == 0) ? ' ' : (char) (grid[x][y][z] + '0')) << " ";
            }
            cout << endl;
        }
        cout << endl << endl;
    }
}

void simulation(char ***grid, long long N, int generations) {
    char ***new_grid = alloc_grid(N), ***temp;

    for (int i = 0; i < generations; i++) {
        for (int x = 0; x < N; x++) {
            for (int y = 0; y < N; y++) {
                for (int z = 0; z < N; z++) {
                    new_grid[x][y][z] = next_state(x, y, z, grid, N);
                }
            }
        }
        
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
