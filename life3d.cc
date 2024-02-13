#include <iostream>
using namespace std;

#include <omp.h>
#include <iomanip> // setprecision function
#include <unistd.h> // sleep function

#include "grid.hh"

void simulation() {
    return;
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

    simulation();

    exec_time += omp_get_wtime();

    cerr << fixed << setprecision(1) << exec_time << "s\n";

    print_result(grid, N);
    return 0;
}
