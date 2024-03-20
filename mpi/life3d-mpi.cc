#include <iostream>
using namespace std;

#include <mpi.h>
// #include <omp.h>

#include <cstring>	// memcpy function
#include <iomanip>	// setprecision function

#include "grid.hh"

#define N_SPECIES 9

int local_cells[N_SPECIES + 1];
int max_cells[N_SPECIES + 1];
int generation[N_SPECIES + 1];

char next_state(int x, int y, int z, char ***grid, long long N) {
	static int ax, ay, az;
	int sum = 0;
	int n_species[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	for (int i = -1; i <= 1; i++) {
		ax = x + i;
		for (int j = -1; j <= 1; j++) {
			ay = (y + j) % N;
			for (int k = -1; k <= 1; k++) {
				if (i == 0 && j == 0 && k == 0) continue;
				az = (z + k) % N;
				sum += ((int)grid[ax][ay][az] != 0);
				n_species[(int)grid[ax][ay][az]]++;
			}
		}
	}

	x %= N;
	y %= N;
	z %= N;

	if ((int)grid[x][y][z]) {
		if (sum > 4 && sum < 14)
			return (char)grid[x][y][z];
		else
			return (char)0;
	} else {
		if (sum > 6 && sum < 11) {
			int max = 1;
			for (int i = 2; i <= N_SPECIES; i++)
				if (n_species[i] > n_species[max]) max = i;
			return (char)max;
		} else
			return (char)0;
	}
}

void print_result(char ***grid, long long N) {
	for (int i = 1; i <= N_SPECIES; i++) {
		cout << i << " " << max_cells[i] << " " << generation[i] << endl;
	}
}

void init_max_cells() {
	for (int i = 1; i <= N_SPECIES; i++) {
		generation[i] = 0;
		max_cells[i] = 0;
	}
}

void count_cells(char ***grid, long long N, long long block_size) {
	for (int x = 1; x <= block_size; x++) {
		for (int y = 0; y < N; y++) {
			for (int z = 0; z < N; z++) {
				local_cells[(int)grid[x][y][z]]++;
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

void simulation(char ***grid, long long N, int generations, int rank, int size) {
	long long block_size = BLOCK_SIZE(rank, size, N);
	char ***new_grid = alloc_grid(N, block_size + 2), ***temp;
	long long N2 = N * N;

	MPI_Request reqs[4];

	if (!rank) init_max_cells();
	count_cells(grid, N, block_size);
	MPI_Reduce(local_cells, max_cells, N_SPECIES + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	for (int i = 0; i < generations; i++) {
		int cells[N_SPECIES + 1] = {0};

		MPI_Irecv(new_grid[0][0], N2, MPI_CHAR, PREV_BLOCK(rank, size), 1, MPI_COMM_WORLD, reqs + 2);
		MPI_Irecv(new_grid[block_size + 1][0], N2, MPI_CHAR, NEXT_BLOCK(rank, size), 0, MPI_COMM_WORLD, reqs + 3);

		for (int y = 0; y < N; y++) {
			for (int z = 0; z < N; z++) {
				new_grid[1][y][z] = next_state(1, y + N, z + N, grid, N);
				cells[(int)new_grid[1][y][z]]++;
			}
		}
		MPI_Isend(new_grid[1][0], N2, MPI_CHAR, PREV_BLOCK(rank, size), 0, MPI_COMM_WORLD, reqs);

		for (int y = 0; y < N; y++) {
			for (int z = 0; z < N; z++) {
				new_grid[block_size][y][z] = next_state(block_size, y + N, z + N, grid, N);
				cells[(int)new_grid[block_size][y][z]]++;
			}
		}
		MPI_Isend(new_grid[block_size][0], N2, MPI_CHAR, NEXT_BLOCK(rank, size), 1, MPI_COMM_WORLD, reqs + 1);

		for (int x = 2; x < block_size; x++) {
			for (int y = 0; y < N; y++) {
				for (int z = 0; z < N; z++) {
					new_grid[x][y][z] = next_state(x, y + N, z + N, grid, N);
					cells[(int)new_grid[x][y][z]]++;
				}
			}
		}

		MPI_Reduce(cells, local_cells, N_SPECIES + 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (!rank) get_max(local_cells, i + 1);

		MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

		temp = grid;
		grid = new_grid;
		new_grid = temp;
	}
}

int main(int argc, char *argv[]) {
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (argc != 5 || atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atof(argv[3]) < 0 || atof(argv[3]) > 1) {
		cerr << "Usage: " << argv[0] << " <generations> <N> <density> <seed>\n";
		return 1;
	}

	int generations = atoi(argv[1]);
	long long N = atoi(argv[2]);
	float density = atof(argv[3]);
	int seed = atoi(argv[4]);

	double exec_time;

	char ***grid = gen_initial_grid(N, density, seed, rank, size);

	MPI_Barrier(MPI_COMM_WORLD);
	exec_time = -MPI_Wtime();

	simulation(grid, N, generations, rank, size);

	MPI_Barrier(MPI_COMM_WORLD);
	exec_time += MPI_Wtime();

	if (!rank) cerr << fixed << setprecision(1) << exec_time << "s\n";

	if (!rank) print_result(grid, N);

	MPI_Finalize();
}
