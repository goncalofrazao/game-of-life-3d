#include <iostream>
using namespace std;

#include <mpi.h>

#include "grid.hh"

#define N_SPECIES 9

unsigned int seed;

void init_r4uni(int input_seed) { seed = input_seed + 987654321; }

float r4_uni() {
	int seed_in = seed;

	seed ^= (seed << 13);
	seed ^= (seed >> 17);
	seed ^= (seed << 5);

	return 0.5 + 0.2328306e-09 * (seed_in + (int)seed);
}

char ***alloc_grid(long long N, long long block_size) {
	int x, y, z;

	char ***grid = (char ***)malloc(block_size * sizeof(char **));
	if (grid == NULL) {
		cout << "Memory allocation failed" << endl;
		exit(1);
	}

	for (x = 0; x < block_size; x++) {
		grid[x] = (char **)malloc(N * sizeof(char *));
		if (grid[x] == NULL) {
			cout << "Memory allocation failed" << endl;
			exit(1);
		}

		grid[x][0] = (char *)calloc(N * N, sizeof(char));
		if (grid[x][0] == NULL) {
			cout << "Memory allocation failed" << endl;
			exit(1);
		}

		for (y = 1; y < N; y++) grid[x][y] = grid[x][0] + y * N;
	}

	return grid;
}

char ***gen_initial_grid(long long N, float density, int input_seed, int rank, int size) {
	int x, y, z;
	char aux;

	MPI_Request reqs[4];

	long long block_low = BLOCK_LOW(rank, size, N);
	long long block_high = BLOCK_HIGH(rank, size, N);
	long long block_size = block_high - block_low + 1;

	char ***grid = alloc_grid(N, block_size + 2);

	init_r4uni(input_seed);

	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				if (r4_uni() < density) {
					aux = (int)(r4_uni() * N_SPECIES) + 1;
					if (BLOCK_OWNER(x, size, N) == rank) grid[x - block_low + 1][y][z] = aux;
				}

	MPI_Isend(grid[1][0], N * N, MPI_CHAR, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, reqs);
	MPI_Isend(grid[block_size][0], N * N, MPI_CHAR, (rank + 1) % size, 0, MPI_COMM_WORLD, reqs + 1);
	MPI_Irecv(grid[0][0], N * N, MPI_CHAR, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, reqs + 2);
	MPI_Irecv(grid[block_size + 1][0], N * N, MPI_CHAR, (rank + 1) % size, 0, MPI_COMM_WORLD, reqs + 3);

	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

	return grid;
}
