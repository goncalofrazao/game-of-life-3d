#include <iostream>
using namespace std;

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

	long long block_low = BLOCK_LOW(rank, size, N) - 1;
	long long block_high = BLOCK_HIGH(rank, size, N) + 1;
	long long block_size = block_high - block_low + 1;

	long long last_col = block_size - 1;
	long long first_col = 0;

	char ***grid = alloc_grid(N, block_size);

	init_r4uni(input_seed);

	if (rank == 0) {
		block_high++;
		for (x = 1; x < N; x++)
			for (y = 0; y < N; y++)
				for (z = 0; z < N; z++)
					if (r4_uni() < density) {
						aux = (int)(r4_uni() * N_SPECIES) + 1;
						if (x <= block_high) grid[x][y][z] = aux;
					}

		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				if (r4_uni() < density) grid[0][y][z] = (int)(r4_uni() * N_SPECIES) + 1;

	} else if (rank == size - 1) {
		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				if (r4_uni() < density) grid[last_col][y][z] = (int)(r4_uni() * N_SPECIES) + 1;

		for (x = 1; x < N; x++)
			for (y = 0; y < N; y++)
				for (z = 0; z < N; z++)
					if (r4_uni() < density) {
						aux = (int)(r4_uni() * N_SPECIES) + 1;
						if (x >= block_low) grid[x - block_low][y][z] = aux;
					}

	} else {
		for (x = 0; x < N; x++)
			for (y = 0; y < N; y++)
				for (z = 0; z < N; z++)
					if (r4_uni() < density) {
						aux = (int)(r4_uni() * N_SPECIES) + 1;
						if (x >= block_low && x <= block_high) grid[x - block_low][y][z] = aux;
					}
	}

	return grid;
}
