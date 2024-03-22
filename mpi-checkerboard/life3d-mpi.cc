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

char next_state(int x, int y, int z, char ***grid) {
	static int ax, ay, az;
	int sum = 0;
	int n_species[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	for (int i = -1; i <= 1; i++) {
		ax = x + i;
		for (int j = -1; j <= 1; j++) {
			ay = y + j;
			for (int k = -1; k <= 1; k++) {
				if (i == 0 && j == 0 && k == 0) continue;
				az = z + k;
				sum += ((int)grid[ax][ay][az] != 0);
				n_species[(int)grid[ax][ay][az]]++;
			}
		}
	}

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

void count_cells(char ***grid, long long local_N[]) {
	for (int x = 1; x <= local_N[0]; x++) {
		for (int y = 1; y <= local_N[1]; y++) {
			for (int z = 1; z <= local_N[2]; z++) {
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

void simulation(char ***grid, long long local_N[], int generations, int rank, MPI_Comm grid_comm) {
	char ***new_grid = alloc_grid(local_N[0] + 2, local_N[1] + 2, local_N[2] + 2), ***temp;

	MPI_Request reqs[4];

	if (!rank) init_max_cells();
	count_cells(grid, local_N);
	MPI_Reduce(local_cells, max_cells, N_SPECIES + 1, MPI_INT, MPI_SUM, 0, grid_comm);

	for (int i = 0; i < generations; i++) {
		int cells[N_SPECIES + 1] = {0};

		for (int x = 1; x <= local_N[0]; x++) {
			for (int y = 1; y <= local_N[1]; y++) {
				for (int z = 1; z <= local_N[2]; z++) {
					new_grid[x][y][z] = next_state(x, y, z, grid);
					cells[(int)new_grid[x][y][z]]++;
				}
			}
		}

		MPI_Reduce(cells, local_cells, N_SPECIES + 1, MPI_INT, MPI_SUM, 0, grid_comm);
		if (!rank) get_max(local_cells, i + 1);

		exchange_boudaries(new_grid, local_N, grid_comm);

		temp = grid;
		grid = new_grid;
		new_grid = temp;
	}
}

void print_grid(char ***grid, long long N, int id[], int p[]) {
	MPI_Request reqs[4];
	long long block_sizes[3] = {BLOCK_SIZE(id[0], p[0], N), BLOCK_SIZE(id[1], p[1], N), BLOCK_SIZE(id[2], p[2], N)};

	int l[3] = {p[0] - 1, p[1] - 1, p[2] - 1};
	int global_size = DIMS_TO_INDEX(l, p) + 1;

	MPI_Barrier(MPI_COMM_WORLD);

	int a = 1;
	for (int i = 0; i < global_size; i++) {
		if (DIMS_TO_INDEX(id, p) == i) {
			cout << "Rank " << i << endl;
			for (int x = 1 - a; x <= block_sizes[0] + a; x++) {
				cout << "x = " << x + BLOCK_LOW(id[0], p[0], N) - 1 << endl;
				for (int y = 1 - a; y <= block_sizes[1] + a; y++) {
					for (int z = 1 - a; z <= block_sizes[2] + a; z++) {
						cout << (int)grid[x][y][z] << " ";
					}
					cout << endl;
				}
				cout << endl;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void mpi_checkerboard_3d(int N, int generations, float density, int seed) {
	int rank, p;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	int dims[3] = {0, 0, 0};
	MPI_Dims_create(p, 3, dims);

	int periods[3] = {1, 1, 1};
	MPI_Comm grid_comm;
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &grid_comm);

	int coords[3];
	MPI_Cart_coords(grid_comm, rank, 3, coords);

	long long local_N[3] = {BLOCK_SIZE(coords[0], dims[0], N), BLOCK_SIZE(coords[1], dims[1], N), BLOCK_SIZE(coords[2], dims[2], N)};

	int up, down, left, right, front, back;
	MPI_Cart_shift(grid_comm, 0, 1, &down, &up);
	MPI_Cart_shift(grid_comm, 1, 1, &left, &right);
	MPI_Cart_shift(grid_comm, 2, 1, &back, &front);

	double exec_time;

	char ***grid = gen_initial_grid(N, density, seed, coords, dims);

	MPI_Barrier(MPI_COMM_WORLD);
	exec_time = -MPI_Wtime();

	simulation(grid, local_N, generations, rank, grid_comm);
	// print_grid(grid, N, coords, dims);

	MPI_Barrier(MPI_COMM_WORLD);
	exec_time += MPI_Wtime();

	if (!rank) cerr << fixed << setprecision(1) << exec_time << "s\n";

	if (!rank) print_result(grid, N);
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

	mpi_checkerboard_3d(N, generations, density, seed);

	MPI_Finalize();
}
