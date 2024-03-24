#include <cmath>
#include <iostream>
#include <mpi.h>
#include "grid.hh"

#define N_SPECIES 9

using namespace std;

unsigned int seed;

void init_r4uni(int input_seed) { seed = input_seed + 987654321; }

float r4_uni() {
	int seed_in = seed;

	seed ^= (seed << 13);
	seed ^= (seed >> 17);
	seed ^= (seed << 5);

	return 0.5 + 0.2328306e-09 * (seed_in + (int)seed);
}

// Function to allocate memory for the 3D grid, with specific dimensions
char ***alloc_grid(long long x_dim, long long y_dim, long long z_dim) {
	// Initialize variables of iteration
	long long x, y, z;

	//Allocate in memory a single contiguous block of size x_dim * y_dim * z_dim (total number of elements)
	long long total_size = x_dim * y_dim * z_dim;
	char *block = (char *)calloc(total_size, sizeof(char));
	if (block == NULL) {
		cout << "Memory allocation failed" << endl;
		exit(1);
	}

	// Allocate memory for the grid (first dimension)
	char ***grid = (char ***)malloc(x_dim * sizeof(char **));
	if (grid == NULL) {
		cout << "Memory allocation failed" << endl;
		exit(1);
	}

	// Allocate memory for the grid (second dimension)
	for (x = 0; x < x_dim; x++) {
		grid[x] = (char **)malloc(y_dim * sizeof(char *));
		if (grid[x] == NULL) {
			cout << "Memory allocation failed" << endl;
			exit(1);
		}
	}

	// Point each element of the grid to the corresponding element of the block
	for (x = 0; x < x_dim; x++) {
		for (y = 0; y < y_dim; y++) {
			grid[x][y] = block + x * y_dim * z_dim + y * z_dim;
		}
	}

	return grid;
}

// Function to exchange the boundaries between neighbor nodes in the 3D grid
void exchange_boudaries(char ***grid, long long local_N[], MPI_Comm grid_comm) {

	// Compute the neighbors of the current node in the 3D grid
	int up, down, left, right, front, back;
	MPI_Cart_shift(grid_comm, 0, 1, &down, &up);
	MPI_Cart_shift(grid_comm, 1, 1, &left, &right);
	MPI_Cart_shift(grid_comm, 2, 1, &back, &front);

	// Initiliaze variables for the MPI communication
	MPI_Request reqs[4];
	MPI_Datatype x_plane, y_plane, z_plane;

	// Create the MPI vector for the x_plane
	//  2D plane with dimensions (y*z) -> y blocks of z elements each
	//  stride is z + 2 because of ghost cells - see alloc_grid function
	MPI_Type_vector(local_N[1], local_N[2], local_N[2] + 2, MPI_CHAR, &x_plane);

	// Create the MPI strided vector for the y_plane
	// plane with dimensions (x*z), but we want to send vertices received in x_plane exchange, so (x+2) blocks of z elements each
	// stride is (z+2)*(y+2) - see alloc_grid function
	MPI_Type_vector(local_N[0] + 2, local_N[2], (local_N[1] + 2) * (local_N[2] + 2), MPI_CHAR, &y_plane);

	// Create the MPI strided vector for the z_plane
	//  2D plane with dimensions (x*y), but we want to send vertices received in x_plane and y_plane exchange, so (x+2)*(y+2) blocks with one element each
	//  stride is z+2 - see alloc_grid function
	MPI_Type_vector((local_N[0] + 2) * (local_N[1] + 2), 1, local_N[2] + 2, MPI_CHAR, &z_plane);

	// Commit the MPI datatypes
	MPI_Type_commit(&x_plane);
	MPI_Type_commit(&y_plane);
	MPI_Type_commit(&z_plane);

	// Exchange x_plane boudaries with up and down neighbors (send and receive 2 planes), and wait for the communication to finish
	MPI_Isend(&grid[1][1][1], 1, x_plane, down, 1, grid_comm, reqs);
	MPI_Isend(&grid[local_N[0]][1][1], 1, x_plane, up, 0, grid_comm, reqs + 1);
	MPI_Irecv(&grid[0][1][1], 1, x_plane, down, 0, grid_comm, reqs + 2);
	MPI_Irecv(&grid[local_N[0] + 1][1][1], 1, x_plane, up, 1, grid_comm, reqs + 3);
	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

	// Exchange y_plane boudaries with left and right neighbors (send and receive 2 planes), and wait for the communication to finish
	MPI_Isend(&grid[0][local_N[1]][1], 1, y_plane, right, 0, grid_comm, reqs);
	MPI_Isend(&grid[0][1][1], 1, y_plane, right, 1, grid_comm, reqs + 1);
	MPI_Irecv(&grid[0][0][1], 1, y_plane, left, 0, grid_comm, reqs + 2);
	MPI_Irecv(&grid[0][local_N[1] + 1][1], 1, y_plane, left, 1, grid_comm, reqs + 3);
	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

	// Exchange z_plane boudaries with front and back neighbors (send and receive 2 planes), and wait for the communication to finish
	MPI_Isend(&grid[0][0][local_N[2]], 1, z_plane, front, 0, grid_comm, reqs);
	MPI_Isend(&grid[0][0][1], 1, z_plane, back, 1, grid_comm, reqs + 1);
	MPI_Irecv(&grid[0][0][0], 1, z_plane, back, 0, grid_comm, reqs + 2);
	MPI_Irecv(&grid[0][0][local_N[2] + 1], 1, z_plane, front, 1, grid_comm, reqs + 3);
	MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

	// Free the MPI datatypes
	MPI_Type_free(&x_plane);
	MPI_Type_free(&y_plane);
	MPI_Type_free(&z_plane);
}

// Function to generate the initial grid for the simulation
char ***gen_initial_grid(long long N, float density, int input_seed, int coords[], int dims[]) {

	// Initialize variables 
	int x, y, z;
	char aux;

	// Calculate the size of the local grid for the current node (3D block)
	long long local_N[3] = {BLOCK_SIZE(coords[0], dims[0], N), BLOCK_SIZE(coords[1], dims[1], N), BLOCK_SIZE(coords[2], dims[2], N)};

	// Calculate the lower bounds of the local grid for the current node (3D block)
	long long x_low = BLOCK_LOW(coords[0], dims[0], N);
	long long y_low = BLOCK_LOW(coords[1], dims[1], N);
	long long z_low = BLOCK_LOW(coords[2], dims[2], N);

	// Allocate memory for the local grid with 2 extra cells in each dimension for the boundaries
	char ***grid = alloc_grid(local_N[0] + 2, local_N[1] + 2, local_N[2] + 2);

	init_r4uni(input_seed);

	// Fill the grid with random values based on the density
	// Although this process is the same for every node, each one only stores in memory its corresponding block
	for (x = 0; x < N; x++)
		for (y = 0; y < N; y++)
			for (z = 0; z < N; z++)
				if (r4_uni() < density) {
					aux = (int)(r4_uni() * N_SPECIES) + 1;
					if (CUBE_OWNER(x, y, z, dims, coords, N)) {
						grid[x - x_low + 1][y - y_low + 1][z - z_low + 1] = aux;
					}
				}

	// Create a 3D grid communicator with periodic boundaries (wraps around)
	int periods[3] = {1, 1, 1};
	MPI_Comm grid_comm;
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &grid_comm);

	// Exchange boundaries between nodes (to ensure that each node has the needed data to start the simulation)
	exchange_boudaries(grid, local_N, grid_comm);

	return grid;
}
