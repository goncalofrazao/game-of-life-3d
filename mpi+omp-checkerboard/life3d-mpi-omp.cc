#include <mpi.h>
#include <omp.h>

#include <cstring>	// memcpy function
#include <iomanip>	// setprecision function
#include <iostream>

#include "grid.hh"

#define N_SPECIES 9

using namespace std;

// Global variables to store the number of cells of each species in the local grid, the maximum number of cells of each species and the generation in which it was reached
int local_cells[N_SPECIES + 1];
int max_cells[N_SPECIES + 1];
int generation[N_SPECIES + 1];

// MPI data types and requests
MPI_Datatype x_plane, y_plane, z_plane;
MPI_Request send_x[2], send_y[2], send_z[2], recv_x[2], recv_y[2], recv_z[2];
int _up, _down, _left, _right, _front, _back;
int faces_recv;

// Function to calculate the next state of a cell, based on its neighbors states
char next_state(int x, int y, int z, char ***grid) {
	// Initialize variables
	static int ax, ay, az;
	int sum = 0;

	// Array to store the number of neighbors of each species
	int n_species[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	// Loop through the neighbors of the cell
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

	// If current cell is alive
	if ((int)grid[x][y][z]) {
		if (sum > 4 && sum < 14)
			return (char)grid[x][y][z];
		else
			return (char)0;

		// If current cell is dead
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

// Function to print the results of the simulation
void print_result(char ***grid, long long N) {
	for (int i = 1; i <= N_SPECIES; i++) {
		cout << i << " " << max_cells[i] << " " << generation[i] << endl;
	}
}

// Function to initialize the global variables
void init_max_cells() {
	for (int i = 1; i <= N_SPECIES; i++) {
		generation[i] = 0;
		max_cells[i] = 0;
	}
}

// Function to count the number of cells of each species in the 'local' grid (used only for generation 0)
void count_cells(char ***grid, long long local_N[]) {
	#pragma omp parallel for reduction(+:local_cells)
	for (int x = 1; x <= local_N[0]; x++) {
		for (int y = 1; y <= local_N[1]; y++) {
			for (int z = 1; z <= local_N[2]; z++) {
				local_cells[(int)grid[x][y][z]]++;
			}
		}
	}
}

// Function to update the global variables max_cells and generation, if the current generation has more cells of a species than the previous ones
void get_max(int *cells, int gen) {
	for (int i = 1; i <= N_SPECIES; i++) {
		if (cells[i] > max_cells[i]) {
			max_cells[i] = cells[i];
			generation[i] = gen;
		}
	}
}

// Function assynchronously send boundaries and test if they have been received - returns 1 if so
int test_recv_data(MPI_Comm grid_comm, char ***new_grid, long long local_N[]) {
	int flag;

	switch (faces_recv) {
		case 0: 
			// Send x-plane boundaries assynchronously
			faces_recv++;
			MPI_Isend(&new_grid[1][1][1], 1, x_plane, _down, 1, grid_comm, send_x);
			MPI_Isend(&new_grid[local_N[0]][1][1], 1, x_plane, _up, 0, grid_comm, send_x + 1);
			break;

		case 1: 
			// Test if x-plane boundaries have been received
			MPI_Testall(2, recv_x, &flag, MPI_STATUSES_IGNORE);
			
			//If so, send y-plane boundaries assynchronously
			if (flag) {
				faces_recv++;  // faces_recv to 2 if x faces have been received
				MPI_Isend(&new_grid[0][1][1], 1, y_plane, _left, 1, grid_comm, send_y);
				MPI_Isend(&new_grid[0][local_N[1]][1], 1, y_plane, _right, 0, grid_comm, send_y + 1);
			}
			break;

		case 2: 
			// Test if y-plane boundaries have been received
			MPI_Testall(2, recv_y, &flag, MPI_STATUSES_IGNORE);

			// If so, send z-plane boundaries assynchronously
			if (flag) {
				faces_recv++; // faces_recv to 3 if y faces have been received
				MPI_Isend(&new_grid[0][0][1], 1, z_plane, _back, 1, grid_comm, send_z);
				MPI_Isend(&new_grid[0][0][local_N[2]], 1, z_plane, _front, 0, grid_comm, send_z + 1);
				return 1;
			}
			break;

		case 4:
			// Test if all boundaries have been sent and received
			// Works as a barrier to synchronize all nodes
			MPI_Waitall(2, send_x, MPI_STATUSES_IGNORE);
			MPI_Waitall(2, send_y, MPI_STATUSES_IGNORE);
			MPI_Waitall(2, send_z, MPI_STATUSES_IGNORE);
			MPI_Waitall(2, recv_z, MPI_STATUSES_IGNORE);
			return 1;
			break;

		default:
			return 1;
	}
	return 0;
}

// Function to receive the boundaries from the neighbor nodes assynchronously
void recv_all(MPI_Comm grid_comm, char ***new_grid, long long local_N[]) {
	// Receive the x-plane boundaries from the up and down neighbors
	MPI_Irecv(&new_grid[0][1][1], 1, x_plane, _down, 0, grid_comm, recv_x);
	MPI_Irecv(&new_grid[local_N[0] + 1][1][1], 1, x_plane, _up, 1, grid_comm, recv_x + 1);

	// Receive the y-plane boundaries from the left and right neighbors
	MPI_Irecv(&new_grid[0][0][1], 1, y_plane, _right, 0, grid_comm, recv_y);
	MPI_Irecv(&new_grid[0][local_N[1] + 1][1], 1, y_plane, _left, 1, grid_comm, recv_y + 1);

	// Receive the z-plane boundaries from the front and back neighbors
	MPI_Irecv(&new_grid[0][0][0], 1, z_plane, _back, 0, grid_comm, recv_z);
	MPI_Irecv(&new_grid[0][0][local_N[2] + 1], 1, z_plane, _front, 1, grid_comm, recv_z + 1);
}

// Function to simulate the game of life
void simulation(char ***grid, char*** new_grid, long long local_N[], int generations, int rank, MPI_Comm grid_comm) {
	// Declare a temporary grid pointer
	char ***temp;
	
	// Initialize the global variables in root node
	if (!rank) init_max_cells();

	// Count cells in the initial grid, locally and then reducing the results to the root node
	count_cells(grid, local_N);
	MPI_Reduce(local_cells, max_cells, N_SPECIES + 1, MPI_INT, MPI_SUM, 0, grid_comm);

	// Run the simulation for the specified number of generations
	for (int i = 0; i < generations; i++) {
		// Initialize the local cells array (also counts dead cells: N_SPECIES+1)
		int cells[N_SPECIES + 1] = {0};
		faces_recv = 0;

		// Receive the boundaries from the neighbor nodes assynchronously (Irecv)
		recv_all(grid_comm, new_grid, local_N);

		// Calculate the next state of the boundary cells (x = 1, x = local_N[0])
		for (int y = 1; y <= local_N[1]; y++) {
			for (int z = 1; z <= local_N[2]; z++) {
				new_grid[1][y][z] = next_state(1, y, z, grid); // x = 1
				cells[(int)new_grid[1][y][z]]++;
				new_grid[local_N[0]][y][z] = next_state(local_N[0], y, z, grid); // x = local_N[0]
				cells[(int)new_grid[local_N[0]][y][z]]++;
			}
		}
		// Send x plane boundaries
		test_recv_data(grid_comm, new_grid, local_N);

		// Calculate the next state of the boundary cells (y = 1, y = local_N[1])
		#pragma omp parallel for reduction(+:cells)
		for (int x = 2; x < local_N[0]; x++) {
			for (int z = 1; z <= local_N[2]; z++) {
				new_grid[x][1][z] = next_state(x, 1, z, grid);
				cells[(int)new_grid[x][1][z]]++;
				new_grid[x][local_N[1]][z] = next_state(x, local_N[1], z, grid);
				cells[(int)new_grid[x][local_N[1]][z]]++;
			}
		}
		// Test if x plane boundaries have been received and if so send y plane boundaries
		test_recv_data(grid_comm, new_grid, local_N);

		// Calculate the next state of the boundary cells (z = 1, z = local_N[2])
		#pragma omp parallel for reduction(+:cells)
		for (int x = 2; x < local_N[0]; x++) {
			for (int y = 2; y < local_N[1]; y++) {
				new_grid[x][y][1] = next_state(x, y, 1, grid);
				cells[(int)new_grid[x][y][1]]++;
				new_grid[x][y][local_N[2]] = next_state(x, y, local_N[2], grid);
				cells[(int)new_grid[x][y][local_N[2]]]++;
			}
		}
		// Test if y and x plane boundaries have been received and if so send z plane boundaries
		test_recv_data(grid_comm, new_grid, local_N);

		// Calculate the next state for the inner cells of the local grid
		#pragma omp parallel for reduction(+:cells)
		for (int x = 2; x < local_N[0]; x++) {
			for (int y = 2; y < local_N[1]; y++) {
				for (int z = 2; z < local_N[2]; z++) {
					new_grid[x][y][z] = next_state(x, y, z, grid);
					cells[(int)new_grid[x][y][z]]++;
				}
			}
		}
		// Test if z plane boundaries have been received
		test_recv_data(grid_comm, new_grid, local_N);

		// Reduce the local cells array counter to the root node
		MPI_Reduce(cells, local_cells, N_SPECIES + 1, MPI_INT, MPI_SUM, 0, grid_comm);

		// Update the global variables max_cells and generation (only in root node)
		if (!rank) get_max(local_cells, i + 1);

		// Syncronize all nodes before the next generation: make sure boundaries were exchanged
		while (!test_recv_data(grid_comm, new_grid, local_N));
		faces_recv++; // faces_recv to 4 if z faces have been received
		test_recv_data(grid_comm, new_grid, local_N);

		// Swap the pointers of the grid and the new_grid
		temp = grid;
		grid = new_grid;
		new_grid = temp;
	}
}

// Function to run the simulation using the checkerboard method
void mpi_checkerboard_3d(int N, int generations, float density, int seed) {
	// Store the node rank and the total number of nodes
	int rank, p;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Create a 3D grid communicator with periodic boundaries (wraps around)
	int dims[3] = {0, 0, 0};
	int periods[3] = {1, 1, 1};
	MPI_Comm grid_comm;
	MPI_Dims_create(p, 3, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &grid_comm);

	// Get the coordinates of the current node in the 3D grid
	int coords[3];
	MPI_Cart_coords(grid_comm, rank, 3, coords);

	// Calculate the size of the local grid for the current node (3D block)
	long long local_N[3] = {BLOCK_SIZE(coords[0], dims[0], N), BLOCK_SIZE(coords[1], dims[1], N), BLOCK_SIZE(coords[2], dims[2], N)};

	// Compute node neighbors in the 3D grid
	MPI_Cart_shift(grid_comm, 0, 1, &_down, &_up);
	MPI_Cart_shift(grid_comm, 1, 1, &_left, &_right);
	MPI_Cart_shift(grid_comm, 2, 1, &_back, &_front);

	// Declare MPI datatype for the cube planes
	MPI_Type_vector(local_N[1], local_N[2], local_N[2] + 2, MPI_CHAR, &x_plane);
	MPI_Type_vector(local_N[0] + 2, local_N[2], (local_N[1] + 2) * (local_N[2] + 2), MPI_CHAR, &y_plane);
	MPI_Type_vector((local_N[0] + 2) * (local_N[1] + 2), 1, local_N[2] + 2, MPI_CHAR, &z_plane);
	MPI_Type_commit(&x_plane);
	MPI_Type_commit(&y_plane);
	MPI_Type_commit(&z_plane);

	// Generate the initial grid and a 'next generation' grid (with ghost cells)
	char ***grid = gen_initial_grid(N, density, seed, coords, dims);
	char ***new_grid = alloc_grid(local_N[0] + 2, local_N[1] + 2, local_N[2] + 2);

	// Mark the initial time instant before the simulation starts (using MPI_barrier to synchronize all nodes)
	MPI_Barrier(MPI_COMM_WORLD);
	double exec_time = -MPI_Wtime();

	// Run the simulation
	simulation(grid, new_grid, local_N, generations, rank, grid_comm);

	// Mark the final time instant after the simulation ends (using MPI_barrier to synchronize all nodes)
	MPI_Barrier(MPI_COMM_WORLD);
	exec_time += MPI_Wtime();

	// Print the execution time of the simulation to stderr (only by the master node)
	if (!rank) cerr << fixed << setprecision(1) << exec_time << "s\n";

	// Print the results of the simulation (only by the master node)
	if (!rank) print_result(grid, N);

	MPI_Type_free(&x_plane);
	MPI_Type_free(&y_plane);
	MPI_Type_free(&z_plane);
}

/* MAIN PROGRAM */
int main(int argc, char *argv[]) {
	// Check if the input is valid
	if (argc != 5 || atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atof(argv[3]) < 0 || atof(argv[3]) > 1) {
		cerr << "Usage: " << argv[0] << " <generations> <N> <density> <seed>\n";
		return 0;
	}

	// Parse the input
	int generations = atoi(argv[1]);
	long long N = atoi(argv[2]);
	float density = atof(argv[3]);
	int seed = atoi(argv[4]);

	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Set the number of threads for OpenMP
	omp_set_num_threads(4);

	// Run the simulation - checkerboard method
	mpi_checkerboard_3d(N, generations, density, seed);

	// Finalize MPI
	MPI_Finalize();

	// Exit
	return 0;
}
