#ifndef _GRID_HH_
#define _GRID_HH_

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#define BLOCK_OWNER(index, p, n) (((p) * ((index) + 1) - 1) / (n))

#define CUBE_OWNER(x, y, z, p, id, n) (BLOCK_OWNER(x, p[0], n) == id[0] && BLOCK_OWNER(y, p[1], n) == id[1] && BLOCK_OWNER(z, p[2], n) == id[2])

#define INDEX_TO_DIMS(i, p) \
	{ i / (p[1] * p[2]), (i / p[2]) % p[1], i % p[2] }
#define DIMS_TO_INDEX(id, p) (id[0] * p[1] * p[2] + id[1] * p[2] + id[2])

#define NEXT_BLOCK(id, p) (((id) + 1) % p)
#define PREV_BLOCK(id, p) (((id) - 1 + p) % p)

char ***alloc_grid(long long x_dim, long long y_dim, long long z_dim);
void exchange_boudaries(char ***grid, long long local_N[], MPI_Comm grid_comm);
char ***gen_initial_grid(long long N, float density, int input_seed, int coords[], int dims[]);

#endif	// _GRID_HH_
