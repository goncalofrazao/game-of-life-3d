#ifndef _GRID_HH_
#define _GRID_HH_

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#define BLOCK_OWNER(index, p, n) (((p) * ((index) + 1) - 1) / (n))

char ***alloc_grid(long long N, long long block_size);
char ***gen_initial_grid(long long N, float density, int input_seed, int rank, int size);

#endif	// _GRID_HH_
