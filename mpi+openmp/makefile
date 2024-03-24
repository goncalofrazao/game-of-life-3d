all: life3d-mpi-omp

life3d-mpi-omp: life3d-mpi-omp.cc grid.cc
	mpic++ -fopenmp -O2 -o life3d-mpi-omp life3d-mpi-omp.cc grid.cc

o3: life3d-mpi-omp.cc grid.cc
	mpic++ -fopenmp -O3 -o life3d-mpi-omp-omp life3d-mpi-omp.cc grid.cc

clean:
	rm -f life3d-mpi-omp