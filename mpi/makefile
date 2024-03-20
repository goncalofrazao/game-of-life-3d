
all: life3d-mpi

life3d-mpi: life3d-mpi.cc grid.cc
	mpic++ -fopenmp -O2 -o life3d-mpi life3d-mpi.cc grid.cc

o3: life3d-mpi.cc grid.cc
	mpic++ -fopenmp -O3 -o life3d-mpi life3d-mpi.cc grid.cc

clean:
	rm -f life3d-mpi