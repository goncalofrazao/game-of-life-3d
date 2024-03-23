all: life3d-omp

life3d-omp: life3d-omp.cc grid.cc
	g++ -fopenmp -O2 -o life3d-omp life3d-omp.cc grid.cc

o3: life3d-omp.cc grid.cc
	g++ -fopenmp -O3 -o life3d-omp life3d-omp.cc grid.cc

my-pc: life3d-omp.cc grid.cc
	g++-13 -fopenmp -O3 -o life3d-omp life3d-omp.cc grid.cc

clean:
	rm -f life3d-omp