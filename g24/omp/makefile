all: life3d-omp

life3d-omp: life3d-omp.cc grid.cc
	g++ -fopenmp -O2 -o life3d-omp life3d-omp.cc grid.cc

clean:
	rm -f life3d-omp