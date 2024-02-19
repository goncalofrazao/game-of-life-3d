all: life3d

life3d: life3d.cc grid.cc
	g++-13 -fopenmp -O3 -o life3d life3d.cc grid.cc

clean:
	rm -f life3d