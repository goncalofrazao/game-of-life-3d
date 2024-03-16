all: life3d

life3d: life3d.cc grid.cc
	g++ -fopenmp -O2 -o life3d life3d.cc grid.cc

clean:
	rm -f life3d