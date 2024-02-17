all: life3d

life3d: life3d.cc grid.cc
	/opt/homebrew/bin/g++-13 -o life3d life3d.cc grid.cc -fopenmp -O3

clean:
	rm -f life3d