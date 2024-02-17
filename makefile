all: life3d life3d2 life3d3

LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
CPPFLAGS="-I/opt/homebrew/opt/libomp/include"

life3d: life3d.cc grid.cc
	/opt/homebrew/bin/g++-13 -o life3d life3d.cc grid.cc $(LDFLAGS) $(CPPFLAGS) -fopenmp -O3

clean:
	rm -f life3d