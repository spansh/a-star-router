.SUFFIXES:
MAKEFLAGS += -r
SHELL := /bin/bash
.DELETE_ON_ERROR:
.PHONY: all test clean


all: AStarRouter

%: %.cpp
	g++ -std=c++17 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -o $@ $^ -lz -lstdc++fs -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu/  -lspatialindex_c -lspatialindex

clean: AStarRouter
	rm -rf $^

production:
	g++ -std=c++17 -O3 -msseregparm -flto -msse -msse2 -msse3 -msse4.2 -minline-all-stringops -mfpmath=sse -march=native -Wall -Wextra -pedantic -o AStarRouter  AStarRouter.cpp -lz -lstdc++fs -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lspatialindex_c -lspatialindex -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -ltcmalloc

profile:
	g++ -std=c++17 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -o AStarRouter  AStarRouter.cpp -lz -lstdc++fs -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu/  -lspatialindex_c -lspatialindex -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -ltcmalloc -Wl,--no-as-needed,-lprofiler,--as-needed
