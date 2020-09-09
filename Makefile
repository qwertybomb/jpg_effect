cc = clang++
flags = -march=native -Ofast -Wall -Wextra -std=gnu++20
linker_flags = -lpng

all: jpg.cc
	$(cc) $(flags) jpg.cc -o main $(linker_flags)
