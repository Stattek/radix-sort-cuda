.PHONY: run build debug build_mpi

run:
	@echo "---Sequential Radix Sort---"
	@./build/src/RadixSort 5 3
	@echo "---MPI Radix Sort---"
	@mpiexec -n 1 radix_sort_mpi.out 5 3

build:
	@mkdir -p build/; cd build/; cmake -DCMAKE_BUILD_TYPE=Release ../ -G Ninja; ninja;
	@mpicc -g -Wall -o radix_sort_mpi.out radix_sort_mpi.c -lm

debug:
	@mkdir -p build/; cd build/; cmake -DCMAKE_BUILD_TYPE=Debug ../ -G Ninja; ninja;
