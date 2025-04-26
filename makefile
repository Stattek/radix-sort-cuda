.PHONY: build run test

build:
	mpic++ -Wall -g radix_sort_mpi.cpp -o rsort_mpi.out -fopenmp

run:
	mpirun -np 1 ./rsort_mpi.out input.txt 1

define fail
	(echo "FAIL";\
	exit 1;)
endef

define pass
	echo "PASS";
endef

test:
	mpirun -np 1 ./rsort_mpi.out test_data/input.txt 1 || $(fail)
	mpirun -np 2 ./rsort_mpi.out test_data/input.txt 1 || $(fail)
	mpirun -np 3 ./rsort_mpi.out test_data/input.txt 1 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/input.txt 1 || $(fail)
	mpirun -np 5 ./rsort_mpi.out test_data/input.txt 1 || $(fail)
	mpirun -np 6 ./rsort_mpi.out test_data/input.txt 1 || $(fail)
	mpirun -np 7 ./rsort_mpi.out test_data/input.txt 1 || $(fail)
	mpirun -np 8 ./rsort_mpi.out test_data/input.txt 1 || $(fail)
	mpirun -np 1 ./rsort_mpi.out test_data/input.txt 2 || $(fail)
	mpirun -np 2 ./rsort_mpi.out test_data/input.txt 2 || $(fail)
	mpirun -np 3 ./rsort_mpi.out test_data/input.txt 2 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/input.txt 2 || $(fail)
	mpirun -np 5 ./rsort_mpi.out test_data/input.txt 2 || $(fail)
	mpirun -np 6 ./rsort_mpi.out test_data/input.txt 2 || $(fail)
	mpirun -np 7 ./rsort_mpi.out test_data/input.txt 2 || $(fail)
	mpirun -np 8 ./rsort_mpi.out test_data/input.txt 2 || $(fail)
	mpirun -np 1 ./rsort_mpi.out test_data/input.txt 3 || $(fail)
	mpirun -np 2 ./rsort_mpi.out test_data/input.txt 3 || $(fail)
	mpirun -np 3 ./rsort_mpi.out test_data/input.txt 3 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/input.txt 3 || $(fail)
	mpirun -np 5 ./rsort_mpi.out test_data/input.txt 3 || $(fail)
	mpirun -np 6 ./rsort_mpi.out test_data/input.txt 3 || $(fail)
	mpirun -np 7 ./rsort_mpi.out test_data/input.txt 3 || $(fail)
	mpirun -np 8 ./rsort_mpi.out test_data/input.txt 3 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/random_data_1000000_values_1_digits.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/random_data_1000000_values_2_digits.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/random_data_1000000_values_3_digits.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/random_data_1000000_values_7_digits.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/random_data_1000000_values_8_digits.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/random_data_1000000_values_9_digits.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/random_data_100000_values_8_digits.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/random_data_100000_values_9_digits.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/rust_test_input.txt 4 || $(fail)
	mpirun -np 4 ./rsort_mpi.out test_data/rust_test_input2.txt 4 || $(fail)
	$(pass)
