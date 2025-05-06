# README

## OpenMP and OpenMPI

### Overview

Some important functions are:

- `readIntFromFile()`
  - Reads integers from file and saves them to an array.
- `computeOffsets()`
  - Computes the offsets from the count matrix.
- `computeLocalOffsets()`
  - Computes offsets for the local part of the array using the global offset matrix.
- `placeValuesFromOffset()`
  - Moves values from the local array to their correct positions in the global array using the computed offsets.
- `updateCountMatrix()`
  - Updates the count matrix for the local array by counting occurrences of values.
- `flipSignBits()`
  - Flips the sign bits of all values in the array. Allows for sorting negative numbers as if they were positive.
- `main()`
  - Where the main logic for the Radix Sort occurs. This includes all of the MPI calls.

### Building and Running

To compile and run the OpenMP and OpenMPI version of the program, `hybrid_radix_sort.cpp`, the makefile can be used.

To build, use the command:

```sh
make build
```

To run locally, use the command:

```sh
make run
```

## CUDA

The CUDA version has much of the general logic in place for performing Radix Sort, but,
unfortunately, suffers from bugs that results in incorrect sort order.

NOTE: the CUDA version has many debug prints which are disabled with `#if` blocks.
To enable these, set the define `DO_CUDA_DEBUG` at the top of the file to `1`.

### Overview

Some important functions are:

- `readIntFromFile()`
  - Reads integers from file and saves them to an array.
- `computeOffsets()`
  - Computes the offsets from the count matrix.
- `computeLocalOffsets()`
  - Computes offsets for the local part of the array using the global offset matrix.
- `placeValuesFromOffset()`
  - Moves values from the local array to their correct positions in the global array using the computed offsets.
- `updateCountMatrix()`
  - Updates the count matrix for the local array by counting occurrences of values.
- `flipSignBits()`
  - Flips the sign bits of all values in the array. Allows for sorting negative numbers as if they were positive.
- `sumOffsetTable()`
  - Finds the sum for a CUDA block by iterating through all values and saving the sum in a global offset array.
- `prefixSum()`
  - Finds a partial prefix sum for the current iteration in the parallel scan algorithm.
- `parallelScan()`
  - Performs the logic for the parallel scan algorithm.
- `reorderElements()`
  - Reorders elements into their final, sorted positions based on the local offsets and global offsets.
- `main()`
  - Where the main logic for the Radix Sort occurs. This includes much of the CUDA function calls.

### Building and Running

To compile and run the CUDA version of the program, `cuda_radix_sort.cu`, the makefile can be used.

To build, use the command:

```sh
make cuda
```

To run locally, use the command:

```sh
make runcuda
```
