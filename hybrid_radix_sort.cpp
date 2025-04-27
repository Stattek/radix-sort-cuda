/*
 * Compile: mpic++ -Wall hybrid_radix_sort.cpp -o hybrid_radix_sort.out -fopenmp
 * Run: mpirun -np <num_processes> ./hybrid_radix_sort.out <input_file_name> <omp_num_threads>
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#define ARRAY_PRINT_THRESHOLD 20

#define NUM_BASE 256
#define COUNT_ARRAY_SIZE NUM_BASE // the count array will always hold the same number of values as the number of digits
#define INITIAL_ARRAY_SIZE 20

/**
 * @brief Finds the maxmimum value in an array at a digit and outputs it.
 *
 * @param array The array to search.
 * @param arrayLen The length of the array.
 * @param digit The digit to find the maximum value for.
 * @param output The output maximum value.
 *
 * @returns `true` on failure, `false` on success.
 */
static bool getMax(const uint *array, const uint arrayLen, uint *output)
{
    if (!array || !output)
    {
        return true; // fail
    }

    uint maxValue = array[0];
// find the maximum
// NOTE: why was the OMP reduction not working here?
#pragma omp parallel for
    for (uint i = 1; i < arrayLen; i++)
    {
#pragma omp critical
        if (maxValue < array[i])
        {
            maxValue = array[i];
        }
    }

    *output = maxValue;
    return false; // success
}

static unsigned long long myPow(uint first, uint exponent)
{
    unsigned long long sum = 1;
    for (uint i = 0; i < exponent; i++)
    {
        sum *= first;
    }

    return sum;
}

/**
 * @brief Gets the number of digits that this value has.
 *
 * @param value The value to find the number of digits of.
 * @return The number of digits.
 */
static uint getNumDigits(uint value)
{
    uint numDigits = 1; // we should start at one
    value /= NUM_BASE;
    while (value > 0)
    {
        numDigits++;
        value /= NUM_BASE;
    }

    return numDigits;
}

/**
 * @brief Prints the array.
 *
 * @param name The name of the array.
 * @param array The array to print.
 * @param arrayLen The length of the array.
 */
static void printArray(const char *name, uint *array, uint arrayLen)
{
    printf("\nArray %s, length %d:", name, arrayLen);
    if (arrayLen <= ARRAY_PRINT_THRESHOLD)
    {
        printf(" [");
        for (uint i = 0; i < arrayLen; i++)
        {
            printf("%d", array[i]);

            if (i != arrayLen - 1)
            {
                printf(", ");
            }
        }
        printf("]\n");
    }
    else
    {
        printf(" Array above printing threshold %d\n", ARRAY_PRINT_THRESHOLD);
    }
}

/**
 * @brief Reads an integer array from a file.
 *
 * @param fileName The name of the file to read from.
 * @param outputNumElements The length of the output array.
 * @note The output array is allocated on the heap and MUST be deallocated by the user.
 *
 * @returns True on failure, false on success.
 */
static uint *readIntArrayFromFile(const char *fileName, uint &outputNumElements)
{
    if (!fileName)
    {
        // bad argument
        return NULL;
    }

    std::vector<uint> output;

    // since this can fail when reading from input
    try
    {
        std::ifstream inputFile(fileName, std::ios_base::in);
        int curInt;
        while (inputFile >> curInt)
        {
            output.push_back((uint)curInt);
        }
    }
    catch (...)
    {
        return NULL; // fail
    }

    if (output.empty())
    {
        // either no integers or a bad file name
        return NULL;
    }

    uint *outputPointer = new uint[output.size()];
    for (uint i = 0; i < (uint)output.size(); i++)
    {
        outputPointer[i] = output[i];
    }
    outputNumElements = output.size();

    return outputPointer; // success
}

/**
 * @brief Checks if the array is sorted.
 *
 * @param array The array to check.
 * @param arrayLen The array length.
 *
 * @returns `true` if the array is sorted, `false` otherwise.
 */
static bool isSorted(int *array, uint arrayLen)
{
    bool output = true;
#pragma omp parallel for
    for (uint i = 1; i < arrayLen; i++)
    {
        if (array[i - 1] > array[i])
        {
            output = false;
        }
    }

    return output;
}

/**
 * @brief Computes offsets from the countmatrix.
 *
 * @param countMatrix The count matrix.
 * @param numSections The number of sections that the array has been split into.
 * @param numValues The number of values per section.
 * @param offsetMatrix The offset matrix.
 */
static void computeOffsets(uint **countMatrix, uint numSections, uint numValues, uint **offsetMatrix)
{
// Initialize the offset matrix
#pragma omp parallel for
    for (uint m = 0; m < numSections; m++)
    {
        for (uint n = 0; n < numValues; n++)
        {
            offsetMatrix[m][n] = 0;

            // First term: Sum of counts for all values < n in all sections
            for (uint i = 0; i < n; i++)
            {
                for (uint j = 0; j < numSections; j++)
                {
                    offsetMatrix[m][n] += countMatrix[j][i];
                }
            }

            // Second term: Sum of counts for value n in all sections < m
            for (uint j = 0; j < m; j++)
            {
                offsetMatrix[m][n] += countMatrix[j][n];
            }
#if 0
            // Debug print to check offsetMatrix values
            // printf("OffsetMatrix[%d][%d] = %d\n", m, n, offsetMatrix[m][n]);
#endif
        }
    }
}

/**
 * @brief Computes offsets for the local part of the array using the global offset matrix.
 *
 * @param localArray The local array to read from
 * @param localArraySize The size of the local array.
 * @param offsetMatrix The offsetMatrix.
 * @param numValues The number of values to save in the local offsets.
 * @param rank The current processor rank.
 * @param offsets The offsets array to modify.
 * @param digit The current digit to compute for.
 */
static void computeLocalOffsets(const uint *localArray, const uint localArraySize, uint **offsetMatrix,
                                const uint numValues, const uint rank, uint *offsets, const uint digit)
{
    // Create a local copy of the offsets for the current process to track updates
    uint *localOffsets = new uint[numValues];
#pragma omp parallel for
    for (uint i = 0; i < numValues; i++)
    {
        localOffsets[i] = offsetMatrix[rank][i];
    }

    // Assign indexes to each value in the local array
    for (uint i = 0; i < localArraySize; i++)
    {
        uint value = (localArray[i] / digit) % NUM_BASE;
        offsets[i] = localOffsets[value]; // Assign the current offset for the value
        localOffsets[value]++;            // Increment the offset for the next occurrence

#if 0
        // Debug print to check offsets
        printf("Process %d: Value %d assigned offset %d\n", rank, value, offsets[i]);
#endif
    }

    delete[] localOffsets;
    localOffsets = NULL;
}

/**
 * @brief Moves values from the local array to their correct positions in the global array
 * using the computed offsets.
 *
 * @param localArray The local array to read from.
 * @param localArraySize The length of the local array.
 * @param offsets The offsets array.
 * @param globalArray The global array to place values into their correct positions.
 */
static void placeValuesFromOffset(uint *localArray, uint localArraySize, uint *offsets, uint *globalArray)
{
#pragma omp parallel for
    for (uint i = 0; i < localArraySize; i++)
    {
        globalArray[offsets[i]] = localArray[i];
    }
}

/**
 * @brief Updates the count matrix for the local array.
 *
 * @param countMatrix The countmatrix to modify.
 * @param localArray The local array to read from.
 * @param localArraySize The local array size.
 * @param digit The current digit to update the count matrix on.
 */
static void updateCountMatrix(uint *countMatrix, const uint *localArray, const uint localArraySize, const uint digit)
{
    // Reset the count matrix for the current process
    (void)memset(countMatrix, 0, sizeof(uint) * COUNT_ARRAY_SIZE);

    // Count the occurrences of each digit at the current place value in the local array
    // NOTE: parallelizing this loop causes great slowdown
    for (uint i = 0; i < localArraySize; i++)
    {
        uint digitValue = (localArray[i] / digit) % NUM_BASE;
        // NOTE: this makes it really slow
        {
            countMatrix[digitValue]++;
        }
    }
}

static void flipSignBits(int *array, uint arrayLength)
{
    // create mask for flipping the sign bit
    uint mask = -1;
    mask -= 1;
    mask = (mask >> 1); // shift one to the right so we flip only the sign bit
    mask = ~mask;

    // flip all the sign bits
#pragma omp parallel for
    for (uint i = 0; i < arrayLength; i++)
    {
        array[i] ^= mask;
    }
}

int main(int argc, char *argv[])
{
    // initialize MPI
    MPI_Init(&argc, &argv);

    // setup ranks and nproc
    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm comm = MPI_COMM_WORLD;

    /* shared variables */
    uint inputArraySize = 0, maxDigit = 0;
    unsigned long long maxPossibleValue = 0;
    uint *inputArray = NULL;
    uint *outputArray = NULL;

    // temp arrays for saving results
    int *tempSendRecvCounts = NULL;
    int *tempDisplacements = NULL;

    /* radix sort setup */
    if (rank == 0)
    {
        if (argc != 3)
        {
            printf("Usage: %s <input_file_name> <omp_num_threads>\n", argv[0]);
            MPI_Abort(comm, 1);
        }

        const char *inputFileName = argv[1];

#ifdef _OPENMP // in case the compiler doesn't have openmp
        int omp_threads = atoi(argv[2]);
        omp_set_num_threads(omp_threads);
#endif

        inputArray = readIntArrayFromFile(inputFileName, inputArraySize);
        if (!inputArray)
        {
            fprintf(stderr, "Could not read array from file %s. File may be empty or not exist\n", inputFileName);
            MPI_Abort(comm, 1);
        }

        // initial info
        printf("\nInput file = \"%s\" OpenMPI numProcessors = %d", inputFileName, nproc);
#ifdef _OPENMP
        printf(", OpenMP numThreads = %d\n", omp_threads);
#else
        printf("\n");
#endif

        // print the array
        printArray("Initial", inputArray, inputArraySize);

        // flip bits, then do the rest of the setup
        flipSignBits((int *)inputArray, inputArraySize);

        outputArray = new uint[inputArraySize];
        (void)memset(outputArray, 0, sizeof(uint) * inputArraySize);

        uint maxValue;
        if (getMax(inputArray, inputArraySize, &maxValue))
        {
            fprintf(stderr, "Could not get the maximum value in the input array\n");
            MPI_Abort(comm, 1);
        }

        // find out the number of digits in this maximum value
        maxDigit = getNumDigits(maxValue);

        // FIXME: evil max possible value, don't like this
        maxPossibleValue = myPow(NUM_BASE, maxDigit);
    }

    MPI_Bcast(&inputArraySize, 1, MPI_UNSIGNED, 0, comm);
    MPI_Bcast(&maxPossibleValue, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);

    // the base size of our local array
    uint baseLocalArraySize = inputArraySize / nproc;
    uint remainder = inputArraySize % nproc; // remainder for the rest of the array

    // since remainder is less than nproc, if the rank is less than the remainder, add 1 to it
    uint localArraySize = (rank < (int)remainder) ? baseLocalArraySize + 1 : baseLocalArraySize;

    // Create a local count array to store the count of each digit for the current process
    uint *localCountArray = new uint[COUNT_ARRAY_SIZE];

    // Create a local offset array to store the offsets for the current process
    uint *localOffsetArray = new uint[localArraySize];

    // local array to sort
    uint *localArray = new uint[localArraySize];

    // TODO: use 2d representation of matrices instead of doing this weird stuff to work with MPI
    // Output Array
    // Allocate countMatrix and offsetMatrix as contiguous blocks
    uint *flatCountMatrix = new uint[nproc * COUNT_ARRAY_SIZE];
    uint **countMatrix = new uint *[nproc];
#pragma omp parallel for
    for (int i = 0; i < nproc; i++)
    {
        countMatrix[i] = &flatCountMatrix[i * COUNT_ARRAY_SIZE];
    }

    uint *flatOffsetMatrix = new uint[nproc * COUNT_ARRAY_SIZE];
    uint **offsetMatrix = new uint *[nproc];
#pragma omp parallel for
    for (int i = 0; i < nproc; i++)
    {
        offsetMatrix[i] = &flatOffsetMatrix[i * COUNT_ARRAY_SIZE];
    }

    // set up the temporary arrays to be able to do scatters and gathers
    // set up the temp arrays
    tempSendRecvCounts = new int[nproc];
    tempDisplacements = new int[nproc];

    int curDisplacement = 0; // the current displacement
    // NOTE: cannot parallelize with OMP due to a dependency on the last iteration
    // iterate through process ranks
    for (int i = 0; i < nproc; i++)
    {
        // set up the temp receive counts array and calculate their values if
        // there is a remainder
        if (i < (int)remainder) // current process rank is less than remainder
        {
            tempSendRecvCounts[i] = (int)baseLocalArraySize + 1;
        }
        else
        {
            tempSendRecvCounts[i] = baseLocalArraySize;
        }

        // add the current displacement to the temp displacement array
        tempDisplacements[i] = curDisplacement;
        curDisplacement += tempSendRecvCounts[i];
    }

    double startTime = MPI_Wtime();

    // NOTE: do not parallelize anything with MP that calls MPI functions
    for (unsigned long long digit = 1; digit < maxPossibleValue; digit *= NUM_BASE)
    {
        // scatter the input array into local arrays
        MPI_Scatterv(inputArray, tempSendRecvCounts, tempDisplacements, MPI_UNSIGNED, localArray, (int)localArraySize, MPI_UNSIGNED, 0, comm);

        // update the local count array as the matrix
        updateCountMatrix(localCountArray, localArray, localArraySize, digit);

        // Gather localCountArray into flatCountMatrix
        MPI_Gather(localCountArray, COUNT_ARRAY_SIZE, MPI_UNSIGNED,
                   flatCountMatrix, COUNT_ARRAY_SIZE, MPI_UNSIGNED, 0, comm);

        MPI_Bcast(flatCountMatrix, nproc * COUNT_ARRAY_SIZE, MPI_UNSIGNED, 0, comm);

        // compute offsets
        computeOffsets(countMatrix, nproc, COUNT_ARRAY_SIZE, offsetMatrix);

        // gather offsetMatrix
        MPI_Gather(offsetMatrix[rank], COUNT_ARRAY_SIZE, MPI_UNSIGNED,
                   flatOffsetMatrix, COUNT_ARRAY_SIZE, MPI_UNSIGNED, 0, comm);
        MPI_Bcast(flatOffsetMatrix, nproc * COUNT_ARRAY_SIZE, MPI_UNSIGNED, 0, comm);

        // compute local offsets
        computeLocalOffsets(localArray, localArraySize, offsetMatrix,
                            COUNT_ARRAY_SIZE, rank, localOffsetArray, digit);

        uint *tempOffsetArray = new uint[inputArraySize];
        MPI_Gatherv(localOffsetArray, localArraySize, MPI_UNSIGNED, tempOffsetArray,
                    tempSendRecvCounts, tempDisplacements, MPI_UNSIGNED, 0, comm);

        // do the move values
        if (rank == 0)
        {
            placeValuesFromOffset(inputArray, inputArraySize, tempOffsetArray, outputArray);
            // Copy the output array back to the input array for the next iteration
            for (uint i = 0; i < inputArraySize; i++)
            {
                inputArray[i] = outputArray[i];
            }
        }

        delete[] tempOffsetArray;
        tempOffsetArray = NULL;
    }

    MPI_Barrier(comm); // Ensure all processes are done

    // save time
    double elapsedTime = MPI_Wtime() - startTime;

    if (rank == 0)
    {

        flipSignBits((int *)outputArray, inputArraySize);

        printArray("Final", outputArray, inputArraySize);

        if (isSorted((int *)outputArray, inputArraySize))
        {
            printf("\n\nThe array is sorted in %lf second(s).\n", elapsedTime);
        }
        else
        {
            printf("\n\nThe array is not sorted.\n");
            MPI_Abort(comm, 1);
        }

        // delete values for rank 0
        delete[] inputArray;
        inputArray = NULL;
        delete[] outputArray;
        outputArray = NULL;
    }

    // delete shared values
    delete[] tempSendRecvCounts;
    tempSendRecvCounts = NULL;
    delete[] tempDisplacements;
    tempDisplacements = NULL;
    delete[] localCountArray;
    localCountArray = NULL;
    delete[] localOffsetArray;
    localOffsetArray = NULL;
    delete[] flatCountMatrix;
    flatCountMatrix = NULL;
    delete[] countMatrix;
    countMatrix = NULL;
    delete[] flatOffsetMatrix;
    flatOffsetMatrix = NULL;
    delete[] offsetMatrix;
    offsetMatrix = NULL;
    delete[] localArray;
    localArray = NULL;

    MPI_Finalize();
    return 0;
}