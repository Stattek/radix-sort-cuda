/*
 * Compile: mpic++ -Wall radix_sort_mpi.cpp -o rsort_mpi.out -fopenmp
 * Run: mpirun -np <num_processes> ./rsort_mpi.out <sizeOfArray> <#Digits> <Omp_threads>
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#define ARRAY_PRINT_THRESHOLD 20

#define NUM_DIGITS 10
#define COUNT_ARRAY_SIZE NUM_DIGITS // the count array will always hold the same number of values as the number of digits
#define INITIAL_ARRAY_SIZE 20

/**
 * @brief Reads an integer array from a file.
 *
 * @param fileName The name of the file to read from.
 * @param output The output array.
 * @param outputLength The length of the output array.
 * @note The output array is allocated on the heap and MUST be deallocated by the user.
 *
 * @returns True on failure, false on success.
 */
static inline bool
readIntArrayFromFile(const char *fileName, int *output, int &outputLength)
{
    if (!fileName)
    {
        // bad argument
        return true;
    }
    output = (int *)malloc(sizeof(int) * INITIAL_ARRAY_SIZE);
    if (!output)
    {
        // could not allocate
        return true; // fail
    }
    outputLength = INITIAL_ARRAY_SIZE; // contiguous size on heap
    int outputNumElements = 0;

    // since this can fail when reading from input
    try
    {
        std::ifstream inputFile(fileName, std::ios_base::in);
        int curInt;
        while (inputFile >> curInt)
        {
            if (outputNumElements >= outputLength)
            {
                // we have to resize the array
                outputLength *= 2; // double the length of the array
                output = (int *)realloc(output, sizeof(int) * outputLength);
                if (!output)
                {
                    // could not allocate
                    return true;
                }
            }

            output[outputLength] = curInt;
            outputNumElements++;
        }
    }
    catch (...)
    {
        return true; // fail
    }

    return false; // success
}

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
static bool getMax(const int *array, const int arrayLen, int *output)
{
    if (!array || !output)
    {
        return true; // fail
    }

    int maxValue = 0;
    // find the maximum

#pragma omp parallel for reduction(max : maxValue)
    for (int i = 0; i < arrayLen; i++)
    {
        maxValue = array[i];
    }

    *output = maxValue;
    return false; // success
}

/**
 * @brief Gets the number of digits that this value has.
 *
 * @param value The value to find the number of digits of.
 * @return The number of digits.
 */
static int getNumDigits(int value)
{
    int numDigits = 0;
    value /= NUM_DIGITS;
    while (value > 0)
    {
        numDigits++;
        value /= NUM_DIGITS;
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
static void printArray(const char *name, int *array, int arrayLen)
{
    if (arrayLen <= ARRAY_PRINT_THRESHOLD)
    {
        printf("Array %s: [", name);
        for (int i = 0; i < arrayLen; i++)
        {
            printf("%d", array[i]);

            if (i != arrayLen - 1)
            {
                printf(", ");
            }
        }
        printf("]\n");
    }
}

/**
 * @brief Checks if the array is sorted.
 *
 * @param array The array to check.
 * @param arrayLen The array length.
 *
 * @returns `true` if the array is sorted, `false` otherwise.
 */
static bool isSorted(int *array, int arrayLen)
{
    bool output = true;
    for (int i = 1; i < arrayLen; i++)
    {
        if (array[i - 1] > array[i])
        {
            output = false;
            break;
        }

        if (!output)
        {
            break;
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
static void computeOffsets(int **countMatrix, int numSections, int numValues, int **offsetMatrix)
{
// Initialize the offset matrix
#pragma omp parallel for
    for (int m = 0; m < numSections; m++)
    {
        for (int n = 0; n < numValues; n++)
        {
            offsetMatrix[m][n] = 0;

            // First term: Sum of counts for all values < n in all sections
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < numSections; j++)
                {
                    offsetMatrix[m][n] += countMatrix[j][i];
                }
            }

            // Second term: Sum of counts for value n in all sections < m
            for (int j = 0; j < m; j++)
            {
                offsetMatrix[m][n] += countMatrix[j][n];
            }
            // Debug print to check offsetMatrix values
            // printf("OffsetMatrix[%d][%d] = %d\n", m, n, offsetMatrix[m][n]);
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
static void computeLocalOffsets(const int *localArray, const int localArraySize, int **offsetMatrix,
                                const int numValues, const int rank, int *offsets, const int digit)
{
    // Create a local copy of the offsets for the current process to track updates
    int *localOffsets = new int[numValues];
#pragma omp parallel for
    for (int i = 0; i < numValues; i++)
    {
        localOffsets[i] = offsetMatrix[rank][i];
    }

    // Assign indexes to each value in the local array
    for (int i = 0; i < localArraySize; i++)
    {
        int value = (localArray[i] / digit) % NUM_DIGITS;
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
static void placeValuesFromOffset(int *localArray, int localArraySize, int *offsets, int *globalArray)
{
#pragma omp parallel for
    for (int i = 0; i < localArraySize; i++)
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
static void updateCountMatrix(int *countMatrix, const int *localArray, const int localArraySize, const int digit)
{
    // Reset the count matrix for the current process
    (void)memset(countMatrix, 0, sizeof(int) * COUNT_ARRAY_SIZE);

// Count the occurrences of each digit at the current place value in the local array
#pragma omp parallel for
    for (int i = 0; i < localArraySize; i++)
    {
        int digitValue = (localArray[i] / digit) % NUM_DIGITS;
#pragma omp critical
        {
            countMatrix[digitValue]++;
        }
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
    int inputArraySize = 0, maxDigit = 0, maxPossibleValue = 0;
    int *inputArray = NULL;
    int *outputArray = NULL;

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
        int omp_threads = atoi(argv[3]);
        omp_set_num_threads(omp_threads);
        printf("DEBUG: num threads: %d\n", omp_threads);
#endif

        outputArray = new int[inputArraySize];
        for (int i = 0; i < inputArraySize; i++)
        {
            outputArray[i] = 0; // Initialize with a default value
        }

        if (readIntArrayFromFile(inputFileName, inputArray, inputArraySize))
        {
            fprintf(stderr, "Could not read array from file %s\n", inputFileName);
            MPI_Abort(comm, 1);
        }

        int maxValue;
        if (getMax(inputArray, inputArraySize, &maxValue))
        {
            fprintf(stderr, "Could not get the maximum value in the input array\n");
            MPI_Abort(comm, 1);
        }

        // find out the number of digits in this maximum value
        maxDigit = getNumDigits(maxValue);

        printArray("Initial", inputArray, inputArraySize);

        // FIXME: evil max possible value, don't like this
        maxPossibleValue = (maxDigit > 9) ? __INT_MAX__ : (int)(pow(NUM_DIGITS, maxDigit));
    }

    MPI_Bcast(&inputArraySize, 1, MPI_INT, 0, comm);

    // the base size of our local array
    int baseLocalArraySize = inputArraySize / nproc;
    int remainder = inputArraySize % nproc; // remainder for the rest of the array

    // since remainder is less than nproc, if the rank is less than the remainder, add 1 to it
    int localArraySize = (rank < remainder) ? baseLocalArraySize + 1 : baseLocalArraySize;

    // Create a local count array to store the count of each digit for the current process
    int *localCountArray = new int[COUNT_ARRAY_SIZE];

    // Create a local offset array to store the offsets for the current process
    int *localOffsetArray = new int[localArraySize];

    // local array to sort
    int *localArray = new int[localArraySize];

    // TODO: use 2d representation of matrices instead of doing this weird stuff to work with MPI

    // Output Array
    // Allocate countMatrix and offsetMatrix as contiguous blocks
    int *flatCountMatrix = new int[nproc * COUNT_ARRAY_SIZE];
    int **countMatrix = new int *[nproc];
#pragma omp parallel for
    for (int i = 0; i < nproc; i++)
    {
        countMatrix[i] = &flatCountMatrix[i * COUNT_ARRAY_SIZE];
    }

    int *flatOffsetMatrix = new int[nproc * COUNT_ARRAY_SIZE];
    int **offsetMatrix = new int *[nproc];
#pragma omp parallel for
    for (int i = 0; i < nproc; i++)
    {
        offsetMatrix[i] = &flatOffsetMatrix[i * COUNT_ARRAY_SIZE];
    }

    // set up the temporary arrays to be able to do scatters and gathers
    // set up the temp arrays
    tempSendRecvCounts = new int[nproc];
    tempDisplacements = new int[nproc];

    // NOTE: cannot parallelize due to a dependency on the last iteration
    int curDisplacement = 0; // the current displacement
    for (int i = 0; i < nproc; i++)
    {
        // set up the temp receive counts array and calculate their values if
        // there is a remainder
        if (rank < remainder)
        {
            tempSendRecvCounts[i] = baseLocalArraySize + 1;
        }
        else
        {
            tempSendRecvCounts[i] = baseLocalArraySize;
        }

        // add the current displacement to the temp displacement array
        tempDisplacements[i] = curDisplacement;
        curDisplacement += tempSendRecvCounts[i];
    }

    // NOTE: do not parallelize anything that calls MPI functions
    for (int digit = 1; digit < maxPossibleValue; digit *= NUM_DIGITS)
    {
        // scatter the input array into local arrays
        MPI_Scatterv(inputArray, tempSendRecvCounts, tempDisplacements, MPI_INT, localArray, localArraySize, MPI_INT, 0, comm);

        // update the local count array as the matrix
        updateCountMatrix(localCountArray, localArray, localArraySize, digit);

        // Gather localCountArray into flatCountMatrix
        MPI_Gather(localCountArray, COUNT_ARRAY_SIZE, MPI_INT, flatCountMatrix, COUNT_ARRAY_SIZE, MPI_INT, 0, comm);
        MPI_Bcast(flatCountMatrix, nproc * COUNT_ARRAY_SIZE, MPI_INT, 0, comm);

        // compute offsets
        computeOffsets(countMatrix, nproc, COUNT_ARRAY_SIZE, offsetMatrix);

        // gather offsetMatrix
        MPI_Gather(offsetMatrix[rank], COUNT_ARRAY_SIZE, MPI_INT, flatOffsetMatrix, COUNT_ARRAY_SIZE, MPI_INT, 0, comm);
        MPI_Bcast(flatOffsetMatrix, nproc * COUNT_ARRAY_SIZE, MPI_INT, 0, comm);

        // compute local offsets
        computeLocalOffsets(localArray, localArraySize, offsetMatrix, COUNT_ARRAY_SIZE, rank, localOffsetArray, digit);

        int *tempOffsetArray = new int[inputArraySize];
        MPI_Gatherv(localOffsetArray, localArraySize, MPI_INT, tempOffsetArray, tempSendRecvCounts, tempDisplacements, MPI_INT, 0, comm);

        // do the move values
        if (rank == 0)
        {
            placeValuesFromOffset(inputArray, inputArraySize, tempOffsetArray, outputArray);
            // Copy the output array back to the input array for the next iteration
            for (int i = 0; i < inputArraySize; i++)
            {
                inputArray[i] = outputArray[i];
            }
        }

        delete[] tempOffsetArray;
        tempOffsetArray = NULL;
    }

    MPI_Barrier(comm); // Ensure all processes are done

    if (rank == 0)
    {

        printArray("Final", outputArray, inputArraySize);

        if (isSorted(outputArray, inputArraySize))
        {
            printf("The array is sorted.\n");
        }
        else
        {
            printf("The array is not sorted.\n");
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