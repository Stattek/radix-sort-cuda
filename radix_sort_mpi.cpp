/*
 * Compile: mpic++ -Wall radix_sort_mpi.cpp -o rsort_mpi
 * Run: mpirun -np <num_processes> ./rsort_mpi <sizeOfArray> <#Digits>
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
using namespace std;

inline void generateRandomArray(int *array, int arraySize, int maxVal)
{
    srand(time(0));
    for (int i = 0; i < arraySize; i++)
    {
        array[i] = rand() % maxVal;
    }
}

inline void printArray(int *array, int arraySize)
{
    for (int i = 0; i < arraySize; i++)
    {
        printf("%d ", array[i]);
    }
    printf("\n");
}

void countSort(int *array, int arraySize, int digit)
{
    int *output = new int[arraySize];
    int count[10] = {0};

    // Count the occurrences of each digit at the current place value
    for (int i = 0; i < arraySize; i++)
    {
        count[(array[i] / digit) % 10]++;
    }
    // Update count[i] to store the position of the next occurrence of the digit in the output array
    for (int i = 1; i < 10; i++)
    {
        count[i] += count[i - 1];
    }

    // Build the output array by placing elements in their correct sorted position
    for (int i = arraySize - 1; i >= 0; i--)
    {
        output[--count[(array[i] / digit) % 10]] = array[i];
    }

    // Copy the sorted elements from the output array back to the original array
    for (int i = 0; i < arraySize; i++)
    {
        array[i] = output[i];
    }

    delete[] output;
}

inline bool isSorted(int *array, int arraySize)
{
    bool check = true;
    for (int i = 1; i < arraySize; i++)
    {
        if (array[i - 1] > array[i])
        {
            printf("Error at: array[%d] = %d, array[%d] = %d\n", i - 1, array[i - 1], i, array[i]);
            check = false;
        }
    }
    return check;
}

/**
 */
void computeOffsets(int **countMatrix, int numSections, int numValues, int **offsetMatrix)
{

    // Initialize the offset matrix
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
 * Computes offsets for the local part of the array using the global offset matrix.
 */
void computeLocalOffsets(int *localArray, int localArraySize, int **offsetMatrix, int numValues, int rank, int *offsets, int digit)
{
    // Create a local copy of the offsets for the current process to track updates
    int *localOffsets = new int[numValues];
    for (int i = 0; i < numValues; i++)
    {
        localOffsets[i] = offsetMatrix[rank][i];
    }

    // Assign indexes to each value in the local array
    for (int i = 0; i < localArraySize; i++)
    {
        int value = (localArray[i] / digit) % 10;
        offsets[i] = localOffsets[value]; // Assign the current offset for the value
        localOffsets[value]++;            // Increment the offset for the next occurrence
        // Debug print to check offsets
        printf("Process %d: Value %d assigned offset %d\n", rank, value, offsets[i]);
    }

    delete[] localOffsets;
}
/**
 * Moves values from the local array to their correct positions in the global array
 * using the computed offsets.
 */
void moveValues(int *localArray, int localArraySize, int *offsets, int *globalArray, int globalArraySize)
{
    for (int i = 0; i < localArraySize; i++)
    {
        // Debug print to check offsets and values
        // printf("Process: Writing value %d at index %d\n", localArray[i], offsets[i]);

        // Check for out-of-bounds access
        if (offsets[i] < 0 || offsets[i] >= globalArraySize)
        {
            printf("Error: Offset %d is out of bounds for globalArray\n", offsets[i]);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        globalArray[offsets[i]] = localArray[i];
    }
}
/**
 * Updates the count matrix for the local array.
 */
void updateCountMatrix(int *countMatrix, int *localArray, int localArraySize, int digit, int rank)
{
    // Reset the count matrix for the current process
    for (int i = 0; i < 10; i++)
    {
        countMatrix[i] = 0;
    }

    // Count the occurrences of each digit at the current place value in the local array
    for (int i = 0; i < localArraySize; i++)
    {
        int digitValue = (localArray[i] / digit) % 10;
        countMatrix[digitValue]++;
    }
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Usage: %s <sizeOfArray> <#Digits>\n", argv[0]);
        return -1;
    }

    int inputArraySize = atoi(argv[1]);
    int maxDigit = atoi(argv[2]);
    int maxPossibleValue = (maxDigit > 9) ? __INT_MAX__ : (int)(pow(10, maxDigit));
    int *outputArray = new int[inputArraySize];
    for (int i = 0; i < inputArraySize; i++)
    {
        outputArray[i] = 0; // Initialize with a default value
    }
    MPI_Init(&argc, &argv);

    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int *inputArray = nullptr;
    int baseLocalArraySize = inputArraySize / nproc;
    int remainder = inputArraySize % nproc;
    int localArraySize = (rank < remainder) ? baseLocalArraySize + 1 : baseLocalArraySize;
    // Create a local count array to store the count of each digit for the current process
    int *localCountArray = new int[10]();
    // Create a local offset array to store the offsets for the current process
    int *localOffsetArray = new int[localArraySize]();

    int *localArray = new int[localArraySize];

    // temp array for saving results
    int *tempRecvCounts = NULL;
    int *tempDispls = NULL;

    // Output Array

    // Allocate countMatrix and offsetMatrix as contiguous blocks
    int *flatCountMatrix = new int[nproc * 10]();
    int **countMatrix = new int *[nproc];
    for (int i = 0; i < nproc; i++)
    {
        countMatrix[i] = &flatCountMatrix[i * 10];
    }

    int *flatOffsetMatrix = new int[nproc * 10]();
    int **offsetMatrix = new int *[nproc];
    for (int i = 0; i < nproc; i++)
    {
        offsetMatrix[i] = &flatOffsetMatrix[i * 10];
    }

    if (rank == 0)
    {
        inputArray = new int[inputArraySize];
        generateRandomArray(inputArray, inputArraySize, maxPossibleValue);
        printf("Initial input:\n");
        printArray(inputArray, inputArraySize);

        // set up the temp array
        tempRecvCounts = new int[nproc];
        tempDispls = new int[nproc];
        int sum = 0;
        for (int i = 0; i < nproc; i++)
        {
            if (rank < remainder)
            {
                tempRecvCounts[i] = baseLocalArraySize + 1;
            }
            else
            {
                tempRecvCounts[i] = baseLocalArraySize;
            }

            tempDispls[i] = sum;
            sum += tempRecvCounts[i];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Calculate displacements and send counts for scattering
    int *sendCounts = new int[nproc];
    int *displacements = new int[nproc];
    for (int i = 0; i < nproc; i++)
    {
        sendCounts[i] = (i < remainder) ? baseLocalArraySize + 1 : baseLocalArraySize;
        displacements[i] = (i == 0) ? 0 : displacements[i - 1] + sendCounts[i - 1];
    }

    for (int digit = 1; digit < maxPossibleValue; digit *= 10)
    {
        MPI_Scatterv(inputArray, sendCounts, displacements, MPI_INT, localArray, localArraySize, MPI_INT, 0, MPI_COMM_WORLD);

        printf("Process %d received local array:\n", rank);
        printArray(localArray, localArraySize);

        // OLD: countSort(localArray, localArraySize, digit);
        updateCountMatrix(localCountArray, localArray, localArraySize, digit, rank);
        MPI_Barrier(MPI_COMM_WORLD);
        // Gather localCountArray into flatCountMatrix
        MPI_Gather(localCountArray, 10, MPI_INT, flatCountMatrix, 10, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(flatCountMatrix, nproc * 10, MPI_INT, 0, MPI_COMM_WORLD);
        // compute offsets
        computeOffsets(countMatrix, nproc, 10, offsetMatrix);
        // gather offsetMatrix
        MPI_Gather(offsetMatrix[rank], 10, MPI_INT, flatOffsetMatrix, 10, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(flatOffsetMatrix, nproc * 10, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        // compute local offsets
        computeLocalOffsets(localArray, localArraySize, offsetMatrix, 10, rank, localOffsetArray, digit);
        MPI_Barrier(MPI_COMM_WORLD);
#if 0
        // move values
        moveValues(localArray, localArraySize, localOffsetArray, outputArray, inputArraySize);
#endif

        int *tempOffsetArray = new int[inputArraySize];
        MPI_Gatherv(localOffsetArray, localArraySize, MPI_INT, tempOffsetArray, tempRecvCounts, tempDispls, MPI_INT, 0, MPI_COMM_WORLD);
        printf("tempOffsetArray\n");
        printArray(tempOffsetArray, inputArraySize);

        // do the move values
        if (rank == 0)
        {
            moveValues(inputArray, localArraySize, localOffsetArray, outputArray, inputArraySize);
        }

        printf("proc %d - Array after digit %d:\n", rank, digit);
        printArray(outputArray, inputArraySize);
        int * temp = inputArray;
        inputArray = outputArray;
        outputArray = temp;
        
        MPI_Barrier(MPI_COMM_WORLD);
#if 0
        // Reduce the outputArray across all processes to gather the final sorted array
        MPI_Reduce(rank == 0 ? MPI_IN_PLACE : outputArray, outputArray, inputArraySize, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

        MPI_Bcast(inputArray, inputArraySize, MPI_INT, 0, MPI_COMM_WORLD);
#endif

        // do something here
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (rank == 0)
    {
        printf("Offset Matrix:\n");
        for (int i = 0; i < nproc; i++)
        {
            for (int j = 0; j < 10; j++)
            {
                printf("%d ", offsetMatrix[i][j]);
            }
            printf("\n");
        }

        printf("Count Matrix:\n");
        for (int i = 0; i < nproc; i++)
        {
            for (int j = 0; j < 10; j++)
            {
                printf("%d ", countMatrix[i][j]);
            }
            printf("\n");
        }
    }

    // MPI_Gatherv(localArray, localArraySize, MPI_INT, inputArray, sendCounts, displacements, MPI_INT, 0, MPI_COMM_WORLD);
    //  Print the count matrix for debugging

    MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes synchronize
    if (rank == 0)
    {

        printf("Final Output:\n");
        printArray(outputArray, inputArraySize);

        if (isSorted(outputArray, inputArraySize))
        {
            printf("The array is sorted.\n");
        }
        else
        {
            printf("The array is not sorted.\n");
        }

        delete[] inputArray;
    }
    delete[] outputArray;
    delete[] flatCountMatrix;
    delete[] countMatrix;
    delete[] flatOffsetMatrix;
    delete[] offsetMatrix;
    delete[] localArray;
    delete[] sendCounts;
    delete[] displacements;

    MPI_Finalize();
    return 0;
}