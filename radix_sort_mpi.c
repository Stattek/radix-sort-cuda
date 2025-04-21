/**
 * @file radix_sort_mpi.cpp
 * @brief Radix sort implementation in OpenMPI
 * @version 0.1
 * @date 2025-04-21
 *
 * @note Compile with mpicc -g -Wall -o radix_sort_mpi.o radix_sort_mpi.c -lm
 * @note Run with mpiexec -n 1 radix_sort_mpi.o 5 3
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <stdbool.h>
#include <time.h>

bool getMax(const int *array, const int arrayLen, const int digit, int *output);
bool countingSort(int *array, const int arrayLen, int *outputArray, const int digit);

bool radixSort(int **array, const int arrayLen, const int maxDigit);

#define DIGIT_AT(value, digit) ((value % (int)pow(10, digit)) / (int)pow(10, digit - 1))

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
bool getMax(const int *array, const int arrayLen, const int digit, int *output)
{
    if (!array || !output)
    {
        return true; // fail
    }

    int max = DIGIT_AT(array[0], digit);
    // find the maximum
    for (int i = 1; i < arrayLen; i++)
    {
        if (DIGIT_AT(array[i], digit) > max)
        {
            max = DIGIT_AT(array[i], digit);
        }
    }

    *output = max;
    return false;
}

/**
 * @brief Performs counting sort on an array on the specified digit.
 *
 * @param array The array to do the counting sort on.
 * @param arrayLen The length of the input array.
 * @param outputArray The output array. Must be at least as long as `array`.
 * @param digit The significant digit to sort by.
 *
 * @returns `true` on failure, `false` on success.
 */
bool countingSort(int *array, const int arrayLen, int *outputArray, const int digit)
{
    if (!array || !outputArray)
    {
        return true; // failure
    }
    int maxElement = 0;
    if (getMax(array, arrayLen, digit, &maxElement))
    {
        return true;
    }

    // create the count array
    int countArrayLen = maxElement + 1;
    int *countArray = (int *)calloc(countArrayLen, sizeof(int));

    // count the values and put them in the countArray
    for (int i = 0; i < arrayLen; i++)
    {
        countArray[DIGIT_AT(array[i], digit)]++;
    }

    // prefix sum
    for (int i = 1; i < countArrayLen; i++)
    {
        countArray[i] = countArray[i - 1] + countArray[i];
    }

    // put values in place
    for (int i = arrayLen - 1; i >= 0; i--)
    {
        outputArray[--countArray[DIGIT_AT(array[i], digit)]] = array[i];
    }

    free(countArray);
    countArray = NULL;

    return false; // success
}

/**
 * @brief Performs radix sort on an array.
 *
 * @param inputArray The array to do the counting sort on and modify.
 * @param arrayLen The length of the input array.
 * @param maxDigit The largest significant digit to sort by.
 *
 * @returns `true` on failure, `false` on success.
 */
bool radixSort(int **array, const int arrayLen, const int maxDigit)
{
    if (!array)
    {
        return true; // error
    }

    int *outputArray = (int *)malloc(sizeof(int) * arrayLen);
    if (!outputArray)
    {
        return true;
    }

    bool isInOrder = true;
    for (int i = 1; i <= maxDigit; i++)
    {
        if (countingSort(*array, arrayLen, outputArray, i))
        {
            fprintf(stderr, "ERROR: Could not perform counting sort on digit %d\n", i);
            return true; // error
        }

        // swap the output array and input array
        int *temp = *array;
        *array = outputArray;
        outputArray = temp;
        isInOrder = !isInOrder;
    }

    // since we were swapping the pointers
    if (!isInOrder)
    {
        // swap the output array and input array
        int *temp = *array;
        *array = outputArray;
        outputArray = temp;
    }

    // free and set to the output array
    free(*array);
    *array = outputArray;

    return false;
}

/**
 * @brief Generates a random array for sorting.
 * @note The array MUST BE FREED by the user of this function.
 *
 * @param arrayLen The length of the array.
 * @param maxDigit The maximum digit to generate.
 *
 * @returns The randomly generated array or `NULL` upon error.
 */
static int *generateRandomArray(int arrayLen, int maxDigit)
{
    int *array = (int *)malloc(sizeof(int) * arrayLen);

    srand(time(0));
    for (int i = 0; i < arrayLen; i++)
    {
        array[i] = rand() % (int)pow(10, maxDigit);
    }

    return array;
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

int main(int argc, char *argv[])
{
    int nproc = 0, rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int arrayLen = 0;
    int maxDigit = 0;

    int *array = NULL;

    if (rank == 0)
    {
        if (argc != 3)
        {
            fprintf(stderr, "Usage: %s <arrayLen> <maxDigit#>\n\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        arrayLen = atoi(argv[1]);
        maxDigit = atoi(argv[2]);

        array = generateRandomArray(arrayLen, maxDigit);
        if (!array)
        {
            fprintf(stderr, "ERROR: Could not generate starting array\n");
            return 1;
        }
    }

    MPI_Bcast(&maxDigit, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int receiveNum = arrayLen / nproc;

    // scatter both arrays
    int *subarray = malloc(sizeof(int) * receiveNum);
    MPI_Scatter(array, receiveNum, MPI_INT, subarray, receiveNum, MPI_INT, 0, MPI_COMM_WORLD);

    free(subarray);

    if (rank == 0)
    {
        printArray("before", array, arrayLen);
        radixSort(&array, arrayLen, maxDigit);
        printArray("after", array, arrayLen);

        if (isSorted(array, arrayLen))
        {
            printf("The array is sorted successfully.\n");
        }
        else
        {
            fprintf(stderr, "The array is not sorted.\n");
            free(array);
            array = NULL;
            return 1;
        }

        free(array);
        array = NULL;
    }

    MPI_Finalize();
    return 0;
}