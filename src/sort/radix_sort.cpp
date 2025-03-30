/**
 * Summary: Radix sort implementation.
 */
#include "radix_sort.h"
#include "counting_sort.h"
#include <stdlib.h>
#include <stdio.h>

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