#include "counting_sort.h"
#include <stdlib.h>
#include <string.h>
#include <cmath>

/**
 * @brief Finds the maxmimum value in an array and outputs it.
 *
 * @param array The array to search.
 * @param arrayLen The length of the array.
 * @param output The output maximum value.
 *
 * @returns `true` on failure, `false` on success.
 */
bool getMax(int *array, int arrayLen, int *output)
{
    if (!array || !output)
    {
        return true; // fail
    }

    int max = array[0];
    // find the maximum
    for (int i = 1; i < arrayLen; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
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
    if (getMax(array, arrayLen, &maxElement))
    {
        return true;
    }

    // create the count array
    int countArrayLen = maxElement + 1;
    int *countArray = (int *)calloc(countArrayLen, sizeof(int));

    // count the values and put them in the countArray
    for (int i = 0; i < arrayLen; i++)
    {
        int tempIndex = array[i] % (int)pow(10, digit);
        tempIndex /= (int)pow(10, digit - 1);
        countArray[(array[i] % (int)pow(10, digit))]++;
    }

    for (int i = 1; i < countArrayLen; i++)
    {
        countArray[i] = countArray[i - 1] + countArray[i];
    }

    for (int i = arrayLen - 1; i >= 0; i--)
    {
        outputArray[--countArray[array[i]]] = array[i];
    }

    free(countArray);
    countArray = NULL;

    return false; // success
}