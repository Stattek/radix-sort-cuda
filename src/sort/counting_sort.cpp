#include "counting_sort.h"
#include <stdlib.h>
#include <string.h>

bool getMax(int *array, int arrayLen, int *output)
{
    if (!array || !output)
    {
        return true;
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

bool countingSort(int *array, int arrayLen, int *outputArray)
{
    if (!array)
    {
        return true; // failure
    }
    int maxElement = 0;
    if (getMax(array, arrayLen, &maxElement))
    {
        return true;
    }

    int *countArray = (int *)calloc(maxElement + 1, sizeof(int));
    int digit = 1; // the digit we are focused on

    for (int i = 0; i < arrayLen; i++)
    {
        countArray[array[i] % (10 * digit)]++;
    }

    for (int i = 0; i < arrayLen; i++)
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