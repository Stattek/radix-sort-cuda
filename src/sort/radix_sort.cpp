/**
 * Summary: Radix sort implementation.
 */
#include "radix_sort.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>

bool getMax(int *array, int n, int *output)
{
    if (!array || !output)
    {
        return true;
    }

    int max = array[0];
    // find the maximum
    for (int i = 1; i < n; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
        }
    }

    *output = max;
    return false;
}

bool generateRandomArray(int *array, int arraySize, int maxDigit)
{
    if (!array)
    {
        return true; // error
    }

    srand(time(0));
    for (int i = 0; i < arraySize; i++)
    {
        array[i] = rand() % (int)(pow(10, maxDigit));
    }

    return false; // success
}