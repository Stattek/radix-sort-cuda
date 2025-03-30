/**
 * Summary: Radix sort implementation.
 */
#include "radix_sort.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>

int getMax(int *array, int n)
{
    int max = array[0];
    for (int i = 1; i < n; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
        }
    }
    return max;
}

void generateRandomArray(int *array, int arraySize, int maxDigit)
{
    srand(time(0));
    for (int i = 0; i < arraySize; i++)
    {
        array[i] = rand() % (int)(pow(10, maxDigit));
    }
}