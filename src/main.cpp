#include <stdio.h>
#include <iostream>
#include "radix_sort.h"
#include "counting_sort.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>

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
        array[i] = rand() % (maxDigit * 10);
    }

    return array;
}

int main(int argc, char *argv[])
{

    if (argc != 3)
    {
        printf("Usage: %s <arrayLen> <maxDigit#>\n\n", argv[0]);
        return 1;
    }

    int arrayLen = atoi(argv[1]);
    int maxDigit = atoi(argv[2]);

    int *array = generateRandomArray(arrayLen, maxDigit);
    if (!array)
    {
        fprintf(stderr, "ERROR: Could not generate starting array\n");
        return 1;
    }
    printf("helloworld\n");

    int *outputArray = (int *)malloc(sizeof(int) * arrayLen);
    if (!outputArray)
    {
        fprintf(stderr, "ERROR: Could not allocate output array\n");
        return 1;
    }

    if (countingSort(array, arrayLen, outputArray, 1))
    {
        fprintf(stderr, "ERROR: Could not perform counting sort\n");
        return 1;
    }
}
