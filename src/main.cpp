#include <stdio.h>
#include <iostream>
#include "radix_sort.h"
#include "counting_sort.h"
#include <stdlib.h>
#include <cmath>
#include <iostream>

#define PRINT_THRESHOLD 50

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
    printf("Array %s (Length %d):\n", name, arrayLen);

    if (arrayLen <= PRINT_THRESHOLD)
    {
        printf("[");
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
    return 0;
}
