#include <stdio.h>
#include <iostream>
#include "radix_sort.h"

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        printf("Usage: %s <sizeOffArray> <maxDigit#>\n", argv[0]);
        return -1;
    }

    int arraySize = atoi(argv[1]);
    int maxDigit = atoi(argv[2]) - 1;
    int array[maxDigit];

    generateRandomArray(array, arraySize, maxDigit);

    int max = getMax(array, arraySize);
}
