#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
using namespace std;

int getMax(int array[], int n)
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
        array[i] = rand() % (int)(pow(10,maxDigit)); 
    }
}

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


