/*
 * Compile: g++ -Wall radix_sort.cpp -o rsort
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
using namespace std;

inline void generateRandomArray(int *array, int arraySize, int maxVal)
{
    srand(time(0));
    for (int i = 0; i < arraySize; i++)
    {
        array[i] = rand() % maxVal;
    }
}

inline void printArray(int *array, int arraySize)
{
    for (int i = 0; i < arraySize; i++)
    {
        printf("%d ", array[i]);
    }
    printf("\n");
}

inline void countSort(int *inputArray, int inputArraySize, int *outArray, int digit)
{
    int countArray[10];
    for (int i = 0; i < 10; i++)
    {
        countArray[i] = 0;
    }

    for (int i = 0; i < inputArraySize; i++)
    {

        countArray[(inputArray[i] / digit) % 10]++;
    }
    int sum = 0;
    for (int i = 0; i < 10; i++)
    {
        sum += countArray[i];
        countArray[i] = sum;
    }

    for (int i = inputArraySize - 1; i >= 0; i--)
    {

        countArray[(inputArray[i] / digit) % 10]--;
        int newIndex = countArray[(inputArray[i] / digit) % 10];
        outArray[newIndex] = inputArray[i];
    }
    int *temp = inputArray;
    inputArray = outArray;
    outArray = temp;
}

int main(int argc, char *argv[])
{

    if (argc != 3)
    {
        printf("Usage: %s <sizeOfArray> <#Digits>\n", argv[0]);
        return -1;
    }

    int inputArraySize = atoi(argv[1]);
    int maxDigit = atoi(argv[2]);
    int maxPossibleValue;

    if (maxDigit > 9)
    {
        printf("Max Value = 2147483647. Integer cap reached, all digit values 10 and above are set to the integer cap.\n");
        maxPossibleValue = 2147483647;
    }
    else
    {
        maxPossibleValue = (int)(pow(10, maxDigit));
        printf("Max Value: %d\n", maxPossibleValue);
    }

    int *inputArray = new int[inputArraySize];
    int *outArray = new int[inputArraySize];

    generateRandomArray(inputArray, inputArraySize, maxPossibleValue);

    printf("Initial input:\n");
    printArray(inputArray, inputArraySize);

    for (long int i = 1; i < maxPossibleValue; i *= 10)
    {
        countSort(inputArray, inputArraySize, outArray, i);
    }

    printf("Final Output:\n");
    printArray(outArray, inputArraySize);

    delete[] inputArray;
    delete[] outArray;
}
