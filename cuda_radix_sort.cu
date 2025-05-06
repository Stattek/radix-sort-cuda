#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cuda.h>
#include <vector>
#include <fstream>
#include <omp.h>
#include <math.h>

#define DO_CUDA_DEBUG 0 // enables debug prints, disables flipping bits, and uses base 10 for running

#if DO_CUDA_DEBUG
#define ARRAY_PRINT_THRESHOLD 30
#else
#define ARRAY_PRINT_THRESHOLD 20
#endif

#define NUM_THREADS 4 // number of CUDA threads

#if DO_CUDA_DEBUG
#define NUM_BASE 10
#else
#define NUM_BASE 256
#endif
#define COUNT_ARRAY_SIZE NUM_BASE // the count array will always hold the same number of values as the number of digits

/**
 * @brief Prints the array.
 *
 * @param name The name of the array.
 * @param array The array to print.
 * @param arrayLen The length of the array.
 */
static void printArray(const char *name, uint *array, uint arrayLen)
{
    printf("\nArray %s, length %d:", name, arrayLen);
    if (arrayLen <= ARRAY_PRINT_THRESHOLD)
    {
        printf(" [");
        for (uint i = 0; i < arrayLen; i++)
        {
            printf("%d", array[i]);

            if (i != arrayLen - 1)
            {
                printf(", ");
            }
        }
        printf("]\n");
    }
    else
    {
        printf(" Array above printing threshold %d\n", ARRAY_PRINT_THRESHOLD);
    }
}

/**
 * @brief Reads an integer array from a file.
 *
 * @param fileName The name of the file to read from.
 * @param outputNumElements The length of the output array.
 * @note The output array is allocated on the heap and MUST be deallocated by the user.
 *
 * @returns True on failure, false on success.
 */
static uint *readIntArrayFromFile(const char *fileName, uint &outputNumElements)
{
    if (!fileName)
    {
        // bad argument
        return NULL;
    }

    std::vector<uint> output;

    // since this can fail when reading from input
    try
    {
        std::ifstream inputFile(fileName, std::ios_base::in);
        int curInt;
        while (inputFile >> curInt)
        {
            output.push_back((uint)curInt);
        }
    }
    catch (...)
    {
        return NULL; // fail
    }

    if (output.empty())
    {
        // either no integers or a bad file name
        return NULL;
    }

    uint *outputPointer = new uint[output.size()];
    for (uint i = 0; i < (uint)output.size(); i++)
    {
        outputPointer[i] = output[i];
    }
    outputNumElements = output.size();

    return outputPointer; // success
}

/**
 * @brief Flips all the sign bits in an array.
 *
 * @param array The array to flip sign bits in.
 * @param arrayLength The lengh of the arrray.
 */
static void flipSignBits(int *array, uint arrayLength)
{
    // create mask for flipping the sign bit
    uint mask = -1;
    mask -= 1;
    mask = (mask >> 1); // shift one to the right so we flip only the sign bit
    mask = ~mask;

    // flip all the sign bits
    for (uint i = 0; i < arrayLength; i++)
    {
        array[i] ^= mask;
    }
}

/**
 * @brief Finds the maxmimum value in an array at a digit and outputs it.
 *
 * @param array The array to search.
 * @param arrayLen The length of the array.
 * @param digit The digit to find the maximum value for.
 * @param output The output maximum value.
 *
 * @returns `true` on failure, `false` on success.
 */
static bool getMax(const uint *array, const uint arrayLen, uint *output)
{
    if (!array || !output)
    {
        return true; // fail
    }

    uint maxValue = array[0];
    // find the maximum
    for (uint i = 1; i < arrayLen; i++)
    {
        if (maxValue < array[i])
        {
            maxValue = array[i];
        }
    }

    *output = maxValue;
    return false; // success
}

/**
 * @brief Gets the number of digits that this value has.
 *
 * @param value The value to find the number of digits of.
 * @return The number of digits.
 */
static uint getNumDigits(uint value)
{
    uint numDigits = 1; // we should start at one
    value /= NUM_BASE;
    while (value > 0)
    {
        numDigits++;
        value /= NUM_BASE;
    }

    return numDigits;
}

/**
 * @brief Finds the power of a value.
 *
 * @param value The value to find power of.
 * @param exponent The exponent.
 * @return The result of value to the power of exponent.
 */
static unsigned long long myPowUint(uint value, uint exponent)
{
    unsigned long long sum = 1;
#pragma omp parallel for reduction(* : sum)
    for (uint i = 0; i < exponent; i++)
    {
        sum *= value;
    }

    return sum;
}

/**
 * @brief Calculates elapsed time from two times.
 *
 * @param initialTime The initial time.
 * @param finalTime The final time.
 * @returns The elapsed time.
 */
static double calculateElapsedTime(timespec &initialTime, timespec &finalTime)
{
    double output = (finalTime.tv_sec - initialTime.tv_sec) + (finalTime.tv_nsec - initialTime.tv_nsec);
    if (output < 0)
    {
        return 0; // bad time
    }

    return output;
}

/**
 * @brief Checks if the array is sorted.
 *
 * @param array The array to check.
 * @param arrayLen The array length.
 *
 * @returns `true` if the array is sorted, `false` otherwise.
 */
static bool isSorted(int *array, uint arrayLen)
{
    bool output = true;
#pragma omp parallel for
    for (uint i = 1; i < arrayLen; i++)
    {
        if (array[i - 1] > array[i])
        {
            output = false;
        }
    }

    return output;
}

/**
 * @brief Updates the count matrix for the local array.
 *
 * @param countMatrix The countmatrix to modify.
 * @param localArray The local array to read from.
 * @param localArraySize The local array size.
 * @param digit The current digit to update the count matrix on.
 */
__global__ void updateCountMatrix(uint *countMatrix, const uint *localArray,
                                  const uint localArraySize, const uint digit)
{
    // Count the occurrences of each digit at the current place value in the local array
    int index = (blockIdx.x * blockDim.x) + threadIdx.x;
#if DO_CUDA_DEBUG // DEBUG: debug print
    printf("DEBUG: index=%d\n", index);
#endif
    if (index < localArraySize)
    {
        uint digitValue = (localArray[index] / digit) % NUM_BASE;
        countMatrix[(blockIdx.x * COUNT_ARRAY_SIZE) + digitValue]++;

#if DO_CUDA_DEBUG // DEBUG: debug print
        printf("DEBUG: index=%d, blockIdx.x=%d,blockDim.x=%d, threadIdx.x=%d, countmatrix[%d]=%u, digitValue=%u\n",
               index, blockIdx.x,
               blockDim.x, threadIdx.x, (blockIdx.x * COUNT_ARRAY_SIZE) + digitValue,
               countMatrix[(blockIdx.x * COUNT_ARRAY_SIZE) + digitValue], digitValue);
#endif
    }
}

/**
 * @brief Finds the offset table from the block sums.
 *
 * @param deviceGlobalOffset The device global offset array.
 * @param countMatrix The count matrix.
 */
__global__ void sumOffsetTable(uint *deviceGlobalOffset, uint *countMatrix)
{
    int sum = 0;

    // sum this count array
    for (int i = 0; i < COUNT_ARRAY_SIZE; i++)
    {
        sum += countMatrix[(blockIdx.x * COUNT_ARRAY_SIZE) + i];
#if DO_CUDA_DEBUG // DEBUG: debug print
        printf("DEBUG: sum=%d, i=%d, blockSumsIdx=%d, blockIdx=%d deviceGlobalOffset[blockIdx.x]=%d\n",
               sum, i, blockIdx.x, blockIdx, deviceGlobalOffset[blockIdx.x]);
#endif
    }
    deviceGlobalOffset[blockIdx.x] = sum;
#if DO_CUDA_DEBUG // DEBUG: debug print
    printf("DEBUG: blockIdx.x=%d, deviceGlobalOffset[blockIdx.x]=%d\n", blockIdx.x, deviceGlobalOffset[blockIdx.x]);
#endif
}

__global__ void shiftOffsetTable(uint *deviceGlobalOffset, uint *newBlockSums, int numBlocks, int iterationNum)
{
#if DO_CUDA_DEBUG // DEBUG: debug print
    printf("DEBUG: shiftOffsetTable iteration=%d\n", iterationNum);
#endif
    if (iterationNum == 0)
    {
        newBlockSums[0] = 0;
        int index = blockIdx.x + 1;
        if (index < numBlocks && blockIdx.x < numBlocks)
        {
            newBlockSums[index] = deviceGlobalOffset[blockIdx.x];
        }
    }
    else
    {
        int powResult = (int)pow(2, iterationNum - 1);
        if (blockIdx.x >= powResult)
        {
            int first = blockIdx.x - powResult;
            int second = blockIdx.x;
#if DO_CUDA_DEBUG // DEBUG: debug print
            printf("DEBUG: blockIdx.x=%d first=%d, second=%d, iteration=%d, pow(2, iterationNum - 1)=%d\n",
                   blockIdx.x, first, second, iterationNum, powResult);
#endif

            newBlockSums[blockIdx.x] = deviceGlobalOffset[first] + deviceGlobalOffset[second];
        }
    }
}

/**
 * @brief Performs the parallel scan algorithm to find prefix sum.
 *
 * @param numBlocks The number of blocks.
 * @param deviceGlobalOffset The global offset.
 */
static void parallelScan(int numBlocks, uint *deviceGlobalOffset)
{
    // run until log2(numBlocks)
    for (int i = 0; i <= (int)log2(numBlocks); i++)
    {
#if DO_CUDA_DEBUG
        printf("DEBUG: CPU LOOP i=%d\n", i);
#endif
        uint *deviceNewBlockSums = NULL;
        cudaMallocManaged(&deviceNewBlockSums, sizeof(uint) * numBlocks);
        cudaMemset(deviceNewBlockSums, 0, sizeof(uint) * numBlocks);
        shiftOffsetTable<<<numBlocks, 1>>>(deviceGlobalOffset, deviceNewBlockSums, numBlocks, i);
        cudaDeviceSynchronize();
        // copy over the new sums to the block sums
        cudaMemcpy(deviceGlobalOffset, deviceNewBlockSums, sizeof(uint) * numBlocks,
                   cudaMemcpyKind::cudaMemcpyDeviceToDevice);
        cudaDeviceSynchronize();
        cudaFree(deviceNewBlockSums);
        uint *tempCopy = (uint *)malloc(sizeof(uint) * numBlocks);
#if DO_CUDA_DEBUG // DEBUG: debug print
        cudaMemcpy(tempCopy, deviceGlobalOffset, sizeof(uint) * numBlocks,
                   cudaMemcpyKind::cudaMemcpyDeviceToHost);
        printArray("deviceGlobalOffset", tempCopy, numBlocks);
#endif
    }
}

/**
 * @brief Reorders elements into sorted positions.
 *
 * @param countMatrix The count matrix.
 * @param inputArray The input array to read from.
 * @param outputArray The output array to write to.
 * @param deviceGlobalOffset The device global offset array.
 * @param localArraySize The local array size.
 * @param digit The current digit.
 */
__global__ void reorderElements(uint *countMatrix, const uint *inputArray, uint *outputArray,
                                const uint *deviceGlobalOffset, const uint localArraySize,
                                const uint digit)
{
    // Count the occurrences of each digit at the current place value in the local array
    int index = (blockIdx.x * blockDim.x) + threadIdx.x;
#if DO_CUDA_DEBUG // DEBUG: debug print
    printf("DEBUG: index=%d\n", index);
#endif
    if (index < localArraySize)
    {
        uint digitValue = (inputArray[index] / digit) % NUM_BASE;
        int localOffset = countMatrix[(blockIdx.x * COUNT_ARRAY_SIZE) + digitValue];
        int globalOffset = deviceGlobalOffset[blockIdx.x];

        int globalIdx = index - localOffset + globalOffset;

        outputArray[globalIdx] = inputArray[index];

#if DO_CUDA_DEBUG // DEBUG: debug print
        printf("DEBUG: index=%d, blockIdx.x=%d,blockDim.x=%d, threadIdx.x=%d, countmatrix[%d]=%u, digitValue=%u\n",
               index, blockIdx.x,
               blockDim.x, threadIdx.x, (blockIdx.x * COUNT_ARRAY_SIZE) + digitValue,
               countMatrix[(blockIdx.x * COUNT_ARRAY_SIZE) + digitValue], digitValue);
#endif
    }
}

int main(int argc, char *argv[])
{
    /* shared variables */
    uint inputArraySize = 0, maxDigit = 0;
    unsigned long long maxPossibleValue = 0;
    uint *inputArray = NULL;
    uint *outputArray = NULL;

    // temp arrays for saving results
    int *numValues = NULL;
    int *tempDisplacements = NULL;

    /* radix sort setup */

    if (argc != 2)
    {
        printf("Usage: %s <input_file_name>\n", argv[0]);
        return 1;
    }

    const char *inputFileName = argv[1];

    inputArray = readIntArrayFromFile(inputFileName, inputArraySize);
    if (!inputArray)
    {
        fprintf(stderr, "Could not read array from file %s. File may be empty or not exist\n", inputFileName);
        return 1;
    }

    // initial info
    printf("\nInput file = \"%s\"\n", inputFileName);

    // print the array
    printArray("Initial", inputArray, inputArraySize);

    // flip bits, then do the rest of the setup
#ifndef DO_CUDA_DEBUG // TODO: bring back
    flipSignBits((int *)inputArray, inputArraySize);
#endif

    outputArray = new uint[inputArraySize];
    (void)memset(outputArray, 0, sizeof(uint) * inputArraySize);

    uint maxValue;
    if (getMax(inputArray, inputArraySize, &maxValue))
    {
        fprintf(stderr, "Could not get the maximum value in the input array\n");
        return 1;
    }

    // find out the number of digits in this maximum value
    maxDigit = getNumDigits(maxValue);

    maxPossibleValue = myPowUint(NUM_BASE, maxDigit);

    // the number of blocks we need to handle each array
    uint numBlocks = ceil((double)inputArraySize / NUM_THREADS);
#if DO_CUDA_DEBUG
    printf("DEBUG: numBlocks=%d\n", numBlocks);
#endif

    // the base size of our local array
    uint localArraySize = numBlocks * NUM_THREADS;

    // Create a local count array to store the count of each digit for the current process
    uint *localCountArray = new uint[COUNT_ARRAY_SIZE];
    // Create a local offset array to store the offsets for the current process
    uint *localOffsetArray = new uint[localArraySize];
    // local array to sort
    uint *localArray = new uint[localArraySize];

    cudaError_t err = cudaError::cudaSuccess;
    // copy the entire array to the GPU
    uint *deviceInputArray = NULL;
    err = cudaMallocManaged(&deviceInputArray, sizeof(uint) * inputArraySize);
    if (err)
    {
        fprintf(stderr, "Could not malloc device input array");
        return 1;
    }
    err = cudaMemcpy(deviceInputArray, inputArray, sizeof(uint) * inputArraySize, cudaMemcpyKind::cudaMemcpyHostToDevice);
    if (err)
    {
        fprintf(stderr, "Could not copy input array to GPU");
        return 1;
    }

    uint *deviceOutputArray = NULL;
    err = cudaMallocManaged(&deviceOutputArray, sizeof(uint) * inputArraySize);
    if (err)
    {
        fprintf(stderr, "Could not malloc device output array");
        return 1;
    }
    err = cudaMemset(deviceOutputArray, 0, sizeof(uint) * inputArraySize);
    if (err)
    {
        fprintf(stderr, "Could not memset output array");
        return 1;
    }

    // Allocate deviceCountMatrix and offsetMatrix as contiguous blocks
    uint *deviceCountMatrix = NULL;
    err = cudaMallocManaged(&deviceCountMatrix, sizeof(uint) * numBlocks * COUNT_ARRAY_SIZE);
    if (err)
    {
        fprintf(stderr, "Could not malloc count matrix on GPU");
        return 1;
    }
    err = cudaMemset(deviceCountMatrix, 0, sizeof(uint) * numBlocks * COUNT_ARRAY_SIZE);
    if (err)
    {
        fprintf(stderr, "Could not memset count matrix");
        return 1;
    }

    // initialize offset matrix
    uint *offsetMatrix = NULL;
    err = cudaMallocManaged(&offsetMatrix, sizeof(uint) * numBlocks * COUNT_ARRAY_SIZE);
    if (err)
    {
        fprintf(stderr, "Could not malloc offset matrix on GPU");
        return 1;
    }
    err = cudaMemset(offsetMatrix, 0, sizeof(uint) * numBlocks * COUNT_ARRAY_SIZE);
    if (err)
    {
        fprintf(stderr, "Could not memset offset matrix");
        return 1;
    }

    timespec startTime, finalTime;
    clock_gettime(CLOCK_REALTIME, &startTime);

    for (unsigned long long digit = 1; digit < maxPossibleValue; digit *= NUM_BASE)
    {
        /* GPU */

        // for every block, make sure that it gets values copied from CPU to the GPU
        // NOTE: we do not need to copy the input array value here because we have the whole input array on device

        // perform counts on the array for this digit
        updateCountMatrix<<<numBlocks, NUM_THREADS>>>(deviceCountMatrix,
                                                      deviceInputArray, inputArraySize, digit);
        uint *debugCountMatrix = (uint *)malloc(sizeof(uint) * numBlocks * COUNT_ARRAY_SIZE);
        cudaMemcpy(debugCountMatrix, deviceCountMatrix, sizeof(uint) * numBlocks * COUNT_ARRAY_SIZE,
                   cudaMemcpyKind::cudaMemcpyDeviceToHost);
        printArray("debugCountMatrix", debugCountMatrix, numBlocks * COUNT_ARRAY_SIZE);

        uint *deviceGlobalOffset = NULL;
        cudaMallocManaged(&deviceGlobalOffset, sizeof(uint) * numBlocks);
        // perform scan algorithm

        // find initial prefix sums
        sumOffsetTable<<<numBlocks, 1>>>(deviceGlobalOffset, deviceCountMatrix);
        cudaDeviceSynchronize();

        // perform parallel scan algorithm
        parallelScan(numBlocks, deviceGlobalOffset);

        // reorder elements back into the original array
        reorderElements<<<numBlocks, NUM_THREADS>>>(deviceCountMatrix, deviceInputArray, deviceOutputArray, deviceGlobalOffset,
                                                    localArraySize, digit);

        // swap the pointers
        uint *temp = deviceInputArray;
        deviceInputArray = deviceOutputArray;
        deviceOutputArray = temp;

        cudaFree(deviceGlobalOffset);
    }

    cudaDeviceSynchronize(); // DEBUG: remove eventually?
    clock_gettime(CLOCK_REALTIME, &finalTime);
    double elapsedTime = calculateElapsedTime(startTime, finalTime);

    // swap the pointers again :)
    if (maxPossibleValue > 1)
    {
        uint *temp = deviceInputArray;
        deviceInputArray = deviceOutputArray;
        deviceOutputArray = temp;
    }

    // copy back from device
    err = cudaMemcpy(outputArray, deviceOutputArray, sizeof(uint) * inputArraySize, cudaMemcpyKind::cudaMemcpyDeviceToHost);
    if (err)
    {
        fprintf(stderr, "Could not copy device output array to host");
        return 1;
    }

// save time
#ifndef DO_CUDA_DEBUG // TODO: bring back
    flipSignBits((int *)outputArray, inputArraySize);
#endif

    printArray("Final", outputArray, inputArraySize);

    if (isSorted((int *)outputArray, inputArraySize))
    {
        printf("\n\nThe array is sorted in %lf second(s).\n", elapsedTime);
    }
    else
    {
        printf("\n\nThe array is not sorted.\n");
        return 1;
    }

    // delete shared values
    cudaFree(deviceInputArray);
    delete[] inputArray;
    inputArray = NULL;
    delete[] outputArray;
    outputArray = NULL;
    delete[] numValues;
    numValues = NULL;
    delete[] tempDisplacements;
    tempDisplacements = NULL;
    delete[] localCountArray;
    localCountArray = NULL;
    delete[] localOffsetArray;
    localOffsetArray = NULL;
    delete[] deviceCountMatrix;
    deviceCountMatrix = NULL;
    delete[] deviceCountMatrix;
    deviceCountMatrix = NULL;
    delete[] offsetMatrix;
    offsetMatrix = NULL;
    delete[] offsetMatrix;
    offsetMatrix = NULL;
    delete[] localArray;
    localArray = NULL;

    return 0;
}
