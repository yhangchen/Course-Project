#ifndef CUDA_LIFE_H
#define CUDA_LIFE_H
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#define min(a, b) (a < b ? a : b)
#define index(i, L) ((0) ? (i + L) : ((i >= L) ? (i - L) : (i)))
typedef unsigned char uchar;
typedef unsigned int uint;

void GameofLife_simpleGPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t numThread);

void GameofLife_tilingGPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t numThread);

// initialize a game-of-life-space.
static inline void init_space(uchar *lifeSpace, size_t spaceSize, size_t SEED)
{
    srand(SEED);
    size_t population = 0;

    for (size_t i = 0; i < spaceSize; i++)
    {
        lifeSpace[i] = rand() % 2;
        if (lifeSpace[i] == 1)
            population++;
    }
    printf("Initial Population is %zu\n", population);
}

// Returns the current time in microseconds
static inline long long start_timer()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000000 + tv.tv_usec;
}

// Prints the time elapsed since the specified time
static inline long long stop_timer(long long start_time, char *name)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    long long end_time = tv.tv_sec * 1000000 + tv.tv_usec;
    printf("%s: %.5f sec\n\n", name, ((float)(end_time - start_time)) / (1000 * 1000));
    return end_time - start_time;
}

static inline void GameofLife_CPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t num_CPU_threads)
{
    size_t xId, yId, xLeft, xRight, yBottom, yAbove;
    size_t numAlive;
    size_t population = 0;
    size_t spaceSize = spaceWidth * spaceHeight;
    // naive implementation

    for (size_t t = 0; t < iteration; t++)
    {
#pragma omp parallel num_threads(num_CPU_threads) shared(lifeSpace_1, lifeSpace_2, spaceWidth, spaceSize) default(none) private(xId, yId, xLeft, xRight, yBottom, yAbove, numAlive)
#pragma omp for
        for (uint cellId = 0; cellId < spaceSize; cellId++)
        {
            xId = cellId % spaceWidth;
            yId = cellId - xId;
            xLeft = (xId + spaceWidth - 1) % spaceWidth;
            xRight = (xId + 1) % spaceWidth;
            yBottom = (yId + spaceSize - spaceWidth) % spaceSize;
            yAbove = (yId + spaceWidth) % spaceSize;

            numAlive = lifeSpace_1[xLeft + yId] + lifeSpace_1[xLeft + yAbove] + lifeSpace_1[xLeft + yBottom] + lifeSpace_1[xId + yAbove] + lifeSpace_1[xId + yBottom] + lifeSpace_1[xRight + yId] + lifeSpace_1[xRight + yAbove] + lifeSpace_1[xRight + yBottom];

            lifeSpace_2[cellId] = (numAlive == 3 || (lifeSpace_1[cellId] && numAlive == 2)) ? 1 : 0;
        }
#pragma omp for
        for (size_t cellId = 0; cellId < spaceSize; cellId++)
            lifeSpace_1[cellId] = lifeSpace_2[cellId];
    }

    for (size_t cellId = 0; cellId < spaceSize; cellId++)
    {
        lifeSpace_2[cellId] = lifeSpace_1[cellId];
        if (lifeSpace_2[cellId] == 1)
            population++;
    }
    // report results
    printf("lifeSpace_1 size: %d x %d, total generations: %d\n", spaceWidth, spaceHeight, iteration);
    printf("Final Population is %d\n\n", population);
}
#endif
