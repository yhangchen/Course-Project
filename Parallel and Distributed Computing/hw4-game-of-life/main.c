#include "CudaLife.h"
#include <string.h>

#define SEED 0

int main(int argc, char const *argv[])
{
    size_t spaceWidth = 2000, spaceHeight = 2000, iteration = 100, num_CPU_threads = 8, num_GPU_threads = 16;
    if (argc > 1)
        spaceWidth = atoi(argv[1]);
    if (argc > 2)
        spaceHeight = atoi(argv[2]);
    if (argc > 3)
        iteration = atoi(argv[3]);
    if (argc > 4)
        num_CPU_threads = atoi(argv[4]);
    if (argc > 5)
        num_GPU_threads = atoi(argv[5]);

    size_t spaceSize = spaceWidth * spaceHeight, mallocSize = spaceSize * sizeof(uchar);
    uchar *lifeSpace_0, *lifeSpace_1, *lifeSpace_2, *lifeSpace_3, *lifeSpace_4;

    lifeSpace_0 = (uchar *)malloc(mallocSize);
    lifeSpace_1 = (uchar *)malloc(mallocSize);
    lifeSpace_2 = (uchar *)malloc(mallocSize);
    lifeSpace_3 = (uchar *)malloc(mallocSize);
    lifeSpace_4 = (uchar *)malloc(mallocSize);

    printf("************TEST START************\n\n");

    init_space(lifeSpace_0, spaceSize, SEED);

    memcpy(lifeSpace_1, lifeSpace_0, mallocSize);
    long long cpu_start_time = start_timer();
    GameofLife_CPU(lifeSpace_1, lifeSpace_2, spaceWidth, spaceHeight, iteration, num_CPU_threads);
    printf("  Num CPU Threads: %d\n", num_CPU_threads);
    long long cpu_time = stop_timer(cpu_start_time, "CPU");

    memcpy(lifeSpace_1, lifeSpace_0, mallocSize);
    long long gpu_start_time = start_timer();
    GameofLife_simpleGPU(lifeSpace_1, lifeSpace_3, spaceWidth, spaceHeight, iteration, num_GPU_threads);
    long long gpu_time = stop_timer(gpu_start_time, "simpleGPU");

    memcpy(lifeSpace_1, lifeSpace_0, mallocSize);
    gpu_start_time = start_timer();
    GameofLife_tilingGPU(lifeSpace_1, lifeSpace_4, spaceWidth, spaceHeight, iteration, num_GPU_threads);
    gpu_time = stop_timer(gpu_start_time, "tilingGPU");

    char FLAG = 1;
    for (size_t i = 0; i < spaceSize; i++)
    {
        if (lifeSpace_2[i] != lifeSpace_3[i] || lifeSpace_2[i] != lifeSpace_4[i])
        {
            printf("GPU results are incorrect. The first incorrect position is\n");
            printf("%d, ", i);
            FLAG = 0;
            break;
        }
    }

    if (FLAG == 1)
        printf("GPU results are correct.\n\n************TEST END************\n\n");
    free(lifeSpace_0);
    free(lifeSpace_1);
    free(lifeSpace_2);
    free(lifeSpace_3);
    free(lifeSpace_4);
}
