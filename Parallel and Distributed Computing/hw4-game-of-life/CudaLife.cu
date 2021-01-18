#include <cuda.h>
#include <stdio.h>
#include <assert.h>

typedef unsigned char uchar;
typedef unsigned int uint;


extern "C" void GameofLife_simpleGPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t numThread);
extern "C" void GameofLife_tilingGPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t tileWidth);


__global__ 
void simpleLifeKernel(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight)
{
    size_t spaceSize = spaceWidth * spaceHeight;
    size_t xId, yId, xLeft, xRight, yBottom, yAbove;
    size_t numAlive;
    for (int cellId = blockIdx.x * blockDim.x + threadIdx.x; cellId < spaceSize; cellId += blockDim.x * gridDim.x)
    {
        xId = cellId % spaceWidth;
        yId = cellId - xId;
        xLeft = (xId + spaceWidth - 1) % spaceWidth;
        xRight = (xId + 1) % spaceWidth;
        yBottom = (yId + spaceSize - spaceWidth) % spaceSize;
        yAbove = (yId + spaceWidth) % spaceSize;
        numAlive = lifeSpace_1[xLeft + yId] + lifeSpace_1[xLeft + yAbove] + lifeSpace_1[xLeft + yBottom] + lifeSpace_1[xId + yAbove] + lifeSpace_1[xId + yBottom] + lifeSpace_1[xRight + yId] + lifeSpace_1[xRight + yAbove] + lifeSpace_1[xRight + yBottom];

        lifeSpace_2[cellId] = numAlive == 3 || (numAlive == 2 && lifeSpace_1[cellId]) ? 1 : 0;
    }
}

void GameofLife_simpleGPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t numThread)
{
    size_t spaceSize = spaceHeight * spaceWidth;
    size_t mallocSize = spaceSize * sizeof(uchar);
    unsigned short numBlock = (unsigned short)min(spaceSize / numThread, (size_t)65535);
    uchar *d_lifeSpace_1;
    uchar *d_lifeSpace_2;
    cudaEvent_t start, stop;
    float elapsed_time = 0.0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaMalloc((void **)&d_lifeSpace_1, mallocSize);
    cudaMalloc((void **)&d_lifeSpace_2, mallocSize);
    cudaMemcpy(d_lifeSpace_1, lifeSpace_1, mallocSize, cudaMemcpyHostToDevice);

    dim3 grid_dim(numBlock, 1, 1);
    dim3 block_dim(numThread, 1, 1);
    cudaEventRecord(start, 0);
    for (size_t i = 0; i < iteration; ++i)
    {
        simpleLifeKernel<<<grid_dim, block_dim>>>(d_lifeSpace_1, d_lifeSpace_2, spaceWidth, spaceHeight);
        cudaMemcpy(d_lifeSpace_1, d_lifeSpace_2, mallocSize, cudaMemcpyDeviceToDevice);
        cudaDeviceSynchronize();
    }
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time, start, stop);

    cudaMemcpy(lifeSpace_2, d_lifeSpace_2, mallocSize, cudaMemcpyDeviceToHost);

    cudaFree(d_lifeSpace_1);
    cudaFree(d_lifeSpace_2);
    printf("  grid  dim:  %d, %d, %d.\n", grid_dim.x, grid_dim.y, grid_dim.z);
    printf("  block dim: %d, %d, %d.\n", block_dim.x, block_dim.y, block_dim.z);
    printf("  kernel time: %.5f sec.\n", elapsed_time / 1000.);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

__global__ 
void tilingLifeKernel(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t tileWidth)
{
    extern __shared__ uchar TMP[];

    uint bx = blockIdx.x; 
    uint by = blockIdx.y;
    uint tx = threadIdx.x;
    uint ty = threadIdx.y;
    uint yId = (by * tileWidth + ty) * spaceWidth;
    uint xId = bx * tileWidth + tx;
    size_t spaceSize = spaceHeight * spaceWidth;
    uint xLeft = (xId + spaceWidth - 1) % spaceWidth;
    uint xRight = (xId + 1) % spaceWidth;
    uint yBottom = (yId + spaceSize - spaceWidth) % spaceSize;
    uint yAbove = (yId + spaceWidth) % spaceSize;

    TMP[(ty + 1) * (tileWidth + 2) + tx + 1] = lifeSpace_1[xId + yId];
    if (ty == 0)
        TMP[tx + 1] = lifeSpace_1[xId + yBottom];
    if (ty == tileWidth - 1)
        TMP[(tileWidth + 1) * (tileWidth + 2) + tx + 1] = lifeSpace_1[xId + yAbove];
    if (tx == 0)
        TMP[(ty + 1) * (tileWidth + 2) + 0] = lifeSpace_1[xLeft + yId];
    if (tx == tileWidth - 1)
        TMP[(ty + 1) * (tileWidth + 2) + tileWidth + 1] = lifeSpace_1[xRight + yId];
    if (ty == 0 && tx == 0)
        TMP[0] = lifeSpace_1[xLeft + yBottom];
    if (ty == tileWidth - 1 && tx == 0)
        TMP[(tileWidth + 1) * (tileWidth + 2) + 0] = lifeSpace_1[xLeft + yAbove];
    if (ty == tileWidth - 1 && tx == tileWidth - 1)
        TMP[(tileWidth + 1) * (tileWidth + 2) + tileWidth + 1] = lifeSpace_1[xRight + yAbove];
    if (ty == 0 && tx == tileWidth - 1)
        TMP[tileWidth + 1] = lifeSpace_1[xRight + yBottom];
    __syncthreads();
    uchar numAlive = TMP[ty * (tileWidth + 2) + tx] + TMP[ty * (tileWidth + 2) + tx + 1] + TMP[ty * (tileWidth + 2) + tx + 2] + TMP[(ty + 1) * (tileWidth + 2) + tx] + TMP[(ty + 1) * (tileWidth + 2) + tx + 2] + TMP[(ty + 2) * (tileWidth + 2) + tx] + TMP[(ty + 2) * (tileWidth + 2) + tx + 1] + TMP[(ty + 2) * (tileWidth + 2) + tx + 2];
    
    lifeSpace_2[xId + yId] = numAlive == 3 || (numAlive == 2 && TMP[(ty + 1) * (tileWidth + 2) + tx + 1]) ? 1 : 0;
}

void GameofLife_tilingGPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t numThread)
{
    uint tileWidth = numThread;

    size_t spaceSize = spaceHeight * spaceWidth;
    size_t mallocSize = spaceSize * sizeof(uchar);
    uchar *d_lifeSpace_1;
    uchar *d_lifeSpace_2;

    cudaEvent_t start, stop;
    float elapsed_time = 0.0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaMalloc((void **)&d_lifeSpace_1, mallocSize);
    cudaMalloc((void **)&d_lifeSpace_2, mallocSize);
    cudaMemcpy(d_lifeSpace_1, lifeSpace_1, mallocSize, cudaMemcpyHostToDevice);

    assert((spaceWidth % tileWidth == 0) && (spaceHeight % tileWidth == 0));
    uint widthPhase = spaceWidth / tileWidth;
    uint heightPhase = spaceHeight / tileWidth;

    dim3 grid_dim(widthPhase, heightPhase, 1);
    dim3 block_dim(tileWidth, tileWidth, 1);
    cudaEventRecord(start, 0);

    uint shardSize = (tileWidth + 2) * (tileWidth + 2) * sizeof(uchar);
    for (size_t i = 0; i < iteration; ++i)
    {
        tilingLifeKernel<<<grid_dim, block_dim, shardSize>>>(d_lifeSpace_1, d_lifeSpace_2, spaceWidth, spaceHeight, tileWidth);
        cudaMemcpy(d_lifeSpace_1, d_lifeSpace_2, mallocSize, cudaMemcpyDeviceToDevice);
        cudaDeviceSynchronize();
    }
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed_time, start, stop);

    cudaMemcpy(lifeSpace_2, d_lifeSpace_2, mallocSize, cudaMemcpyDeviceToHost);

    cudaFree(d_lifeSpace_1);
    cudaFree(d_lifeSpace_2);
    printf("  grid  dim:  %d, %d, %d.\n", grid_dim.x, grid_dim.y, grid_dim.z);
    printf("  block dim: %d, %d, %d.\n", block_dim.x, block_dim.y, block_dim.z);
    printf("  kernel time: %.5f sec.\n", elapsed_time / 1000.);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

}


// __global__ void bitEncodeKernel(const uchar *lifeSpace, uchar *encodedLifeSpace, size_t encodedLifeSpaceSize)
// {

//     for (size_t encodedId = blockIdx.x * blockDim.x + threadIdx.x; encodedId < encodedLifeSpaceSize; encodedId += blockDim.x * gridDim.x)
//     {

//         size_t cellId = encodedId << 3;

//         uchar result = lifeSpace[cellId] << 7 | lifeSpace[cellId + 1] << 6 | lifeSpace[cellId + 2] << 5 | lifeSpace[cellId + 3] << 4 | lifeSpace[cellId + 4] << 3 | lifeSpace[cellId + 5] << 2 | lifeSpace[cellId + 6] << 1 | lifeSpace[cellId + 7];

//         encodedLifeSpace[encodedId] = result;
//     }
// }

// void runBitEncodeKernel(const uchar *lifeSpace, uchar *encodedLifeSpace, size_t spaceWidth, size_t spaceHeight)
// {

//     assert(spaceWidth % 8 == 0);
//     size_t encodedSpaceWidth = spaceWidth / 8;
//     size_t encodedLifeSpaceSize = encodedSpaceWidth * spaceHeight;

//     size_t threadsCount = 256;
//     assert(encodedLifeSpaceSize % threadsCount == 0);
//     unsigned short blocksCount = (unsigned short)min(encodedLifeSpaceSize / threadsCount, (size_t)65535);

//     bitEncodeKernel<<<blocksCount, threadsCount>>>(lifeSpace, encodedLifeSpace, encodedLifeSpaceSize);
//     checkCudaErrors(cudaDeviceSynchronize());
// }

// __global__ void bitDecodeKernel(const uchar* encodedLifeSpace, uchar* lifeSpace, size_t encodedLifeSpaceSize) {

//     for (size_t i = blockIdx.x * blockDim.x + threadIdx.x; i < encodedLifeSpaceSize; i += blockDim.x * gridDim.x) 
//     {
//         uint cellId = i << 3;
//         uchar i_state = encodedLifeSpace[i];

//         lifeSpace[cellId] = i_state >> 7;
//         lifeSpace[cellId + 1] = (i_state >> 6) & 1;
//         lifeSpace[cellId + 2] = (i_state >> 5) & 1;
//         lifeSpace[cellId + 3] = (i_state >> 4) & 1;
//         lifeSpace[cellId + 4] = (i_state >> 3) & 1;
//         lifeSpace[cellId + 5] = (i_state >> 2) & 1;
//         lifeSpace[cellId + 6] = (i_state >> 1) & 1;
//         lifeSpace[cellId + 7] = i_state & 1;
//     }
// }

// void runBitDecodeKernel(const uchar *encodedLifeSpace, uchar *lifeSpace, size_t spaceWidth, size_t spaceHeight)
// {
//     assert(spaceWidth % 8 == 0);
//     size_t encodedSpaceWidth = spaceWidth / 8;
//     size_t encodedLifeSpaceSize = encodedSpaceWidth * spaceHeight;

//     size_t threadsCount = 256;
//     assert(encodedLifeSpaceSize % threadsCount == 0);
//     unsigned short blocksCount = (unsigned short)min(encodedLifeSpaceSize / threadsCount, (size_t)65535);

//     bitDecodeKernel<<<blocksCount, threadsCount>>>(encodedLifeSpace, lifeSpace, encodedLifeSpaceSize);
//     checkCudaErrors(cudaDeviceSynchronize());
// }