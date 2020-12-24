#include <cuda.h>
#include <stdio.h>
#include <math.h>

#define BLOCK_WIDTH 16 
#define TILE_WIDTH  BLOCK_WIDTH

extern "C" void gpu_mat_mul(float* h_M, float* h_N, float* h_P, int width);

__global__
void gpu_mat_mul_kernel(float* M, float* N, float* P, int width){

  __shared__ float Mds[TILE_WIDTH][TILE_WIDTH];
  __shared__ float Nds[TILE_WIDTH][TILE_WIDTH];

  int bx = blockIdx.x; 
  int by = blockIdx.y;
  int tx = threadIdx.x;
  int ty = threadIdx.y;

  // Identify the row and column of the P element to work on
  // Each thread works on an element of P
  int Row = by * TILE_WIDTH + ty;
  int Col = bx * TILE_WIDTH + tx;

  float sum = 0;
  int phase_num = width/TILE_WIDTH;

  // Each thread loads 'Row'th row of M and 'Col'th column of N
  for (int ph = 0; ph < phase_num; ph++) {    

    // Collaboratively load data into shared memory
    Mds[ty][tx] = M[Row * width + ph * TILE_WIDTH + tx];   
    Nds[ty][tx] = N[(ph * TILE_WIDTH + ty) * width + Col];

    __syncthreads();
    for (int k = 0; k < TILE_WIDTH; k++) { 
      sum += Mds[ty][k] * Nds[k][tx];
    }
    __syncthreads();
  }

  P[Row * width + Col] = sum;
}

void gpu_mat_mul(float* h_M, float* h_N, float* h_P, int width) {
  float *d_M, *d_N, *d_P;

  size_t size_of_float = sizeof(float);
  size_t size_M = width * width * size_of_float;
  size_t size_N = width * width * size_of_float;
  size_t size_P = width * width * size_of_float;

  cudaMalloc((void**)&d_M, size_M);
  cudaMalloc((void**)&d_N, size_N);
  cudaMalloc((void**)&d_P, size_P);
    
  cudaMemcpy(d_M, h_M, size_M, cudaMemcpyHostToDevice);
  cudaMemcpy(d_N, h_N, size_N, cudaMemcpyHostToDevice);

  cudaEvent_t start, stop;
  float elapsed_time = 0.0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  dim3 grid_dim(width/BLOCK_WIDTH, width/BLOCK_WIDTH, 1);
  dim3 block_dim(BLOCK_WIDTH, BLOCK_WIDTH, 1);
  gpu_mat_mul_kernel<<<grid_dim, block_dim>>>(d_M, d_N, d_P, width);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  cudaMemcpy(h_P, d_P, size_P, cudaMemcpyDeviceToHost);
    
  // Free device memory for M, N, P
  cudaFree(d_M);
  cudaFree(d_N);
  cudaFree(d_P);

  cudaEventElapsedTime(&elapsed_time, start, stop);
    
  printf("  grid  dim:  %d, %d, %d.\n", grid_dim.x, grid_dim.y, grid_dim.z);
  printf("  block dim: %d, %d, %d.\n", block_dim.x, block_dim.y, block_dim.z);
  printf("  kernel time: %.5f sec\n", elapsed_time / 1000);

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}


