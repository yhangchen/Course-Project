/* 
  Foundations of Parallel and Distributed Computing, Falls 2019.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to do vector add in cuda. 
*/

#include <stdio.h>
#include <math.h>
#include <cuda.h>

extern "C" void vec_add(float *h_A, float *h_B, float *h_C, int n);

const int threads_per_block = 256;

// Compute vector sum C = A + B
// Each thread performs one pair-wise addition
__global__
void vec_add_kernel(float *d_A, float *d_B, float *d_C, int n) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i < n) d_C[i] = d_A[i] + d_B[i];
}

void vec_add(float *h_A, float *h_B, float *h_C, int n) {
  int size = n * sizeof(float);
  float *d_A, *d_B, *d_C;
   
  cudaEvent_t start, stop;
  float elapsed_time = 0.0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  //1. Allocate device memory and transfer data from host to device
  cudaMalloc((void **) &d_A, size);
  cudaMalloc((void **) &d_B, size);
  cudaMalloc((void **) &d_C, size);
  cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);
    
  // 2. Configure execution and launch kernel
  dim3 grid_dim(ceil(n/threads_per_block), 1, 1);
  dim3 block_dim(threads_per_block, 1, 1);
  cudaEventRecord(start, 0);
  vec_add_kernel<<<grid_dim, block_dim>>>(d_A, d_B, d_C, n);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_time, start, stop);  
  
  // 3. Transfer result back to host and free device memory 
  cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);
  cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
  
  printf("  grid  dim:  %d, %d, %d.\n", grid_dim.x, grid_dim.y, grid_dim.z);
  printf("  block dim: %d, %d, %d.\n", block_dim.x, block_dim.y, block_dim.z);
  printf("  kernel time: %.5f sec.\n", elapsed_time / 1000.);
  
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}

