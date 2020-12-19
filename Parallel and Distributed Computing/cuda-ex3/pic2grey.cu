#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include "pic2grey.h"

#define BLOCK_DIM 16

__global__
void gpu_pic2grey_kernel(unsigned char * Pout, unsigned char * Pin, int width, int height) {
  int Col = threadIdx.x + blockIdx.x * blockDim.x; 
  int Row = threadIdx.y + blockIdx.y * blockDim.y; 

  if (Col < width && Row < height) {
    int grey_offset = Row * width + Col;
    int rgb_offset = grey_offset * CHANNELS;
    unsigned char r = Pin[rgb_offset + 0]; 
    unsigned char g = Pin[rgb_offset + 1]; 
    unsigned char b = Pin[rgb_offset + 2]; 
    Pout[grey_offset] = 0.21f * r + 0.71f * g + 0.07f * b;
  } 
}

void gpu_pic2grey(unsigned char *h_out, unsigned char *h_in, int width, int height){
  unsigned char* d_out;
  unsigned char* d_in;

  size_t size_in = width * height * CHANNELS * sizeof(unsigned char);
  size_t size_out = width * height * sizeof(unsigned char);

  // Allocates object in the device global memory
  cudaMalloc((void **)&d_out, size_out);
  cudaMalloc((void **)&d_in, size_in);

  // Memory data transfer from host to device
  cudaMemcpy(d_in, h_in, size_in, cudaMemcpyHostToDevice);

  cudaEvent_t start, stop;
  float elapsed_time = 0.0;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);
  
  // Invoke kernel to do the computation on device
  dim3 grid_dim(ceil(width/(float)(BLOCK_DIM)),ceil(height/(float)(BLOCK_DIM)));
  dim3 block_dim(BLOCK_DIM,BLOCK_DIM);
  gpu_pic2grey_kernel<<<grid_dim, block_dim>>>(d_out, d_in, width, height);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  // Transfer back the result from d_Pout to h_Pout
  cudaMemcpy(h_out, d_out, size_out, cudaMemcpyDeviceToHost);
  
  // Free device memory for Pout, Pin
  cudaFree(d_out);
  cudaFree(d_in);
  
  cudaEventElapsedTime(&elapsed_time, start, stop);
  
  printf("  grid  dim:  %d, %d, %d.\n", grid_dim.x, grid_dim.y, grid_dim.z);
  printf("  block dim: %d, %d, %d.\n", block_dim.x, block_dim.y, block_dim.z);
  printf("  kernel time: %.5f sec.\n", elapsed_time / 1000.);

  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}
