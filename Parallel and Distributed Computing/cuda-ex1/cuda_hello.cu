/* 
  Foundations of Parallel and Distributed Computing, Falls 2019.
  Instructor: Prof. Chao Yang @ Peking University.
  This code is a pilot code for cuda. 
*/
#include <cuda.h>
#include <stdio.h>

__global__ void mykernel(void) {
    printf("Hello World from GPU by thread %d!\n", threadIdx.x);
}

int main(void) {
    // Use 1 thread
    mykernel<<<1,1>>>();
    // Flush the output
    cudaDeviceSynchronize();
    // Use 4 threads
    mykernel<<<1,4>>>();
    // Flush the output
    cudaDeviceSynchronize();
    return 0;
}
