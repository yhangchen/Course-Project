/* 
  Foundations of Parallel and Distributed Computing, Falls 2019.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to do vector add C = A + B in cuda. 
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vec_add.h"

void cpu_vec_add(float *h_A, float *h_B, float *h_C, int n) { 
  for (int i = 0; i < n; i++) 
    h_C[i] = h_A[i] + h_B[i];
}

int main(int argc, char **argv){
  
  // Determine the vector length
  int N = 1024000;  // default value
  if (argc > 1) N = atoi(argv[1]); // user-specified value
  printf("\nVector length: %d.\n", N);

  // Memory allocation for vectors A and B
  srand(time(NULL));
  float* h_A = rand_vec(N);
  float* h_B = rand_vec(N);

  // Memory allocation for vector C
  float* cpu_C = raw_vec(N);
  float* gpu_C = raw_vec(N);
    
  long long cpu_start_time = start_timer();
  cpu_vec_add(h_A, h_B, cpu_C, N);
  long long cpu_time = stop_timer(cpu_start_time, "CPU");

  long long gpu_start_time = start_timer();
  vec_add(h_A, h_B, gpu_C, N);
  long long gpu_time = stop_timer(gpu_start_time, "GPU");

  // Check the correctness of the GPU results
  int num_wrong = 0;
  for (int i = 0; i < N; i++) {
    if (fabs(cpu_C[i] - gpu_C[i]) > 0.000001) num_wrong++;
  }
	
  // Report the correctness results
  if (num_wrong) printf("%d / %d values incorrect.\n", num_wrong, N);
  else           printf("All values correct.\n");
}

