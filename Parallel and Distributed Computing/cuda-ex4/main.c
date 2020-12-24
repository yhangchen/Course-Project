/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to compute matrix multiplicaiton. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mat_mul.h"

int main(int argc, char **argv){
  
  int width = 512;  // default value
  if (argc > 1) width = atoi(argv[1]); // user-specified value
  printf("\nMatrix width: %d.\n", width);

  srand(time(NULL));
  float* M = rand_mat(width, width);
  float* N = rand_mat(width, width);

  float* cpu_P = raw_mat(width, width);
  float* gpu_P = raw_mat(width, width);

  long long cpu_start_time = start_timer();
  cpu_mat_mul(M, N, cpu_P, width);
  long long cpu_time = stop_timer(cpu_start_time, "CPU");

  long long gpu_start_time = start_timer();
  gpu_mat_mul(M, N, gpu_P, width);
  long long gpu_time = stop_timer(gpu_start_time, "GPU");

  // Check the correctness of the GPU results
  int num_wrong = 0;
  for (int i = 0; i < width * width; i++) {
    if (fabs(cpu_P[i] - gpu_P[i]) > 0.000001) num_wrong++;
  }
	
  // Report the correctness results
  if (num_wrong) printf("GPU %d / %d values incorrect\n", num_wrong, N);
  else           printf("GPU all values correct\n");

}
