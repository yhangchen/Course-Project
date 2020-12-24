#ifndef MAT_MUL_H
#define MAT_MUL_H

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

void cpu_mat_mul(float* h_A, float* h_B, float* h_C, int width);

void gpu_mat_mul(float* h_A, float* h_B, float* h_C, int width);

// Returns a randomized matrix containing mxn elements
static inline float *rand_mat(int m, int n) {
  float *mat = (float *) malloc(m * n * sizeof(float));
  if (mat == NULL) { 
    printf("Error allocating CPU memory");
    exit(1);
  }
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      mat[i * n + j] = (float)(rand() % 100);
    }
  }
  return mat;
}

// Returns a raw matrix containing n elements
static inline float *raw_mat(int m, int n) {
  float *mat = (float *) malloc(m * n * sizeof(float));
  if (mat == NULL) { 
    printf("Error allocating CPU memory");
    exit(1);
  }
  return mat;
}

// Returns the current time in microseconds
static inline long long start_timer() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000000 + tv.tv_usec;
}

// Prints the time elapsed since the specified time
static inline long long stop_timer(long long start_time, char *name) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  long long end_time = tv.tv_sec * 1000000 + tv.tv_usec;
  printf("%s: %.5f sec\n", name, ((float) (end_time - start_time)) / (1000 * 1000));
  return end_time - start_time;
}

#endif
