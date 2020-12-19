#ifndef VEC_ADD_H
#define VEC_ADD_H

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

void vec_add(float* h_A, float* h_B, float* h_C, int n);

// Returns a randomized vector containing n elements
static inline float *rand_vec(int n) {
  float *vec = (float *) malloc(n * sizeof(float));
  if (vec == NULL) { 
    printf("Error allocating CPU memory");
      exit(1);
  }
  for (int i = 0; i < n; i++) vec[i] = (float) (rand() % 100);
  return vec;
}

// Returns a raw vector containing n elements
static inline float *raw_vec(int n) {
  float *vec = (float *) malloc(n * sizeof(float));
  if (vec == NULL) { 
    printf("Error allocating CPU memory");
      exit(1);
  }
  return vec;
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




