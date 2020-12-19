#ifndef PIC2GREY_H
#define PIC2GREY_H
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#define CHANNELS 3

void cpu_pic2grey(unsigned char * Pout, unsigned char * Pin, int width, int height);

void gpu_pic2grey(unsigned char * Pout, unsigned char * Pin, int width, int height);

// Returns the current time in microseconds
extern "C" inline long long start_timer() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000000 + tv.tv_usec;
}

// Prints the time elapsed since the specified time
inline long long stop_timer(long long start_time, char *name) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  long long end_time = tv.tv_sec * 1000000 + tv.tv_usec;
  printf("%s: %.5f sec\n", name, ((float) (end_time - start_time)) / (1000 * 1000));
  return end_time - start_time;
}


#endif
