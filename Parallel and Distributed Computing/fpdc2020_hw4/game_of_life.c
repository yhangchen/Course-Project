/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This is a naive implement of Conwey's Game of Life.
  Version 0.1 (date: 04/16/2020)
*/

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define L1   2000
#define L2   2000
#define LT   100
#define SEED 0
#define LS   (L1*L2)
#define index(i,L) ((i<0)?(i+L):((i>=L)?(i-L):(i)))

double get_walltime() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return (double) (tp.tv_sec + tp.tv_usec*1e-6); 
}

int main(int argc, char *argv[]) {
  int    i, t, x, y, x0, y0;
  int    ineighbor, neighbors, population0, population;
  int    *world, *result;
  double time;

  // initialize data
  world  = malloc(LS * sizeof(int));
  result = malloc(LS * sizeof(int));
  srand(SEED);
  population0 = 0;
  for (i = 0; i < LS; i++) {
    world[i] = rand()%2;
    if (world[i]==1) population0++;
  }

  // naive implementation
  time = get_walltime();
  for (t = 0; t < LT; t++) {
    population = 0;
    for (i = 0; i < LS; i++) {
      x = i % L1;
      y = i / L1;
      neighbors = 0;
      for (y0 = y - 1; y0 <= y + 1; y0++) {
        for (x0 = x - 1; x0 <= x + 1; x0++) {
          ineighbor = index(y0, L2) * L1 + index(x0, L1);
          if (world[ineighbor] != 0) neighbors++;
        }
      }
      if (world[i] != 0) neighbors--;
      if (neighbors == 3 || (neighbors == 2 && (world[i] != 0))) {
        result[i] = 1;
        population++;
      } else {
        result[i] = 0;
      }
    }
    for (i = 0; i < LS; i++) {
      world[i] = result[i];
    }
  }
  time = get_walltime() - time;

  // report results
  printf("World size: %d x %d, total generations: %d\n", L1, L2, LT);
  printf("Population is changed from %d to %d\n", population0, population);
  printf("Wall time: %f\n", time);

  // cleanup data
  free(world);  world  = NULL;
  free(result); result = NULL;
  
  return 0; 
}
