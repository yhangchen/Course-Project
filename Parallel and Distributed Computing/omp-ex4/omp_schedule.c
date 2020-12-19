/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code examines the usage of schedule clauses. 
*/

#include <omp.h>
#include <stdio.h> 
#include <math.h> 

#define n 12

int main(int argc, char *argv[]){
  int    tid, i, j;
  double sum = 0.;

  printf("tid 0\ttid 1\ttid 2\ttid 3\n");
#pragma omp parallel for num_threads(4) \
  reduction(+:sum) schedule(static)
  for (i = 0; i < n; i++) {
    for (j = 0; j < 100000; j++) {
      sum += sin(i+j);
    }
    tid = omp_get_thread_num();
    if (tid == 0) {
      printf("i = %d\n", i);
    } else if (tid == 1) {
      printf("\ti = %d\n", i);
    } else if (tid == 2) {
      printf("\t\ti = %d\n", i);
    } else {
      printf("\t\t\ti = %d\n", i);
    }
  }

  return 0;
}
