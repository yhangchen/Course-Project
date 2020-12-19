/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code examines the usage of synchronization constructs. 
*/

#include <omp.h>
#include <stdio.h> 

#define n 10

int main(int argc, char *argv[]){
  int  nt, tid;
  int  i, a[n], sum = 0;

#pragma omp parallel private(i, tid)
  {
    nt = omp_get_num_threads();
    tid = omp_get_thread_num();
    for (i = tid; i < n; i += nt) {
      a[i] = i;
    }
#pragma omp barrier
#pragma omp single
    for (i = 0; i < n; i++) {
      printf(" tid = %d/%d, a[%d] = %d\n", tid, nt, i, a[i]);
    }
#pragma omp for
    for (i = 0; i < n; i++) {
#pragma omp critical (summation)
      sum += a[i];
    }
#pragma omp master
    printf(" sum = %d\n", sum);
  }
  return 0;
}
