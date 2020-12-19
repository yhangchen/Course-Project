/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows the most basic usage of OpenMP. 
*/

#include <omp.h>   // omp header file
#include <stdio.h> // standard I/O

int main(int argc, char *argv[])
{
  int nthreads, tid;
  double t0, t1;
  omp_set_num_threads(4);
  t0 = omp_get_wtime();
#pragma omp parallel private(tid)
  {
    nthreads = omp_get_num_threads(); // get num of threads
    tid = omp_get_thread_num();       // get my thread id
    printf("From thread %d out of %d, Hello World!\n", tid, nthreads);
  }
  t1 = omp_get_wtime();
  printf("Time elapsed is %f.\nThat's all, folks!\n", t1 - t0);
  return 0;
}
