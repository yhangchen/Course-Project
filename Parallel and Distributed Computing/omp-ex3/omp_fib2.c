/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  Compute Fibonacci Numbers with the task construct faster.
*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

uint64_t fib(int n);

uint64_t *a;

int main(int argc, char *argv[])
{
  uint64_t res;
  int n = 30;
  double t0, t1;

  a = malloc((n + 1) * sizeof(uint64_t));
  memset(a, 0, (n + 1) * sizeof(uint64_t));

  t0 = omp_get_wtime();
#pragma omp parallel
  {
#pragma omp single
    res = fib(n);
  }
  t1 = omp_get_wtime();

  printf(" Number of threads = %d\n", omp_get_max_threads());
  printf(" The %d-th Fibonacci Number is %lu\n", n, res);
  printf(" Wall clock time = %f\n", t1 - t0);

  free(a);
  return 0;
}

uint64_t fib(int n)
{
  uint64_t x, y, res;

  if (n < 2)
    res = n;
  else
  {
#pragma omp task shared(x)
    x = fib(n - 1);
#pragma omp task shared(y)
    y = fib(n - 2);
#pragma omp taskwait
    res = x + y;
  }
  a[n] = res;

  return a[n];
}
