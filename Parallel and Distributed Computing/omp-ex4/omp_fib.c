/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  Compute Fibonacci Numbers with the task construct. 
*/

#include <omp.h>
#include <stdio.h>
#include <math.h>

int fib(int n);

int main(int argc, char *argv[]){
  int     res, n = 30;
  double  t0, t1;

  t0 = omp_get_wtime();
#pragma omp parallel
  {
#pragma omp single
    res = fib(n);
  }
  t1 = omp_get_wtime();

  printf(" Number of threads = %d\n", omp_get_max_threads());
  printf(" The %d-th Fibonacci Number is %d\n", n, res);
  printf(" Wall clock time = %f\n", t1-t0);

  return 0;
}

int fib(int n)
{
  int x, y;

  if (n < 2)  return n;
#pragma omp task shared(x)
  x = fib(n-1);
#pragma omp task shared(y)
  y = fib(n-2);
#pragma omp taskwait
  return(x+y);
}
