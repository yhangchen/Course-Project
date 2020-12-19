/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to calculate pi with simd construct. 
*/

#include <omp.h>
#include <stdio.h>
#include <math.h>

#define PI25DT 3.141592653589793238462643

#pragma omp declare simd 
double f(double x) {
  return (16.0*(x-1.0)/(x*x*x*x-2.0*x*x*x+4.0*x-4.0));
}

int main(int argc, char *argv[]){
  int    n, i;
  double pi, h, x, t0, t1;

  t0 = omp_get_wtime();
  n = 1024000000;
  h = 1.0 / (double) n;
  pi = 0.0;
#pragma omp parallel for simd private(x) reduction(+:pi) 
  for (i = 0; i < n; i++) {
    x = h * ((double)i + 0.5);
    pi += h * f(x);
  }
  t1 = omp_get_wtime();
  printf(" Number of threads = %d\n", omp_get_max_threads());
  printf(" pi is approximately %.20f\n", pi);
  printf(" Error is %.16f\n", fabs(pi-PI25DT));
  printf(" Wall clock time = %f\n", t1-t0);

  return 0;
}
