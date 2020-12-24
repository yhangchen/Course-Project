#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI25DT 3.141592653589793238462643

int main(int argc, char *argv[])
{
  int nthreads = 8, tid, n, i,
      choice = 1;
  double pi, h, x, t0, t1;
  if (argc > 1)
  {
    nthreads = atoi(argv[1]);
  }
  if (argc > 2)
  {
    choice = atoi(argv[2]);
  }
  t0 = omp_get_wtime();
  n = 11111111;
  if (argc > 3)
  {
    n = atoi(argv[3]);
  }
  h = 1.0 / (double)n;
  pi = 0.0;
  omp_set_num_threads(nthreads);
  switch (choice)
  {
  case 1:
  {
#pragma omp parallel default(shared) private(tid, i, x) reduction(+ \
                                                                  : pi)
    {
      nthreads = omp_get_num_threads();
      tid = omp_get_thread_num();
      for (i = tid + 1; i <= n; i += nthreads)
      {
        x = h * ((double)i - 0.5);
        pi += 4.0 * h * sqrt(1. - x * x);
      }
    }
    break;
  }
  case 2:
  {
#pragma omp parallel for default(shared) private(x) reduction(+ \
                                                              : pi)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  case 3:
  {
#pragma omp parallel default(shared)
#pragma omp for private(x) reduction(+ \
                                     : pi)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  case 4:
  {
#pragma omp parallel for default(shared) private(x) schedule(static, 1) reduction(+ \
                                                                                  : pi)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  case 5:
  {
#pragma omp parallel for default(shared) private(x) schedule(static) reduction(+ \
                                                                               : pi)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  case 6:
  {
#pragma omp parallel for default(shared) private(x) schedule(dynamic) reduction(+ \
                                                                                : pi)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  case 7:
  {
#pragma omp parallel for default(shared) private(x) schedule(guided) reduction(+ \
                                                                               : pi)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  case 8:
  {
#pragma omp parallel for default(shared) private(x) schedule(auto) reduction(+ \
                                                                             : pi)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  case 9:
  {
#pragma omp parallel for default(shared) private(x)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
#pragma omp critical(summation)
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  case 10:
  {
#pragma omp parallel for default(shared) private(x)
    for (i = 1; i <= n; i++)
    {
      x = h * ((double)i - 0.5);
#pragma omp atomic
      pi += 4.0 * h * sqrt(1. - x * x);
    }
    break;
  }
  default:
    printf("Wrong Choice! Enter int from 1 to 10.\n");
    return 0;
  }

  t1 = omp_get_wtime();
  printf(" Number of threads = %d\n", nthreads);
  printf(" Parallel Choice = %d\n", choice);
  printf(" pi is approximately %.16f\n", pi);
  printf(" Error is %.4e\n", fabs(pi - PI25DT));
  printf(" Wall clock time = %f\n", t1 - t0);

  return 0;
}
