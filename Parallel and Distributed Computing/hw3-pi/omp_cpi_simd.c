#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI25DT 3.141592653589793238462643

#pragma omp declare simd
int main(int argc, char *argv[])
{
    int n, i, nthreads = 8;
    double pi, h, x, t0, t1;

    t0 = omp_get_wtime();
    if (argc > 1)
    {
        nthreads = atoi(argv[1]);
    }
    n = 111111;
    if (argc > 2)
    {
        n = atoi(argv[2]);
    }
    h = 1.0 / (double)n;
    pi = 0.0;
    omp_set_num_threads(nthreads);
#pragma omp parallel for simd private(x) reduction(+ \
                                                   : pi)
    for (i = 0; i < n; i++)
    {
        x = h * ((double)i + 0.5);
        pi += 4.0 * h * sqrt(1. - x * x);
    }
    t1 = omp_get_wtime();
    printf(" Number of threads = %d\n", omp_get_max_threads());
    printf(" pi is approximately %.20f\n", pi);
    printf(" Error is %.4e\n", fabs(pi - PI25DT));
    printf(" Wall clock time = %f\n", t1 - t0);

    return 0;
}
