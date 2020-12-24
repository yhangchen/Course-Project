#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI25DT 3.141592653589793238462643

double cpi(int start, int end, int n);

int main(int argc, char *argv[])
{
    int nthreads = 8;
    if (argc > 1)
    {
        nthreads = atoi(argv[1]);
    }
    int n = 111111;
    if (argc > 2)
    {
        n = atoi(argv[2]);
    }    double t0, t1, pi;

    t0 = omp_get_wtime();
    omp_set_num_threads(nthreads);
#pragma omp parallel
    {
#pragma omp single
        pi = cpi(1, n + 1, n);
    }
    t1 = omp_get_wtime();

    printf(" Number of threads = %d\n", omp_get_max_threads());
    printf(" pi is approximately %.16f\n", pi);
    printf(" Error is %.4e\n", fabs(pi - PI25DT));
    printf(" Wall clock time = %f\n", t1 - t0);

    return 0;
}

double cpi(int start, int end, int n)
{
    int middle = (start + end) / 2,
        diff = end - start;
    double h = 1.0 / (double)n;
    double res = 0.0, x, y;
    if (diff < n/1024)
    {
        for (int i = start; i < end; i++)
        {
            x = h * ((double)i - 0.5);
            res += 4.0 * h * sqrt(1. - x * x);
        }
        return res;
    }
    else
    {
#pragma omp task shared(x)
        x = cpi(start, middle, n);
#pragma omp task shared(y)
        y = cpi(middle, end, n);
#pragma omp taskwait
        return x + y;
    }
}