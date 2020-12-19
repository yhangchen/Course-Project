/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  Testing the nested parallelism of OpenMP. 
  The code is taken from the official omp-4.5 example set. 
*/

#include <omp.h>
#include <stdio.h>
int main (void)
{
  omp_set_dynamic(0);
#pragma omp parallel num_threads(2)
  {

    omp_set_nested(1);
#pragma omp parallel num_threads(3)
    {
#pragma omp single
      printf ("Inner: num_thds=%d\n", omp_get_num_threads());
    }

#pragma omp barrier
    omp_set_nested(0);
#pragma omp parallel num_threads(3)
    {
#pragma omp single
      printf ("Inner: num_thds=%d\n", omp_get_num_threads());
    }

#pragma omp barrier
#pragma omp single
    printf ("Outer: num_thds=%d\n", omp_get_num_threads());

  }
  return 0;
}
