/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  Testing the threadprivate data type of OpenMP. 
  The code is modified from the LLNL omp tutorial. 
*/

#include <omp.h>
#include <stdio.h>

int a = 0, b = 0;
#pragma omp threadprivate(a)

int compare(int x)
{
  static int n = 2;
#pragma omp threadprivate(n)
  if (n < x)
    n = x;
  return n;
}

int main(int argc, char *argv[])
{
  int tid, c, d;
  omp_set_dynamic(0);
  printf("1st Parallel Region:\n");
#pragma omp parallel num_threads(4) private(tid, b, c, d)
  {
    tid = omp_get_thread_num();
    a = tid + 1;
    b = tid + 2;
    c = compare(a);
    d = compare(b);
    printf("Thread %d: a,b,c,d = %d %d %d %d\n", tid, a, b, c, d);
  }
  printf("Serial Region: a,b = %d %d\n", a, b);
  printf("2nd Parallel Region:\n");
#pragma omp parallel num_threads(4) private(tid, c, d)
  {
    tid = omp_get_thread_num();
    c = compare(a + 2);
    d = compare(b + 3);
    printf("Thread %d: a,b,c,d = %d %d %d %d\n", tid, a, b, c, d);
  }

  return 0;
}
