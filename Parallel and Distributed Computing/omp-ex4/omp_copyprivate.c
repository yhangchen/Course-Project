/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  Testing the copyin and copyprivate clause. 
*/

#include <omp.h>
#include <stdio.h>

int counter = 0;
#pragma omp threadprivate(counter)

int increment_counter()
{
  return (++counter);
}

int main(int argc, char *argv[])
{
  int tid, c;
  omp_set_dynamic(0);
  omp_set_num_threads(4);
  printf("1st Parallel Region:\n");
#pragma omp parallel private(tid, c)
  {
    tid = omp_get_thread_num();
#pragma omp single copyprivate(counter)
    counter = 50 + tid;
    c = increment_counter();
    printf("ThreadId: %d, count = %d\n", tid, c);
#pragma omp barrier
    counter = 100 + tid;
    c = increment_counter();
    printf("ThreadId: %d, count = %d\n", tid, c);
  }
  printf("2nd Parallel Region:\n");
#pragma omp parallel private(tid, c) copyin(counter)
  {
    tid = omp_get_thread_num();
    c = increment_counter();
    printf("ThreadId: %d, count = %d\n", tid, c);
  }
  return 0;
}
