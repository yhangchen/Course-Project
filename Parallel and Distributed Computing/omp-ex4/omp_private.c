/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code examines the usage of private-like clauses. 
*/

#include <omp.h>
#include <stdio.h> 

int main(int argc, char *argv[]){
  int  i, k;

  printf("Test:\n");
  k = 100;
  printf(" first k = %d\n", k);
#pragma omp parallel for private(k) 
  for (i = 0; i < 10; i++) {
    k = i;
    printf("  private k = %d\n", k);
  }
  printf(" last k = %d\n", k);

  return 0;
}
