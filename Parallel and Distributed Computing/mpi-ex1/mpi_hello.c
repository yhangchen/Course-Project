/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows the most basic usage of MPI. 
*/

#include <mpi.h> // mpi header file
#include <stdio.h> // standard I/O

int main(int argc, char *argv[]){
  int size, rank;
  
  MPI_Init(&argc, &argv); // initialize MPI 
  MPI_Comm_size(MPI_COMM_WORLD, &size); // get num of procs
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get my rank
  printf("From process %d out of %d, Hello World!\n", rank, size);
  if (rank == 0) printf("That's all, folks!\n");
  MPI_Finalize(); // done with MPI
  return 0;
} 

