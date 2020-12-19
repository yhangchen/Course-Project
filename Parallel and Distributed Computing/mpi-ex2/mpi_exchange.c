/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to use MPI send/recv for data exchange. 
*/
#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
  int a, b, size, rank;

  MPI_Init(&argc, &argv); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

  if (rank == 0) {
    a = -1;
    MPI_Send(&a, 1, MPI_INT, 1-rank, 0, MPI_COMM_WORLD);
    printf("Process %d sent token %d to process %d\n", rank, a, 1-rank);
    MPI_Recv(&b, 1, MPI_INT, 1-rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Process %d received token %d from process %d\n", rank, b, 1-rank);
  } else if (rank == 1) {
    a = 1;
    MPI_Recv(&b, 1, MPI_INT, 1-rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&a, 1, MPI_INT, 1-rank, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0;
}
