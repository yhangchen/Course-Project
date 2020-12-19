/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to use MPI_Bcast/MPI_Reduce. 
*/

#include "mpi.h"
#include <stdio.h>

#define MAX_P 16 // Maximum number of processes
#define ROOT   0 // Rank of the root process

int main(int argc, char **argv){
  int rank, size;
  int data, data_reduce;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Check if the number of processes is too large */
  if ( size > MAX_P ) {
    printf("Number of processes %d is too large. About MPI!\n", size);
    MPI_Abort(MPI_COMM_WORLD, 911);
  }

  /* Set the data on root process */
  if (rank == ROOT) {
    data = 100;
    printf("On root proc %d, data to be broadcast: %d\n", rank, data);
  }

  /* Broadcast the data to everybody */
  MPI_Bcast(&data, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  /* Everybody shows the broadcast data */
  printf("On proc %d, data after broadcasting = %d\n", rank, data);
  MPI_Barrier(MPI_COMM_WORLD);

  /* Modify the data to reduce */
  data = data + rank;
  printf("On proc %d, data to be reduced = %d\n", rank, data);

  /* Reduce everybody's data to root */
  MPI_Reduce(&data, &data_reduce, 1, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);

  /* On root, show the reduced data */
  if (rank == ROOT) {
    printf("On root proc %d, data reduced: %d\n", rank, data_reduce);
  }

  MPI_Finalize();
  return 0;
}