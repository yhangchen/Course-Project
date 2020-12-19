/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to use MPI_Scatter/MPI_Gather. 
  It is modified from:
  http://mpi.deino.net/mpi_functions/MPI_Scatter.html
*/

#include "mpi.h"
#include <stdio.h>

#define MAX_P 16 // Maximum number of processes
#define N      2 // Number of data per process
#define ROOT   0 // Rank of the root process

int main(int argc, char **argv){
  int rank, size, i, j;
  int data_root[MAX_P][N], data[N];

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
    printf("On root proc %d, data to be scattered:\n", rank);
    for ( j=0; j<size; j++) {
      for ( i=0; i<N; i++ ) {
        data_root[j][i] = j+i;
        printf("(%d, %d) = %d  ", j, i, data_root[j][i]);
      }
      printf("\n");
    }
  }

  /* Scatter the data to everybody */
  MPI_Scatter(&data_root[0][0], N, MPI_INT, &data[0], N, MPI_INT, ROOT, MPI_COMM_WORLD);

  /* Everybody shows the scatterd data */
  for (i=0; i<N; i++) {
    printf("On proc %d, scattered data[%d] = %d\n", rank, i, data[i]);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /* Gather everybody's data to root */
  MPI_Gather(&data[0], N, MPI_INT, &data_root[0][0], N, MPI_INT, ROOT, MPI_COMM_WORLD);

  /* On root, show the gathered data */
  if (rank == ROOT) {
    printf("On root proc %d, data gathered:\n", rank);
    for ( j=0; j<size; j++) {
      for ( i=0; i<N; i++ ) {
        printf("(%d, %d) = %d  ", j, i, data_root[j][i]);
      }
      printf("\n");
    }
  }

  MPI_Finalize();
  return 0;
}