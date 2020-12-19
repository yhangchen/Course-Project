/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how MPI_Comm_split works. 
  It is modified from:
  https://github.com/wesleykendall/mpitutorial/tree/gh-pages
  /tutorials/introduction-to-groups-and-communicators/code
*/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv){
  int rank, size, color, sub_rank, sub_size;
  MPI_Comm sub_comm;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  color = rank/4;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub_comm);

  MPI_Comm_rank(sub_comm, &sub_rank);
  MPI_Comm_size(sub_comm, &sub_size);

  printf("World rank/size: %d/%d -- Sub rank/size: %d/%d\n", rank, size, sub_rank, sub_size);

  MPI_Comm_free(&sub_comm);
  MPI_Finalize();
  return 0;
}
