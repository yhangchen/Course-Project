/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows the consepts of MPI_Group. 
  It is modified from:
  https://github.com/wesleykendall/mpitutorial/blob/gh-pages
  /tutorials/introduction-to-groups-and-communicators/code/groups.c
*/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv){
  int       rank, size, sub_rank, sub_size;
  MPI_Group group, sub_group;
  MPI_Comm  sub_comm;
  const int ranks[4] = {2, 3, 5, 7};

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Comm_group(MPI_COMM_WORLD, &group);
  MPI_Group_incl(group, 4, ranks, &sub_group);

  MPI_Comm_create(MPI_COMM_WORLD, sub_group, &sub_comm);

  if (sub_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(sub_comm, &sub_rank);
    MPI_Comm_size(sub_comm, &sub_size);
  } else {
    sub_rank = -1;
    sub_size = -1;
  }

  printf("World rank/size: %d/%d --- Sub rank/size: %d/%d\n", rank, size, sub_rank, sub_size);

  MPI_Group_free(&group);
  MPI_Group_free(&sub_group);

  if (sub_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&sub_comm);
  }
  MPI_Finalize();
  return 0;
}
