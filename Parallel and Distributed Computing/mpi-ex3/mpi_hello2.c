/* 
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows the most basic usage of MPI. 
*/

#include "mpi.h" // mpi header file
#include <stdio.h> // standard I/O

int main(int argc, char *argv[]){
  int size, rank, len, ver, sver;
  char name[MPI_MAX_PROCESSOR_NAME];
  double t0, t1;

  MPI_Init(&argc, &argv); // initialize MPI 
  MPI_Comm_size(MPI_COMM_WORLD, &size); // get num of procs
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // get my rank  
  MPI_Get_processor_name(name, &len); // get node name
  MPI_Get_version(&ver, &sver); // get mpi version

  if ( size > 16 ) {
    printf("Number of processes %d is too large. About MPI!\n", size);
    MPI_Abort(MPI_COMM_WORLD, 911); // about mpi for large size
  }

  t0 = MPI_Wtime(); // tick
  printf("On %s, from process %d of %d: Hello World!\n", name, rank, size);
  fflush(stdout); // flush standard output
  MPI_Barrier(MPI_COMM_WORLD); // barrier for the flush
  t1 = MPI_Wtime(); // tock
  
  if (rank == 0) printf("MPI version is %1d.%1d. Time elapsed is %f.\nThat's all, folks!\n", ver, sver, t1-t0);

  MPI_Finalize(); // done with MPI
  return 0;
} 
