#!/bin/bash

#SBATCH -o job_%j.out
#SBATCH --partition=cpu
#SBATCH -J mpi-ex1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 10:00

# you'll need mpich to run this job
mpiexec -n 4 ./bcast_reduce
#mpiexec -n 4 ./scatter_gather
#mpiexec -n 4 ./hello2
#mpiexec -n 1 ./cpi
#mpiexec -n 2 ./cpi
#mpiexec -n 4 ./cpi
#mpiexec -n 8 ./cpi

