#!/bin/bash

#SBATCH -o job_%j.out
#SBATCH --partition=cpu
#SBATCH -J mpi-ex1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH -t 10:00

# you'll need mpich to run this job
mpiexec -n 2 ./exchange
#mpiexec -n 2 ./exchange2
