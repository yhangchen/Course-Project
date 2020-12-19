#!/bin/bash

#SBATCH -o job_%j.out
#SBATCH --partition=cpu
#SBATCH -J mpi-ex4
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 10:00

# you'll need mpich to run this job
mpiexec -n 16 ./split
#mpiexec -n 8 ./group
