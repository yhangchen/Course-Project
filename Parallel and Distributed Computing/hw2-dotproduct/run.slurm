#!/bin/bash

#SBATCH -o job_%j.out
#SBATCH --partition=cpu
#SBATCH -J mpi-ex4
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 10:00

# you'll need mpich to run this job
mpiexec -n 4 ./hw2 1000 0
mpiexec -n 4 ./hw2 1000 1
mpiexec -n 4 ./hw2 1000 2
mpiexec -n 4 ./hw2 1000 3
mpiexec -n 4 ./hw2 1000 4
mpiexec -n 8 ./hw2 1000 0
mpiexec -n 8 ./hw2 1000 1
mpiexec -n 8 ./hw2 1000 2
mpiexec -n 8 ./hw2 1000 3
mpiexec -n 8 ./hw2 1000 4
mpiexec -n 16 ./hw2 1000 0
mpiexec -n 16 ./hw2 1000 1
mpiexec -n 16 ./hw2 1000 2
mpiexec -n 16 ./hw2 1000 3
mpiexec -n 16 ./hw2 1000 4
mpiexec -n 4 ./hw2 729 0
mpiexec -n 4 ./hw2 729 1
mpiexec -n 4 ./hw2 729 2
mpiexec -n 4 ./hw2 729 3
mpiexec -n 4 ./hw2 729 4
mpiexec -n 8 ./hw2 729 0
mpiexec -n 8 ./hw2 729 1
mpiexec -n 8 ./hw2 729 2
mpiexec -n 8 ./hw2 729 3
mpiexec -n 8 ./hw2 729 4
mpiexec -n 16 ./hw2 729 0
mpiexec -n 16 ./hw2 729 1
mpiexec -n 16 ./hw2 729 2
mpiexec -n 16 ./hw2 729 3
mpiexec -n 16 ./hw2 729 4
mpiexec -n 4 ./hw2 10005 0
mpiexec -n 4 ./hw2 10005 1
mpiexec -n 4 ./hw2 10005 2
mpiexec -n 4 ./hw2 10005 3
mpiexec -n 4 ./hw2 10005 4
mpiexec -n 8 ./hw2 10005 0
mpiexec -n 8 ./hw2 10005 1
mpiexec -n 8 ./hw2 10005 2
mpiexec -n 8 ./hw2 10005 3
mpiexec -n 8 ./hw2 10005 4
mpiexec -n 16 ./hw2 10005 0
mpiexec -n 16 ./hw2 10005 1
mpiexec -n 16 ./hw2 10005 2
mpiexec -n 16 ./hw2 10005 3
mpiexec -n 16 ./hw2 10005 4