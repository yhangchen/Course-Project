#!/bin/bash

#SBATCH -o job_%j.out
#SBATCH --partition=cpu
#SBATCH -J omp-ex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 10:00

./cpi 8 1 111111
./cpi 8 2 111111
./cpi 8 3 111111
./cpi 8 4 111111
./cpi 8 5 111111
./cpi 8 6 111111
./cpi 8 7 111111
./cpi 8 8 111111
./cpi 8 9 111111
./cpi 8 10 111111
./cpi_task 8 111111
./cpi_simd 8 111111