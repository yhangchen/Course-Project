#!/bin/bash

#SBATCH -o job_%j_%N.out
#SBATCH --partition=gpu
#SBATCH -J vec_add 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 10:00

./vec_add 102400

./vec_add 409600

./vec_add 819200

./vec_add 1024000
 

