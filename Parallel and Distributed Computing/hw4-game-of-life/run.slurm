#!/bin/bash

#SBATCH -o job_%j_%N.out
#SBATCH --partition=gpu
#SBATCH -J gpu-ex1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -t 10:00

# test of CPU, not we need to modify --ntasks-per-node
# ./game_of_life 400 400 100 1 8
# ./game_of_life 400 2000 100 1 8
# ./game_of_life 2000 2000 100 1 8

# ./game_of_life 400 400 100 4 16
# ./game_of_life 400 2000 100 4 16
# ./game_of_life 2000 2000 100 4 16

# ./game_of_life 400 400 100 8 16
# ./game_of_life 400 2000 100 8 16
# ./game_of_life 2000 2000 100 8 16


# test of GPU
# ./game_of_life 400 400 100 8 8
# ./game_of_life 400 2000 100 8 8
# ./game_of_life 2000 2000 100 8 8

# ./game_of_life 400 400 100 8 16
# ./game_of_life 400 2000 100 8 16
# ./game_of_life 2000 2000 100 8 16

# the following codes cannot be tested on CPU.
# ./game_of_life 40000 40000 100 8 8
# ./game_of_life 40000 40000 100 8 16
# ./game_of_life 40000 40000 100 8 32
# ./game_of_life 40000 40000 100 8 64
