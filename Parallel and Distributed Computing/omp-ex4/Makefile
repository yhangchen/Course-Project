CC = gcc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = hello cpi private schedule sync fib fib2 nested threadprivate copyprivate cpi3_simd cpi3_nosimd

all: ${EXE}

%: omp_%.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm

cpi3_simd: omp_cpi3_simd.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -fopt-info-vec

cpi3_nosimd: omp_cpi3_nosimd.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -fopt-info-vec -fno-tree-vectorize

clean:
	rm -f $(OBJ) $(EXE)

