CC = gcc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = cpi cpi_task cpi_simd

all: ${EXE}


%: omp_%.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm -std=c99
cpi_simd: omp_cpi_simd.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -fopt-info-vec -lm

clean:
	rm -f $(OBJ) $(EXE)
