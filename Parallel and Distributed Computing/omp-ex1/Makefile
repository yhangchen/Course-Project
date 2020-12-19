CC = gcc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = hello cpi

all: ${EXE}

hello: omp_hello.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp

cpi: omp_cpi.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm

clean:
	rm -f $(OBJ) $(EXE)

