CC = gcc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = hello cpi private schedule sync fib fib2 nested 

all: ${EXE}

%: omp_%.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm

clean:
	rm -f $(OBJ) $(EXE)

