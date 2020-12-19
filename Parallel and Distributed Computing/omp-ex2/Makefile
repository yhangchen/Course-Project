CC = gcc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = hello cpi private schedule

all: ${EXE}

hello: omp_hello.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp

cpi: omp_cpi.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm

cpi2: omp_cpi2.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm

cpi3: omp_cpi3.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm

private: omp_private.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp

private2: omp_private2.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp

schedule: omp_schedule.c
	$(CC) -o $@ $^ $(FLAGS) -fopenmp -lm

clean:
	rm -f $(OBJ) $(EXE)

