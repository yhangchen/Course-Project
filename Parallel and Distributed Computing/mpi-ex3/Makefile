default: all

CC = mpicc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = bcast_reduce scatter_gather hello2 cpi

all: ${EXE}

bcast_reduce: mpi_bcast_reduce.c
	$(CC) -o $@ $^ $(FLAGS) 
scatter_gather: mpi_scatter_gather.c
	$(CC) -o $@ $^ $(FLAGS) 
hello2: mpi_hello2.c
	$(CC) -o $@ $^ $(FLAGS)
cpi: mpi_cpi.c
	$(CC) -o $@ $^ $(FLAGS) -lm
clean:
	rm -f $(OBJ) $(EXE)

