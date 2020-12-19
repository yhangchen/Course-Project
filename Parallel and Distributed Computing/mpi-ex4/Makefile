CC = mpicc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = split group

all: ${EXE}

split: mpi_split.c
	$(CC) -o $@ $^ $(FLAGS)
group: mpi_group.c
	$(CC) -o $@ $^ $(FLAGS)
clean:
	rm -f $(OBJ) $(EXE)

