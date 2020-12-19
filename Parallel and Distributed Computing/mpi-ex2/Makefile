CC = mpicc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = exchange exchange2 

all: ${EXE}

exchange: mpi_exchange.c
	$(CC) -o $@ $^ $(FLAGS)
exchange2: mpi_exchange2.c
	$(CC) -o $@ $^ $(FLAGS)
clean:
	rm -f $(OBJ) $(EXE)

