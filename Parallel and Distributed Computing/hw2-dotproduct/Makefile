CC = mpicc
FLAGS = -O3
OBJ = *.o
EXE = hw2

all: ${EXE}

hw2: hw2.c
	$(CC) -o $@ $^ $(FLAGS)
clean:
	rm -f $(OBJ) $(EXE)