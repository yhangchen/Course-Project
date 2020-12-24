CC = gcc
FLAGS = -O3 -Wall
OBJ = *.o
EXE = game_of_life

all: ${EXE}

game_of_life: game_of_life.c makefile
	$(CC) -o $@ game_of_life.c $(FLAGS) 

clean:
	rm -f $(OBJ) $(EXE)
