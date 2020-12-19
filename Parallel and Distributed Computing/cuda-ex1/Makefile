CC = nvcc 
FLAGS =
OBJ = *.o
EXE = hello 

all: ${EXE}

hello: cuda_hello.cu
	$(CC) -o $@ $^
clean:
	rm -f $(OBJ) $(EXE)

