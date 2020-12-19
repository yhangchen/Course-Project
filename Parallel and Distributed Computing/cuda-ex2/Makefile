CC = gcc
CFLAGS = -std=c99

NVCC = nvcc
NVCC_FLAGS = -std=c++11 -O3 -Wno-deprecated-gpu-targets

LIBRARIES = -L/mnt/lustrefs/softwares/cuda/10.2/lib64 -lcudart -lm

default: vec_add

vec_add: main.o vec_add.o 
	$(CC) $^ -o $@ $(LIBRARIES)

main.o: main.c
	$(CC) $(CFLAGS) -c $^ -o $@ 

vec_add.o: vec_add.cu
	$(NVCC) $(NVCC_FLAGS) -c $^ -o $@

clean:
	rm -f *.o vec_add
	
