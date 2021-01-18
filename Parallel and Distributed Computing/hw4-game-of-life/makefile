CC = gcc
CFLAGS = -std=c99 -O3 -fopenmp

NVCC = nvcc
NVCC_FLAGS = --gpu-architecture=sm_50 -std=c++11 -O3 -Wno-deprecated-gpu-targets

LIBRARIES = -L/mnt/lustrefs/softwares/cuda/10.2/lib64 -lcudart -lm -fopenmp

game_of_life: main.o CudaLife.o 
	$(CC) $^ -o $@ $(LIBRARIES)

main.o: main.c
	$(CC) $(CFLAGS) -c $^ -o $@

CudaLife.o: CudaLife.cu
	$(NVCC) $(NVCC_FLAGS) -c $^ -o $@

clean:
	rm -f *.o game_of_life
	
