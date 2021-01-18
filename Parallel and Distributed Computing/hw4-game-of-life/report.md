# CUDA's Implementation of Game of Life

**Author Yihang Chen (陈奕行)**

## 1. Introduction

The Game of Life, also known simply as Life, is a cellular automaton devised by the British mathematician John Horton Conway in 1970.[1] It is a zero-player game, meaning that its evolution is determined by its initial state, requiring no further input. One interacts with the Game of Life by creating an initial configuration and observing how it evolves. 

The universe of the Game of Life is an infinite, two-dimensional orthogonal grid of square cells, each of which is in one of two possible states, live or dead, (or populated and unpopulated, respectively). Every cell interacts with its eight neighbours, which are the cells that are horizontally, vertically, or diagonally adjacent. At each step in time, the following transitions occur:

1. Any live cell with two or three live neighbours survives.
2. Any dead cell with three live neighbours becomes a live cell.
3. All other live cells die in the next generation. Similarly, all other dead cells stay dead.

For simplicity, we consider the cyclic space. This means that left neighbor of the leftmost column is the rightmost column and vice versa. We use the following variables in our programmes.

```c
size_t spaceWidth;
size_t spaceHeight;
size_t spaceSize; 
usigned char *lifeSpace; //unsigned char, size of "spaceSize".
```

Note that "unsigned char" is used to store the input to decrease the storage of the game's space. 

## 2. CPU implementation

The CPU implementation is coded in "CudaLife.h".

```c
static inline void GameofLife_CPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t num_CPU_threads)
```

For lifeSpace_1 to lifeSpace_2, we use the following codes:

```c
#pragma omp for
        for (size_t cellId = 0; cellId < spaceSize; cellId++)
        {
            xId = cellId % spaceWidth;
            yId = cellId - xId;
            xLeft = (xId + spaceWidth - 1) % spaceWidth;
            xRight = (xId + 1) % spaceWidth;
            yBottom = (yId + spaceSize - spaceWidth) % spaceSize;
            yAbove = (yId + spaceWidth) % spaceSize;

            numAlive = lifeSpace_1[xLeft + yId] + lifeSpace_1[xLeft + yAbove] + lifeSpace_1[xLeft + yBottom] + lifeSpace_1[xId + yId] + lifeSpace_1[xId + yAbove] + lifeSpace_1[xId + yBottom] + lifeSpace_1[xRight + yId] + lifeSpace_1[xRight + yAbove] + lifeSpace_1[xRight + yBottom];

            lifeSpace_2[cellId] = (numAlive == 3 || (lifeSpace_1[cellId] && numAlive == 2)) ? 1 : 0;
        }
```

Note yId is not the vertical coordinate. In this way, we can save some floating point arithmetic. Then we use

```c
#pragma omp for
        for (size_t cellId = 0; cellId < spaceSize; cellId++)
            lifeSpace_1[cellId] = lifeSpace_2[cellId];
```

to move the data linked to lifeSpace_2 to lifeSpace_1.

### 

```c
__global__ void simpleLifeKernel(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight)
```
We use OpenMP to perform parallel programming. 

### 3. Simple GPU implementation

We require the number of GPU Blocks fewer than 65535 (maximum of unsigned short).

```c
unsigned short numBlock = (unsigned short)min(spaceSize / numThread, (size_t)65535);
```

Since spaceSize / numThread might be larger than 65535.

In rhe kernel

```c
__global__ 
void simpleLifeKernel(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight)
```

we need to incorporate for loop 

```c
for (int cellId = blockIdx.x * blockDim.x + threadIdx.x; cellId < spaceSize; cellId += blockDim.x * gridDim.x)
```

Even when spaceSize % numThread != 0, this implementation is still correct. Yet different threads might be sightly imbalanced. 

We call the kernel by 

```c
void GameofLife_simpleGPU(uchar *lifeSpace_1, uchar *lifeSpace_2, size_t spaceWidth, size_t spaceHeight, size_t iteration, size_t numThread)
```

### 4. Tiled GPU implementation
We use shared memory to load data for each block to reduce the time consumption on memory access. Specifically, we use ''blockSize = tileWidth''. In each block, we load a (tileWidth + 2) * (tileWidth + 2) array to store the information of each tile and two additional rows and columns around it. Specifially, we use
```c
tilingLifeKernel<<<grid_dim, block_dim, shardSize>>>(d_lifeSpace_1, d_lifeSpace_2, spaceWidth, spaceHeight, tileWidth); // inside function GameofLife_tilingGPU.
...
extern __shared__ uchar TMP[]; // inside function tilingLifeKernel.
``` 
to launch the kernel.

The following "if" clauses are adopted to store the information for each block.
```c
TMP[(ty + 1) * (tileWidth + 2) + tx + 1] = lifeSpace_1[xId + yId];
if (ty == 0)
    TMP[tx + 1] = lifeSpace_1[xId + yBottom];
if (ty == tileWidth - 1)
    TMP[(tileWidth + 1) * (tileWidth + 2) + tx + 1] = lifeSpace_1[xId + yAbove];
if (tx == 0)
    TMP[(ty + 1) * (tileWidth + 2) + 0] = lifeSpace_1[xLeft + yId];
if (tx == tileWidth - 1)
    TMP[(ty + 1) * (tileWidth + 2) + tileWidth + 1] = lifeSpace_1[xRight + yId];
if (ty == 0 && tx == 0)
    TMP[0] = lifeSpace_1[xLeft + yBottom];
if (ty == tileWidth - 1 && tx == 0)
    TMP[(tileWidth + 1) * (tileWidth + 2) + 0] = lifeSpace_1[xLeft + yAbove];
if (ty == tileWidth - 1 && tx == tileWidth - 1)
    TMP[(tileWidth + 1) * (tileWidth + 2) + tileWidth + 1] = lifeSpace_1[xRight + yAbove];
if (ty == 0 && tx == tileWidth - 1)
    TMP[tileWidth + 1] = lifeSpace_1[xRight + yBottom];
__syncthreads();

```
Note that we must have spaceWidth / tileWidth and spaceHeight / tildWidth is integer.  

### 5. Numerical results
The program can be executed by:
```bash
module add cuda
make
sbatch run.slurm
```

Input Size
| 1st param | 2nd param | 3rd param | 4th param | 5th param |
|---|---|---|---|---|
| width | height | iteration | num_CPU_threads | num_GPU_threads |

In the following, we take num_GPU_threads = 1 / 4 / 8. We need to change 
```bash
#SBATCH --ntasks-per-node=8 
```
in run.slurm.

| width | height | iteration | CPU threads | OPENMP TIME |
|---|---|---|---|---|
|400|400|100|1|  0.89570  |
|400|2000|100|1|  4.47279 |
|2000|2000|100|1| 22.39981 |
|400|400|100|4| 0.75427   |
|400|2000|100|4|  3.77121 |
|2000|2000|100|4| 18.76387 |
|400|400|100|8|  0.41985  |
|400|2000|100|8| 1.99657  |
|2000|2000|100|8| 10.02208 |

It acheives approximately 2x accelerations for the number of threads to be 8 compared to serial implementation. Clearly, such acceleration is not satisfactory. 

In the following, we take num_GPU_threads = 8 / 16. The kernel is inside the bracket.
| width | height | iteration | GPU threads| simple GPU time  | tiling GPU time |
|---|---|---|---|---|---|
|400|400|100|8|0.29551(0.00477)|0.00337(0.00208)|
|400|2000|100|8|0.22029(0.01849)|0.00699(0.00567)|
|2000|2000|100|8|0.35675(0.08472)|0.02985(0.02292)|
|400|400|100|16|0.22529(0.00309)|0.00285(0.00212)|
|400|2000|100|16|0.31690(0.01042)|0.01255(0.00573)|
|2000|2000|100|16|0.29920(0.04581)|0.02881(0.02544)|

Clearly, it acheives 400 ~ 500 x acceleration by tiled GPU implementation. We can find that the kernel time for simple GPU consumes approximate twice the time by tiling GPU.

It achieves 10x - 50x acceleration for simple GPU and 500x acceleration for tiling GPU, compared with 8 threads OPENMP, which implied over 1000x acceleration compared with serial codes on CPU.

We also consider a much larger grid, and only test GPU implementation:
| width | height | iteration | GPU threads| simple GPU time  | tiling GPU time |
|---|---|---|---|---|---|
|40000|40000|100|8|28.42349(26.26473)|16.21500(7.81866)|
|40000|40000|100|16|16.39350(14.21215)|16.62610(7.63598)|
|40000|40000|100|32|10.75161(8.59912)|16.24306(7.81866)|
|40000|40000|100|64|7.58502(9.76549)|9.21603(0.74369)|

Note we at most launch 65535 blocks in simple GPU. At this time, the tiling GPU does not have much advatage, since a tile is of size 627x627 B (tileWidth=64), approximately 400kB, which is larger than the default 48kB shared memory. Hence, the memory overflows and the acceleration disappears.