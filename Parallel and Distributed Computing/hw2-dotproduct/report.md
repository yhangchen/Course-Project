# MPI 点乘

## 1.程序设计

首先，我们需要考虑进程个数size无法整除$n$ 的情形，我们对每个进程分配${\rm local\_n}=[n/{\rm size}]$ 个元素，并将剩余的元素在0号进程里进行计算。

### 1.1 Scatter Gather

利用scatter/gather，程序可以很容易地写出。核心部分如下：

```
MPI_Scatter(a_cut, local_n, MPI_DOUBLE, sub_a, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Gather(&sub_res, 1, MPI_DOUBLE, sub_ress, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
```

其中a_cut的大小为local_n。

注意的是Scatter和Gather是在每个进程里运行的。我们将这种方法编号设为1.

### 1.2 Scatterv Gather

利用Scatterv可以避免特殊划分，其官方文档如下：

```c
int MPIAPI MPI_Scatterv(
  _In_  void         *sendbuf,
  _In_  int          *sendcounts,
  _In_  int          *displs,
        MPI_Datatype sendtype,
  _Out_ void         *recvbuf,
        int          recvcount,
        MPI_Datatype recvtype,
        int          root,
        MPI_Comm     comm
);
```

可以看出只需要指定每一部分的长度（sendcounts）即可，然后没有间断或重复（displs可以由sendcounts累加得到），即

```c
      for (int i = 0; i < size; i++)
      {
        sendcounts[i] = local_n;
      }
      sendcounts[0] += n - cut;
      displs[0] = 0;
      for (int i = 1; i < size; i++)
      {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
      }

```

稍微修改就得到结果，我们将这种方法编号设为0.

### 1.2 Send/Irecv, ISend/Recv, ISend/IRecv.

首先建立一个二维数组mat_a/mat_b，将a_cut分为${\rm size}\times n$的二维矩阵。然后，将a_cut[i]从第0号进程发送至第i号进程。第i号进程接收之后，将计算结果发回第0号进程。有以下几点值得注意：

（1）第0号进程不能free(mat_a); free(mat_b)；否则有一些可能还没有传出去。

（2）使用Irecv时一定加Wait，不使用时不加Wait。

特别地，我们将这三种方法分别设为2，3，4.

## 2.数值结果

```
mpiexec -n nprocs ./hw2 n choice
```

其中nproc为进程数目，$n$为向量长度，choice为选择的方法。以下是在数院机器上运行的结果

| 方法\进程个数(n=1000) | 4         | 8         | 16        | 串行      |
| :-------------------: | --------- | --------- | --------- | --------- |
|           0           | 0.007676s | 0.012279s | 0.032876s | 0.000006s |
|           1           | 0.007551s | 0.019448s | 0.045187s | ~         |
|           2           | 0.007929s | 0.019737s | 0.042239s | ~         |
|           3           | 0.007824s | 0.020254s | 0.025717s | ~         |
|           4           | 0.007840s | 0.004252s | 0.046052s | ~         |

| 方法\进程个数(n=729) | 4         | 8         | 16        | 串行      |
| :------------------: | --------- | --------- | --------- | --------- |
|          0           | 0.007570s | 0.020642s | 0.031019s | 0.000004s |
|          1           | 0.010789s | 0.019684s | 0.032153s | ~         |
|          2           | 0.008110s | 0.021254s | 0.043337s | ~         |
|          3           | 0.007645s | 0.011976s | 0.027537s | ~         |
|          4           | 0.007572s | 0.014229s | 0.023909s | ~         |

注意到实际上随着线程增加，反而速度变慢，可以看出主要时间都花在了通信上面。余下其不同方法差距其实不大。