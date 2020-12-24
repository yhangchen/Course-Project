# OPENMP 计算 $\pi$

## 1. 积分公式

$$
\pi = 4\int_0^1  \sqrt{1-x^2} {\rm d}x
$$

由于我们关心不同并行算法的效果，我们采用矩形积分公式并使用等距格点。

## 2. 并行方法

### 2.1 对parallel和for的修改（omp_cpi.c）

1. parallel
2. parallel for
3. parallel + for
4. parallel for + static, 1
5. parallel for + static
6. parallel for + dynamic
7. parallel for + guided
8. parallel for + auto
9. parallel for + critical
10. parallel for + atomic

三个参数分别为，线程数目，并行策略（1~10），以及n大小。

### 2.2 使用Task结构 （omp_cpi_task.c）

具体而言，我们考虑了一段积分double cpi(int start, int end, int n)，其中n代表grid point总数，以上函数返回从第start个节点至第end-1个节点上计算出来的积分。在end-start<n/1024时候，我们直接计算这一段的积分，以减少调度开销。

两个参数分别为，线程数目，以及n大小。

### 2.3 使用SIMD （omp_cpi_simd.c）

两个参数分别为，线程数目，以及n大小。

## 3.数值结果

nthreads = 8, n = 111111.

| method                  | Error      | Time     |
| ----------------------- | ---------- | -------- |
| parallel                | 9.2998e-09 | 0.000289 |
| parallel for            | 9.2998e-09 | 0.000338 |
| parallel + for          | 9.2998e-09 | 0.001098 |
| parallel for + static,1 | 9.2998e-09 | 0.000300 |
| parallel for + static   | 9.2998e-09 | 0.000287 |
| parallel for + dynamic  | 9.2998e-09 | 0.008532 |
| parallel for + guided   | 9.2998e-09 | 0.000303 |
| parallel for + auto     | 9.2998e-09 | 0.000264 |
| parallel for + critical | 9.2998e-09 | 0.040159 |
| parallel for + atomic   | 9.2998e-09 | 0.010112 |
| task                    | 9.2998e-09 | 0.005748 |
| SIMD                    | 9.2998e-09 | 0.000216 |

nthreads = 2, n = 111111.

| method                  | Error      | Time     |
| ----------------------- | ---------- | -------- |
| parallel                | 9.2998e-09 | 0.000394 |
| parallel for            | 9.2998e-09 | 0.000411 |
| parallel + for          | 9.2998e-09 | 0.000414 |
| parallel for + static,1 | 9.2998e-09 | 0.000384 |
| parallel for + static   | 9.2998e-09 | 0.000376 |
| parallel for + dynamic  | 9.2998e-09 | 0.009497 |
| parallel for + guided   | 9.2998e-09 | 0.000365 |
| parallel for + auto     | 9.2998e-09 | 0.000352 |
| parallel for + critical | 9.2998e-09 | 0.016658 |
| parallel for + atomic   | 9.2998e-09 | 0.008754 |
| task                    | 9.2998e-09 | 0.006241 |
| SIMD                    | 9.2998e-09 | 0.000326 |

我们使用更大的n。

nthreads = 8, n = 11111111

| method                  | Error      | Time     |
| ----------------------- | ---------- | -------- |
| parallel                | 9.2868e-12 | 0.009028 |
| parallel for            | 9.2868e-12 | 0.008318 |
| parallel + for          | 9.2868e-12 | 0.007723 |
| parallel for + static,1 | 9.2868e-12 | 0.007232 |
| parallel for + static   | 9.2868e-12 | 0.006810 |
| parallel for + dynamic  | 9.2868e-12 | 0.667579 |
| parallel for + guided   | 9.2868e-12 | 0.004352 |
| parallel for + auto     | 9.2868e-12 | 0.004361 |
| parallel for + critical | 9.2868e-12 | 2.219530 |
| parallel for + atomic   | 9.2868e-12 | 0.830800 |
| task                    | 9.2868e-12 | 0.039022 |
| SIMD                    | 9.2868e-12 | 0.010588 |



不难发现使用guided dynamic最为高效。在chunk为1时，由于调度开销的存在，static比dynamic更为高效。parallel 与 for的不同写法不影响效率。SIMD将循环向量化也较为高效。同时，我们发现atomic比critical更为高效，因为其更针对于简单的相加步骤。而task使用了二叉树不断迭代，而且在首末端较近的时候直接计算（避免最后过多的调度开销），效率明显提高。

因而，由于此问题每一步计算较为简单，如何减少调度开销时加速并行计算的核心问题。