/*
  Foundations of Parallel and Distributed Computing, Falls 2020.
  Instructor: Prof. Chao Yang @ Peking University.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

// Returns a randomized vector containing n elements
double *rand_vec(int n)
{
    double *vec = (double *)malloc(n * sizeof(double));
    if (vec == NULL)
    {
        printf("Error allocating CPU memory");
        exit(1);
    }
    int i;
    for (i = 0; i < n; i++)
    {
        vec[i] = (double)(rand() % 100);
    }
    return vec;
}

// Calculate the result of vector dot
double dot(double *a, double *b, int n)
{
    int i;
    double result = 0.0;
    for (i = 0; i < n; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

int main(int argc, char **argv)
{
    int n = 1024;
    if (argc > 1)
    {
        n = atoi(argv[1]);
    }

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double serial_res = 0.0, parallel_res = 0.0;
    double *a, *b;
    double start, end;
    double serial_start, serial_end;

    // Generate random vector element only in process 0
    if (rank == 0)
    {
        srand(time(NULL));
        a = rand_vec(n);
        b = rand_vec(n);
        serial_start = MPI_Wtime();
        serial_res = dot(a, b, n);
        serial_end = MPI_Wtime();
        start = MPI_Wtime();
    }
    // Write your parallel code below
    int choice = 0;
    int local_n = n / size;
    int cut = local_n * size;
    double a_cut[cut];
    double b_cut[cut];
    int sendcounts[size], displs[size];
    MPI_Request request_scatter[size];
    MPI_Request request_gather[size];
    MPI_Status status;
    if (argc > 2)
    {
        choice = atoi(argv[2]);
    }

    switch (choice)
    {
    case 0:
    {
        if (rank == 0)
        {

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
        }
        int recvcount = local_n;
        if (rank == 0)
        {
            recvcount += n - cut;
        }
        double *sub_a = (double *)malloc(sizeof(double) * recvcount);
        double *sub_b = (double *)malloc(sizeof(double) * recvcount);

        MPI_Scatterv(a, sendcounts, displs, MPI_DOUBLE, sub_a,
                     recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(b, sendcounts, displs, MPI_DOUBLE, sub_b,
                     recvcount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double sub_res = dot(sub_a, sub_b, recvcount);
        double *sub_ress = NULL;
        if (rank == 0)
        {
            sub_ress = (double *)malloc(sizeof(double) * size);
        }
        MPI_Gather(&sub_res, 1, MPI_DOUBLE, sub_ress, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            for (int i = 0; i < size; i++)
            {
                parallel_res += sub_ress[i];
            }

            free(sub_ress);
            free(sub_a);
            free(sub_b);
        }
        break;
    }
    case 1:
    { // Scatter, Gather
        if (rank == 0)
        {

            for (int i = 0; i < cut; i++)
            {
                a_cut[i] = a[i];
                b_cut[i] = b[i];
            }
        }
        double *sub_a = (double *)malloc(sizeof(double) * local_n);
        double *sub_b = (double *)malloc(sizeof(double) * local_n);
        MPI_Scatter(a_cut, local_n, MPI_DOUBLE, sub_a,
                    local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(b_cut, local_n, MPI_DOUBLE, sub_b,
                    local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double sub_res = dot(sub_a, sub_b, local_n);
        double *sub_ress = NULL;
        if (rank == 0)
        {
            sub_ress = (double *)malloc(sizeof(double) * size);
        }
        MPI_Gather(&sub_res, 1, MPI_DOUBLE, sub_ress, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            for (int i = cut; i < n; i++)
            {
                parallel_res += a[i] * b[i];
            }
            for (int i = 0; i < size; i++)
            {
                parallel_res += sub_ress[i];
            }

            free(sub_ress);
            free(sub_a);
            free(sub_b);
        }
        break;
    }
    case 2:
    { // Send, Irecv
        if (rank == 0)
        {
            double mat_a[size][local_n];
            double mat_b[size][local_n];
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < local_n; j++)
                {
                    mat_a[i][j] = a[i * local_n + j];
                    mat_b[i][j] = b[i * local_n + j];
                }
                if (i > 0)
                {
                    MPI_Send(mat_a[i], local_n, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
                    MPI_Send(mat_b[i], local_n, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
                }
            }
            for (int i = cut; i < n; i++)
            {
                parallel_res += a[i] * b[i];
            }
            parallel_res += dot(mat_a[0], mat_b[0], local_n);
            double recv_res;
            for (int i = 1; i < size; i++)
            {
                MPI_Irecv(&recv_res, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request_gather[i]);
                MPI_Wait(&request_gather[i], &status);
                parallel_res += recv_res;
            }
        }
        if (rank != 0)
        {
            double *recv_a = (double *)malloc(sizeof(double) * local_n);
            double *recv_b = (double *)malloc(sizeof(double) * local_n);
            MPI_Irecv(recv_a, local_n, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request_scatter[rank]);
            MPI_Irecv(recv_b, local_n, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request_scatter[rank]);
            MPI_Wait(&request_scatter[rank], &status);
            double send_res = dot(recv_a, recv_b, local_n);
            MPI_Send(&send_res, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
            free(recv_a);
            free(recv_b);
        }
        break;
    }
    case 3:
    { // Isend, Recv.
        if (rank == 0)
        {
            double mat_a[size][local_n];
            double mat_b[size][local_n];
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < local_n; j++)
                {
                    mat_a[i][j] = a[i * local_n + j];
                    mat_b[i][j] = b[i * local_n + j];
                }
                if (i > 0)
                {
                    MPI_Isend(mat_a[i], local_n, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request_scatter[i]);
                    MPI_Isend(mat_b[i], local_n, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request_scatter[i]);
                }
            }
            for (int i = cut; i < n; i++)
            {
                parallel_res += a[i] * b[i];
            }
            parallel_res += dot(mat_a[0], mat_b[0], local_n);
            double recv_res;
            for (int i = 1; i < size; i++)
            {
                MPI_Recv(&recv_res, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                parallel_res += recv_res;
            }
        }
        if (rank != 0)
        {
            double *recv_a = (double *)malloc(sizeof(double) * local_n);
            double *recv_b = (double *)malloc(sizeof(double) * local_n);
            MPI_Recv(recv_a, local_n, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(recv_b, local_n, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            double send_res = dot(recv_a, recv_b, local_n);
            MPI_Isend(&send_res, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request_gather[rank]);
            free(recv_a);
            free(recv_b);
        }
        break;
    }
    case 4:
    { // Isend, Irecv
        if (rank == 0)
        {
            double mat_a[size][local_n];
            double mat_b[size][local_n];
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < local_n; j++)
                {
                    mat_a[i][j] = a[i * local_n + j];
                    mat_b[i][j] = b[i * local_n + j];
                }
                if (i > 0)
                {
                    MPI_Isend(mat_a[i], local_n, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request_scatter[i]);
                    MPI_Isend(mat_b[i], local_n, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request_scatter[i]);
                }
            }
            for (int i = cut; i < n; i++)
            {
                parallel_res += a[i] * b[i];
            }
            parallel_res += dot(mat_a[0], mat_b[0], local_n);
            double recv_res;
            for (int i = 1; i < size; i++)
            {
                MPI_Irecv(&recv_res, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &request_gather[i]);
                MPI_Wait(&request_gather[i], &status);
                parallel_res += recv_res;
            }
        }
        if (rank != 0)
        {
            double *recv_a = (double *)malloc(sizeof(double) * local_n);
            double *recv_b = (double *)malloc(sizeof(double) * local_n);
            MPI_Irecv(recv_a, local_n, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request_scatter[rank]);
            MPI_Irecv(recv_b, local_n, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request_scatter[rank]);
            MPI_Wait(&request_scatter[rank], &status);
            double send_res = dot(recv_a, recv_b, local_n);
            MPI_Isend(&send_res, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &request_gather[rank]);
            free(recv_a);
            free(recv_b);
        }
        break;
    }
    default:
        printf("Wrong Choice! Enter int 0,1,2,3,4.\n");
        return 0;
    }

    // Check correctness and report the execution time
    if (rank == 0)
    {
        end = MPI_Wtime();
        if (fabs(serial_res - parallel_res) > 0.000001)
        {
            printf("Incorrect Parallel Implementation.\n");
        }
        else
        {
            printf(
                "Correct, the result is %f.\n Parallel calculation using %lfs.\n "
                "Serial calculation using %lfs.\n",
                parallel_res, (end - start), (serial_end - serial_start));
        }
        free(a);
        free(b);
    }
    MPI_Finalize();
    return 0;
}