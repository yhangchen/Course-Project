#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>



const double pi = acos(-1.0);



int reverse_bit(int i, int n) {
    int ans = 0;
    for (int j = 0; j < n; ++j)
    {
        ans = (ans << 1) | (i & 1);
        i >>= 1;
    }
    return ans;
}

double f(double t) {
    return exp(-t*t/10) * (sin(2*t) + 2*cos(4*t) + 0.4*sin(t)*sin(50*t));
}

void FFT(double complex *x, double complex *y, int n, int choice) {
    int N = 1 << n;
    for (int i = 0; i < N; ++i){
        int ind = reverse_bit(i, n);
        y[i] = x[ind];
    }
    for (int i = 1; i <= n; ++i) {
        int j = 1 << i;
        for (int s = 0; s < N; s += j) {
            for (int k = 0; k < j / 2; ++k) {
                double complex w = cos(2 * pi * k / j) - choice * sin(2 * pi * k / j) * _Complex_I;
                double complex u = y[s + k];
                double complex v = y[s + k + j / 2];
                y[s + k] = u + w * v;
                y[s + k + j / 2] = u - w * v;
            }
        }
    }
    if (choice < 0){
        for (int i = 0; i <= N; ++i){
            y[i] = y[i] / N;
        }
    }
}

int main(void){
    const int n = 8, m = 6;
    int N = 1 << n;
    double complex z[N];
    double complex z1[N];
    double complex z2[N];
    for (int k = 0; k < N; ++k){
        z[k] = f(2 * k * pi / N);
    }
    FFT(z, z1, n, 1);
    // char st0[100];
    // sprintf(st0, "freq.csv", m);
    // frequency of signal
    // FILE *fp0 = NULL;
    // fp0 = fopen(st0, "w");

    // for (int k = 0; k < N; ++k){
    //     fprintf(fp0, "%g, %g\n", creal(z1[k]), cimag(z1[k]));
    // }
    // fclose(fp0);
    // fp0 = NULL;

    for (int k = m; k < N - m; ++k){
        z1[k] = 0;
    }

    FFT(z1, z2, n, -1);

    char st[100];
    sprintf(st, "out_%d.csv", m);

    FILE *fp = NULL;
    fp = fopen(st, "w");

    for (int k = 0; k < N; ++k){
        fprintf(fp, "%g, %g, %g, %g\n", creal(z[k]), cimag(z[k]), creal(z2[k]), cimag(z2[k]));
    }
    fclose(fp);
    fp = NULL;

    return 0;
}


