/*
 * $Smake: g++ -o %F %f
 *
 * THIS PROGRAM CONTAINS ERRORS
 *
 * This program contains several errors that are here to be found using
 * the gdb debugger and valgrind memory checker.
 *
 * This code is supposed to create a small square matrix A and vector x
 * and then form the product y = Ax.
 */

#include <cstdio>

int main(int argc, char *argv[])
{
    const int N = 5;

    // allocate contiguous memory for A

    double** a = new double* [N];   // memory for array of row pointers
    a[0] = new double [N * N];      // memory for array
    for (int i = 1; i < N; i++)
    {
        a[i] = &a[0][i * N];        // assign row pointers
    }

    // allocate memory for x and y vectors

    double* x = new double [N];
    double* y = new double [N];

    // initalize array and vector

    for (int i = 0; i <= N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            a[i][j] = 10.0 * (i + 1) + j;
        }
        x[i] = 5.0 * (i + 1);
    }

    // compute product y = Ax

    for (int i = 0; i < N; i++)
    {
        y[i] = 0.0;
        for (int j = 0; j <= N; j++)
        {
            y[i] += a[i][j] * x[j];
        }
    }

    // display results

    for (int i = 0; i < N; i++)
    {
        printf("[");
        for (int j = 0; j < N; j++)
        {
            printf(" %4.1f", a[i][j]);
        }
        printf("] [ %4.1f] ", x[i]);
        printf(i == N / 2 ? "=" : " ");
        printf(" [ %6.1f]\n", y[i]);
    }

    // release memory

    delete [] a;
    delete [] a[0];

    return 0;
}
