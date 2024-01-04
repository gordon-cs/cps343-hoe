/*
 * $Smake: g++ -DN=500 -Wall -O3 -funroll-loops -o %F %f
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark ijk, jki, and ikj matrix-matrix products.
 *
 * "Matrix with 1D arrays" version.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>

#if !defined(N)
# define N 500  // default matrix dimension
#endif

const double EPSILON = 1.0E-14;
#define IDX(i,j,stride) ((i)*(stride)+j) // row major
//#define IDX(i,j,stride) ((i)+(stride)*j) // columm major

//----------------------------------------------------------------------------
// Returns the number of seconds since some fixed arbitrary time in the past

double wtime(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) (ts.tv_sec + ts.tv_nsec / 1.0e9);
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using ijk loop order

void matmat_ijk(double* c, double* a, double* b, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            c[IDX(i,j,n)] = 0.0;
            for (int k = 0; k < n; k++)
            {
                c[IDX(i,j,n)] += a[IDX(i,k,n)] * b[IDX(k,j,n)];
            }
        }
    }
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using jki loop order

void matmat_jki(double* c, double* a, double* b, int n)
{
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < n; i++)
        {
            c[IDX(i,j,n)] = 0.0;
        }
        for (int k = 0; k < n; k++)
        {
            for (int i = 0; i < n; i++)
            {
                c[IDX(i,j,n)] += a[IDX(i,k,n)] * b[IDX(k,j,n)];
            }
        }
    }
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using ikj loop order

void matmat_ikj(double* c, double* a, double* b, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            c[IDX(i,j,n)] = 0.0;
        }
        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < n; j++)
            {
                c[IDX(i,j,n)] += a[IDX(i,k,n)] * b[IDX(k,j,n)];
            }
        }
    }
}

//----------------------------------------------------------------------------
// Verify product

int verify(double* d, double* c, int n)
{
    int status = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (fabs(c[IDX(i,j,n)] != d[IDX(i,j,n)]) > EPSILON)
            {
                status++;
                printf("[%d][%d]: c = %f; d = %f\n", i, j,
                        c[IDX(i,j,n)], d[IDX(i,j,n)]);
            }
        }
    }
    return status;
}

//----------------------------------------------------------------------------
// Main program

int main(int argc, char* argv[])
{
    double t0, t1;
    double ijk_time, jki_time, ikj_time;

    const double gflop_count = 2.0 * N * N * N / 1.0e9;

    // allocate memory for matrices use 1-D linear arrays
    double* a = (double*) new double [N * N];
    double* b = (double*) new double [N * N];
    double* c_ijk = (double*) new double [N * N];
    double* c_jki = (double*) new double [N * N];
    double* c_ikj = (double*) new double [N * N];

    // initialize matrices A and B with random values
    srandom((unsigned int) time(NULL));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            a[IDX(i,j,N)] = (double) random() / RAND_MAX;
            b[IDX(i,j,N)] = (double) random() / RAND_MAX;
        }
    }

    // compute ijk product
    t0 = wtime();
    matmat_ijk(c_ijk, a, b, N);
    t1 = wtime();
    ijk_time = t1 - t0;

    // compute jki product
    t0 = wtime();
    matmat_jki(c_jki, a, b, N);
    t1 = wtime();
    jki_time = t1 - t0;

    // compute ikj product
    t0 = wtime();
    matmat_ikj(c_ikj, a, b, N);
    t1 = wtime();
    ikj_time = t1 - t0;

    // output results
    printf("%-15s (%d) ", "1D array (ROW)", N);
    printf("ijk: %6.3f gflops, jki: %6.3f gflops, ikj: %6.3f gflops\n",
    gflop_count / ijk_time, gflop_count / jki_time, gflop_count / ikj_time);

    // verify products
    if (verify(c_ijk, c_jki, N))
    {
        printf("Verification error: c_ijk != c_jki\n");
        exit(1);
    }
    if (verify(c_ijk, c_ikj, N))
    {
        printf("Verification error: c_ijk != c_ikj\n");
        exit(1);
    }
    if (verify(c_jki, c_ikj, N))
    {
        printf("Verification error: c_jki != c_ikj\n");
        exit(1);
    }

    // all done
    delete [] a;
    delete [] b;
    delete [] c_ijk;
    delete [] c_jki;
    delete [] c_ikj;

    return 0;
}
