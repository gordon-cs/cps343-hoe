/*
 * $Smake: g++ -DN=500 -Wall -O3 -funroll-loops -o %F %f
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark various ijk-form matrix-matrix products.
 *
 * "Matrix as array of pointers" version.
 */

#include <cstdio>
#include <cmath>
#include <time.h>

#if !defined(N)
# define N 500  // default matrix dimension
#endif

#define IDX(i,j,stride) ((i)*(stride)+(j)) // row-major
//#define IDX(i,j,stride) ((i)+(stride)*(j)) // column-major

//----------------------------------------------------------------------------
// Returns the number of seconds since some fixed arbitrary time in the past

double wtime(void)
{
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) (ts.tv_sec + ts.tv_nsec / 1.0e9);
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using ijk loop order

void matmat_ijk(double* c, double* a, double* b, int n)
{
    // initialize result matrix to zero
    for (int i = 0; i < n * n; i++)
        c[i] = 0.0;

    // compute product
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                c[IDX(i,j,n)] += a[IDX(i,k,n)] * b[IDX(k,j,n)];
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using ikj loop order

void matmat_ikj(double* c, double* a, double* b, int n)
{
    // initialize result matrix to zero
    for (int i = 0; i < n * n; i++)
        c[i] = 0.0;

    // REPLACE THIS COMMENT WITH CODE TO COMPUTE PRODUCT
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using jik loop order

void matmat_jik(double* c, double* a, double* b, int n)
{
    // initialize result matrix to zero
    for (int i = 0; i < n * n; i++)
        c[i] = 0.0;

    // REPLACE THIS COMMENT WITH CODE TO COMPUTE PRODUCT
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using jki loop order

void matmat_jki(double* c, double* a, double* b, int n)
{
    // initialize result matrix to zero
    for (int i = 0; i < n * n; i++)
        c[i] = 0.0;
 
    // REPLACE THIS COMMENT WITH CODE TO COMPUTE PRODUCT
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using kij loop order

void matmat_kij(double* c, double* a, double* b, int n)
{
    // initialize result matrix to zero
    for (int i = 0; i < n * n; i++)
        c[i] = 0.0;

    // REPLACE THIS COMMENT WITH CODE TO COMPUTE PRODUCT
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product using kji loop order

void matmat_kji(double* c, double* a, double* b, int n)
{
    // initialize result matrix to zero
    for (int i = 0; i < n * n; i++)
        c[i] = 0.0;

    // REPLACE THIS COMMENT WITH CODE TO COMPUTE PRODUCT
}

//----------------------------------------------------------------------------
// Verify product

double verify(double* c, int n)
{
    double checksum = 0.0;
    for (int i = 0; i < n * n; i++)
        checksum += c[i];
    return checksum;
}

//----------------------------------------------------------------------------
// Main program

int main(int argc, char* argv[])
{
    double t0, t1;
    double etime;
    const double mflop_count = 2.0 * N * N * N / 1.0e6;

    printf("Matrix-Matrix multiply (1D Array w/macro): Matrices are %dx%d\n",
            N, N);

    // allocate memory for matrices.
    double* a = new double [N * N];
    double* b = new double [N * N];
    double* c = new double [N * N];
    
    // initialize matrices
    srandom((unsigned int) time(NULL));
    for (int i = 0; i < N * N; i++)
    {
        a[i] = (double) random() / RAND_MAX;
        b[i] = (double) random() / RAND_MAX;
    }

    // compute ijk product
    t0 = wtime();
    matmat_ijk(c, a, b, N);
    t1 = wtime();
    etime = t1 - t0;
    printf("ijk: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify(c, N));

    // compute ikj product
    t0 = wtime();
    matmat_ikj(c, a, b, N);
    t1 = wtime();
    etime = t1 - t0;
    printf("ikj: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify(c, N));

    // compute jik product
    t0 = wtime();
    matmat_jik(c, a, b, N);
    t1 = wtime();
    etime = t1 - t0;
    printf("jik: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify(c, N));

    // compute jki product
    t0 = wtime();
    matmat_jki(c, a, b, N);
    t1 = wtime();
    etime = t1 - t0;
    printf("jki: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify(c, N));

    // compute kij product
    t0 = wtime();
    matmat_kij(c, a, b, N);
    t1 = wtime();
    etime = t1 - t0;
    printf("kij: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify(c, N));

    // compute kji product
    t0 = wtime();
    matmat_kji(c, a, b, N);
    t1 = wtime();
    etime = t1 - t0;
    printf("kji: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify(c, N));

    // all done
    delete [] a;
    delete [] b;
    delete [] c;

    return 0;
}
