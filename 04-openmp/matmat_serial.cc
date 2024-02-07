/*
 * $Smake: g++ -Wall -O2 -o %F %f
 *
 * Benchmark serial and parallel matrix-matrix products.  Parallel version
 * uses OpenMP and parallelizes the outer loop so that rows are assigned to
 * threads.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <time.h>

#if !defined(N)
# define N 1000 // default matrix dimension
#endif
#define IDX(i,j,stride) ((i)*(stride)+(j)) // row major (C/C++)

//----------------------------------------------------------------------------
// Returns the number of seconds since some fixed arbitrary time in the past.

double wtime(void)
{
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return double(ts.tv_sec + ts.tv_nsec / 1.0e9);
}

//----------------------------------------------------------------------------
// Compute matrix-matrix product -- this is baseline serial version

void doSerialProduct(double* c, double* a, double* b, int n)
{
    for (int i = 0; i < n * n; i++)
    {
        c[i] = 0.0;
    }
    for (int i = 0; i < n; i++)
    {
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
// Compute matrix-matrix product -- parallel version
// **** THIS FUNCTION SHOULD BE MODIFIED TO RUN IN PARALLEL ****

void doParallelProduct(double* c, double* a, double* b, int n)
{
    for (int i = 0; i < n * n; i++)
    {
        c[i] = 0.0;
    }

    for (int i = 0; i < n; i++)
    {
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
// verify products; no output means the products match

void verifyProduct(double* a, double* b, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (a[IDX(i,j,n)] != b[IDX(i,j,n)])
            {
                printf("** ERROR ** location [%d][%d]: %f != %f\n",
                        i, j, a[IDX(i,j,n)], b[IDX(i,j,n)]);
            }
        }
    }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    double t1, t2;
    double serial_time;
    double parallel_time;

    // set matrix dimension from command line or use default value

    int n = (argc > 1 ? atoi(argv[1]) : N);
    if (n <= 0)
    {
        fprintf(stderr, "N must be positive, got N = %d\n", n);
        return EXIT_FAILURE;
    }

    printf("Matrix-Matrix multiply: Matrices are %d x %d\n", n, n);

    // allocate and initialize matrices

    double* a = new double [n * n];  // matrix A
    double* b = new double [n * n];  // matrix B
    double* c = new double [n * n];  // matrix C = A * B (serial product)
    double* d = new double [n * n];  // matrix D = A * B (parallel product)

    srandom((unsigned int) time(NULL));
    for (int i = 0; i < n * n; i++)
    {
        a[i] = double(random()) / RAND_MAX;
        b[i] = double(random()) / RAND_MAX;
    }

    // compute serial product

    t1 = wtime();
    doSerialProduct(c, a, b, n);
    t2 = wtime();
    serial_time = t2 - t1;

    // compute product in parallel

    t1 = wtime();
    doParallelProduct(d, a, b, n);
    t2 = wtime();
    parallel_time = t2 - t1;

    // verify product

    verifyProduct(c, d, n);

    // report

    printf("   Serial Time        Parallel Time\n");
    printf("    (seconds)           (seconds)       Speedup\n");
    printf("------------------  ------------------  -------\n");
    printf("%12.6f        %12.6f    %10.3f\n", serial_time,
            parallel_time, serial_time / parallel_time);

    // all done

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
    return 0;
}
