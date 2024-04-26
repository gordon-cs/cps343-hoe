/*
 * $Smake: nvc++ -fast -o %F %f readMatrixSeq.cc wtime.c -lhdf5
 *
 * Power Method
 */

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <unistd.h>
#include <hdf5.h>
#include "wtime.h"
#include "readMatrixSeq.h"

// Compute index into single linear array for matrix element (i,j)
#define IDX(i,j,stride) ((i)*(stride)+(j)) // row major

//----------------------------------------------------------------------------

// Display usage string
void usage(
    char* program_name // in  - name of program
    )
{
    fprintf(stderr, "Usage: %s [-q] [-v] [-e tol] [-m maxiter] filename\n",
            program_name);
}

//----------------------------------------------------------------------------

// Display contents of matrix (stored as 1-D array)
void dumpMatrix(
    double* a,      // in  - address of matrix data
    int rows,       // in  - number of rows in matrix
    int cols,       // in  - number of cols in matrix
    int stride      // in  - row or col length (C/C++ or FORTRAN ordering)
    )
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf(" %8.2f", a[IDX(i,j,stride)]);
        }
        printf("\n");
    }
    printf("\n");
    fflush(stdout);
}

//----------------------------------------------------------------------------

// Display contents of vector
void dumpVector(
    double* x,          // in  - vector (array) of data
    int n               // in  - vector length
    )
{
    for (int i = 0; i < n; i++) printf(" %9.6f\n", x[i]);
    fflush(stdout);
}

//----------------------------------------------------------------------------

// Copy contents of vector x into vector y
void copyVector(
    double* y,          // out - vector (array) of data
    double* x,          // in  - vector (array) of data
    int n               // in  - vector length
    )
{
    for (int i = 0; i < n; i++) y[i] = x[i];
}

//----------------------------------------------------------------------------

// Scale vector
void scaleVector(
    double sf,          // in     - scale factor
    double* x,          // in/out - vector (array) of data
    int n               // in     - vector length
    )
{
    for (int i = 0; i < n; i++) x[i] *= sf;
}

//----------------------------------------------------------------------------

// Compute inner (dot) product of two vectors
double dotProduct(
    double* x,          // in  - vector (array)
    double* y,          // in  - vector (array)
    int n               // in  - vector length
    )
{
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += x[i] * y[i];
    return sum;
}

//----------------------------------------------------------------------------

// Compute L2 norm of vector
double vectorNorm(
    double* x,          // in  - vector (array)
    int n               // in  - vector length
    )
{
    return sqrt(dotProduct(x, x, n));
}

//----------------------------------------------------------------------------

// Matrix-vector product
void matrixVectorProduct(
    double* y,          // out - m-element result vector
    double* a,          // in  - mxn matrix (stored as 1-D array)
    double* x,          // in  - n-element vector
    int m,              // in  - row dimension of matrix, length of y
    int n               // in  - col dimension of matrix, length of x
    )
{
    for (int i = 0; i < m; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
        {
            sum += a[IDX(i,j,n)] * x[j];
        }
        y[i] = sum;
    }
}

//----------------------------------------------------------------------------
// Main program
//----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    double tol = 1.0e-6;                // allowable diff in e-value estimate
    bool quiet = false;                 // true for minimal output
    int verbosity = 0;                  // verbosity output level
    int max_iter = 1000;                // maximum number of iterations
    int num_iter = 0;                   // number of iterations done
    double* A = NULL;                   // matrix
    double* x = NULL;                   // current estimate of eigenvector
    double* y = NULL;                   // new estimate of eigenvector
    int nrow, ncol, n;                  // number of rows/cols in matrix
    double lambda_prev, lambda;         // old and new estimates of e-value
    double delta;                       // difference in eigenvector estimates
    double vector_mag;                  // vector magnitude (L2-norm)
    double t1, t2;                      // timing variables
    double read_time, compute_time;     // elapsed read and compute times

    // process command line options
    int opt;
    while ((opt = getopt(argc, argv, "e:m:qv")) != -1)
    {
        switch (opt)
        {
            case 'e':
                tol = atof(optarg);
                break;
            case 'm':
                max_iter = atoi(optarg);
                break;
            case 'q':
                quiet = true;
                break;
            case 'v':
                verbosity++;
                break;
            default:
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // validate parameters that may have been changed on command line
    if (tol <= 0.0)
    {
        fprintf(stderr, "tolerance must be positive\n");
        exit(EXIT_FAILURE);
    }
    if (max_iter <= 0)
    {
        fprintf(stderr, "maximum iterations must be positive\n");
        exit(EXIT_FAILURE);
    }

    // remove command line options from argument list
    argv[optind - 1] = argv[0];
    argc -= (optind - 1);
    argv += (optind - 1);

    // make sure we've got a filename argument
    if (argc != 2)
    {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    // read matrix data from HDF5 file
    t1 = wtime();
    readMatrix(argv[1], "/A/value", &A, &nrow, &ncol);
    t2 = wtime();
    read_time = t2 - t1;
    n = nrow;

    // make sure matrix is square
    if (nrow != ncol)
    {
        fprintf(stderr, "matrix must be square: rows = %d, cols = %d\n",
                nrow, ncol);
        return EXIT_FAILURE;
    }
    if (!quiet && verbosity > 1) dumpMatrix(A, n, n, n);

    // allocate memory for current and previous eigenvector vector estimates
    x = new double [n];
    y = new double [n];

    // initialize estimate of eigenvector
    for (int i = 0; i < n; i++) x[i] = 1.0 / sqrt(n);

    // main power method loop
    lambda = 0.0;
    lambda_prev = lambda + 2.0 * tol;
    num_iter = 0;
    delta = fabs(lambda - lambda_prev);
    t1 = wtime();
    while (delta > tol && num_iter++ < max_iter)
    {
        // compute new eigenvector estimate y := A*x
        matrixVectorProduct(y, A, x, n, n);

        // update estimates of eigenvalue and unnormalized eigenvector
        lambda_prev = lambda;
        lambda = dotProduct(x, y, n);
        delta = fabs(lambda - lambda_prev);

        // update and normalize eigenvector
        copyVector(x, y, n);
        vector_mag = vectorNorm(x, n);
        scaleVector(1.0 / vector_mag, x, n);

        // display status if requested
        if (!quiet)
        {
            if (verbosity > 1) dumpVector(x, n);
            if (verbosity > 0)
            {
                printf("k = %d: lambda = %f, diff = %e\n",
                       num_iter, lambda, delta);
            }
        }
    }
    t2 = wtime();
    compute_time = t2 - t1;

    // report
    if (num_iter > max_iter)
    {
        fprintf(stderr, "*** WARNING ****: ");
        fprintf(stderr, "maximum number of iterations exceeded\n");
    }

    if (quiet)
    {
        printf("%5d x %5d %5d %10.6f %10.6f %10.6f\n",
               n, n, num_iter, lambda, read_time, compute_time);
    }
    else
    {
        printf("matrix dimensions: %d x %d; tolerance: %8.4e; ", n, n, tol );
        printf("max iterations: %d\n", max_iter);
        printf("elapsed HDF5 read time = %10.6f seconds\n", read_time);
        printf("elapsed compute time   = %10.6f seconds\n", compute_time);
        printf("eigenvalue = %f found in %d iterations\n", lambda, num_iter);
    }

    // all done, cleanup and quit
    destroyMatrix(A);
    delete [] x;
    delete [] y;

    return 0;
}
