/*
 * $Smake: pgc++ -fast -o %F %f wtime.c readMatrixSeq.cc -L$LIBRARY_PATH -lhdf5
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

// Compute index into single linear array for matrix element (i,j)
#define IDX(i,j,stride) ((i)*(stride)+(j)) // row major (c/c++ ordering)

//----------------------------------------------------------------------------

// Prototype of function declared in external file
void readMatrix( const char* fname, const char* path, double** a,
                 int* rows, int* cols );

//----------------------------------------------------------------------------

// Display usage string
void usage(
    char* program_name // in  - name of program
    )
{
    fprintf( stderr, "Usage: %s [-v] [-e tol] [-m maxiter] filename\n",
             program_name );
}

//----------------------------------------------------------------------------

// Display contents of matrix (stored as 1-D array)
void dumpMatrix(
    double* a,      // in  - address of matrix data
    int rows,       // in  - number of rows in matrix
    int cols,       // in  - number of cols in matrix
    int stride      // in  - row length in memory (assuming C/C++ storage)
    )
{
    for ( int i = 0; i < rows; i++ )
    {
        for ( int j = 0; j < cols; j++ )
        {
            printf( " %8.2f", a[IDX(i,j,stride)] );
        }
        printf( "\n" );
    }
    printf( "\n" );
    fflush( stdout );
}

//----------------------------------------------------------------------------

// Display contents of vector
void dumpVector(
    double* x,       // in  - vector (array) of data
    int n            // in  - vector length
    )
{
    for ( int i = 0; i < n; i++ )
    {
        printf( " %8.2f\n", x[i] );
    }
    fflush( stdout );
}

//----------------------------------------------------------------------------
// Main program
//----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    double tol = 1.0e-6;                     // allowable diff in e-value est
    int maxIter = 500;                       // max number of iterations
    int verbosity = 0;                       // verbosity output level
    double* a;                               // matrix
    double* x;                               // current est. of eigenvector
    double* y;                               // new est. of eigenvector
    int n;                                   // number of rows/cols in matrix
    double lambda_old, lambda;               // old and new est. of e-value
    double t1, t2, read_time, compute_time;  // timing variables
    int numIter;                             // iteration counter

    // Process command line
    int opt;
    while ( ( opt = getopt( argc, argv, "ve:m:" ) ) != -1 )
    {
        switch ( opt )
        {
            case 'e':
                tol = atof( optarg );
                break;
            case 'm':
                maxIter = atoi( optarg );
                break;
            case 'v':
                verbosity++;
                break;
            default:
                usage( argv[0] );
                return EXIT_FAILURE;
        }
    }
    argv[optind - 1] = argv[0];
    argc -= ( optind - 1 );
    argv += ( optind - 1 );

    // Validate parameters that may have been changed on command line
    if ( tol <= 0.0 )
    {
        fprintf( stderr, "tolerance must be positive\n" );
        exit( EXIT_FAILURE );
    }
    if ( maxIter <= 0 )
    {
        fprintf( stderr, "maximum iterations must be positive\n" );
        exit( EXIT_FAILURE );
    }

    // Make sure we've got a filename argument
    if ( argc != 2 )
    {
        usage( argv[0] );
        return EXIT_FAILURE;
    }

    if ( verbosity > 0 )
    {
        printf( "tolerance = %e\n", tol );
        printf( "maximum number of iterations = %d\n", maxIter );
    }

    // Read matrix data
    {
        int nrow, ncol;
        t1 = wtime();
        readMatrix( argv[1], "/A/value", &a, &nrow, &ncol );
        t2 = wtime();
        read_time = t2 - t1;

        if ( nrow != ncol )
        {
            fprintf( stderr, "matrix must be square: rows = %d, cols = %d\n",
                     nrow, ncol );
            exit( EXIT_FAILURE );
        }
        n = nrow;
        if ( verbosity > 2 )
        {
            printf( "Matrix A:\n" );
            dumpMatrix( a, nrow, ncol, ncol );
        }
    }

    // Allocate memory for previous and current eigenvector vector estimates
    x = new double [n];
    y = new double [n];

    // Initialize estimate of eigenvector
    for ( int i = 0; i < n; i++ ) x[i] = 1.0;

    // Main power method loop
    lambda = 0.0;
    lambda_old = lambda + 2.0 * tol;
    numIter = 0;
    t1 = wtime();
    while ( fabs( lambda - lambda_old ) > tol && numIter <= maxIter )
    {
        // normalize x
        double xscale = 0.0;
        for ( int i = 0; i < n; i++ ) xscale += x[i] * x[i];
        xscale = sqrt( xscale );
        for ( int i = 0; i < n; i++ ) x[i] /= xscale;

        // compute y := A*x
        for ( int i = 0; i < n; i++ )
        {
            double sum = 0.0;
            for ( int j = 0; j < n; j++ )
            {
                sum += x[j] * a[IDX(i,j,n)];
            }
            y[i] = sum;
        }

        // compute estimate of eigenvalue and copy vector y in to x
        lambda_old = lambda;
        lambda = 0.0;
        for ( int i = 0; i < n; i++ )
        {
            lambda += y[i] * x[i];
            x[i] = y[i];
        }

        numIter++;

        if ( verbosity > 2 ) dumpVector( x, n );

        if ( verbosity > 1 )
        {
            printf( "k = %d: lambda = %f, diff = %e\n",
                    numIter, lambda, fabs( lambda - lambda_old ) );
        }
    }
    t2 = wtime();
    compute_time = t2 - t1;

    // Report
    if ( numIter > maxIter )
    {
        fprintf( stderr, "*** WARNING ****: " );
        fprintf( stderr, "maximum number of iterations exceeded\n" );
    }

    printf( "magnitude of eigenvalue = %f found in %d iterations\n",
            lambda, numIter );
    printf( "elapsed read time    = %10.6f seconds\n", read_time );
    printf( "elapsed compute time = %10.6f seconds\n", compute_time );

    // All done, cleanup and quit
    delete [] x;
    delete [] y;
    delete [] a;

    return 0;
}
