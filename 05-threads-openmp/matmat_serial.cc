/*
 * $Smake: g++ -DN=1000 -Wall -O2 -o %F %f -lrt
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark serial and parallel matrix-matrix products (after modification).
 * The parallel version should use OpenMP.
 */

#include <cstdio>
#include <cstdlib>
#include <time.h>

using namespace std;

#if !defined(N)
# define N 1000  /* default matrix dimension */
#endif

// Using static memory allocation.
// Arrays are global so they will not be allocated on stack AND so that
// references to them do not have to be passed to thread functions.

double a[N][N];  /* matrix A */
double b[N][N];  /* matrix B */
double c[N][N];  /* matrix C = A * B (serial product) */
double d[N][N];  /* matrix D = A * B (parallel product */

//----------------------------------------------------------------------------
// Returns the number of seconds since some fixed arbitrary time in the past.

double wtime( void )
{
    timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return double( ts.tv_sec + ts.tv_nsec / 1.0e9 );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    double t1, t2;
    double serial_time;
    double parallel_time;

    srandom( (unsigned int) time( NULL ) );

    printf( "Matrix-Matrix multiply: Matrices are %d x %d\n", N, N );

    // initialize matrices

    for ( int i = 0; i < N; i++ )
    {
        for ( int j = 0; j < N; j++ )
        {
            a[i][j] = double( random() ) / RAND_MAX;
            b[i][j] = double( random() ) / RAND_MAX;
        }
    }

    // begin serial product

    t1 = wtime();
    for ( int i = 0; i < N; i++ )
    {
        for ( int j = 0; j < N; j++ )
        {
            c[i][j] = 0.0;
        }
        for ( int k = 0; k < N; k++ )
        {
            for ( int j = 0; j < N; j++ )
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    t2 = wtime();
    serial_time = t2 - t1;

    //********************************************************************
    //* Modify nested loops below so the computation is done in parallel *
    //********************************************************************

    t1 = wtime();
    for ( int i = 0; i < N; i++ )
    {
        for ( int j = 0; j < N; j++ )
        {
            d[i][j] = 0.0;
        }
        for ( int k = 0; k < N; k++ )
        {
            for ( int j = 0; j < N; j++ )
            {
                d[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    t2 = wtime();
    parallel_time = t2 - t1;

    printf( "   Serial Time        Parallel Time\n" );
    printf( "    (seconds)           (seconds)       Speedup\n" );
    printf( "------------------  ------------------  -------\n" );
    printf( "%12.6f        %12.6f    %10.3f\n", serial_time,
            parallel_time, serial_time / parallel_time );

    // verify products; no output means the products match

    for ( int i = 0; i < N; i++ )
    {
        for ( int j = 0; j < N; j++ )
        {
            if ( c[i][j] != d[i][j] )
            {
                printf( "[%d][%d]: c = %f; d = %f\n", i, j, c[i][j], d[i][j] );
            }
        }
    }

    // all done

    return 0;
}
