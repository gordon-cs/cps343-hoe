/*
 * $Smake: gcc -DN=500 -Wall -O3 -o %F %f -lrt
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark ijk, jki, and ikj matrix-matrix products.
 *
 * "Matrix with 1D arrays" version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#if !defined(N)
# define N 500  /* default matrix dimension */
#endif

#define IDX(i,j,rowlen) ((i)*(rowlen)+j)

/*----------------------------------------------------------------------------
 * Returns the number of seconds since some fixed arbitrary time in the past
 */

double wtime( void )
{
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return (double) ( ts.tv_sec + ts.tv_nsec / 1.0e9 );
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using ijk loop order
 */
void matmat_ijk( double* c, double* a, double* b, int n )
{
    int i, j, k;
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            c[IDX(i,j,n)] = 0.0;
            for ( k = 0; k < n; k++ )
            {
                c[IDX(i,j,n)] += a[IDX(i,k,n)] * b[IDX(k,j,n)];
            }
        }
    }
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using jki loop order
 */
void matmat_jki( double* c, double* a, double* b, int n )
{
    int i, j, k;
    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < n; i++ )
        {
            c[IDX(i,j,n)] = 0.0;
        }
        for ( k = 0; k < n; k++ )
        {
            for ( i = 0; i < n; i++ )
            {
                c[IDX(i,j,n)] += a[IDX(i,k,n)] * b[IDX(k,j,n)];
            }
        }
    }
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using ikj loop order
 */
void matmat_ikj( double* c, double* a, double* b, int n )
{
    int i, j, k;
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            c[IDX(i,j,n)] = 0.0;
        }
        for ( k = 0; k < n; k++ )
        {
            for ( j = 0; j < n; j++ )
            {
                c[IDX(i,j,n)] += a[IDX(i,k,n)] * b[IDX(k,j,n)];
            }
        }
    }
}

/*----------------------------------------------------------------------------
 * Verify product
 */
int verify( double* d, double* c, int n )
{
    int i, j;
    int status = 0;
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            if ( c[IDX(i,j,n)] != d[IDX(i,j,n)] )
            {
                status++;
                printf( "[%d][%d]: c = %f; d = %f\n", i, j,
                        c[IDX(i,j,n)], d[IDX(i,j,n)] );
            }
        }
    }
    return status;
}

/*----------------------------------------------------------------------------
 * Main program
 */
int main( int argc, char* argv[] )
{
    double t0, t1;
    double ijk_time, jki_time, ikj_time;
    int i, j;

    const double mflop_count = 2.0 * N * N * N / 1.0e6;

    double* a;   /* matrix A */
    double* b;   /* matrix B */
    double* c1;   /* matrix C = A * B (computed via ijk loop order) */
    double* c2;   /* matrix C = A * B (computed via ikj loop order) */
    double* c3;   /* matrix D = A * B (computed via jki loop order) */

    printf( "Matrix-Matrix multiply (1D arrays): Matrices are %dx%d\n", N, N );

    /*
     * allocate memory for matrices.
     * in this case, matrices are linear arrays (row by row).
     */
    a = (double*) malloc( N * N * sizeof( double ) );
    b = (double*) malloc( N * N * sizeof( double ) );
    c1 = (double*) malloc( N * N * sizeof( double ) );
    c2 = (double*) malloc( N * N * sizeof( double ) );
    c3 = (double*) malloc( N * N * sizeof( double ) );
    
    /*
     * initialize matrices
     */
    srandom( (unsigned int) time( NULL ) );
    for ( i = 0; i < N; i++ )
    {
        for ( j = 0; j < N; j++ )
        {
            a[IDX(i,j,N)] = (double) random() / RAND_MAX;
            b[IDX(i,j,N)] = (double) random() / RAND_MAX;
        }
    }

    /*
     * ijk product
     */
    t0 = wtime();
    matmat_ijk( c1, a, b, N );
    t1 = wtime();
    ijk_time = t1 - t0;

    /*
     * jki product
     */
    t0 = wtime();
    matmat_jki( c2, a, b, N );
    t1 = wtime();
    jki_time = t1 - t0;

    /*
     * ikj product
     */
    t0 = wtime();
    matmat_ikj( c3, a, b, N );
    t1 = wtime();
    ikj_time = t1 - t0;

    /*
     * output results
     */
    printf( "        ijk                 jki                ikj\n" );
    printf( "------------------  ------------------  ------------------\n" );
    printf( "%10.6g sec%16.6g sec%16.6g sec\n", ijk_time, jki_time, ikj_time );
    printf( "%10.2f mflops %12.2f mflops %12.2f mflops\n",
            mflop_count / ijk_time, mflop_count / jki_time,
            mflop_count / ikj_time );

    /*
     * verify products
     */
    if ( verify( c1, c2, N ) )
    {
        printf( "Verification error: c1 != c2\n" );
        exit( 1 );
    }
    if ( verify( c1, c3, N ) )
    {
        printf( "Verification error: c1 != c3\n" );
        exit( 1 );
    }
    if ( verify( c2, c3, N ) )
    {
        printf( "Verification error: c2 != c3\n" );
        exit( 1 );
    }

    /*
     * all done
     */
    free( a );
    free( b );
    free( c1 );
    free( c2 );
    free( c3 );

    return 0;
}
