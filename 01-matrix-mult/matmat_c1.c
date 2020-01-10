/*
 * $Smake: gcc -DN=500 -Wall -O3 -o %F %f
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark ijk, jki, and ikj matrix-matrix products.
 *
 * "Static Matrix" version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#if !defined(N)
# define N 500  /* default matrix dimension */
#endif

const double EPSILON = 1.0E-14;

/*
 * using static memory allocation --
 * arrays are global so they will not be allocated on stack
 */
double a[N][N];    /* matrix A */
double b[N][N];    /* matrix B */
double c1[N][N];   /* matrix C = A * B (computed via ijk loop order) */
double c2[N][N];   /* matrix C = A * B (computed via ikj loop order) */
double c3[N][N];   /* matrix D = A * B (computed via jki loop order) */

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
void matmat_ijk( double c[N][N], double a[N][N], double b[N][N], int n )
{
    int i, j, k;
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            c[i][j] = 0.0;
            for ( k = 0; k < n; k++ )
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using jki loop order
 */
void matmat_jki( double c[N][N], double a[N][N], double b[N][N], int n )
{
    int i, j, k;
    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < n; i++ )
        {
            c[i][j] = 0.0;
        }
        for ( k = 0; k < n; k++ )
        {
            for ( i = 0; i < n; i++ )
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using ikj loop order
 */
void matmat_ikj( double c[N][N], double a[N][N], double b[N][N], int n )
{
    int i, j, k;
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            c[i][j] = 0.0;
        }
        for ( k = 0; k < n; k++ )
        {
            for ( j = 0; j < n; j++ )
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

/*----------------------------------------------------------------------------
 * Verify product
 */
int verify( double c[N][N], double d[N][N], int n )
{
    int i, j;
    int status = 0;
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            if ( fabs( c[i][j] - d[i][j] ) > EPSILON )
            {
                status++;
                printf( "[%d][%d]: c = %f; d = %f\n", i, j, c[i][j], d[i][j] );
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

    const double gflop_count = 2.0 * N * N * N / 1.0e9;


    /*
     * initialize matrices
     */
    srandom( (unsigned int) time( NULL ) );
    for ( i = 0; i < N; i++ )
    {
        for ( j = 0; j < N; j++ )
        {
            a[i][j] = (double) random() / RAND_MAX;
            b[i][j] = (double) random() / RAND_MAX;
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
    printf( "%-15s (%d) ", "2D static array", N );
    printf( "ijk: %6.3f gflops, jki: %6.3f gflops, ikj: %6.3f gflops\n",
    gflop_count / ijk_time, gflop_count / jki_time, gflop_count / ikj_time );

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
    return 0;
}
