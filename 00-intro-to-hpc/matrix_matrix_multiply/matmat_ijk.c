/*
 * $Smake: gcc -DN=500 -Wall -O3 -o %F %f -lrt
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark various ijk-form matrix-matrix products.
 *
 * "Matrix as array of pointers" version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#if !defined(N)
# define N 500  /* default matrix dimension */
#endif

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
void matmat_ijk( double** c, double** a, double** b, int n )
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
 * Compute matrix-matrix product using ikj loop order
 */
void matmat_ikj( double** c, double** a, double** b, int n )
{
    /*
     * REPLACE LINES BELOW WITH BODY OF FUNCTION
     */
    int i, j;
    for ( i = 0; i < n; i++ )
        for ( j = 0; j < n; j++ )
            c[i][j] = 0.0;
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using jik loop order
 */
void matmat_jik( double** c, double** a, double** b, int n )
{
    /*
     * REPLACE LINES BELOW WITH BODY OF FUNCTION
     */
    int i, j;
    for ( i = 0; i < n; i++ )
        for ( j = 0; j < n; j++ )
            c[i][j] = 0.0;
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using jki loop order
 */
void matmat_jki( double** c, double** a, double** b, int n )
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
 * Compute matrix-matrix product using kij loop order
 */
void matmat_kij( double** c, double** a, double** b, int n )
{
    /*
     * REPLACE LINES BELOW WITH BODY OF FUNCTION
     */
    int i, j;
    for ( i = 0; i < n; i++ )
        for ( j = 0; j < n; j++ )
            c[i][j] = 0.0;
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using kji loop order
 */
void matmat_kji( double** c, double** a, double** b, int n )
{
    /*
     * REPLACE LINES BELOW WITH BODY OF FUNCTION
     */
    int i, j;
    for ( i = 0; i < n; i++ )
        for ( j = 0; j < n; j++ )
            c[i][j] = 0.0;
}

/*----------------------------------------------------------------------------
 * Verify product
 */
double verify( double** c, int n )
{
    int i, j;
    double checksum = 0.0;
    for ( i = 0; i < n; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            checksum += c[i][j];
        }
    }
    return checksum;
}

/*----------------------------------------------------------------------------
 * Main program
 */
int main( int argc, char* argv[] )
{
    double t0, t1;
    double etime;
    const double mflop_count = 2.0 * N * N * N / 1.0e6;
    int i, j;
    double** a;   /* matrix A */
    double** b;   /* matrix B */
    double** c;   /* matrix C = A * B */

    printf( "Matrix-Matrix multiply (array of pointers): Matrices are %dx%d\n",
            N, N );

    /*
     * allocate memory for matrices.
     * in this case, matrices are constructed as an array of row pointers
     * each pointing to a row in a contiguous memory block.
     */
    a = (double**) malloc( N * sizeof( double* ) );
    b = (double**) malloc( N * sizeof( double* ) );
    c = (double**) malloc( N * sizeof( double* ) );
    a[0] = (double*) malloc( N * N * sizeof( double ) );
    b[0] = (double*) malloc( N * N * sizeof( double ) );
    c[0] = (double*) malloc( N * N * sizeof( double ) );
    for ( i = 1; i < N; i++ )
    {
        a[i] = &a[0][i * N];
        b[i] = &b[0][i * N];
        c[i] = &c[0][i * N];
    }
    
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
    matmat_ijk( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "ijk: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify( c, N ) );

    /*
     * ikj product
     */
    t0 = wtime();
    matmat_ikj( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "ikj: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify( c, N ) );

    /*
     * jik product
     */
    t0 = wtime();
    matmat_jik( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "jik: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify( c, N ) );

    /*
     * jki product
     */
    t0 = wtime();
    matmat_jki( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "jki: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify( c, N ) );

    /*
     * kij product
     */
    t0 = wtime();
    matmat_kij( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "kij: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify( c, N ) );

    /*
     * kji product
     */
    t0 = wtime();
    matmat_kji( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "kji: %10.6f sec,%12.2f mflops,   checksum = %18.6f\n",
            etime, mflop_count / etime, verify( c, N ) );

    /*
     * all done
     */
    free( a[0] );
    free( b[0] );
    free( c[0] );
    free( a );
    free( b );
    free( c );

    return 0;
}
