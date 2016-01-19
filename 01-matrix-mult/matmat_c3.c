/*
 * $Smake: gcc -DN=500 -Wall -O3 -o %F %f -lrt
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark ijk, jki, and ikj matrix-matrix products.
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
 * Verify product
 */
double verify( double** c, double** d, int n )
{
    int i, j;
    double checksum = 0.0;
    for ( i = 0; i < n; i++ )
    {
	for ( j = 0; j < n; j++ )
	{
	    checksum += c[i][j];
	    if ( c[i][j] != d[i][j] )
	    {
		printf( "[%d][%d]: c = %f; d = %f\n", i, j, c[i][j], d[i][j] );
	    }
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
    double ijk_time;
    double ikj_time;
    double jki_time;
    double cs1, cs2, cs3;
    const double mflop_count = 2.0 * N * N * N / 1.0e6;
    int i, j;
    double** a;   /* matrix A */
    double** b;   /* matrix B */
    double** c1;   /* matrix C = A * B (computed via ijk loop order) */
    double** c2;   /* matrix C = A * B (computed via ikj loop order) */
    double** c3;   /* matrix D = A * B (computed via jki loop order) */

    printf( "Matrix-Matrix multiply (array of pointers): Matrices are %dx%d\n",
	    N, N );

    /*
     * allocate memory for matrices.
     * in this case, matrices are constructed as an array of row pointers
     * each pointing to a row in a contiguous memory block.
     */
    a = (double**) malloc( N * sizeof( double* ) );
    b = (double**) malloc( N * sizeof( double* ) );
    c1 = (double**) malloc( N * sizeof( double* ) );
    c2 = (double**) malloc( N * sizeof( double* ) );
    c3 = (double**) malloc( N * sizeof( double* ) );
    a[0] = (double*) malloc( N * N * sizeof( double ) );
    b[0] = (double*) malloc( N * N * sizeof( double ) );
    c1[0] = (double*) malloc( N * N * sizeof( double ) );
    c2[0] = (double*) malloc( N * N * sizeof( double ) );
    c3[0] = (double*) malloc( N * N * sizeof( double ) );
    for ( i = 1; i < N; i++ )
    {
	a[i] = &a[0][i * N];
	b[i] = &b[0][i * N];
	c1[i] = &c1[0][i * N];
	c2[i] = &c2[0][i * N];
	c3[i] = &c3[0][i * N];
    }
    
    /*
     * initialize matrices
     */
    srandom( (unsigned int) time( NULL ) );
    for ( i = 0; i < N; i++ )
    {
	for ( j = 0; j < N; j++ )
	{
/*	    a[i][j] = 0.1 * ( ( 2 * (i+1) + 5 * (j+1) ) % 10 );
	    b[i][j] = 0.1 * ( ( 4 * (i+1) + 3 * (j+1) ) % 10 );*/
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
     * ikj product
     */
    t0 = wtime();
    matmat_ikj( c2, a, b, N );
    t1 = wtime();
    ikj_time = t1 - t0;

    /*
     * jki product
     */
    t0 = wtime();
    matmat_jki( c3, a, b, N );
    t1 = wtime();
    jki_time = t1 - t0;

    /*
     * output results
     */
    printf( "        ijk                 jki                ikj\n" );
    printf( "------------------  ------------------  ------------------\n" );
    printf( "%10.6g sec%16.6g sec%16.6g sec\n",
	    ijk_time, jki_time, ikj_time );
    printf( "%10.2f mflops %12.2f mflops %12.2f mflops\n",
	    mflop_count / ijk_time, mflop_count / jki_time,
	    mflop_count / ikj_time );

    /*
     * verify products
     */
    cs1 = verify( c1, c2, N );
    cs2 = verify( c1, c3, N );
    cs3 = verify( c2, c3, N );

    if ( cs1 != cs2 || cs1 != cs3 || cs2 != cs3 )
    {
	printf( "checksum error: %18.9f %18.9f %18.9f\n", cs1, cs2, cs3 );
    }

    /*
     * all done
     */
    free( a[0] );
    free( b[0] );
    free( c1[0] );
    free( c2[0] );
    free( c3[0] );
    free( a );
    free( b );
    free( c1 );
    free( c2 );
    free( c3 );

    return 0;
}
