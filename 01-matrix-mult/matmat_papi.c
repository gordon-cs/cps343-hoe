/*
 * $Smake: gcc -DN=500 -Wall -O3 -o %F %f -lrt -lpapi
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark various ijk-form matrix-matrix products.
 *
 * PAPI version: be sure to load PAPI module before compiling and running:
 *     $ module load papi
 *
 * "Matrix as array of pointers" version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <papi.h>

#if !defined(N)
# define N 500  /* default matrix dimension */
#endif

int event_set = PAPI_NULL;
long long L1_misses[6];
long long L2_misses[6];

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
 *
 */
void handle_papi_error( int ret )
{
    fprintf( stderr, "PAPI error %d: %s\n", ret, PAPI_strerror( ret ) );
    exit( EXIT_FAILURE );
}
    
/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using ijk loop order
 */
void matmat_ijk( double** c, double** a, double** b, int n )
{
    int i, j, k;
    long long values[3];
    int ret;
    ret = PAPI_reset( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    ret = PAPI_start( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
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
    ret = PAPI_stop( event_set, values );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    L1_misses[0] = values[0];
    L2_misses[0] = values[1];
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using ikj loop order
 */
void matmat_ikj( double** c, double** a, double** b, int n )
{
    int i, j, k;
    long long values[3];
    int ret;
    ret = PAPI_reset( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    ret = PAPI_start( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
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
    ret = PAPI_stop( event_set, values );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    L1_misses[1] = values[0];
    L2_misses[1] = values[1];
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using jik loop order
 */
void matmat_jik( double** c, double** a, double** b, int n )
{
    int i, j, k;
    long long values[3];
    int ret;
    ret = PAPI_reset( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    ret = PAPI_start( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    for ( j = 0; j < n; j++ )
    {
	for ( i = 0; i < n; i++ )
	{
	    c[i][j] = 0.0;
	    for ( k = 0; k < n; k++ )
	    {
		c[i][j] += a[i][k] * b[k][j];
	    }
	}
    }
    ret = PAPI_stop( event_set, values );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    L1_misses[2] = values[0];
    L2_misses[2] = values[1];
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using jki loop order
 */
void matmat_jki( double** c, double** a, double** b, int n )
{
    int i, j, k;
    long long values[3];
    int ret;
    ret = PAPI_reset( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    ret = PAPI_start( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
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
    ret = PAPI_stop( event_set, values );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    L1_misses[3] = values[0];
    L2_misses[3] = values[1];
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using kij loop order
 */
void matmat_kij( double** c, double** a, double** b, int n )
{
    int i, j, k;
    long long values[3];
    int ret;
    ret = PAPI_reset( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    ret = PAPI_start( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    for ( i = 0; i < n; i++ )
    {
	for ( j = 0; j < n; j++ )
	{
	    c[i][j] = 0.0;
	}
    }
    for ( k = 0; k < n; k++ )
    {
	for ( i = 0; i < n; i++ )
	{
	    for ( j = 0; j < n; j++ )
	    {
		c[i][j] += a[i][k] * b[k][j];
	    }
	}
    }
    ret = PAPI_stop( event_set, values );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    L1_misses[4] = values[0];
    L2_misses[4] = values[1];
}

/*----------------------------------------------------------------------------
 * Compute matrix-matrix product using kji loop order
 */
void matmat_kji( double** c, double** a, double** b, int n )
{
    int i, j, k;
    long long values[3];
    int ret;
    ret = PAPI_reset( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    ret = PAPI_start( event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    for ( i = 0; i < n; i++ )
    {
	for ( j = 0; j < n; j++ )
	{
	    c[i][j] = 0.0;
	}
    }
    for ( k = 0; k < n; k++ )
    {
	for ( j = 0; j < n; j++ )
	{
	    for ( i = 0; i < n; i++ )
	    {
		c[i][j] += a[i][k] * b[k][j];
	    }
	}
    }
    ret = PAPI_stop( event_set, values );
    if ( ret != PAPI_OK ) handle_papi_error( ret );
    L1_misses[5] = values[0];
    L2_misses[5] = values[1];
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
    int event_code;
    int ret;

    printf( "Matrix-Matrix multiply (array of pointers): Matrices are %dx%d\n",
	    N, N );

    if ( PAPI_library_init( PAPI_VER_CURRENT ) != PAPI_VER_CURRENT )
    {
	fprintf( stderr, "PAPI Library mismatch; rebuild executable\n" );
	exit( EXIT_FAILURE );
    }

    ret = PAPI_create_eventset( &event_set );
    if ( ret != PAPI_OK ) handle_papi_error( ret );

    event_code = PAPI_L1_DCM;
    ret = PAPI_add_event( event_set, event_code );
    if ( ret != PAPI_OK ) handle_papi_error( ret );

    event_code = PAPI_L2_DCM;
    ret = PAPI_add_event( event_set, event_code );
    if ( ret != PAPI_OK ) handle_papi_error( ret );

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
    for ( i = 0; i < N; i++ )
    {
	for ( j = 0; j < N; j++ )
	{
	    a[i][j] = 0.1 * ( ( 2 * (i+1) + 5 * (j+1) ) % 10 );
	    b[i][j] = 0.1 * ( ( 4 * (i+1) + 3 * (j+1) ) % 10 );
	}
    }
    printf(
	"         seconds       MFLOPS    L1 misses    L2 misses     checksum\n"
	);
    printf(
	"--- ------------ ------------ ------------ ------------ ------------\n"
	);

    /*
     * ijk product
     */
    t0 = wtime();
    matmat_ijk( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "%s %12.6f %12.2f %12lld %12lld %12.2f\n", "ijk", etime,
	    mflop_count / etime, L1_misses[0], L2_misses[0], verify( c, N ) );

    /*
     * ikj product
     */
    t0 = wtime();
    matmat_ikj( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "%s %12.6f %12.2f %12lld %12lld %12.2f\n", "ikj", etime,
	    mflop_count / etime, L1_misses[1], L2_misses[1], verify( c, N ) );

    /*
     * jik product
     */
    t0 = wtime();
    matmat_jik( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "%s %12.6f %12.2f %12lld %12lld %12.2f\n", "jik", etime,
	    mflop_count / etime, L1_misses[2], L2_misses[2], verify( c, N ) );

    /*
     * jki product
     */
    t0 = wtime();
    matmat_jki( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "%s %12.6f %12.2f %12lld %12lld %12.2f\n", "jki", etime,
	    mflop_count / etime, L1_misses[3], L2_misses[3], verify( c, N ) );

    /*
     * kij product
     */
    t0 = wtime();
    matmat_kij( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "%s %12.6f %12.2f %12lld %12lld %12.2f\n", "kij", etime,
	    mflop_count / etime, L1_misses[4], L2_misses[4], verify( c, N ) );

    /*
     * kji product
     */
    t0 = wtime();
    matmat_kji( c, a, b, N );
    t1 = wtime();
    etime = t1 - t0;
    printf( "%s %12.6f %12.2f %12lld %12lld %12.2f\n", "kji", etime,
	    mflop_count / etime, L1_misses[5], L2_misses[5], verify( c, N ) );

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
