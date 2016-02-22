/*
 * $Smake: g++ -DN=500 -Wall -O2 -o %F %f -lpthread -lrt
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Benchmark serial and parallel matrix-matrix products.  The parallel
 * version uses pthreads and statically assigns a block of N/num_thread
 * rows to each thread.
 */

#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <pthread.h>

using namespace std;

#define MAX_THREADS 512

#if !defined(N)
# define N 500  /* default matrix dimension */
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
// Structure to hold the start and ending rows in a block of rows

class Rows
{
public:
    int start;
    int end;
};

//----------------------------------------------------------------------------
// Compute the product of a range of rows in A with the matrix B

void* rowsTimesCols( void* arg )
{
    Rows* row = (Rows*) arg;
    for ( int i = row->start; i < row->end; i++ )
    {
        for ( int j = 0; j < N; j++ )
        {
            d[i][j] = 0.0;
            for ( int k = 0; k < N; k++ )
            {
                d[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    pthread_exit( NULL );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    double t1, t2;
    double serial_time;
    double multithreaded_time;
    int numberOfThreads;

    if ( argc != 2 )
    {
        printf( "Usage: %s NUMBER_OF_THREADS\n", argv[0] );
        return 0;
    }
    else
    {
        numberOfThreads = atoi( argv[1] );
    }
    if ( numberOfThreads <= 0 || numberOfThreads > MAX_THREADS )
    {
        fprintf( stderr,
                 "Number of threads must be greater than 0 and less than %d\n",
                 MAX_THREADS );
        exit( EXIT_FAILURE );
    }

    srandom( (unsigned int) time( NULL ) );

    printf( "Matrix-Matrix multiply: Matrices are %d x %d\n", N, N );
    printf( "Static thread assignment of %d threads\n", numberOfThreads );

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
            for ( int k = 0; k < N; k++ )
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    t2 = wtime();
    serial_time = t2 - t1;

    // begin multithreaded product

    pthread_attr_t attr;
    pthread_t threadID[numberOfThreads];
    Rows row[numberOfThreads];

    // set thread attributes

    pthread_attr_init( &attr );
    pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );

    // start threads

    t1 = wtime();
    for ( int n = 0; n < numberOfThreads; n++ )
    {
        row[n].start =     n     * N / numberOfThreads;
        row[n].end   = ( n + 1 ) * N / numberOfThreads;
        pthread_create( &threadID[n], NULL, rowsTimesCols, (void*) &row[n] );
    }

    // wait for threads to finish

    for ( int n = 0; n < numberOfThreads; n++ )
    {
        pthread_join( threadID[n], NULL );
    }
    t2 = wtime();
    multithreaded_time = t2 - t1;

    pthread_attr_destroy( &attr );

    printf( "   Serial Time     Multithreaded Time\n" );
    printf( "    (seconds)           (seconds)       Speedup\n" );
    printf( "------------------  ------------------  -------\n" );
    printf( "%12.6f        %12.6f    %10.3f\n", serial_time,
            multithreaded_time, serial_time / multithreaded_time );

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
