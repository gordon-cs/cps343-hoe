/*
 * $Smake: gcc -O2 -o %F_2D %f -lrt; g++ -DUSE_MACRO -O2 -o %F_1D %f -lrt
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * Written to test timing differences between using pointers for 2D array
 * vs. using linear access into a 1D array.
 *
 * Two-dimensional arrays declared in C or C++ cannot be passed to functions
 * if the second dimension is not known and specified at compile time.  There
 * are two ways to work around this:
 * 1. Use a single, linear array to store the data and compute the offsets
 *    in the array manually, or
 * 2. Allocate an array of pointers to pointers to rows in the array, allocate
 *    contiguous linear data for the array, and set the pointers to the start
 *    of each row.
 * The second approach has the advantage that that the code can address array
 * elements using standard subscript notation but it has memory overhead for
 * the row pointers.  The first approach uses a minimum of memory, but the
 * allocated row length must be explicitly known (or computed) any time
 * an array location is computed.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#if !defined(N)
# define N 1000
#endif

#if defined(USE_MACRO)
# define IDX(row,col,stride) ((row)*(stride)+(col))
#endif

/*
 *----------------------------------------------------------------------------
 * Returns the number of seconds since some fixed arbitrary time in the past
 */

double wtime( void )
{
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return (double) ( ts.tv_sec + ts.tv_nsec / 1.0e9 );
}

/*
 *----------------------------------------------------------------------------
 *----------------------------------------------------------------------------
 * Main program
 */

int main( int argc, char *argv[] )
{
    const int n = N;

    /*
     * allocate contiguous memory for matrices
     */
#if defined(USE_MACRO)
    double* a = (double*) malloc( n * n * sizeof( double ) );
    double* b = (double*) malloc( n * n * sizeof( double ) );
    double* c = (double*) malloc( n * n * sizeof( double ) );
#else
    double** a = (double**) malloc( n * sizeof( double* ) );
    double** b = (double**) malloc( n * sizeof( double* ) );
    double** c = (double**) malloc( n * sizeof( double* ) );
    a[0] = (double*) malloc( n * n * sizeof( double ) );
    b[0] = (double*) malloc( n * n * sizeof( double ) );
    c[0] = (double*) malloc( n * n * sizeof( double ) );
    for ( int i = 1; i < n; i++ )
    {
        a[i] = &a[0][i * n];
        b[i] = &b[0][i * n];
        c[i] = &c[0][i * n];
    }
#endif

    /*
     * initalize array and vector
     */
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < n; j++ )
        {
#if defined(USE_MACRO)
            a[IDX(i,j,n)] = (double) random() / RAND_MAX;
            b[IDX(i,j,n)] = (double) random() / RAND_MAX;
            c[IDX(i,j,n)] = 0.0;
#else
            a[i][j] = (double) random() / RAND_MAX;
            b[i][j] = (double) random() / RAND_MAX;
            c[i][j] = 0.0;
#endif
        }
    }

    /*
     * compute product using ikj loop ordering
     */
    double t1 = wtime();
    for ( int i = 0; i < n; i++ )
    {
        for ( int k = 0; k < n; k++ )
        {
            for ( int j = 0; j < n; j++ )
            {
#if defined(USE_MACRO)
                c[IDX(i,j,n)] += a[IDX(i,k,n)] * b[IDX(k,j,n)];
#else
                c[i][j] += a[i][k] * b[k][j];
#endif
            }
        }
    }
    double t2 = wtime();

    /*
     * report results
     */
#if defined(USE_MACRO)
    printf( "1D ARRAY ACCESS:       " );
#else
    printf( "2D ARRAY (W/POINTERS): " );
#endif
    printf( "%f seconds\n", t2 - t1 );

    /*
     * release memory
     */
#if defined(USE_MACRO)
    free( a );
    free( b );
    free( c );
#else
    free( a[0] ); free( a );
    free( b[0] ); free( b );
    free( c[0] ); free( c );
#endif

    return 0;
}
