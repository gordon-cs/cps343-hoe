/*
 * $Smake: g++ -O2 -o %F_2D %f -lrt; g++ -DUSE_MACRO -O2 -o %F_1D %f -lrt
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

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>

using namespace std;

#if !defined(N)
# define N 1000
#endif

#if defined(USE_MACRO)
# define idx(i,j,rowlen) ((i)*(rowlen)+j)
#endif

//----------------------------------------------------------------------------
// Returns the number of seconds since some fixed arbitrary time in the past

double wtime( void )
{
    timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return double( ts.tv_sec + ts.tv_nsec / 1.0e9 );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc, char *argv[] )
{
    const int n = N;

    // allocate contiguous memory for matrices

#if defined(USE_MACRO)
    double* a = new double [n * n];
    double* b = new double [n * n];
    double* c = new double [n * n];
#else
    double** a = new double* [n];
    double** b = new double* [n];
    double** c = new double* [n];
    a[0] = new double [n * n];
    b[0] = new double [n * n];
    c[0] = new double [n * n];
    for ( int i = 1; i < n; i++ )
    {
	a[i] = &a[0][i * n];
	b[i] = &b[0][i * n];
	c[i] = &c[0][i * n];
    }
#endif

    // initalize array and vector

    for ( int i = 0; i < n; i++ )
    {
	for ( int j = 0; j < n; j++ )
	{
#if defined(USE_MACRO)
	    a[idx(i,j,n)] = (double) random() / RAND_MAX;
	    b[idx(i,j,n)] = (double) random() / RAND_MAX;
	    c[idx(i,j,n)] = 0.0;
#else
	    a[i][j] = (double) random() / RAND_MAX;
	    b[i][j] = (double) random() / RAND_MAX;
	    c[i][j] = 0.0;
#endif
	}
    }

    // compute product

    double t1 = wtime();
    for ( int i = 0; i < n; i++ )
    {
	for ( int k = 0; k < n; k++ )
	{
	    for ( int j = 0; j < n; j++ )
	    {
#if defined(USE_MACRO)
		c[idx(i,j,n)] += a[idx(i,k,n)] * b[idx(k,j,n)];
#else
		c[i][j] += a[i][k] * b[k][j];
#endif
	    }
	}
    }
    double t2 = wtime();
#if defined(USE_MACRO)
    cout << "1D ARRAY ACCESS:       ";
#else
    cout << "2D ARRAY (W/POINTERS): ";
#endif
    cout << t2 - t1 << " seconds" << endl;

    // release memory

#if defined(USE_MACRO)
    delete [] a;
    delete [] b;
    delete [] c;
#else
    delete [] a[0]; delete [] a;
    delete [] b[0]; delete [] b;
    delete [] c[0]; delete [] c;
#endif

    return 0;
}
