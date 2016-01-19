/*
 * Uses OpenMP: dynamic load balancing version
 *
 * Compute approximation of pi by evaluating the integral of 4/(1+x^2) on
 * the interval [0,1] using the midpoint rule.
 *
 * $Smake: g++ -Wall -O2 -fopenmp -o %F %f
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * This program is released into the public domain.
 */

#include <iostream>
#include <omp.h>

//----------------------------------------------------------------------------

// Evalute the integral of 4/(1+x^2) on [a,b] using the midpoint rule.

double integrate( double a, double b, long n )
{
    const double dx = ( b - a ) / double( n );
    double sum = 0.0;

    for ( long i = 0; i < n; i++ )
    {
        double x = a + ( i + 0.5 ) * dx;
	sum += 4.0 / ( 1.0 + x * x );
    }
    return sum * dx;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    using namespace std;
    const long num_intervals = 400000000L;
    const int numSections = 100;
    const double a = 0.0;
    const double b = 1.0;
    const double subIntervalLength = ( b - a ) / double( numSections );

    double t1, t2;
    double pi = 0.0;
    int section = 1;

    // compute pi

    t1 = omp_get_wtime();
    #pragma omp parallel default(shared) reduction (+:pi)
    {
	double x0, x1;
	int i;
	#pragma omp critical
	    i = section++;
	while ( i <= numSections )
	{
	    x0 = a + subIntervalLength * ( i - 1 );
	    x1 = a + subIntervalLength * i;
	    pi += integrate( x0, x1, num_intervals / numSections );
	    #pragma omp critical
		i = section++;
	}
    }
    t2 = omp_get_wtime();

    // display result

    const int opsPerStep = 7;
    double gflops = opsPerStep * ( num_intervals / ( t2 - t1 ) ) / 1.0e+9;

    cout.precision( 16 );
    cout << "pi = " << pi;
    cout.precision( 4 );
    cout << " computed in " << t2 - t1 << " seconds; rate = "
	 << gflops << " GFLOPS" << endl;

    // all done

    return 0;
}
