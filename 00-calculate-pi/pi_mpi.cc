/*
 * Uses MPI
 *
 * Compute approximation of pi by evaluating the integral of 4/(1+x^2) on
 * the interval [0,1] using the midpoint rule.
 *
 * $Smake: mpic++ -Wall -O2 -o %F %f
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * This program is released into the public domain.
 */

#include <iostream>
#include <mpi.h>

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
    const double a = 0.0;
    const double b = 1.0;

    int numberOfProcesses;
    int rank;
    double mySum;
    double t1, t2;
    double pi;

    // initialize MPI

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &numberOfProcesses );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // determine the portion of the interval to compute

    double delta = ( b - a ) / double( numberOfProcesses );
    double start = rank * delta;
    double end = ( rank + 1 ) * delta;

    // compute pi by computing local sums and reducing to single value

    t1 = MPI_Wtime();
    mySum = integrate( start, end, num_intervals / numberOfProcesses );
    MPI_Reduce( &mySum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    t2 = MPI_Wtime();

    // display result

    if ( rank == 0 )
    {
        const int opsPerStep = 7;
        double gflops = opsPerStep * ( num_intervals / ( t2 - t1 ) ) / 1.0e+9;

        cout.precision( 16 );
        cout << "pi = " << pi;
        cout.precision( 4 );
        cout << " computed in " << t2 - t1 << " seconds; rate = "
             << gflops << " GFLOPS" << endl;
    }

    // all done

    MPI_Finalize();

    return 0;
}
