/*
 * Solve the Helmholtz equation
 *                Laplacian(u) - 0.04u = 0
 * on a rectangular domain using finite differences and Jacobi iterations.
 *
 * $Smake: mpic++ -O2 -Wall -o %F %f
 *
 * To profile add "-pg" to compiler command line and use "-O" instead of "-O2":
 *     g++ -pg -O -Wall -o %F %f
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 * Original Version: February 14, 2013 (based on helmholtz.cc)
 *
 */

#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <mpi.h>

//----------------------------------------------------------------------------
// Compute solution of Helmholtz equation Laplacian(u)-0.04u = 0
//
// Input:
//    double x:    x coordinate, 0.0 <= x <= 1.0
//    double y:    y coordinate, 0.0 <= y <= 1.0 (when domain is square)
//
// Returns:
//    double:      value of solution at (x,y)

double solution( double x, double y )
{
    return cosh( x / 5.0 ) + cosh( y / 5.0 );
}

//----------------------------------------------------------------------------
// Display program usage statement.
//
// Input:
//    char* program:  program name string

void usage( char* program )
{
    printf( "Usage: %s NX NY [DIMX DIMY [ITERATIONS_BETWEEN_CHECKS]]\n\n",
             program );
    printf( "- NX and NY are the grid dimensions.\n" );
    printf( "- DIMX and DIMY are the process grid dimensions.\n" );
    printf( "- ITERATIONS_BETWEEN_CHECKS is number of SOR iterations to\n" );
    printf( "    perform between convergence checks; default is 1.\n\n" );
}

//----------------------------------------------------------------------------
// Dump data array to stdout.  This should only be used when nx and ny are
// small, e.g. less than 20.
//
// Input:
//    double** v:  data array
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction

void showGrid( double** v, int nx, int ny )
{
    printf( "------------------------------------------------------------\n" ); 
    for ( int j = ny - 1; j >= 0; j-- )
    {
        for ( int i = 0; i < nx; i++ )
        {
            printf( " %6.4f", v[i][j] );
        }
        printf( "\n" );
    }
    printf( "------------------------------------------------------------\n" ); 
}

//----------------------------------------------------------------------------

void decompose1D( int n, int nproc, int rank, int& s, int& e )
// This function is based on the FORTRAN subroutine MPE_DECOMP1D in the
// file UsingMPI/intermediate/decomp.f supplied with the book Using MPI
// by Gropp et al.  It has been adapted to use 0-based indexing.
{
    int nlocal  = n / nproc;
    int deficit = n % nproc;

    s = rank * nlocal + ( rank < deficit ? rank : deficit );
    e = s + nlocal - ( rank < deficit ? 0 : 1 );

    if ( ( e >= n ) || ( rank == nproc - 1 ) ) e = n - 1;
}

//----------------------------------------------------------------------------

void decompose2D( MPI_Comm comm, int nx, int ny,
                  int& sx, int& ex, int& sy, int& ey )
{
    int dims[2];
    int periods[2];
    int coords[2];
    
    MPI_Cart_get( comm, 2, dims, periods, coords );
    decompose1D( nx, dims[0], coords[0], sx, ex );
    decompose1D( ny, dims[1], coords[1], sy, ey );
}

//----------------------------------------------------------------------------

void exchangeSlices( double** u, int nx, int ny, 
                     int up, int down, int left, int right, 
                     MPI_Datatype xSlice, MPI_Datatype ySlice, MPI_Comm comm )
{
    enum { T1, T2, T3, T4 };

    // Exchange x-slices with my top and bottom neighbors
    MPI_Sendrecv( &u[0][ny-2], 1, xSlice, up,   T1,
                  &u[0][0],    1, xSlice, down, T1,
                  comm, MPI_STATUS_IGNORE );
    MPI_Sendrecv( &u[0][1],    1, xSlice, down, T2,
                  &u[0][ny-1], 1, xSlice, up,   T2,
                  comm, MPI_STATUS_IGNORE );

    // Exchange y-slices with my left and right neighbors
    MPI_Sendrecv( &u[nx-2][0], 1, ySlice, right, T3,
                  &u[0][0],    1, ySlice, left,  T3,
                  comm, MPI_STATUS_IGNORE );
    MPI_Sendrecv( &u[1][0],    1, ySlice, left,  T4,
                  &u[nx-1][0], 1, ySlice, right, T4,
                  comm, MPI_STATUS_IGNORE );
}

//----------------------------------------------------------------------------
// Initialize subdomain boundary of data array using boundary data and
// interior of data array by linear interpolation of boundary data.
//
// Input:
//    double h:    grid spacing
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction
//
// Output:
//    double** u:  data array

void initializeDomain( double** u, double x0, double xn, double y0, double yn,
                       double h, int nx, int ny, int sx, int sy,
                       int up, int down, int left, int right )
{
    // start by initializing everything to zero

    for ( int i = 0; i < nx; i++ )
    {
        for ( int j = 0; j < ny; j++ )
        {
            u[i][j] = 0.0;
        }
    }

    // if down, up, left, or right are less than zero then the lower, upper,
    // left, or right edge respectively is a boundary edge and should have
    // boundary values set.

    if ( down < 0 )
    {
        for ( int i = 0; i < nx; i++ )
        {
            u[i][0] = solution( ( sx + i ) * h, y0 );
        }
    }
    if ( up < 0 )
    {
        for ( int i = 0; i < nx; i++ )
        {
            u[i][ny-1] = solution( ( sx + i ) * h, yn );
        }
    }
    if ( left < 0 )
    {
        for ( int j = 0; j < ny; j++ )
        {
            u[0][j] = solution( x0, ( sy + j ) * h );
        }
    }
    if ( right < 0 )
    {
        for ( int j = 0; j < ny; j++ )
        {
            u[nx-1][j] = solution( xn, ( sy + j ) * h );
        }
    }
}

//----------------------------------------------------------------------------
// Perform Jacobi sweep for Helmholtz equation
// Laplacian(u) + f * u = g
//
// Input:
//    double** u:  data array
//    double f:    coefficient of u
//    double g:    right-hand-side
//    double h:    grid spacing
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction
//
// Output:
//    double** v:  updated data array

void jacobiSweep( double** v, double** u, double f, double g, double h,
                  int nx, int ny )
{
    const double A = h * h * g;
    const double B = 1.0 / ( 4.0 - h * h * f );

    for ( int i = 1; i < nx - 1; i++ )
    {
        for ( int j = 1; j < ny - 1; j++ )
        {
            v[i][j] = B * ( u[i][j-1] + u[i-1][j] + u[i][j+1] + u[i+1][j] - A );
        }
    }
}

//----------------------------------------------------------------------------
// Computes actual L-infinity error norm between data in u and true solution
// Of course, this is only possible if the true solution is known...
//
// Input:
//    double** u:  data array
//    double h:    grid spacing
//    int nx, ny:  number of grid points in x and y directions
//    int sx, sy;  initial grid point indicies in x and y direction
//
// Returns:
//    double:      sum_{i,j} | u(i,j)-uhat(i,j) |,  uhat is true solution

double errorNorm( double** u, double h, int nx, int ny, int sx, int sy )
{
    double sum = 0.0;
    for ( int i = 1; i < nx - 1; i++ )
    {
        double x = ( sx + i ) * h;
        for ( int j = 1; j < ny - 1; j++ )
        {
            double y = ( sy + j ) * h;
            sum += fabs( u[i][j] - solution( x, y ) );
        }
    }
    return sum;
}

//----------------------------------------------------------------------------
// Computes L-infinity norm between data in u and v.
//
// Input:
//    double** u:  data array
//    double** v:  second data array
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction
//
// Returns:
//    double:      sum_{i,j} | u(i,j)-v(i,j) |

double diffNorm( double** u, double** v, int nx, int ny )
{
    double sum = 0.0;
    for ( int i = 1; i < nx - 1; i++ )
    {
        for ( int j = 1; j < ny - 1; j++ )
        {
            sum += fabs( u[i][j] - v[i][j] );
        }
    }
    return sum;
}

//----------------------------------------------------------------------------
// Copy data from one array to another
//
// Input:
//    double** v:  data array
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction
//
// Output:
//    double** u:  copy of data in v

void gridCopy( double** u, double** v, int nx, int ny )
{
    for ( int i = 0; i < nx; i++ )
    {
        for ( int j = 0; j < ny; j++ )
        {
            u[i][j] = v[i][j];
        }
    }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main( int argc, char *argv[] )
{
    int verbosity = 1;
    int numberOfProcesses;
    int rank;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &numberOfProcesses );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // Set up cartesian grid of processors.  A new communicator is created
    // and rank is updated.

    MPI_Comm c2d;
    int dims[2]    = { 0, 0 };  // Allow MPI to determine grid dimensions
    int periods[2] = { 0, 0 };  // Non-periodic
    int reorder    = 1;         // Allow processes to be re-ranked

    // Process command line: two arguments (plus program name) are required,
    // third and fourth are optional but must both appear if either does.
    // A fifth argument is used if present; all others are ignored.

    if ( argc < 3 || argc == 4 )
    {
        if ( rank == 0 ) usage( argv[0] );
        MPI_Finalize();
        exit( EXIT_FAILURE );
    }

    int NX = atoi( argv[1] );
    int NY = atoi( argv[2] );
    if ( NX <= 0 || NY <= 0 )
    {
        if ( rank == 0 )
        {
            fprintf( stderr, "Error: both NX and NY must be positive.\n" );
        }
        MPI_Finalize();
        exit( EXIT_FAILURE );
    }

    if ( argc > 4 )
    {
        dims[0] = atoi( argv[3] );
        dims[1] = atoi( argv[4] );
    }

    int iterStep = 1;
    if ( argc > 5 )
    {
        iterStep = atoi( argv[5] );
    }
    if ( iterStep <= 0 )
    {
        if ( rank == 0 )
        {
            fprintf( stderr,
                     "Warning: ITERATIONS_BETWEEN_CHECKS must be positive\n" );
            fprintf( stderr, "Resetting to 1\n" );
        }
        iterStep = 1;
    }

    // Set up Cartesian grid

    if ( verbosity > 0 && rank == 0 )
    {
        printf( "Iterations between checks = %d\n", iterStep );
        printf( "Requested dims[] = %d %d\n", dims[0], dims[1] );
    }

    MPI_Dims_create( numberOfProcesses, 2, dims );
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, reorder, &c2d );
    MPI_Comm_rank( c2d, &rank );

    if ( verbosity > 0 )
    {
        int coords[2];
        MPI_Cart_get( c2d, 2, dims, periods, coords );
        if ( rank == 0 )
        {
            printf( "Actual dims[]    = %d %d\n", dims[0], dims[1] );
        }
        printf( "Rank %2d: coords[]   = %d %d\n", rank, coords[0], coords[1] );
    }

    // Figure out who my neighbors are

    int up, down, left, right;
    MPI_Cart_shift( c2d, 0, 1, &left, &right );
    MPI_Cart_shift( c2d, 1, 1, &down, &up );

    // Figure out the size of my portion of the grid.  Note that we adjust
    // the starting and ending grid indices along boundaries with neighboring
    // subdomains; this provides room to hold the internal "boundary" data
    // that will be received from the neighbor.

    int sx, ex, sy, ey;
    decompose2D( c2d, NX, NY, sx, ex, sy, ey );
    if ( left  >= 0 ) sx--;
    if ( right >= 0 ) ex++;
    if ( down  >= 0 ) sy--;
    if ( up    >= 0 ) ey++;

    int nx = ex - sx + 1;
    int ny = ey - sy + 1;

    double h = 1.0 / ( NX - 1 );

    // Create datatypes for exchanging x and y slices

    MPI_Datatype xSlice;
    MPI_Type_vector( nx, 1, ny, MPI_DOUBLE, &xSlice );
    MPI_Type_commit( &xSlice );

    MPI_Datatype ySlice;
    MPI_Type_vector( ny, 1, 1, MPI_DOUBLE, &ySlice );
    MPI_Type_commit( &ySlice );

    // Create my portion of the grid.  For the exchange to work properly
    // we must have a constant stride in each dimension.  This is
    // accomplished by allocating an array of pointers then allocating
    // the full data array to the first pointer.  The remaining pointers
    // are set to point to the start of each "row" of contiguous data
    // in the single linear array.

    double** u = new double* [nx];
    double** v = new double* [nx];
    u[0] = new double [nx * ny];
    v[0] = new double [nx * ny];
    for ( int i = 1; i < nx; i++ )
    {
        u[i] = &u[0][i * ny];
        v[i] = &v[0][i * ny];
    }

    // Set boundary values and fill interior of domain with initial estimate.
    // We need both u and v to have the boundary values so we copy u to v.

    initializeDomain( u, 0.0, 1.0, 0.0, ( NY - 1 ) * h,
                      h, nx, ny, sx, sy, up, down, left, right );

//    gridCopy( v, u, nx, ny );

    // Set iteration parameters

    const long maxIter = 10 * NX * NY; // max number of Jacobi iterations
    const double tolerance = 1e-6;     // bound on inf-norm of consecutive solns

    // Perform Jacobi iterations

    double t1 = MPI_Wtime();

    long k = 0;
    double norm = 2 * tolerance;

    while ( norm > tolerance && k++ < maxIter )
    {
        exchangeSlices( u, nx, ny, up, down, left, right, xSlice, ySlice, c2d );
        
        // Do Jacobi iteration, updating v using values in u
        jacobiSweep( v, u, -0.04, 0.0, h, nx, ny );

        if ( k % iterStep == 0 )
        {
            // compute norm between old and new data values
            double my_norm = diffNorm( u, v, nx, ny );
            MPI_Allreduce( &my_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, c2d );
            if ( verbosity > 1 && rank == 0 ) printf( "%6ld: %e\n", k, norm );
        }

        gridCopy( u, v, nx, ny );
    }

    double t2 = MPI_Wtime();

    // done iterating; print out results

    double err;
    double myErr = errorNorm( u, h, nx, ny, sx, sy );
    MPI_Allreduce( &myErr, &err, 1, MPI_DOUBLE, MPI_SUM, c2d );

#ifdef DEBUG
    int junk = 0;    
    if ( rank > 0 )
    {
        // wait for next lower rank process to send to me
        MPI_Recv( &junk, 1, MPI_INT, rank - 1, 99, c2d, MPI_STATUS_IGNORE );
    }
    showGrid( u, nx, ny );
    if ( rank < numberOfProcesses - 1 )
    {
        // send to next higher rank process
        MPI_Send( &junk, 1, MPI_INT, rank + 1, 99, c2d );
    }
#endif

    if ( rank == numberOfProcesses - 1 )
    {
        printf( "%ld iterations done in %g seconds\n", k, t2 - t1 );
        printf( "difference norm = %e\n", norm );
        printf( "error norm      = %e\n", err );
    }

    // Release memory and quit

    delete [] u[0];
    delete [] v[0];
    delete [] u;
    delete [] v;

    MPI_Type_free( &xSlice );
    MPI_Type_free( &ySlice );

    MPI_Finalize();

    return 0;
}
