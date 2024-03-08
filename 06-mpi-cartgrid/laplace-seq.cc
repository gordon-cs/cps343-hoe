/*
 * $Smake: g++ -Wall -O3 -funroll-loops -o %F %f
 *
 * Solves the Laplace Equation uxx + uyy = 0 on the unit square with
 * boundary conditions u(x,0)=0, u(x,1)=1+x(1-x), u(0,y)=y, u(1,y)=y^2.
 * The Jacobi Method with second order centered finite differences are
 * used.
 * 
 *            | u(x,1)=1+x(1-x)
 *      (0,1) |---------------- (1,1)
 *            |               |
 *            |               |
 *            |               |
 * u(0,y) = y |               | u(1,y) = y^2
 *            |               |
 *            |               |
 *            |               |
 *            ---------------------
 *          (0,0)           (1,0)
 *                 u(x,0)=0
 *
 * A one-dimensional array is used to store the two-dimensional domain
 * and the macro IDX is used to index the array using two indices.
 */

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cmath>
#include "wtime.c"

const int    DEFAULT_DIMENSION = 200;
const int    DEFAULT_ITERATIONS = 500000;
const int    DEFAULT_ITERATION_STRIDE = 100;
const double DEFAULT_TOLERANCE = 1.0e-6;

#define IDX(i,j) ((i)*(ny)+(j)) // presumes variable ny contains grid height

//---------------------------------------------------------------------------

// Initialize 2D grid with zeros
void init_grid(
    int nx,     // number of x grid points 
    int ny,     // number of y grid points
    double* u   // grid data
    )
{
    for (int i = 0; i < nx * ny; i++) u[i] = 0.0;
}

//---------------------------------------------------------------------------

// Set boundary values on 2D grid
void impose_boundary_conditions(
    double x0, double xn,  // left and right endpoints of domain
    double y0, double yn,  // bottom and top endpoints of domain
    int nx,                // number of x grid points 
    int ny,                // number of y grid points
    double* u              // grid data
    )
{
    // set top and bottom boundary values
    for (int i = 0; i < nx; i++)
    {
        double x = x0 + (xn - x0) * i / (nx - 1);
        u[IDX(i,0)]    = 0.0;                    // f(x,y0)
        u[IDX(i,ny-1)] = 1.0 + x * (1.0 - x);  // f(x,y1)
    }

    // set left and right boundary values
    for (int j = 0; j < ny; j++)
    {
        double y = y0 + (yn - y0) * j / (ny - 1);
        u[IDX(0,j)]    = y;       // f(x0,y)
        u[IDX(nx-1,j)] = y * y;   // f(x1,y)
    }
}

//---------------------------------------------------------------------------

// Perform a single Jacobi Sweep on grid u storing results in grid v
void jacobi_sweep(
    double* v,  // updated grid data
    double* u,  // grid data
    int nx,     // number of x grid points 
    int ny      // number of y grid points
    )
{
    for (int i = 1; i < nx-1; i++)
    {
        for (int j = 1; j < ny-1; j++)
        {
            v[IDX(i,j)] = (u[IDX(i-1,j)] + u[IDX(i+1,j)]
                           + u[IDX(i,j-1)] + u[IDX(i,j+1)]) / 4.0;
        }
    }
}

//---------------------------------------------------------------------------

// Copy grid data
void copy_grid(
    double* v,  // destination grid data
    double* u,  // source grid data
    int nx,     // number of x grid points 
    int ny      // number of y grid points
    )
{
    for (int i = 0; i < nx * ny; i++) v[i] = u[i];
}

//---------------------------------------------------------------------------

// Compute L-infinity norm between interior values of u and v
double norm(
    double* v,  // updated grid data
    double* u,  // original grid data
    int nx,     // number of x grid points 
    int ny      // number of y grid points
    )
{
    double s = 0.0;
    for (int i = 1; i < nx-1; i++)
    {
        for (int j = 1; j < ny-1; j++)
        {
            s += fabs(v[IDX(i,j)] - u[IDX(i,j)]);
        }
    }
    return s;
}

//----------------------------------------------------------------------------

// Display program usage statement
void usage(
    char* program   //program name string
    )
{
    printf("Usage: %s ", program);
    printf("[-n N] [-e TOL] [-m MAXITER] [-s ITERATION_STRIDE] [-v]\n");
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    int nx = DEFAULT_DIMENSION;
    int ny = DEFAULT_DIMENSION;
    int max_iter = DEFAULT_ITERATIONS;
    double tolerance = DEFAULT_TOLERANCE;
    int iterations_between_checks = DEFAULT_ITERATION_STRIDE;
    int verbosity = 0;

    // Process command line
    int c;
    while ((c = getopt(argc, argv, "e:hm:n:s:v")) != -1)
    {
        switch (c)
        {
            case 'e':
                tolerance = atof(optarg);
                if (tolerance <= 0.0) tolerance = DEFAULT_TOLERANCE;
                break;
            case 'n':
                nx = ny = atoi(optarg);
                if (nx <= 0) nx = ny = DEFAULT_DIMENSION;
                break;
            case 'm':
                max_iter = atoi(optarg);
                if (max_iter <= 0) max_iter = DEFAULT_ITERATIONS;
                break;
            case 's':
                iterations_between_checks = atoi(optarg);
                if (iterations_between_checks <= 0)
                    iterations_between_checks = DEFAULT_ITERATION_STRIDE;
                break;
            case 'v':
                verbosity++;
                break;
            case 'h':
            default:
                usage(argv[0]);
                return 0;
        }
    }

    // Allocate memory for grid and copy of grid
    double* u = new double [nx * ny];
    double* v = new double [nx * ny];

    // Prepare grid
    init_grid(nx, ny, u);
    impose_boundary_conditions(0.0, 1.0, 0.0, 1.0, nx, ny, u);
    copy_grid(v, u, nx, ny);

    // Do Jacobi iterations until convergence or too many iterations
    double t0 = wtime();
    int k = 0;
    double alpha = 2 * tolerance;
    while (k++ < max_iter && alpha > tolerance)
    {
        jacobi_sweep(v, u, nx, ny);

        // check to see how much we've changed from prior estimate
        if (k % iterations_between_checks == 0)
        {
            alpha = norm(u, v, nx, ny);
            if (verbosity > 0) printf("%6d %e\n", k, alpha);
        }
        copy_grid(u, v, nx, ny);
    }
    double t1 = wtime();

    // Report results
    printf("Iterations: %d  Difference Norm: %12.6e  ", --k, alpha);
    printf("Time: %f seconds\n", t1 - t0);

    // All done - clean up and exit
    delete [] u;
    delete [] v;

    return 0;
}
