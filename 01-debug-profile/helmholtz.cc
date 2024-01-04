/*
 * Solve the Helmholtz equation
 *                Laplacian(u) - 0.04u = 0
 * on a rectangular domain using finite differences and a red-black SOR
 * iteration.
 *
 * $Smake: g++ -O2 -Wall -o %F %f
 *
 * To profile add "-pg" to compiler command line and use "-O" instead of "-O2":
 *     g++ -pg -O -Wall -o %F %f
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 * Original Version: October 21, 2008
 * Current Version: September 6, 2010
 *
 */

#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <time.h>

//----------------------------------------------------------------------------
// Returns the number of seconds since some fixed arbitrary time in the past.

double wtime(void)
{
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return double(ts.tv_sec + ts.tv_nsec / 1.0e9);
}

//----------------------------------------------------------------------------
// Dump data array to stdout.  This should only be used when nx and ny are
// small, e.g. less than 20.
//
// Input:
//    double** v:  data array
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction

void showGrid(double** v, int nx, int ny)
{
    printf("------------------------------------------------------------\n"); 
    for (int j = ny - 1; j >= 0; j--)
    {
        for (int i = 0; i < nx; i++)
        {
            printf(" %6.4f", v[i][j]);
        }
        printf("\n");
    }
    printf("------------------------------------------------------------\n"); 
}

//----------------------------------------------------------------------------
// Compute solution of Helmholtz equation Laplacian(u)-0.04u = 0
//
// Input:
//    double x:    x coordinate, 0.0 <= x <= 1.0
//    double y:    y coordinate, 0.0 <= y <= 1.0 (when domain is square)
//
// Returns:
//    double:      value of solution at (x,y)

double solution(double x, double y)
{
    return cosh(x / 5.0) + cosh(y / 5.0);
}

//----------------------------------------------------------------------------
// Initialize boundary of data array using boundary data and interior of
// data array by linear interpolation of boundary data.
//
// Input:
//    double h:    grid spacing
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction
//
// Output:
//    double** u:  data array

void initializeDomain(double** u, double h, int nx, int ny)
{
    const int nxm1 = nx - 1;
    const int nym1 = ny - 1;

    double x, y;

    // top and bottom boundaries
    y = nym1 * h;
    for (int i = 0; i < nx; i++) 
    {
        x = i * h;
        u[i][0] = solution(x, 0.0);
        u[i][nym1] = solution(x, y);
    }

    // left and right boundaries
    x = 1.0;
    for (int j = 0; j < ny; j++)
    {
        y = j * h;
        u[0][j] = solution(0.0, y);
        u[nxm1][j] = solution(x, y);
    }

    // linearly interpolate boundary data across interior of domain;
    // this gives us a head start toward convergence...
    for (int i = 1; i < nxm1; i++)
    {
        for (int j = 1; j < nym1; j++)
        {
            double a = (u[0][j] * (nxm1 - i) + u[nxm1][j] * i) / nxm1;
            double b = (u[i][0] * (nym1 - j) + u[i][nym1] * j) / nym1;
            u[i][j] = 0.5 * (a + b);
        }
    }
}

//----------------------------------------------------------------------------
// Perform red or black SOR sweep for Helmholtz equation
// Laplacian(u) + f * u = g
//
// Input:
//    int color:   either 0 or 1 indicating color of sweep
//    double** u:  data array
//    double f:    coefficient of u
//    double g:    right-hand-side
//    double w:    SOR parameter omega
//    double h:    grid spacing
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction
//
// Output:
//    double** u:  interior red or black values updated

void SORHelmholtz(int color, double** u, double f, double g, double w,
                   double h, int nx, int ny)
{
    const double A = h * h * g;
    const double wB = w / (4.0 - h * h * f);

    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1 + (i + color) % 2; j < ny - 1; j += 2)
        {
            u[i][j] = (1.0 - w) * u[i][j]
                + wB * (u[i][j-1] + u[i-1][j] + u[i][j+1] + u[i+1][j] - A);
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
//    int nx:      number of grid points in x direction
//    int ny:      number of grid points in y direction
//
// Returns:
//    double:      sum_{i,j} | u(i,j)-uhat(i,j) |,  uhat is true solution

double errorNorm(double** u, double h, int nx, int ny)
{
    double sum = 0.0;
    for (int i = 1; i < nx - 1; i++)
    {
        double x = i * h;
        for (int j = 1; j < ny - 1; j++)
        {
            double y = j * h;
            sum += fabs(u[i][j] - solution(x, y));
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

double diffNorm(double** u, double** v, int nx, int ny)
{
    double sum = 0.0;
    for (int i = 1; i < nx - 1; i++)
    {
        for (int j = 1; j < ny - 1; j++)
        {
            sum += fabs(u[i][j] - v[i][j]);
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

void gridCopy(double** u, double** v, int nx, int ny)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            u[i][j] = v[i][j];
        }
    }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    // Get grid dimensions from command line

    if (argc < 3 || argc > 4)
    {
        printf("Usage %s NX NY [ITERATIONS_BETWEEN_CHECKS]\n", argv[0]);
        printf("\nNX and NY are the grid dimensions.\n");
        printf("ITERATIONS_BETWEEN_CHECKS is the number of Jacobi ");
        printf("iterations to\n");
        printf("  perform between convergence checks; default is 1.\n\n");
        exit(EXIT_FAILURE);
    }
    int nx = atoi(argv[1]);    // number of grid points in x direction
    int ny = atoi(argv[2]);    // number of grid points in y direction
    if (nx <= 0 || ny <= 0)
    {
        fprintf(stderr, "Error: both NX and NY must be positive\n");
        exit(EXIT_FAILURE);
    }

    // Get optional command line parameter (default is 1) for the 
    // number of iterations between checks for convergence

    int iterationsBetweenChecks = 1; // default is to check every iteration
    if (argc == 4)
    {
        iterationsBetweenChecks = atoi(argv[3]);
    }
    if (iterationsBetweenChecks < 0)
    {
        fprintf(stderr,
                 "Error: ITERATIONS_BETWEEN_CHECKS must be positive\n");
        exit(EXIT_FAILURE);
    }
    
    double h = 1.0 / (nx - 1); // use horizontal spacing as grid spacing

    // Create the grid.  We want a a constant stride in each dimension
    // This is accomplished by allocating an array of pointers then
    // allocating the full data array to the first pointer.  The
    // remaining pointers are set to point to the start of each "row"
    // of contiguous data in the single linear array.

    double** u = new double* [nx];
    double** v = new double* [nx];
    u[0] = new double [nx * ny];
    v[0] = new double [nx * ny];
    for (int i = 1; i < nx; i++)
    {
        u[i] = &u[0][i * ny];
        v[i] = &v[0][i * ny];
    }

    // Set boundary values and fill interior of domain with initial estimate

    initializeDomain(u, h, nx, ny);

    // Set iteration parameters

    const int maxIter = 10 * nx * ny; // bound on number of SOR iterations
    const double tolerance = 1e-6;    // bound on inf-norm of consecutive solns
    const int red = 0;
    const int black = (red + 1) % 2;

    // Estimate of optimal SOR parameter omega

    const double w = 2.0
        / (1.0 + 0.5 * (sin(M_PI / (nx-1)) + sin(M_PI / (ny-1))));

    // Perform SOR iterations

    double t1 = wtime();

    int k = 0;
    double norm = 2 * tolerance;
    while (norm > tolerance && k++ < maxIter)
    {
        if (k % iterationsBetweenChecks == 0)
        {
            // we're going to a convergence check after this iteration,
            // need to save a copy of the data array
            gridCopy(v, u, nx, ny);
        }

        // perform both red and black SOR sweeps

        SORHelmholtz(red,   u, -0.04, 0.0, w, h, nx, ny);
        SORHelmholtz(black, u, -0.04, 0.0, w, h, nx, ny);

        if (k % iterationsBetweenChecks == 0)
        {
            // compute norm between old and new data values
            norm = diffNorm(u, v, nx, ny);
        }
    }

    double t2 = wtime();

    // done iterating; print out results

    double err = errorNorm(u, h, nx, ny);

#ifdef DEBUG
    showGrid(u, nx, ny);
#endif

    printf("%d iterations done in %e seconds\n", k, t2 - t1);
    printf("difference norm = %e\n", norm);
    printf("error norm      = %e\n", err);

    // Release memory and quit

    delete [] u[0];
    delete [] v[0];
    delete [] u;
    delete [] v;

    return 0;
}
