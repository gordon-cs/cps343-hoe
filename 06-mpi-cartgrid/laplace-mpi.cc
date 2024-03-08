/*
 * $Smake: mpic++ -Wall -O3 -funroll-loops -o %F %f cartesian-mpi.o
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
#include <unistd.h>
#include <mpi.h>
#include "cartesian-mpi.h"

const int    DEFAULT_DIMENSION = 200;
const int    DEFAULT_ITERATIONS = 500000;
const int    DEFAULT_ITERATION_STRIDE = 100;
const double DEFAULT_TOLERANCE = 1.0e-6;

#define IDX(i,j,stride) ((i)*(stride)+(j)) // row major

//---------------------------------------------------------------------------

// Print out 2D grid data
void show_grid(
    Cartesian_Block* grid,    // parameters for local grid block
    double* v                 // local grid data
    )
{
    const int nx = grid->nx;
    const int ny = grid->ny;

    printf("--------------------------------------------------------\n");
    for (int j = ny - 1; j >= 0; j--)
    {
        for (int i = 0; i < nx; i++)
        {
            printf(" %6.2f", v[IDX(i,j,ny)]);
        }
        printf("\n");
    }
    printf("--------------------------------------------------------\n");
    fflush(stdout);
}

//---------------------------------------------------------------------------

// Display local grid data in process rank order
void dump_grid(
    int my_rank,              // rank of this process
    int num_proc,             // number of processes
    Cartesian_Block* grid,    // parameters for local grid block
    double* u,                // local grid data
    MPI_Comm comm2d           // Cartesian communicator
    )
{
    for (int i = 0; i < num_proc; i++)
    {
        if (my_rank == i)
        {
            printf("Rank: %d\n", my_rank);
            show_grid(grid, u);
        }
        // barrier and delay so output comes out in process rank order
        MPI_Barrier(comm2d);
        usleep(10000);
    }
}

//---------------------------------------------------------------------------

// Initialize 2D grid with zeros
void init_grid(
    Cartesian_Block* grid,    // parameters for local grid block
    double* u                 // local grid data
    )
{
    for (int i = 0; i < grid->nx * grid->ny; i++) u[i] = 0.0;
}

//---------------------------------------------------------------------------

// Set boundary values on 2D grid
void impose_boundary_conditions(
    double x0, double xn,     // left and right endpoints of domain
    double y0, double yn,     // bottom and top endpoints of domain
    int NX, int NY,           // dimensions of full grid
    Cartesian_Block* grid,    // parameters for local grid block
    double* u                 // local grid data
    )
{
    const int nx = grid->nx;
    const int ny = grid->ny;

    // set top and bottom boundary values
    for (int i = 0; i < nx; i++)
    {
        double x = x0 + (xn - x0) * (grid->x0 + i) / (NX - 1);
        if (grid->below_neighbor < 0) u[IDX(i,0,ny)]    = 0.0;
        if (grid->above_neighbor < 0) u[IDX(i,ny-1,ny)] = 1.0 + x * (1.0 - x);
    }

    // set left and right boundary values
    for (int j = 0; j < ny; j++)
    {
        double y = y0 + (yn - y0) * (grid->y0 + j) / (NY - 1);
        if (grid->left_neighbor  < 0) u[IDX(0,j,ny)]    = y;
        if (grid->right_neighbor < 0) u[IDX(nx-1,j,ny)] = y * y;
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
            v[IDX(i,j,ny)] = (u[IDX(i-1,j,ny)] + u[IDX(i+1,j,ny)]
                              + u[IDX(i,j-1,ny)] + u[IDX(i,j+1,ny)]) / 4.0;
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
            s += fabs(v[IDX(i,j,ny)] - u[IDX(i,j,ny)]);
        }
    }
    return s;
}

/*----------------------------------------------------------------------------
 * Exchange halo boundary data between 2D grid blocks.
 *
 * Input:
 *   double* u             - 2-D array holding single block of grid data
 *   Cartesian_Block* grid - pointer to grid block data
 *   MPI_Datatype x_slice  - type representing horizontal row of data
 *   MPI_Datatype y_slice  - type representing vertical column of data
 *   MPI_Comm comm         - communicator
 *
 * Output:
 *   double* u             - 2-D array holding grid block with update halo
 */
void exchange_halo_data(
    double* u,               // local grid data
    Cartesian_Block* grid,   // parameters for local grid block
    MPI_Datatype x_slice,    // horizontal slice of local grid
    MPI_Datatype y_slice,    // vertical slice of local grid
    MPI_Comm comm            // Cartesian communicator
    )
{
    const int tag = 0;
    const int nx = grid->nx;
    const int ny = grid->ny;
    
    // Send top row of my data to bottom halo of neighbor above me and
    // receive bottom row of same neighbor's data into my top halo
    MPI_Sendrecv(&u[IDX(0,ny-2,ny)], 1, x_slice, grid->above_neighbor, tag,
                 &u[IDX(0,0,ny)],    1, x_slice, grid->below_neighbor, tag,
                 comm, MPI_STATUS_IGNORE);
    
    // Send bottom row of my data to top halo of neighbor below me and
    // receive top row of same neighbor's data into my bottom halo
    MPI_Sendrecv(&u[IDX(0,1,ny)],    1, x_slice, grid->below_neighbor, tag,
                 &u[IDX(0,ny-1,ny)], 1, x_slice, grid->above_neighbor, tag,
                 comm, MPI_STATUS_IGNORE);

    // Send right column of my data to left halo of neighbor to my right
    // and receive left column of same neighbor's data into my right halo
    MPI_Sendrecv(&u[IDX(nx-2,0,ny)], 1, y_slice, grid->right_neighbor, tag,
                 &u[IDX(0,0,ny)],    1, y_slice, grid->left_neighbor,  tag,
                 comm, MPI_STATUS_IGNORE);

    // Send left column of my data to right halo of neighbor to my left
    // and receive right column of same neighbor's data into my left halo
    MPI_Sendrecv(&u[IDX(1,0,ny)],    1, y_slice, grid->left_neighbor,  tag,
                 &u[IDX(nx-1,0,ny)], 1, y_slice, grid->right_neighbor, tag,
                 comm, MPI_STATUS_IGNORE);
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

    int num_proc;              // number of processes
    int my_rank;               // rank of this process
    int dims[2] = {0,0};       // Cartesian block array dimensions
    int periodic[2] = {0,0};   // zero values indicate boundaries not periodic
    int may_rerank = 1;        // nonzero to indicate process may be reranked
    Cartesian_Block halo_grid; // parameters for local grid block w/halo
    Cartesian_Block orig_grid; // parameters for local grid block wo/halo
    MPI_Datatype x_slice;      // type representing horizontal row of data
    MPI_Datatype y_slice;      // type representing vertical row of data
    MPI_Comm comm2d;           // Cartesian communicator

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // Process command line
    int c;
    opterr = 0;  // suppress getopt error messages
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
                if (my_rank == 0) usage(argv[0]);
                MPI_Finalize();
                return 0;
        }
    }

    // Set up Cartesian communicator
    comm2d = mpi_cart_setup(num_proc, nx, ny, may_rerank, &my_rank, dims,
                            periodic, &x_slice, &y_slice, &halo_grid,
                            &orig_grid);

    // Allocate memory for grid and copy of grid
    double* u = new double [halo_grid.nx * halo_grid.ny];
    double* v = new double [halo_grid.nx * halo_grid.ny];

    // Prepare grid
    init_grid(&halo_grid, u);
    impose_boundary_conditions(0.0, 1.0, 0.0, 1.0, nx, ny, &halo_grid, u);
    copy_grid(v, u, halo_grid.nx, halo_grid.ny);

    if (verbosity > 1) dump_grid(my_rank, num_proc, &halo_grid, u, comm2d);

    // Do Jacobi iterations until convergence or too many iterations
    double t0 = MPI_Wtime();
    int k = 0;
    double alpha = 2 * tolerance;
    while (k++ < max_iter && alpha > tolerance)
    {
        exchange_halo_data(u, &halo_grid, x_slice, y_slice, comm2d);
        jacobi_sweep(v, u, halo_grid.nx, halo_grid.ny);

        // check to see how much we've changed from prior estimate
        if (k % iterations_between_checks == 0)
        {
            double my_alpha = norm(u, v, halo_grid.nx, halo_grid.ny);
            MPI_Allreduce(&my_alpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, comm2d);
            if (verbosity > 0 && my_rank == 0) printf("%6d %e\n", k, alpha);
            if (verbosity > 1)
                dump_grid(my_rank, num_proc, &halo_grid, u, comm2d);
        }
        copy_grid(u, v, halo_grid.nx, halo_grid.ny);
    }
    double t1 = MPI_Wtime();

    // Report results
    if (my_rank == 0)
    {
        printf("Iterations: %d Difference Norm: %12.6e  ", --k, alpha);
        printf("Time: %f seconds\n", t1 - t0);
    }

    // All done - clean up and exit
    delete [] u;
    delete [] v;

    MPI_Finalize();
    return 0;
}
