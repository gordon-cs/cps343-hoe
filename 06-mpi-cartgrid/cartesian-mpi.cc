/*
 * Setup MPI Cartesian communicator
 */

#include "cartesian-mpi.h"   // defines Cartesian_Block structure

/*----------------------------------------------------------------------------
 * Computes the starting and ending displacements for the ith
 * subinterval in an n-element array given that there are m
 * subintervals of approximately equal size.
 *
 * Input:
 *    int n    - length of array (array indexed [0]..[n-1])
 *    int m    - number of subintervals
 *    int i    - subinterval number
 *
 * Output:
 *    int* s   - location to store subinterval starting index
 *    int* e   - location to store subinterval ending index
 *
 * Suppose we want to partition a 100-element array into 3
 * subintervals of roughly the same size.  The following three
 * pairs of calls find the starting and ending indices of each
 * subinterval:
 *   decompose1d(100, 3, 0, &s, &e);  (now s =  0, e = 33)
 *   decompose1d(100, 3, 1, &s, &e);  (now s = 34, e = 66)
 *   decompose1d(100, 3, 2, &s, &e);  (now s = 67, e = 99)
 *
 * The subinterval length can be computed with e - s + 1.
 *
 * Based on the FORTRAN subroutine MPE_DECOMP1D in the file
 * UsingMPI/intermediate/decomp.f supplied with the book
 * "Using MPI" by Gropp et al.  It has been adapted to use
 * 0-based indexing.
 */
static void decompose1d(int n, int m, int i, int* s, int* e)
{
    const int length  = n / m;
    const int deficit = n % m;
    *s =  i * length + (i < deficit ? i : deficit);
    *e = *s + length - (i < deficit ? 0 : 1);
    if ((*e >= n) || (i == m - 1)) *e = n - 1;
}

/*----------------------------------------------------------------------------
 * Create and initalize Cartesian communicator
 *
 * Input:
 *   int num_proc    - number of processes
 *   int NX, NY      - number of grid points in x and y directions
 *   int may_rerank  - nonzero if processes can be reranked
 *
 * Input/Output:
 *   int* rank       - process rank
 *   int* dims       - process grid dimensions (zeros to allow MPI to choose)
 *   int* periodic   - zeros if grid is periodic
 *
 * Output:
 *   MPI_Datatype* x_slice      - type representing horizontal row of data
 *   MPI_Datatype* y_slize      - type representing vertical column of data
 *   Cartesian_Block* halo_grid - parameters for grid block including halo
 *   Cartesian_Block* orig_grid - parameters for grid block without halo
 *
 * Returns
 *   MPI_Comm        - MPI communicator for Cartesian grid
 */
MPI_Comm mpi_cart_setup(int num_proc, int NX, int NY, int may_rerank,
                        int* rank, int* dims, int* periodic,
                        MPI_Datatype* x_slice, MPI_Datatype* y_slice,
                        Cartesian_Block* halo_grid,
                        Cartesian_Block* orig_grid)
{
    MPI_Comm comm2d;    // communicator for Cartesian grid
    int coords[2];      // coordinates of block within Cartesian grid

    // Set up Cartesian grid of processors.  A new communicator is
    // created we get our rank within it.

    MPI_Dims_create(num_proc, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, may_rerank, &comm2d);
    MPI_Cart_get(comm2d, 2, dims, periodic, coords);
    MPI_Comm_rank(comm2d, rank);

    // Figure out who my neighbors are.  The values of left_neighbor,
    // right_neighbor, below_neighbor, and above_neighbor will be set
    // to the rank of the process responsible for the corresponding
    // block relative to the position of the block we are responsible
    // for.  If there is no neighbor in a particular direction (i.e.,
    // we are on a boundary) the returned rank will be MPI_PROC_NULL
    // calls to MPI_sendrecv() will be silently ignored.

    MPI_Cart_shift(comm2d, 0, 1, &(orig_grid->left_neighbor),
                   &(orig_grid->right_neighbor));
    MPI_Cart_shift(comm2d, 1, 1, &(orig_grid->below_neighbor),
                   &(orig_grid->above_neighbor));

    // Figure out the extents and size of my portion of the grid.

    decompose1d(NX, dims[0], coords[0], &(orig_grid->x0), &(orig_grid->x1));
    decompose1d(NY, dims[1], coords[1], &(orig_grid->y0), &(orig_grid->y1));
    orig_grid->nx = orig_grid->x1 - orig_grid->x0 + 1;
    orig_grid->ny = orig_grid->y1 - orig_grid->y0 + 1;

    // Adjust domain parameters to account for inter-domain halo boundary
    // data.  If we have a neighbor in a given direction (rank of neighbor
    // is non-negative) then we need to adjust the starting or ending index.

    halo_grid->left_neighbor  = orig_grid->left_neighbor;
    halo_grid->right_neighbor = orig_grid->right_neighbor;
    halo_grid->above_neighbor = orig_grid->above_neighbor;
    halo_grid->below_neighbor = orig_grid->below_neighbor;
    halo_grid->x0 = orig_grid->x0 - (halo_grid->left_neighbor  >= 0 ? 1 : 0);
    halo_grid->x1 = orig_grid->x1 + (halo_grid->right_neighbor >= 0 ? 1 : 0);
    halo_grid->y0 = orig_grid->y0 - (halo_grid->below_neighbor >= 0 ? 1 : 0);
    halo_grid->y1 = orig_grid->y1 + (halo_grid->above_neighbor >= 0 ? 1 : 0);
    halo_grid->nx = halo_grid->x1 - halo_grid->x0 + 1;
    halo_grid->ny = halo_grid->y1 - halo_grid->y0 + 1;

    // Create datatypes for exchanging x and y slices

    MPI_Type_vector(halo_grid->nx, 1, halo_grid->ny, MPI_DOUBLE, x_slice);
    MPI_Type_commit(x_slice);
    MPI_Type_vector(halo_grid->ny, 1, 1, MPI_DOUBLE, y_slice);
    MPI_Type_commit(y_slice);

    return comm2d;
}
