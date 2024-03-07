/*
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 * Written: August 2012.
 * Revised: April 2013, March 2016, March 2018, March 2020.
 *
 * $Smake: mpic++ -Wall -O2 -o %F %f -lhdf5
 *
 * This program demonstrates how to set up and use a Cartesian grid when data
 * along the boundary between neighboring regions must be exchanged.  It also
 * illustrates the use of MPI_Gather() and how to write a data file in
 * parallel using HDF5.
 *
 * Suppose a rectangular 10x10 two-dimensional domain is to be partitioned
 * into four subdomains.  Often computation on the domain requires the
 * stencil
 *                 *
 *                 |
 *             *---X---*
 *                 |
 *                 *
 *
 * when updating grid points that are on the boundary of a subdomain
 * that are adjacent to another subdomain.  This can be facilitated by
 * "adding" additional rows and columns to subdomains that will hold
 * copies of the interior boundary grid points from adjacent subdomains.
 * These are often called "halo" or "ghost" boundaries.
 *
 * Here is the unpartitioned grid, showing the boundary grid points ("*")
 * whose data is used but unchanged, and the grid points on the domain
 * interior ("O") that are updated by computation.
 *
 *      0   1   2   3   4   5   6   7   8   9
 *
 *   0  *---*---*---*---*---*---*---*---*---*
 *      |   |   |   |   |   |   |   |   |   |
 *   1  *---O---O---O---O---O---O---O---O---*
 *      |   |   |   |   |   |   |   |   |   |
 *   2  *---O---O---O---O---O---O---O---O---*
 *      |   |   |   |   |   |   |   |   |   |
 *   3  *---O---O---O---O---O---O---O---O---*    "*" is boundary node
 *      |   |   |   |   |   |   |   |   |   |
 *   4  *---O---O---O---O---O---O---O---O---*    "O" is interior node
 *      |   |   |   |   |   |   |   |   |   |
 *   5  *---O---O---O---O---O---O---O---O---*
 *      |   |   |   |   |   |   |   |   |   |
 *   6  *---O---O---O---O---O---O---O---O---*
 *      |   |   |   |   |   |   |   |   |   |
 *   7  *---O---O---O---O---O---O---O---O---*
 *      |   |   |   |   |   |   |   |   |   |
 *   8  *---O---O---O---O---O---O---O---O---*
 *      |   |   |   |   |   |   |   |   |   |
 *   9  *---*---*---*---*---*---*---*---*---*
 *
 * Below are the four partitioned grid parts.  In each grid subdomain,
 * "O,A,B,C,D" represent grid points that are computed.  The points
 * labeled "A,B,C,D" are the grid points adjacent to an interior
 * boundary so they must be copied to added positions in adjacent
 * subdomains, indicated by "a,b,c,d" (the halo/ghost boundaries).
 *
 *      0   1   2   3   4  5   4  5   6   7   8   9
 *
 *   0  *---*---*---*---*--*   *--*---*---*---*---*
 *      |   |   |   |   |         |   |   |   |   |
 *   1  *---O---O---O---A--b   a--B---O---O---O---*
 *      |   |   |   |   |         |   |   |   |   |
 *   2  *---O---O---O---A--b   a--B---O---O---O---*
 *      |   |   |   |   |         |   |   |   |   |
 *   3  *---O---O---O---A--b   a--B---O---O---O---*   "*" is boundary node
 *      |   |   |   |   |         |   |   |   |   |
 *   4  *---A---A---A---A--b   a--B---B---B---B---*   "O" is interior node
 *      |   |   |   |   |         |   |   |   |   |
 *   5  *   d   d   d   d         c   c   c   c   *   "A,B,C,D" are computed
 *                                                    interior boundary nodes
 *   4  *   a   a   a   a         b   b   b   b   *
 *      |   |   |   |   |         |   |   |   |   |   "a,b,c,d" are copies of
 *   5  *---D---D---D---D--c   d--C---C---C---C---*   interior boundary nodes
 *      |   |   |   |   |         |   |   |   |   |
 *   6  *---O---O---O---D--c   d--C---O---O---O---*
 *      |   |   |   |   |         |   |   |   |   |
 *   7  *---O---O---O---D--c   d -C---O---O---O---*
 *      |   |   |   |   |         |   |   |   |   |
 *   8  *---O---O---O---D--c   d--C---O---O---O---*
 *      |   |   |   |   |         |   |   |   |   |
 *   9  *---*---*---*---*--*   *--*---*---*---*---*
 *
 */

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <mpi.h>
#include <hdf5.h>

#include "cartesian-mpi.cc"

#define IDX(i,j,stride) ((i)*(stride)+(j)) // row major (C/C++)

/*----------------------------------------------------------------------------
 * Print out 2D grid data
 *
 * Input:
 *   double* v       - two-dimensional array holding grid data
 *   int nx, ny      - dimensions of grid
 *
 * Output:
 *   None, other than output to stdout.
 */
void show_grid(double* v, int nx, int ny)
{
    printf("--------------------------------------------------------\n");
    for (int j = ny - 1; j >= 0; j--)
    {
        for (int i = 0; i < nx; i++) printf(" %6.4f", v[IDX(i,j,ny)]);
        printf("\n");
    }
    printf("--------------------------------------------------------\n");
    fflush(stdout);
}

/*----------------------------------------------------------------------------
 * Print out 2D grid data in rank-order
 *
 * Input:
 *   double* u                   - two-dimensional array holding grid data
 *   Cartesian_Block* halo_grid  - grid parameters
 *   int num_proc                - number of processes
 *   int rank                    - process rank
 *   MPI_Comm comm               - MPI Cartesian communicator
 *
 * Output:
 *   None, other than output to stdout.
 */
void dump_grid_rank_order(double* u, Cartesian_Block* halo_grid,
                          int num_proc, int rank, MPI_Comm comm)
{
    for (int i = 0; i < num_proc; i++)
    {
        if (rank == i)
        {
            printf("Rank: %d\n", rank);
            show_grid(u, halo_grid->nx, halo_grid->ny);
        }
        // Use barrier and delay so output comes out in process rank order
        MPI_Barrier(comm);
        usleep(10000); // 0.01 sec
    }
    
}

/*----------------------------------------------------------------------------
 * Exchange halo boundary data between 2D grid blocks.
 *
 * Input:
 *   double** u            - 2-D array holding single block of grid data
 *   Cartesian_Block* grid - pointer to grid block data
 *   MPI_Datatype x_slice  - type representing horizontal row of data
 *   MPI_Datatype y_slice  - type representing vertical column of data
 *   MPI_Comm comm         - communicator
 *
 * Output:
 *   double* u             - 2-D array holding grid block with update halo
 */
void exchange_halo_data(double* u, Cartesian_Block* grid,
                        MPI_Datatype x_slice, MPI_Datatype y_slice,
                        MPI_Comm comm)
{
    const int tag = 0;
    const int nx = grid->nx;
    const int ny = grid->ny;
    const int m = grid->ny;

    // Send top row of my data to bottom halo of neighbor above me and
    // receive bottom row of same neighbor's data into my top halo
    MPI_Sendrecv(&u[IDX(0,ny-2,m)], 1, x_slice, grid->above_neighbor, tag,
                 &u[IDX(0,0,m)],    1, x_slice, grid->below_neighbor, tag,
                 comm, MPI_STATUS_IGNORE);

    // Send bottom row of my data to top halo of neighbor below me and
    // receive top row of same neighbor's data into my bottom halo
    MPI_Sendrecv(&u[IDX(0,1,m)],    1, x_slice, grid->below_neighbor, tag,
                 &u[IDX(0,ny-1,m)], 1, x_slice, grid->above_neighbor, tag,
                 comm, MPI_STATUS_IGNORE);

    // Send right column of my data to left halo of neighbor to my right
    // and receive left column of same neighbor's data into my right halo
    MPI_Sendrecv(&u[IDX(nx-2,0,m)], 1, y_slice, grid->right_neighbor, tag,
                 &u[IDX(0,0,m)],    1, y_slice, grid->left_neighbor,  tag,
                 comm, MPI_STATUS_IGNORE);

    // Send left column of my data to right halo of neighbor to my left
    // and receive right column of same neighbor's data into my left halo
    MPI_Sendrecv(&u[IDX(1,0,m)],    1, y_slice, grid->left_neighbor,  tag,
                 &u[IDX(nx-1,0,m)], 1, y_slice, grid->right_neighbor, tag,
                 comm, MPI_STATUS_IGNORE);
}

/*----------------------------------------------------------------------------
 * Check return values from HDF5 routines
 */
#define CHKERR(status,name) if (status < 0) \
        fprintf(stderr, "Error: nonzero status (%d) in %s\n", status, name)

/*----------------------------------------------------------------------------
 * Create HDF5 file and store grid data (excluding halo) in it.
 *
 * Input:
 *   char* fname           - name of file to create
 *   char* dname           - name of dataset in file
 *   double** u            - local portion of array to write to file
 *   int NX, NY            - dimensions of overall array
 *   Cartesian_Block* orig_grid - grid block parameters (data only, no halo)
 *   Cartesian_BLock* halo_grid - grid block parameters (including halo data)
 *   MPI_Comm comm         - MPI communicator used by procs sharing in write
 *
 * Output:
 *   None, other than created file
 */
void writeFile(const char* fname, const char* dname, double* u,
               int NX, int NY, Cartesian_Block* orig_grid,
               Cartesian_Block* halo_grid, MPI_Comm comm)
{
    hid_t plist_id;       // property list
    hid_t file_id;        // output file
    hid_t dataspace_id;   // dataspace of output file
    hid_t dataset_id;     // dataset in output file
    hid_t memspace_id;    // dataspace of local memory
    hsize_t dimsf[2];     // dimensions of dataspace in file
    hsize_t dimsm[2];     // dimensions of dataspace in memory
    hsize_t count[2];     // number of blocks
    hsize_t offset[2];    // start values of data block locations
    herr_t status;

    // Create property list for output file using data from MPI communicator
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if (plist_id < 0) fprintf(stderr, "Could not create property list\n");
    status = H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);
    CHKERR(status, "H5Pset_fapl_mpio()");

    // Create the output file - property list not needed after this
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    if (file_id < 0) fprintf(stderr, "Could not output file\n");
    CHKERR(H5Pclose(plist_id), "H5Pclose()");

    // Create output file dataspace and dataset that will hold data from all
    // processes.  Use dimensions for entire grid.
    dimsf[0] = NX;
    dimsf[1] = NY;
    dataspace_id = H5Screate_simple(2, dimsf, NULL);
    if (dataspace_id < 0) fprintf(stderr, "Could not create dataspace\n");

    dataset_id = H5Dcreate(file_id, dname, H5T_IEEE_F64LE, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) fprintf(stderr, "Could not create dataset\n");

    // Create local dataspace corresponding to process's portion of the
    // global dataspace.  Use the actual dimensions for the portion of
    // the data this process is responsible for, including halo.
    dimsm[0] = halo_grid->nx;
    dimsm[1] = halo_grid->ny;
    memspace_id = H5Screate_simple(2, dimsm, NULL);
    if (memspace_id < 0) fprintf(stderr, "Could not create memspace\n");

    // Define the hyperslab in the global dataspace corresponding to the
    // data to be written by this process.
    offset[0] = orig_grid->x0;
    offset[1] = orig_grid->y0;
    count[0]  = orig_grid->nx;
    count[1]  = orig_grid->ny;
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL,
                                 count, NULL);
    CHKERR(status, "H5Sselect_hyperslab() dataspace");

    // Define the hyperslab in the local dataspace corresponding to the
    // data to be written by this process.  Offsets are nonzero
    // if there is a ghost/halo boundary to ignore on the left or bottom.
    offset[0] = orig_grid->x0 - halo_grid->x0;
    offset[1] = orig_grid->y0 - halo_grid->y0;
    status = H5Sselect_hyperslab(memspace_id, H5S_SELECT_SET, offset, NULL,
                                 count, NULL);
    CHKERR(status, "H5Sselect_hyperslab() memspace");

    // Set transfer mode
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if (plist_id < 0) fprintf(stderr, "Could not create property list\n");
    //status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    CHKERR(status, "H5Pset_dxpl_mpio()");

    // Finally - we can write the data!
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id,
                      dataspace_id, plist_id, u);
    CHKERR(status, "H5Dwrite()");

    // All done -- release all remaining open dataspaces, datasets, etc.
    CHKERR(H5Pclose(plist_id), "H5Pclose()");
    CHKERR(H5Sclose(memspace_id), "H5Sclose()");
    CHKERR(H5Dclose(dataset_id), "H5Dclose()");
    CHKERR(H5Sclose(dataspace_id), "H5Sclose()");
    CHKERR(H5Fclose(file_id), "H5Fclose()");
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    int NX = 10;                // number of grid points in x direction
    int NY = 10;                // number of grid points in y direction
    int num_proc;               // number of participating processes
    int rank;                   // process rank within communicator
    int dims[2] = {0, 0};       // allow MPI to choose grid block dimensions
    int periodic[2] = {0, 0};   // domain is non-periodic
    int may_rerank = 1;         // allow processes to be re-ranked
    Cartesian_Block halo_grid;  // parameters for grid block including halo
    Cartesian_Block orig_grid;  // parameters for grid block without halo
    double* u = NULL;           // array to hold grid block data
    MPI_Comm comm2d;            // Cartesian communicator
    MPI_Datatype x_slice;       // datatype for horizontal slice (single row)
    MPI_Datatype y_slice;       // datatype for vertical slice (single column)

    // Initialize MPI

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Process command line.  All arguments are optional but first and
    // second must occur as a pair, if at all, and third and fourth
    // must occur as a pair, if at all.

    if (argc == 2 || argc == 4 || argc > 5)
    {
        if (rank == 0)
        {
            printf("Usage: %s [NX NY [DIMX DIMY]]\n\n", argv[0]);
            printf("- NX and NY are the grid dimensions.\n");
            printf("- DIMX and DIMY are the process grid dimensions.\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Override default grid dimensions if requested.  Make sure both
    // values are positive.

    if (argc >= 3)
    {
        NX = atoi(argv[1]);
        NY = atoi(argv[2]);
    }
    if (NX <= 0 || NY <= 0)
    {
        if (rank == 0)
        {
            fprintf(stderr, "Error: both NX and NY must be positive.\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Override default process grid dimensions if requested.  Make
    // sure that if the user specified block grid dimensions they are
    // consistent with the number of processes.

    if (argc >= 5)
    {
        dims[0] = atoi(argv[3]);
        dims[1] = atoi(argv[4]);
        if (dims[0] * dims[1] != num_proc)
        {
            if (rank == 0)
            {
                fprintf(stderr, "Product of grid block dimensions must ");
                fprintf(stderr, "match the number of processes\n");
            }
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }

    // Set up Cartesian grid of processors.  A new communicator is
    // created and our rank may be adjusted relative to this new
    // communicator.

    comm2d = mpi_cart_setup(num_proc, NX, NY, may_rerank, &rank, dims,
                            periodic, &x_slice, &y_slice, &halo_grid,
                            &orig_grid);

    // Create my portion of the grid.  For the exchange to work
    // properly we must have a constant stride in each dimension.

    u = new double [halo_grid.nx * halo_grid.ny];

    // Since this is a demonstration program, here we initialize the
    // local portion of the grid with values that indicate their
    // original position in the grid.  Assuming nx and ny are each
    // less than 100 these values have the form R.XXYY where: R is the
    // rank of of the process that created the data, XX is the x
    // coordinate in the grid (0 is at left), and YY is the y
    // coordinate in the grid (0 is at bottom)

    for (int j = 0; j < halo_grid.ny; j++)
    {
        for (int i = 0; i < halo_grid.nx; i++)
        {
            u[IDX(i,j,halo_grid.ny)] = rank + 0.01 * (i + halo_grid.x0)
                + 0.0001 * (j + halo_grid.y0);
        }
    }

    // Ready to start exchanging data.  Wait for my turn and then
    // display my portion of the grid.

    if (rank == 0)
    {
        printf("\n");
        printf("********************************************************\n");
        printf("**** Local grids BEFORE exchange ***********************\n");
        printf("********************************************************\n");
    }
    dump_grid_rank_order(u, &halo_grid, num_proc, rank, comm2d);

    /*******************************************************************
     * Now we're done setting things up.  In a "real" program, now     *
     * we'd begin doing whatever we need to do with the data.  In this *
     * case, we merely exchange the slices to copy halo boundary data  *
     * to adjacent processes.                                          *
     *******************************************************************/

    exchange_halo_data(u, &halo_grid, x_slice, y_slice, comm2d);

    /*******************************************************************
     * Here we can do work that can be done once the exchange is done. *
     *******************************************************************/

    // Exchange cycle is complete.  Wait for my turn and then display
    // my portion of the grid.

    if (rank == 0)
    {
        printf("\n");
        printf("********************************************************\n");
        printf("**** Local grids AFTER exchange ************************\n");
        printf("********************************************************\n");
    }
    dump_grid_rank_order(u, &halo_grid, num_proc, rank, comm2d);

    // Write grid data to HDF5 file

    writeFile("cart-hdf5.dat", "/grid", u, NX, NY, &orig_grid, &halo_grid,
              comm2d);

    // Release memory and datatypes and then quit

    delete [] u;
    MPI_Type_free(&x_slice);
    MPI_Type_free(&y_slice);

    MPI_Finalize();
    return 0;
}
