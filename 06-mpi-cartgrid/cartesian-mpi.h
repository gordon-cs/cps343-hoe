#ifndef CARTESIAN_MPI_H
#define CARTESIAN_MPI_H

#include <mpi.h>

// Define structure to hold Cartesian communicator parameters.  Each process
// has its own grid that corresponds to a portion of the domain.  These
// parameters apply locally and not to the entire domain.

typedef struct Cartesian_Block
{
    int x0, x1;         // extreme left and right indices of grid region
    int y0, y1;         // extreme bottom and top indices of grid region
    int nx, ny;         // grid extent in horizontal and vertical directions
    int above_neighbor; // rank of process handling block above
    int below_neighbor; // rank of process handling block below
    int left_neighbor;  // rank of process handling block to left
    int right_neighbor; // rank of process handling block to right
} Cartesian_Block;

// Main function prototype

MPI_Comm mpi_cart_setup(int num_proc, int NX, int NY, int may_rerank,
                        int* rank, int* dims, int* periodic,
                        MPI_Datatype* x_slice, MPI_Datatype* y_slice,
                        Cartesian_Block* halo_grid,
                        Cartesian_Block* orig_grid);
#endif
