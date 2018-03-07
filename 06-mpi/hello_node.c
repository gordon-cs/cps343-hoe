/*
 * $Smake: mpicc -Wall -g -O2 -o %F %f
 */

#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

#define STR_LENGTH 256

int main ( int argc, char *argv[] )
{
    int rank;
    int nproc;
    char nodeName[STR_LENGTH];

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    gethostname( nodeName, STR_LENGTH );
    printf( "hello from %-12s: process %d of %d\n", nodeName, rank, nproc );

    MPI_Finalize();

    return 0;
}
