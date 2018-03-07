// $Smake: mpic++ -Wall -O2 -o %F %f
// Demonstrate basic MPI set up and simple message passing

#include <iostream>
#include <mpi.h>

int main( int argc, char* argv[] )
{
    using namespace std;

    int pRank;          // process rank
    int nProc;          // number of processes
    int msg;            // message buffer
    int tag = 99;       // message tag
    MPI_Status status;  // message receive status

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &pRank );
    MPI_Comm_size( MPI_COMM_WORLD, &nProc );

    if ( pRank == 0 )
    {
	cout << "There are " << nProc << " processes " << endl;
	// send a message to each non-rank 0 process
	for ( int i = 1; i < nProc; i++ )
	{
	    msg = 100 * i;
	    MPI_Send( &msg, 1, MPI_INT, i, tag, MPI_COMM_WORLD );
	}
    }
    else
    {
	// receive message from rank 0 process
	MPI_Recv( &msg, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status );
	cout << "Rank " << pRank << " process received " << msg << endl;
    }

    MPI_Finalize();

    return 0;
}
