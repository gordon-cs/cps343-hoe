// $Smake: mpic++ -Wall -O2 -o %F %f

#include <cstdio>
#include <mpi.h>

int main( int argc, char* argv[] )
{
    int my_rank;
    int num_proc;
    int msg;
    const int tag = 42; // the answer to the ultimate question
    MPI_Status status;

    // Initalize MPI
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &num_proc );

    // Set who I'll be receive from (if anyone) and sending to (if anyone)
    const int prev = my_rank - 1;
    const int next = my_rank + 1;

    if ( my_rank == 0 )
    {
	// Rank 0 process is the master; it starts the token by sending it
	// to the rank 1 process.
	msg = 1000;
        printf( "Process %d sending %d to process %d\n", my_rank, msg, next );
	MPI_Send( &msg, 1, MPI_INT, next, tag, MPI_COMM_WORLD );
    }
    else
    {
	// We're not the rank 0 process so we wait for the token to arrive
	// from our predecessor, increment it, and send it along to our
	// successor, if any.
	MPI_Recv( &msg, 1, MPI_INT, prev, tag, MPI_COMM_WORLD, &status );
        printf( "Process %d received %d\n", my_rank, msg );

        // Increment message so we now it's been here
	msg++;

        // If we're not the last process, send the message on...
        if ( my_rank < num_proc - 1 )
        {
            printf( "Process %d sending %d to process %d\n",
                    my_rank, msg, next );
            MPI_Send( &msg, 1, MPI_INT, next, tag, MPI_COMM_WORLD );
        }
    }

    // All done, time to clean up
    MPI_Finalize();

    return 0;
}
