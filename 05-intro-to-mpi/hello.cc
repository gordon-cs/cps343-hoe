// $Smake: mpic++ -Wall -O2 -o %F %f
// Demonstrate basic MPI set up and simple message passing

#include <cstdio>
#include <mpi.h>

int main(int argc, char* argv[])
{
    using namespace std;

    int my_rank;        // process rank
    int num_proc;       // number of processes
    int msg;            // message buffer
    int tag = 99;       // message tag
    MPI_Status status;  // message receive status

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    if (my_rank == 0)
    {
        printf("Hello from Rank %2d: There are %d processes\n",
                my_rank, num_proc);
	// send a message to each non-rank 0 process
	for (int i = 1; i < num_proc; i++)
	{
	    msg = 100 * i;
	    MPI_Send(&msg, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
	}
    }
    else
    {
	// receive message from rank 0 process
	MPI_Recv(&msg, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        printf("Hello from Rank %2d: Received %d\n", my_rank, msg);
    }

    MPI_Finalize();

    return 0;
}
