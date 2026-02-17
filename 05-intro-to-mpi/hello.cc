// $Smake: mpic++ -Wall -O3 -o %F %f
// Demonstrate basic MPI set up and simple message passing

#include <iostream>
#include <iomanip>
#include <mpi.h>

int main(int argc, char* argv[])
{
    int my_rank;        // process rank
    int num_proc;       // number of processes
    int msg;            // message buffer
    int tag = 42;       // message tag (can really be any positive integer)
    MPI_Status status;  // message receive status

    // Initalize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

    if (my_rank == 0)
    {
        // Rank 0 process announces itself and then sends a message to
        // all other processes.  The messages are just 100 times the
        // destination process number.
        std::cout << "Hello from process " << my_rank
                  << " of " << num_proc << std::endl;
	for (int i = 1; i < num_proc; i++)
	{
	    msg = 100 * i;
	    MPI_Send(&msg, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
	}
    }
    else
    {
        // Processes with rank greater than 0 wait to receive a message
        // from the rank 0 process and then display it.
	MPI_Recv(&msg, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        std::cout << "Process " << std::setw(2) << my_rank
                  << " received message " << std::setw(4) << msg << std::endl;
    }

    // Clean up MPI
    MPI_Finalize();

    return 0;
}
