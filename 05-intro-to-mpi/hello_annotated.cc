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
    int tag = 42;       // message tag (answer to the ultimate question)
    MPI_Status status;  // message receive status

    // Initalize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(
        MPI_COMM_WORLD, // communicator
        &my_rank        // address of location to store the rank for this
                        //   process out of all processes connected through
                        //   MPI_COMM_WORLD
        );
    MPI_Comm_size(
        MPI_COMM_WORLD, // communicator
        &num_proc       // address of location to store number of processes
                        //   that are connected through MPI_COMM_WORLD
        );

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
	    MPI_Send(
                &msg,          // address of start of data to send
                1,             // number of elements to send
                MPI_INT,       // type of elements being sent
                i,             // rank of process to send data to
                tag,           // tag for messsage (must match receiving tag)
                MPI_COMM_WORLD // communicator we are connected to/through
                );
	}
    }
    else
    {
        // Processes with rank greater than 0 wait to receive a message
        // from the rank 0 process and then display it.
	MPI_Recv(
            &msg,              // address to store received data
            1,                 // number of elements to receive
            MPI_INT,           // type of elements being received
            0,                 // rank of process to receive data from
            tag,               // tag for message (must match sending tag)
            MPI_COMM_WORLD,    // communictor we are connected to/through
            &status            // address to store status block
            );
        std::cout << "Process " << std::setw(2) << my_rank
                  << " received message " << std::setw(4) << msg << std::endl;
    }

    // Clean up MPI
    MPI_Finalize();

    return 0;
}
