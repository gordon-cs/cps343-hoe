// $Smake: g++ -std=c++11 -Wall -O2 -pthread -o %F %f

/*
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Gordon College
 * February 4, 2002
 * February 15, 2018 -- updated to use C++11 threads
 *
 * This program demonstrates a very simple use of C++11 threads.
 * A class for counting semaphores based on POSIX semaphores is created.
 *
 * This program starts a single child thread that increments a global variable.
 * Two semaphores are used for synchronization - one to tell the child thread
 * that an increment is desired and the other to tell the parent that the
 * the increment has completed.
 */

#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <semaphore.h>
#include <thread>

//============================================================================
// Semaphore Class
//============================================================================

// Define a class based on POSIX semaphore functions

// The semaphores defined here are counting semaphores.  By default the
// semaphores are appropriate for a MUTEX - the first time the wait function
// is called the caller is allowed to pass.  Instantiating a semaphore with
// a value other than 1 changes this behavior; in particular if the value is
// 0 then the first caller of wait will block until a signal arrives (unless
// a signal has already arrived, in which case it will pass immediately).

class Semaphore
{
  private:
    sem_t _semaphore;

  public:
    Semaphore( int value = 1 )	{ (void) sem_init( &_semaphore, 0, value ); };
   ~Semaphore( void )		{ (void) sem_destroy( &_semaphore ); };
    void wait( void )		{ (void) sem_wait( &_semaphore ); };
    void post( void )		{ (void) sem_post( &_semaphore ); };
};

//============================================================================
// Global (Shared) Variables
//============================================================================

// Define some global variables.  These must be global because they are shared
// by both the parent and the child thread.

Semaphore s1( 0 );	// Semaphore will block on first call to wait()
Semaphore s2( 0 );	// ditto
int counter = 0;	// Initalize shared counter

//============================================================================
// Thread run method (body of thread)
//============================================================================

void my_run( int n )
// Thread run method.  After starting this waits for a signal from the parent
// indicating that the shared variable "counter" should be incremented.  The
// increment is completed and a signal is sent back to the parent indicating
// that fact.  The main loop is executed until the counter equals or exceeds
// the integer pointed to by the argument to this method.
{
    std::cout << "child  thread: running" << std::endl;
    std::cout << "child  thread: argument is " << n << std::endl;
    std::cout << "child  thread: counter = " << counter << std::endl;

    while ( counter < n )
    {
	s1.wait();	// wait for signal to increment
	counter++;
	s2.post();	// signal that increment is complete
    }

    std::cout << "child  thread: done" << std::endl;
}

//============================================================================
// Main program
//============================================================================

int main( void )
{
    int numberOfIncrements = 10;	// Desired number of steps

    std::cout << "parent thread: running" << std::endl;

    // Instantiate and start the child thread
    std::thread child( my_run, numberOfIncrements );
    
    // Main loop
    while ( counter < numberOfIncrements )
    {
	sleep( 1 );
	s1.post();	// signal we are ready for increment
	s2.wait();	// wait for signal that increment is complete
	std::cout << "parent thread: counter = " << counter << std::endl;
    }

    // Wait for child thread to finish
    child.join();

    return 0;
}
