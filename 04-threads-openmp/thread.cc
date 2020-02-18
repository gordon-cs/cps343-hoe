/*
 * Jonathan Senning <senning@gordon.edu>
 * Gordon College
 * Written: February  4, 2002
 * Revised: February 20, 2016
 * 
 * $Smake: g++ -O2 -Wall -o %F %f -lpthread
 *
 * This program demonstrates a very simple use of pthreads in a C++ program.
 * Two simple classes are created; one for counting semaphores and the other
 * for simple threads.  Pthreads are used for threads, and the class defined
 * for threads here is not sophisticated at all.
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
#include <pthread.h>

using namespace std;

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
    sem_t mySemaphore;

public:
    Semaphore( int value = 1 ) { (void) sem_init( &mySemaphore, 0, value ); };
   ~Semaphore( void )          { (void) sem_destroy( &mySemaphore ); };
    void wait( void )          { (void) sem_wait( &mySemaphore ); };
    void post( void )          { (void) sem_post( &mySemaphore ); };
};

//============================================================================
// Thread Class
//============================================================================

// Define a class for simple threads

// In this class all the constructor does is record the run method for the
// thread being created.  The thread is not actually created and started
// until the start() method is invoked.  The single parameter to the start
// method is passed directly to the run method for the thread, and should
// be a void pointer.  The thread continues to exist until the run method
// terminates.

class Thread
{
private:
    pthread_t threadID;             // thread identifier
    void* (*myRunMethod)(void *);   // pointer to thread's run method

public:
    Thread( void* (*run)(void *) ) : myRunMethod( run ) {};
   ~Thread( void ) {};
    void start( void* arg );
    void join( void ) { pthread_join( threadID, NULL ); };
};

void Thread::start( void* arg )
// Start the thread
{
    if ( pthread_create( &threadID, NULL, myRunMethod, arg ) != 0 )
    {
        cerr << "pthread_create: could not create thread" << endl;
        exit( EXIT_FAILURE );
    }
}

//============================================================================
// Global (Shared) Variables
//============================================================================

// Define some global variables.  These must be global because they are shared
// by both the parent and the child thread.

Semaphore s1( 0 );      // Semaphore will block on first call to wait()
Semaphore s2( 0 );      // ditto

int counter = 0;        // Initalize shared counter

//============================================================================
// Thread run method (body of thread)
//============================================================================

void my_run( int* arg )
// Thread run method.  After starting this waits for a signal from the parent
// indicating that the shared variable "counter" should be incremented.  The
// increment is completed and a signal is sent back to the parent indicating
// that fact.  The main loop is executed until the counter equals or exceeds
// the integer pointed to by the argument to this method.
{
    cout << "child  thread: running" << endl;
    cout << "child  thread: argument is " << *arg << endl;
    cout << "child  thread: counter = " << counter << endl;

    while ( counter < *arg )
    {
        s1.wait();      // wait for signal to increment
        counter++;
        s2.post();      // signal that increment is complete
    }
    
    cout << "child  thread: done" << endl;
}

//============================================================================
// Main program
//============================================================================

int main( void )
{
    int numberOfIncrements = 10;   // Desired number of steps

    // Instantiate the child thread
    Thread child( (void* (*)(void *)) &my_run );
    
    cout << "parent thread: running" << endl;

    // Start the child thread running
    child.start( &numberOfIncrements );

    // Main loop
    while ( counter < numberOfIncrements )
    {
        sleep( 1 );
        s1.post();      // signal we are ready for increment
        s2.wait();      // wait for signal that increment is complete
        cout << "parent thread: counter = " << counter << endl;
    }

    // Wait for child thread to finish
    child.join();

    return 0;
}
