/*
 * Jonathan Senning <senning@gordon.edu>
 * Gordon College
 * Written: February 26, 2002
 * Revised: February 20, 2016
 *
 * $Smake: g++ -Wall -O -o %F %f -lpthread
 *
 * This program implements a solution of the Producer-Consumer problem using
 * a bounded buffer.  Semaphores are used to to control access to the buffer
 * and to ensure that an empty buffer is not read from and a full buffer is
 * not written to.
 */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <semaphore.h>
#include <pthread.h>

using namespace std;

const unsigned int max_produce_time = 100000;
const unsigned int max_consume_time = 150000;

//============================================================================
// Program Classes and supporting methods
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

// Define a class for a bounded Buffer

// This class implements a bounded buffer using a circular array of integers.

class BoundedBuffer
{
private:
    int *myBuffer;
    int  mySize;
    int  myHead, myTail;
    bool myIsFull;

public:
    BoundedBuffer( int n )
        : mySize( n ), myHead( 0 ), myTail( 0 ), myIsFull( false )
        { myBuffer = new int [mySize]; };
   ~BoundedBuffer( void ) { delete [] myBuffer; };
    void insert( int x );
    int  remove( void );
};

void BoundedBuffer::insert( int x )
// Insert x into the buffer
{
    if ( !myIsFull )
    {
        myBuffer[myTail] = x;
        myTail = ( myTail + 1 ) % mySize;
    }
    else
    {
        cerr << "**** Attempt to insert into a full buffer ****" << endl;
    }
    myIsFull = ( myHead == myTail );
}

int BoundedBuffer::remove( void )
// Return next value from the buffer
{
    int x = -1;
    if ( myIsFull || myHead != myTail )
    {
        x = myBuffer[myHead];
        myHead = ( myHead + 1 ) % mySize;
    }
    else
    {
        cerr << "**** Attempt to read from an empty buffer ****" << endl;
    }
    myIsFull = false;
    return x;
}

//============================================================================
// Global Variables
//============================================================================

const int bufsize = 4;           // Size of bounded buffer

BoundedBuffer buffer( bufsize ); // The shared buffer

Semaphore buflock;               // Mutual exclusion for buffer access
Semaphore output;                // Mutual exclusion for output access
Semaphore notfull( bufsize );    // Number of empty slots in buffer
Semaphore notempty( 0 );         // Number of full slots in buffer

//============================================================================
// Program Methods
//============================================================================

void producerRun( int* arg )
// Run method for the producer thread -- do what the producer does...
{
    int numberOfItems = *arg;
    int data = 0;

    struct random_data prng;
    char state[32];
    initstate_r( 0, state, sizeof( state ), &prng );
    srandom_r( (unsigned int) time( NULL ), &prng );

    while ( numberOfItems-- > 0 )
    {
        int rnd;
        random_r( &prng, &rnd );
        usleep( rnd % max_produce_time );  // simulate data production
        data++;

        output.wait();
        cout << "<-- producer: sending " << data << endl << flush;
        output.post();

        notfull.wait();         // wait until buffer has room

        buflock.wait();         // enter critical section
        buffer.insert( data );  // store data in buffer
        notempty.post();        // signal buffer has data
        buflock.post();         // leave critical section
    }
}

//============================================================================

void consumerRun( int* arg )
// Run method for consumer thread -- do what the consumer does...
{
    int numberOfItems = *arg;
    int data;
    int rnd;

    struct random_data prng;
    char state[32];
    initstate_r( 0, state, sizeof( state ), &prng );
    srandom_r( (unsigned int) time( NULL ) + 100000, &prng );

    while ( numberOfItems-- > 0 )
    {
        notempty.wait();        // wait until buffer has some data

        buflock.wait();         // enter critcal section
        data = buffer.remove(); // get data from buffer 
        notfull.post();         // signal buffer has at least one space
        buflock.post();         // leave critical section

        output.wait();
        cout << "--> consumer: received " << data << endl << flush;
        output.post();

        random_r( &prng, &rnd );
        usleep( rnd % max_consume_time ); // simulate using data
    }
}

//============================================================================

int main( int argc, char* argv[] )
{
    int numberOfItems = 10;
    if ( argc > 1 ) numberOfItems = atoi( argv[1] );

    // Create threads for producer and consumer

    Thread producer( (void* (*)(void *)) &producerRun );
    Thread consumer( (void* (*)(void *)) &consumerRun );

    // Start the threads

    producer.start( &numberOfItems );
    consumer.start( &numberOfItems );

    // Wait for threads to finish

    producer.join();
    consumer.join();

    return 0;
}
