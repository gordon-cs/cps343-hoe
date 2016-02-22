/*
 * Jonathan Senning <senning@gordon.edu>
 * Gordon College
 * Written: Spring 1998
 * Revised: February 20, 2016
 *
 * $Smake: g++ -Wall -O -o %F %f -lpthread
 *
 * Simulate operation of a bank.  The bank has a single waiting line served
 * by a number of tellers.  Customers arrive with a exponentially-distributed
 * random arrival time and each requires a exponentially-distributed random
 * service time.  POSIX semaphores are used for synchronization.
 */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <list>
#include <unistd.h>
#include <pthread.h>
#include <semaphore.h>
#include "randomInt.h"

using namespace std;

const unsigned int clockTick = 10000; // microseconds

//============================================================================
// Classes
//============================================================================

// Semaphore class based on POSIX semaphore functions

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

// Customer class

class Customer
{
private:
    pthread_t threadID;                // thread identifier
    void*     (*myRunMethod)(void *);  // pointer to thread's run method
    int       myID;                    // customer identifier
    int       myServiceTime;           // time to serivce customer

public:
    Customer( void* (*)(void *) = NULL, int custID = -1, int serviceTime = 0 );
   ~Customer( void ) {};
    void start( void );
    void wait( void ) { pthread_join( threadID, NULL ); };
    int getNumber( void ) { return myID; };
    int getServiceTime( void ) { return myServiceTime; };
};

Customer::Customer( void* (*run)(void *), int custID, int serviceTime )
    : myRunMethod( run ), myID( custID ), myServiceTime( serviceTime ) {}

void Customer::start( void )
// Start the thread representing the customer
{
    if ( ! myRunMethod )
    {
        cerr << "Customer::start(): Customer not initialized" << endl;
        exit( EXIT_FAILURE );
    }
    if ( pthread_create( &threadID, NULL, myRunMethod, this ) != 0 )
    {
        cerr << "could not create thread for Customer" << endl;
        exit( EXIT_FAILURE );
    }
}

//============================================================================
// Global (shared) variables -- shared between all threads
//============================================================================

int numberCustomersServed = 0;
int customerTimeInLine = 0;

Semaphore tellersFree;  // Counting semaphore for number of free tellers
Semaphore outputLock;   // mutex so only one thread at a time is doing output
Semaphore updateLock;   // mutex for shared variable updates

//============================================================================
// Functions
//============================================================================

int bankClock( void )
// Returns the number of clock ticks (virtual seconds) since the first call
// to this function.
//
// THIS IS NOT THREAD SAFE - use of static variable here assumes a single
// thread calls this function once before any chance of multiple threads
// calling it.
{
    static int startTime = 0;

    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    int t = int( 1.0e6 * ( ts.tv_sec + ts.tv_nsec / 1.0e9 ) / clockTick );

    if ( startTime == 0 ) startTime = t;
    return t - startTime;
}

//============================================================================

void* CustomerRun( Customer* customer )
// Simulates a customer entering bank, waiting in line, proceeding to a
// teller window, doing business, and exiting the bank.
{
    // Enter bank

    int enterBankTime = bankClock();

    outputLock.wait();
    cout << "Time " << setw( 3 ) << enterBankTime
         << ": customer " << setw( 2 ) << customer->getNumber()
         << " arrives in line" << endl << flush;
    outputLock.post();

    // Wait for a teller to be free

    tellersFree.wait();

    // Move to free teller

    int beginServiceTime = bankClock();
    int waitTime = beginServiceTime - enterBankTime;

    outputLock.wait();
    cout << "Time " << setw( 3 ) << beginServiceTime
         << ": customer " << setw( 2 ) << customer->getNumber()
         << "  starts being served";
    if ( waitTime > 0 ) cout << "       (wait time was " << waitTime << ")";
    cout << endl << flush;
    outputLock.post();

    // Simulate customer interacting with teller

    usleep( customer->getServiceTime() * clockTick );

    // Done with teller, exit bank

    int endServiceTime = bankClock();

    outputLock.wait();
    cout << "Time " << setw( 3 ) << endServiceTime
         << ": customer " << setw( 2 ) << customer->getNumber()
         << "   finishes         "
         << "     (service time was " << endServiceTime - beginServiceTime
         << ")" << endl << flush;
    outputLock.post();

    tellersFree.post();  // Teller service this customer is now free

    // Update global statistics

    updateLock.wait();
    numberCustomersServed++;
    customerTimeInLine += waitTime;
    updateLock.post();

    // terminate thread

    pthread_exit( NULL );
    return NULL;
}

//============================================================================
// Main program
//============================================================================

int main( int argc, char *argv[] )
{
    // Read command line arguments

    if ( argc != 5 )
    {
        cerr << "usage:" << endl;
        cerr << argv[0] << " <# of tellers> <arrival spacing>";
        cerr << " <service time> <simulation length>" << endl;
        exit( EXIT_FAILURE );
    }

    int numberOfTellers         = atoi( argv[1] );
    int meanTimeBetweenArrivals = atoi( argv[2] );
    int meanServiceTime         = atoi( argv[3] );
    int simulationDuration      = atoi( argv[4] );

    cout << "Number of tellers = " << numberOfTellers << endl;
    cout << "Arrival spacing   = " << meanTimeBetweenArrivals << endl;
    cout << "Service time      = " << meanServiceTime << endl;
    cout << "simulation length = " << simulationDuration << endl;

    // Create list for customers and initialize the waiting line (teller)
    // semaphore so that the first "numberOfTellers" customers pass
    // through immediately

    list<Customer*> customers;
    tellersFree = Semaphore( numberOfTellers );

    // Run the simulation

    outputLock.wait();
    cout << "=======================================================" << endl;
    cout << "Time " << setw( 3 ) << bankClock() << ": Bank is open"   << endl;
    cout << "=======================================================" << endl;
    cout << flush;
    outputLock.post();

    int customerNumber = 1;
    int timeToNextArrival = getRandomInteger( meanTimeBetweenArrivals );
    while ( bankClock() < simulationDuration )
    {
        while ( timeToNextArrival == 0 )
        {
            // A new customer arrives at the bank
            int serviceTime = getRandomInteger( meanServiceTime );
            Customer* newCustomer = new Customer( 
                (void* (*)(void *)) &CustomerRun, customerNumber++, serviceTime
                );
            customers.push_back( newCustomer );
            newCustomer->start();

            timeToNextArrival = getRandomInteger( meanTimeBetweenArrivals );
        }
        timeToNextArrival--;
        usleep( clockTick );
    }

    outputLock.wait();
    cout << "=======================================================" << endl;
    cout << "Time " << setw( 3 ) << bankClock() << ": Bank is closed" << endl;
    cout << "=======================================================" << endl;
    cout << flush;
    outputLock.post();

    // Wait for all customers to finish

    while ( !customers.empty() )
    {
        customers.front()->wait();
        delete customers.front();
        customers.pop_front();
    }

    outputLock.wait();
    cout << "=======================================================" << endl;
    cout << "All customers have left the bank" << endl;
    cout << "=======================================================" << endl;
    cout << flush;
    outputLock.post();

    // Produce the final report

    cout << endl;
    cout << "Number of customers served:   "
         << setw( 6 ) << numberCustomersServed << endl;
    cout << "Total time in line:           "
         << setw( 6 ) << customerTimeInLine << endl;
    cout << "Average time waiting in line: "
         << setw( 6 )  << setprecision( 4 )
         << double( customerTimeInLine ) / double( numberCustomersServed )
         << endl;

    return 0;
}
