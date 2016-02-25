/*
 * Returns a random integer, exponentially distributed about mean with range
 * of values truncated to lie between 1 and 5 * mean.  The first time this
 * routine is called it initialized the random number generate with 
 * 
 * See http://www.lter.umn.edu/tools/utility/rexpon.d for a justification of
 * the formula used here.
 * 
 * Written by R. Bjork
 * Modified and Extended by J. Senning 2/22/98
 */

#include <cstdlib>
#include <cmath>
#include <ctime>
#include "randomInt.h"

//----------------------------------------------------------------------------

// Static variable to record the number of times the generator has been
// initialized.

static bool initialized = false;

//----------------------------------------------------------------------------

void initializeRandomInteger( int seed )
// Initializes the random number generator.  If seed is negative then the
// system clock is used to initialize the generator.
{
    if ( seed < 0 )
    {
        seed = (unsigned int) time( (time_t *) 0 );
    }
    srandom( seed );
    initialized = true;
}

//----------------------------------------------------------------------------

int getRandomInteger( int mean )
// Computes a random integer from an exponetial distribution with a
// specified mean.  If this routine is called and the generator has not
// yet been initialized, it initializes it using the system clock.
{    
    if ( !initialized )
    {
        initializeRandomInteger( -1 );
    }

    // Random number from uniform distibution on [0,1)

    double u = (double) ( ( random() & RAND_MAX ) / ( 1.0 + RAND_MAX ) );

    // Find number from exponentially distribution with specified mean
    // and round to an integer.

    int value = (int) ( 0.5 - mean * log( u ) );

    return value;
}
