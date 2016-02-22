#ifndef __RANDOM_INT_H__
#define __RANDOM_INT_H__

//----------------------------------------------------------------------------

// Initializes the random number generator.  If seed is negative then the
// system clock is used to initialize the generator.

void initializeRandomInteger( int seed );

//----------------------------------------------------------------------------

// Computes a random integer from an exponetial distribution with a
// specified mean.  If this routine is called and the generator has not
// yet been initialized, it initializes it using the system clock.

int getRandomInteger( int mean );

//----------------------------------------------------------------------------

#endif /* __RANDOM_INT_H__ */
