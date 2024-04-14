/*
 * Sequential version
 *
 * $Smake: g++ -Wall -O2 -o %F %f wtime.c -lgsl -lgslcblas
 *
 * Usage:
 *    pi_sequential [-n NUM_SAMPLES] [-q]
 *
 * where the [] indicate an optional parameter.  Default values for the 
 * optional parameters are provided.
 *
 * This program is released into the public domain.
 *
 * ==========================================================================
 *
 * Estimate Pi using Monte Carlo sampling.  This is done using a 1x1 square
 * and a quarter circle of radius 1.  
 *
 *                    1 ----------------------------
 *                      |* * * *                   |
 *                      |        * *               |
 *                      |            * *           |
 *                      |                 *        |
 *                      |                   *      |
 *                      |                     *    |
 *                      |                     *    |
 *                      |                       *  |
 *                      |                       *  |
 *                      |                         *|
 *                      |                         *|
 *                      |                         *|
 *                      |                         *|
 *                    0 ----------------------------
 *                      0                          1
 *
 * The area of the square is 1 and the area the quarter circle is Pi/4.  The
 * ratio of these values should be the same as the ratio of samples inside
 * the quarter circle to the total number of samples, so
 *
 *             samples inside quarter circle     Pi/4
 *             ----------------------------- =  ----  = Pi/4
 *                total number of samples         1
 *
 * Thus, Pi can be estimated as
 *
 *                     4 (samples inside quarter circle)
 *               Pi ~= ---------------------------------
 *                        (total numeber of samples)
 *
 * ==========================================================================
 */

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <unistd.h>
#include "wtime.h"
#include <gsl/gsl_rng.h>

//----------------------------------------------------------------------------
// Generate n random (x,y) samples in a 1x1 box and return the number that
// fall inside a quarter circle of radius 1.
//
// Input:   long n               - number of samples
//          gsl_rng* rng         - pointer to rng instance
//
// Returns: long                 - number of samples inside quarter circle

long genSamples(const long n, const gsl_rng* rng)
{
    long count = 0L;

    for (long i = 0L; i < n; i++)
    {
        const double x = gsl_rng_uniform(rng);
        const double y = gsl_rng_uniform(rng);
        count += (x * x + y * y < 1.0 ? 1 : 0);
    }

    return count;
}

//----------------------------------------------------------------------------
// Set up pseudo-random number generator (PRNG), generate n points in unit
// square, and return count of inside unit circle.
//
// Input:   long numSamples  - number of samples
//
// Returns: long             - number of samples inside quarter circle

long genRandomSamples(const long numSamples)
{
    // Set up random number generator using PID and system time as seed
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 100 * getpid() + time(NULL));

    // Generate the samples and count those inside the circle
    long count = genSamples(numSamples, rng);

    // Clean up
    gsl_rng_free(rng);
    return count;
}

//----------------------------------------------------------------------------
// Display results with or without labels
//
// Input:   long numSamples  - number of samples
//          double estimate  - estimate of Pi
//          double wtime     - elapsed wall-clock time
//          bool noLabels    - don't show labels if true
//
// Returns: nothing

void displayResults(long numSamples, double estimate, double wtime,
                    bool noLabels)
{
    const char* format = 
        (noLabels ?
         "%12.10f %10.3e %10.6f %ld\n" :
         "Pi: %12.10f, error: %10.3e, seconds: %g, samples: %ld\n");
    printf(format, estimate, M_PI - estimate, wtime, numSamples);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    long numSamples = 50000000L;
    bool quiet = false;

    // Process command line
    int c;
    while ((c = getopt(argc, argv, "n:q")) != -1)
    {
        switch(c)
        {
            case 'n':
                numSamples = atol(optarg);
                if (numSamples <= 0)
                {
                    fprintf(stderr, "number of samples must be positive\n");
                    fprintf(stderr, "got: %ld\n", numSamples);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'q':
                quiet = true;
                break;
            default:
                fprintf(stderr, "usage: %s [-n NUM_SAMPLES] [-q]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // Get samples and compute estimate of Pi
    double t1 = wtime();
    long count = genRandomSamples(numSamples);
    double result = 4.0 * count / numSamples;
    double t2 = wtime();

    // Display result
    displayResults(numSamples, result, t2 - t1, quiet);

    // All done
    return 0;
}
