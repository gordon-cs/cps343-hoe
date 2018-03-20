/*
 * Returns the number of seconds since some fixed arbitrary time in
 * the past.  Resolution varies depending on what functions are
 * available during compilation.  If POSIX timers are available timer
 * should have an accuracy of about 1 nanosecond, otherwise the
 * resolution will only be about 10 milliseconds.
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 */

#include <unistd.h>
#include <time.h>
#include "wtime.h"

#if defined(_POSIX_TIMERS)

/*---- POSIX timer version --------------------------------------------------*/

double wtime( void )
{
    struct timespec ts;

    // prefer monotonic clock if it's available
# if defined(_POSIX_MONOTONIC_CLOCK)
    clock_gettime( CLOCK_MONOTONIC, &ts );
# else
    clock_gettime( CLOCK_REALTIME, &ts );
# endif

    return (double) ( ts.tv_sec + ts.tv_nsec / 1.0e9 );
}

#else  // defined(_POSIX_TIMERS)

/*---- Non-POSIX timer version ----------------------------------------------*/

#include <sys/times.h>
double wtime( void )
{
    struct tms buf;
    return double( times( &buf ) ) / double( sysconf( _SC_CLK_TCK ) );
}

#endif // defined(_POSIX_TIMERS)
