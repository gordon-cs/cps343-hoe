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
#include "wtime.h"

#if defined(__MACH__)
# include <mach/clock.h>
# include <mach/mach.h>
#elif defined(_POSIX_TIMERS)
# include <time.h>
#else
# include <sys/time.h>
#endif

double wtime( void )
{
#if defined(__MACH__)
    /* Apple OS X */
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service( mach_host_self(), CALENDAR_CLOCK, &cclock );
    clock_get_time( cclock, &mts );
    mach_port_deallocate( mach_task_self(), cclock );
    return (double) ( mts.tv_sec + mts.tv_nsec / 1.0e9 );
#elif defined(_POSIX_TIMERS)
    /* Linux */
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts ); /* could be CLOCK_REALTIME */
    return (double) ( ts.tv_sec + ts.tv_nsec / 1.0e9 );
#else
    /* Other */
    struct tms buf;
    return (double) times( &buf ) / (double) sysconf( _SC_CLK_TCK ) ;
#endif
}
