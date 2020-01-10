#ifndef __WTIME_H__
#define __WTIME_H__

/*
 * Prototype function that returns the number of seconds since a fixed
 * arbitrary time in the past.  Resolution varies depending on what
 * functions are available during compilation.  If MACH or POSIX timers
 * are available timer should have an accuracy of about 1 nanosecond,
 * otherwise the resolution will only be about 10 milliseconds.
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 */

#ifdef __cplusplus
extern "C" {
#endif

double wtime( void );

#ifdef __cplusplus
}
#endif

#endif /* __WTIME_H__ */
