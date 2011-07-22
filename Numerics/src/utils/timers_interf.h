/* Siconos-Numerics, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#ifndef TIMERS_INTERF_H
#define TIMERS_INTERF_H

#include "NumericsConfig.h"

/*!\file timers.h
  \brief A common interface to different timers, see test/testop3x3.c for examples
*/

#ifdef __GNUC__
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif

#define HAVE_FFTW_CYCLE_H
#define _BSD_SOURCE

#ifdef WITH_TIMERS
#if(                               \
  !defined(TIMER_CLOCK_GETTIME) && \
  !defined(TIMER_ATL_WALLTIME) &&  \
  !defined(TIMER_ATL_CPUTIME) &&   \
  !defined(TIMER_CLOCK) &&         \
  !defined(TIMER_FFTW_CYCLE) &&    \
  !defined(TIMER_SYSTIMES))
/* the default timer */
#ifdef HAVE_TIME_H
#define TIMER_CLOCK_GETTIME 1
#else
#undef WITH_TIMERS
#endif
#endif
#endif

#ifdef HAVE_TIME_H
#include <time.h>
#endif

#ifdef HAVE_SYSTIMES_H
#include <sys/times.h>
#endif

#ifdef HAVE_TIME_H
#define DECL_TIMER_CLOCK_GETTIME(T)                                     \
  struct timespec __timer__##T##1[1],__timer__##T##2[1]; long T##delta;
#define START_TIMER_CLOCK_GETTIME(T)                          \
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,__timer__##T##1)
#define STOP_TIMER_CLOCK_GETTIME(T)                           \
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,__timer__##T##2)
#define ELAPSED_CLOCK_GETTIME(T) \
  ({T##delta= __timer__##T##2->tv_nsec-__timer__##T##1->tv_nsec;\
       if (T##delta < 0) T##delta += 1000000000;\
       (T##delta)*1e-9;})
#define TIMER_TICK_CLOCK_GETTIME long
#endif

#ifdef HAVE_ATLAS_AUX_H
#include <atlas_aux.h>
#define DECL_TIMER_ATL_WALLTIME(T)                  \
  double __timer__##T##1,__timer__##T##2
#define START_TIMER_ATL_WALLTIME(T)             \
  __timer__##T##1 = ATL_walltime()
#define STOP_TIMER_ATL_WALLTIME(T)              \
  __timer__##T##2 = ATL_walltime()
#define ELAPSED_ATL_WALLTIME(T)                 \
  __timer__##T##2-__timer__##T##1
#define TIMER_TICK_ATL_WALLTIME double


#define DECL_TIMER_ATL_CPUTIME(T)                   \
  double __timer__##T##1,__timer__##T##2
#define START_TIMER_ATL_CPUTIME(T)              \
  __timer__##T##1 = ATL_cputime()
#define STOP_TIMER_ATL_CPUTIME(T)               \
  __timer__##T##2 = ATL_cputime()
#define ELAPSED_ATL_CPUTIME(T)                  \
  __timer__##T##2-__timer__##T##1
#define TIMER_TICK_ATL_CPUTIME double

#endif

#ifdef HAVE_TIME_H
#define DECL_TIMER_CLOCK(T)                           \
  clock_t __timer__##T##1,__timer__##T##2
#define START_TIMER_CLOCK(T)                    \
  __timer__##T##1 = clock()
#define STOP_TIMER_CLOCK(T)                     \
  __timer__##T##2 = clock()
#define ELAPSED_CLOCK(T)                            \
  (__timer__##T##2-__timer__##T##1)*1e-6
#define TIMER_TICK_CLOCK clock_t
#endif

#ifdef HAVE_FFTW_CYCLE_H

#define asm __asm
#include "fftw_cycle.h"
#define DECL_TIMER_FFTW_CYCLE(T)                              \
  ticks __timer__##T##1, __timer__##T##2; double __timer__##T##total MAYBE_UNUSED =0.
#define START_TIMER_FFTW_CYCLE(T)               \
  __timer__##T##1 = getticks()
#define STOP_TIMER_FFTW_CYCLE(T)                \
  __timer__##T##2 = getticks()
#define ELAPSED_FFTW_CYCLE(T)                                           \
  ({__timer__##T##total+=elapsed(__timer__##T##2,__timer__##T##1); __timer__##T##total ;})
#define TIMER_TICK_FFTW_CYCLE double

#endif

#ifdef HAVE_SYSTIMES_H
#define DECL_TIMER_SYSTIMES(T)                               \
  struct tms __timer__##T##1[1], __timer__##T##2[1]
#define START_TIMER_SYSTIMES(T)                \
  times(__timer__##T##1)
#define STOP_TIMER_SYSTIMES(T)                 \
  times(__timer__##T##2)
#define ELAPSED_SYSTIMES(T)                                        \
  __timer__##T##2->tms_utime - __timer__##T##1->tms_utime
#define TIMER_TICK_SYSTIMES clock_t

#endif

#if defined(WITH_TIMERS)
#if defined(TIMER_CLOCK_GETTIME) && defined(_POSIX_TIMERS)
#define DECL_TIMER DECL_TIMER_CLOCK_GETTIME
#define START_TIMER START_TIMER_CLOCK_GETTIME
#define STOP_TIMER STOP_TIMER_CLOCK_GETTIME
#define ELAPSED ELAPSED_CLOCK_GETTIME
#define TIMER_TICK TIMER_TICK_CLOCK_GETTIME
#define TIMER_NAME "clock_gettime"
#elif defined(TIMER_ATL_WALLTIME)
#define DECL_TIMER DECL_TIMER_ATL_WALLTIME
#define START_TIMER START_TIMER_ATL_WALLTIME
#define STOP_TIMER STOP_TIMER_ATL_WALLTIME
#define ELAPSED ELAPSED_ATL_WALLTIME
#define TIMER_TICK TIMER_TICK_ATL_WALLTIME
#define TIMER_NAME "atl_walltime"
#elif defined(TIMER_ATL_CPUTIME)
#define DECL_TIMER DECL_TIMER_ATL_CPUTIME
#define START_TIMER START_TIMER_ATL_CPUTIME
#define STOP_TIMER STOP_TIMER_ATL_CPUTIME
#define ELAPSED ELAPSED_ATL_CPUTIME
#define TIMER_TICK TIMER_TICK_ATL_CPUTIME
#define TIMER_NAME "atl_cputime"
#elif defined(TIMER_CLOCK)
#define DECL_TIMER DECL_TIMER_CLOCK
#define START_TIMER START_TIMER_CLOCK
#define STOP_TIMER STOP_TIMER_CLOCK
#define ELAPSED ELAPSED_CLOCK
#define TIMER_TICK TIMER_TICK_CLOCK
#define TIMER_NAME "clock"
#elif defined(TIMER_FFTW_CYCLE)
#define DECL_TIMER DECL_TIMER_FFTW_CYCLE
#define START_TIMER START_TIMER_FFTW_CYCLE
#define STOP_TIMER STOP_TIMER_FFTW_CYCLE
#define ELAPSED ELAPSED_FFTW_CYCLE
#define TIMER_TICK TIMER_TICK_FFTW_CYCLE
#define TIMER_NAME "fftw_cycle"
#elif defined(TIMER_SYSTIMES)
#define DECL_TIMER DECL_TIMER_SYSTIMES
#define START_TIMER START_TIMER_SYSTIMES
#define STOP_TIMER STOP_TIMER_SYSTIMES
#define ELAPSED ELAPSED_SYSTIMES
#define TIMER_TICK TIMER_TICK_SYSTIMES
#define TIMER_NAME "systimes"
#else
#undef WITH_TIMERS
#endif
#else
#define DECL_TIMER(X)
#define START_TIMER(X)
#define STOP_TIMER(X)
#define ELAPSED(X) 0.
#endif

#if defined(TIMER_TICK)
#define DECL_TIMER_TICK(X) TIMER_TICK X
#define GET_ELAPSED(TIMER,X) do {X = ELAPSED(TIMER);} while(0)
#define PRINT_ELAPSED(TIMER) do {printf("%s %s:%g\n",TIMER_NAME,#TIMER,(double) ELAPSED(TIMER));} while(0)
#else
#define DECL_TIMER_TICK(X)
#define GET_ELAPSED(TIMER,X)
#define PRINT_ELAPSED(X)
#endif

#endif



