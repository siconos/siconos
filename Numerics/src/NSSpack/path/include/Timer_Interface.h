/*****************************************************************************/
/* Timer_Interface.h                                                         */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Interface used to change the behaviour of timing routines.              */
/*                                                                           */
/*   User defined timers are communicated through a Timer_Interface, which   */
/*   should be set using the Timer_SetInterface( ) function.                 */
/*****************************************************************************/

#ifndef TIMER_INTERFACE_H
#define TIMER_INTERFACE_H

#include "Types.h"

/*****************************************************************************/
/* Timer_Interface declaration.                                              */
/*****************************************************************************/
/*                                                                           */
/* - timer_data is a user defined piece of data passed as the first argument */
/*   to all of the interface functions.                                      */
/*                                                                           */
/* - create callback function should return a newly allocated timer.         */
/*                                                                           */
/* - destroy callback function should deallocate the indicated timer.        */
/*                                                                           */
/* - start callback function should initialize the time for the indicated    */
/*   timer to the current time.                                              */
/*                                                                           */
/* - query callback function shound return the number of seconds that have   */
/*   elapsed since the inticated timer was started.                          */
/*                                                                           */
/*****************************************************************************/

typedef struct
{
  Void *timer_data;

  Void  *(CB_FPTR create)(Void *data);
  Void(CB_FPTR destroy)(Void *data, Void *v);

  Void(CB_FPTR start)(Void *data, Void *v);
  Double(CB_FPTR query)(Void *data, Void *v);
} Timer_Interface;

/*****************************************************************************/
/* Interface functions.                                                      */
/*****************************************************************************/
/*                                                                           */
/* Timer_Default      - reset the Timer_Interface to the default             */
/*                      implemetation.                                       */
/*                                                                           */
/* Timer_SetInterface - set the Timer_Interface.  This function should be    */
/*                      called during the system setup before calling any    */
/*                      option, creation, or algorithm routines.             */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(void) Timer_Default(Void);
FUN_DECL(void) Timer_SetInterface(Timer_Interface *p);

#endif

