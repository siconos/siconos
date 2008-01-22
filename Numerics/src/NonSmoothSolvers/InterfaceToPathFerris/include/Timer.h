/*****************************************************************************/
/* Timer.h                                                                   */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Functions used obtain the time spent in a part of the code.             */
/*****************************************************************************/

#ifndef TIMER_H
#define TIMER_H

#include "Types.h"

/*****************************************************************************/
/* Timer allocation and deallocation routines.                               */
/*****************************************************************************/
/*                                                                           */
/* - Timer_Create  - allocate a new timer instance and return it.            */
/*                                                                           */
/* - Timer_Destroy - deallocate the previously allocated timer instance      */
/*                   indicated.                                              */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void *) Timer_Create(Void);
FUN_DECL(Void) Timer_Destroy(Void *v);

/*****************************************************************************/
/* Timer routines.                                                           */
/*****************************************************************************/
/*                                                                           */
/* - Timer_Start - set the initial time for the indicated allocated timer to */
/*                 the current time.                                         */
/*                                                                           */
/* - Timer_Query - obtain the number of seconds that have elapsed since the  */
/*                 indicated, allocated timer has been started.              */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void) Timer_Start(Void *v);
FUN_DECL(Double) Timer_Query(Void *v);

#endif

