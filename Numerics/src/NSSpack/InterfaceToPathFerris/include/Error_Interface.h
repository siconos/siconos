/*****************************************************************************/
/* Error_Interface.h                                                         */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Interface used to change the behaviour of error and warning reporting.  */
/*                                                                           */
/*   User defined error and warning functions are communicated through an    */
/*   Error_Interface, which should be set using the Error_SetInterface( )    */
/*   function.                                                               */
/*****************************************************************************/

#ifndef ERROR_INTERFACE_H
#define ERROR_INTERFACE_H

#include "Types.h"

/*****************************************************************************/
/* Error_Interface declaration.                                              */
/*****************************************************************************/
/*                                                                           */
/* - error_data is a user defined piece of data passed as the first argument */
/*   to all of the interface functions.                                      */
/*                                                                           */
/* - error callback function should not return.  The msg argument indicates  */
/*   the problem.  Before the handler is invoked, the message will already   */
/*   have been output to the log and status files.                           */
/*                                                                           */
/* - warning callback function should return.  The msg argument indicates    */
/*   the warning.  Before the handler is invoked, the message will already   */
/*   have been output to the log and status files.                           */
/*                                                                           */
/*****************************************************************************/

typedef struct
{
  Void *error_data;

  Void(CB_FPTR error)(Void *data, Char *msg);
  Void(CB_FPTR warning)(Void *data, Char *msg);
} Error_Interface;

/*****************************************************************************/
/* Interface functions.                                                      */
/*****************************************************************************/
/*                                                                           */
/* Error_Default      - reset the Error_Interface to the default             */
/*                      implemetation.                                       */
/*                                                                           */
/* Error_SetInterface - set the Error_Interface.  This function should be    */
/*                      called during the system setup before calling any    */
/*                      option, creation, or algorithm routines.             */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void) Error_Default(Void);
FUN_DECL(Void) Error_SetInterface(Error_Interface *i);

#endif

