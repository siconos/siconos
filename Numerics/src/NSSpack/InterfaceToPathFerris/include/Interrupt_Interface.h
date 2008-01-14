/*****************************************************************************/
/* Interrupt_Interface.h                                                     */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Interface used to change the behaviour of interrupts that occur within  */
/*   the algorithm call.                                                     */
/*                                                                           */
/*   User defined interrupt functions are communicated through an            */
/*   Interrupt_Interface, which should be set using the                      */
/*   Interrupt_SetInterface( ) function.                                     */
/*****************************************************************************/

#ifndef INTERRUPT_INTERFACE_H
#define INTERRUPT_INTERFACE_H

#include "Types.h"

/*****************************************************************************/
/* Interrupt_Interface declaration.                                          */
/*****************************************************************************/
/*                                                                           */
/* - interrupt_data is a user defined piece of data passed as the first      */
/*   argument to all of the interface functions.                             */
/*                                                                           */
/* - set callback function should set the interrupt handler.  Typically,     */
/*   this should only look for SIGINT interrupts.                            */
/*                                                                           */
/* - restore callback function should restore the interrupt handler to the   */
/*   previous state.  We guarantee the a set will be called before each      */
/*   restore.                                                                */
/*                                                                           */
/* - check callback function should return a positive number if a SIGINT     */
/*   has been encountered and zero otherwise.  We guarantee that this        */
/*   function will only be called between a set and restore.                 */
/*                                                                           */
/****************************************************************************/

typedef struct
{
  Void *interrupt_data;

  Void(CB_FPTR set)(Void *data);
  Void(CB_FPTR restore)(Void *data);
  Int(CB_FPTR check)(Void *data);
} Interrupt_Interface;

/*****************************************************************************/
/* Interface functions.                                                      */
/*****************************************************************************/
/*                                                                           */
/* Interface_Default      - reset the Interrupt_Interface to the default     */
/*                          implemetation.                                   */
/*                                                                           */
/* Interrupt_SetInterface - set the Interrupt_Interface.  This function      */
/*                          should be called during the system setup before  */
/*                          calling any option, creation, or algorithm       */
/*                          routines.                                        */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void) Interrupt_Default(Void);
FUN_DECL(Void) Interrupt_SetInterface(Interrupt_Interface *i);

#endif

