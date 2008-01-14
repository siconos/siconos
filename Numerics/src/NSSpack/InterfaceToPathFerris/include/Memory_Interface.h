/*****************************************************************************/
/* Memory_Interface.h                                                        */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Interface used to change the behaviour of memory allocation and         */
/*   deallocation.                                                           */
/*                                                                           */
/*   User defined allocate and free functions are communicated through a     */
/*   Memory_Interface, which should be set using the Memory_SetInterface( )  */
/*   function.                                                               */
/*****************************************************************************/

#ifndef MEMORY_INTERFACE_H
#define MEMORY_INTERFACE_H

#include "Types.h"

/*****************************************************************************/
/* Memory_Interface declaration.                                             */
/*****************************************************************************/
/*                                                                           */
/* - memory_data is a user defined piece of data passed as the first         */
/*   argument to all of the interface functions.                             */
/*                                                                           */
/* - allocate callback function should allocate n bytes of general memory    */
/*   and return the newly allocated space.  Return NULL if the space cannot  */
/*   be allocated.                                                           */
/*                                                                           */
/* - allocate_factors callback function should allocate n bytes of storage   */
/*   for use in matrix factorizations and return the newly allocated space.  */
/*   Return NULL if the space cannot be allocated.                           */
/*                                                                           */
/* - free callback function should deallocate the general memory that was    */
/*   previously allocated thorugh the allocate callback.                     */
/*                                                                           */
/* - free_factors callback function should deallocate the factorization      */
/*   storage previously allocated through the allocate_factors callback.     */
/*                                                                           */
/* NOTE: PATH will only ever have one set of factors allocated at any given  */
/* time.  However, the allocate_factors and free_factors may be called       */
/* several times in succession.                                              */
/*                                                                           */
/*****************************************************************************/

typedef struct
{
  Void *memory_data;

  Void *(CB_FPTR allocate)(Void *data, Unsigned Long n);
  Void *(CB_FPTR allocate_factors)(Void *data, Unsigned Long n);
  Void(CB_FPTR free)(Void *data, Void *v);
  Void(CB_FPTR free_factors)(Void *data, Void *v);
} Memory_Interface;

/*****************************************************************************/
/* Interface functions.                                                      */
/*****************************************************************************/
/*                                                                           */
/* Memory_Default      - reset the Memory_Interface to the default           */
/*                       implemetation.                                      */
/*                                                                           */
/* Memory_SetInterface - set the Memory_Interface.  This function should be  */
/*                       called during the system setup before calling any   */
/*                       option, creation, or algorithm routines.            */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(void) Memory_Default(Void);
FUN_DECL(void) Memory_SetInterface(Memory_Interface *i);

#endif

