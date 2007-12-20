/*****************************************************************************/
/* Memory.h                                                                  */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Functions used to allocate and deallocate memory.                       */
/*****************************************************************************/

#ifndef MEMORY_H
#define MEMORY_H

#include "Types.h"

/*****************************************************************************/
/* Memory routines.                                                          */
/*****************************************************************************/
/*                                                                           */
/* - Memory_Allocate        - allocate n bytes of general memory.  This      */
/*                            function uses the current Memory_Interface.    */
/*                            The function does not return if the memory     */
/*                            cannot be allocated.                           */
/*                                                                           */
/* - Memory_AllocateFactors - allocate n bytes for use in matrix             */
/*                            factorizations.  This function uses the        */
/*                            current Memory_Interface.  The function does   */
/*                            not return if the memory cannot be allocated.  */
/*                                                                           */
/* - Memory_Free            - free the general memory allocated with the     */
/*                            Memory_Allocate( ) function.  This function    */
/*                            uses the current Memory_Interface.             */
/*                                                                           */
/* - Memory_FreeFactors     - free the memory allocated with the             */
/*                            Memory_AllocateFactors( ) function.  This      */
/*                            function uses the current Memory_Interface.    */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void *) Memory_Allocate(Unsigned Long n);
FUN_DECL(Void *) Memory_AllocateFactors(Unsigned Long n);
FUN_DECL(Void) Memory_Free(Void *v);
FUN_DECL(Void) Memory_FreeFactors(Void *v);

#endif

