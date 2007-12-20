/*****************************************************************************/
/* Error.h                                                                   */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Functions used to report errors and warnings.                           */
/*****************************************************************************/

#ifndef ERROR_H
#define ERROR_H

#include "Types.h"

/*****************************************************************************/
/* Error routines.                                                           */
/*****************************************************************************/
/*                                                                           */
/* - Error   - report an error; the function does not return.  The input     */
/*             arguments are the same as is used for the C standard printf   */
/*             routine.  This function uses the current Error_Interface      */
/*             and Output_Interface.                                         */
/*                                                                           */
/* - Warning - report a warning to the used; the function returns control    */
/*             back to the calling program.  The input arguments are the     */
/*             same as is used for the C standard printf routine.  This      */
/*             function uses the current Error_Interface and                 */
/*             Output_Interface.                                             */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void) Error(const Char *fmt, ...);
FUN_DECL(Void) Warning(const Char *fmt, ...);

#endif
