/*****************************************************************************/
/* PathOptions.h                                                             */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Function used to add the PathOptions to the option set.                 */
/*****************************************************************************/

#ifndef PATHOPTIONS_H
#define PATHOPTIONS_H

#include "Types.h"
#include "Options.h"

/*****************************************************************************/
/* Main functionality.                                                       */
/*****************************************************************************/
/*                                                                           */
/* Path_AddOptions - add the options for Path to the allocated option        */
/*                   interface.  This should be done right after allocating  */
/*                   the option set and before using the Option_Default( )   */
/*                   function.                                               */
/*                                                                           */
/*****************************************************************************/

FUN_DECL(Void) Path_AddOptions(Options_Interface *i);

#endif

