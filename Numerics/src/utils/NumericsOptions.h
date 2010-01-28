/* Siconos-Numerics, Copyright INRIA 2005-2010.
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

#ifndef NumericsOptions_H
#define NumericsOptions_H

/*!\file NumericsOptions.h
  \brief General options for Numerics functions, structures and so on (mainly used to send information from Kernel to Numerics).
  \author Franck Perignon
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Structure used to set general options of Numerics functions, structures and so on.
    \param verbose mode (0: off, 1: on)
*/
typedef struct
{
  int verboseMode;
} NumericsOptions;


/* Verbose mode */
extern int verbose;

#ifdef __cplusplus
extern "C" {
#endif

  /* Set verbose mode in numerics
     \param newVerboseMode, 0 no verbose, 1 verbose.
   */
  void setNumericsVerbose(int newVerboseMode);

  /* Set global option for numerics
     \param opt, a NumericsOptions structure
   */
  void setNumericsOptions(NumericsOptions* opt);

  /* message output and exit with error
     \param name of the function where error occurs
     \param output message
  */
  void numericsError(char* functionName, char* message);

  /* message output without exit
     \param name of the function where warning occurs
     \param output message
  */
  void numericsWarning(char* functionName, char* message);

#ifdef __cplusplus
}
#endif


#endif
