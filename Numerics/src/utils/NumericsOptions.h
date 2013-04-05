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

#ifndef NumericsOptions_H
#define NumericsOptions_H

/*!\file NumericsOptions.h
  \brief General options for Numerics functions, structures and so on 
         (mainly used to send information from Kernel to Numerics).
  \author Franck Perignon
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OUTPUT_ON_ERROR 4

/** Structure used to set general options of Numerics functions,
  * structures and so on.
  *
  *  \param verboseMode (0: off, 1: on)
  *  \param outputMode  (0: none, 1: c file, 2: dat file, 3: fclib file) 
  *                      with OUTPUT_ON_ERROR flag : only on solver error
  *  \param fileName    output filename

  *  \param description experience description
  *  \param counter     a counter for file numeration
  */

/* note : swig is ok with char[n] (wrapper checks string length) */
typedef struct
{
  int verboseMode;
  int outputMode;
  char fileName[256];   
  char title[256];
  char description[1024];
  char math_info[1024];
  int counter;
} NumericsOptions;




/* Verbose mode */
extern int verbose;

#ifdef __cplusplus
extern "C"
{
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


  /* set default values for NumericsOptions */
  void setDefaultNumericsOptions(NumericsOptions* opt);

#ifdef __cplusplus
}
#endif


#endif
