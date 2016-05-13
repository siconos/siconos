/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
#include "misc.h"

#define OUTPUT_ON_ERROR 4
#define OUTPUT_ON_ITER  8

/* note : swig is ok with char[n] (wrapper checks string length) */


/** Structure used to set general options of Numerics functions,
*/
typedef struct
{
  int verboseMode; /**< 0: off, 1: on */
} NumericsOptions;




/* Verbose mode */
extern int verbose;

#ifdef __cplusplus
extern "C"
{
#endif

  /* Set verbose mode in numerics
     \param newVerboseMode 0 no verbose, 1 verbose.
   */
  void setNumericsVerbose(int newVerboseMode);

  /* Set global option for numerics
     \param opt a NumericsOptions structure
   */
  void setNumericsOptions(NumericsOptions* opt);

  /* message output and exit with error
     \param functionName name of the function where error occurs
     \param message output message
  */
  void numericsError(char* functionName, char* message) NO_RETURN;

  /* message output without exit
     \param functionName name of the function where warning occurs
     \param message output message
  */
  void numericsWarning(char* functionName, char* message);


  /* set default values for NumericsOptions
   * \param opt a NumericsOptions structure
   */
  void setDefaultNumericsOptions(NumericsOptions* opt);

#ifdef __cplusplus
}
#endif


#endif
