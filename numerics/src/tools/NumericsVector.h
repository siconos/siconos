/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#ifndef NumericsVector_H
#define NumericsVector_H

/*!\file NumericsVector.h
  \brief Structure definition and functions related to vector storage in Numerics
  \author Franck Perignon
*/
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include "NumericsFwd.h"
#include "SiconosConfig.h"



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Screen display of the vector content stored as a double * array
      \param m the vector to be displayed
      \param nRow the number of rows
   */
  void NV_display(double * m, int nRow);

  void NV_write_in_file_python(double * m,  int nRow, FILE* file);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
