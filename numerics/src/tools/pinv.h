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

#ifndef Pinv_H
#define Pinv_H

#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** Compute the pseudo-inverse of dense matrix with column major storage
  \param A the matrix to be inversed
  \param n the number of rows of A
  \param m the number of columns of A
  \param tolerance threshold used to validate the computation: if the error is less than this value, the computation is considered valid
  */
  double pinv(double * A, int n, int m, double tolerance);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif
#endif
