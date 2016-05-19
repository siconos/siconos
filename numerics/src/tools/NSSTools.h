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

#ifndef NSSTOOLS_H
#define NSSTOOLS_H

/*!\file NSSTools.h
  Header to collect basic tools, structures definition or any usefull things for NSSpack

*/

#include "SiconosConfig.h"

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Search for the max. element of a vector
      \param[in] x the vector
      \param[in,out] sol the  solution, value of the greatest element of x
      \param[in] n  size of x
  */
  void max_part(double* x, double* sol, int n);

  /** compare two double a and b, and return the max.
   *  \param a  double*
   *  \param b  double*
   *  \param c  double*, the max
   */
  void maxf(double* a, double* b , double* c);

  /** Search for the min. element of a vector
      \param[in] x the vector
      \param[in,out] sol solution, value of the smallest element of x
      \param[in] n size of x
  */
  void min_part(double* x,  double* sol , int n);

  /** compare two double a and b, and return the min.
   *  \param a double*
   *  \param b double*
   *  \param c double*, the min
   */
  void minf(double* a, double* b, double* c);

  /** Positive part values of the components of a vector
      \param[in] n size of x
      \param[in] x the vector
      \param[out] x_plus solution vector of positive part values of x components
  */
  void pos_part(unsigned n, double* x, double* x_plus);

  /** Absolute values of the components of a vector
      \param[in] x the vector
      \param[in,out] sol solution, vector of absolute values of x components
      \param[in] n size of x
  */
  void abs_part(double* x, double* sol, int n);

  /**
      Input na, a, nb, b
      Output nc, c
      a and b: interger vectors in increasing order
      c : vector of integers of a that are not in b.
      \author Nineb Sheherazade & Dureisseix David.
  */
  void diffns(int *na, int *a, int *nb, int * b, int *nc, int *c);

  /** */
  void sortsn_(int*ddl_i, int *sort, int *n);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
