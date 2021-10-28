/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#ifndef NONSMOOTHNEWTONNEIGH_H
#define NONSMOOTHNEWTONNEIGH_H

#include "NonSmoothNewton.h"  // for NewtonFunctionPtr
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

/*!\file NonSmoothNewtonNeighbour.h
  Typedef and functions declarations related to non-smooth Newton solver

  Solve \f$ \phi(z) = 0 \f$ using a Newton method.

  The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
  for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

 */



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  double* nonSmoothNewtonNeighInitMemory(int n, double * dWork, int * iWork);

  int nonSmoothNewtonNeigh(int, double*, NewtonFunctionPtr*, NewtonFunctionPtr*, int*, double*);

  int nonSmoothNewtonNeigh_getNbIWork(int n, int m);
  int nonSmoothNewtonNeigh_getNbDWork(int n, int m);

  /*only for debug*/
  void NSNN_thisIsTheSolution(int n, double * z);
  void  NSNN_reset(void);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
