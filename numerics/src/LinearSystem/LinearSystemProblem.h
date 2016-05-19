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
#ifndef LINEARSYSTEM_PROBLEM_H
#define LINEARSYSTEM_PROBLEM_H

/*!\file LinearSystemProblem.h
  \brief Structure used to define a Linear Complementarity Problem

  \author O. Bonnefon
*/

/* \page Linear System problems

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ z + q = w =0\\
  \end{array}
  \right.
  \f$

  \f$ w, z, q\f$ are vectors of size n and \f$ M \f$ is a nXn matrix.


  Use the generic function LinearSystem_driver() to call one the the specific solvers listed below:



*/

#include "NumericsMatrix.h"
extern int SICONOS_LS_0;
/** Linear Complementarity Problem elements
 */
typedef struct
{
  int size; /**< dim of the problem */
  NumericsMatrix* M; /**< matrix of the linear system */
  double * q; /**< vector */
} LinearSystemProblem;




#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  double LinearSystem_computeError(LinearSystemProblem* problem, double *z);
  int LinearSystem_getNbDwork(LinearSystemProblem* problem, SolverOptions* options);
  int LinearSystem_getNbIwork(LinearSystemProblem* problem, SolverOptions* options);
  int LinearSystem_setDefaultSolverOptions(LinearSystemProblem* problem, SolverOptions* options, int solverId);
  void LinearSystem_freeProblem(LinearSystemProblem *problem);
  int LinearSystem_newFromFile(LinearSystemProblem* problem, FILE* file);
  void LinearSystem_display(LinearSystemProblem* p);
  int LinearSystem_driver(LinearSystemProblem* problem, double *z , double *w, SolverOptions* options);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

