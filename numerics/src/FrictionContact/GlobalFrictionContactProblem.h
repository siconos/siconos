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
#ifndef GLOBALFRICTIONCONTACTPROBLEM_H
#define GLOBALFRICTIONCONTACTPROBLEM_H

/*!\file GlobalFrictionContactProblem.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
  \author Vincent Acary.
*/

#include "NumericsFwd.h"
#include "SiconosConfig.h"
#include <stdio.h>

/** \struct GlobalFrictionContactProblem GlobalFrictionContactProblem.h
 * The structure that defines a Friction-Contact (3D or 2D ) problem \f$\mathrm{PFC}(M,H,q,b,\mu)\f$  such that
 * \f{eqnarray*}{
 * \begin{cases}
 *  M v =  q +  H r \\
 *  u = H^\top v + b \\
 *    \hat u = u +\left[
 *      \left[\begin{array}{c}
 *          \mu^\alpha \|u^\alpha_{T}\|\\
 *         0 \\
 *         0
 *        \end{array}\right]^T, \alpha = 1 \ldots n_c
 *    \right]^T \\ \                                \
 *    C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu}
 * \end{cases}
 * \f}
 * and the set \f$C^{\alpha,\star}_{\mu^\alpha}\f$ is its dual.

*/
struct GlobalFrictionContactProblem
{
  /** dimension \f$d=2\f$ or \f$d=3\f$ of the contact space (3D or 2D ) */
  int dimension;
  /** the number of contacts \f$ n_c \f$ */
  int numberOfContacts;
  /** \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
      a matrix with \f$ n\f$ stored in NumericsMatrix structure */
  NumericsMatrix* M;
  /**  \f${H} \in {{\mathrm{I\!R}}}^{n \times m} \f$,
      a matrix with \f$ m = d  n_c\f$ stored in NumericsMatrix structure */
  NumericsMatrix* H;
  /** \f${q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double* q;
  /** \f${b} \in {{\mathrm{I\!R}}}^{m} \f$ */
  double* b;
  /** mu \f${\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of friction coefficients
      (\f$ n_c =\f$ numberOfContacts) */
  double* mu;
  /** opaque environment, solver specific */
  void* env; 
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void globalFrictionContact_display(GlobalFrictionContactProblem*  problem);


  int globalFrictionContact_printInFile(GlobalFrictionContactProblem*  problem, FILE* file);

  int globalFrictionContact_newFromFile(GlobalFrictionContactProblem*  problem, FILE* file);

  static inline void globalFrictionContact_null(GlobalFrictionContactProblem*  problem)
  {
    problem->M = NULL;
    problem->H = NULL;
    problem->q = NULL;
    problem->b = NULL;
    problem->mu = NULL;
    problem->env = NULL;
    problem->numberOfContacts = 0;
    problem->dimension = 0;

  }

 void freeGlobalFrictionContactProblem(GlobalFrictionContactProblem* problem);


  

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif
#endif
