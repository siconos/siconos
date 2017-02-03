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

/*!\file VariationalInequality.h
  \brief Definition of a structure to handle Variational Inequalities (VI).
*/

/** \page viProblem Variational Inequality (VI)
 *
 * \section viIntro Problem statement
 *  Given
 * <ul>
 *   <li> an integer \f$n\f$, the dimension of the ambient space,</li>
 *   <li> a mapping \f$ F\colon \mathrm{I\!R}^n \rightarrow \mathrm{I\!R}^n\f$</li>
 *   <li> a set  \f$ {X} \in {{\mathrm{I\!R}}}^n\f$</li>
 * </ul>
 * the variational inequality problem is to find a vector \f$z\in{{\mathrm{I\!R}}}^n\f$,
 * \f{equation*}{
 * F(z)^T(y-z) \geq 0,\quad \text{ for all } y \in X
 * \f}
 * or equivalently,
 * \f{equation*}{
 * - F(z) \in \mathcal{N}_X(z)
 * \f}
 * where \f$\mathcal{N}_X\f$ is the normal cone to \f$X\f$ at \f$z\f$.
 *
 * Reference
 *
 * Facchinei, Francisco; Pang, Jong-Shi (2003), 
 * <i>Finite Dimensional Variational Inequalities and Complementarity Problems</i>, Vol. 1 & 2,
 * Springer Series in Operations Research, Berlin-Heidelberg-New York: Springer-Verlag.
 *
 * The problem is stored and given to the solver in Siconos/Numerics thanks to
 *  a C structure VariationalProblem.
 *
 *  \section viSolversList Available solvers for Variational Inequality
 * Use the generic function variationalInequality_driver() to call one the the specific solvers listed below:
 * - variationalInequality_ExtraGradient() : Extra gradient solver.
 *      SolverId: SICONOS_VI_EG = 1000.
 * - variationalInequality_FixedPointProjection() : Fixed-point solver.
 *      SolverId: SICONOS_VI_EG = 1001.
 * - variationalInequality_HyperplaneProjection() : Hyperplane Projection
 *      based Solver. SolverId: SICONOS_VI_HP_STR = 1002.
 * - variationalInequality_box_newton_QiLSA() : Solver using the merit
 * function proposed by Qi for box-constrained VI. SolverId:
 * SICONOS_VI_BOX_QI_STR = 1003
 *
 * (see the functions/solvers list in VariationalInequality_Solvers.h)
 *
 */

#ifndef VARIATIONALINEQUALITY_H
#define VARIATIONALINEQUALITY_H

#include "NumericsFwd.h"
#include <stdio.h>
#include "SiconosConfig.h"

typedef void * (FVIPtr)(void*, double *, double *);
typedef void (*ptrFunctionVI)(void *self, int n, double* x, double* fx);
typedef void (*ptrFunctionVI_nabla)(void *self, int n, double* x, NumericsMatrix* nabla_F);


/** \struct VariationalInequality VariationalInequality.h
 * 
 */
struct VariationalInequality
{
  int size; /**< size of the VI \f$ n \f$ */
  void *env; /**< pointer onto env object (which is self is the simplest case)*/
  ptrFunctionVI F; /**< Function of the VI */
  ptrFunctionVI_nabla compute_nabla_F; /**< Function to compute the jacobian of F */
  void (*ProjectionOnX)(void *self, double *x, double * PX); /**< Projection on X of the VI */
  double normVI; /**< Norm of the VI problem to compute relative solution */
  int istheNormVIset; /**< Boolean to know if the norm is set 
   * If not (istheNormVIset=0) it will be computed in the first call of variationalInequality_computeError
   * By default, set istheNormVIset =0 */
  void* set; /**< opaque struct that represent the set K (possibly empty) */
  NumericsMatrix* nabla_F; /**< storage for \f$\nabla_x F\f$*/
};


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** display a VariationalInequalityProblem
   * \param problem the problem to display
   */
  void variationalInequality_display(VariationalInequality*  problem);

  /** print a VariationalInequalityProblem in a file (numerics .dat format)
   * \param problem the problem to print out
   * \param file the dest file
   * \return ok if successfull
   */
  int variationalInequality_printInFile(VariationalInequality*  problem, FILE* file);

  /** read a VariationalInequalityProblem in a file (numerics .dat format)
   * \param problem the problem to read
   * \param file the target file
   * \return ok if successfull
   */
  int variationalInequality_newFromFile(VariationalInequality*  problem, FILE* file);

  /** free a VariationalInequalityProblem
   * \param problem the problem to free
   */
  void freeVariationalInequalityProblem(VariationalInequality* problem);

  /** Clear VariationalInequality structure: set all pointeurs to NULL, double and int to 0.
   * \param vi the problem to clear
   */
  void variationalInequality_clear(VariationalInequality* vi);

  /** new VariationalInequality problem
    * \param size size of the ambient space for the VI
    * \return a initialized VariationalInequality struct
    */
  VariationalInequality* variationalInequality_new(int size);

  /** get the environment from the struct
   * \param problem a VariationalInequality problem
   * \return the environment from the struct
   */
  static inline void* VI_get_env(void* problem)
  {
    return ((VariationalInequality*) problem)->env;
  }


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

