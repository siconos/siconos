/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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
 * the variational inequality problem is to find a vector \f$z\in{{\mathrm{I\!R}}}^n \in X\f$,
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
 * <ul>
 *
 * <li> variationalInequality_ExtraGradient() : Extra gradient solver.
 *      SolverId: SICONOS_VI_EG = 1000. </li>
 * <li> variationalInequality_FixedPointProjection() : Fixed-point solver.
 *      SolverId: SICONOS_VI_EG = 1001. </li>
 * <li> variationalInequality_HyperplaneProjection() : Hyperplane Projection
 *      based Solver. SolverId: SICONOS_VI_HP_STR = 1002. </li>
 * <li> variationalInequality_box_newton_QiLSA : Solver using the merit
 * function proposed by Qi for box-constrained VI. SolverId:
 * SICONOS_VI_BOX_QI_STR = 1003
 *
 * </ul>
 * (see the functions/solvers list in VariationalInequality_Solvers.h)
 *
 *
 */

#ifndef VARIATIONALINEQUALITY_H
#define VARIATIONALINEQUALITY_H

#include "NumericsMatrix.h"

typedef void * (FVIPtr)(void*, double *, double *);
typedef void (*ptrFunctionVI)(void *self, int n, double* x ,double* fx);
enum VI_SET_TYPE { VI_SET_IS_BOX, VI_SET_IS_POLYHEDRON };


/** \struct VariationalInequality
 *
 */
typedef struct VariationalInequality_
{
  int size; /**< size of the VI \f$ n \f$ */
  void *env; /**< pointer onto env object (which is self is the simplest case)*/
  ptrFunctionVI F; /**< Function of the VI */
  ptrFunctionVI compute_nabla_F; /**< Function to compute the jacobian of F */
  void (*ProjectionOnX)(void *self, double *x, double * PX); /**< Projection on X of the VI */
  double normVI; /**< Norm of the VI problem to compute relative solution */
  int istheNormVIset; /**< Boolean to know if the norm is set 
   * If not (istheNormVIset=0) it will be computed in the first call of variationalInequality_computeError
   * By default, set istheNormVIset =0 */
  void* set; /**< opaque struct that represent the set K (possibly empty) */
  double* nabla_F; /**< storage for \f$\nabla_x F\f$*/
} VariationalInequality;


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


  /* /\** new VariationalInequality  from minimal set of data */
  /*  * */
  /*  *\/ */
  /* VariationalInequality* variationalInequality_new(int size, CallbackVI * callback); */



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

