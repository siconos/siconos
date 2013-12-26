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
  \brief Definition of a structure to handle with variational inequality problems.
*/

/** \page viProblem Variational Inequality (VI)
 *
 * \section viIntro Problem statement
 *  Given
 * <ul>
 *   <li> an integer \f$n\f$, the dimension of the ambient space,</li>
 *   <li> a mapping \f$ F:{\mathrm{I\!R}^n \rightarrow {\mathrm{I\!R}^n, \f$</li>
 *   <li> a set  \f$ {X} \in {{\mathrm{I\!R}}}^n\f$</li>
 * </ul>
 * the variational inequality problem  is to find a vector \f$z\in{{\mathrm{I\!R}}}^n \in X\f$,
 * \f{eqnarray*}{
 * F(z)(y-z) \geq 0,\quad \text{ for all } y \in X
 * \f}
 * or equivalently,
 * \f{eqnarray*}{
 * - F(z) \in N_X(z)
 * \f}
 * where \f$N_X\f$ is the normal cone to \f$X\f$.
 *
 * Reference
 *
 * Facchinei, Francisco; Pang, Jong-Shi (2003), 
 * Finite Dimensional Variational Inequalities and Complementarity Problems, Vol. 1 & 2,
 * Springer Series in Operations Research, Berlin-Heidelberg-New York: Springer-Verlag,
 *
 * The problem is stored and given to the solver in Siconos/Numerics thanks to
 *  a C structure VariationalProblem .
 *
 *  \section viSolversList Available solvers for Variational Inequality (see VI_cst.h)
 * Use the generic function variationalInequality_driver() to call one the the specific solvers listed below:
 *
 * <ul>
 *
 * <li> variationalInequality_ExtraGradient() :Extra gradient solver.
 *       SolverId : SICONOS_VI_EG =1000, </li>
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




/** \struct VariationalInequality
 *
 */
typedef struct
{
  /** size of the VI \f$ n \f$ */
  int size;
  
  /** pointer onto env object (which is self is the simplest case)*/
  void *env;
  
  /** Function of the VI */
  void (*F)(void *self, double * x ,double *fx);  
  
  /** Projection on X of the VI */
  void (*ProjectionOnX)(void *self, double *x, double *projectionOnX); 

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
   */
  int variationalInequality_printInFile(VariationalInequality*  problem, FILE* file);

  /** read a VariationalInequalityProblem in a file (numerics .dat format)
   * \param problem the problem to read
   * \param file the target file
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

