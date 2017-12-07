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

/*!\file ConvexQP.h
  \brief Definition of a structure to handle Convex Quadratic Problem.
*/

/** \page convexqpProblem ConvexQP
 *
 * \section convexQPintro Problem statement
 * Given
 * <ul>
 *   <li> an integer \f$n\f$, the dimension of the ambient space,</li>
 *   <li> a SDP matrix \f$ M \in  \mathrm{I\!R}^{n \times n}\f$</li>
 *   <li> a vector \f$ q \in  \mathrm{I\!R}^n\f$</li>
 *   <li> a matrix \f$ A \in  \mathrm{I\!R}^{m times n}\f$ of constraints</li>
 *   <li> a vector \f$ b \in  \mathrm{I\!R}^m\f$</li>
 *   <li> a convex set \f$ {C} \in {{\mathrm{I\!R}}}^m\f$</li>
 * </ul>
 * the convex QP problem is to find a vector \f$z\in{{\mathrm{I\!R}}}^n\f$,
 * \f{equation*}{
 *   \begin{array}{lcl}
 *     \min & & \frac{1}{2} z^T M z + z^T q \\
 *      s.t  & & A z + b  \in C \\
 *   \end{array}
 * \f}
 * and is most simple example is when \f$ b= 0 A =I\f$ and we obtain
 *    \f{equation*}{
 *  \begin{array}{lcl}
 *    \min & & \frac{1}{2} z^T M z + Z^T q \\
 *    s.t  & &  z  \in C \\
 *  \end{array}
 * \f}
 */

#ifndef CONVEXQP_H
#define CONVEXQP_H

#include "NumericsFwd.h"
#include <stdio.h>
#include "SiconosConfig.h"


/** \struct ConvexQP ConvexQP.h
 *
 */
struct ConvexQP
{
  int size; /**< size  \f$ n \f$ */
  int m; /**< m \f$ m \f$ */
  void *env; /**< pointer onto env object (which is self is the simplest case)*/
  NumericsMatrix *M; /**< Matrix M that defines the quadratic term in the cost function. **/
  double *q; /**< vector q that defines the linear term in the cost function. **/
  NumericsMatrix *A; /**< Matrix A that defines the constraints. If it is NULL, we assume that A is the identity matrix **/
  double *b;  /**< vector b that defines the constant term in the constraints. **/
  void (*ProjectionOnC)(void *self, double *x, double * PX); /**< Projection on C  */
  double normConvexQP; /**< Norm of the  problem to compute relative solution */
  int istheNormConvexQPset; /**< Boolean to know if the norm is set 
                             * If not (istheNormConvexQPset=0) it will be computed in the first call of convexQP_computeError
                             * By default, set istheNormConvexQPset =0 */
  void* set; /**< opaque struct that represent the set C (possibly empty) */
};


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** display a ConvexQPProblem
   * \param convexQP the problem to display
   */
  void convexQP_display(ConvexQP*  cqp);

  /** print a ConvexQPProblem in a file (numerics .dat format)
   * \param cqp the problem to print out
   * \param file the dest file
   * \return ok if successfull
   */
  int convexQP_printInFile(ConvexQP*  cqp, FILE* file);

  /** read a ConvexQPProblem in a file (numerics .dat format)
   * \param cqp the problem to read
   * \param file the target file
   * \return ok if successfull
   */
  int convexQP_newFromFile(ConvexQP*  cqp, FILE* file);

  /** free a ConvexQPProblem
   * \param cqp the problem to free
   */
  void freeConvexQPProblem(ConvexQP* cqp);

  /** Clear ConvexQP structure: set all pointeurs to NULL, double and int to 0.
   * \param cqp the problem to clear
   */
  void convexQP_clear(ConvexQP* cqp);

  /** new ConvexQP problem
    * \param size size of the ambient space for the CQP
    * \return a initialized ConvexQP struct
    */
  ConvexQP* convexQP_new(int size);

  /** new ConvexQP problem
    * \return an empty CQP
    */
  ConvexQP* newCQP(void);

  /** get the environment from the struct
   * \param cqp a ConvexQP problem
   * \return the environment from the struct
   */
  static inline void* convexQP_get_env(void* cqp)
  {
    return ((ConvexQP*) cqp)->env;
  }


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif

