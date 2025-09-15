/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#ifndef CONVEXQP_H
#define CONVEXQP_H
#include <stdio.h>  // for FILE

#include "NumericsFwd.h"    // for ConvexQP, NumericsMatrix
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

/** \struct ConvexQP ConvexQP.h
 *
 */
struct ConvexQP {
  int size;  /**< size  \f$ n \f$ */
  int m;     /**< m \f$ m \f$ */
  void* env; /**< pointer onto env object (which is self is the simplest case)*/
  NumericsMatrix*
      M;     /**< a n x n matrix M that defines the quadratic term in the cost function. **/
  double* q; /**< a vector q of size n that defines the linear term in the cost function. **/
  NumericsMatrix* A; /**< a m x n matrix A that defines the constraints. If it is NULL, we
                        assume that A is the identity matrix **/
  double* b; /**< a vector b of size m that defines the constant term in the constraints. **/
  void (*ProjectionOnC)(void* self, double* x, double* PX); /**< Projection on C  */
  double normConvexQP;      /**< Norm of the  problem to compute relative solution */
  int istheNormConvexQPset; /**< Boolean to know if the norm is set
                             * If not (istheNormConvexQPset=0) it will be computed in the first
                             * call of convexQP_compute_error By default, set
                             * istheNormConvexQPset =0 */
  void* set;                /**< opaque struct that represent the set C (possibly empty) */
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif
/** display a ConvexQPProblem
 *
 *  \param cqp the problem to display
 */
void convexQP_display(ConvexQP* cqp);

/** print a ConvexQPProblem in a file (numerics .dat format)
 *
 *  \param cqp the problem to print out
 *  \param file the dest file
 *  \return ok if successfull
 */
int convexQP_printInFile(ConvexQP* cqp, FILE* file);

/** read a ConvexQPProblem in a file (numerics .dat format)
 *
 *  \param cqp the problem to read
 *  \param file the target file
 *  \return ok if successfull
 */
int convexQP_newFromFile(ConvexQP* cqp, FILE* file);

/** free a ConvexQPProblem
 *
 *  \param cqp the problem to free
 */
void convexQP_free(ConvexQP* cqp);

/** Clear ConvexQP structure: set all pointeurs to NULL, double and int to 0.
 *
 *  \param cqp the problem to clear
 */
void convexQP_clear(ConvexQP* cqp);

/** new ConvexQP problem
 *
 *  \param size size of the ambient space for the CQP
 *  \return a initialized ConvexQP struct
 */
ConvexQP* convexQP_new(int size);

/** get the environment from the struct
 *
 *  \param cqp a ConvexQP problem
 *  \return the environment from the struct
 */
static inline void* convexQP_get_env(void* cqp) { return ((ConvexQP*)cqp)->env; }

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
