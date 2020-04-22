/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*! \file FirstOrderNonLinearR.hpp
\brief General interface for relations.
 */

#ifndef FirstOrderNonLinearR_H
#define FirstOrderNonLinearR_H

#include "FirstOrderR.hpp"
#include "Interaction.hpp"

/** Pointer to function for plug-in for operators related to input and its gradients.*/

/** FirstOrder Non Linear Relation.
 *
 *  This is the most generic relation for First Order Dynamical Systems, with:
 * \rststar
 *
 * .. math::
 *
 *     y &=& h(X,t,\lambda,Z) \\
 *     R &=& g(X,t,\lambda,Z)
 *
 * \endrststar
 *
 *  where X, Z, and R corresponds to DynamicalSystem variables.
 *  If more than 2 DynamicalSystem are involved in the Interaction, then X = [x1 x2], Z=[z1 z2] R = [r1 r2].
 *
 *  \f$ y \f$ and \f$ \lambda \f$ are specific variables of the Interaction (see this class for more details).
 *
 * Let us define the following jacobians:
 * \rststar
 *
 * .. math::
 *
 *   C &=& \nabla_x h\\
 *   B &=& \nabla_{lambda} g\\
 *   D &=& \nabla_{lambda} h\\
 *   K &=& \nabla_x g.
 * \endrststar
 *
 *  There are 2 ways to define this relation:
 *  - by using the plugin mechanism of calling C functions
 *  - by using the inheritance mechanism (of C++ or Python) and overloading methods.
 *
 *  For the plugins, the following definitions are mandatory:
 *  - A function to compute \f$h\f$ with signature
 *  @code (double time, unsigned x_size, double *x, unsigned size_lambda, double* lambda, double *y, unsigned z_size, double *z) @endcode
 *  - A function to compute \f$g\f$ with signature
 *  @code (double time, unsigned x_size, double *x, unsigned size_lambda, double* lambda, double *r, unsigned z_size, double *z) @endcode
 *
 *  Note that the size of \f$y\f$ is the same as \f$ \lambda \f$, and the size of \f$R\f$ is the same as \f$X\f$.
 *  Thus those are not specified in the plugin function signatures.
 *
 *  For the various jacobians, there are two possibilities:
 *  If one is constant, the value may directly be set: for instance, if \f$C\f$ is constant, then one can use setCPtr to fix the value.
 *  A word of cautions: whenever a jacobian matrix is fixed using this call, then the corresponding C++ function (and not plugin) is not called anymore.
 *  A small example: if \f$ C\f$ is fixed via setCPtr, then computeJachx is never called again.
 *
 *  The other option is to use the plugin mechanism. They all share the same signature:
 *  @code (double time, unsigned x_size, double *x, unsigned size_lambda, double* lambda, double *mat, unsigned z_size, double *z) @endcode
 *  where mat is the pointer to the array of values for each Jacobian. This implies that only dense matrix are supported.
 *
 *
 */

class FirstOrderNonLinearR : public FirstOrderR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderNonLinearR);




public:

  /** basic constructor
  */
  FirstOrderNonLinearR(): FirstOrderR(RELATION::NonLinearR) {}

  /** destructor
  */
  virtual ~FirstOrderNonLinearR() {};

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction using this relation
   */
  void initialize(Interaction& inter);

  /** check sizes of the relation specific operators.
   * \param inter an Interaction using this relation
   */
  virtual void checkSize(Interaction& inter);

  /** default function to compute \f$h\f$
   * \param time    current time
   * \param x       current state variables
   * \param lambda  current nonsmooth variables
   * \param z       current auxiliary variable
   * \param[out] y  output value
   */
  virtual void computeh(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SiconosVector& y);

  /** default function to compute \f$g\f$
   * \param time    current time
   * \param x       current state variables
   * \param lambda  current nonsmooth variables
   * \param z       current auxiliary variable
   * \param[out] r  input value
  */
  virtual void computeg(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, BlockVector& r);

  /** default function to compute \f$ C = \nabla_x h \f$
   * \param time    current time
   * \param x       current state variables
   * \param lambda  current nonsmooth variables
   * \param z       current auxiliary variable
   * \param[out] C  jacobian matrix
   */
  virtual void computeJachx(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SimpleMatrix& C);

  /** default function to compute \f$ D = \nabla_{\lambda} h \f$
   * \param time    current time
   * \param x       current state variables
   * \param lambda  current nonsmooth variables
   * \param z       current auxiliary variables
   * \param[out] D  jacobian matrix
   */
  virtual void computeJachlambda(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SimpleMatrix& D);

  virtual void computeJach(double time, Interaction& inter);

  /** default function to compute \f$ B = \nabla_{\lambda}g \f$
  * \param time current time
  * \param x       current state variables
  * \param lambda  current nonsmooth variables
  * \param z       current auxiliary variables
  * \param[out] B  jacobian matrix
  */
  virtual void computeJacglambda(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SimpleMatrix& B);

  /** default function to compute \f$ K = \nabla_{\lambda}g \f$
   * \param time    current time
   * \param x       current state variables
   * \param lambda  current nonsmooth variables
   * \param z       current auxiliary variables
   * \param[out] K  jacobian matrix
   */
  virtual void computeJacgx(double time, const BlockVector& x, const SiconosVector& lambda, BlockVector& z, SimpleMatrix& K);

  virtual void computeJacg(double time, Interaction& inter);

  /** default function to compute y, using the data from the Interaction and DS
  *  \param time current time (not used)
  *  \param inter Interaction using this Relation
  *  \param level not used
  */
  virtual void computeOutput(double time, Interaction& inter, unsigned int level = 0);

  /** default function to compute r, using the data from the Interaction and DS
   *  \param time current time (not used)
   *  \param inter Interaction using this Relation
   *  \param level not used
   */
  virtual void computeInput(double time, Interaction& inter, unsigned int level = 0);

  /** return true if the relation requires the computation of residu
   * \return true if residu are required, false otherwise
   */
  virtual bool requireResidu()
  {
    return true;
  }
 
};
TYPEDEF_SPTR(FirstOrderNonLinearR)

#endif
