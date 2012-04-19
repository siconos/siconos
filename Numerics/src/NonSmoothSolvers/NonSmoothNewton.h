/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef NONSMOOTHNEWTON_H
#define NONSMOOTHNEWTON_H

/*!\file NonSmoothNewton.h
  Typedef and functions declarations related to non-smooth Newton solver

  Solve \f$ \phi(z) = 0 \f$ using a Newton method.

  The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth Newton-type methods
  for nonlinear complementarity problems", in Computational Optimization and Applications, 11, 227-251 (1998).

  \author Houari Khenous - Franck Perignon

 */

/* Pointer to function that corresponds to the function \f$ \phi \f$ */
typedef void (*NewtonFunctionPtr)(int, double*, double*, int);


#ifdef __cplusplus
extern "C"
{
#endif
  /** Armijo Linesearch
   *  \param n  size of the vector z
   *  \param z unknown vector
   *  \param dir  search direction
   *  \param psi_k initial value of the merit function
   *  \param descentCondition descent condition
   *  \param phi pointer to function used to compute phi(z)
   */
  void linesearch_Armijo(int n, double *z, double* dir, double psi_k,
                         double descentCondition, NewtonFunctionPtr* phi);


  /** Newton solver with line Search
   \param n size of the vector z
   \param z the vector z, unknown vector, in-out argument
   \param phi pointer to \f$ \phi \f$ function
   \param jacobianPhi pointer to \f$ \nabla_z \phi(z) \f$ function
   \param iparam vector of int parameters:\n
     - [0] : max. number of iterations
     - [1] : number of iterations processed
   \param dparam vector of double parameters:\n
     - [0]: tolerance
     - [1]: error
  */
  int nonSmoothNewton(int n, double* z, NewtonFunctionPtr* phi,
                      NewtonFunctionPtr* jacobianPhi,
                      int* iparam, double* dparam);

  /** Newton solver without line Search
  \param n size of the vector z
  \param z the vector z, unknown vector, in-out argument
  \param phi pointer to \f$ \phi \f$ function
  \param jacobianPhi pointer to \f$ \nabla_z \phi(z) \f$ function
  \param iparam vector of int parameters:\n
   - [0] : max. number of iterations
   - [1] : number of iterations processed
  \param dparam vector of double parameters:\n
   - [0]: tolerance
   - [1]: error
  */
  int nonSmoothDirectNewton(int n, double* z, NewtonFunctionPtr* phi,
                            NewtonFunctionPtr* jacobianPhi,
                            int* iparam, double* dparam)


#ifdef __cplusplus
}
#endif

#endif
