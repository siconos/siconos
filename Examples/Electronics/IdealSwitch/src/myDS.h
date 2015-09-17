/* Siconos-Kernel , Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

/*! \file MyDSDS.h
  First Order Non Linear Dynamical Systems
*/

#ifndef MYDSDS_H
#define MYDSDS_H

#include "SiconosKernel.hpp"

/**  General First Order Non Linear Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) April 29, 2004
 *
 * This class defines and computes a generic n-dimensional
 * dynamical system of the form :
 * \f[
 * M \dot x = f(x,t,z) + r,
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state.
 *    - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.
 *  For example, z may be used to set some perturbation parameters, or to control the system (z will be set by some actuators) or anything else.
 *
 *  with \f$ f : R^{n} \times R  \mapsto  R^{n}   \f$ .
 *  and M a nXn matrix.
 *
 * By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
 * and the initial conditions are given by
 *  \f[
 *  x(t_0)=x_0
 * \f]
 * To define a boundary Value Problem, the pointer on  a BoundaryCondition must be set.
 *
 * \f$ f(x,t) \f$ is a plug-in function, and can be computed using computeF(t).
 * Its Jacobian according to x is denoted jacobianfx, and computed thanks to computeJacobianXF(t).
 * f and jacobianfx can be plugged to external functions thanks to setComputeFFunction/setComputeJacobianXFFunction.
 *
 * Right-hand side of the equation is computed thanks to computeRhs(t).
 *
 * \f[
 *    \dot x =  M^{-1}(f(x,t,z)+ r)
 * \f]
 *
 * Its Jacobian according to x is jacobianXRhs:
 *
 *  \f[
 *   jacobianXRhs = \nabla_xrhs(x,t,z) = M^{-1}\nabla_xf(x,t,z)
 *  \f]
 *
 * At the time:
 *  - M is considered to be constant. (ie no plug-in, no jacobian ...)
 *  - M is not allocated by default. The only way to use M is setM or setMPtr.
 *
 */
class MyDS : public FirstOrderNonLinearDS
{

public:
  /** default constructor
   * \param the type of the system
   */
  MyDS(SP::SiconosVector x0);


  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~MyDS() {};


  /** Default function to compute \f$ f: (x,t)\f$
   * \param double time : current time
   */
  virtual void computeF(double);

  /** function to compute \f$ f: (x,t)\f$ with x different from current saved state.
   * \param double time : current time
   * \param SP::SiconosVector
   */
  virtual void computeF(double, SP::SiconosVector);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeJacobianfx(double, bool  = false);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$ with x different from current saved state.
   *  \param double time : current time
   *  \param SP::SiconosVector
   */
  virtual void computeJacobianfx(double, const SiconosVector&);

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeRhs(double, bool  = false);
  virtual void resetNonSmoothPart(unsigned int level);

};

TYPEDEF_SPTR(MyDS);

#endif


