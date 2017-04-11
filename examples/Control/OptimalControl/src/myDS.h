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

protected :
  SimpleMatrix * Q;
  SimpleMatrix * K1;
  //SimpleMatrix * K1T;

public:
  /** default constructor
   * \param the type of the system
   */
  MyDS(SP::SiconosVector x0);


  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~MyDS();

  /** function to compute \f$ f: (x,t)\f$ with x different from current saved state.
   * \param double time : current time
   * \param SP::SiconosVector
   */
  virtual void computef(double, SP::SiconosVector);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$ with x different from current saved state.
   *  \param double time : current time
   *  \param SP::SiconosVector
   */
  virtual void computeJacobianfx(double, SP::SiconosVector v);

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */



  void alpha(double t, SP::SiconosVector xvalue, SP::SiconosVector alpha);
  void JacobianXalpha(double t, SP::SiconosVector xvalue, SP::SiconosMatrix JacXalpha);



};

TYPEDEF_SPTR(MyDS);

#endif


