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
/*! \file LagrangianCompliantR.hpp

 */
#ifndef LagrangianCompliantR_H
#define LagrangianCompliantR_H

#include "LagrangianR.hpp"

/** Lagrangian Compliant Relation: Scleronomous, Non-Holonomic (function of lambda).
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date Apr 27, 2004
 *
 * \f[
 * Y[0] = y = h(q,\lambda(t),z)
 * \f]
 *
 * \f[
 * Y[1] = \dot y = G0(q,\lambda(t),z)\dot q + G1((q,\lambda(t),z)\dot\lambda(t)
 * \f]
 *
 * \f[
 * p = G0^t(q,\lambda(t),z)\lambda(t)
 * \f]
 *
 * with
 * \f[
 * G0(q,\lambda(t),z) = \nabla_q h(q,\lambda(t),z)
 * \f]
 * \f[
 * G1(q,\lambda(t),z) = \nabla_{\lambda}h(q,\lambda(t),z)
 * \f]
 *
 * h, G0 and G1 are connected to user-defined functions.
 *
 */
class LagrangianCompliantR : public LagrangianR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianCompliantR);


  /** LagrangianR plug-in to compute h(q,lambda,z)
  * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
  * @param q : pointer to the first element of q
  * @param sizeY : size of vector y (ie of lambda and of the interaction)
  * @param lambda : pointer to lambda of the interaction
  * @param[in,out] y : pointer to the first element of y
  * @param sizeZ : size of vector z.
  * @param[in,out] z : a vector of user-defined parameters
  */
  SP::PluggedObject _pluginJachlambda;

  /** default constructor
  */
  LagrangianCompliantR() : LagrangianR(RELATION::CompliantR) {  };

  /** initialize G matrices or components specific to derived classes.
   * \param inter : the Interaction
   * \param DSlink : block vectors from dynamical systems
   * \param workV : work vectors
   * \param workM : work vectors
  */
  void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);
  void _zeroPlugin();

public:

  /** constructor from a set of data
  *  \param pluginh the name of the plugin to compute h
  *  \param pluginJacobianhq the name of the plugin to compute the gradient of h w.r.t q
  *  \param pluginJacobianhlambda the name of the plugin to compute the gradient of h w.r.t \f$\lambda\f$
  */
  LagrangianCompliantR(const std::string& pluginh, const std::string& pluginJacobianhq, const std::string& pluginJacobianhlambda);

  /** destructor
  */
  virtual ~LagrangianCompliantR() {};

  /** compute y = h(q,v,t) using plug-in mechanism
  * \param time the current time
  * \param q vector of coordinates
  * \param lambda vector for \f[ \lambda \f]
  * \param z parameter vector
  * \param y vector for y
  */
  virtual void computeh(double time, SiconosVector& q, SiconosVector& lambda, SiconosVector& z, SiconosVector& y);

  /** compute the jacobian of h w.r.t. q using plug-in mechanism
  * \param time current time
  * \param q vector of coordinates
  * \param lambda vector for \f[ \lambda \f]
  * \param z parameter vector
  */
  virtual void computeJachq(double time, SiconosVector& q, SiconosVector& lambda, SiconosVector& z);

  /** compute the jacobian of h w.r.t. \f$\lambda\f$ using plug-in mechanism
  * \param time  current time
  * \param q0 coordinates
  * \param lambda vector for \f[ \lambda \f]
  * \param z vector of parameters
  */
  virtual void computeJachlambda(double time, SiconosVector& q0, SiconosVector& lambda, SiconosVector& z);

  /** to compute output
  *  \param time the current time
  *  \param inter the Interaction owning y
  *  \param interProp Interaction properties
  *  \param derivativeNumber the number of the derivative to compute,
  *  optional, default = 0.
  */
  void computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber = 0);

  /** to compute the input
  *  \param time the current time
  *  \param inter the Interaction owning lambda
  *  \param interProp Interaction properties
  *  \param level "derivative" order of lambda used to compute input
  */
  void computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);

  /* compute all the H Jacobian */
  void computeJach(double time, Interaction& inter, InteractionProperties& interProp);

  /* compute all the G Jacobian */
  void computeJacg(double time, Interaction& inter, InteractionProperties& interProp)
  {
    ;
  }

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianCompliantR)

#endif
