/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
  */
  void initComponents(Interaction& inter);
  void zeroPlugin();

public:

  /** constructor from a set of data
  *  \param pluginh the name of the plugin to compute h
  *  \param computeJacobianhq the name of the plugin to compute the gradient of h w.r.t q
  *  \param computeJacobianhlambda the name of the plugin to compute the gradient of h w.r.t \f$\lambda\f$
  */
  LagrangianCompliantR(const std::string& pluginh, const std::string& pluginJacobianhq, const std::string& pluginJacobianhlambda);

  /** destructor
  */
  virtual ~LagrangianCompliantR() {};

  /** compute y = h(q,v,t) using plug-in mechanism
  * \param time current time
  * \param inter the Interaction
  */
  void computeh(double time, Interaction& inter);

  /** compute the jacobian of h w.r.t. q using plug-in mechanism
  * \param time current time
  * \param inter the Interaction
  */
  void computeJachq(double time, Interaction& inter);

  /** compute the jacobian of h w.r.t. \f$\lambda\f$ using plug-in mechanism
  * \param time current time
  * \param inter the Interaction
  */
  void computeJachlambda(double time, Interaction& inter);

  const std::string getJachlambdaName() const;
  const std::string getJachqName() const;


  /** to compute output
  *  \param time the current time
  *  \param inter the Interaction owning y
  *  \param level number of the derivative to compute, optional, default = 0.
  */
  void computeOutput(double time, Interaction& inter, unsigned int level = 0);

  /** to compute the input
  *  \param time the current time
  *  \param inter the Interaction owning lambda
  *  \param level "derivative" order of lambda used to compute input
  */
  void computeInput(double time, Interaction& inter, unsigned int level = 0);

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianCompliantR)

#endif
