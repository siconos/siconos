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
/*! \file LagrangianScleronomousR.hpp

 */
#ifndef LagrangianScleronomousR_H
#define LagrangianScleronomousR_H

#include "LagrangianR.hpp"
/** \class LagrangianScleronomousR
 * \brief   Scleronomic Lagrangian (Non Linear) Relations
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date Apr 27, 2004
 *
 * Scleronomic Relation (constraint) :
 * \f[
 * y = h(q,z) \\
 * \f]
 *
 * \f[
 * \dot y = \nabla^\top_q h(q,z) \dot q
 * \f]
 * or more generally
 * \f[
 * \dot y = H(q,z) \dot q
 * \f]
 *
 * and by duality
 *
 * \f[
 * p = \nabla_q h(q,z)\lambda
 * \f]
 * or more generally
 * \f[
 * p = H^\top(q,z)\lambda
 * \f]
 *
 * with
 *
 * \f[
 * H^\top(q,z) = \nabla_q h(q,z)
 * \f]
 *
 * is the pure Lagrangian setting.
 *
 *  y (or its discrete approximation) is stored in y[0]
 * \f$ \dot y \f$ (or its discrete approximation) is  stored in y[1]
 *  higher level y[i] can be used for storing higher levels of derivatives.
 *
 * Jacobians and h are connected to plug-in functions.\n
 * The plugin function to compute h(q,z) needs the following parameters:\n
 * --> sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in the interaction\n
 * --> q : pointer to the first element of q \n
 * --> sizeY : size of vector y (ie of the interaction) \n
 * --> [in,out] y : pointer to the first element of y \n
 * --> sizeZ : size of vector z \n
 * --> [in,out] z: pointer to z vector(s) from DS. \n
 * Its signature must be "void plugin(unsigned int, double*, unsigned int, double*, unsigned int, double*)"\n\n
 * The plugin function to compute G0(q,z), gradient of h according to q, needs the following parameters: \n
 *--> sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in the interaction  \n
 *--> q : pointer to the first element of q  \n
 *--> sizeY : size of vector y (ie of the intercation) \n
 *--> [in,out] H : pointer to the first element of H (sizeY X sizeDS matrix)\n
 * --> sizeZ : size of vector z \n
 * -->[in,out] z: pointer to z vector(s) from DS.\n
 * Its signature must be "void plugin(unsigned int, double*, unsigned int, double*, unsigned int, double*)"\n
 *
 */
class LagrangianScleronomousR : public LagrangianR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianScleronomousR);

  /** LagrangianScleronomousR plug-in to compute G0(q,z), gradient of h according to q
  * @param sizeQ size of q = sum of the sizes of all the DynamicalSystems involved in the interaction
  * @param q  pointer to the first element of q
  * @param sizeY  size of vector y (ie of the intercation)
  * @param[in,out] G0 : pointer to the first element of G0 (sizeY X sizeDS matrix)
  * @param sizeZ size of vector z
  * @param[in,out] z: pointer to z vector(s) from DS.
  */
  /** Plugin object for the time--derivative of Jacobian i.e.
   * \f[\frac{d}{dt} \nabla^T_{q} h(t,q,\dot q,\ldots).\f]
   * stored in _dotjachq
   */
  SP::PluggedObject _plugindotjacqh;

  /** Product of the time--derivative of Jacobian with the velocity qdot */
  SP::SiconosVector _dotjacqhXqdot;

  /** reset all plugins */
  virtual void _zeroPlugin();

  /** basic constructor */
  LagrangianScleronomousR(): LagrangianR(RELATION::ScleronomousR)
  {
    _zeroPlugin();
  }

  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

public:

  /** constructor from a set of data
  *  \param pluginh the name of the plugin to compute h(q,z).
  * The signature  of the plugged function must be:
  *  "void pluginH(unsigned int, double*, unsigned int, double*, unsigned int, double*)"
  *  \param pluginJacobianhq the name of the plugin to compute jacobian h according to q.\n
  * The signature  of the plugged function must be:
  *  "void pluginG0(unsigned int, double*, unsigned int, double*, unsigned int, double*)"
  *
  */
  LagrangianScleronomousR(const std::string& pluginh, const std::string& pluginJacobianhq);

  /** constructor from a set of data used for EventDriven Scheme
  *  \param pluginh the name of the plugin to compute h(q,z).
  * The signature  of the plugged function must be:
  *  "void pluginH(unsigned int, double*, unsigned int, double*, unsigned int, double*)"
  *  \param pluginJacobianhq the name of the plugin to compute jacobian h according to q.\n
  * The signature  of the plugged function must be:
  *  "void pluginG0(unsigned int, double*, unsigned int, double*, unsigned int, double*)"
  * \param pluginDotJacobianhq the name of the plugin to compute the derivative of H Jacobian with respect to time
  * The signature of the plugged function must be:
  * "void pluginS0(unsigned int, double*,unsigned int, double*, unsigned int, double*, unsigned int, double*)"
  *
  */
  LagrangianScleronomousR(const std::string& pluginh, const std::string& pluginJacobianhq, const std::string& pluginDotJacobianhq);

  /** destructor
  */
  virtual ~LagrangianScleronomousR() {};

  /** \return the product of  the time--derivative of Jacobian with the velocity qdot */
  inline SP::SiconosVector dotjacqhXqdot()
  {
    return _dotjacqhXqdot;
  };

  /** to compute y = h(q,z) using plug-in mechanism
   * \param q the BlockVector of coordinates
   * \param z the BlockVector of parameters
   * \param y the output
   */
  virtual void computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y);

  /** to compute the jacobian of h using plug-in mechanism.
   * Index shows which jacobian is computed
   * \param q the BlockVector of coordinates
   * \param z the BlockVector of parameters
   */
  virtual void computeJachq(SiconosVector& q, SiconosVector& z);

  /** to compute the product of  the time--derivative of Jacobian with the velocity qdot
   * \param time double, current time
   * \param inter interaction that owns the relation
   * \param DSlink
   */
  void computedotjacqhXqdot(double time, Interaction& inter, VectorOfBlockVectors& DSlink);

  /* compute all the H Jacobian
   * \param time double, current time
   * \param inter interaction that owns the relation
   * \param interProp
   */
  void computeJach(double time, Interaction& inter, InteractionProperties& interProp);

  /* compute all the G Jacobian
   * \param time double, current time
   * \param inter interaction that owns the relation
   * \param interProp
   */
  void computeJacg(double time, Interaction& inter, InteractionProperties& interProp)
  {
    ;
  }

  /** to compute the time derivative of the Jacobian with respect to time using plug-in mechanism
   * \param q the BlockVector of coordinates
   * \param z the BlockVector of parameters
   * \param qDot q the BlockVector of derivative of coordinates
   */
  virtual void computeDotJachq(SiconosVector& q, SiconosVector& z, SiconosVector& qDot);

  /** to compute output
   * \param time the current time
   * \param inter interaction that owns the relation
   * \param interProp the InteractionProperties of this Interaction
   * \param derivativeNumber number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double time, Interaction& inter, InteractionProperties& interProp,
                             unsigned int derivativeNumber = 0);

  /** to compute p
   * \param time the current time
   * \param inter interaction that owns the relation
   * \param interProp the InteractionProperties of this Interaction
   * \param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction& inter, InteractionProperties& interProp,
                    unsigned int level = 0);


  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianScleronomousR)

#endif // LAGRANGIANRELATION_H
