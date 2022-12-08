/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include "SimpleMatrixFriends.hpp"
/** 
    Scleronomic Lagrangian (Non Linear) Relations

    \f[
    y = h(q,z)
    \f] 

    \f[
    \dot y = \nabla^\top_q h(q,z) \dot q
    \f]

    or more generally

    \f[
     \dot y = H(q,z) \dot q
    \f]

    and by duality

    \f[
    p = \nabla_q h(q,z)\lambda
    \f]
    
    or more generally

    \f[
    p = H^\top(q,z)\lambda
    \f]

    with

    \f[
    H^\top(q,z) = \nabla_q h(q,z)
    \f]
    
    is the pure Lagrangian setting.
    
    y (or its discrete approximation) is stored in y[0]
    \f$ \dot y \f$ (or its discrete approximation) is  stored in y[1]
    higher level y[i] can be used for storing higher levels of derivatives.

    Jacobians and h are connected to plug-in functions.
    
    The plugin function to compute h(q,z) needs the following parameters:

    --> sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved
    in the interaction

    --> q : pointer to the first element of q

    --> sizeY : size of vector y (ie of the interaction)

    --> [in,out] y : pointer to the first element of y

    --> sizeZ : size of vector z

    --> [in,out] z: pointer to z vector(s) from DS.

    Its signature must be "void plugin(unsigned int, double*, unsigned int,
    double*, unsigned int, double*)"
    
    The plugin function to compute G0(q,z),
    gradient of h according to q, needs the following parameters:

    --> sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in
    the interaction

    --> q : pointer to the first element of q

    --> sizeY : size of vector y (ie of the intercation)

    --> [in,out] H : pointer to the first element of H (sizeY X sizeDS matrix)

    --> sizeZ : size of vector z

    -->[in,out] z: pointer to z vector(s) from DS.

    Its signature must be "void plugin(unsigned int, double*, unsigned int,
    double*, unsigned int, double*)"

 */
class LagrangianScleronomousR : public LagrangianR {

protected:
  
  ACCEPT_SERIALIZATION(LagrangianScleronomousR);

  /* LagrangianScleronomousR plug-in to compute G0(q,z), gradient of h
   *  according to q
   *
   *  @param sizeQ size of q = sum of the sizes of all the DynamicalSystems
   *  involved in the interaction
   *  @param q  pointer to the first element of q
   *  @param sizeY  size of vector y (ie of the intercation)
   *  @param[in,out] G0 : pointer to the first element of G0 (sizeY X sizeDS
   *  matrix)
   *  @param sizeZ size of vector z
   *  @param[in,out] z: pointer to z vector(s) from DS.
   */
  
  /** Plugin object for the time--derivative of Jacobian i.e.
   *  \f$ \frac{d}{dt} \nabla^T_{q} h(t,q,\dot q,\ldots). \f$ 
   * stored in _dotjachq
   */
  SP::PluggedObject _plugindotjacqh{nullptr};

  /** Product of the time--derivative of Jacobian with the velocity qdot */
  SP::SiconosVector _dotjacqhXqdot{nullptr};

  /** reset all plugins */
  void _zeroPlugin() override;

  /** basic constructor */
  LagrangianScleronomousR() : LagrangianR(RELATION::ScleronomousR) {
    _zeroPlugin();
  }

public:
  /** constructor from a set of data
   *
   *  \param pluginh the name of the plugin to compute h(q,z).
   *  The signature  of the plugged function must be:
   *  "void pluginH(unsigned int, double*, unsigned int, double*, unsigned int,
   *  double*)" 
   *  \param pluginJacobianhq the name of the plugin to compute
   *  jacobian h according to q. The signature  of the plugged function must
   *  be: "void pluginG0(unsigned int, double*, unsigned int, double*, unsigned
   *  int, double*)"
   *
   */
  LagrangianScleronomousR(const std::string &pluginh,
                          const std::string &pluginJacobianhq);

  /** constructor from a set of data used for EventDriven Scheme
   *
   *  \param pluginh the name of the plugin to compute h(q,z).
   *  The signature  of the plugged function must be:
   *  "void pluginH(unsigned int, double*, unsigned int, double*, unsigned int,
   *  double*)" 
   *  \param pluginJacobianhq the name of the plugin to compute
   *  jacobian h according to q. The signature  of the plugged function must
   *  be: "void pluginG0(unsigned int, double*, unsigned int, double*, unsigned
   *  int, double*)" 
   *  \param pluginDotJacobianhq the name of the plugin to compute
   *  the derivative of H Jacobian with respect to time The signature of the
   *  plugged function must be: "void pluginS0(unsigned int, double*,unsigned
   *  int, double*, unsigned int, double*, unsigned int, double*)"
   *
   */
  LagrangianScleronomousR(const std::string &pluginh,
                          const std::string &pluginJacobianhq,
                          const std::string &pluginDotJacobianhq);

  /** destructor
   */
  virtual ~LagrangianScleronomousR() noexcept = default;

  void initialize(Interaction &inter) override;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  void checkSize(Interaction &inter) override;

  /** \return the product of  the time--derivative of Jacobian with the velocity
   * qdot */
  inline SP::SiconosVector dotjacqhXqdot() { return _dotjacqhXqdot; };

  /** 
      to compute the output y = h(q,z) of the Relation
      
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
      \param y the resulting vector
  */
  virtual void computeh(const BlockVector &q, BlockVector &z, SiconosVector &y);

  /** 
      to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())
      
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
  */
  virtual void computeJachq(const BlockVector &q, BlockVector &z);

  /**
     to compute the time derivative of the Jacobian. Result in _dotjachq
     (access: dotjachq()) 
     
     \param q coordinates of the dynamical systems involved in the relation
     \param z user defined parameters (optional)
     \param time derivatives of q
  */
  virtual void computeDotJachq(const BlockVector &q, BlockVector &z,
                               const BlockVector &qDot);

  /** to compute the product of  the time--derivative of Jacobian with the
   *  velocity qdot
   *
   *  \param time double, current time
   *  \param inter interaction 
   *  \param DSlink
   */
  void computedotjacqhXqdot(double time, Interaction &inter,
                            VectorOfBlockVectors &DSlink);

  /** compute all the H Jacobian
   *
   *  \param time double, current time
   *  \param inter interaction that owns the relation
   *  \param interProp
   */
  void computeJach(double time, Interaction &inter) override;

  /** compute all the G Jacobian
   *  
   *  \param time double, current time
   *  \param inter interaction that owns the relation
   *  \param interProp
   */
  void computeJacg(double time, Interaction &inter) override {}

  /** to compute output
   * 
   *  \param time the current time
   *  \param inter interaction that owns the relation
   *  \param derivativeNumber number of the derivative to compute, optional,
   *  default = 0.
   */
  void computeOutput(double time, Interaction &inter,
                     unsigned int derivativeNumber = 0) override;

  /** to compute p
   * 
   *  \param time the current time
   *  \param inter interaction that owns the relation
   *  \param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction &inter,
                    unsigned int level = 0) override;

  ACCEPT_STD_VISITORS();
};

TYPEDEF_SPTR(LagrangianScleronomousR)

#endif // LAGRANGIANRELATION_H
