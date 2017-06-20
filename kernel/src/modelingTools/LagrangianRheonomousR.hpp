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
/*! \file LagrangianRheonomousR.hpp

 */
#ifndef LagrangianRheonomousR_H
#define LagrangianRheonomousR_H

#include "LagrangianR.hpp"

/** Lagrangian (Non Linear) Rheonomous Relation 
    
    \author SICONOS Development Team - copyright INRIA
    \date February 28, 2007
    
    This class provides tools to describe non linear relation of the type:
    
    \f[
    y &= h(q,t,z) \\
    \dot y &=  \nabla^\top_q(q,t,z)\dot q + \frac{\partial }{\partial t}h(q,t,z) \\
    \f]

    or more generally
    
    \f[
    \dot y =  H(q,t,z)\dot q + \frac{\partial }{\partial t}h(q,t,z)
    \f]
    
    and by duality

    \f[
    p = H^\top(q,t,z)\lambda
    \f]
    
    The following operators (and their jacobians) can be plugged, in the usual way (see User Guide, 'User-defined plugins')
    
    - \f$ h(q,t,z)\f$
    - \f$ \nabla_q h(q,t,z)\f$
    - \f$ \dot h(q,t,z)\f$

    The plugin functions must fit with the following signature (FPtr4):
    
    void func(unsigned int qsize, double* q, double time, unsigned int ysize, double* buffer , unsigned int sizez, double* z)
    
    buffer being either \f$y\f$, \f$\dot h\f$ or \f$\nabla_qh\f$.
 */
class LagrangianRheonomousR : public LagrangianR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianRheonomousR);

  /** plugged vector used to compute hDot */
  SP::SiconosVector _hDot;

  /** LagrangianRheonomousR plug-in to compute hDot(q,t,z)
  * @param sizeDS sum of the sizes of all the DynamicalSystems involved in the interaction
  * @param q pointer to the first element of q
  * @param time current time
  * @param sizeY size of vector hDot (ie of the intercation)
  * @param[in,out] pointer to the first element of hDot
  * @param sizeZ size of vector z
  * @param[in,out] z a vector of user-defined parameters
  */
  SP::PluggedObject _pluginhDot;

  /** initialize G matrices or components specific to derived classes.
   * \param inter the Interaction
   * \param DSlink block vectors from dynamical systems
   * \param workV work vectors
   * \param workM work vectors
   */
  void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

  /** default constructor
  */
  LagrangianRheonomousR(): LagrangianR(RELATION::RheonomousR)
  {
    _zeroPlugin();
  };
  void _zeroPlugin();
public:

  /** constructor from a set of data
  *  \param pluginh name of the plugin to compute h.\n
  * Its signature must be "void userPluginH(unsigned int, double*, double, unsigned int, double*, unsigned int, double*)"
  *  \param pluginJacobianhq name of the plugin  to compute jacobian h according to q.\n
  * Its signature must be "void userPluginG0(unsigned int, double*, double, unsigned int, double*, unsigned int, double*)"
  *  \param pluginDoth name of the plugin to compute hDot. \n
  * Its signature must be "void userPluginHDot(unsigned int, double*, double, unsigned int, double*, unsigned int, double*)
  */
  LagrangianRheonomousR(const std::string& pluginh, const std::string& pluginJacobianhq, const std::string& pluginDoth);

  /** destructor
  */
  virtual ~LagrangianRheonomousR() {};

  // -- hDot --

  /** get vector hDot
  *  \return a SiconosVector

  inline const SiconosVector gethDot() const { return *hDot; }
  */

  /** get a pointer on vector hDot
  *  \return a smart pointer on a SiconosVector
  */
  inline SP::SiconosVector hDot() const
  {
    return _hDot;
  }

  /** to set a specified function to compute function hDot
  *  \param pluginpath the complete path to the plugin
  *  \param name the name of the function to use in this plugin
  */
  void setComputehDotFunction(const std::string& pluginpath, const std::string& name);

  /** to compute y = h(t,q,z) using plug-in mechanism
  * \param time current time
  * \param q the vector of coordinates
  * \param z the vector of parameters
  * \param y the y vector
  */
  virtual void computeh(double time, SiconosVector& q, SiconosVector& z, SiconosVector& y);

  /** to compute hDot using plug-in mechanism
   * \param time current time
   * \param q the vector of coordinates
   * \param z the vector of parameters
   */
  virtual void computehDot(double time, SiconosVector& q, SiconosVector& z);

  /** to compute the jacobian of h using plug-in mechanism. Index shows which jacobian is computed
  * \param time double, current time
  * \param q the coordinates vector
  * \param z the parameters vector
  */
  virtual void computeJachq(double time, SiconosVector& q, SiconosVector& z);


  /* compute all the H Jacobian */
  void computeJach(double time, Interaction& inter, InteractionProperties& interProp);
  /* compute all the G Jacobian */
  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp)
  {
    ;
  }


  /** to compute output
   * \param time current time
   * \param inter the Interaction
   * \param interProp the Interaction properties
  *  \param derivativeNumber number of the derivative to compute, optional, default = 0.
  */
  virtual void computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber = 0);

  /** to compute p
   * \param time current time
   * \param inter the Interaction
   * \param interProp the Interactions properties
   * \param level "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianRheonomousR)

#endif
