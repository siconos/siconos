/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

    This class provides tools to describe non linear relation of the type:

    \rst
    .. math::


        y = h(q,t,z) \\
        \\dot y =  \nabla^\top_q(q,t,z)\\dot q + \frac{\partial }{\partial
   t}h(q,t,z) \\ \endrst

    or more generally

    \rst
    .. math::

        \\dot y =  H(q,t,z)\\dot q + \frac{\partial }{\partial t}h(q,t,z)

    \endrst

    and by duality

    \rst

    .. math::

        p = H^\top(q,t,z)\lambda

    \endrst

    The following operators (and their jacobians) can be plugged, in the usual
   way (see User Guide, 'User-defined plugins')

    - \f$ h(q,t,z)\f$
    - \f$ \nabla_q h(q,t,z)\f$
    - \f$ \dot h(q,t,z)\f$

    The plugin functions must fit with the following signature (FPtr4):

    void func(unsigned int qsize, double* q, double time, unsigned int ysize,
   double* buffer , unsigned int sizez, double* z)

    buffer being either \f$y\f$, \f$\dot h\f$ or \f$\nabla_qh\f$.
 */
class LagrangianRheonomousR : public LagrangianR {

protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(LagrangianRheonomousR);

  /** plugged vector used to compute hDot */
  SP::SiconosVector _hDot;

  /** LagrangianRheonomousR plug-in to compute hDot(q,t,z)
   */
  SP::PluggedObject _pluginhDot;

  /** default constructor
   */
  LagrangianRheonomousR() : LagrangianR(RELATION::RheonomousR)
  {
    _zeroPlugin();
  };
  void _zeroPlugin() override;

public:
  /** constructor from a set of data
   *  \param pluginh name of the plugin to compute h.
   * Its signature must be "void userPluginH(unsigned int, double*, double,
   * unsigned int, double*, unsigned int, double*)" \param pluginJacobianhq name
   * of the plugin  to compute jacobian h according to q. Its signature must be
   * "void userPluginG0(unsigned int, double*, double, unsigned int, double*,
   * unsigned int, double*)" \param pluginDoth name of the plugin to compute
   * hDot. Its signature must be "void userPluginHDot(unsigned int, double*,
   * double, unsigned int, double*, unsigned int, double*)
   */
  LagrangianRheonomousR(const std::string &pluginh,
                        const std::string &pluginJacobianhq,
                        const std::string &pluginDoth);

  /** destructor
   */
  virtual ~LagrangianRheonomousR(){};

  /** initialize G matrices or components specific to derived classes.
   * \param inter the Interaction
   */
  void initialize(Interaction &inter) override;

  /** check sizes of the relation specific operators.
   * \param inter an Interaction using this relation
   */
  void checkSize(Interaction &inter) override;

  // -- hDot --

  /** get a pointer on vector hDot
   *  \return a smart pointer on a SiconosVector
   */
  inline SP::SiconosVector hDot() const { return _hDot; }

  /** to set a specified function to compute function hDot
   *  \param pluginpath the complete path to the plugin
   *  \param name the name of the function to use in this plugin
   */
  void setComputehDotFunction(const std::string &pluginpath,
                              const std::string &name);

  /** to compute the output y = h(t,q,z) of the Relation
      \param time current time value
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
      \param y the resulting vector
  */
  virtual void computeh(double time, const BlockVector &q, BlockVector &z,
                        SiconosVector &y);

  /** to compute the time-derivative of the output y = h(t,q,z), saved in
     attribute _hDot (access: hDot()) \param time current time value \param q
     coordinates of the dynamical systems involved in the relation \param z user
     defined parameters (optional)
  */
  virtual void computehDot(double time, const BlockVector &q, BlockVector &z);

  /** to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())
      \param time current time value
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
  */
  virtual void computeJachq(double time, const BlockVector &q, BlockVector &z);

  /* compute all the H Jacobian */
  void computeJach(double time, Interaction &inter) override;
  /* compute all the G Jacobian */
  void computeJacg(double time, Interaction &inter) override {}

  /** to compute output
   * \param time current time
   * \param inter the Interaction
   *  \param derivativeNumber number of the derivative to compute, optional,
   * default = 0.
   */
  void computeOutput(double time, Interaction &inter,
                     unsigned int derivativeNumber = 0) override;

  /** to compute p
   * \param time current time
   * \param inter the Interaction
   * \param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction &inter,
                    unsigned int level = 0) override;

  ACCEPT_STD_VISITORS();
};

TYPEDEF_SPTR(LagrangianRheonomousR)

#endif
