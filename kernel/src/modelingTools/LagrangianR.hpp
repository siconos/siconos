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
/*! \file LagrangianR.hpp

 */
#ifndef LAGRANGIANRELATION_H
#define LAGRANGIANRELATION_H

#include "Relation.hpp"
#include "Interaction.hpp"
#include "PluginTypes.hpp"

/**
   Lagrangian Non Linear Relation (generic interface)
   
   This class is an interface for specific Lagrangian Relations used for Lagrangian dynamical systems.
   
   -  \f$ y = h(t,q,\dot q,\ldots) \f$  describes the constraint (the relation)
   
   - The Jacobian of the constraints with respect to the coodinates  \f$ q \f$ 
   i.e.  \f$ \nabla^T_q h(t,q,\dot q,\ldots) \f$ , is accessed with jachq().
   
   This Jacobian is mainly used for Newton linearization and to compute the time-derivative of the constraint,
   
   \f$ y = h(q,\ldots) \f$  that is  \f$ \dot y (t) = \nabla^T_q h(t,q,\dot q,\ldots) (q) \dot q +\ldots \f$  
   
   This object can also store more general linearized part of the gap function.
   If  \f$ y=h(q) \f$  models a gap function, then the time derivative
   can be generically  written as
   
   \f$ \dot y (t) = H(q,\ldots) \dot q  +\ldots.  \f$ 
   The matrix  \f$ H(q,\ldots) \f$ can also be accessed using jachq().
   
   - The Jacobian of the constraints with respect to the generalized velocities  \f$ \dot q \f$
   
   i.e.  \f$ \nabla^\top_{\dot q} h(t,q,\dot q,\ldots) \f$  is accessed using jachqDot().
   
   - The time-derivative of Jacobian of the constraints with respect to the generalized coordinates \f$ q \f$ 
   
   i.e.  \f$ \frac{d}{dt} \nabla^\top_{q} h(t,q,\dot q,\ldots). \f$ , is accessed using dotJachq().
   
   This value is useful to compute the second-order time--derivative of the constraints with respect to time.
   
   All these operators can be defined with user-defined plugins.

 */


class LagrangianR : public Relation
{
public:
  enum LagrangianRDS  {z, q0, q1, q2, p0, p1, p2, DSlinkSize};

protected:

  ACCEPT_SERIALIZATION(LagrangianR);

  /** Jacobian matrices of  \f$ y = h(t,q,\dot q,\ldots) \f$  */
  SP::SimpleMatrix _jachlambda{nullptr};

  /** The Jacobian of the constraints with respect to the generalized coodinates   \f$ q \f$ 
   *  i.e.  \f$ \nabla^\top_q h(t,q,\dot q,\ldots) \f$ 
   */
  SP::SimpleMatrix _jachq{nullptr};

  /**The Jacobian of the constraints with respect to the generalized velocities   \f$ \dot q \f$ 
   *  i.e.  \f$ \nabla^\top_{\dot q} h(t,q,\dot q,\ldots) \f$ 
   */
  SP::SimpleMatrix _jachqDot{nullptr};

  /** The time-derivative of Jacobian of the constraints with respect
   *  to the generalized coordinates   \f$  q \f$ 
   *  i.e.  \f$ \frac{d}{dt} \nabla^\top_{ q} h(t,q,\dot q,\ldots). \f$ 
   *  This value is useful to compute the second-order
   *  time--derivative of the constraints with respect to time.
   */
  SP::SimpleMatrix _dotjachq{nullptr};

  SP::PluggedObject _pluginJachq{nullptr};

  /** basic constructor
   *
   *  \param lagType the sub-type of the relation
   */
  LagrangianR(RELATION::SUBTYPES lagType): Relation(RELATION::Lagrangian, lagType) {}


  void _zeroPlugin() override;

public:

  /** destructor
   */
  virtual ~LagrangianR() noexcept = default;


  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction using this relation
   */
  inline void initialize(Interaction &inter) override {};

  // -- Jach --

  /** get a pointer on matrix Jach[index]
   *
   *  \return a pointer on a SimpleMatrix
   */
  inline SP::SimpleMatrix jachq() const
  {
    return _jachq;
  }
  inline SP::SimpleMatrix jachqDot() const
  {
    return _jachqDot;
  }
  inline SP::SimpleMatrix dotJachq() const
  {
    return _dotjachq;
  }
  inline SP::SimpleMatrix jachlambda() const
  {
    return _jachlambda;
  }

  /** set Jach[index] to pointer newPtr (pointer link)
   *
   *  \param newPtr the new matrix
   */
  inline void setJachqPtr(SP::SimpleMatrix newPtr)
  {
    _jachq = newPtr ;
  }

  inline SP::SimpleMatrix C() const override
  {
    return _jachq;
  }
  
  inline SP::SimpleMatrix H() const override
  {
    return _jachq;
  }

  /** main relation members display
  */
  void display() const override;

};
TYPEDEF_SPTR(LagrangianR)
#endif // LAGRANGIANRELATION_H
