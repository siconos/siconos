/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/** Lagrangian (Non Linear) Relation (generic interface)
 *
 * Relations for Lagrangian Dynamical Systems.
 * This class is only an interface for specific (Linear, scleronomic, rheomomic  ...)
 * Lagrangian Relations (see derived classes).
 *
 *  Class name = type+subType.
 *
 * If \f$y = h(t,q,\dot q,\ldots)\f$ describes the constraint (the relation) , all the gradients of h
 * are handled by the following SimpleMatrix and SiconosVector objects.
 *
 * - The Jacobian of the constraints with respect to the coodinates  \f$q\f$
 * 
 *  i.e. \f$\nabla^T_q h(t,q,\dot q,\ldots)\f$  is stored in  SP::SimpleMatrix _jachq .
 *
 * This Jacobian is mainly used for Newton linearization and to compute the time-derivative of the constraint \f$y = h(q,\ldots)\f$ that is
 *  \f$\dot y (t) = \nabla^T_q h(t,q,\dot q,\ldots) (q) \dot q +\ldots\f$
 * This object can also store
 * more general linearized part of the gap function. If \f$y=h(q)\f$ models a gap function, then the time--derivative
 * can be generically  written as
 * \f$\dot y (t) = H(q,\ldots) \dot q  +\ldots. \f$
 * The matrix \f$H(q,\ldots) \f$ is also stored in   SP::SimpleMatrix _jachq </li>
 *
 * - The Jacobian of the constraints with respect to the generalized velocities  \f$\dot q\f$
 *  i.e. \f$\nabla^\top_{\dot q} h(t,q,\dot q,\ldots)\f$ is stored in  SP::SimpleMatrix _jachqDot </li>
 *
 * - The time-derivative of Jacobian of the constraints with respect to the generalized coordinates  \f$ q\f$
 *  i.e. \f$\frac{d}{dt} \nabla^\top_{q} h(t,q,\dot q,\ldots).\f$. This value is useful to compute the second-order
 * time--derivative of the constraints with respect to time.</li>
 *
 * 
 *
 * In corresponding derived classes, h and Jacobians are connected to plug-in functions (user-defined).
 *
 */


class LagrangianR : public Relation
{
public:
  enum LagrangianRDS  {z, q0, q1, q2, p0, p1, p2, DSlinkSize};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianR);

  /** Jacobian matrices of \f$y = h(t,q,\dot q,\ldots)\f$ */

  SP::SimpleMatrix _jachlambda;

  /**The Jacobian of the constraints with respect to the generalized coodinates  \f$q\f$
   *  i.e. \f$\nabla^\top_q h(t,q,\dot q,\ldots)\f$
   */
  SP::SimpleMatrix _jachq;

  /**The Jacobian of the constraints with respect to the generalized velocities  \f$\dot q\f$
   *  i.e. \f$\nabla^\top_{\dot q} h(t,q,\dot q,\ldots)\f$
   */
  SP::SimpleMatrix _jachqDot;

  /**The time-derivative of Jacobian of the constraints with respect
     to the generalized coordinates  \f$ q\f$
   * i.e. \f$\frac{d}{dt} \nabla^\top_{ q} h(t,q,\dot q,\ldots).\f$
   * This value is useful to compute the second-order
   * time--derivative of the constraints with respect to time.
   */
  SP::SimpleMatrix _dotjachq;

  SP::PluggedObject _pluginJachq;

  /** basic constructor
   * \param lagType the sub-type of the relation
   */
  LagrangianR(RELATION::SUBTYPES lagType): Relation(RELATION::Lagrangian, lagType) {}


  virtual void _zeroPlugin();

public:

  /** destructor
  */
  virtual ~LagrangianR() {};

  // -- Jach --

  /** get a pointer on matrix Jach[index]
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
  *  \param newPtr the new matrix
  */
  inline void setJachqPtr(SP::SimpleMatrix newPtr)
  {
    _jachq = newPtr ;
  }

  inline SP::SimpleMatrix C() const
  {
    return _jachq;
  }

  /** initialize components specific to derived classes.
   * \param inter the interaction using this relation
   */
  virtual void initialize(Interaction& inter) {};
  
  /** check sizes of the relation specific operators.
   * \param inter an Interaction using this relation
   */
  virtual void checkSize(Interaction& inter) = 0;

  /* compute all the H Jacobian 
   * \param time
   * \param inter
   */
  virtual void computeJach(double time, Interaction& inter) = 0 ;
  /* compute all the G Jacobian
   * \param time
   * \param inter
   */
  virtual void computeJacg(double time, Interaction& inter) = 0 ;

  /** main relation members display
  */
  void display() const;

};
TYPEDEF_SPTR(LagrangianR)
#endif // LAGRANGIANRELATION_H
