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
/*! \file NewtonEulerR.hpp

 */
#ifndef NEWTONEULERRELATION_H
#define NEWTONEULERRELATION_H

#include "Relation.hpp"
#include "PluginTypes.hpp"

/** NewtonEuler (Non Linear) Relation (generic interface)
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date Apr 27, 2004
 *
 * \class NewtonEulerR
 * Relations for NewtonEuler Dynamical Systems. This class is only an
 * interface for specific (Linear, Scleronomous ...)  NewtonEuler
 * Relations (see derived classes).
 *
 *  Class name = type+subType.
 *
 * If y = h(...), all the gradients of are handled by G object.
 * For example, G[0] = \f$ \nabla_q h(q,...) \f$.
 *
 * In corresponding derived classes, h and Gi are connected to plug-in functions (user-defined).
 * For more details, see the DevNotes.pdf, chapter NewtonEuler.
 */

class NewtonEulerR : public Relation
{
public:
// add deltaq ??? -- xhub 30/03/2014
  enum NewtonEulerRDS  {xfree, z, q0, velocity, dotq, p0, p1, p2, DSlinkSize};
  // enum NewtonEulerRVec {xfree, z, q0, dotq, p0, p1, p2, workVecSize};
  // enum NewtonEulerRMat {C, D, F, workMatSize};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEulerR);

  /** Jacobian matrices of H */
  /** Jacobian matrices of \f$y = h(t,q,\dot q,\ldots)\f$ */

  /**The Jacobian of the constraints with respect to the generalized coodinates  \f$q\f$
   *  i.e. \f[\nabla^T_q h(t,q,\dot q,\ldots)\f]
   */
  SP::SimpleMatrix _jachq;

  /**The Jacobian of the constraints with respect to the generalized velocities  \f$\dot q\f$
   *  i.e. \f[\nabla^T_{\dot q} h(t,q,\dot q,\ldots)\f]
   */
  SP::SimpleMatrix _jachqDot;

  /**The time-derivative of Jacobian of the constraints with respect
     to the generalized coordinates  \f$ q\f$
   * i.e. \f[\frac{d}{dt} \nabla^T_{q} h(t,q,\dot q,\ldots).\f]
   * This value is useful to compute the second-order
   * time--derivative of the constraints with respect to time.
   */
  SP::SimpleMatrix _dotjachq;

  SP::SimpleMatrix _jachlambda;
  SP::SimpleMatrix _jacglambda;

  /**vector e*/
  SP::SiconosVector _e;
  /*Used for the projection formulation*/

  /**vector of contact forces, ie: _contactForce = B lambda. Useful for the end user.*/
  SP::SiconosVector _contactForce;

  /**updated in computeJachqT:
  In the case of the bilateral constrains, it is _jachq._T.
  In the case of a local frame, _jachqT is built from the geometrical datas(local frame, point of contact).*/
  SP::SimpleMatrix _jachqT;

  /** local storage of _T as working vector to compute JachqT from q */
  SP::SimpleMatrix _T;

  /** basic constructor
  \param lagType the sub-type of the relation
  */
  NewtonEulerR(RELATION::SUBTYPES lagType): Relation(RELATION::NewtonEuler, lagType) {}



public:
  NewtonEulerR(): Relation(RELATION::NewtonEuler, RELATION::NonLinearR) {}

  /** destructor
  */
  virtual ~NewtonEulerR() {};

  // -- Jach --

  /** get a pointer on matrix Jach[index]
  *  \return a pointer on a SimpleMatrix
  */
  inline SP::SimpleMatrix jachq() const
  {
    return _jachq;
  }

  //proj_with_q  inline SP::SimpleMatrix jachqProj() const { return _jachqProj; }
  void setJachq(SP::SimpleMatrix newJachq);

  inline SP::SimpleMatrix jachqDot() const
  {
    return _jachqDot;
  }
  inline SP::SimpleMatrix dotJachq() const
  {
    assert(_dotjachq);
    return _dotjachq;
  }

  inline SP::SiconosVector secondOrderTimeDerivativeTerms()
  {
    assert(_secondOrderTimeDerivativeTerms);
    return _secondOrderTimeDerivativeTerms;
  };

  inline SP::SimpleMatrix jachlambda() const
  {
    return _jachlambda;
  }
  inline SP::SimpleMatrix jacglambda() const
  {
    return _jacglambda;
  }
  inline void setE(SP::SiconosVector newE)
  {
    _e = newE;
  }

  inline SP::SimpleMatrix jachqT() const
  {
    return _jachqT;
  }
  inline void setJachqT(SP::SimpleMatrix newJachqT)
  {
    _jachqT = newJachqT;
  }

  /** set Jach[index] to pointer newPtr (pointer link)
  *  \param newPtr the new matrix
  */
  void setJachqPtr(SP::SimpleMatrix newPtr);

  /** Plugin object for the time--derivative of Jacobian i.e.
  * \f[\frac{d}{dt} \nabla^T_{q} h(t,q,\dot q,\ldots).\f]
  * stored in _dotjachq
  */
  SP::PluggedObject _plugindotjacqh;

  /**  the additional  terms of the second order time derivative of y
   *
   *    \f$ \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$
   *
   */
  SP::SiconosVector _secondOrderTimeDerivativeTerms;

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction using this relation
   * \param DSlink the container of the link to DynamicalSystem attributes
  */
  void initializeDSLink(Interaction& inter, VectorOfBlockVectors& DSlink);

  /** initialize components specific to derived classes.
   * \param inter  Interaction associated with the Relation
   * \param DSlink
   * \param workV
   * \param workM
   */
  virtual void initializeWorkVectorsAndMatrices(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

  /** to compute y = h(q,v,t) using plug-in mechanism
  * \param time current time
  * \param q0 the position
  * \param y the output
  */
  virtual void computeh(double time, BlockVector& q0, SiconosVector& y);

  /** default function to compute jacobianH
   * \param time current time
   * \param inter the interaction using this relation
   */

  virtual void computeJachlambda(double time, Interaction& inter)
  {
    ;
  }
  /** compute the jacobian of h w.r.t. q
   * \param time current time
   * \param inter the interaction using this relation
   * \param q0  the container of the block vector to the dynamical system
   */
  virtual void computeJachq(double time, Interaction& inter, SP::BlockVector q0)
  {
    ;
  }

  /** compute the jacobian of h w.r.t. \f$\dot{q}\f$
   * \param time current time
   * \param inter the interaction using this relation
   */
  virtual void computeJachqDot(double time, Interaction& inter)
  {
    /* \warning. This method should never be called, since we are only considering
     * holonomic NewtonEulerR up to now
     */
    assert(0) ;
  }
  virtual void computeDotJachq(double time, BlockVector& workQ, BlockVector& workZ, BlockVector& workQdot);


  /** compute the jacobian of h w.r.t. \f$\dot{q}\f$
   * \param time current time
   * \param inter the interaction using this relation
   */
  virtual void computeJacglambda(double time, Interaction& inter)
  {
    ;
  }
  /** compute the jacobian of h w.r.t. \f$\dot{q}\f$
   * \param time current time
   * \param inter the interaction using this relation
   */
  virtual void computeJacgq(double time, Interaction& inter)
  {
    ;
  }
  /** compute the jacobian of h w.r.t. \f$\dot{q}\f$
   * \param time current time
   * \param inter the interaction using this relation
   */
  virtual void computeJacgqDot(double time, Interaction& inter)
  {
    ;
  }

  /* default implementation consists in multiplying jachq and T
   * in this implementation we use _T which is consitent which directly
   * computed with computeT(q) when q is given
   * this one in more consistent with the notion of function of q
   *
   *  \param inter interaction that owns the relation
   *  \param q0  the block vector to the dynamical system position
  */
  virtual void computeJachqT(Interaction& inter, SP::BlockVector q0);


  /** compute all the jacobian of h
   * \param time current time
   * \param inter the interaction using this relation
   * \param interProp Interaction properties
   */
  virtual void computeJach(double time, Interaction& inter, InteractionProperties& interProp);

  /** compute all the jacobian of g
   * \param time current time
   * \param inter the interaction using this relation
   * \param interProp Interaction properties
   */
  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp)
  {
    computeJacgq(time, inter);
    computeJacgqDot(time, inter);
    computeJacglambda(time, inter);
  }

  /** To compute the terms of the second order time derivative of y
      \f$ \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$
      \param time  current time
      \param inter interaction that owns the relation
      \param DSlink the container of the link to DynamicalSystem attributes
      \param ds1 dynamical system linked to this interaction (source)
      \param ds2 second ds linked to this interaction (target). If there is
      only one ds in the inter, call this function with ..., ds, ds)
   */
  void computeSecondOrderTimeDerivativeTerms(double time, Interaction& inter, VectorOfBlockVectors& DSlink, SP::DynamicalSystem ds1, SP::DynamicalSystem ds2);

  /** to compute output
   * \param time current time
   * \param inter the interaction using this relation
   * \param derivativeNumber number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double time, Interaction& inter, unsigned int derivativeNumber = 0);

  /** to compute the input
   * \param time current time
   * \param inter the interaction using this relation
   * \param level number of the derivative to compute, optional, default = 0.
   */
  virtual void computeInput(double time, Interaction& inter, unsigned int level = 0);

  /**
  * return a SP on the C matrix.
  * The matrix C in the linear case, else it returns Jacobian of the output with respect to x.
  * \return SP::SimpleMatrix
  */
  virtual inline SP::SimpleMatrix C() const
  {
    return _jachq;
  }
  /**
  * return a SP on the D matrix.
  * The matrix D in the linear case, else it returns Jacobian of the output with respect to lambda.
  * \return SP::SimpleMatrix
  */
  virtual inline SP::SimpleMatrix D() const
  {
    return _jachlambda;
  }
  /**
  * return a SP on the B matrix.
  * The matrix B in the linear case, else it returns Jacobian of the input with respect to lambda.
  * \return SP::SimpleMatrix
  */
  virtual inline SP::SimpleMatrix B() const
  {
    return _jacglambda;
  }
  /** A buffer containing the forces due to this.
  It is an output unused for the computation.
  Fix : is it usefull ?
  \return SP::SiconosVector
  */
  inline SP::SiconosVector contactForce() const
  {
    return _contactForce;
  };

  ACCEPT_STD_VISITORS();

};
#endif // NEWTONEULERRELATION_H
