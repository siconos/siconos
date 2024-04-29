/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include <cassert>

#include "Relation.hpp"

namespace siconos::modeling {

class DynamicalSystem;

/**
    NewtonEuler (Non Linear) Relation (generic interface)

    Relations for NewtonEuler Dynamical Systems. This class is only an
    interface for specific (Linear, Scleronomous ...)  NewtonEuler
    Relations (see derived classes).

    Class name = type+subType.

    If y = h(...), all the gradients of are handled by G object.
    For example, G[0] = \f$ \nabla_q h(q,...) \f$.

    In corresponding derived classes, h and Gi are connected to plug-in functions
    (user-defined). For more details, see the DevNotes.pdf, chapter NewtonEuler.
*/

class NewtonEulerR : public Relation {
 public:
  // add deltaq ??? -- xhub 30/03/2014
  enum NewtonEulerRDS { z, q0, velocity, dotq, p0, p1, p2, DSlinkSize };

 protected:
  ACCEPT_SERIALIZATION(NewtonEulerR);

  /* Jacobian matrices of H */
  /* Jacobian matrices of  \f$ y = h(t,q,\dot q,\ldots) \f$  */

  /** The Jacobian of the constraints with respect to the generalized coodinates
   *  \f$ q \f$  i.e. \f[\nabla^T_q h(t,q,\dot q,\ldots)\f]
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _jachq{nullptr};

  /** The Jacobian of the constraints with respect to the generalized velocities
   *  \f$ \dot q \f$  i.e. \f[\nabla^T_{\dot q} h(t,q,\dot q,\ldots)\f]
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _jachqDot{nullptr};

  /** The time-derivative of Jacobian of the constraints with respect
   *  to the generalized coordinates  \f$ q \f$
   *  i.e. \f[\frac{d}{dt} \nabla^T_{q} h(t,q,\dot q,\ldots).\f]
   *  This value is useful to compute the second-order
   *  time--derivative of the constraints with respect to time.
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _dotjachq{nullptr};

  std::shared_ptr<siconos::algebra::SiconosMatrix> _jachlambda{nullptr};
  std::shared_ptr<siconos::algebra::SiconosMatrix> _jacglambda{nullptr};

  /** vector e*/
  std::shared_ptr<siconos::algebra::SiconosVector> _e{nullptr};
  /*Used for the projection formulation*/

  /** vector of contact forces, ie: _contactForce = B lambda. Useful for the end user.*/
  std::shared_ptr<siconos::algebra::SiconosVector> _contactForce{nullptr};

  /**
     updated in computeJachqT:
     In the case of the bilateral constrains, it is _jachq._T.
     In the case of a local frame, _jachqT is built from the geometrical
     datas(local frame, point of contact).*/
  std::shared_ptr<siconos::algebra::SiconosMatrix> _jachqT{nullptr};

  /** local storage of _T as working vector to compute JachqT from q */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _T{nullptr};

  // /** basic constructor
  //  *
  //  *  \param lagType the sub-type of the relation
  //  */
  // NewtonEulerR(RelationSubType lagType) : Relation(RelationType::NewtonEuler, lagType) {}

 public:
  /** Default constructor */
  NewtonEulerR() : Relation(RelationType::NewtonEuler, RelationSubType::NonLinearR){};

  /** destructor
   */
  virtual ~NewtonEulerR() noexcept = default;

  // -- Jach --

  /** get a pointer on matrix Jach[index]
   *
   *  \return a pointer on a SiconosMatrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> jachq() const { return _jachq; }

  // proj_with_q  inline std::shared_ptr<siconos::algebra::SiconosMatrix> jachqProj() const {
  // return _jachqProj;
  // }
  void setJachq(std::shared_ptr<siconos::algebra::SiconosMatrix> newJachq);

  inline std::shared_ptr<siconos::algebra::SiconosMatrix> jachqDot() const { return _jachqDot; }
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> dotJachq() const
  {
    assert(_dotjachq);
    return _dotjachq;
  }

  inline std::shared_ptr<siconos::algebra::SiconosVector> secondOrderTimeDerivativeTerms()
  {
    assert(_secondOrderTimeDerivativeTerms);
    return _secondOrderTimeDerivativeTerms;
  };

  inline std::shared_ptr<siconos::algebra::SiconosMatrix> jachlambda() const
  {
    return _jachlambda;
  }
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> jacglambda() const
  {
    return _jacglambda;
  }
  inline void setE(std::shared_ptr<siconos::algebra::SiconosVector> newE) { _e = newE; }

  inline std::shared_ptr<siconos::algebra::SiconosMatrix> jachqT() const { return _jachqT; }
  inline void setJachqT(std::shared_ptr<siconos::algebra::SiconosMatrix> newJachqT)
  {
    _jachqT = newJachqT;
  }
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> H() const override { return _jachqT; }

  /** set Jach[index] to pointer newPtr (pointer link)
   *
   *  \param newPtr the new matrix
   */
  void setJachqPtr(std::shared_ptr<siconos::algebra::SiconosMatrix> newPtr);

  /** Plugin object for the time--derivative of Jacobian i.e.
   *  \f$ \frac{d}{dt} \nabla^T_{q} h(t,q,\dot q,\ldots) \f$
   *  stored in _dotjachq
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _plugindotjacqh;

  /**  the additional  terms of the second order time derivative of y
   *
   *   \f$ \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$
   *
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _secondOrderTimeDerivativeTerms;

  /** initialize components specific to derived classes.
   *
   *  \param inter  Interaction associated with the Relation
   */
  virtual void initialize(Interaction &inter) override;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  virtual void checkSize(Interaction &inter) override;

  /**
      to compute the output y = h(t,q,z) of the Relation

      \param time current time value
      \param q coordinates of the dynamical systems involved in the relation
      \param y the resulting vector
  */
  virtual void computeh(double time, const siconos::algebra::BlockVector &q0,
                        siconos::algebra::SiconosVector &y);

  /** default function to compute jacobianH
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */

  virtual void computeJachlambda(double time, Interaction &inter) { ; }
  /** compute the jacobian of h w.r.t. q
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  the container of the block vector to the dynamical system
   */
  virtual void computeJachq(double time, Interaction &inter,
                            std::shared_ptr<siconos::algebra::BlockVector> q0)
  {
    ;
  }

  /** compute the jacobian of h w.r.t.  \f$ \dot{q} \f$
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJachqDot(double time, Interaction &inter)
  {
    /* \warning. This method should never be called, since we are only
     * considering holonomic NewtonEulerR up to now
     */
    assert(0);
  }
  virtual void computeDotJachq(double time, const siconos::algebra::BlockVector &workQ,
                               siconos::algebra::BlockVector &workZ,
                               const siconos::algebra::BlockVector &workQdot);

  /** compute the jacobian of h w.r.t.  \f$ \dot{q} \f$
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJacglambda(double time, Interaction &inter) { ; }
  /** compute the jacobian of h w.r.t.  \f$ \dot{q} \f$
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJacgq(double time, Interaction &inter) { ; }
  /** compute the jacobian of h w.r.t.  \f$ \dot{q} \f$
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJacgqDot(double time, Interaction &inter) { ; }

  /** default implementation consists in multiplying jachq and T
   *  in this implementation we use _T which is consitent which directly
   *  computed with computeT(q) when q is given
   *  this one in more consistent with the notion of function of q
   *
   *  \param inter interaction that owns the relation
   *  \param q0  the block vector to the dynamical system position
   */
  virtual void computeJachqT(Interaction &inter,
                             std::shared_ptr<siconos::algebra::BlockVector> q0);

  /** compute all the jacobian of h
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJach(double time, Interaction &inter) override;

  /** compute all the jacobian of g
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJacg(double time, Interaction &inter) override
  {
    computeJacgq(time, inter);
    computeJacgqDot(time, inter);
    computeJacglambda(time, inter);
  }

  /**
      To compute the terms of the second order time derivative of y
      \f$  \nabla_q h(q) \dot T v + \frac{d}{dt}(\nabla_q h(q) ) T v \f$

      \param time  current time
      \param inter interaction that owns the relation
      \param ds1 dynamical system linked to this interaction (source)
      \param ds2 second ds linked to this interaction (target). If there is
      only one ds in the inter, call this function with ..., ds, ds)
   */
  void computeSecondOrderTimeDerivativeTerms(double time, Interaction &inter,
                                             std::shared_ptr<DynamicalSystem> ds1,
                                             std::shared_ptr<DynamicalSystem> ds2);

  /** to compute output
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param derivativeNumber number of the derivative to compute, optional,
   *  default = 0.
   */
  virtual void computeOutput(double time, Interaction &inter,
                     unsigned int derivativeNumber = 0) override;

  /** to compute the input
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param level number of the derivative to compute, optional, default = 0.
   */
  virtual void computeInput(double time, Interaction &inter, unsigned int level = 0) override;

  /** return a SP on the C matrix.
   *  The matrix C in the linear case, else it returns Jacobian of the output
   *  with respect to x.
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> C() const override { return _jachq; }
  /** return a SP on the D matrix.
   *  The matrix D in the linear case, else it returns Jacobian of the output
   *  with respect to lambda.
   */
  virtual inline std::shared_ptr<siconos::algebra::SiconosMatrix> D() const
  {
    return _jachlambda;
  }
  /** return a SP on the B matrix.
   *  The matrix B in the linear case, else it returns Jacobian of the input with
   *  respect to lambda.
   */
  virtual inline std::shared_ptr<siconos::algebra::SiconosMatrix> B() const
  {
    return _jacglambda;
  }

  /**
      A buffer containing the forces due to this.
      It is an output unused for the computation.
      Fix : is it usefull ?

      \return std::shared_ptr<siconos::algebra::SiconosVector>
  */
  inline std::shared_ptr<siconos::algebra::SiconosVector> contactForce() const
  {
    return _contactForce;
  };

  void display() const override {}
 
  // Visitors stuff
  void accept(std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const override;

  
};
}  // namespace siconos::modeling
#endif  // NEWTONEULERRELATION_H
