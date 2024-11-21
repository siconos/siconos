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
/*! \file LagrangianR.hpp

 */
#ifndef LAGRANGIANRELATION_H
#define LAGRANGIANRELATION_H

#include "Relation.hpp"

namespace siconos::modeling {
/**
   Lagrangian Non Linear Relation (generic interface)

   This class is the interface to relations used for Lagrangian (2nd order) systems.

   -  \f$ y = h(t,q,\dot q,\ldots) \f$  describes the constraint

    h may be \f$ h(q)\f$  (scleronomous), \f$ h(q,t)\f$  (rheonomous), \f$ h(q, \lambda)\f$
  (compliant)... depending on the chosen derived class.

  The Jacobian of the constraints with respect to the coodinates  \f$ q \f$
  i.e.  \f$ \nabla^T_q h(t,q,\dot q,\ldots) \f$ is always defined.
  Other jacobians are defined only when required in the derived classes, check the API if
  required.


  This Jacobians are mainly used for Newton linearization and to compute the time-derivative of
  the constraint, \f$ y = h(q,\ldots) \f$

  that is  \f$ \dot y (t) = \nabla^T_q h(t,q,\dot q,\ldots) (q) \dot q +\ldots \f$

  This object can also store more general linearized part of the gap function.
  If \f$ y=h(q) \f$ models a gap function, then the time derivative can be generically  written
  as

  \f$ \dot y (t) = H(q,\ldots) \dot q  +\ldots.  \f$

  All operators of Lagrangian relations can be plugged to user-defined functions.
  The signature of the plugged functions depends on which variables are taken into
  account in h and its jacobians. Here again, check in the proper derived class API to find
  which functions are available.
 */

class LagrangianR : public Relation {
 public:
  enum LagrangianRDS { z, q0, q1, q2, p0, p1, p2, DSlinkSize };

 protected:
  ACCEPT_SERIALIZATION(LagrangianR);

  /** The Jacobian of the constraints with respect to the generalized coordinates   \f$ q \f$
   *  i.e.  \f$ \nabla^\top_q h(t,q,\dot q,\ldots) \f$
   */
  std::shared_ptr<siconos::algebra::MapType> jacobianhOver_q_view_{nullptr};

  /** internal (optional) storage used for \f$ \nabla^\top_q h(t,q,\dot q,\ldots) \f$ */
  std::unique_ptr<std::vector<double>> jacobianhOver_q_internal_storage_{nullptr};

  /** True if \f$ \nabla^\top_q h(t,q,\dot q,\ldots) \f$ is a constant matrix */
  bool hasConstantJacobianhOver_q_{false};

  /** basic constructor
   *
   *  \param lagType the sub-type of the relation
   */
  LagrangianR(RelationSubType lagType) : Relation(RelationType::Lagrangian, lagType) {}

 public:
  /** destructor
   */
  virtual ~LagrangianR() noexcept = default;

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction using this relation
   */
  inline void initialize(Interaction &inter) override {};

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  virtual void checkSize(Interaction &inter) override {};
  // Does nothing by default. Reimplement if required.

  /** \return a read-only view on \f$ \nabla^\top_q h(q, \ldots) \f$ matrix */
  inline const auto jacobianhOver_q() const {
    return siconos::algebra::ConstMapType(jacobianhOver_q_view_->data(),
                                          jacobianhOver_q_view_->rows(),
                                          jacobianhOver_q_view_->cols());
  }

  /** main relation members display */
  void display() const override;

  // Visitors stuff
  void accept(std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const override;
};

}  // namespace siconos::modeling

#endif  // LAGRANGIANRELATION_H
