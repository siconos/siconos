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
/*! \file LagrangianCompliantLinearTIR.hpp
 */
#ifndef LAGRANGIANCOMPLIANTLINEARRELATION_H
#define LAGRANGIANCOMPLIANTLINEARRELATION_H

#include "LagrangianR.hpp"

namespace siconos::modeling {
/**  Lagrangian Compliant Linear Relation.

This class provides tools to describe Lagrangian (2nd order) linear relations, with:

\f$
y= Cq + e + D\lambda
\f$

\f$
p = C^t \lambda
\f$

Minimal setup requires C and D matrices.

Note: considering the notations of LagrangianR or LagrangianCompliantR, we have

    \f$
    C = \nabla_q h(q,\lambda(t))
    \f$
and
    \f$
    D = \nabla_{\lambda}h(q,\lambda(t))
    \f$

*/
class LagrangianCompliantLinearTIR : public LagrangianR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianCompliantLinearTIR);

  /** D matrix */
  std::shared_ptr<siconos::algebra::MapType> DMatrix_view_{nullptr};

  /** e*/
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

 public:
  /** Default constructor
   */
  LagrangianCompliantLinearTIR() : LagrangianR(RelationSubType::CompliantLinearTIR) {};

  /** create the Relation from a set of data
   *  \param C the matrix C
   *  \param D the matrix D
   */
  LagrangianCompliantLinearTIR(Eigen::Ref<siconos::algebra::SiconosMatrix> C,
                               Eigen::Ref<siconos::algebra::SiconosMatrix> D);

  /** destructor */
  virtual ~LagrangianCompliantLinearTIR() noexcept = default;

  /** initialize LagrangianCompliantLinearTIR specific operators.
   * \param inter an Interaction using this relation
   */
  void initialize(Interaction &inter) override;

  /** check sizes LagrangianCompliantLinearTIR specific operators.
   * \param inter an Interaction using this relation
   */
  void checkSize(Interaction &inter) override;

  /** \return a read-only view on the C matrix */
  inline const auto CMatrix() const { return jacobianhOver_q(); }

  /** \return a read-only view on the D matrix */
  inline const auto DMatrix() const {
    return siconos::algebra::ConstMapType(DMatrix_view_->data(), DMatrix_view_->rows(),
                                          DMatrix_view_->cols());
  }

  /** \return  a read-only view on e vector */
  inline const auto eVector() const {
    return siconos::algebra::ConstMapVectorType(eVector_view_->data(), eVector_view_->size());
  }

  /** True if e is defined */
  bool haseVector() const { return eVector_view_ != nullptr; }

  /** set a constant e vector. Warning: shared memory with input
   *
   *  \param newe e value
   */
  void seteVector(Eigen::Ref<siconos::algebra::SiconosVector> newe);

  /** default function to compute y
   *  \param time dummy parameter for this kind of relation
   *  \param inter the Interaction we want to update
   *  \param derivativeNumber the derivative of y we want to compute
   */
  void computeOutput(double time, Interaction &inter,
                     unsigned int derivativeNumber = 0) override;

  /** default function to compute r
   *  \param time dummy parameter for this kind of relation
   *  \param inter the Interaction we want to update
   *  \param level the derivative of lambda we want to compute
   */
  void computeInput(double time, Interaction &inter, unsigned int level = 0) override;

  /* compute all the H Jacobian
   *  \param time dummy parameter for this kind of relation
   *  \param inter the Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJach(double time, Interaction &inter) override {}

  /** print the data to the screen
   */
  void display() const override;

  /**
   * \return true if the relation is linear.
   */

  bool isLinear() override { return true; }
};
}  // namespace siconos::modeling
#endif
