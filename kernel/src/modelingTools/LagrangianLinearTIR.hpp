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
/*! \file LagrangianLinearTIR.hpp

 */
#ifndef LAGRANGIANLINEARRELATION_H
#define LAGRANGIANLINEARRELATION_H

#include "LagrangianR.hpp"

namespace siconos::modeling {
/**
    Lagrangian Linear Relation.

    Lagrangian Relation with:

    \f$ y= Cq + e \f$

    \f$ p = C^t \lambda \f$

    C is the only required input to built a LagrangianLinearTIR.
    C storage is \f$ \nabla^\top_q h(q, ...)\f$ from
    LagrangianR

*/
class LagrangianLinearTIR : public LagrangianR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianLinearTIR);

  /** e*/
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

 public:
  // /** Default constructor */
  // LagrangianLinearTIR() : LagrangianR(RelationSubType::LinearTIR) {};

  /** Build a relation with only C matrix defined. Use set.. if e is required.
   *
   *  \param C the matrix C
   */
  LagrangianLinearTIR(Eigen::Ref<siconos::algebra::SiconosMatrix> C);

  /** destructor */
  virtual ~LagrangianLinearTIR() noexcept = default;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  void checkSize(const Interaction &inter) const override;
  ;

  /** \return a read-only view on the C matrix */
  inline const auto CMatrix() const { return jacobianhOver_q(); }

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
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param derivativeNumber the derivative of y we want to compute
   */
  void computeOutput(double time, Interaction &inter,
                     unsigned int derivativeNumber = 0) override;

  /** default function to compute r
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param level the derivative of lambda we want to compute
   */
  void computeInput(double time, Interaction &inter, unsigned int level = 0) override;

  /** compute all the H Jacobian
   *
   *  \param time not used
   *  \param inter the Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJach(double time, Interaction &inter) override {}

  /** print the data to the screen
   */
  void display() const override;

  /** \return true if the relation is linear.
   */

  bool isLinear() override  // final would be better but swig does not like it
  {
    return true;
  }
};

}  // namespace siconos::modeling
#endif
