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
/*! \file FirstOrderLinearTIR.hpp

 */
#ifndef FirstOrderLinearTIR_H
#define FirstOrderLinearTIR_H

#include "FirstOrderR.hpp"

namespace siconos::modeling {
/**
  First Order Linear Relation with time-invariant coeff.

    \f[
    y &=& Cx(t) + D\lambda + e \\
    R &=& B\lambda
    \f]

    Storage for:
    - B: \f$ \nabla_{\lambda} g(x,t,\lambda) \f$ in FirstOrderR
    - C: \f$ \nabla_{x} h(x,t,\lambda) \f$ in FirstOrderR
    - D: \f$ \nabla_{\lambda} h(x,t,\lambda) \f$ in FirstOrderR
    - e: attribute in this class

*/
class FirstOrderLinearTIR : public FirstOrderR {
 protected:
  ACCEPT_SERIALIZATION(FirstOrderLinearTIR);

  /** e */
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

  /** True if e(t) is taken into account */
  bool haseVector_{false};

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction that owns this relation
   */
  void initialize(Interaction &inter) override;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  void checkSize(Interaction &inter) override;

 public:
  /** minimal constuctor. Use setXXX functions to fix B,C, D ... */
  FirstOrderLinearTIR() : FirstOrderR(RelationSubType::LinearTIR) {};

  /** Build a time-invariant coeff. linear relation
   *
   *  \param C matrix coeff C
   *  \param B matrix coeff B
   */
  FirstOrderLinearTIR(Eigen::Ref<siconos::algebra::SiconosMatrix> C,
                      Eigen::Ref<siconos::algebra::SiconosMatrix> B);

  /** destructor */
  virtual ~FirstOrderLinearTIR() noexcept = default;

  /** \return a read-only view on B(t) matrix */
  inline const auto B() const { return jacobiangOver_lambda(); }

  /** Set a constant B matrix for the system
   *
   *  \param newValue B matrix
   *
   */
  void setConstantB(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** \return a read-only view on C(t) matrix */
  inline const auto C() const { return jacobianhOver_state(); }

  /** Set a constant C matrix for the system
   *
   *  \param newValue C matrix
   *
   */
  void setConstantC(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** \return a read-only view on D(t) matrix */
  inline const auto D() const { return jacobianhOver_lambda(); }

  /** Set a constant D matrix for the system
   *
   *  \param newValue D matrix
   *
   */
  void setConstantD(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** \return  a read-only view on e(t) */
  inline const auto eVector() const {
    return siconos::algebra::ConstMapVectorType(eVector_view_->data(), eVector_view_->size());
  }

  /** set a constant e vector
   *
   *  \param neweVector e vector
   */
  void setConstanteVector(Eigen::Ref<siconos::algebra::SiconosVector> neweVector);

  /** True if e(t) is taken into account */
  bool haseVector() const { return haseVector_; }

  // /** default function to compute h = y = Cx(t) + Fz + Dlambda + e
  //  *
  //  *  \param x
  //  *  \param lambda
  //  *  \param z
  //  *  \param y the resulting vector
  //  */
  // void computeh(const siconos::algebra::BlockVector &x,
  //               const siconos::algebra::SiconosVector &lambda,
  //               siconos::algebra::BlockVector &z, siconos::algebra::SiconosVector &y);

  // /** default function to compute g = Blambda
  //  *
  //  *  \param lambda
  //  *  \param r non-smooth input
  //  */
  // void computeg(const siconos::algebra::SiconosVector &lambda,
  //               siconos::algebra::BlockVector &r);

  /** default function to compute y
   *
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param level
   */
  void computeOutput(double time, Interaction &inter, unsigned int level = 0) override;

  /** default function to compute r
   *
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param level
   */
  void computeInput(double time, Interaction &inter, unsigned int level = 0) override;

  /** print the data to the screen
   */
  void display() const override;

  /** determine if the Relation is linear
   *
   *  \return true if the relation is linear.
   */
  inline bool isLinear() override { return true; }

  // Jacobians: required to fullfill base abstract class API but do nothing.
  // Note FP: final would be better than override but swig cannot handle it.
  void computeJach(double time, Interaction &inter) override {};
};
}  // namespace siconos::modeling

#endif
