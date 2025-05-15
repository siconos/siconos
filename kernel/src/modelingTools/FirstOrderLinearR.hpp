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
/*! \file FirstOrderLinearR.hpp

 */
#ifndef FirstOrderLinearR_H
#define FirstOrderLinearR_H

#include "FirstOrderR.hpp"

namespace siconos::modeling {

/**
   First Order Linear Relation with time-dependant operators:

   \f[

   y &=& C(t)x(t) + D(t)\lambda + e(t) \\
   R &=& B(t) \lambda
   \f]

   The following operators can be can be set by user-defined functions:
   B(t), C(t), D(t), e(t)

   Storage for:
    - B: \f$ \nabla_{\lambda} g(x,t,\lambda) \f$ in FirstOrderR
    - C: \f$ \nabla_{x} h(x,t,\lambda) \f$ in FirstOrderR
    - D: \f$ \nabla_{\lambda} h(x,t,\lambda) \f$ in FirstOrderR
    - e: attribute in this class

*/
class FirstOrderLinearR : public FirstOrderR {
 protected:
  ACCEPT_SERIALIZATION(FirstOrderLinearR);

  /** e(t) operator */
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

  /** internal (optional) storage used for e(t) */
  std::unique_ptr<std::vector<double>> eVector_internal_storage_{nullptr};

  /** function wrapper used to compute e(t) */
  siconos::modeling::func_prototypes::FunctionS_V computeeVector_{nullptr};

  /** True if e(t) is taken into account and constant */
  bool hasConstanteVector_{false};

  /** True if e(t) is taken into account */
  bool haseVector_{false};

  /** function wrapper used to compute  B(t) */
  siconos::modeling::func_prototypes::FunctionS_M computeB_{nullptr};

  /** function wrapper used to compute  C(t) */
  siconos::modeling::func_prototypes::FunctionS_M computeC_{nullptr};

  /** function wrapper used to compute D(t) */
  siconos::modeling::func_prototypes::FunctionS_M computeD_{nullptr};

  /** initialize the relation (check sizes, memory allocation in workV and workM ...)
   *
   *  \param inter Interaction using this Relation
   */
  void initialize(Interaction &inter) override;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  void checkSize(const Interaction &inter) const override;

 public:
  /** default (and only) constructor */
  FirstOrderLinearR() : FirstOrderR(RelationSubType::LinearR) {};

  /** destructor
   */
  virtual ~FirstOrderLinearR() noexcept = default;

  /** \return a read-only view on B(t) matrix */
  inline auto B() const { return jacobiangOver_lambda(); }

  /** Set a constant B matrix for the system
   *
   *  \param newValue B matrix
   *
   */
  void setConstantB(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue);

  /** set a user-defined function to compute B(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeBFunction(const siconos::modeling::func_prototypes::FunctionS_M &fct);

  /** Computes B(t)
   *  \param time current time value
   */
  void computeB(double time);

  /** \return a read-only view on C(t) matrix */
  inline auto C() const { return jacobianhOver_state(); }

  /** Set a constant C matrix for the system
   *
   *  \param newValue C matrix
   *
   */
  void setConstantC(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue);

  /** set a user-defined function to compute C(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeCFunction(const siconos::modeling::func_prototypes::FunctionS_M &fct);

  /** Computes C(t)
   *  \param time current time value
   */
  void computeC(double time);

  /** \return a read-only view on D(t) matrix */
  inline auto D() const { return jacobianhOver_lambda(); }

  /** Set a constant D matrix for the system
   *
   *  \param newValue D matrix
   *
   */
  void setConstantD(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue);

  /** set a user-defined function to compute D(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeDFunction(const siconos::modeling::func_prototypes::FunctionS_M &fct);

  /** Computes D(t)
   *  \param time current time value
   */
  void computeD(double time);

  /** \return  a read-only view on e(t) */
  inline auto eVector() const {
    return siconos::algebra::ConstMapVectorType(eVector_view_->data(), eVector_view_->size());
  }

  /** set a constant e vector
   *
   *  \param neweVector e vector
   */
  void setConstanteVector(Eigen::Ref<siconos::algebra::SiconosVector> neweVector);

  /** True if e(t) is taken into account */
  bool haseVector() const { return haseVector_; }

  /** set a user-defined function to compute external forces
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeeVectorFunction(const siconos::modeling::func_prototypes::FunctionS_V &fct);

  /** Computes e(t)
   *  \param time current time value
   */
  void computee(double time);

  // /**
  //    to compute the output y = h(t,x,...) of the Relation

  //    \param time current time value
  //    \param x coordinates of the dynamical systems involved in the relation
  //    \param lambda interaction \f$\lambda\f$ vector
  //    \param z user defined parameters (optional)
  //    \param y the resulting vector
  // */
  // virtual void computeh(double time, const siconos::algebra::BlockVector &x,
  //                       const siconos::algebra::SiconosVector &lambda,
  //                       siconos::algebra::BlockVector &z, siconos::algebra::SiconosVector
  //                       &y);

  // /**
  //    to compute the nonsmooth input r = g(t,x,...) of the Relation

  //    \param time current time value
  //    \param lambda interaction \f$\lambda\f$ vector
  //    \param z user defined parameters (optional)
  //    \param r the resulting vector
  // */
  // virtual void computeg(double time, const siconos::algebra::SiconosVector &lambda,
  //                       siconos::algebra::BlockVector &z, siconos::algebra::BlockVector
  //                       &r);

  /** default function to compute y
   *
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param level dummy parameter, always=0
   */
  void computeOutput(double time, Interaction &inter, unsigned int level = 0) override;

  /** default function to compute r
   *
   *  \param time current time
   *  \param inter Interaction using this Relation
   *  \param level dummy parameter, always=0
   */
  void computeInput(double time, Interaction &inter, unsigned int level = 0) override;

  /** print the data to the screen
   */
  void display() const override;

  /** determines if the Relation is linear
   *
   *  \return true if the relation is linear.
   */
  bool isLinear() const override { return true; }

  // Jacobians: required to fullfill base abstract class API but do nothing.
  // Note FP: final would be better than override but swig cannot handle it.
  void computeJach(double time, Interaction &inter) override {};
};
}  // namespace siconos::modeling

#endif
