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

/*! \file FirstOrderLinearDS.hpp
 */
#ifndef FirstOrderLinearDS_H
#define FirstOrderLinearDS_H

#include "FirstOrderNonLinearDS.hpp"

namespace siconos::modeling {

/**
   First Order Linear Systems - \f$ M(t) \dot x = A(t)x(t)+ b(t) + r, \quad
   x(t_0)=x_0 \f$.

   This class represents first order linear systems of the form:

   \f[
   M(t) \dot x = A(t)x(t)+ b(t) + r,
   x(t_0)=x_0
   \f]
   where

   - \f$ x \in R^{n}  \f$ is the state,
   - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
   - \f$ M \in R^{n\times n} \f$ is an invertible matrix
   - \f$ A \in R^{n\times n} \f$
   - \f$ b \in R^{n} \f$


   The following operators can be plugged, in the usual way (see User Guide)

   - \f$ A(t) \f$
   - \f$ b(t) \f$
   - \f$ M(t) \f$

   Storage:

    - \f$ A(t) \f$ in jacobianRhsOver_x (\f$ \nabla_x f(x,t) \f$) of FirstOrderNonLinearDS
    - \f$ b(t) \f$ local
    - \f$ M(t) \f$ in MMatrix FirstOrderNonLinearDS


*/
class FirstOrderLinearDS : public FirstOrderNonLinearDS {
 protected:
  ACCEPT_SERIALIZATION(FirstOrderLinearDS);

  /** function wrapper used to compute  A(t) */  // to replace computeJacobianfOver_x_
  siconos::modeling::func_prototypes::FunctionS_M computeA_{nullptr};

  /** b(t) operator */
  std::shared_ptr<siconos::algebra::MapVectorType> bVector_view_{nullptr};

  /** internal (optional) storage used for b(t) */
  std::unique_ptr<std::vector<double>> bVector_internal_storage_{nullptr};

  /** function wrapper used to compute b(t) */
  siconos::modeling::func_prototypes::FunctionS_V computebVector_{nullptr};

  /** True if b(t) is taken into account and constant */
  bool hasConstantbVector_{false};

  /** True if b(t) is taken into account */
  bool hasbVector_{false};

  // /** default constructor */
  // FirstOrderLinearDS() = default;

 public:
  /** constructor from initial state
   *
   *  \param newX0 the initial state of this DynamicalSystem
   */
  FirstOrderLinearDS(Eigen::Ref<siconos::algebra::SiconosVector> x0)
      : FirstOrderNonLinearDS(x0) {};

  /** Build a time-invariant coeff. linear DS
   *
   *  \param newX0 the initial state of this DynamicalSystem
   *  \param A matrix coeff A (in \f$ M\dot x = Ax+b \f$)
   *  \param b vector coeff b
   */
  FirstOrderLinearDS(Eigen::Ref<siconos::algebra::SiconosVector> x0,
                     Eigen::Ref<siconos::algebra::SiconosMatrix> A,
                     Eigen::Ref<siconos::algebra::SiconosVector> b);

  /** Copy constructor
   *
   *  \param FOLDS the original FirstOrderLinearDS we want to copy
   */
  FirstOrderLinearDS(const FirstOrderLinearDS &FOLDS);

  /** destructor */
  virtual ~FirstOrderLinearDS() noexcept = default;

  /** \return a read-only view on A(t) matrix */
  inline  auto A() const { return jacobianfOver_x(); }

  /** \return true if A is taken into account */
  auto hasA() { return hasJacobianfOver_x_; }

  /** Set a constant A matrix for the system
   *
   *  \param newValue A matrix
   *
   */
  void setConstantA(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** set a user-defined function to compute A(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeAFunction(const siconos::modeling::func_prototypes::FunctionS_M &fct);

  /** Computes A(t)
   *  \param time current time value
   */
  void computeA(double time);

  /** \return  a read-only view on b(t) */
  inline  auto bVector() const {
    return siconos::algebra::ConstMapVectorType(bVector_view_->data(), bVector_view_->size());
  }

  /** set a constant e vector
   *
   *  \param newbVector e vector
   */
  void setConstantbVector(Eigen::Ref<siconos::algebra::SiconosVector> newbVector);

  /** True if b(t) is taken into account */
  bool hasbVector() const { return hasbVector_; }

  /** set a user-defined function to compute external forces
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputebVectorFunction(const siconos::modeling::func_prototypes::FunctionS_V &fct);

  /** Computes b(t)
   *  \param time current time value
   */
  void computeb(double time);

  /** update right-hand side for the current state
   *
   *  \param time of interest
   */
  void computeRhs(double time) override;

  /** update \f$ \nabla_x rhs \f$ for the current state
   *
   *  \param time of interest
   */
  void computeJacobianRhsOver_x(double time) override;

  /**
      Call all user-defined operators to update the DS components

      \param time value
  */
  void updatePlugins(double time) override;

  /** data display on screen  */
  void display(bool brief = true) const override;

  /** \return always true */
  bool isLinear() const override { return true; }
};
}  // namespace siconos::modeling
#endif  // FirstOrderLinearDS_H
