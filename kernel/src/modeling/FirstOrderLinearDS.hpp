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
#include "SiconosMatrix.hpp"

namespace siconos::modeling {

/**
   First Order Linear Systems

   \f$ M(t) \dot x = A(t)x(t)+ b(t) + r, \quad x(t_0)=x_0 \f$

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
  siconos::algebra::DenseVectorStorage bVector_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_bVector(F&& f) const {
    return siconos::algebra::visitStorage(bVector_storage_, std::forward<F>(f),
                                          "bVector_storage_");
  }

  /** function wrapper used to compute b(t) */
  siconos::modeling::func_prototypes::FunctionS_V computebVector_{nullptr};

  /** True if b(t) is taken into account and constant */
  bool hasConstantbVector_{false};

  /** True if b(t) is taken into account */
  bool hasbVector_{false};

  // /** default constructor */
  // FirstOrderLinearDS() = default;

 public:
  using FirstOrderNonLinearDS::FirstOrderNonLinearDS;

  /** Constructor from initial state
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means that for initial and A
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   *  @param[in] x0 initial state
   *  @param[in] A matrix coeff A (in \f$ M\dot x = Ax+b \f$)
   *  @param[in] b vector coeff b (in \f$ M\dot x = Ax+b \f$)
   *  @param tag Pass siconos::algebra::alias_t to select this overload
   * (rather than copy version)
   */
  FirstOrderLinearDS(Eigen::Ref<siconos::algebra::SiconosVector> x0,
                     Eigen::Ref<siconos::algebra::SiconosDenseMatrix> A,
                     Eigen::Ref<siconos::algebra::SiconosVector> b,
                     siconos::algebra::AliasTag);

  /** Constructor from initial state
   *
   *  initial state and A attributes will be initialised (copied)
   *  from the input vectors/matrices
   *
   *  @param[in] x0 initial state
   *  @param[in] A matrix coeff A (in \f$ M\dot x = Ax+b \f$)
   *  @param[in] b vector coeff b (in \f$ M\dot x = Ax+b \f$)
   *  @param tag Pass siconos::algebra::copy_t to select this overload
   * (rather than alias version)
   */
  FirstOrderLinearDS(const siconos::algebra::SiconosVector& x0,
                     const siconos::algebra::SiconosDenseMatrix& A,
                     const siconos::algebra::SiconosVector& b, siconos::algebra::CopyTag);

  /** Copy constructor
   *
   *  \param FOLDS the original FirstOrderLinearDS we want to copy
   */
  FirstOrderLinearDS(const FirstOrderLinearDS& FOLDS);

  /** destructor */
  virtual ~FirstOrderLinearDS() noexcept = default;

  /** \return a read-only view on A(t) matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> A() const {
    return jacobianfOver_x();
  }

  /** \return true if A is taken into account */
  auto hasA() const { return hasJacobianfOver_x_; }

  /** \return true if A is taken into account and is time independant */
  auto hasConstantA() const { return hasConstantJacobianfOver_x_; }

  /** @brief set a constant A matrix for the system
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   * @param newValue matrix to be copied. Its shape must match dimension() x dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  void setConstantA(const siconos::algebra::SiconosDenseMatrix& newValue,
                    siconos::algebra::CopyTag tag);

  /** @brief set a constant A matrix for the system
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue external force vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setConstantA(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
                    siconos::algebra::AliasTag tag);

  /** set a user-defined function to compute A(t)
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeAFunction(const siconos::modeling::func_prototypes::FunctionS_M& fct);

  /** Computes A(t)
   *  \param time current time value
   */
  void computeA(double time);

  /** \return  a read-only view on f(x,t)  */
  Eigen::Ref<const siconos::algebra::SiconosVector> bVector() const {
    return use_bVector(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector>(v); });
  }

  /** @brief set a constant b(t)  vector
   *
   * Warning : deep copy of the provided vector into internal attribute
   *
   * @param newValue vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   * version)
   *
   */
  void setConstantbVector(const siconos::algebra::SiconosVector& newValue,
                          siconos::algebra::CopyTag tag);

  /** @brief set a constant b(t) vector
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setConstantbVector(Eigen::Ref<siconos::algebra::SiconosVector> newValue,
                          siconos::algebra::AliasTag tag);

  /** True if b(t) is taken into account */
  bool hasbVector() const { return hasbVector_; }

  /** set a user-defined function to compute external forces
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputebVectorFunction(const siconos::modeling::func_prototypes::FunctionS_V& fct);

  /** Remove any plugin connection for b(t)
   * Seems to be useful only for MatrixIntegrator and control stuff. Use it with care ...
   */
  void clearbVector();

  /** Computes b(t)
   *  \param time current time value
   */
  void computeb(double time);

  bool isbVectorPlugged() const { return computebVector_ != nullptr; }

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
