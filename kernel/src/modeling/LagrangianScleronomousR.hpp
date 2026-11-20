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

/*! @file LagrangianScleronomousR.hpp Lagrangian scleronomous ( h(q) ) relations (dense
 * storage)*/
#ifndef LagrangianScleronomousR_H
#define LagrangianScleronomousR_H

#include "LagrangianR.hpp"

namespace siconos::modeling {

/**
  @brief Lagrangian (2nd order) non linear relations, with y=h(q) (dense storage)

  \f[
  \begin{aligned}
    y &= h(q), \\
    \dot y &= H(q)\, \dot q, \qquad
      H(q) = \nabla_q^\top h(q)
      = \left[ \frac{\partial h_i}{\partial q_j} \right]_{1\le i\le m,\;1\le j\le n}, \\
    \ddot y &= H(q)\, \ddot q + \dot H(q, \dot q)\, \dot q, \\
    p &= H(q)^\top \lambda.
  \end{aligned}
  \f]

  \f$ H \f$ is the contact configuration matrix (see LagrangianR). The following operators
  can be set by user-defined functions, or computed by overriding the virtual functions:

  - \f$ h(q) \f$: setComputehFunction, or override computeh;
  - \f$ H(q) \f$: setConstantContactMatrix, setComputeContactMatrixFunction, or override
  computeContactMatrix;
  - \f$ \dot H(q, \dot q) \f$: setComputeContactMatrixDotFunction.

  <b>Former names</b> (deprecated, kept for one release): setConstantJacobianhOver_q,
  setComputeJacobianhOver_qFunction, computeJacobianhOver_q (a class overriding it is still
  called), setComputejacobianhOver_q_dotFunction, computeJacobianhOver_q_dot,
  jacobianhOver_q_dot_X_qdot, computeJacobianhOver_q_dot_X_qdot.
*/
class LagrangianScleronomousR : public LagrangianR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianScleronomousR);

  /** function wrapper used to compute \f$ h(q) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_V computeh_{nullptr};

  /** function wrapper used to compute \f$ H(q) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_M computeContactMatrix_{nullptr};

  /** \f$ \dot H = \frac{d}{dt} H(q) \f$, used to compute the second-order time derivative
   *  of \f$ y \f$. */
  std::shared_ptr<siconos::algebra::SiconosMatrix> contactMatrixDot_{nullptr};

  /** function wrapper used to compute \f$ \dot H(q, \dot q) \f$ */
  siconos::modeling::func_prototypes::FunctionBVBV_M computeContactMatrixDot_{nullptr};

  /** True if \f$ \dot H \f$ is taken into account */
  bool hasContactMatrixDot_{false};

  /** \f$ \dot H \dot q \f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> contactMatrixDotV_{nullptr};

  /** Deprecated names of the members above, kept for one release. */
  [[deprecated("use computeContactMatrix_")]] siconos::modeling::func_prototypes::FunctionBV_M&
      computejacobianhOver_q_{computeContactMatrix_};
  [[deprecated("use contactMatrixDot_")]] std::shared_ptr<siconos::algebra::SiconosMatrix>&
      jacobianhOver_q_dot_{contactMatrixDot_};
  [[deprecated("use hasContactMatrixDot_")]] bool& hasJacobianhOver_q_dot_{
      hasContactMatrixDot_};
  [[deprecated("use contactMatrixDotV_")]] std::shared_ptr<siconos::algebra::SiconosVector>&
      jacobianhOver_q_dot_X_qdot_{contactMatrixDotV_};

 private:
  /** flag of the dispatch to the deprecated virtual name (see RelationCompatibility.hpp) */
  bool computeContactMatrixLegacy_{false};

 public:
  /** basic constructor */
  LagrangianScleronomousR() : LagrangianR(RelationSubType::ScleronomousR) {}

  /** destructor */
  virtual ~LagrangianScleronomousR() noexcept = default;

  /** initialize matrices or components
   *  @param inter the Interaction
   */
  void initialize(Interaction& inter) override;

  /** set a user-defined function to compute \f$ h(q) \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehFunction(const siconos::modeling::func_prototypes::FunctionBV_V& fct);

  /**
    to compute the output y = h(q) of the Relation

    @param q coordinates of the dynamical systems involved in the relation
    @param y the resulting vector
  */
  virtual void computeh(const siconos::algebra::BlockVector& q,
                        Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** Set a constant \f$ H \f$ (warning: shared memory, newValue must outlive the relation)
   *
   *  @param newValue the constant matrix
   */
  void setConstantContactMatrix(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** set a user-defined function to compute \f$ H(q) \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...), called as fct(q, H)
   */
  void setComputeContactMatrixFunction(
      const siconos::modeling::func_prototypes::FunctionBV_M& fct);

  /** Computes \f$ H(q) \f$ (contactMatrix_view_). The default calls the user-defined function,
   * if any. Override this function to compute \f$ H \f$ in a derived class.
   *
   *  @param q coordinates of the dynamical systems involved in the relation
   */
  virtual void computeContactMatrix(const siconos::algebra::BlockVector& q);

  /** set a user-defined function to compute \f$ \dot H(q, \dot q) \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...), called as
   *  fct(q, qdot, contactMatrixDot)
   */
  void setComputeContactMatrixDotFunction(
      const siconos::modeling::func_prototypes::FunctionBVBV_M& fct);

  /** Update \f$ \dot H(q, \dot q) \f$
   *  @param q 'list' of state vectors (for all ds involved in the interaction)
   *  @param qdot 'list' of state vectors (for all ds involved in the interaction)
   */
  void computeContactMatrixDot(const siconos::algebra::BlockVector& q,
                               const siconos::algebra::BlockVector& qdot);

  /** @return a read-only view on \f$ \dot H \dot q \f$ */
  inline auto contactMatrixDotV() const {
    return siconos::algebra::ConstMapVectorType(contactMatrixDotV_->data(),
                                                contactMatrixDotV_->size());
  }

  /** compute \f$ \dot H \dot q \f$, the terms of \f$ \ddot y \f$ that do not depend on
   *  \f$ \ddot q \f$
   *
   *  @param time double, current time
   *  @param inter interaction
   */
  void computeContactMatrixDotV(double time, Interaction& inter);

  /** compute all the jacobians of h
   *
   *  @param time double, current time
   *  @param inter interaction that owns the relation
   */
  void computeJach(double time, Interaction& inter) override;

  /** to compute output
   *
   *  @param time current time
   *  @param inter interaction that owns the relation
   *  @param derivativeNumber number of the derivative to compute, optional,
   *  default = 0.
   */
  void computeOutput(double time, Interaction& inter,
                     siconos::algebra::blocks::size_type derivativeNumber = 0) override;

  /** to compute p
   *
   *  @param time current time
   *  @param inter interaction that owns the relation
   *  @param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction& inter,
                    siconos::algebra::blocks::size_type level = 0) override;

  virtual void accept(relations::Visitor& tourist) const override { tourist.visit(*this); }

  /** \name Deprecated names, kept for one release */
  ///@{
  [[deprecated("use setConstantContactMatrix")]] void setConstantJacobianhOver_q(
      Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
    setConstantContactMatrix(newValue);
  }
  [[deprecated("use setComputeContactMatrixFunction")]] void setComputeJacobianhOver_qFunction(
      const siconos::modeling::func_prototypes::FunctionBV_M& fct) {
    setComputeContactMatrixFunction(fct);
  }
  /** \deprecated former name of computeContactMatrix. Still called for classes that override
   * it. */
  [[deprecated("override/call computeContactMatrix")]] virtual void computeJacobianhOver_q(
      const siconos::algebra::BlockVector& q);
  [[deprecated("use setComputeContactMatrixDotFunction")]] void
  setComputejacobianhOver_q_dotFunction(
      const siconos::modeling::func_prototypes::FunctionBVBV_M& fct) {
    setComputeContactMatrixDotFunction(fct);
  }
  [[deprecated("use computeContactMatrixDot")]] void computeJacobianhOver_q_dot(
      const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot) {
    computeContactMatrixDot(q, qdot);
  }
  [[deprecated("use contactMatrixDotV()")]] auto jacobianhOver_q_dot_X_qdot() const {
    return contactMatrixDotV();
  }
  [[deprecated("use computeContactMatrixDotV")]] void computeJacobianhOver_q_dot_X_qdot(
      double time, Interaction& inter) {
    computeContactMatrixDotV(time, inter);
  }
  ///@}
};
}  // namespace siconos::modeling
#endif
