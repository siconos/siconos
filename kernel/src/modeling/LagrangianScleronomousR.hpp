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
/*! \file LagrangianScleronomousR.hpp

 */
#ifndef LagrangianScleronomousR_H
#define LagrangianScleronomousR_H

#include "LagrangianR.hpp"

namespace siconos::modeling {
/**
  Scleronomic Lagrangian (Non Linear) Relations

  This class provides tools to describe Lagrangian (2nd order) non linear relations, with:

    \f[
    y = h(q)
    \f]

    \f[
    \dot y = \nabla^\top_q h(q) \dot q
    \f]

    or more generally

    \f[
     \dot y = H(q) \dot q
    \f]

    and by duality

    \f[
    p = \nabla_q h(q)\lambda
    \f]

    or more generally

    \f[
    p = H^\top(q)\lambda
    \f]

    with

    \f[
    H^\top(q) = \nabla_q h(q)
    \f]

    The following operators can be set by user-defined functions:

    - \f$ h(q) \f$
    - \f$ \nabla_q h(q) \f$

   */
class LagrangianScleronomousR : public LagrangianR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianScleronomousR);

  /** function wrapper used to compute \f$ h(q) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_V computeh_{nullptr};

  /** function wrapper used to compute  \f$ \nabla^\top_q h(q) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_M computejacobianhOver_q_{nullptr};

  /** \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q)) \f$
   *  This value is useful to compute the second-order
   *  derivative of the constraints with respect to time.
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> jacobianhOver_q_dot_{nullptr};

  /** function wrapper used to compute \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q)) \f$
   */
  siconos::modeling::func_prototypes::FunctionBVBV_M computejacobianhOver_q_dot_{nullptr};

  /** True if  \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q)) \f$ is taken into account */
  bool hasJacobianhOver_q_dot_{false};

  /** \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q)).\dot q\f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> jacobianhOver_q_dot_X_qdot_{nullptr};

 public:
  /** basic constructor */
  LagrangianScleronomousR() : LagrangianR(RelationSubType::ScleronomousR) {}

  /** destructor */
  virtual ~LagrangianScleronomousR() noexcept = default;

  void initialize(Interaction &inter) override;

  /** set a user-defined function to compute \f$ h(q) \f$  \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehFunction(const siconos::modeling::func_prototypes::FunctionBV_V &fct);

  /**
    to compute the output y = h(q) of the Relation

    \param q coordinates of the dynamical systems involved in the relation
    \param y the resulting vector
  */
  virtual void computeh(const siconos::algebra::BlockVector &q,
                        Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** Set a constant  \f$ \nabla^\top_q h(q) \f$
   *
   *  \param newValue the constant matrix
   *
   */
  void setConstantJacobianhOver_q(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** set a user-defined function to compute \f$ \nabla^\top_q h(q) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianhOver_qFunction(
      const siconos::modeling::func_prototypes::FunctionBV_M &fct);

  /** Computes \f$ \nabla^\top_q h(q) \f$
   * \param q coordinates of the dynamical systems involved in the relation
   */
  virtual void computeJacobianhOver_q(const siconos::algebra::BlockVector &q);

  /** set a user-defined function to compute
   *  \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q))\f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputejacobianhOver_q_dotFunction(
      const siconos::modeling::func_prototypes::FunctionBVBV_M &fct);

  /** Update \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q))\f$
   *  \param q 'list' of state vectors (for all ds involved in the interaction)
   *  \param qdot 'list' of state vectors (for all ds involved in the interaction)
   */
  void computejacobianhOver_q_dot(const siconos::algebra::BlockVector &q,
                                  const siconos::algebra::BlockVector &qdot);

  /** \return a read-only view on \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q)).\dot q\f$
   * vector */
  inline auto jacobianhOver_q_dot_X_qdot() const {
    return siconos::algebra::ConstMapVectorType(jacobianhOver_q_dot_X_qdot_->data(),
                                                jacobianhOver_q_dot_X_qdot_->size());
  }

  /** to compute \f$ \frac{\partial}{\partial t}(\nabla^T_{q} h(q)).\dot q\f$
   *
   *  \param time double, current time
   *  \param inter interaction
   */
  void computeJacobianhOver_q_dot_X_qdot(double time, Interaction &inter);

  /** compute all the H Jacobian
   *
   *  \param time double, current time
   *  \param inter interaction that owns the relation
   *  \param interProp
   */
  void computeJach(double time, Interaction &inter) override;

  /** to compute output
   *
   *  \param time the current time
   *  \param inter interaction that owns the relation
   *  \param derivativeNumber number of the derivative to compute, optional,
   *  default = 0.
   */
  void computeOutput(double time, Interaction &inter,
                     siconos::algebra::blocks::size_type derivativeNumber = 0) override;

  /** to compute p
   *
   *  \param time the current time
   *  \param inter interaction that owns the relation
   *  \param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction &inter,
                    siconos::algebra::blocks::size_type level = 0) override;

  virtual void accept(relations::Visitor &tourist) const override { tourist.visit(*this); }
};
}  // namespace siconos::modeling

#endif  // LAGRANGIANRELATION_H
