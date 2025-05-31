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
/*! \file LagrangianCompliantR.hpp

 */
#ifndef LagrangianCompliantR_H
#define LagrangianCompliantR_H

#include "LagrangianR.hpp"

namespace siconos::modeling {
/**
    Lagrangian Compliant Relation: Scleronomous, Non-Holonomic

    This class provides tools to describe Lagrangian (2nd order) non linear relations, with:

    \f$
    y = h(q,\lambda(t))
    \f$

    \f$
    \dot y = G0(q,\lambda(t))\dot q + G1((q,\lambda(t))\dot\lambda(t)
    \f$

    \f$
    p = G0^t(q,\lambda(t))\lambda(t)
    \f$

    with
    \f$
    G0(q,\lambda(t)) = \nabla_q h(q,\lambda(t))
    \f$

    \f$
    G1(q,\lambda(t)) = \nabla_{\lambda}h(q,\lambda(t))
    \f$

    The following operators can be set by user-defined functions:

    - \f$ h(q, \lambda) \f$
    - \f$ G0(q, \lambda) \f$
    - \f$ G1(q, \lambda) \f$

*/
class LagrangianCompliantR : public LagrangianR {
 protected:
  /** function wrapper used to compute \f$ h(q, \lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionBVV_V computeh_{nullptr};

  /** function wrapper used to compute  \f$ G_0=\nabla^\top_q h(q, \lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionBVV_M computejacobianhOver_q_{nullptr};

  /** jacobian matrix \f$G_1 = \nabla^\top_{\lambda}h(q, \lambda) \f$  */
  std::shared_ptr<siconos::algebra::SiconosMatrix> jacobianhOver_lambda_{nullptr};

  /** function wrapper used to compute  \f$ G_1=\nabla^\top_{\lambda} h(q, \lambda) \f$ */
  siconos::modeling::func_prototypes::FunctionBVV_M computejacobianhOver_lambda_{nullptr};

  ACCEPT_SERIALIZATION(LagrangianCompliantR);

  /** initialize G matrices or components specific to derived classes
   *
   *  \param inter : the Interaction
   */
  void initialize(Interaction &inter) override;

 public:
  /** default and only constructor */
  LagrangianCompliantR() : LagrangianR(RelationSubType::CompliantR) {};

  /** destructor */
  virtual ~LagrangianCompliantR() noexcept = default;

  /** set a user-defined function to compute \f$ h(q, \lambda) \f$  \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehFunction(const siconos::modeling::func_prototypes::FunctionBVV_V &fct);

  /**
    to compute the output y = h(q, \lambda) of the Relation

    \param q coordinates of the dynamical systems involved in the relation
    \param lambda interaction  \f$ \lambda \f$  vector
    \param y the resulting vector
  */
  virtual void computeh(const siconos::algebra::BlockVector &q,
                        const Eigen::Ref<const siconos::algebra::SiconosVector> &lambda,
                        Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** set a user-defined function to compute \f$ \nabla^\top_q h(q, \lambda) \f$ \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianhOver_qFunction(
      const siconos::modeling::func_prototypes::FunctionBVV_M &fct);

  /** Computes \f$ \nabla^\top_q h(q, \lambda) \f$
   *  \param q coordinates of the dynamical systems involved in the relation
   *  \param lambda interaction  \f$ \lambda \f$  vector
   */
  virtual void computeJacobianhOver_q(
      const siconos::algebra::BlockVector &q,
      const Eigen::Ref<const siconos::algebra::SiconosVector> &lambda);

  /** set a user-defined function to compute \f$ \nabla^\top_{\lambda} h(q, \lambda) \f$ \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianhOver_lambdaFunction(
      const siconos::modeling::func_prototypes::FunctionBVV_M &fct);

  /** Computes \f$ \nabla^\top_{\lambda} h(q, \lambda) \f$
   *  \param q coordinates of the dynamical systems involved in the relation
   *  \param lambda interaction  \f$ \lambda \f$  vector
   */
  virtual void computeJacobianhOver_lambda(
      const siconos::algebra::BlockVector &q,
      const Eigen::Ref<const siconos::algebra::SiconosVector> &lambda);

  /*  \return a read-only view on the matrix \f$ \nabla^\top_{\lambda}h(q,\lambda) \f$*/
  inline const siconos::algebra::ConstMapType jacobianhOver_lambda() const override {
    return siconos::algebra::ConstMapType(jacobianhOver_lambda_->data(),
                                          jacobianhOver_lambda_->rows(),
                                          jacobianhOver_lambda_->cols());
  }

  /** to compute output
   *
   *  \param time the current time
   *  \param inter the Interaction owning y
   *  \param derivativeNumber the number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(double time, Interaction &inter,
                     unsigned int derivativeNumber = 0) override;
  /** to compute the input
   *
   *  \param time the current time
   *  \param inter the Interaction owning lambda
   *  \param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction &inter, unsigned int level = 0) override;

  /** compute all the H Jacobian */
  void computeJach(double time, Interaction &inter) override;
};
}  // namespace siconos::modeling

#endif
