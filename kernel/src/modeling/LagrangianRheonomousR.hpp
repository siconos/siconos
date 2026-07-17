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

/** @file LagrangianRheonomousR.hpp  Rhenonomous (h(q,t)) Lagrangian (2nd order) non linear
 * relations */

#ifndef LagrangianRheonomousR_H
#define LagrangianRheonomousR_H

#include "LagrangianR.hpp"

namespace siconos::modeling {

/**
 * @brief describe Lagrangian (2nd order) non linear relations, with h(q,t)
 *
 *
 * \f[
 *
 * \dot y =  H(q,t)\dot q + \frac{\partial }{\partial t}h(q,t)
 *  \f]
 *
 *   and by duality
 *
 *   \f[
 *
 *        p = H^\top(q,t)\lambda
 *   \f]
 *
 *  where
 *
 * \f[
 * H(q,t)
 * =
 * \frac{\partial h}{\partial q}(q,t)
 * =
 * \left[
 * \frac{\partial h_i}{\partial q_j}
 * \right]_{1\le i\le m,\;1\le j\le n}.
 * \f]
 *
 * is the Jacobian matrix of the relation.
 *
 *   The following operators can be set by user-defined functions:
 *
 *   - \f$ h(q,t) \f$
 *   - \f$ H(q,t) \f$
 *   - \f$ \dot h(q,t) \f$
 *
 */
class LagrangianRheonomousR : public LagrangianR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianRheonomousR);

  /** function wrapper used to compute \f$ h(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_V computeh_{nullptr};

  /** function wrapper used to compute  \f$ H(q, t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_M computejacobianhOver_q_{nullptr};

  /** \f$ \frac{\partial }{\partial t}h(q,t) \f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> hdot_{nullptr};

  /** function wrapper used to compute \f$ \frac{\partial }{\partial t}h(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_V computehdot_{nullptr};

  /** True if \f$ \frac{\partial }{\partial t}h(q,t) \f$ is taken into account */
  bool hashdot_{false};

 public:
  /** default and only constructor */
  LagrangianRheonomousR() : LagrangianR(RelationSubType::RheonomousR) {}

  /** destructor
   */
  virtual ~LagrangianRheonomousR() noexcept = default;

  /** initialize G matrices or components specific to derived classes.
   *
   *  \param inter the Interaction
   */
  void initialize(Interaction& inter) override;

  /** \return  a read-only view on \f$ \frac{\partial }{\partial t}h(q,t) \f$  */
  inline auto hdot() const {
    return siconos::algebra::ConstMapVectorType(hdot_->data(), hdot_->size());
  }

  /** set a user-defined function to compute \f$ \frac{\partial }{\partial t}h(q,t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehdotFunction(const siconos::modeling::func_prototypes::FunctionBVS_V& fct);

  /** Update \f$ \frac{\partial }{\partial t}h(q,t) \f$
   *  \param position 'list' of state vectors (for all ds involved in the interaction)
   *  \param time the current time
   */
  virtual void computehdot(const siconos::algebra::BlockVector& position, double time);

  /** set a user-defined function to compute \f$ h(q,t) \f$  \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehFunction(const siconos::modeling::func_prototypes::FunctionBVS_V& fct);

  /**
    to compute the output y = h(q, t) of the Relation

    \param q coordinates of the dynamical systems involved in the relation
    \param time current time value
    \param y the resulting vector
  */
  virtual void computeh(const siconos::algebra::BlockVector& q, double time,
                        Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** set a user-defined function to compute \f$ \nabla^\top_q h(q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianhOver_qFunction(
      const siconos::modeling::func_prototypes::FunctionBVS_M& fct);

  /** Computes \f$ \nabla^\top_q h(q, t) \f$
   *  \param q coordinates of the dynamical systems involved in the relation
   *  \param time current time value
   */
  virtual void computeJacobianhOver_q(const siconos::algebra::BlockVector& q, double time);

  /** compute all the H Jacobian */
  void computeJach(double time, Interaction& inter) override;

  /** to compute output
   *
   *  \param time current time
   *  \param inter the Interaction
   *  \param derivativeNumber number of the derivative to compute, optional,
   *  default = 0.
   */
  void computeOutput(double time, Interaction& inter,
                     siconos::algebra::blocks::size_type derivativeNumber = 0) override;

  /** to compute p
   *
   *  \param time current time
   *  \param inter the Interaction
   *  \param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction& inter,
                    siconos::algebra::blocks::size_type level = 0) override;
};
}  // namespace siconos::modeling
#endif
