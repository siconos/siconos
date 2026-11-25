// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/** @file LagrangianRheonomousR.hpp  Rhenonomous (h(q,t)) Lagrangian (2nd order) non linear
 * relations */

#ifndef LagrangianRheonomousR_H
#define LagrangianRheonomousR_H

#include "LagrangianR.hpp"

namespace siconos::modeling {

/**
 * @brief Lagrangian (2nd order) non linear relations, with y=h(q,t) (dense storage)
 *
 * \f[
 * \begin{aligned}
 *   y &= h(q,t), \\
 *   \dot y &= H(q,t)\, \dot q + \frac{\partial h}{\partial t}(q,t), \qquad
 *     H(q,t) = \nabla_q^\top h(q,t)
 *     = \left[ \frac{\partial h_i}{\partial q_j} \right]_{1\le i\le m,\;1\le j\le n}, \\
 *   p &= H(q,t)^\top \lambda.
 * \end{aligned}
 * \f]
 *
 * \f$ H \f$ is the contact configuration matrix (see LagrangianR). The following operators
 * can be set by user-defined functions, or computed by overriding the virtual functions:
 *
 *   - \f$ h(q,t) \f$: setComputehFunction, or override computeh;
 *   - \f$ H(q,t) \f$: setComputeContactMatrixFunction, or override computeContactMatrix;
 *   - \f$ \frac{\partial h}{\partial t}(q,t) \f$: setComputedhdtFunction, or override
 *     computedhdt.
 *
 * <b>Former names</b> (deprecated, kept for one release): setComputeJacobianhOver_qFunction,
 * computeJacobianhOver_q, hdot, setComputehdotFunction, computehdot (a class overriding
 * computeJacobianhOver_q or computehdot is still called).
 */
class LagrangianRheonomousR : public LagrangianR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianRheonomousR);

  /** function wrapper used to compute \f$ h(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_V computeh_{nullptr};

  /** function wrapper used to compute \f$ H(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_M computeContactMatrix_{nullptr};

  /** \f$ \frac{\partial h}{\partial t}(q,t) \f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> dhdt_{nullptr};

  /** function wrapper used to compute \f$ \frac{\partial h}{\partial t}(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_V computedhdt_{nullptr};

  /** True if \f$ \frac{\partial h}{\partial t} \f$ is taken into account */
  bool hasdhdt_{false};

 private:
  /** flags of the dispatch to the deprecated virtual names (see RelationCompatibility.hpp) */
  bool computeContactMatrixLegacy_{false};
  bool computedhdtLegacy_{false};

 public:
  /** default and only constructor */
  LagrangianRheonomousR() : LagrangianR(RelationSubType::RheonomousR) {}

  /** destructor */
  virtual ~LagrangianRheonomousR() noexcept = default;

  /** initialize matrices or components
   *  @param inter the Interaction
   */
  void initialize(Interaction& inter) override;

  /** set a user-defined function to compute \f$ h(q,t) \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehFunction(const siconos::modeling::func_prototypes::FunctionBVS_V& fct);

  /**
    to compute the output y = h(q, t) of the Relation

    @param q coordinates of the dynamical systems involved in the relation
    @param time current time value
    @param y the resulting vector
  */
  virtual void computeh(const siconos::algebra::BlockVector& q, double time,
                        Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** set a user-defined function to compute \f$ H(q, t) \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...), called as
   *  fct(q, t, H)
   */
  void setComputeContactMatrixFunction(
      const siconos::modeling::func_prototypes::FunctionBVS_M& fct);

  /** Computes \f$ H(q, t) \f$ (contactMatrix_view_). The default calls the user-defined
   * function, if any. Override this function to compute \f$ H \f$ in a derived class.
   *
   *  @param q coordinates of the dynamical systems involved in the relation
   *  @param time current time value
   */
  virtual void computeContactMatrix(const siconos::algebra::BlockVector& q, double time);

  /** @return a read-only view on \f$ \frac{\partial h}{\partial t}(q,t) \f$ */
  inline auto dhdt() const {
    return siconos::algebra::ConstMapVectorType(dhdt_->data(), dhdt_->size());
  }

  /** set a user-defined function to compute \f$ \frac{\partial h}{\partial t}(q,t) \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  void setComputedhdtFunction(const siconos::modeling::func_prototypes::FunctionBVS_V& fct);

  /** Update \f$ \frac{\partial h}{\partial t}(q,t) \f$. The default calls the user-defined
   *  function, if any. Override this function to compute it in a derived class.
   *
   *  @param position 'list' of state vectors (for all ds involved in the interaction)
   *  @param time the current time
   */
  virtual void computedhdt(const siconos::algebra::BlockVector& position, double time);

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
  [[deprecated("use setComputeContactMatrixFunction")]] void setComputeJacobianhOver_qFunction(
      const siconos::modeling::func_prototypes::FunctionBVS_M& fct) {
    setComputeContactMatrixFunction(fct);
  }
  /** \deprecated former name of computeContactMatrix. Still called for classes that override
   * it. */
  [[deprecated("override/call computeContactMatrix")]] virtual void computeJacobianhOver_q(
      const siconos::algebra::BlockVector& q, double time);
  [[deprecated("use dhdt()")]] auto hdot() const { return dhdt(); }
  [[deprecated("use setComputedhdtFunction")]] void setComputehdotFunction(
      const siconos::modeling::func_prototypes::FunctionBVS_V& fct) {
    setComputedhdtFunction(fct);
  }
  /** \deprecated former name of computedhdt. Still called for classes that override it. */
  [[deprecated("override/call computedhdt")]] virtual void computehdot(
      const siconos::algebra::BlockVector& position, double time);
  ///@}
};
}  // namespace siconos::modeling
#endif
