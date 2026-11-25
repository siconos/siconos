/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/** @file LagrangianSparseRheonomousR.hpp  Rhenonomous (h(q,t)) Lagrangian (2nd order) non
 * linear relations, sparse storage */

#ifndef LagrangianSparseRheonomousR_H
#define LagrangianSparseRheonomousR_H

#include "BlockVector.hpp"
#include "LagrangianSparseR.hpp"
#include "StorageTools.hpp"

namespace siconos::modeling {

/**
 * @brief Lagrangian (2nd order) non linear relations, with y=h(q,t) (sparse storage)
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
 * Sparse-storage counterpart of LagrangianRheonomousR (same notation, see LagrangianR).
 * The following operators can be set by user-defined functions:
 *
 *   - \f$ h(q,t) \f$ (setComputehFunction);
 *   - \f$ H(q,t) \f$ (setConstantContactMatrix, setComputeContactMatrixFunction);
 *   - \f$ \frac{\partial h}{\partial t}(q,t) \f$ (setComputedhdtFunction, or override
 *     computedhdt).
 *
 * Former names (deprecated, kept for one release): *JacobianhOver_q* = *H*, hdot = dhdt,
 * setComputehdotFunction = setComputedhdtFunction, computehdot = computedhdt (a class
 * overriding computehdot is still called).
 */
class LagrangianSparseRheonomousR : public LagrangianSparseR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianSparseRheonomousR);

  /** function wrapper used to compute \f$ h(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_V computeh_{nullptr};

  /** function wrapper used to compute  \f$ H(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_Ms computeContactMatrix_{nullptr};
  siconos::modeling::func_prototypes::RFunctionBVS_Ms computeContactMatrix_python_{nullptr};

  /** \f$ \frac{\partial }{\partial t}h(q,t) \f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> dhdt_{nullptr};

  /** function wrapper used to compute \f$ \frac{\partial }{\partial t}h(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBVS_V computedhdt_{nullptr};

  /** True if \f$ \frac{\partial }{\partial t}h(q,t) \f$ is taken into account */
  bool hasdhdt_{false};

 public:
  /** default and only constructor */
  LagrangianSparseRheonomousR() : LagrangianSparseR(RelationSubType::RheonomousR) {}

  /** destructor */
  virtual ~LagrangianSparseRheonomousR() noexcept = default;

  /** initialize G matrices or components specific to derived classes.
   *
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

  /** @brief Set a constant \f$ H \f$
   *  The input matrix is copied.
   *
   *  @param newValue H matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setConstantContactMatrix(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    // Warning: we can't check dimensions --> they depend on the interaction
    contactMatrix_storage_ =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
    hasContactMatrix_ = true;
    hasConstantContactMatrix_ = true;
    computeContactMatrix_ = nullptr;
    computeContactMatrix_python_ = nullptr;
  }

  /** @brief set a constant \f$ H \f$
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   * @param jac new stiffness matrix
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   */
  void setConstantContactMatrix(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac,
                                siconos::algebra::AliasTag);

  /** set a user-defined function to compute \f$ H \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  template <typename Func>
  void setComputeContactMatrixFunction(Func&& f) {
    using Ret = std::invoke_result_t<Func, const siconos::algebra::BlockVector&, double,
                                     siconos::algebra::SiconosSparseMatrix&>;
    // if (!std::holds_alternative<siconos::algebra::OwnedSparse>(contactMatrix_storage_)) {
    ////// contactMatrix_storage_ =
    ////////  std::make_unique<siconos::algebra::SiconosSparseMatrix>(ndof_, ndof_);
    // We can't check size of the Jacobian, it depends on the relation
    //}

    hasContactMatrix_ = true;
    hasConstantContactMatrix_ = false;

    if constexpr (std::is_same_v<Ret, void>) {
      computeContactMatrix_ = std::forward<Func>(f);
      computeContactMatrix_python_ = nullptr;
    } else {
      computeContactMatrix_python_ = std::forward<Func>(f);
      computeContactMatrix_ = nullptr;
    }
  }

  /** to compute  \f$ H(q, t)  \f$
   *  @param positions q vectors
   *  @param time the current time
   */
  void computeContactMatrix(const siconos::algebra::BlockVector& positions, double time);

  /** @return  a read-only view on \f$ \frac{\partial }{\partial t}h(q,t) \f$  */
  inline auto dhdt() const {
    return siconos::algebra::ConstMapVectorType(dhdt_->data(), dhdt_->size());
  }

  /** set a user-defined function to compute \f$ \frac{\partial }{\partial t}h(q,t) \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  void setComputedhdtFunction(const siconos::modeling::func_prototypes::FunctionBVS_V& fct);

  /** Update \f$ \frac{\partial }{\partial t}h(q,t) \f$
   *  @param position 'list' of state vectors (for all ds involved in the interaction)
   *  @param time the current time
   */
  virtual void computedhdt(const siconos::algebra::BlockVector& position, double time);

 private:
  /** flag of the dispatch to the deprecated virtual name (see RelationCompatibility.hpp) */
  bool computedhdtLegacy_{false};

 public:
  /** compute all the jacobians of h */
  void computeJach(double time, Interaction& inter) override;

  /** to compute output
   *
   *  @param time current time
   *  @param inter the Interaction
   *  @param derivativeNumber number of the derivative to compute, optional,
   *  default = 0.
   */
  void computeOutput(double time, Interaction& inter,
                     siconos::algebra::blocks::size_type derivativeNumber = 0) override;

  /** to compute p
   *
   *  @param time current time
   *  @param inter the Interaction
   *  @param level "derivative" order of lambda used to compute input
   */
  void computeInput(double time, Interaction& inter,
                    siconos::algebra::blocks::size_type level = 0) override;

  /** \name Deprecated names, kept for one release */
  ///@{
  template <typename T>
  [[deprecated("use setConstantContactMatrix")]] void setConstantJacobianhOver_q(
      T&& newValue, siconos::algebra::CopyTag tag) {
    setConstantContactMatrix(std::forward<T>(newValue), tag);
  }
  [[deprecated("use setConstantContactMatrix")]] void setConstantJacobianhOver_q(
      Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac, siconos::algebra::AliasTag tag) {
    setConstantContactMatrix(jac, tag);
  }
  template <typename Func>
  [[deprecated("use setComputeContactMatrixFunction")]] void setComputeJacobianhOver_qFunction(
      Func&& f) {
    setComputeContactMatrixFunction(std::forward<Func>(f));
  }
  [[deprecated("use computeContactMatrix")]] void computeJacobianhOver_q(
      const siconos::algebra::BlockVector& positions, double time) {
    computeContactMatrix(positions, time);
  }
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
