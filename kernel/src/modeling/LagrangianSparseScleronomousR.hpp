/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

/** @file LagrangianSparseScleronomousR.hpp Lagrangian scleronomous ( h(q) ) relations, sparse
 * storage */

#ifndef LagrangianSparseScleronomousR_H
#define LagrangianSparseScleronomousR_H

#include "BlockVector.hpp"
#include "LagrangianSparseR.hpp"
#include "StorageTools.hpp"

namespace siconos::modeling {

/**
 * @brief Lagrangian (2nd order) non linear relations, with y=h(q) (sparse storage)
 *
 *
 * \f[
 *
 * y =  h(q) \\
 *
 * \dot y =  J^h_q(q)\dot q
 *  \f]
 *
 *   and by duality
 *
 *   \f[
 *
 *        p = {J^h_q}^\top(q)\lambda
 *   \f]
 *
 *  where
 *
 * \f[
 * J^h_q(q)
 * =
 * \frac{\partial h}{\partial q}(q)
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
 *   - \f$ h(q) \f$
 *   - \f$ J^h_q(q) \f$
 *   - \f$  \frac{\partial}{\partial t}(J^h_q) \f$
 *
 */
class LagrangianSparseScleronomousR : public LagrangianSparseR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianSparseScleronomousR);

  /** function wrapper used to compute \f$ h(q) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_V computeh_{nullptr};

  /** function wrapper used to compute  \f$ J^h_q(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_Ms computejacobianhOver_q_{nullptr};
  siconos::modeling::func_prototypes::RFunctionBV_Ms computejacobianhOver_q_python_{nullptr};

  /** \f$ \frac{\partial}{\partial t}(J^h_q) \f$
   *  This value is useful to compute the second-order
   *  derivative of the constraints with respect to time.
   */
  siconos::algebra::SparseStorage jacobianhOver_q_dot_storage_{std::monostate{}};

  /** true if \f$ \frac{\partial}{\partial t}(J^h_q) \f$ is not set or set and constant */
  bool hasConstantJacobianhOver_q_dot_{true};

  /** True if \f$ \frac{\partial}{\partial t}(J^h_q) \f$ is required */
  bool hasJacobianhOver_q_dot_{false};

  /**
   * @brief Utility function providing uniform access to the \f$ \frac{\partial}{\partial
   * t}(J^h_q) \f$ matrix.
   *
   * generic and efficient access to the underlying `SiconosSparseMatrix` object, avoiding any
   * copy.
   *
   * @tparam F a callable type (e.g., lambda) taking a reference to
   *           `SiconosSparseMatrix` or `const SiconosSparseMatrix&`.
   *
   * @param f the functor to apply to the matrix.
   * @return whatever is returned by @p f.
   *
   * @throws std::runtime_error if no matrix (owned or external) is set.
   *
   * @note This function does not copy the matrix; it forwards a reference
   *       to the actual internal or external object.
   */
  template <typename F>
  decltype(auto) useJacobianhOver_q_dot(F&& f) const {
    return siconos::algebra::visitStorage(jacobianhOver_q_dot_storage_, std::forward<F>(f),
                                          "jacobianhOver_q_dot_storage_");
  }

  template <typename F>
  decltype(auto) useJacobianhOver_q_dot(F&& f) {
    return siconos::algebra::visitStorage(jacobianhOver_q_dot_storage_, std::forward<F>(f),
                                          "jacobianhOver_q_dot_storage_");
  }

  /** function wrapper used to compute \f$ \frac{\partial}{\partial t}(J^h_q) \f$ */
  siconos::modeling::func_prototypes::FunctionBVBV_Ms computejacobianhOver_q_dot_{nullptr};
  siconos::modeling::func_prototypes::RFunctionBVBV_Ms computejacobianhOver_q_dot_python_{
      nullptr};

  /** \f$ \frac{\partial}{\partial t}(J^h_q) .\dot q\f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> jacobianhOver_q_dot_X_qdot_{nullptr};

 public:
  /** basic constructor */
  LagrangianSparseScleronomousR() : LagrangianSparseR(RelationSubType::ScleronomousR) {}

  /** destructor */
  virtual ~LagrangianSparseScleronomousR() noexcept = default;

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

  /** @brief Set a constant \f$ J^h_q \f$
   *  The input matrix is copied.
   *
   *  @param newValue jacobianhOver_q matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setConstantJacobianhOver_q(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    // Warning: we can't check dimensions --> they depend on the interaction
    jacobianhOver_q_storage_ =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
    hasJacobianhOver_q_ = true;
    hasConstantJacobianhOver_q_ = true;
    computejacobianhOver_q_ = nullptr;
    computejacobianhOver_q_python_ = nullptr;
  }

  /** @brief set a constant \f$ J^h_q \f$
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
  void setConstantJacobianhOver_q(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac,
                                  siconos::algebra::AliasTag);

  /** set a user-defined function to compute \f$ J^h_q \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  template <typename Func>
  void setComputeJacobianhOver_qFunction(Func&& f) {
    using Ret = std::invoke_result_t<Func, const siconos::algebra::BlockVector&,
                                     siconos::algebra::SiconosSparseMatrix&>;
    hasJacobianhOver_q_ = true;
    hasConstantJacobianhOver_q_ = false;

    if constexpr (std::is_same_v<Ret, void>) {
      computejacobianhOver_q_ = std::forward<Func>(f);
      computejacobianhOver_q_python_ = nullptr;
    } else {
      computejacobianhOver_q_python_ = std::forward<Func>(f);
      computejacobianhOver_q_ = nullptr;
    }
  }

  /** to compute  \f$ J^h_q(q)  \f$
   *  @param positions q vectors
   *  @param time the current time
   */
  void computeJacobianhOver_q(const siconos::algebra::BlockVector& positions);

  /** @brief Set a constant \f$ \frac{\partial}{\partial t}(J^h_q) \f$  \f$
   *  The input matrix is copied.
   *
   *  @param newValue jacobianhOver_q matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setConstantJacobianhOver_q_dot(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    // Warning: we can't check dimensions --> they depend on the interaction
    jacobianhOver_q_dot_storage_ =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
    hasJacobianhOver_q_dot_ = true;
    hasConstantJacobianhOver_q_dot_ = true;
    computejacobianhOver_q_dot_ = nullptr;
    computejacobianhOver_q_dot_python_ = nullptr;
  }

  /** @brief set a constant \f$ \frac{\partial}{\partial t}(J^h_q) \f$
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
  void setConstantJacobianhOver_q_dot(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac,
                                      siconos::algebra::AliasTag);

  /** set a user-defined function to compute \f$ \frac{\partial}{\partial t}(J^h_q) \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  template <typename Func>
  void setComputeJacobianhOver_q_dotFunction(Func&& f) {
    using Ret = std::invoke_result_t<Func, const siconos::algebra::BlockVector&,
                                     const siconos::algebra::BlockVector&,
                                     siconos::algebra::SiconosSparseMatrix&>;
    hasJacobianhOver_q_dot_ = true;
    hasConstantJacobianhOver_q_dot_ = false;

    if constexpr (std::is_same_v<Ret, void>) {
      computejacobianhOver_q_dot_ = std::forward<Func>(f);
      computejacobianhOver_q_dot_python_ = nullptr;
    } else {
      computejacobianhOver_q_dot_python_ = std::forward<Func>(f);
      computejacobianhOver_q_dot_ = nullptr;
    }
  }

  /** to compute  \f$ \frac{\partial}{\partial t}(J^h_q) \f$
   *  @param positions q vectors
   *  @param qdot time derivatives of the state vector
   */
  void computeJacobianhOver_q_dot(const siconos::algebra::BlockVector& positions,
                                  const siconos::algebra::BlockVector& qdot);

  /** @return a read-only view on \f$ \frac{\partial}{\partial t}( J^h_q(q)).\dot q\f$
   * vector */
  inline auto jacobianhOver_q_dot_X_qdot() const {
    return siconos::algebra::ConstMapVectorType(jacobianhOver_q_dot_X_qdot_->data(),
                                                jacobianhOver_q_dot_X_qdot_->size());
  }

  /** @return a read-only reference on \f$ \frac{\partial}{\partial t}(J^h_q) \f$ matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> jacobianhOver_q_dot() const {
    return useJacobianhOver_q_dot([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>(M);
    });
  }

  /*  @return a reference to the \f$ \frac{\partial}{\partial t}(J^h_q) \f$ matrix
   *
   *   pybind11 use only!
   */
  const siconos::algebra::SiconosSparseMatrix& jacobianhOver_q_dot_py() const {
    return siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
        jacobianhOver_q_dot_storage_,
        [](const auto& M) -> const siconos::algebra::SiconosSparseMatrix& { return M; },
        "jacobianhOver_q_dot_storage_");
  }

  /** @return True if  \f$ \frac{\partial}{\partial t}(J^h_q) \f$ matrix has been set */
  bool hasJacobianhOver_q_dot() const { return hasJacobianhOver_q_dot_; }

  /** to compute \f$ \frac{\partial}{\partial t}( J^h_q(q)).\dot q\f$
   *
   *  @param time double, current time
   *  @param inter interaction
   */
  void computeJacobianhOver_q_dot_X_qdot(double time, Interaction& inter);

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
};
}  // namespace siconos::modeling
#endif
