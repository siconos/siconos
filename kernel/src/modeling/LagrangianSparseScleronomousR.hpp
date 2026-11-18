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
 * \f[
 * \begin{aligned}
 *   y &= h(q), \\
 *   \dot y &= H(q)\, \dot q, \qquad
 *     H(q) = \nabla_q^\top h(q)
 *     = \left[ \frac{\partial h_i}{\partial q_j} \right]_{1\le i\le m,\;1\le j\le n}, \\
 *   p &= H(q)^\top \lambda.
 * \end{aligned}
 * \f]
 *
 * Sparse-storage counterpart of LagrangianScleronomousR (same notation, see LagrangianR).
 * The following operators can be set by user-defined functions:
 *
 *   - \f$ h(q) \f$ (setComputehFunction);
 *   - \f$ H(q) \f$ (setConstantH, setComputeHFunction);
 *   - \f$ \dot H \f$ (setConstantHdot, setComputeHdotFunction).
 *
 * Former names (deprecated, kept for one release): *JacobianhOver_q* = *H*,
 * *JacobianhOver_q_dot* = *Hdot*, jacobianhOver_q_dot_X_qdot = Hdotv.
 */
class LagrangianSparseScleronomousR : public LagrangianSparseR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianSparseScleronomousR);

  /** function wrapper used to compute \f$ h(q) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_V computeh_{nullptr};

  /** function wrapper used to compute  \f$ H(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_Ms computeH_{nullptr};
  siconos::modeling::func_prototypes::RFunctionBV_Ms computeH_python_{nullptr};

  /** \f$ \dot H \f$
   *  This value is useful to compute the second-order
   *  derivative of the constraints with respect to time.
   */
  siconos::algebra::SparseStorage Hdot_storage_{std::monostate{}};

  /** true if \f$ \dot H \f$ is not set or set and constant */
  bool hasConstantHdot_{true};

  /** True if \f$ \dot H \f$ is required */
  bool hasHdot_{false};

  /**
   * @brief Utility function providing uniform access to the \f$ \dot H \f$ matrix.
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
  decltype(auto) useHdot(F&& f) const {
    return siconos::algebra::visitStorage(Hdot_storage_, std::forward<F>(f), "Hdot_storage_");
  }

  template <typename F>
  decltype(auto) useHdot(F&& f) {
    return siconos::algebra::visitStorage(Hdot_storage_, std::forward<F>(f), "Hdot_storage_");
  }

  /** function wrapper used to compute \f$ \dot H \f$ */
  siconos::modeling::func_prototypes::FunctionBVBV_Ms computeHdot_{nullptr};
  siconos::modeling::func_prototypes::RFunctionBVBV_Ms computeHdot_python_{nullptr};

  /** \f$ \dot H .\dot q\f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> Hdotv_{nullptr};

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

  /** @brief Set a constant \f$ H \f$
   *  The input matrix is copied.
   *
   *  @param newValue H matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setConstantH(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    // Warning: we can't check dimensions --> they depend on the interaction
    H_storage_ =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
    hasH_ = true;
    hasConstantH_ = true;
    computeH_ = nullptr;
    computeH_python_ = nullptr;
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
  void setConstantH(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac,
                    siconos::algebra::AliasTag);

  /** set a user-defined function to compute \f$ H \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  template <typename Func>
  void setComputeHFunction(Func&& f) {
    using Ret = std::invoke_result_t<Func, const siconos::algebra::BlockVector&,
                                     siconos::algebra::SiconosSparseMatrix&>;
    hasH_ = true;
    hasConstantH_ = false;

    if constexpr (std::is_same_v<Ret, void>) {
      computeH_ = std::forward<Func>(f);
      computeH_python_ = nullptr;
    } else {
      computeH_python_ = std::forward<Func>(f);
      computeH_ = nullptr;
    }
  }

  /** to compute  \f$ H(q)  \f$
   *  @param positions q vectors
   *  @param time the current time
   */
  void computeH(const siconos::algebra::BlockVector& positions);

  /** @brief Set a constant \f$ \dot H \f$
   *  The input matrix is copied.
   *
   *  @param newValue H matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setConstantHdot(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    // Warning: we can't check dimensions --> they depend on the interaction
    Hdot_storage_ =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
    hasHdot_ = true;
    hasConstantHdot_ = true;
    computeHdot_ = nullptr;
    computeHdot_python_ = nullptr;
  }

  /** @brief set a constant \f$ \dot H \f$
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
  void setConstantHdot(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac,
                       siconos::algebra::AliasTag);

  /** set a user-defined function to compute \f$ \dot H \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  template <typename Func>
  void setComputeHdotFunction(Func&& f) {
    using Ret = std::invoke_result_t<Func, const siconos::algebra::BlockVector&,
                                     const siconos::algebra::BlockVector&,
                                     siconos::algebra::SiconosSparseMatrix&>;
    hasHdot_ = true;
    hasConstantHdot_ = false;

    if constexpr (std::is_same_v<Ret, void>) {
      computeHdot_ = std::forward<Func>(f);
      computeHdot_python_ = nullptr;
    } else {
      computeHdot_python_ = std::forward<Func>(f);
      computeHdot_ = nullptr;
    }
  }

  /** to compute  \f$ \dot H \f$
   *  @param positions q vectors
   *  @param qdot time derivatives of the state vector
   */
  void computeHdot(const siconos::algebra::BlockVector& positions,
                   const siconos::algebra::BlockVector& qdot);

  /** @return a read-only view on \f$ \dot H \dot q\f$
   * vector */
  inline auto Hdotv() const {
    return siconos::algebra::ConstMapVectorType(Hdotv_->data(), Hdotv_->size());
  }

  /** @return a read-only reference on \f$ \dot H \f$ matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> Hdot() const {
    return useHdot([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>(M);
    });
  }

  /*  @return a reference to the \f$ \dot H \f$ matrix
   *
   *   pybind11 use only!
   */
  const siconos::algebra::SiconosSparseMatrix& Hdot_py() const {
    return siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
        Hdot_storage_,
        [](const auto& M) -> const siconos::algebra::SiconosSparseMatrix& { return M; },
        "Hdot_storage_");
  }

  /** @return True if  \f$ \dot H \f$ matrix has been set */
  bool hasHdot() const { return hasHdot_; }

  /** to compute \f$ \dot H \dot q\f$
   *
   *  @param time double, current time
   *  @param inter interaction
   */
  void computeHdotv(double time, Interaction& inter);

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
  [[deprecated("use setConstantH")]] void setConstantJacobianhOver_q(
      T&& newValue, siconos::algebra::CopyTag tag) {
    setConstantH(std::forward<T>(newValue), tag);
  }
  [[deprecated("use setConstantH")]] void setConstantJacobianhOver_q(
      Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac, siconos::algebra::AliasTag tag) {
    setConstantH(jac, tag);
  }
  template <typename Func>
  [[deprecated("use setComputeHFunction")]] void setComputeJacobianhOver_qFunction(Func&& f) {
    setComputeHFunction(std::forward<Func>(f));
  }
  [[deprecated("use computeH")]] void computeJacobianhOver_q(
      const siconos::algebra::BlockVector& positions) {
    computeH(positions);
  }
  template <typename T>
  [[deprecated("use setConstantHdot")]] void setConstantJacobianhOver_q_dot(
      T&& newValue, siconos::algebra::CopyTag tag) {
    setConstantHdot(std::forward<T>(newValue), tag);
  }
  [[deprecated("use setConstantHdot")]] void setConstantJacobianhOver_q_dot(
      Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac, siconos::algebra::AliasTag tag) {
    setConstantHdot(jac, tag);
  }
  template <typename Func>
  [[deprecated("use setComputeHdotFunction")]] void setComputeJacobianhOver_q_dotFunction(
      Func&& f) {
    setComputeHdotFunction(std::forward<Func>(f));
  }
  [[deprecated("use computeHdot")]] void computeJacobianhOver_q_dot(
      const siconos::algebra::BlockVector& positions,
      const siconos::algebra::BlockVector& qdot) {
    computeHdot(positions, qdot);
  }
  [[deprecated("use Hdotv()")]] auto jacobianhOver_q_dot_X_qdot() const { return Hdotv(); }
  [[deprecated("use Hdot()")]] Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>
  jacobianhOver_q_dot() const {
    return Hdot();
  }
  [[deprecated("use Hdot_py()")]] const siconos::algebra::SiconosSparseMatrix&
  jacobianhOver_q_dot_py() const {
    return Hdot_py();
  }
  [[deprecated("use hasHdot()")]] bool hasJacobianhOver_q_dot() const { return hasHdot(); }
  [[deprecated("use computeHdotv")]] void computeJacobianhOver_q_dot_X_qdot(
      double time, Interaction& inter) {
    computeHdotv(time, inter);
  }
  ///@}
};
}  // namespace siconos::modeling
#endif
