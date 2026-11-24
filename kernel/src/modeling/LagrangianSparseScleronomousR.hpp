// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
 *   - \f$ H(q) \f$ (setConstantContactMatrix, setComputeContactMatrixFunction);
 *   - \f$ \dot H \f$ (setConstantContactMatrixDot, setComputeContactMatrixDotFunction).
 *
 * Former names (deprecated, kept for one release): *JacobianhOver_q* = *ContactMatrix*,
 * *JacobianhOver_q_dot* = *contactMatrixDot*, jacobianhOver_q_dot_X_qdot = contactMatrixDotV.
 */
class LagrangianSparseScleronomousR : public LagrangianSparseR {
 protected:
  ACCEPT_SERIALIZATION(LagrangianSparseScleronomousR);

  /** function wrapper used to compute \f$ h(q) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_V computeh_{nullptr};

  /** function wrapper used to compute  \f$ H(q,t) \f$ */
  siconos::modeling::func_prototypes::FunctionBV_Ms computeContactMatrix_{nullptr};
  siconos::modeling::func_prototypes::RFunctionBV_Ms computeContactMatrix_python_{nullptr};

  /** \f$ \dot H \f$
   *  This value is useful to compute the second-order
   *  derivative of the constraints with respect to time.
   */
  siconos::algebra::SparseStorage contactMatrixDot_storage_{std::monostate{}};

  /** true if \f$ \dot H \f$ is not set or set and constant */
  bool hasConstantContactMatrixDot_{true};

  /** True if \f$ \dot H \f$ is required */
  bool hasContactMatrixDot_{false};

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
  decltype(auto) useContactMatrixDot(F&& f) const {
    return siconos::algebra::visitStorage(contactMatrixDot_storage_, std::forward<F>(f),
                                          "contactMatrixDot_storage_");
  }

  template <typename F>
  decltype(auto) useContactMatrixDot(F&& f) {
    return siconos::algebra::visitStorage(contactMatrixDot_storage_, std::forward<F>(f),
                                          "contactMatrixDot_storage_");
  }

  /** function wrapper used to compute \f$ \dot H \f$ */
  siconos::modeling::func_prototypes::FunctionBVBV_Ms computeContactMatrixDot_{nullptr};
  siconos::modeling::func_prototypes::RFunctionBVBV_Ms computeContactMatrixDot_python_{
      nullptr};

  /** \f$ \dot H .\dot q\f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> contactMatrixDotV_{nullptr};

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
    using Ret = std::invoke_result_t<Func, const siconos::algebra::BlockVector&,
                                     siconos::algebra::SiconosSparseMatrix&>;
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

  /** to compute  \f$ H(q)  \f$
   *  @param positions q vectors
   *  @param time the current time
   */
  void computeContactMatrix(const siconos::algebra::BlockVector& positions);

  /** @brief Set a constant \f$ \dot H \f$
   *  The input matrix is copied.
   *
   *  @param newValue H matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setConstantContactMatrixDot(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    // Warning: we can't check dimensions --> they depend on the interaction
    contactMatrixDot_storage_ =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
    hasContactMatrixDot_ = true;
    hasConstantContactMatrixDot_ = true;
    computeContactMatrixDot_ = nullptr;
    computeContactMatrixDot_python_ = nullptr;
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
  void setConstantContactMatrixDot(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac,
                                   siconos::algebra::AliasTag);

  /** set a user-defined function to compute \f$ \dot H \f$
   *
   *  @param fct the user-defined function (std::function, lambda ...)
   */
  template <typename Func>
  void setComputeContactMatrixDotFunction(Func&& f) {
    using Ret = std::invoke_result_t<Func, const siconos::algebra::BlockVector&,
                                     const siconos::algebra::BlockVector&,
                                     siconos::algebra::SiconosSparseMatrix&>;
    hasContactMatrixDot_ = true;
    hasConstantContactMatrixDot_ = false;

    if constexpr (std::is_same_v<Ret, void>) {
      computeContactMatrixDot_ = std::forward<Func>(f);
      computeContactMatrixDot_python_ = nullptr;
    } else {
      computeContactMatrixDot_python_ = std::forward<Func>(f);
      computeContactMatrixDot_ = nullptr;
    }
  }

  /** to compute  \f$ \dot H \f$
   *  @param positions q vectors
   *  @param qdot time derivatives of the state vector
   */
  void computeContactMatrixDot(const siconos::algebra::BlockVector& positions,
                               const siconos::algebra::BlockVector& qdot);

  /** @return a read-only view on \f$ \dot H \dot q\f$
   * vector */
  inline auto contactMatrixDotV() const {
    return siconos::algebra::ConstMapVectorType(contactMatrixDotV_->data(),
                                                contactMatrixDotV_->size());
  }

  /** @return a read-only reference on \f$ \dot H \f$ matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> contactMatrixDot() const {
    return useContactMatrixDot([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>(M);
    });
  }

  /*  @return a reference to the \f$ \dot H \f$ matrix
   *
   *   pybind11 use only!
   */
  const siconos::algebra::SiconosSparseMatrix& contactMatrixDot_py() const {
    return siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
        contactMatrixDot_storage_,
        [](const auto& M) -> const siconos::algebra::SiconosSparseMatrix& { return M; },
        "contactMatrixDot_storage_");
  }

  /** @return True if  \f$ \dot H \f$ matrix has been set */
  bool hasContactMatrixDot() const { return hasContactMatrixDot_; }

  /** to compute \f$ \dot H \dot q\f$
   *
   *  @param time double, current time
   *  @param inter interaction
   */
  void computeContactMatrixDotV(double time, Interaction& inter);

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
      const siconos::algebra::BlockVector& positions) {
    computeContactMatrix(positions);
  }
  template <typename T>
  [[deprecated("use setConstantContactMatrixDot")]] void setConstantJacobianhOver_q_dot(
      T&& newValue, siconos::algebra::CopyTag tag) {
    setConstantContactMatrixDot(std::forward<T>(newValue), tag);
  }
  [[deprecated("use setConstantContactMatrixDot")]] void setConstantJacobianhOver_q_dot(
      Eigen::Map<siconos::algebra::SiconosSparseMatrix>& jac, siconos::algebra::AliasTag tag) {
    setConstantContactMatrixDot(jac, tag);
  }
  template <typename Func>
  [[deprecated("use setComputeContactMatrixDotFunction")]] void
  setComputeJacobianhOver_q_dotFunction(Func&& f) {
    setComputeContactMatrixDotFunction(std::forward<Func>(f));
  }
  [[deprecated("use computeContactMatrixDot")]] void computeJacobianhOver_q_dot(
      const siconos::algebra::BlockVector& positions,
      const siconos::algebra::BlockVector& qdot) {
    computeContactMatrixDot(positions, qdot);
  }
  [[deprecated("use contactMatrixDotV()")]] auto jacobianhOver_q_dot_X_qdot() const {
    return contactMatrixDotV();
  }
  [[deprecated(
      "use contactMatrixDot()")]] Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>
  jacobianhOver_q_dot() const {
    return contactMatrixDot();
  }
  [[deprecated("use contactMatrixDot_py()")]] const siconos::algebra::SiconosSparseMatrix&
  jacobianhOver_q_dot_py() const {
    return contactMatrixDot_py();
  }
  [[deprecated("use hasContactMatrixDot()")]] bool hasJacobianhOver_q_dot() const {
    return hasContactMatrixDot();
  }
  [[deprecated("use computeContactMatrixDotV")]] void computeJacobianhOver_q_dot_X_qdot(
      double time, Interaction& inter) {
    computeContactMatrixDotV(time, inter);
  }
  ///@}
};
}  // namespace siconos::modeling
#endif
