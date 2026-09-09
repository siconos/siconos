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

/** @file LagrangianSparseR.hpp
 */

#ifndef LAGRANGIANSPARSERELATION_H
#define LAGRANGIANSPARSERELATION_H

#include "Relation.hpp"
#include "StorageTools.hpp"

namespace siconos::modeling {
/**
   @brief interface to 2nd order Lagrangian nonlinear relations (sparse storage)

   \f$ y = h(t,q,\dot q,\ldots) \f$  describes the constraint

  - Sparse storage for matrix-like attributes (see LagrangianSparseR for sparse storage)
  - h may be \f$ h(q)\f$  (scleronomous), \f$ h(q,t)\f$  (rheonomous), \f$ h(q, \lambda)\f$
  (compliant)... Choose the proper derivated class, according to your needs.

  - The Jacobian of the constraints with respect to the coordinates  \f$ q \f$ is always
  defined. Other jacobians are defined only when required in the derived classes.

  Notations:

  \f[
    J^h_q(q, ...)
    = \frac{\partial h}{\partial q} \in \mathbb{R}^{M \times N}
    =
   \left[
   \frac{\partial h_i}{\partial q_j}
   \right]_{1\le i\le M,\;1\le j\le N}.
  \f]

  with \f$ q \in \mathbb{R}^N \f$ and \f$ h \in \mathbb{R}^M \f$

  Same for other jacobians (\f$ J^h_\lambda(q, ...) \f$ ) when they exist.

  This Jacobians are mainly used for Newton linearization and to compute the time-derivative of
  the constraint, \f$ y = h(q,\ldots) \f$

  that is  \f$ \dot y (t) = J^h_q(q, ...)\dot q +\ldots \f$

  This object can also store more general linearized part of the gap function.
  If \f$ y=h(q) \f$ models a gap function, then the time derivative can be generically  written
  as

  \f$ \dot y (t) = H(q,\ldots) \dot q  +\ldots.  \f$

  All operators of Lagrangian relations can be set by user-defined functions.
  The signature of the plugged functions depends on which variables are taken into
  account in h and its jacobians. Here again, check in the proper derived class API to find
  which functions are available.


 */
class LagrangianSparseR : public Relation {
 public:
  /** Names and positions of ds variables that might be used (read only) in compute input and
   * output.
   * Last param, size, gives the size of the BlockVector to be allocated.
   */
  enum class ds_var : std::size_t { q0, q1, q2, p0, p1, p2, size };

 protected:
  ACCEPT_SERIALIZATION(LagrangianSparseR);

  /** The Jacobian of the constraints with respect to the generalized coordinates   \f$ q \f$
   */
  siconos::algebra::SparseStorage jacobianhOver_q_storage_{std::monostate{}};

  /** true if jacobianhOver_q is not set or set and constant */
  bool hasConstantJacobianhOver_q_{true};

  /** True if jacobianhOver_q is required */
  bool hasJacobianhOver_q_{false};

  /**
   * @brief Utility function providing uniform access to the  \f$ J^h_q \f$ matrix.
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
  decltype(auto) useJacobianhOver_q(F&& f) const {
    return siconos::algebra::visitStorage(jacobianhOver_q_storage_, std::forward<F>(f),
                                          "jacobianhOver_q_storage_");
  }

  template <typename F>
  decltype(auto) useJacobianhOver_q(F&& f) {
    return siconos::algebra::visitStorage(jacobianhOver_q_storage_, std::forward<F>(f),
                                          "jacobianhOver_q_storage_");
  }

  /** @brief minimal (protected) constructor
   *
   *  @param lagType the sub-type of the relation
   */
  LagrangianSparseR(RelationSubType lagType)
      : Relation(RelationType::LagrangianSparse, lagType) {}

 public:
  /** @brief destructor */
  virtual ~LagrangianSparseR() noexcept = default;

  /**
   *  @brief initialize the relation (check sizes, memory allocation ...)
   *
   *  @param inter the interaction using this relation
   */
  inline void initialize(Interaction& inter) override {};

  /** check sizes of the relation specific operators.
   *
   *  @param inter an Interaction using this relation
   */
  virtual void checkSize(const Interaction& inter) const override {};

  /** @return a read-only reference on \f$ J^h_q \f$ matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> jacobianhOver_q() const {
    return useJacobianhOver_q([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>(M);
    });
  }

  /*  @return a reference to the \f$ J^h_q  \f$ matrix
   *
   *   pybind11 use only!
   */
  const siconos::algebra::SiconosSparseMatrix& jacobianhOver_q_py() const {
    return siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
        jacobianhOver_q_storage_,
        [](const auto& M) -> const siconos::algebra::SiconosSparseMatrix& { return M; },
        "jacobianhOver_q_storage_");
  }

  /** @return True if  \f$ J^h_q \f$ matrix has been set */
  bool hasJacobianhOver_q() const { return hasJacobianhOver_q_; }

  /** main relation members display */
  void display() const override;

  /** @brief allocation of memory space for relations in the graph
   *
   *      Warning: internal use only (called from Topology)
   *  @param ds_vars a container of vectors (pointers), from the parent interaction
   *  @param ds1 first ds concerned by the relation
   *  @param ds2 second ds concerned by the relation
   */
  virtual void allocate_read_dynamical_systems_var_vectors(
      std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
      modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const override;

  virtual void accept(relations::Visitor& tourist) const override { tourist.visit(*this); }
};

}  // namespace siconos::modeling

#endif  // LAGRANGIANSPARSERELATION_H
