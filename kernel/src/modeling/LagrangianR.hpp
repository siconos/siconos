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

/** @file LagrangianR.hpp
 *  @brief Interface to the Lagrangian (2nd order) relations (dense storage).
 */
#ifndef LAGRANGIANRELATION_H
#define LAGRANGIANRELATION_H

#include "Relation.hpp"
#include "RelationCompatibility.hpp"

namespace siconos::modeling {
/**
  @brief Interface to 2nd order Lagrangian relations (dense storage)

  The dynamical systems are Lagrangian, with coordinates \f$ q \in \mathbb{R}^n \f$
  (stacked for the dynamical systems of the interaction) and velocities
  \f$ v = \dot q \f$. The relation, of size \f$ m \f$, reads
  \f[
  \begin{aligned}
    y &= h(q, \ldots) \in \mathbb{R}^m, \\
    \dot y &= H(q, \ldots)\, v + \ldots, \\
    p &= H(q, \ldots)^\top \lambda,
  \end{aligned}
  \f]
  where \f$ H \in \mathbb{R}^{m \times n} \f$ (H()) is the <em>contact configuration
  matrix</em>: it maps the generalized velocities to the local (relative) velocities of the
  relation. When the relation is given by a function \f$ h \f$, \f$ H \f$ is its Jacobian,
  \f[
    H(q, \ldots) = G(q, \ldots) = \nabla_q^\top h(q, \ldots)
    = \left[ \frac{\partial h_i}{\partial q_j} \right]_{1\le i\le m,\;1\le j\le n},
  \f]
  since \f$ \dot q = v \f$ (\f$ T = I \f$ with the notation of NewtonEulerR); G() returns
  the same matrix as H(). For contact relations with friction or rolling resistance (e.g.
  Lagrangian2dVelocityR), the rows of \f$ H \f$ beyond the normal one are not the Jacobian of
  a function: they contain the change to the local contact frame and the relative velocity
  between the bodies.

  - h may be \f$ h(q) \f$ (scleronomic, LagrangianScleronomousR), \f$ h(q,t) \f$
    (rheonomic, LagrangianRheonomousR), \f$ h(q, \lambda) \f$ (compliant,
    LagrangianCompliantR), or linear (LagrangianLinearTIR, LagrangianCompliantLinearTIR).
  - Dense storage for matrix-like attributes (see LagrangianSparseR for sparse storage).
  - All operators can be set by user-defined functions; their signatures depend on the
    derived class.

  <b>Former names</b> (deprecated, kept for one release): jacobianhOver_q() = H(), and in
  the derived classes, *JacobianhOver_q* = *H*, *jacobianhOver_q_dot* = *Hdot*,
  jacobianhOver_q_dot_X_qdot = Hdotv.
 */
class LagrangianR : public Relation {
 public:
  /** Names and positions of ds variables that might be used (read only) in compute input and
   * output.
   * Last param, size, gives the size of the BlockVector to be allocated.
   */
  enum class ds_var : std::size_t { q0, q1, q2, p0, p1, p2, size };

 protected:
  ACCEPT_SERIALIZATION(LagrangianR);

  /** View onto the contact configuration matrix \f$ H \in \mathbb{R}^{m \times n} \f$
   *  (either H_storage_ or external memory) */
  std::shared_ptr<siconos::algebra::MapType> H_view_{nullptr};

  /** internal (optional) storage used for \f$ H \f$ */
  std::unique_ptr<std::vector<double>> H_storage_{nullptr};

  /** True if \f$ H \f$ is a constant matrix */
  bool hasConstantH_{false};

  /** Deprecated names of the members above, kept for one release. */
  [[deprecated(
      "use H_view_")]] std::shared_ptr<siconos::algebra::MapType>& jacobianhOver_q_view_{
      H_view_};
  [[deprecated("use H_storage_")]] std::unique_ptr<std::vector<double>>&
      jacobianhOver_q_internal_storage_{H_storage_};
  [[deprecated("use hasConstantH_")]] bool& hasConstantJacobianhOver_q_{hasConstantH_};

  /** @brief minimal (protected) constructor
   *
   *  @param lagType the sub-type of the relation
   */
  LagrangianR(RelationSubType lagType) : Relation(RelationType::Lagrangian, lagType) {}

 public:
  /** @brief destructor */
  virtual ~LagrangianR() noexcept = default;

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

  /** @return a read-only view on the contact configuration matrix
   *  \f$ H \in \mathbb{R}^{m \times n} \f$, \f$ \dot y = H v + \ldots \f$ */
  inline auto H() const {
    return siconos::algebra::ConstMapType(H_view_->data(), H_view_->rows(), H_view_->cols());
  }

  /** @return a read-only view on \f$ G = \nabla_q^\top h \f$, which is the same matrix as
   *  H() for Lagrangian relations (\f$ \dot q = v \f$) */
  inline auto G() const { return H(); }

  /** @deprecated use H() */
  [[deprecated("use H()")]] inline auto jacobianhOver_q() const { return H(); }

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

#endif  // LAGRANGIANRELATION_H
