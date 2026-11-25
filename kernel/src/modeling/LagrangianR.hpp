// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

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
  where \f$ H \in \mathbb{R}^{m \times n} \f$ (contactMatrix()) is the <em>contact
  configuration matrix</em>: it maps the generalized velocities to the local (relative)
  velocities of the relation. When the relation is given by a function \f$ h \f$, \f$ H \f$ is
  its Jacobian,
  \f[
    H(q, \ldots) = G(q, \ldots) = \nabla_q^\top h(q, \ldots)
    = \left[ \frac{\partial h_i}{\partial q_j} \right]_{1\le i\le m,\;1\le j\le n},
  \f]
  since \f$ \dot q = v \f$ (\f$ T = I \f$ with the notation of NewtonEulerR);
  constraintJacobian() returns the same matrix as contactMatrix(). For contact relations with
  friction or rolling resistance (e.g. Lagrangian2dVelocityR), the rows of \f$ H \f$ beyond the
  normal one are not the Jacobian of a function: they contain the change to the local contact
  frame and the relative velocity between the bodies.

  - h may be \f$ h(q) \f$ (scleronomic, LagrangianScleronomousR), \f$ h(q,t) \f$
    (rheonomic, LagrangianRheonomousR), \f$ h(q, \lambda) \f$ (compliant,
    LagrangianCompliantR), or linear (LagrangianLinearTIR, LagrangianCompliantLinearTIR).
  - Dense storage for matrix-like attributes (see LagrangianSparseR for sparse storage).
  - All operators can be set by user-defined functions; their signatures depend on the
    derived class.

  <b>Names in the code.</b> The matrices are named after what they are, not after their
  letter: \f$ H \f$ is contactMatrix() (setConstantContactMatrix, computeContactMatrix,
  setComputeContactMatrixFunction, ...), \f$ G \f$ is constraintJacobian(),
  \f$ \dot H \f$ is contactMatrixDot() and \f$ \dot H v \f$ is contactMatrixDotV().

  <b>Former names</b> (deprecated, kept for one release): jacobianhOver_q() = contactMatrix(),
  and in the derived classes, *JacobianhOver_q* = *ContactMatrix*, *jacobianhOver_q_dot* =
  *contactMatrixDot*, jacobianhOver_q_dot_X_qdot = contactMatrixDotV.
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
   *  (either contactMatrix_storage_ or external memory) */
  std::shared_ptr<siconos::algebra::MapType> contactMatrix_view_{nullptr};

  /** internal (optional) storage used for \f$ H \f$ */
  std::unique_ptr<std::vector<double>> contactMatrix_storage_{nullptr};

  /** True if \f$ H \f$ is a constant matrix */
  bool hasConstantContactMatrix_{false};

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
  inline auto contactMatrix() const {
    return siconos::algebra::ConstMapType(
        contactMatrix_view_->data(), contactMatrix_view_->rows(), contactMatrix_view_->cols());
  }

  /** @return a read-only view on \f$ G = \nabla_q^\top h \f$, which is the same matrix as
   *  contactMatrix() for Lagrangian relations (\f$ \dot q = v \f$) */
  inline auto constraintJacobian() const { return contactMatrix(); }

  /** @deprecated use contactMatrix() */
  [[deprecated("use contactMatrix()")]] inline auto jacobianhOver_q() const {
    return contactMatrix();
  }

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
