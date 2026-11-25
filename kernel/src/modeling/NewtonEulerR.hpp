/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/*! \file NewtonEulerR.hpp
 *  \brief Generic relation for Newton-Euler dynamical systems.
 */
#ifndef NEWTONEULERR_HPP
#define NEWTONEULERR_HPP

#include <cassert>
#include <optional>

#include "Relation.hpp"
#include "RelationCompatibility.hpp"

namespace siconos::modeling {

class DynamicalSystem;

/** Generic relation for Newton-Euler dynamical systems (one or two bodies).

  <b>Kinematics.</b> Body \f$ i \in \{1, 2\} \f$ is described by its coordinates
  \f$ q_i = (x_{G_i}, p_i) \in \mathbb{R}^7 \f$ (position of the centre of mass and unit
  quaternion) and its twist \f$ v_i = (v_{G_i}, \Omega_i) \in \mathbb{R}^6 \f$, where
  \f$ v_{G_i} \f$ is expressed in the inertial frame and the angular velocity
  \f$ \Omega_i \f$ in the body frame. They are linked by
  \f[
    \dot q_i = T(q_i)\, v_i, \qquad
    T(q_i) = \begin{bmatrix} I_3 & 0 \\ 0 & \Psi(p_i) \end{bmatrix} \in \mathbb{R}^{7\times 6},
  \f]
  with \f$ \dot p_i = \Psi(p_i)\,\Omega_i = \frac{1}{2}\, p_i \otimes (0, \Omega_i) \f$
  (see newton_euler::computeT). For an interaction with \f$ n \in \{1, 2\} \f$ bodies,
  \f$ q = (q_1, q_2) \in \mathbb{R}^{7n} \f$ and \f$ v = (v_1, v_2) \in \mathbb{R}^{6n} \f$.

  <b>Relation.</b> With \f$ m \f$ the size of the interaction,
  \f[
  \begin{aligned}
    y &= h(q) \in \mathbb{R}^m, \\
    \dot y &= H(q)\, v, \\
    \ddot y &= H(q)\, \dot v + \dot H(q)\, v,
  \end{aligned}
  \f]
  where \f$ H(q) \in \mathbb{R}^{m \times 6n} \f$ (contactMatrix()) is the <em>contact
  configuration matrix</em>: it maps the twists to the local (relative) velocities of the
  relation. For a relation given by a function \f$ h(q) \f$, it is
  \f[
    H(q) = G(q)\, T(q), \qquad G(q) = \nabla_q^\top h(q) \in \mathbb{R}^{m \times 7n}
    \quad (\text{constraintJacobian()}),
  \f]
  and \f$ \dot H v = (\dot G\, T + G\, \dot T)\, v = \dot G\, \dot q + G\, \dot T v \f$
  (contactMatrixDotV()). For contact relations with friction or rolling resistance, the rows of
  \f$ H \f$ beyond the normal one are not the Jacobian of a function: they also contain the
  change to the local contact frame and the relative velocity between the bodies (see
  NewtonEulerNormalVelocityR). \f$ G \f$ is then only meaningful for the normal row.

  The generalized forces (impulses) applied to the bodies are, for the input level
  \f$ k \f$,
  \f[
    p_k = H(q)^\top \lambda_k \quad (k = 1, 2), \qquad
    p_0 = G(q)^\top \lambda_0 \quad \text{(projection onto the constraints)}.
  \f]
  Only scleronomic relations are handled: \f$ h \f$ does not depend explicitly on time.

  <b>Default implementation.</b> The relation is linear,
  \f$ h(q) = G\, q + e \f$, with a constant \f$ G \f$ (setConstantConstraintJacobian) and an
  optional constant \f$ e \f$ (setConstanteVector). \f$ H = G\, T(q) \f$ is recomputed from
  \f$ T(q) \f$ in computeContactMatrix.

  <b>Derived classes</b> override:
   - computeh() to compute \f$ y = h(q) \f$;
   - computeConstraintJacobian() to update \f$ G(q) \f$ (only used when the matrix is not
  constant);
   - computeContactMatrix() when \f$ H(q) \f$ is computed directly, without forming \f$ G \f$
  (e.g. contact relations written in a local frame, see NewtonEulerNormalVelocityR).

  \f$ \dot G \f$, needed for \f$ \ddot y \f$, can be given by a user-defined function
  (setComputeConstraintJacobianDotFunction).

  <b>Storage.</b>
   - \f$ G \f$: a shared view on external memory after setConstantConstraintJacobian, internal
  storage otherwise (allocated in initialize, or by the derived class).
   - \f$ H \f$ (contactMatrix_) and the buffer for \f$ T \f$: always allocated in initialize.
   - \f$ e \f$: a shared view on external memory, only after setConstanteVector.
   - \f$ \dot G \f$: allocated in initialize if a function has been set, otherwise at the
     first call of computeConstraintJacobianDot or computeContactMatrixDotV.

  <b>Names in the code.</b> \f$ H \f$ is contactMatrix() (setConstantContactMatrix is not
  available: H depends on T(q)), \f$ G \f$ is constraintJacobian()
  (setConstantConstraintJacobian, computeConstraintJacobian), \f$ \dot G \f$ is
  constraintJacobianDot (computeConstraintJacobianDot, setComputeConstraintJacobianDotFunction)
  and \f$ \dot H v \f$ is contactMatrixDotV().

  <b>Former names</b> (deprecated, kept for one release): H_NE() = constraintJacobian(),
  H_NE_prod_T() = contactMatrix(), setConstantH_NE = setConstantConstraintJacobian,
  computeH_NE_ = computeConstraintJacobian, computeH_NE_prod_T = computeContactMatrix,
  setComputeH_NE_dotFunction = setComputeConstraintJacobianDotFunction, computeH_NE_dot =
  computeConstraintJacobianDot, secondOrderTimeDerivativeTerms() = contactMatrixDotV(),
  computeSecondOrderTimeDerivativeTerms = computeContactMatrixDotV. A derived class that still
  overrides computeH_NE_ or computeH_NE_prod_T keeps being called.

  \see DevNotes, chapter "NewtonEuler Dynamical Systems", section "Mechanical systems
  with bilateral and unilateral constraints".
 */
class NewtonEulerR : public Relation {
 public:
  /** Names and positions of the ds variables read in computeInput and computeOutput:
   *  q0 = \f$ q \f$, velocity = \f$ v \f$ (twist), dotq = \f$ \dot q \f$, and
   *  p0, p1, p2 = \f$ p_k \f$. The last value, size, is the number of variables.
   */
  enum class ds_var : std::size_t { q0, velocity, dotq, p0, p1, p2, size };

 protected:
  ACCEPT_SERIALIZATION(NewtonEulerR);

  /** View onto \f$ G = \nabla_q^\top h \in \mathbb{R}^{m \times 7n} \f$ (either
   *  constraintJacobian_storage_ or external memory, see setConstantConstraintJacobian) */
  std::shared_ptr<siconos::algebra::MapType> constraintJacobian_view_{nullptr};

  /** internal (optional) storage used for \f$ G \f$ */
  siconos::algebra::SiconosMatrix constraintJacobian_storage_{0, 0};

  /** True if \f$ G \f$ is a constant matrix (set with setConstantConstraintJacobian) */
  bool hasConstantConstraintJacobian_{false};

  /** \f$ \dot G = \frac{d}{dt} G(q) \in \mathbb{R}^{m \times 7n} \f$, used to compute
   *  \f$ \dot H v \f$. */
  std::shared_ptr<siconos::algebra::SiconosMatrix> constraintJacobianDot_{nullptr};

  /** user-defined function computing \f$ \dot G(q, \dot q) \f$ */
  siconos::modeling::func_prototypes::FunctionBVBV_M computeConstraintJacobianDot_{nullptr};

  /** View onto the constant vector \f$ e \f$ of the default linear relation */
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

  /** True if \f$ e \f$ is taken into account */
  bool haseVector_{false};

  /** Generalized forces (or impulses) applied to the bodies by the relation,
   *  \f$ H^\top \lambda_k \in \mathbb{R}^{6n} \f$, updated by computeInput for
   *  \f$ k = 1, 2 \f$. Useful for the end user. */
  siconos::algebra::SiconosVector contactForce_;

  /** \f$ H \in \mathbb{R}^{m \times 6n} \f$, the contact configuration matrix:
   *  \f$ \dot y = H v \f$. Computed as \f$ G T \f$ by default, or directly from the contact
   *  geometry by the contact relations. */
  std::shared_ptr<siconos::algebra::SiconosMatrix> contactMatrix_{nullptr};

  /** buffer for \f$ T(q_i) = \begin{bmatrix} I_3 & 0 \\ 0 & \Psi(p_i) \end{bmatrix} \f$
   *  of one body */
  siconos::algebra::SiconosMatrix76 T_buffer_;

  /** \f$ \dot H v = \dot G\, \dot q + G\, \dot T v \f$: the terms of \f$ \ddot y \f$ that
   *  do not depend on \f$ \dot v \f$. */
  siconos::algebra::SiconosVector contactMatrixDotV_;

  /** allocate \f$ \dot G \f$ (constraintJacobianDot_), with the size of \f$ G \f$, if not done
   */
  void allocateConstraintJacobianDot_();

  /** update \f$ G(q) \f$ (constraintJacobian_view_). Called by computeJach unless the matrix
   * is constant. The default does nothing (constant \f$ G \f$). Override this function.
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  coordinates \f$ q \f$ of the bodies of the interaction
   */
  virtual void computeConstraintJacobian(double time, Interaction& inter,
                                         const siconos::algebra::BlockVector& q0);

  /** \deprecated former name of computeConstraintJacobian. Still called for classes that
   * override it. */
  [[deprecated("override computeConstraintJacobian")]] virtual void computeH_NE_(
      double time, Interaction& inter, const siconos::algebra::BlockVector& q0);

  /** \f$ H = G T \f$, body by body: the default of computeContactMatrix. */
  void computeContactMatrixFromConstraintJacobian_(const siconos::algebra::BlockVector& q0);

 private:
  /** flags of the dispatch to the deprecated virtual names (see RelationCompatibility.hpp) */
  bool computeConstraintJacobianLegacy_{false};
  bool computeContactMatrixLegacy_{false};
  bool computeConstraintJacobianDotLegacy_{false};

 public:
  /** Default and only constructor */
  NewtonEulerR() : Relation(RelationType::NewtonEuler, RelationSubType::NonLinearR) {};

  /** destructor */
  virtual ~NewtonEulerR() noexcept = default;

  /** initialize components specific to derived classes.
   *
   *  \param inter  Interaction associated with the Relation
   */
  virtual void initialize(Interaction& inter) override;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  virtual void checkSize(const Interaction& inter) const override;

  /** \return a read-only view on \f$ G = \nabla_q^\top h \in \mathbb{R}^{m \times 7n} \f$ */
  inline auto constraintJacobian() const {
    return siconos::algebra::ConstMapType(constraintJacobian_view_->data(),
                                          constraintJacobian_view_->rows(),
                                          constraintJacobian_view_->cols());
  }

  /** Set a constant \f$ G \f$ for the default linear relation \f$ h(q) = G q + e \f$
   *  (warning: shared memory, newValue must outlive the relation).
   *
   *  \param newValue \f$ G \in \mathbb{R}^{m \times 7n} \f$
   */
  void setConstantConstraintJacobian(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** \return a read-only view on \f$ e \f$. Requires haseVector(). */
  inline auto eVector() const {
    assert(eVector_view_ && "NewtonEulerR::eVector(): no e vector, call setConstanteVector");
    return siconos::algebra::ConstMapVectorType(eVector_view_->data(), eVector_view_->size());
  }

  /** set a constant \f$ e \f$ for the default linear relation \f$ h(q) = G q + e \f$
   *  (warning: shared memory, neweVector must outlive the relation).
   *
   *  \param neweVector \f$ e \in \mathbb{R}^m \f$
   */
  void setConstanteVector(Eigen::Ref<siconos::algebra::SiconosVector> neweVector);

  /** True if \f$ e \f$ is taken into account */
  bool haseVector() const { return haseVector_; }

  /** \return a read-only view on the contact configuration matrix
   *  \f$ H \in \mathbb{R}^{m \times 6n} \f$, \f$ \dot y = H v \f$ */
  inline auto contactMatrix() const {
    return siconos::algebra::ConstMapType(contactMatrix_->data(), contactMatrix_->rows(),
                                          contactMatrix_->cols());
  }

  /** set a user-defined function computing \f$ \dot G(q, \dot q) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...), called as
   *  fct(q, qdot, constraintJacobianDot)
   */
  void setComputeConstraintJacobianDotFunction(
      const siconos::modeling::func_prototypes::FunctionBVBV_M& fct);

  /** Update \f$ \dot G(q, \dot q) \f$ with the user-defined function, if any.
   *
   *  \param q coordinates \f$ q \f$ of the bodies of the interaction
   *  \param qdot time derivatives \f$ \dot q \f$ of the coordinates
   */
  virtual void computeConstraintJacobianDot(const siconos::algebra::BlockVector& q,
                                            const siconos::algebra::BlockVector& qdot);

  /** \return a read-only view on \f$ \dot H v = (\dot G T + G \dot T)\, v \f$, the terms of
   *  \f$ \ddot y \f$ that do not depend on \f$ \dot v \f$ (see computeContactMatrixDotV)
   */
  inline auto contactMatrixDotV() const {
    return siconos::algebra::ConstMapVectorType(contactMatrixDotV_.data(),
                                                contactMatrixDotV_.size());
  }

  /** Compute \f$ \dot H v = \dot G\, \dot q + G\, \dot T v \f$. Requires a function for
   *  \f$ \dot G \f$ (setComputeConstraintJacobianDotFunction) and up-to-date \f$ \dot q \f$.
   *
   *  \param time  current time
   *  \param inter interaction that owns the relation
   *  \param ds1 first body of the interaction (a NewtonEulerDS)
   *  \param ds2 second body (a NewtonEulerDS); ignored if the interaction has one body,
   *  in which case ds1 may be passed twice
   */
  void computeContactMatrixDotV(double time, Interaction& inter,
                                std::shared_ptr<DynamicalSystem> ds1,
                                std::shared_ptr<DynamicalSystem> ds2);

  /** \return a read-only view on the generalized forces \f$ H^\top \lambda_k \f$ applied
   *  to the bodies by this relation (updated by computeInput, levels 1 and 2)
   */
  inline auto contactForce() const {
    return siconos::algebra::ConstMapVectorType(contactForce_.data(), contactForce_.size());
  }

  /** compute the output \f$ y = h(q) \f$. The default implementation is the linear
   *  relation \f$ y = G q + e \f$, using the current constraintJacobian_view_.
   *
   *  \param[in] q1 coordinates of the first body
   *  \param[in] q2 coordinates of the second body (std::nullopt if there is only one body)
   *  \param[in,out] y the resulting vector
   */
  virtual void computeh(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** compute the contact configuration matrix \f$ H(q) \f$ (contactMatrix_). The default
   * implementation forms the product \f$ G\, T(q) \f$ body by body. Override this function to
   * compute
   *  \f$ H \f$ directly.
   *
   *  \param inter interaction that owns the relation
   *  \param q0  coordinates \f$ q \f$ of the bodies of the interaction
   */
  virtual void computeContactMatrix(const Interaction& inter,
                                    const siconos::algebra::BlockVector& q0);

  /** Update \f$ G \f$ (unless constant, computeConstraintJacobian) and then \f$ H \f$
   * (computeContactMatrix), at the current \f$ q \f$.
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJach(double time, Interaction& inter) override;

  /** compute \f$ y = h(q) \f$ (derivativeNumber = 0) or \f$ \dot y = H v \f$
   *  (derivativeNumber = 1), with \f$ H \f$ as last computed by computeJach.
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param derivativeNumber 0 or 1 (default = 0)
   */
  virtual void computeOutput(
      double time, Interaction& inter,
      siconos::algebra::blocks::size_type derivativeNumber = 0) override;

  /** compute the input \f$ p_k = H^\top \lambda_k \f$ (level k = 1 or 2, also stored in
   *  contactForce) or \f$ p_0 = G^\top \lambda_0 \f$ (level 0, projection onto the
   *  constraints).
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param level 0, 1 or 2 (default = 0)
   */
  virtual void computeInput(double time, Interaction& inter,
                            siconos::algebra::blocks::size_type level = 0) override;

  void display() const override {}

  /** @brief allocation of memory space for relations in the graph
   *      Warning: internal use only (called from Topology)
   *  @param dslink a container of vectors (pointers), from the parent interaction
   *  @param ds1 first ds concerned by the relation
   *  @param ds2 second ds concerned by the relation
   */
  virtual void allocate_read_dynamical_systems_var_vectors(
      std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
      modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const override;

  virtual void accept(relations::Visitor& tourist) const override { tourist.visit(*this); }

  /** \name Deprecated names, kept for one release */
  ///@{
  [[deprecated("use constraintJacobian()")]] auto H_NE() const { return constraintJacobian(); }
  [[deprecated("use setConstantConstraintJacobian")]] void setConstantH_NE(
      Eigen::Ref<siconos::algebra::SiconosMatrix> newValue) {
    setConstantConstraintJacobian(newValue);
  }
  [[deprecated("use contactMatrix()")]] auto H_NE_prod_T() const { return contactMatrix(); }
  [[deprecated("use setComputeConstraintJacobianDotFunction")]] void
  setComputeH_NE_dotFunction(const siconos::modeling::func_prototypes::FunctionBVBV_M& fct) {
    setComputeConstraintJacobianDotFunction(fct);
  }
  [[deprecated("use contactMatrixDotV()")]] auto secondOrderTimeDerivativeTerms() const {
    return contactMatrixDotV();
  }
  [[deprecated("use computeContactMatrixDotV")]] void computeSecondOrderTimeDerivativeTerms(
      double time, Interaction& inter, std::shared_ptr<DynamicalSystem> ds1,
      std::shared_ptr<DynamicalSystem> ds2) {
    computeContactMatrixDotV(time, inter, ds1, ds2);
  }
  /** \deprecated former name of computeContactMatrix. Still called for classes that override
   * it. */
  [[deprecated("override/call computeContactMatrix")]] virtual void computeH_NE_prod_T(
      const Interaction& inter, const siconos::algebra::BlockVector& q0);
  /** \deprecated former name of computeConstraintJacobianDot. Still called for classes that
   * override it. */
  [[deprecated("override/call computeConstraintJacobianDot")]] virtual void computeH_NE_dot(
      const siconos::algebra::BlockVector& q, const siconos::algebra::BlockVector& qdot);
  ///@}
};
}  // namespace siconos::modeling
#endif  // NEWTONEULERR_HPP
