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
/*! \file NewtonEulerR.hpp
 *  \brief Generic relation for Newton-Euler dynamical systems.
 */
#ifndef NEWTONEULERR_HPP
#define NEWTONEULERR_HPP

#include <cassert>
#include <optional>

#include "Relation.hpp"

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
    \dot y &= J_h(q)\, \dot q = J_h(q)\, T(q)\, v = H(q)\, v, \\
    \ddot y &= H(q)\, \dot v + \bigl( \dot J_h(q)\, T(q) + J_h(q)\, \dot T(q) \bigr)\, v,
  \end{aligned}
  \f]
  with the two Jacobians
  \f[
    J_h(q) = \nabla_q^\top h(q) \in \mathbb{R}^{m \times 7n}
    \quad (\text{stored as } H_{NE}), \qquad
    H(q) = J_h(q)\, T(q) \in \mathbb{R}^{m \times 6n}
    \quad (\text{stored as } H_{NE}T).
  \f]
  The generalized forces (impulses) applied to the bodies are, for the input level
  \f$ k \f$,
  \f[
    p_k = H(q)^\top \lambda_k \quad (k = 1, 2), \qquad
    p_0 = J_h(q)^\top \lambda_0 \quad \text{(projection onto the constraints)}.
  \f]
  Only scleronomic relations are handled: \f$ h \f$ does not depend explicitly on time.

  <b>Default implementation.</b> The relation is linear,
  \f$ h(q) = J_h\, q + e \f$, with a constant \f$ J_h \f$ (setConstantH_NE) and an
  optional constant \f$ e \f$ (setConstanteVector). \f$ H = J_h T(q) \f$ is recomputed
  from \f$ T(q) \f$ in computeH_NE_prod_T.

  <b>Derived classes</b> override:
   - computeh() to compute \f$ y = h(q) \f$;
   - computeH_NE_() to update \f$ J_h(q) \f$ (only used when the matrix is not constant);
   - computeH_NE_prod_T() when \f$ H(q) \f$ is computed directly, without forming
     \f$ J_h \f$ (e.g. contact relations written in a local frame, see NewtonEuler1DR).

  \f$ \dot J_h \f$, needed for \f$ \ddot y \f$, can be given by a user-defined function
  (setComputeH_NE_dotFunction).

  <b>Storage.</b>
   - \f$ J_h \f$: a shared view on external memory after setConstantH_NE, internal
     storage otherwise (allocated in initialize, or by the derived class).
   - \f$ H \f$ (H_NE_prod_T_) and the buffer for \f$ T \f$: always allocated in initialize.
   - \f$ e \f$: a shared view on external memory, only after setConstanteVector.
   - \f$ \dot J_h \f$: allocated in initialize if a function has been set, otherwise at the
     first call of computeH_NE_dot or computeSecondOrderTimeDerivativeTerms.

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

  /** View onto \f$ J_h = \nabla_q^\top h \in \mathbb{R}^{m \times 7n} \f$ (either
   *  H_NE_internal_storage_ or external memory, see setConstantH_NE) */
  std::shared_ptr<siconos::algebra::MapType> H_NE_view_{nullptr};

  /** internal (optional) storage used for \f$ J_h \f$ */
  siconos::algebra::SiconosMatrix H_NE_internal_storage_{0, 0};

  /** True if \f$ J_h \f$ is a constant matrix (set with setConstantH_NE) */
  bool hasConstantH_NE_{false};

  /** \f$ \dot J_h = \frac{d}{dt} J_h(q) \in \mathbb{R}^{m \times 7n} \f$, used to
   *  compute the second-order time derivative of \f$ y \f$.
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> H_NE_dot_{nullptr};

  /** user-defined function computing \f$ \dot J_h(q, \dot q) \f$ */
  siconos::modeling::func_prototypes::FunctionBVBV_M computeH_NE_dot_{nullptr};

  /** View onto the constant vector \f$ e \f$ of the default linear relation */
  std::shared_ptr<siconos::algebra::MapVectorType> eVector_view_{nullptr};

  /** True if \f$ e \f$ is taken into account */
  bool haseVector_{false};

  /** Generalized forces (or impulses) applied to the bodies by the relation,
   *  \f$ H^\top \lambda_k \in \mathbb{R}^{6n} \f$, updated by computeInput for
   *  \f$ k = 1, 2 \f$. Useful for the end user. */
  siconos::algebra::SiconosVector contactForce_;

  /** \f$ H = J_h T \in \mathbb{R}^{m \times 6n} \f$, the Jacobian with respect to the
   *  twist: \f$ \dot y = H v \f$. Computed as the product \f$ J_h T \f$ by default, or
   *  directly from the contact geometry by the contact relations. */
  std::shared_ptr<siconos::algebra::SiconosMatrix> H_NE_prod_T_{nullptr};

  /** buffer for \f$ T(q_i) = \begin{bmatrix} I_3 & 0 \\ 0 & \Psi(p_i) \end{bmatrix} \f$
   *  of one body */
  siconos::algebra::SiconosMatrix76 T_buffer_;

  /** The terms of \f$ \ddot y \f$ that do not depend on \f$ \dot v \f$:
   *  \f$ \dot J_h\, \dot q + J_h\, \dot T v = (\dot J_h T + J_h \dot T)\, v \f$.
   */
  siconos::algebra::SiconosVector secondOrderTimeDerivativeTerms_;

  /** allocate \f$ \dot J_h \f$ (H_NE_dot_), with the size of \f$ J_h \f$, if not done */
  void allocateH_NE_dot_();

  /** update \f$ J_h(q) \f$ (H_NE_view_). Called by computeJach unless the matrix is
   *  constant. The default does nothing.
   *
   *  \param time current time
   *  \param inter the interaction using this relation
   *  \param q0  coordinates \f$ q \f$ of the bodies of the interaction
   */
  virtual void computeH_NE_(double time, Interaction& inter,
                            const siconos::algebra::BlockVector& q0) {}

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

  /** \return a read-only view on \f$ J_h = \nabla_q^\top h \in \mathbb{R}^{m \times 7n}
   *  \f$ */
  inline auto H_NE() const {
    return siconos::algebra::ConstMapType(H_NE_view_->data(), H_NE_view_->rows(),
                                          H_NE_view_->cols());
  }

  /** Set a constant \f$ J_h \f$ for the default linear relation
   *  \f$ h(q) = J_h q + e \f$ (warning: shared memory, newValue must outlive the
   *  relation).
   *
   *  \param newValue \f$ J_h \in \mathbb{R}^{m \times 7n} \f$
   */
  void setConstantH_NE(Eigen::Ref<siconos::algebra::SiconosMatrix> newValue);

  /** \return a read-only view on \f$ e \f$. Requires haseVector(). */
  inline auto eVector() const {
    assert(eVector_view_ && "NewtonEulerR::eVector(): no e vector, call setConstanteVector");
    return siconos::algebra::ConstMapVectorType(eVector_view_->data(), eVector_view_->size());
  }

  /** set a constant \f$ e \f$ for the default linear relation
   *  \f$ h(q) = J_h q + e \f$ (warning: shared memory, neweVector must outlive the
   *  relation).
   *
   *  \param neweVector \f$ e \in \mathbb{R}^m \f$
   */
  void setConstanteVector(Eigen::Ref<siconos::algebra::SiconosVector> neweVector);

  /** True if \f$ e \f$ is taken into account */
  bool haseVector() const { return haseVector_; }

  /** \return a read-only view on \f$ H = J_h T \in \mathbb{R}^{m \times 6n} \f$, so that
   *  \f$ \dot y = H v \f$ */
  inline auto H_NE_prod_T() const {
    return siconos::algebra::ConstMapType(H_NE_prod_T_->data(), H_NE_prod_T_->rows(),
                                          H_NE_prod_T_->cols());
  }

  /** set a user-defined function computing \f$ \dot J_h(q, \dot q) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...), called as
   *  fct(q, qdot, Jdot)
   */
  void setComputeH_NE_dotFunction(
      const siconos::modeling::func_prototypes::FunctionBVBV_M& fct);

  /** Update \f$ \dot J_h(q, \dot q) \f$ with the user-defined function, if any.
   *
   *  \param q coordinates \f$ q \f$ of the bodies of the interaction
   *  \param qdot time derivatives \f$ \dot q \f$ of the coordinates
   */
  virtual void computeH_NE_dot(const siconos::algebra::BlockVector& q,
                               const siconos::algebra::BlockVector& qdot);

  /** \return a read-only view on the terms of \f$ \ddot y \f$ that do not depend on
   *  \f$ \dot v \f$, \f$ (\dot J_h T + J_h \dot T)\, v \f$
   *  (see computeSecondOrderTimeDerivativeTerms)
   */
  inline auto secondOrderTimeDerivativeTerms() const {
    return siconos::algebra::ConstMapVectorType(secondOrderTimeDerivativeTerms_.data(),
                                                secondOrderTimeDerivativeTerms_.size());
  }

  /** Compute the terms of \f$ \ddot y \f$ that do not depend on \f$ \dot v \f$,
   *  \f$ \dot J_h\, \dot q + J_h\, \dot T v \f$. Requires a function for
   *  \f$ \dot J_h \f$ (setComputeH_NE_dotFunction) and up-to-date \f$ \dot q \f$.
   *
   *  \param time  current time
   *  \param inter interaction that owns the relation
   *  \param ds1 first body of the interaction (a NewtonEulerDS)
   *  \param ds2 second body (a NewtonEulerDS); ignored if the interaction has one body,
   *  in which case ds1 may be passed twice
   */
  void computeSecondOrderTimeDerivativeTerms(double time, Interaction& inter,
                                             std::shared_ptr<DynamicalSystem> ds1,
                                             std::shared_ptr<DynamicalSystem> ds2);

  /** \return a read-only view on the generalized forces \f$ H^\top \lambda_k \f$ applied
   *  to the bodies by this relation (updated by computeInput, levels 1 and 2)
   */
  inline auto contactForce() const {
    return siconos::algebra::ConstMapVectorType(contactForce_.data(), contactForce_.size());
  }

  /** compute the output \f$ y = h(q) \f$. The default implementation is the linear
   *  relation \f$ y = J_h q + e \f$, using the current H_NE_view_.
   *
   *  \param[in] q1 coordinates of the first body
   *  \param[in] q2 coordinates of the second body (std::nullopt if there is only one body)
   *  \param[in,out] y the resulting vector
   */
  virtual void computeh(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y);

  /** compute \f$ H(q) \f$ (H_NE_prod_T_). The default implementation forms the product
   *  \f$ J_h\, T(q) \f$ body by body, with \f$ T(q_i) \f$ computed from \f$ q_i \f$.
   *
   *  \param inter interaction that owns the relation
   *  \param q0  coordinates \f$ q \f$ of the bodies of the interaction
   */
  virtual void computeH_NE_prod_T(const Interaction& inter,
                                  const siconos::algebra::BlockVector& q0);

  /** Update \f$ J_h \f$ (unless constant, computeH_NE_) and then \f$ H \f$
   *  (computeH_NE_prod_T), at the current \f$ q \f$.
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
   *  contactForce) or \f$ p_0 = J_h^\top \lambda_0 \f$ (level 0, projection onto the
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
};
}  // namespace siconos::modeling
#endif  // NEWTONEULERR_HPP
