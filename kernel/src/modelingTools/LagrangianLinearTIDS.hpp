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

/*! \file LagrangianLinearTIDS.hpp */
#ifndef LAGRANGIANTIDS_H
#define LAGRANGIANTIDS_H

#include "LagrangianDS.hpp"
#include "SiconosPointers.hpp"

namespace siconos::modeling {

/**
    Lagrangian Linear Systems with time invariant coefficients -

    \f$ M\dot v + Cv + Kq = F_{ext}(t) + p \f$

    where
    -  \f$ q \in R^{ndof} \f$ is the set of the generalized coordinates,
    - \f$ v \in R^{ndof} \f$  the velocity, i. e. the time derivative of
    the  generalized coordinates.
    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time
    derivative of the  generalized coordinates.
    - \f$ p \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In
    particular case of Non Smooth evolution, the variable p contains the impulse
    and not the force.
    -  \f$ M \in  R^{ndof \times ndof} \f$ is the Mass matrix (access : mass()
    method).
    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix (access : K()
    method).
    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix (access : C()
    method).

    The equation of motion is also shortly denoted as (as in LagrangianDS):

    \f[

    M(q) \dot v = F_{total}(v, q, t) + p

    \f]

    where
    -  \f$ F_{total}(v, q, t) \in R^{ndof} \f$ collects the total forces
    acting on the system, that is
    \f$ F_{total}(v, q, t) =  F_{ext}(t) -  Cv - Kq \f$.

    This vector is saved and may be accessed using totalForces() method.

    If required (e.g. for Event-Driven like simulation), reformulation as a
    first-order system is also available, and writes:

    - \f$ n= 2 ndof \f$
    - \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right] \f$
    - rhs given by:

    \f[

        rhs(x,t) = \left[\begin{array}{c}
        \dot q  \\
        \ddot q = M^{-1}\left[F_{ext}(t) - C \dot q - K q  + p \right] \\
        \end{array}\right]

    \f]

    Its jacobian is:

    \f[

      \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
      0   & I \\
      -M^{-1}K & -M^{-1}C \\
      \end{array}\right]

    \f]

    with the input due to the non smooth law:

    \f[
      r = \left[\begin{array}{c}0 \\ p \end{array}\right]

    \f]
*/
class LagrangianLinearTIDS : public LagrangianDS {
 protected:
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);

  /** stiffness matrix */
  std::shared_ptr<siconos::algebra::MapType> stiffnessMatrix_view_{nullptr};

  /** damping matrix */
  std::shared_ptr<siconos::algebra::MapType> dampingMatrix_view_{nullptr};

  /** default constructor */
  LagrangianLinearTIDS() = default;  // Used in FiniteElementLinearTIDS

 public:
  /** constructor from initial state and mass matrix only. Leads to \f$ M\dot v
   *  = F_{ext}(t) + p \f$ .
   *
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param M mass matrix
   */
  LagrangianLinearTIDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
                       Eigen::Ref<siconos::algebra::SiconosVector> v0,
                       Eigen::Ref<siconos::algebra::SiconosMatrix> M);
  /** destructor */
  ~LagrangianLinearTIDS() noexcept = default;

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  \param t time of initialization
   */
  void initRhs(double t) override;

  /** set the stiffness matrix. Warning: shared memory with input
   *
   *  \param K new stiffness matrix
   */
  void setStiffnessMatrix(Eigen::Ref<siconos::algebra::SiconosMatrix> K);

  /** set the damping matrix. Warning: shared memory with input
   *
   *  \param C new damping matrix
   */
  void setDampingMatrix(Eigen::Ref<siconos::algebra::SiconosMatrix> C);

  /** \return a read-only view on the stiffness matrix */
  inline const auto stiffnessMatrix() const {
    return siconos::algebra::ConstMapType(stiffnessMatrix_view_->data(),
                                          stiffnessMatrix_view_->rows(),
                                          stiffnessMatrix_view_->cols());
  }

  /** True if stiffness matrix is defined */
  bool hasStiffnessMatrix() const { return stiffnessMatrix_view_ != nullptr; }

  /** True if stiffness matrix is defined */
  bool hasDampingMatrix() const { return dampingMatrix_view_ != nullptr; }

  /** \return a read-only view on the damping matrix */
  inline const auto dampingMatrix() const {
    return siconos::algebra::ConstMapType(
        dampingMatrix_view_->data(), dampingMatrix_view_->rows(), dampingMatrix_view_->cols());
  }

  /** Compute  \f$ F_{total}(v,q,t) = -Kq - Cv + f_{ext}(t)\f$
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  void computeTotalForces(const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
                          const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
                          double time) override;

  /** \return true if the Dynamical system is linear.
   */
  bool isLinear() const override { return true; }

  /** print the data onto the screen
   */
  void display(bool brief = true) const override;

  Type acceptType(types::FindType &ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // LAGRANGIANTIDS_H
