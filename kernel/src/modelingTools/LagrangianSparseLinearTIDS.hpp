/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
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

/*! \file LagrangianSparseLinearTIDS.hpp */
#ifndef LAGRANGIANSPARSETIDS_H
#define LAGRANGIANSPARSETIDS_H

#include <siconos/kernel/SiconosMatrix.hpp>

#include "LagrangianSparseDS.hpp"
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
class LagrangianSparseLinearTIDS : public LagrangianSparseDS {
 protected:
  ACCEPT_SERIALIZATION(LagrangianSparseLinearTIDS);

  /** stiffness matrix */
  std::shared_ptr<siconos::algebra::SiconosSparseMatrix> stiffnessMatrix_{nullptr};

  /** damping matrix */
  std::shared_ptr<siconos::algebra::SiconosSparseMatrix> dampingMatrix_{nullptr};

  /** default constructor */
  LagrangianSparseLinearTIDS() = default;  // Used in FiniteElementLinearTIDS

  /** constructor from initial state - Used for RigidBodies
   *
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   */
  LagrangianSparseLinearTIDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
                             Eigen::Ref<siconos::algebra::SiconosVector> v0)
      : LagrangianSparseDS{q0, v0} {}

 public:
  /** constructor from initial state and mass matrix only. Leads to \f$ M\dot v
   *  = F_{ext}(t) + p \f$ .
   * warning: M will be copied into mass attribute.
   *
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param M mass matrix
   */
  LagrangianSparseLinearTIDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
                             Eigen::Ref<siconos::algebra::SiconosVector> v0,
                             const siconos::algebra::SiconosSparseMatrix &M);

  /** destructor */
  ~LagrangianSparseLinearTIDS() noexcept = default;

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  \param t time of initialization
   */
  void initRhs(double t) override;

  /** set the stiffness matrix. Warning: shared memory with input
   *
   *  \param K new stiffness matrix
   */
  void setStiffnessMatrix(siconos::algebra::SiconosSparseMatrix &K);

  /** Set a constant stiffness matrix for the system
   *  The input matrix is copied into the internal stiffness.
   *
   *  \param newValue stiffness matrix
   *
   */
  void setStiffnessMatrixWithCopy(const siconos::algebra::SiconosSparseMatrix &newValue);

  /** set the damping matrix. Warning: shared memory with input
   *
   *  \param C new damping matrix
   */
  void setDampingMatrix(siconos::algebra::SiconosSparseMatrix &C);

  /** Set a constant damping matrix for the system
   *  The input matrix is copied into the internal damping.
   *
   *  \param newValue damping matrix
   *
   */
  void setDampingMatrixWithCopy(const siconos::algebra::SiconosSparseMatrix &newValue);

  /** \return a read-only ref onto the stiffness matrix */
  inline const siconos::algebra::SiconosSparseMatrix &stiffnessMatrix() const {
    return *stiffnessMatrix_;
  }

  /** True if stiffness matrix is defined */
  bool hasStiffnessMatrix() const { return stiffnessMatrix_ != nullptr; }

  /** True if stiffness matrix is defined */
  bool hasDampingMatrix() const { return dampingMatrix_ != nullptr; }

  /** \return a read-only ref onto the damping matrix */
  inline const siconos::algebra::SiconosSparseMatrix &dampingMatrix() const {
    return *dampingMatrix_;
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
