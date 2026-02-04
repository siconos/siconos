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
#include "StorageTools.hpp"

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
  siconos::algebra::DenseStorage stiffnessMatrix_storage_{std::monostate{}};
  template <typename F>
  decltype(auto) useStiffness(F&& f) const {
    return siconos::algebra::visitStorage(stiffnessMatrix_storage_, std::forward<F>(f),
                                          "stiffnessMatrix_storage_");
  }

  /** damping matrix */
  siconos::algebra::DenseStorage dampingMatrix_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) useDamping(F&& f) const {
    return siconos::algebra::visitStorage(dampingMatrix_storage_, std::forward<F>(f),
                                          "dampingMatrix_storage_");
  }

  /** default constructor */
  LagrangianLinearTIDS() = default;  // Used in FiniteElementLinearTIDS

  /** constructor from initial state - Used for RigidBodies
   *
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   */
  LagrangianLinearTIDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
                       Eigen::Ref<siconos::algebra::SiconosVector> v0,
                       siconos::algebra::AliasTag)
      : LagrangianDS{q0, v0, siconos::algebra::alias_t} {}
  LagrangianLinearTIDS(const siconos::algebra::SiconosVector& q0,
                       const siconos::algebra::SiconosVector& v0, siconos::algebra::CopyTag)
      : LagrangianDS{q0, v0, siconos::algebra::copy_t} {}

 public:
  /** constructor from initial state and mass matrix only. Leads to \f$ M\dot v
   *  = F_{ext}(t) + p \f$ .
   *  initial state, velocity and mass attributes will be initialised (copied)
   *  from the input vectors/matrices
   *  @param[in] position initial coordinates
   *  @param[in] velocity initial velocity
   *  @param[in] M mass matrix
   *  @param[in] tag pass siconos::algebra::copy_t to select this overload
   * (rather than alias version)

   */
  LagrangianLinearTIDS(const siconos::algebra::SiconosVector& position,
                       const siconos::algebra::SiconosVector& velocity,
                       const siconos::algebra::SiconosDenseMatrix& M,
                       siconos::algebra::CopyTag tag);

  /** constructor from initial state and mass matrix only. Leads to \f$ M\dot v
   *  = F_{ext}(t) + p \f$ .
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means that for initial position and velocity
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   *  @param[in] position initial coordinates
   *  @param[in] velocity initial velocity
   *  @param[in] M mass matrix
   *  @param tag Pass siconos::algebra::alias_t to select this overload
   * (rather than copy version)
   */
  LagrangianLinearTIDS(Eigen::Ref<siconos::algebra::SiconosVector> position,
                       Eigen::Ref<siconos::algebra::SiconosVector> velocity,
                       Eigen::Ref<siconos::algebra::SiconosDenseMatrix> M,
                       siconos::algebra::AliasTag tag);

  /** destructor */
  ~LagrangianLinearTIDS() noexcept = default;

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  \param t time of initialization
   */
  void initRhs(double t) override;

  /** @brief set a constant stiffness matrix for the system
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   * @param newValue matrix to be copied. Its shape must match dimension() x dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  void setStiffnessMatrix(const siconos::algebra::SiconosDenseMatrix& newValue,
                          siconos::algebra::CopyTag tag);

  /** @brief set a constant stiffness matrix for the system
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue external force vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setStiffnessMatrix(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
                          siconos::algebra::AliasTag tag);

  /** @brief set a constant damping matrix for the system
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   * @param newValue matrix to be copied. Its shape must match dimension() x dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  void setDampingMatrix(const siconos::algebra::SiconosDenseMatrix& newValue,
                        siconos::algebra::CopyTag tag);

  /** @brief set a constant damping matrix for the system
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue external force vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setDampingMatrix(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
                        siconos::algebra::AliasTag tag);

  /** \return a read-only view on the stiffness matrix */
  Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> stiffnessMatrix() const {
    return useStiffness([](auto const& K) {
      return Eigen::Ref<const siconos::algebra::SiconosDenseMatrix>(K);
    });
  }

  /** True if stiffness matrix is defined */
  bool hasStiffnessMatrix() const {
    return !std::holds_alternative<std::monostate>(stiffnessMatrix_storage_);
  }

  /** True if stiffness matrix is defined */
  bool hasDampingMatrix() const {
    return !std::holds_alternative<std::monostate>(dampingMatrix_storage_);
  }

  /** \return a read-only view on the damping matrix */
  Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> dampingMatrix() const {
    return useDamping([](auto const& C) {
      return Eigen::Ref<const siconos::algebra::SiconosDenseMatrix>(C);
    });
  }

  /** Compute  \f$ F_{total}(v,q,t) = -Kq - Cv + f_{ext}(t)\f$
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  void computeTotalForces(const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                          const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
                          double time) override;

  /** \return true if the Dynamical system is linear.
   */
  bool isLinear() const override { return true; }

  /** print the data onto the screen
   */
  void display(bool brief = true) const override;

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // LAGRANGIANTIDS_H
