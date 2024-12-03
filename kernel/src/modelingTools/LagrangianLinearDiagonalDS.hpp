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

/*! \file LagrangianLinearDiagonalDS.hpp */
#ifndef LAGRANGIANLINEARDIAGONALDS_H
#define LAGRANGIANLINEARDIAGONALDS_H

#include "LagrangianDS.hpp"

namespace siconos::modeling {
/**

    Lagrangian Linear Systems with time invariant and diagonal coefficients -  \f$ M\dot v + Cv
   + Kq = F_{ext}(t,z) + p \f$


    where

    -  \f$ q \in R^{ndof} \f$ is the set of the generalized coordinates,
    - \f$ \dot q = v \in R^{ndof} \f$  the velocity, i. e. the time derivative of the
   generalized coordinates.
    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the
   generalized coordinates.
    - \f$ p  \in R^{ndof} \f$  the forces due to the nonsmooth interaction. In the particular
   case of a nonsmooth evolution, the variable p contains the impulse and not the force.
    -  \f$ M \in  R^{ndof \times ndof} \f$ is the mass matrix (access : mass() method).
    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix (access : stiffness()
   method).
    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix (access : damping() method).
    -  \f$ z \in R^{zSize} \f$ is a vector of arbitrary algebraic variables, some sort of
   discret state.

    Remind that the specificity of this class is that all matrices are diagonal (and hence only
   diagonal coefficients are saved in memory).

    For details about dynamical systems in Siconos, please read user's guide.
*/
class LagrangianLinearDiagonalDS : public LagrangianDS {
 protected:
  ACCEPT_SERIALIZATION(LagrangianLinearDiagonalDS);

  /** stiffness matrix */
  std::shared_ptr<siconos::algebra::MapVectorType> stiffnessMatrix_view_{nullptr};

  /** damping matrix */
  std::shared_ptr<siconos::algebra::MapVectorType> dampingMatrix_view_{nullptr};

  /** mass matrix */
  std::shared_ptr<siconos::algebra::MapVectorType> massMatrix_view_{nullptr};

  /** default constructor */
  LagrangianLinearDiagonalDS() = delete;

 public:
  /** constructor from initial state and all operators.
   *
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param stiffness diagonal of the stiffness matrix
   */
  LagrangianLinearDiagonalDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
                             Eigen::Ref<siconos::algebra::SiconosVector> v0,
                             Eigen::Ref<siconos::algebra::SiconosVector> stiffness_diag);

  /* destructor */
  ~LagrangianLinearDiagonalDS() noexcept = default;

  /** set the damping matrix. Warning: shared memory with input
   *
   *  \param C diagonal of the damping matrix
   */
  void setDampingMatrix(Eigen::Ref<siconos::algebra::SiconosVector> C);

  /** set the mass matrix. Warning: shared memory with input
   *
   *  \param M diagonal of the mass matrix
   */
  void setMassMatrix(Eigen::Ref<siconos::algebra::SiconosVector> M);

  /** \return a read-only view on the stiffness matrix */
  inline const auto stiffnessMatrix() const {
    return siconos::algebra::ConstMapVectorType(stiffnessMatrix_view_->data(),
                                                stiffnessMatrix_view_->size());
  }

  /** \return a read-only view on the damping matrix */
  inline const auto dampingMatrix() const {
    return siconos::algebra::ConstMapVectorType(dampingMatrix_view_->data(),
                                                dampingMatrix_view_->size());
  }

  /** \return a read-only view on the damping matrix */
  inline const auto massMatrix() const {
    return siconos::algebra::ConstMapVectorType(massMatrix_view_->data(),
                                                massMatrix_view_->size());
  }

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  \param t time of initialization
   */
  void initRhs(double t) override;

  /** Compute  \f$ F_{total}(v,q,t) = -Kq - Cv + f_{ext}(t)\f$
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  void computeTotalForces(const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
                          const Eigen::Ref<const siconos::algebra::SiconosVector> &q,
                          double time) override;

  /** Compute  \f$ \nabla_qF_{total}(v,q,t) \f$. Nothing done for this kind of system.
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  void computeJacobianTotalForcesOver_q(
      const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
      const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time) override {
    THROW_EXCEPTION("diagonal DS, missing implementation ...");
  };

  /** Compute  \f$ \nabla_{\dot q}F_{total}(v,q,t) \f$. Nothing done for this kind of system.
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  void computeJacobianTotalForcesOver_velocity(
      const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
      const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time) override {
    THROW_EXCEPTION("diagonal DS, missing implementation ...");
  };

  /** True if stiffness matrix is defined */
  bool hasStiffnessMatrix() const { return stiffnessMatrix_view_ != nullptr; }

  /** True if stiffness matrix is defined */
  bool hasDampingMatrix() const { return dampingMatrix_view_ != nullptr; }

  /** True if mass matrix is defined */
  bool hasMassMatrix() const { return massMatrix_view_ != nullptr; }

  /**\return true if the Dynamical system is linear. */
  bool isLinear() override { return true; }

  /** print the data of the dynamical system on the standard output
   */
  void display(bool brief = true) const override;

  Type acceptType(types::FindType &ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // LAGRANGIANLINEARDIAGONALDS_H
