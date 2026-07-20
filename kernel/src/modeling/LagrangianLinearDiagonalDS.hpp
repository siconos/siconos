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

#include <variant>

#include "LagrangianDS.hpp"
#include "StorageTools.hpp"

namespace siconos::modeling {
/**

    Lagrangian Linear Systems with time invariant and diagonal coefficients
    
    \f$ M\dot v + Cv + Kq = F_{ext}(t,z) + p \f$

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

  /** stiffness */
  siconos::algebra::DenseVectorStorage stiffnessMatrix_storage_{std::monostate{}};

  /** damping  */
  siconos::algebra::DenseVectorStorage dampingMatrix_storage_{std::monostate{}};

  /** mass  */
  siconos::algebra::DenseVectorStorage massMatrix_storage_{std::monostate{}};

  /** default constructor */
  LagrangianLinearDiagonalDS() = delete;
  template <typename F>
  decltype(auto) use_stiffnessMatrix(F&& f) {
    return siconos::algebra::visitStorage(stiffnessMatrix_storage_, std::forward<F>(f),
                                          "stiffnessMatrix_storage_");
  }
  template <typename F>
  decltype(auto) use_stiffnessMatrix(F&& f) const {
    return siconos::algebra::visitStorage(stiffnessMatrix_storage_, std::forward<F>(f),
                                          "stiffnessMatrix_storage_");
  }

  template <typename F>
  decltype(auto) use_dampingMatrix(F&& f) {
    return siconos::algebra::visitStorage(dampingMatrix_storage_, std::forward<F>(f),
                                          "dampingMatrix_storage_");
  }
  template <typename F>
  decltype(auto) use_dampingMatrix(F&& f) const {
    return siconos::algebra::visitStorage(dampingMatrix_storage_, std::forward<F>(f),
                                          "dampingMatrix_storage_");
  }

  template <typename F>
  decltype(auto) use_massMatrix(F&& f) {
    return siconos::algebra::visitStorage(massMatrix_storage_, std::forward<F>(f),
                                          "massMatrix_storage_");
  }

  template <typename F>
  decltype(auto) use_massMatrix(F&& f) const {
    return siconos::algebra::visitStorage(massMatrix_storage_, std::forward<F>(f),
                                          "massMatrix_storage_");
  }

 public:
  /** constructor from initial state and stiffness matrix only.
   *  initial state, velocity and mass attributes will be initialised (copied)
   *  from the input vectors/matrices
   *  @param[in] q0 initial coordinates
   *  @param[in] v0 initial velocity
   *  @param[in] stiffness diagonal of the stiffness matrix
   *  @param[in] tag pass siconos::algebra::copy_t to select this overload
   * (rather than alias version)
   */
  LagrangianLinearDiagonalDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
                             Eigen::Ref<siconos::algebra::SiconosVector> v0,
                             Eigen::Ref<siconos::algebra::SiconosVector> stiffness_diag,
                             siconos::algebra::AliasTag tag);

  /** constructor from initial state and stiffness matrix only.
   *  initial state, velocity and mass attributes will be initialised (copied)
   *  from the input vectors/matrices
   *  @param[in] q0 initial coordinates
   *  @param[in] v0 initial velocity
   *  @param[in] stiffness diagonal of the stiffness matrix
   *  @param[in] tag pass siconos::algebra::copy_t to select this overload
   * (rather than alias version)
   */
  LagrangianLinearDiagonalDS(const siconos::algebra::SiconosVector& q0,
                             const siconos::algebra::SiconosVector& v0,
                             const siconos::algebra::SiconosVector& stiffness_diag,
                             siconos::algebra::CopyTag tag);

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
  Eigen::Ref<const siconos::algebra::SiconosVector> stiffnessMatrix() const {
    return use_stiffnessMatrix(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector>(v); });
  }

  /** @brief set a constant stiffness
   *
   * Warning : deep copy of the provided vector into internal attribute
   *
   * @param newValue external vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
  version)
   *
   */
  void setStiffnessMatrix(const siconos::algebra::SiconosVector& newValue,
                          siconos::algebra::CopyTag tag);

  /** @brief set a constant stiffness
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue external vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setStiffnessMatrix(Eigen::Ref<siconos::algebra::SiconosVector> newValue,
                          siconos::algebra::AliasTag tag);

  /** \return a read-only view on the stiffness matrix */
  Eigen::Ref<const siconos::algebra::SiconosVector> dampingMatrix() const {
    return use_dampingMatrix(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector>(v); });
  }

  /** @brief set a constant damping
   *
   * Warning : deep copy of the provided vector into internal attribute
   *
   * @param newValue external vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
  version)
   *
   */
  void setDampingMatrix(const siconos::algebra::SiconosVector& newValue,
                        siconos::algebra::CopyTag tag);

  /** @brief set a constant damping
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue external vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setDampingMatrix(Eigen::Ref<siconos::algebra::SiconosVector> newValue,
                        siconos::algebra::AliasTag tag);

  /** \return a read-only view on the mass */
  Eigen::Ref<const siconos::algebra::SiconosVector> massMatrix() const {
    return use_massMatrix(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector>(v); });
  }

  /** @brief set a constant mass
   *
   * Warning : deep copy of the provided vector into internal attribute
   *
   * @param newValue external vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
  version)
   *
   */
  void setMassMatrix(const siconos::algebra::SiconosVector& newValue,
                     siconos::algebra::CopyTag tag);

  /** @brief set a constant mass
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue external vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setMassMatrix(Eigen::Ref<siconos::algebra::SiconosVector> newValue,
                     siconos::algebra::AliasTag tag);

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
  void computeTotalForces(const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                          const Eigen::Ref<const siconos::algebra::SiconosVector>& q,
                          double time) override;

  /** Compute  \f$ \nabla_qF_{total}(v,q,t) \f$. Nothing done for this kind of system.
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  void computeJacobianTotalForcesOver_q(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time) override {
    THROW_EXCEPTION("diagonal DS, missing implementation ...");
  };

  /** Compute  \f$ \nabla_{\dot q}F_{total}(v,q,t) \f$. Nothing done for this kind of system.
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  void computeJacobianTotalForcesOver_velocity(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time) override {
    THROW_EXCEPTION("diagonal DS, missing implementation ...");
  };

  /** True if stiffness matrix is defined */
  bool hasStiffnessMatrix() const {
    return !std::holds_alternative<std::monostate>(stiffnessMatrix_storage_);
  }

  /** True if stiffness matrix is defined */
  bool hasDampingMatrix() const {
    return !std::holds_alternative<std::monostate>(dampingMatrix_storage_);
  }

  /** True if mass matrix is defined */
  bool hasMassMatrix() const {
    return !std::holds_alternative<std::monostate>(massMatrix_storage_);
  }

  /**\return true if the Dynamical system is linear. */
  bool isLinear() const override { return true; }

  /** print the data of the dynamical system on the standard output
   */
  void display(bool brief = true) const override;

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // LAGRANGIANLINEARDIAGONALDS_H
