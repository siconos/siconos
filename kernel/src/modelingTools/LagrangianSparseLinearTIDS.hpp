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

#include "LagrangianSparseDS.hpp"
#include "SiconosMatrix.hpp"
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
class LagrangianSparseLinearTIDS : public LagrangianSparseDS {
 protected:
  ACCEPT_SERIALIZATION(LagrangianSparseLinearTIDS);

  /** stiffness matrix */
  siconos::algebra::SparseStorage stiffnessMatrix_storage{std::monostate{}};

  /** damping matrix */
  siconos::algebra::SparseStorage dampingMatrix_storage{std::monostate{}};

  /** default constructor */
  LagrangianSparseLinearTIDS() = default;  // Used in FiniteElementLinearTIDS

  /**
   * @brief Utility function providing uniform access to the stiffness matrix.
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
  decltype(auto) useStiffness(F&& f) const {
    return siconos::algebra::visitStorage(stiffnessMatrix_storage, std::forward<F>(f),
                                          "stiffnessMatrix_storage");
  }

  /**
   * @brief Utility function providing uniform access to the damping matrix.
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
  decltype(auto) useDamping(F&& f) const {
    return siconos::algebra::visitStorage(dampingMatrix_storage, std::forward<F>(f),
                                          "dampingMatrix_storage");
  }

  // Non const versions. Not needed but kept for the record ...
  // template <typename F>
  // decltype(auto) useStiffness(F&& f) {
  //   if (owned_stiffnessMatrix_) return f(*owned_stiffnessMatrix_);
  //   if (external_stiffnessMatrix_) return f(*external_mass_mat_);
  //   throw std::runtime_error("mass not set");
  // }
  // template <typename F>
  // decltype(auto) useDamping(F&& f) {
  //   if (owned_dampingMatrix_) return f(*owned_dampingMatrix_);
  //   if (external_dampingMatrix_) return f(*external_dampingMatrix_);
  //   throw std::runtime_error("mass not set");
  // }

 public:
  /** constructor from initial state and velocity with memory alias
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means that for initial position and velocity
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   *  @param[in] position initial coordinates
   *  @param[in] velocity initial velocity
   *  @param tag Pass siconos::algebra::alias_t to select this overload
   * (rather than copy version)
   */
  LagrangianSparseLinearTIDS(Eigen::Ref<siconos::algebra::SiconosVector> position,
                             Eigen::Ref<siconos::algebra::SiconosVector> velocity,
                             siconos::algebra::AliasTag)
      : LagrangianSparseDS(position, velocity, siconos::algebra::alias_t) {}

  /** constructor from initial state and velocity with copy
   *
   *  initial state and velocity attributes will be initialised (copied)
   *  from the input vectors
   *
   *  @param[in] position initial coordinates
   *  @param[in] velocity initial velocity
   *  @param tag Pass siconos::algebra::copy_t to select this overload
   * (rather than alias version)
   */
  LagrangianSparseLinearTIDS(const siconos::algebra::SiconosVector& position,
                             const siconos::algebra::SiconosVector& velocity,
                             siconos::algebra::CopyTag)
      : LagrangianSparseDS(position, velocity, siconos::algebra::copy_t) {}

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
  LagrangianSparseLinearTIDS(const siconos::algebra::SiconosVector& position,
                             const siconos::algebra::SiconosVector& velocity,
                             const siconos::algebra::SiconosSparseMatrix& M,
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
  LagrangianSparseLinearTIDS(Eigen::Ref<siconos::algebra::SiconosVector> position,
                             Eigen::Ref<siconos::algebra::SiconosVector> velocity,
                             Eigen::Map<siconos::algebra::SiconosSparseMatrix>& M,
                             siconos::algebra::AliasTag tag);

  /** destructor */
  ~LagrangianSparseLinearTIDS() noexcept = default;

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  @param t time of initialization
   */
  void initRhs(double t) override;

  /** @brief set the stiffness matrix.
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   * @param K new stiffness matrix
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   */
  void setStiffnessMatrix(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& K,
                          siconos::algebra::AliasTag);

  /** @brief Set a constant stiffness matrix for the system
   *
   *  Warning : deep copy of the provided object into internal attribute
   *
   *  @param newValue stiffness matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setStiffnessMatrix(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    if (newValue.rows() != ndof_ || newValue.cols() != ndof_)
      throw std::invalid_argument("Input stiffness matrix has wrong dimensions");

    stiffnessMatrix_storage =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
  }

  /** @brief set the damping matrix.
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param V new damping matrix
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   */
  void setDampingMatrix(Eigen::Map<siconos::algebra::SiconosSparseMatrix>& C,
                        siconos::algebra::AliasTag);

  /** @brief Set a constant damping matrix for the system
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   *  @param newValue damping matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  template <typename T>
  void setDampingMatrix(T&& newValue, siconos::algebra::CopyTag) {
    static_assert(std::is_same_v<std::decay_t<T>, siconos::algebra::SiconosSparseMatrix>,
                  "Type must be SiconosSparseMatrix");

    if (newValue.rows() != ndof_ || newValue.cols() != ndof_)
      throw std::invalid_argument("Input damping matrix has wrong dimensions");

    dampingMatrix_storage =
        std::make_unique<siconos::algebra::SiconosSparseMatrix>(std::forward<T>(newValue));
  }

  /** \return a read-only ref onto the stiffness matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> stiffnessMatrix() const {
    return useStiffness([](Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> M) {
      return Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>(M);
    });
  }

  /*  \return a reference to the stiffness matrix

    pybind11 use only!
  */
  const siconos::algebra::SiconosSparseMatrix& stiffnessMatrix_py() const {
    return siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
        stiffnessMatrix_storage,
        [](const auto& M) -> const siconos::algebra::SiconosSparseMatrix& { return M; },
        "stiffnessMatrix_storage");
  }

  /** True if stiffness matrix is defined */
  bool hasStiffnessMatrix() const {
    return !std::holds_alternative<std::monostate>(stiffnessMatrix_storage);
  }

  /** True if damping matrix is defined */
  bool hasDampingMatrix() const {
    return !std::holds_alternative<std::monostate>(dampingMatrix_storage);
  }

  /** \return a read-only ref onto the stiffness matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> dampingMatrix() const {
    return useDamping([](Eigen::Ref<const siconos::algebra::SiconosSparseMatrix> M) {
      return Eigen::Ref<const siconos::algebra::SiconosSparseMatrix>(M);
    });
  }

  /*  \return a reference to the damping matrix

    pybind11 use only!
  */
  const siconos::algebra::SiconosSparseMatrix& dampingMatrix_py() const {
    return siconos::algebra::visitStorage<siconos::algebra::AccessMode::OwnedOnly>(
        dampingMatrix_storage,
        [](const auto& M) -> const siconos::algebra::SiconosSparseMatrix& { return M; },
        "dampingMatrix_storage");
  }

  /** Compute  \f$ F_{total}(v,q,t) = -Kq - Cv + f_{ext}(t)\f$
   *
   *  @param velocity vector
   *  @param q state
   *  @param time the current time
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
