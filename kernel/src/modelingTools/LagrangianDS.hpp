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

/*! \file LagrangianDS.hpp
  LagrangianDS class - Second Order Non Linear Dynamical Systems.
*/

#ifndef LagrangianDS_H
#define LagrangianDS_H

#include "FunctionTypes.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"

namespace siconos::modeling {

/**
   Lagrangian non linear dynamical systems -  \f$ M(q) \dot v = F_{total}(v, q, t) + p \f$

   This class defines and computes a generic ndof-dimensional Lagrangian Non Linear Dynamical
   System of the form :

   \f[
     M(q) \dot v + F_{gyr}(v, q) + F_{int}(v , q , t) = F_{ext}(t) + p \\
     \dot q = v
   \f]

   where

   -  \f$ q \in R^{ndof} \f$ is the set of the generalized coordinates,

   - \f$ \dot q =v \in R^{ndof} \f$ the velocity, i. e. the time derivative of the generalized
   coordinates (Lagrangian systems).

   - \f$ \ddot q =\\dot v \in R^{ndof} \f$ the acceleration, i. e. the second time derivative
   of the generalized coordinates.

   - \f$ p \in R^{ndof} \f$ the reaction forces due to the Non Smooth Interaction.

   - \f$ M(q) \in R^{ndof \times ndof} \f$ is the inertia term (access : mass() method).

   - \f$ F_{gyr}(\dot q, q) \in R^{ndof} \f$ is the non linear inertia term
   method).

   - \f$ F_{int}(\dot q , q , t) \in R^{ndof} \f$ are the internal forces
   method).

   - \f$ F_{ext}(t) \in R^{ndof} \f$ are the external forces

   The equation of motion is also shortly denoted as  \f$ M(q) \dot v = F_{total}(v,q, t) + p
   \f$

   where  \f$ F_{total}(v, q, t) \in R^{ndof} \f$ collects the total forces acting on
   the system, that is \f$ F_{total}(v, q, t) =  F_{ext}(t) -  F_{gyr}(v, q) - F_{int}(v, q, t)
   \f$

   This vector is saved and may be accessed using totalForces() method.

   q[i] is the derivative number i of q.
   Thus: q[0]=\f$ q \f$, global coordinates, q[1]= \f$  \dot q \f$ , velocity,
   q[2]= \f$  \ddot q \f$, acceleration.

   The following operators (and their jacobians) can be plugged, in the usual
   way (see User Guide, 'User-defined plugins')

   -  \f$ M(q) \f$  (computeMass(q))
   -  \f$ F_{gyr}(v, q) \f$  (computeFgyr())
   -  \f$ F_{int}(v , q , t) \f$  (computeFint())
   -  \f$ F_{ext}(t) \f$  (computeFext())

   If required (e.g. for Event-Driven like simulation), formulation as a
   first-order system is also available, and writes:

   - \f$ n= 2 ndof \f$
   - \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right] \f$
   - rhs given by:

   \f[

      \dot x = \left[\begin{array}{c}
      \dot q\\
      \ddot q = M^{-1}(q)\left[F_{total}(v, q , t) + p \right]\\
      \end{array}\right]

   \f]   \endrst

   - jacobian of the rhs, with respect to x

   \f[
       \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
       0  & I \\
       \nabla_{q}(M^{-1}(q)F_{total}(v, q , t)) &  \nabla_{\dot q}(M^{-1}(q)F_{total}(v, q ,t))
   \\ \end{array}\right] \f]

   with the input due to the non smooth law:

   \f[

      \left[\begin{array}{c}
      0 \\
      p \end{array}\right]
   \f]

   In that case, use the following methods:
   - initRhs() to allocate/initialize memory for these new operators,
   - rhs() to get the rhs vector
   - computeRhs(), computeJacobianRhsOver_x() ..., to update the content of rhs, its
   jacobians ...

*/
class LagrangianDS : public SecondOrderDS {
 protected:
  ACCEPT_SERIALIZATION(LagrangianDS);

  /** state of the system. See details on top of page. */
  std::vector<std::shared_ptr<siconos::algebra::SiconosVector>> state_q_ = {nullptr, nullptr,
                                                                            nullptr};

  /** Initial position of the system */
  siconos::algebra::DenseVectorStorage q0_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_q0(F&& f) const {
    return siconos::algebra::visitStorage(q0_storage_, std::forward<F>(f), "q0_storage_");
  }

  /** initial velocity of the system */
  siconos::algebra::DenseVectorStorage velocity0_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_velocity0(F&& f) const {
    return siconos::algebra::visitStorage(velocity0_storage_, std::forward<F>(f),
                                          "v0_storage_");
  }

  /** memory of previous coordinates of the system */
  siconos::algebra::SiconosMemory qMemory_;

  /** memory of previous velocities of the system */
  siconos::algebra::SiconosMemory velocityMemory_;

  /** mass matrix of the system (as a view onto memory)*/
  siconos::algebra::DenseStorage mass_storage_{std::monostate{}};

  /** function wrapper used to compute mass */
  siconos::modeling::func_prototypes::FunctionV_M computemass_{nullptr};

  /** true if mass is required/set and constant */
  bool hasConstantMass_{false};

  /** True if mass is required */
  bool hasMass_{false};

  /** inverse or factorization of the mass of the system */
  std::shared_ptr<siconos::algebra::SiconosDenseLUMatrix> LUMass_{nullptr};

  /** internal forces (\f$ F_{int}(v , q , t) \f$) applied to the system */
  std::unique_ptr<siconos::algebra::SiconosVector> fint_{nullptr};

  /** function wrapper used to compute internal forces */
  siconos::modeling::func_prototypes::FunctionVVS_V computefint_{nullptr};

  /** True if internal forces are  taken into account */
  bool hasFint_{false};

  /** jacobianFintOver_q matrix */
  siconos::algebra::DenseStorage jacobianFintOver_q_storage_{std::monostate{}};

  /** function wrapper used to compute jacobianFintOver_q */
  siconos::modeling::func_prototypes::FunctionVVS_M computejacobianFintOver_q_{nullptr};

  /** true if jacobianFintOver_q is not set or set and constant */
  bool hasConstantJacobianFintOver_q_{true};

  /** True if jacobianFintOver_q is required */
  bool hasJacobianFintOver_q_{false};

  /** jacobianFintOver_velocity matrix */
  siconos::algebra::DenseStorage jacobianFintOver_velocity_storage_{std::monostate{}};

  /** function wrapper used to compute jacobianFintOver_velocity */
  siconos::modeling::func_prototypes::FunctionVVS_M computejacobianFintOver_velocity_{nullptr};

  /** true if jacobianFintOver_velocity is not set or required/set and constant */
  bool hasConstantJacobianFintOver_velocity_{true};

  /** True if jacobianFintOver_velocity is required */
  bool hasJacobianFintOver_velocity_{false};

  /** external forces applied to the system */
  siconos::algebra::DenseVectorStorage fext_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_fext(F&& f) const {
    return siconos::algebra::visitStorage(fext_storage_, std::forward<F>(f), "fext_storage_");
  }

  /** function wrapper used to compute external forces */
  siconos::modeling::func_prototypes::FunctionS_V computefext_{nullptr};

  /** True if external forces are taken into account and constant */
  bool hasConstantFext_{false};

  /** True if external forces are taken into account */
  bool hasFext_{false};

  /** non-linear inertia term (\f$ F_{gyr}(v, q) \f$ ) applied to the system */
  std::unique_ptr<siconos::algebra::SiconosVector> fgyr_{nullptr};

  /** function wrapper used to compute non-linear inertia term */
  siconos::modeling::func_prototypes::FunctionVV_V computefgyr_{nullptr};

  /** True if non-linear inertia terms are  taken into account */
  bool hasFgyr_{false};

  /** jacobianFgyrOver_q matrix */
  siconos::algebra::DenseStorage jacobianFgyrOver_q_storage_{std::monostate{}};

  /** function wrapper used to compute jacobianFgyrOver_q */
  siconos::modeling::func_prototypes::FunctionVV_M computejacobianFgyrOver_q_{nullptr};

  /** true if jacobianFgyrOver_q is not set or required/set and constant */
  bool hasConstantJacobianFgyrOver_q_{true};

  /** True if jacobianFgyrOver_q is required */
  bool hasJacobianFgyrOver_q_{false};

  /** jacobianFgyrOver_velocity matrix */
  siconos::algebra::DenseStorage jacobianFgyrOver_velocity_storage_{std::monostate{}};

  /** function wrapper used to compute jacobianFgyrOver_velocity */
  siconos::modeling::func_prototypes::FunctionVV_M computejacobianFgyrOver_velocity_{nullptr};

  /** true if jacobianFgyrOver_velocity is not set or required/set and constant */
  bool hasConstantJacobianFgyrOver_velocity_{true};

  /** True if jacobianFgyrOver_velocity is required */
  bool hasJacobianFgyrOver_velocity_{false};

  /** forces(q[0],q[1],t)= fExt - fInt -FGyr */
  std::shared_ptr<siconos::algebra::SiconosVector> totalForces_{nullptr};

  /** jacobian_q forces*/
  std::unique_ptr<siconos::algebra::SiconosDenseMatrix> jacobianTotalForcesOver_q_{nullptr};

  /** jacobian_{velocity} forces*/
  std::unique_ptr<siconos::algebra::SiconosDenseMatrix> jacobianTotalForcesOver_velocity_{
      nullptr};

  /** memory of previous forces of the system */
  siconos::algebra::SiconosMemory totalForcesMemory_;

  /** memory of previous generalized forces due to constraints */
  std::vector<siconos::algebra::SiconosMemory> pMemory_;

  /** Default constructor */
  LagrangianDS() = default;

  /**
   * @brief Utility function providing uniform access to the  \f$ \nabla_qF_{int} \f$ matrix.
   *
   * generic and efficient access to the underlying `SiconosDenseMatrix` object, avoiding any
   * copy.
   *
   * @tparam F a callable type (e.g., lambda) taking a reference to
   *           `SiconosDenseMatrix` or `const SiconosDenseMatrix&`.
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
  decltype(auto) useJacobianFintOver_q(F&& f) const {
    return siconos::algebra::visitStorage(jacobianFintOver_q_storage_, std::forward<F>(f),
                                          "jacobianFintOver_q_storage_");
  }

  template <typename F>
  decltype(auto) useJacobianFintOver_q(F&& f) {
    return siconos::algebra::visitStorage(jacobianFintOver_q_storage_, std::forward<F>(f),
                                          "jacobianFintOver_q_storage_");
  }

  /**
   * @brief Utility function providing uniform access to the  \f$ \nabla_vF_{int} \f$ matrix.
   *
   * generic and efficient access to the underlying `SiconosDenseMatrix` object, avoiding any
   * copy.
   *
   * @tparam F a callable type (e.g., lambda) taking a reference to
   *           `SiconosDenseMatrix` or `const SiconosDenseMatrix&`.
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
  decltype(auto) useJacobianFintOver_velocity(F&& f) const {
    return siconos::algebra::visitStorage(jacobianFintOver_velocity_storage_,
                                          std::forward<F>(f),
                                          "jacobianFintOver_velocity_storage_");
  }

  template <typename F>
  decltype(auto) useJacobianFintOver_velocity(F&& f) {
    return siconos::algebra::visitStorage(jacobianFintOver_velocity_storage_,
                                          std::forward<F>(f),
                                          "jacobianFintOver_velocity_storage_");
  }

  /**
   * @brief Utility function providing uniform access to the  \f$ \nabla_qF_{gyr} \f$ matrix.
   *
   * generic and efficient access to the underlying `SiconosDenseMatrix` object, avoiding any
   * copy.
   *
   * @tparam F a callable type (e.g., lambda) taking a reference to
   *           `SiconosDenseMatrix` or `const SiconosDenseMatrix&`.
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
  decltype(auto) useJacobianFgyrOver_q(F&& f) const {
    return siconos::algebra::visitStorage(jacobianFgyrOver_q_storage_, std::forward<F>(f),
                                          "jacobianFgyrOver_q_storage_");
  }

  template <typename F>
  decltype(auto) useJacobianFgyrOver_q(F&& f) {
    return siconos::algebra::visitStorage(jacobianFgyrOver_q_storage_, std::forward<F>(f),
                                          "jacobianFgyrOver_q_storage_");
  }

  /**
   * @brief Utility function providing uniform access to the  \f$ \nabla_qF_{gyr} \f$ matrix.
   *
   * generic and efficient access to the underlying `SiconosDenseMatrix` object, avoiding any
   * copy.
   *
   * @tparam F a callable type (e.g., lambda) taking a reference to
   *           `SiconosDenseMatrix` or `const SiconosDenseMatrix&`.
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
  decltype(auto) useJacobianFgyrOver_velocity(F&& f) const {
    return siconos::algebra::visitStorage(jacobianFgyrOver_velocity_storage_,
                                          std::forward<F>(f),
                                          "jacobianFgyrOver_velocity_storage_");
  }

  template <typename F>
  decltype(auto) useJacobianFgyrOver_velocity(F&& f) {
    return siconos::algebra::visitStorage(jacobianFgyrOver_velocity_storage_,
                                          std::forward<F>(f),
                                          "jacobianFgyrOver_velocity_storage_");
  }

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
  LagrangianDS(Eigen::Ref<siconos::algebra::SiconosVector> position,
               Eigen::Ref<siconos::algebra::SiconosVector> velocity,
               siconos::algebra::AliasTag);

  /** constructor from initial state and velocity with copy
   *
   *  all attributes will be initialised (deep copy)
   *  from the input vectors/matrices
   *
   *  @param[in] position initial coordinates
   *  @param[in] velocity initial velocity
   *  @param tag Pass siconos::algebra::copy_t to select this overload
   * (rather than alias version)
   */
  LagrangianDS(const siconos::algebra::SiconosVector& position,
               const siconos::algebra::SiconosVector& velocity, siconos::algebra::CopyTag);

  /** destructor */
  virtual ~LagrangianDS() noexcept = default;

  /** reset the state to the initial state */
  void resetToInitialState() override;

  /** allocate (if needed)  and compute rhs and its jacobian.
   *
   *  \param time of initialization
   */
  void initRhs(double time) override;

  /** set nonsmooth input to zero
   *
   *  \param level input-level to be initialized.
   */
  void initializeNonSmoothInput(unsigned int level) override;

  /** update right-hand side for the current state
   *
   *  \param time of interest
   */
  void computeRhs(double time) override;

  /** update  \f$ \nabla_x rhs \f$  for the current state
   *
   *  \param time of interest
   */
  void computeJacobianRhsOver_x(double time) override;

  /** reset non-smooth part of the rhs (i.e. p), for all 'levels' */
  void resetAllNonSmoothParts() override;

  /** set nonsmooth part of the rhs (i.e. p) to zero for a given level
   *
   *  \param level
   */
  void resetNonSmoothPart(unsigned int level) override;

  /** Compute  \f$ F_{total}(v,q,t) \f$
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  virtual void computeTotalForces(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time);

  /** Compute  \f$ \nabla_qF_{total}(v,q,t) \f$
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  virtual void computeJacobianTotalForcesOver_q(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time);

  /** Compute  \f$ \nabla_{\dot q}F_{total}(v,q,t) \f$
   *
   *  \param velocity vector
   *  \param q state
   *  \param time the current time
   */
  virtual void computeJacobianTotalForcesOver_velocity(
      const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
      const Eigen::Ref<const siconos::algebra::SiconosVector>& q, double time);

  /** \return a read-only view on the generalized coordinates vector (size=dimension()) */
  inline const siconos::algebra::ConstMapVectorType q_read() const override {
    return siconos::algebra::ConstMapVectorType(state_q_[0]->data(), state_q_[0]->size());
  }

  /** \return the generalized coordinates of the system (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> q() const override { return state_q_[0]; }

  // /** generalized coordinates of the system (vector of size dimension())
  //  *
  //  *  \return pointer on a siconos::algebra::SiconosVector
  //  */
  inline siconos::algebra::SiconosVector& q_python() const { return *(state_q_[0]); }

  /** \return  a view on initial position vector q0 */
  inline const siconos::algebra::MapVectorType q0_python()  {
    return *(q0_view_);
  }

  /** \return  a read-only view on velocity vector */
  inline const siconos::algebra::ConstMapVectorType velocity_read() const override {
    return siconos::algebra::ConstMapVectorType(state_q_[1]->data(), state_q_[1]->size());
  }

  /** \return the velocity vector (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> velocity() const { return state_q_[1]; }

  // /** \return  a read-only view on velocity vector */
  inline siconos::algebra::SiconosVector& velocity_python() const { return *(state_q_[1]); }

  /**   \return a read - only reference on the initial state vector */
  Eigen::Ref<const siconos::algebra::SiconosVector> q0() const {
    return use_q0(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector>(v); });
  }

  /** \return a read - only reference on the initial velocity */
  Eigen::Ref<const siconos::algebra::SiconosVector> velocity0() const {
    return use_velocity0(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector>(v); });
  }

  /** \return a  view on the initial velocity vector */  
  inline const siconos::algebra::MapVectorType velocity0_python()  {
    return *(velocity0_view_);
  }

  /** \return a read-only view on acceleration vector */
  inline const siconos::algebra::ConstMapVectorType acceleration_read() const override {
    return siconos::algebra::ConstMapVectorType(state_q_[2]->data(), state_q_[2]->size());
  }

  /** \return the acceleration vector (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> acceleration() const override {
    return state_q_[2];
  }

  /**
   * @brief Utility function providing uniform access to the mass matrix.
   *
   * generic and efficient access to the underlying `SiconosDenseMatrix` object, avoiding any
   * copy.
   *
   * @tparam F a callable type (e.g., lambda) taking a reference to
   *           `SiconosDenseMatrix` or `const SiconosDenseMatrix&`.
   *
   * @param f the functor to apply to the mass matrix.
   * @return whatever is returned by @p f.
   *
   * @throws std::runtime_error if no mass matrix (owned or external) is set.
   *
   * @note This function does not copy the matrix; it forwards a reference
   *       to the actual internal or external object.
   */

  template <typename F>
  decltype(auto) use_mass(F&& f) {
    return siconos::algebra::visitStorage(mass_storage_, std::forward<F>(f), "mass_storage_");
  }

  template <typename F>
  decltype(auto) use_mass(F&& f) const {
    return siconos::algebra::visitStorage(mass_storage_, std::forward<F>(f), "mass_storage_");
  }

  /*  \return a read-only reference on the mass matrix */
  Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> mass() const {
    return use_mass([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosDenseMatrix>(M);
    });
  }

  /** \return LU-factorization of the mass (pointer link) */
  inline std::shared_ptr<siconos::algebra::SiconosDenseLUMatrix> LUMass() const {
    return LUMass_;
  }

  /** @brief set a constant mass matrix for the system
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   * @param newValue matrix to be copied. Its shape must match dimension() x dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  void setConstantMass(const siconos::algebra::SiconosDenseMatrix& newValue,
                       siconos::algebra::CopyTag tag);

  /** @brief set a constant mass matrix for the system
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
  void setConstantMass(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
                       siconos::algebra::AliasTag tag);

  /** \return True if the mass matrix has been set (i.e. different from identity) */
  bool hasMass() const { return hasMass_; }

  /** \return True if the mass matrix has been set (i.e. different from identity) and is
   * constant */
  bool hasConstantMass() const { return hasConstantMass_; }

  /** set a user-defined function to compute the mass matrix
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeMassFunction(const siconos::modeling::func_prototypes::FunctionV_M& fct);

  /** to compute the mass matrix operator \f$ M(q) \f$
   *
   *  \param position q vector
   */
  virtual void computeMass(const Eigen::Ref<const siconos::algebra::SiconosVector>& position);

  /** \return  a read-only view on \f$ F_{int}(\dot q, q, t) \f$  */
  inline auto fint() const {
    return siconos::algebra::ConstMapVectorType(fint_->data(), fint_->size());
  }

  /** True if internal forces are taken into account */
  bool hasInternalForces() const { return hasFint_; }

  /** set a user-defined function to compute \f$ F_{int}(\dot q, q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeFintFunction(const siconos::modeling::func_prototypes::FunctionVVS_V& fct);

  /** Update \f$ F_{int}(\dot q, q, t) \f$
   *  \param velocity \f$ \dot q \f$ vector
   *  \param position q vector
   *  \param time the current time
   */
  void computeFint(const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                   const Eigen::Ref<const siconos::algebra::SiconosVector>& position,
                   double time);

  /*  \return a read-only reference on \f$ \nabla_qF_{int} \f$ matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> jacobianFintOver_q() const {
    return useJacobianFintOver_q([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosDenseMatrix>(M);
    });
  }

  /** @brief set a constant \f$ \nabla_qF_{int} \f$
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   * @param newValue matrix to be copied. Its shape must match dimension() x dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  void setConstantJacobianFintOver_q(const siconos::algebra::SiconosDenseMatrix& newValue,
                                     siconos::algebra::CopyTag tag);

  /** @brief set a constant \f$ \nabla_qF_{int} \f$
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
  void setConstantJacobianFintOver_q(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
                                     siconos::algebra::AliasTag tag);

  /** \return True if  \f$ \nabla_qF_{int} \f$ matrix has been set */
  bool hasJacobianFintOver_q() const { return hasJacobianFintOver_q_; }

  /** set a user-defined function to compute \f$ \nabla_qF_{int}(\dot q, q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianFintOver_qFunction(
      const siconos::modeling::func_prototypes::FunctionVVS_M& fct);

  /** to compute  \f$ \nabla_qF_{int}(\dot q, q, t) \f$
   *  \param velocity \f$ \dot q \f$ vector
   *  \param position q vector
   *  \param time the current time
   */
  void computeJacobianFintOver_q(const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
                                 const Eigen::Ref<siconos::algebra::SiconosVector>& position,
                                 double time);

  /*  \return a read-only reference on \f$ \nabla_{\dot q}F_{int} \f$ matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> jacobianFintOver_velocity()
      const {
    return useJacobianFintOver_velocity([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosDenseMatrix>(M);
    });
  }

  /** @brief set a constant \f$ \nabla_{\dot q}F_{int} \f$
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   * @param newValue matrix to be copied. Its shape must match dimension() x dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  void setConstantJacobianFintOver_velocity(
      const siconos::algebra::SiconosDenseMatrix& newValue, siconos::algebra::CopyTag tag);

  /** @brief set a constant \f$ \nabla_{\dot q}F_{int} \f$
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
  void setConstantJacobianFintOver_velocity(
      Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
      siconos::algebra::AliasTag tag);

  /** \return True if \f$ \nabla_{\dot q}F_{int} \f$ matrix has been set */
  bool hasJacobianFintOver_velocity() const { return hasJacobianFintOver_velocity_; }

  /** set a user-defined function to compute \f$ \nabla_{\dot q}F_{int} \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianFintOver_velocityFunction(
      const siconos::modeling::func_prototypes::FunctionVVS_M& fct);

  /** to compute \f$ \nabla_{\dot q}F_{int} \f$
   *
   *  \param velocity \f$ \dot q \f$ vector
   *  \param position q vector
   *  \param time the current time
   */
  void computeJacobianFintOver_velocity(
      const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
      const Eigen::Ref<siconos::algebra::SiconosVector>& position, double time);

  /** \return  a read-only view on \f$ F_{gyr}(\dot q, q) \f$ */
  inline auto fgyr() const {
    return siconos::algebra::ConstMapVectorType(fgyr_->data(), fgyr_->size());
  }

  // No need to add this method: to set a constant fgyr, just use fext
  //  /** set a constant \f$ F_{gyr}\f$ vector
  //  *
  //  *  \param newFgyr gyroscopic forces vector
  //  */
  // void setConstantFgyr(Eigen::Ref<siconos::algebra::SiconosVector> newFgyr);

  /** True if gyroscopic forces are taken into account */
  bool hasGyroscopicForces() const { return hasFgyr_; }

  /** set a user-defined function to compute \f$ F_{gyr}(\dot q, q) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeFgyrFunction(const siconos::modeling::func_prototypes::FunctionVV_V& fct);

  /** Update \f$ F_{gyr}(\dot q, q) \f$
   *  \param velocity \f$ \dot q \f$ vector
   *  \param position q vector
   */
  void computeFgyr(const Eigen::Ref<const siconos::algebra::SiconosVector>& velocity,
                   const Eigen::Ref<const siconos::algebra::SiconosVector>& position);

  /*  \return a read-only reference on  \f$ \nabla_qF_{gyr} \f$ matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> jacobianFgyrOver_q() const {
    return useJacobianFgyrOver_q([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosDenseMatrix>(M);
    });
  }

  /** @brief set a constant \f$ \nabla_qF_{gyr} \f$
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   * @param newValue matrix to be copied. Its shape must match dimension() x dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  void setConstantJacobianFgyrOver_q(const siconos::algebra::SiconosDenseMatrix& newValue,
                                     siconos::algebra::CopyTag tag);

  /** @brief set a constant  \f$ \nabla_qF_{gyr} \f$
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
  void setConstantJacobianFgyrOver_q(Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
                                     siconos::algebra::AliasTag tag);

  /** \return True if  \f$ \nabla_qF_{gyr} \f$ matrix has been set */
  bool hasJacobianFgyrOver_q() const { return hasJacobianFgyrOver_q_; }

  /** set a user-defined function to compute \f$ \nabla_qF_{gyr}(\dot q, q) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianFgyrOver_qFunction(
      const siconos::modeling::func_prototypes::FunctionVV_M& fct);

  /** to compute  \f$ \nabla_qF_{gyr}(\dot q, q, t) \f$
   *  \param velocity \f$ \dot q \f$ vector
   *  \param position q vector
   */
  void computeJacobianFgyrOver_q(const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
                                 const Eigen::Ref<siconos::algebra::SiconosVector>& position);

  /*  \return a read-only reference on \f$ \nabla_{\dot q}F_{gyr} \f$ matrix */
  inline Eigen::Ref<const siconos::algebra::SiconosDenseMatrix> jacobianFgyrOver_velocity()
      const {
    return useJacobianFgyrOver_velocity([](auto const& M) {
      return Eigen::Ref<const siconos::algebra::SiconosDenseMatrix>(M);
    });
  }

  /** @brief set a constant \f$ \nabla_{\dot q}F_{gyr} \f$
   *
   * Warning : deep copy of the provided object into internal attribute
   *
   * @param newValue matrix to be copied. Its shape must match dimension() x dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
   *
   */
  void setConstantJacobianFgyrOver_velocity(
      const siconos::algebra::SiconosDenseMatrix& newValue, siconos::algebra::CopyTag tag);

  /** @brief set a constant  \f$ \nabla_{\dot q}F_{gyr} \f$
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
  void setConstantJacobianFgyrOver_velocity(
      Eigen::Ref<siconos::algebra::SiconosDenseMatrix> newValue,
      siconos::algebra::AliasTag tag);

  /** \return True if \f$ \nabla_{\dot q}F_{gyr} \f$ matrix has been set */
  bool hasJacobianFgyrOver_velocity() const { return hasJacobianFgyrOver_velocity_; }

  /** set a user-defined function to compute \f$ \nabla_{\dot q}F_{gyr} \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianFgyrOver_velocityFunction(
      const siconos::modeling::func_prototypes::FunctionVV_M& fct);

  /** to compute \f$ \nabla_{\dot q}F_{gyr} \f$
   *
   *  \param velocity \f$ \dot q \f$ vector
   *  \param position q vector
   */
  void computeJacobianFgyrOver_velocity(
      const Eigen::Ref<siconos::algebra::SiconosVector>& velocity,
      const Eigen::Ref<siconos::algebra::SiconosVector>& position);

  /*  \return a read-only reference on the mass matrix */
  Eigen::Ref<const siconos::algebra::SiconosVector> fext() const {
    return use_fext(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector>(v); });
  }

  /** @brief set a constant external forces vector
   *
   * Warning : deep copy of the provided vector into internal attribute
   *
   * @param newValue vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
  version)
   *
   */
  void setConstantFext(const siconos::algebra::SiconosVector& newValue,
                       siconos::algebra::CopyTag tag);

  /** @brief set a constant external forces vector
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means:
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   * @param newValue vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   */
  void setConstantFext(Eigen::Ref<siconos::algebra::SiconosVector> newValue,
                       siconos::algebra::AliasTag tag);

  /** True if external forces are taken into account */
  bool hasExternalForces() const { return hasFext_; }

  /** set a user-defined function to compute external forces
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeFextFunction(const siconos::modeling::func_prototypes::FunctionS_V& fct);

  /** Update external forces values
   *
   *  \param time the current time
   */
  void computeFext(double time);

  /**  \return  a read-only view on \f$  F_{total}(v,q,t) \f$ */
  inline auto totalForces() const {
    return siconos::algebra::ConstMapVectorType(totalForces_->data(), totalForces_->size());
  }

  // FP: A pointer version, required only for hem5OSI - TO BE REVIEWED
  /**  \return  \f$  F_{total}(v,q,t) (pointer link)\f$ */
  inline auto totalForcesPtr() const { return totalForces_; }

  /** \return True if total forces are defined */
  bool hasTotalForces() const { return totalForces_ != nullptr; }

  /** \return True if  \f$ \nabla_q F_{total}\f$ is taken into account */
  bool hasJacobianTotalForcesOver_q() const { return jacobianTotalForcesOver_q_ != nullptr; }

  /** \return True if  \f$ \nabla_q F_{total}\f$ is taken into account and is constant */
  bool hasConstantJacobianTotalForcesOver_q() const {
    return hasConstantJacobianFgyrOver_q_ && hasConstantJacobianFintOver_q_;
  }

  /** \return True if  \f$ \nabla_{velocity} F_{total}\f$ is taken into account */
  bool hasJacobianTotalForcesOver_velocity() const {
    return jacobianTotalForcesOver_velocity_ != nullptr;
  }

  /** \return True if  \f$ \nabla_q F_{total}\f$ is taken into account and is constant */
  bool hasConstantJacobianTotalForcesOver_velocity() const {
    return hasConstantJacobianFgyrOver_velocity_ && hasConstantJacobianFintOver_velocity_;
  }

  /** \return a read-only view on \f$  \nabla_qF_{total}(v,q,t) \f$ */
  inline virtual const siconos::algebra::ConstMapType jacobianTotalForcesOver_q() const {
    return siconos::algebra::ConstMapType(jacobianTotalForcesOver_q_->data(),
                                          jacobianTotalForcesOver_q_->rows(),
                                          jacobianTotalForcesOver_q_->cols());
  }

  /** \return a read-only view on \f$ \nabla_{\dot q}F_{total}(v,q,t) \f$ */
  inline virtual const siconos::algebra::ConstMapType jacobianTotalForcesOver_velocity()
      const {
    return siconos::algebra::ConstMapType(jacobianTotalForcesOver_velocity_->data(),
                                          jacobianTotalForcesOver_velocity_->rows(),
                                          jacobianTotalForcesOver_velocity_->cols());
  }

  /** \return last saved (memory) values of the state vector*/
  inline const siconos::algebra::SiconosMemory& qMemory() override { return qMemory_; }

  /** \return last saved (memory) values of the velocity vector*/
  inline const siconos::algebra::SiconosMemory& velocityMemory() { return velocityMemory_; }

  /** \return last saved (memory) values of p[level] vector
   * \param level required index for p
   */
  inline const siconos::algebra::SiconosMemory& pMemory(unsigned int level) {
    return pMemory_[level];
  }

  /** \return last saved (memory) values of the total forces vector*/
  inline const siconos::algebra::SiconosMemory& forcesMemory() override {
    return totalForcesMemory_;
  }

  /** initialize the siconos::algebra::SiconosMemory objects with a positive
   * size.
   *
   *  \param size the size of the siconos::algebra::SiconosMemory. must be >= 0
   */
  void initMemory(unsigned int size) override;

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   *  \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory() override;

  /** To compute the kinetic energy */
  double computeKineticEnergy();

  /** print the data of the dynamical system on the standard output
   */
  void display(bool brief = true) const override;

  /** Computes post-impact velocity, using pre-impact velocity and impulse (p)
   *  value. Used in EventDriven (LsodarOSI->updateState)
   */
  void computePostImpactVelocity();

  /**
     Allocate memory for q[level], level > 1
     Useful for some integrators that need
     q[2] or other coordinates vectors.

     \param level the required level
   */
  void initMemoryForGeneralizedCoordinates(unsigned int level);

  /**
     Allocate memory and lu-factorize the mass of the system.
     Useful for some integrators with system inversion involving the mass
  */
  void init_lu_mass() override;

  // visitors hook
  virtual void accept(dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }
  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};

}  // namespace siconos::modeling
#endif  // LagrangianDS_H
