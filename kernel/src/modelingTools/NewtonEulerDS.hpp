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

/** \file NewtonEulerDS.hpp
 */

#ifndef NEWTONEULERNLDS_H
#define NEWTONEULERNLDS_H

#include <math.h>

#include "FunctionTypes.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"

namespace siconos::modeling {

/**
    NewtonEuler non linear dynamical systems

    The equations of motion in the Newton-Euler formalism can be stated as

    \f[
    \left\{\begin{array}{rcl}
    M \dot twist =  wrench(twist, q, t) \\
    \dot q &=& T(q) [ v_G, \Omega] \\
    \dot R &=& R \tilde \Omega,\quad R^{-1}=R^T,\quad  \det(R)=1 .
    \end{array}\right.
    \f]

    - \f$ x_G,v_G  \in \RR^3  \f$ position and velocity of the center of mass expressed in the
    inertial frame of reference (world frame)

    - \f$ \Omega \in\RR^3 \f$ angular velocity vector expressed in the body-fixed frame (frame
   attached to the object)

    - \f$ R \f$ rotation matrix from the inertial frame to the
      body-fixed frame \f$ R^{-1}=R^T, \det(R)=1 \f$, i.e \f$ R\in SO^+(3) \f$
      In the current implementation, \f$ R \f$ is parametrized by a unit quaternion.

    - \f$ q =  \left[\begin{array}{c} x_g \\ p \end{array}\right] \in \RR^7\f$
      p being a parametrization, unit quaternion representing the orientation of the solid such
   that

      \f$ R = \Phi(p), \dot p = \psi(p)\Omega \f$

      This unit quaternion encodes the rotation mapping from the inertial
      frame of reference to the body-fixed frame

    - \f$ twist = \left[\begin{array}{c} v_g \\ \Omega \end{array}\right] \in \RR^6\f$


    - \f$  M of size 6x6 \f$ is the total inertia matrix

      \f[ = \left[\begin{array}{cc} mass.I_{3x3}  & 0 \\
               0 &  Inertia \end{array}\right] \f]

      mass being a scalar and Inertia hte constant inertia matrix

    - \f$ wrench(twist, q, t) =

       \left[\begin{array}{c} f_{ext}(t) - f_{int}(twist, q, t) \\
                            -M_{gyr}(twist) + M_{ext}(t) - M_{int}(twist, q, t)
                             \end{array}\right] \in \RR^6\f$


    - \f$ M_{gyr}(twist) = \Omega \wedge Inertia\Omega \f$

    -  \f$ f_{ext}  \in \RR^3  and M_{ext}(t)  \in \RR^3 \f$ are the external applied forces
   and torques

*/
class NewtonEulerDS : public SecondOrderDS {
 protected:
  ACCEPT_SERIALIZATION(NewtonEulerDS);

  /** the twist  \f$= \left[\begin{array}{c} v_g \\ \Omega \end{array}\right] \in \RR^6\f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> twist_{nullptr};

  /** Initial value for the state */
  siconos::algebra::DenseVector7Storage q0_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_q0(F&& f) const {
    return siconos::algebra::visitStorage(q0_storage_, std::forward<F>(f), "q0_storage_");
  }

  /** Initial value for the twist */
  siconos::algebra::DenseVector6Storage twist0_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_twist0(F&& f) const {
    return siconos::algebra::visitStorage(twist0_storage_, std::forward<F>(f),
                                          "twist0_storage_");
  }

  /** vector of coordinates of positions and orientations of the body */
  std::shared_ptr<siconos::algebra::SiconosVector> state_q_{nullptr};

  /** Dimension of state_q_, = 7 because of quaternion parametrization */
  static constexpr siconos::algebra::Index qDim_{7};

  /** The time derivative of  \f$ q \f$ ,  \f$ \dot q(t) \f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> dotq_{nullptr};

  /** time derivative of the twist
   */
  std::shared_ptr<siconos::algebra::SiconosVector> acceleration_{nullptr};

  /** Memory vectors that stores the values within the time--step */
  siconos::algebra::SiconosMemory twistMemory_;
  siconos::algebra::SiconosMemory qMemory_;
  siconos::algebra::SiconosMemory wrenchMemory_;
  siconos::algebra::SiconosMemory dotqMemory_;

  /** Scalar mass of the system */
  double scalarMass_{1.};

  /** inverse or factorization of the mass of the system */
  std::shared_ptr<siconos::algebra::SiconosDenseLUMatrix> LUMass_{nullptr};

  /** total inertia matrix \in \RR\f[ \left[\begin{array}{cc} I_{3x3}  & 0 \\
              0 &  \phi(p) \end{array}\right] \f]

     Matrix depending on the parametrization of the orientation

     \f$ \dot q = T(q) v \∏
  */
  std::shared_ptr<siconos::algebra::SiconosMatrix66> totalInertiaMatrix_{nullptr};

  /** \f[ T(q) = \left[\begin{array}{cc} I_{3x3}  & 0 \\
               0 &  \phi(p) \end{array}\right] \f]

      Matrix depending on the parametrization of the orientation

      \f$ \dot q = T(q) v \
   */
  std::unique_ptr<siconos::algebra::SiconosMatrix76> T_{nullptr};

  /** Time derivative of T(q)
   *
   * \f$ \dot v = \dot T(q) \dot q + T(q) \ddot q \f$
   */
  std::unique_ptr<siconos::algebra::SiconosMatrix76> Tdot_{nullptr};

  /** wrench(twist, q, t)= [ fExt(t) - fInt(twist,q, t) ; mExt(t) - mGyr(twist) -
   * mInt(twist,q,t) ]^T */
  std::shared_ptr<siconos::algebra::SiconosVector6> wrench_{nullptr};

  /** \nabla_{twist} wrench(twist, q, t) */
  std::shared_ptr<siconos::algebra::SiconosMatrix66> jacobianWrenchOver_twist_{nullptr};

  /** \nabla_{q} wrench(twist, q, t) */
  std::shared_ptr<siconos::algebra::SiconosMatrix67> jacobianWrenchOver_q_{nullptr};

  /** external forces applied to the system */
  siconos::algebra::DenseVector3Storage fext_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_fext(F&& f) const {
    return siconos::algebra::visitStorage(fext_storage_, std::forward<F>(f), "fext_storage_");
  }
  template <typename F>
  decltype(auto) use_fext(F&& f) {
    return siconos::algebra::visitStorage(fext_storage_, std::forward<F>(f), "fext_storage_");
  }

  /** function wrapper used to compute external forces \f$f_{ext}(t)\f$ */
  func_prototypes::FunctionS_V computefext_{nullptr};

  /** True if external forces are taken into account and constant */
  bool hasConstantFext_{false};

  /** True if external forces are taken into account */
  bool hasFext_{false};

  /** function wrapper used to compute internal forces \f$f_{int}(twist, q, t)\f$ */
  func_prototypes::FunctionVVS_V computefint_{nullptr};

  /** True if internal forces are  taken into account */
  bool hasFint_{false};

  /** function wrapper used to compute \f$\nabla_{twist}(f_{int}) */
  func_prototypes::FunctionVVS_M computejacobianFintOver_twist_{nullptr};

  /** True if jacobian of fint over twist is required */
  bool hasJacobianFintOver_twist_{false};

  /** function wrapper used to compute \f$\nabla_q(f_{int}) */
  func_prototypes::FunctionVVS_M computejacobianFintOver_q_{nullptr};

  /** True if jacobian of fint over q is required */
  bool hasJacobianFintOver_q_{false};

  /** True to compute \f$\nabla_q(f_{int}) with forward finite differences */
  bool computeJacobianFintOver_q_byFD_{true};

  /** True to compute \f$\nabla_{twist}(f_{int}) with forward finite differences */
  bool computeJacobianFintOver_twist_byFD_{true};

  /** external moments applied to the system */
  siconos::algebra::DenseVector3Storage mext_storage_{std::monostate{}};

  template <typename F>
  decltype(auto) use_mext(F&& f) const {
    return siconos::algebra::visitStorage(mext_storage_, std::forward<F>(f), "mext_storage_");
  }
  template <typename F>
  decltype(auto) use_mext(F&& f) {
    return siconos::algebra::visitStorage(mext_storage_, std::forward<F>(f), "mext_storage_");
  }

  /** function wrapper used to compute external moment \f$m_{ext}(t)\f$ */
  func_prototypes::FunctionS_V computemext_{nullptr};

  /** True if internal forces are taken into account and constant */
  bool hasConstantMext_{false};

  /** True if internal forces are  taken into account */
  bool hasMext_{false};

  /** if true, we assume that mExt is given in inertial frame (default false) */
  bool isMextExpressedInInertialFrame_{false};

  /** if true if gyroscopic forces are taken into account (default true) **/
  bool hasMgyr_{true};

  // Note FP: no used-defined functions for mgyr

  /** function wrapper used to compute internal torques \f$m_{int}(twist, q, t)\f$ */
  func_prototypes::FunctionVVS_V computemint_{nullptr};

  /** True if internal torques are  taken into account */
  bool hasMint_{false};

  /** function wrapper used to compute \f$\nabla_{twist}(m_{int}) */
  func_prototypes::FunctionVVS_M computejacobianMintOver_twist_{nullptr};

  /** function wrapper used to compute \f$\nabla_q(m_{int}) */
  func_prototypes::FunctionVVS_M computejacobianMintOver_q_{nullptr};

  /** True to compute \f$\nabla_q(m_{int}) with forward finite differences */
  bool computeJacobianMintOver_q_byFD_{true};

  /** True to compute \f$\nabla_{twist}(m_{int}) with forward finite differences */
  bool computeJacobianMintOver_twist_byFD_{true};

 public:
  // === CONSTRUCTORS - DESTRUCTOR ===

  /** constructor from a minimum set of data
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * This means that for postion and twist (not inertia: always copied!!)
   *  - ownership stays external
   *  - modifications to the original vector are reflected inside the class
   *
   *  @param[in] position initial coordinates
   *  @param[in] twist initial velocity
   *  @param[in] mass the mass (scalar value)
   *  @param[in] inertia the inertia matrix
   *  @param tag Pass siconos::algebra::alias_t to select this overload
   * (rather than copy version)
   */
  NewtonEulerDS(Eigen::Ref<siconos::algebra::SiconosVector7> position,
                Eigen::Ref<siconos::algebra::SiconosVector6> twist, double mass,
                Eigen::Ref<siconos::algebra::SiconosMatrix33> inertia,
                siconos::algebra::AliasTag tag);
  /** constructor from initial state and velocity with copy
   *
   *  all attributes will be initialised (deep copy)
   *  from the input vectors/matrices
   *
   *  @param[in] position initial coordinates
   *  @param[in] twist initial velocity
   *  @param[in] mass the mass (scalar value)
   *  @param[in] inertia the inertia matrix
   *  @param tag Pass siconos::algebra::copy_t to select this overload
   * (rather than alias version)
   */
  NewtonEulerDS(const siconos::algebra::SiconosVector7& position,
                const siconos::algebra::SiconosVector6& twist, double mass,
                const siconos::algebra::SiconosMatrix33& inertia,
                siconos::algebra::CopyTag tag);

  /** destructor */
  virtual ~NewtonEulerDS() noexcept = default;

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

  /** update \f$ \nabla_x rhs \f$ for the current state
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

  /** \return a read-only view onto the wrench vector */
  inline auto wrench() const {
    return siconos::algebra::ConstMapVector6Type(wrench_->data(), wrench_->size());
  }

  // -- Jacobian Forces w.r.t q --

  /** \return the jacobian matrix of forces, with respect to q */
  inline auto jacobianWrenchOver_q() const {
    return siconos::algebra::ConstMapType(jacobianWrenchOver_q_->data(),
                                          jacobianWrenchOver_q_->rows(),
                                          jacobianWrenchOver_q_->cols());
  }

  /** \return the jacobian matrix  of forces with respect to twist */
  inline auto jacobianWrenchOver_twist() const {
    return siconos::algebra::ConstMapType(jacobianWrenchOver_twist_->data(),
                                          jacobianWrenchOver_twist_->rows(),
                                          jacobianWrenchOver_twist_->cols());
  }

  /** \return dimension of vector q */
  inline auto getqDim() const { return qDim_; }

  /**  \return a read-only view on q state vector */
  inline const siconos::algebra::ConstMapVectorType q_read() const override {
    return siconos::algebra::ConstMapVectorType(state_q_->data(), state_q_->size());
  }

  /** \return the generalized coordinates of the system (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> q() const override { return state_q_; }

  // Only for pybind11
  inline siconos::algebra::SiconosVector& q_python() const { return *(state_q_); }

  /**  \return a read-only view on  $\dot q$ */
  inline auto dotq_read() {
    return siconos::algebra::ConstMapVector7Type(dotq_->data(), dotq_->size());
  }

  /** \return  $\dot q$ (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> dotq() const { return dotq_; }

  /** \return a read-only view on twist = \left[\begin{array}{c} v_g \\ \Omega
   *   \end{array}\right] \in \RR^6\f$ */
  inline auto twist_read() const {
    return siconos::algebra::ConstMapVectorType(twist_->data(), twist_->size());
  }

  /** \return the twist = \left[\begin{array}{c} v_g \\ \Omega \end{array}\right]\f$
   * (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> twist() const { return twist_; }

  // Only for pybind11
  inline siconos::algebra::SiconosVector& twist_python() const { return *(twist_); }

  // FP: override SecondOrderDS. Used only in visitors of MechanicsIO. To be reviewed ...
  inline const siconos::algebra::ConstMapVectorType velocity_read() const override {
    return twist_read();
  }

  /** \return a read-only view onto linear velocity (twist(0:2))*/
  inline auto linearVelocity_read() const {
    return siconos::algebra::ConstMapVector3Type(twist_->data(), 3);
  }

  /** \return a read-only view onto angular velocity (twist(3:6))*/
  inline auto angularVelocity_read() const {
    return siconos::algebra::ConstMapVector3Type(twist_->data() + 3, 3);
  }

  /**   \return a read - only reference on the initial state vector */
  Eigen::Ref<const siconos::algebra::SiconosVector7> q0() const {
    return use_q0(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector7>(v); });
  }

  /** \return a read - only reference on the initial twist vector */
  Eigen::Ref<const siconos::algebra::SiconosVector6> twist0() const {
    return use_twist0(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector6>(v); });
  }

  /** \return a read-only view on acceleration vector */
  inline const siconos::algebra::ConstMapVectorType acceleration_read() const override {
    return siconos::algebra::ConstMapVectorType(acceleration_->data(), acceleration_->size());
  }

  /** \return the acceleration vector (pointer link) */
  std::shared_ptr<siconos::algebra::SiconosVector> acceleration() const override {
    return acceleration_;
  }

  /** \return the linear velocity written in therelative (body) frame */
  siconos::algebra::SiconosVector3 linearVelocityInBodyFrame() const;

  /** \return the angular velocity written in the relative (body) frame */
  siconos::algebra::SiconosVector3 angularVelocityInBodyFrame() const;

  /** \return mass value */
  inline double scalarMass() const { return scalarMass_; };

  /** Modify the scalar mass */
  void setScalarMass(double mass);

  /** \return LU-factorization of the mass (pointer link) */
  inline auto LUMass() const { return LUMass_; }

  /** \return a read-only view on total inertia matrix */
  inline auto totalInertiaMatrix() const {
    return Eigen::Map<const siconos::algebra::SiconosMatrix66>(
        totalInertiaMatrix_->data(), totalInertiaMatrix_->rows(), totalInertiaMatrix_->cols());
  }

  // TEMP FP to be complient with blockCSR
  inline auto totalInertiaMatrixNotCONST() const {
    return Eigen::Map<siconos::algebra::SiconosMatrix66>(
        totalInertiaMatrix_->data(), totalInertiaMatrix_->rows(), totalInertiaMatrix_->cols());
  }

  /** Modify the inertia matrix.

    \param ix x component
    \param iy y component
    \param iz z component
  */
  void setInertia(double ix, double iy, double iz);

  /** \return a read-only view on last computed \f[ T(q) = \left[\begin{array}{cc} I_{3x3}  &
     0
     \\ 0 &  \phi(p) \end{array}\right] \f]T(q)
  */
  inline auto T() const {
    return siconos::algebra::ConstMapType(T_->data(), T_->rows(), T_->cols());
  }

  /** \return a read-only view on last computed \f$ \dot T(q)\f$  */
  inline auto Tdot() {
    return siconos::algebra::ConstMapType(Tdot_->data(), Tdot_->rows(), Tdot_->cols());
  }

  /** \return last saved (memory) values of the state vector*/
  inline const siconos::algebra::SiconosMemory& qMemory() override { return qMemory_; }

  /** \return last saved (memory) values of the twist*/
  inline const siconos::algebra::SiconosMemory& twistMemory() { return twistMemory_; }

  /** \return last saved (memory) values of the wrench*/
  inline const siconos::algebra::SiconosMemory& forcesMemory() override {
    return wrenchMemory_;
  }

  /**   \return a read - only reference on the initial state vector */
  Eigen::Ref<const siconos::algebra::SiconosVector3> fext() const {
    return use_fext(
        [](auto const& v) { return Eigen::Ref<const siconos::algebra::SiconosVector3>(v); });
  }

  /** initialize the siconos::algebra::SiconosMemory objects with a positive size.
   *
   *  \param steps the size of the siconos::algebra::SiconosMemory (i)
   */
  void initMemory(unsigned int steps) override;

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   *  \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory() override;

  inline const siconos::algebra::SiconosMemory& wrenchMemory() { return wrenchMemory_; }

  inline const siconos::algebra::SiconosMemory& dotqMemory() { return dotqMemory_; }

  /** To compute the kinetic energy
   */
  double computeKineticEnergy();

  // --- miscellaneous ---

  /** print the data to the screen
   */
  void display(bool brief = true) const override;

  /** to set which frame is used to write mext
   *  \param value  true when mext is expressed in the inertial frame
   */
  void setIsMextExpressedInInertialFrame(bool value);

  /** To switch on or off Mgyr component
   * \param value true to use Mgyr (which is the default)
   */
  inline void takeIntoAccounMGyr(bool value) { hasMgyr_ = value; }

  void normalizeq();

  /**
      Allocate memory and lu-factorize the mass of the system.
      Useful for some integrators with system inversion involving the mass
  */
  void init_lu_mass() override;

  /** @brief set a constant external forces vector
   *
   * Warning : deep copy of the provided vector into internal attribute
   *
   * @param newValue vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
  version)
   *
   */
  void setConstantFext(const siconos::algebra::SiconosVector3& newValue,
                       siconos::algebra::CopyTag tag);

  /** set a constant external forces vector
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * @param newValue vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   *  \param newFext external forces vector
   */
  void setConstantFext(Eigen::Ref<siconos::algebra::SiconosVector3> newFext,
                       siconos::algebra::AliasTag tag);

  /** True if external forces are taken into account */
  bool hasExternalForces() const { return hasFext_; }

  /** set a user-defined function to compute external forces
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeFextFunction(const func_prototypes::FunctionS_V& fct);

  /** @brief set a constant external moment vector
   *
   * Warning : deep copy of the provided vector into internal attribute
   *
   * @param newValue vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::copy_t to select this overload (rather than alias
  version)
   *
   */
  void setConstantMext(const siconos::algebra::SiconosVector3& newValue,
                       siconos::algebra::CopyTag tag);

  /** set a constant external moment vector
   *
   * Warning : This method does NOT copy the data. Instead, it creates an Eigen::Map
   * pointing directly to the memory provided by the argument.
   *
   * @param newValue vector to be copied. Its size must match dimension()
   * @param tag Pass siconos::algebra::alias_t to select this overload
   *        (rather than copy version)
   *
   *  \param newFext external forces vector
   */
  void setConstantMext(Eigen::Ref<siconos::algebra::SiconosVector3> newFext,
                       siconos::algebra::AliasTag tag);

  /** True if external moments are taken into account */
  bool hasExternalMoment() const { return hasMext_; }

  /** set a user-defined function to compute external moment
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeMextFunction(const func_prototypes::FunctionS_V& fct);

  /** set a user-defined function to compute \f$ f_{int}(twist, q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeFintFunction(const func_prototypes::FunctionVVS_V& fct);

  /** set a user-defined function to compute \f$ \nabla_q f_{int}(twist, q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianFintOver_qFunction(const func_prototypes::FunctionVVS_M& fct);

  /** set a user-defined function to compute \f$ \nabla_{twist} f_{int}(twist, q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianFintOver_twistFunction(const func_prototypes::FunctionVVS_M& fct);

  /** set a user-defined function to compute \f$ m_{int}(twist, q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeMintFunction(const func_prototypes::FunctionVVS_V& fct);

  /** set a user-defined function to compute \f$ \nabla_q m_{int}(twist, q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianMintOver_qFunction(const func_prototypes::FunctionVVS_M& fct);

  /** set a user-defined function to compute \f$ \nabla_{twist} m_{int}(twist, q, t) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianMintOver_twistFunction(const func_prototypes::FunctionVVS_M& fct);

  /** if true, use finite-differences to compute \f$ \nabla_q f_{int}(twist, q, t) \f$  */
  void setComputeJacobianFintOver_q_byFD(bool value);

  /** if true, use finite-differences to compute \f$ \nabla_{twist} f_{int}(twist, q, t) \f$
   */
  void setComputeJacobianFintOver_twist_byFD(bool value);

  /** if true, use finite-differences to compute \f$ \nabla_q m_{int}(twist, q, t) \f$  */
  void setComputeJacobianMintOver_q_byFD(bool value);

  /** if true, use finite-differences to compute \f$ \nabla_{twist} m_{int}(twist, q, t) \f$
   */
  void setComputeJacobianMintOver_twist_byFD(bool value);

  /** True if \nabla_{twist}wrench(v,q,t) is defined */
  bool hasJacobianWrenchOver_twist() const { return jacobianWrenchOver_twist_ != nullptr; }

  /** True if \nabla_{q}wrench(v,q,t) is defined */
  bool hasJacobianWrenchOver_q() const { return jacobianWrenchOver_q_ != nullptr; }

  /**
      compute wrench(twist, q, t) =

       \left[\begin{array}{c} f_{ext}(t) - f_{int}(twist, q, t) \\
                            -M_{gyr}(twist) + M_{ext}(t) - M_{int}(twist, q, t)
                             \end{array}\right] \in \RR^6\f$Compute  \f$ F_{total}(v,q,t) \f$

      \param twist vector
      \param q state
      \param time the current time
  */
  void computeWrench(const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
                     const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time);

  /** Compute  \f$ \nabla_{twist}wrench(v,q,t) \f$
   *
   *  \param twist vector
   *  \param q state
   *  \param time the current time
   */
  void computeJacobianWrenchOver_twist(
      const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time);

  /** Compute  \f$ \nabla_q wrench(v,q,t) \f$
   *
   *  \param twist vector
   *  \param q state
   *  \param time the current time
   */
  void computeJacobianWrenchOver_q(
      const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time);

  /** To compute the jacobian w.r.t q of the internal forces
   *  by forward finite difference
   *
   *  \param time double
   *  \param position const Eigen::Ref<siconos::algebra::SiconosVector> &
   *  \param twist const Eigen::Ref<siconos::algebra::SiconosVector> &
   */
  void computeJacobianFIntqByFD(double time,
                                const Eigen::Ref<siconos::algebra::SiconosVector7>& position,
                                const Eigen::Ref<siconos::algebra::SiconosVector6>& twist);

  /** To compute the jacobian w.r.t v of the internal forces
   *  by forward finite difference
   *
   *  \param time double
   *  \param position const Eigen::Ref<siconos::algebra::SiconosVector> &
   *  \param twist const Eigen::Ref<siconos::algebra::SiconosVector> &
   */
  void computeJacobianFIntvByFD(double time,
                                const Eigen::Ref<siconos::algebra::SiconosVector7>& position,
                                const Eigen::Ref<siconos::algebra::SiconosVector6>& twist);

  void computeT(const Eigen::Ref<const siconos::algebra::SiconosVector7>& q);

  virtual void computeTdot();

  /** Adds a force/torque impulse to a body's FExt and MExt vectors in
   *  either absolute (inertial) or relative (body) frame.  Modifies
   *  contents of _fExt and _mExt! Therefore these must have been set
   *  as constant vectors using setFExtPtr and setMExtPtr prior to
   *  calling this function.  Adjustments to _mExt will take into
   *  account the value of _isMextExpressedInInertialFrame.
   *
   *  \param force A force vector to be added.
   *  \param forceAbsRef If true, force is in inertial frame, otherwise
   *  it is in body frame.
   *  \param pos A position at which force should be applied.
   *  \param posAbsRef If true, pos is in inertial frame, otherwise it
   *  is in body frame.
   */
  void addExtForceAtPos(const Eigen::Ref<siconos::algebra::SiconosVector3>& force,
                        bool forceAbsRef,
                        const Eigen::Ref<siconos::algebra::SiconosVector3>& pos,
                        bool posAbsRef = false);

  virtual void accept(dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};

///////// Free functions /////////

namespace newton_euler {
/** Update the non constant part of \f[ T(q) = \left[\begin{array}{cc} I_{3x3}  & 0 \\
               0 &  \phi(p) \end{array}\right] \f]

    \param[in] q coordinates vector
    \param[in,out] T result

*/
void computeT(const Eigen::Ref<const siconos::algebra::SiconosVector7>& q,
              Eigen::Ref<siconos::algebra::SiconosMatrix76> T);

/** Compute the moment vectors applied to a body with state
 *  q from a force vector at a given position. */

/** To compute the moment vectors applied to a body with state
 *  q from a force vector at a given position.
 *
 *  \param[in] q state vector
 *  \param[in] isMextExpressedInInertialFrame true if mext is defined in inertial frame
 *  \param[in] force a force vector to be added.
 *  \param[in] forceAbsRef if true, force is in inertial frame, otherwise
 *   it is in body frame.
 *  \param[in] pos
 *  \param[in] posAbsRef If true, pos is in inertial frame, otherwise it is in body frame.
 *  \param[in,out] result resulting momemt
 *  \param[in] accumulate if true, mext += result, else mext is reinitialized
 */
void computeMextForceAtPos(const Eigen::Ref<siconos::algebra::SiconosVector7>& q,
                           bool isMextExpressedInInertialFrame,
                           const Eigen::Ref<siconos::algebra::SiconosVector3>& force,
                           bool forceAbsRef,
                           const Eigen::Ref<siconos::algebra::SiconosVector3>& pos,
                           bool posAbsRef, Eigen::Ref<siconos::algebra::SiconosVector3> result,
                           bool accumulate);

/** To add a force impulse to body's forces vector in
 *  either absolute (inertial) or relative (body) frame.
 *
 *  \param[in] q state vector
 *  \param[in] force a force vector to be added.
 *  \param[in] forceAbsRef if true, force is in inertial frame, otherwise
 *   it is in body frame.
 *  \param[in, out] fext resulting force
 *  \param[in] accumulate if true, fext += result, else fext is reinitialized
 */
void computeFextForceAtPos(const Eigen::Ref<siconos::algebra::SiconosVector7>& q,
                           const Eigen::Ref<siconos::algebra::SiconosVector3>& force,
                           bool forceAbsRef, Eigen::Ref<siconos::algebra::MapVector3Type> fExt,
                           bool accumulate);

/** Compute \f$ \nabla_{twist} m_{gyr}(twist) \f$
 *  we compute only the sub-block [:,3:6]
 *  \param[in] twist vector (size 6) = [vg, Omega]^T
 *  \param[in] inertia total inertia matrix, size=6x6
 *  \param[in,out] result the computed jacobian matrix (size 3x6)
 */
void computeJacobianMGyrOver_twist(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosMatrix66>& inertia,
    Eigen::Ref<siconos::algebra::SiconosMatrix36> result);

/** Compute \f$ \nabla_{twist} m_{gyr}(twist) \f$ by forward finite difference
 *  we compute only the sub-block [:,3:6]
 *
 *  \param[in] twist vector (size 6) = [vg, Omega]^T
 *  \param[in] epsilonFD FD step
 *  \param[in] inertia matrix (I), size=3x3
 *  \param[in] computeMgyr function used to compute \f$ m_{gyr}(twist) \f$
 *  \param[in,out] result the computed jacobian matrix (size 3x6)
 */
void computeJacobianMGyrOver_twist_byFD(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist, double epsilonFD,
    const Eigen::Ref<const siconos::algebra::SiconosMatrix66>& inertia,
    const siconos::modeling::func_prototypes::FunctionMV_V& computeMgyr,
    Eigen::Ref<siconos::algebra::SiconosMatrix36> result);

/** Compute \f$ \nabla_{twist} f(twist, q, time) \f$ by forward finite difference
 *   f = fint, mint ...
 *  \param[in] twist vector (size 6) = [vg, Omega]^T
 *  \param[in] q vector q =  \left[\begin{array}{c} x_g \\ p \end{array}\right]
 *  \param[in] time current time
 *  \param[in] epsilonFD FD step
 *  \param[in] computeMint function used to compute \f$ f(twist,q, time) \f$
 *  \param[in,out] result the computed jacobian matrix (size: 3xsize twist)
 */
void computeJacobianFOver_twist_byFD(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time, double epsilonFD,
    const siconos::modeling::func_prototypes::FunctionVVS_V& computeMint,
    Eigen::Ref<siconos::algebra::SiconosMatrix36> result);

/** Compute \f$ \nabla_{q} f(twist, q, time) \f$ by forward finite difference
 *   f = fint, mint ...
 *  \param[in] twist vector (size 6) = [vg, Omega]^T
 *  \param[in] q vector q =  \left[\begin{array}{c} x_g \\ p \end{array}\right]
 *  \param[in] time current time
 *  \param[in] epsilonFD FD step
 *  \param[in] computeMint function used to compute \f$ f(twist,q, time) \f$
 *  \param[in,out] result the computed jacobian matrix (size: 3xsize q)
 */
void computeJacobianFOver_q_byFD(
    const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
    const Eigen::Ref<const siconos::algebra::SiconosVector7>& q, double time, double epsilonFD,
    const siconos::modeling::func_prototypes::FunctionVVS_V& computeMint,
    Eigen::Ref<siconos::algebra::SiconosMatrix37> result);

/** Compute \f$\nabla_q(m_{ext})\f$, required when mext is expressed in the inertial frame.
 *  \param[in] q vector q =  \left[\begin{array}{c} x_g \\ p \end{array}\right]
 *  \param[in] time current time
 *  \param[in] computeMext function used to compute \f$ m_{ext}(time) \f$
 *  \param[in,out] result the computed jacobian matrix (size: 3x7)
 */
void computeJacobianMExtqExpressedInInertialFrame(
    const Eigen::Ref<siconos::algebra::SiconosVector7>& q, double time,
    const siconos::algebra::SiconosVector3& mext,
    Eigen::Ref<siconos::algebra::SiconosMatrix37> result);

/** Compute \f$\nabla_q(m_{ext})\f$, required when mext is expressed in the inertial frame.
 *  \param[in] q vector q =  \left[\begin{array}{c} x_g \\ p \end{array}\right]
 *  \param[in] time current time
 *  \param[in] computeMext function used to compute \f$ m_{ext}(time) \f$
 *  \param[in] isMextExpressedInInertialFrame true if Mext is ... expressed in the inertial
 * frame \param[in] epsilonFD FD step \param[in,out] result the computed jacobian matrix
 * (size: 3x3)
 *
 */
void computeJacobianMExtqExpressedInInertialFrameByFD(
    const Eigen::Ref<siconos::algebra::SiconosVector7>& q, double time,
    const func_prototypes::FunctionS_V& computeMext, bool isMextExpressedInInertialFrame,
    double epsilonFD, Eigen::Ref<siconos::algebra::SiconosMatrix33> result);

/** function to compute gyroscopic forces
 *
 *  \param[in] twist vector (size 6) = [vg, Omega]^T
 *  \param[in] inertia total inertia matrix, size=6x6
 *  \param[in,out] result the computed jacobian matrix (size 3x6)
 */
void computeMgyr(const Eigen::Ref<const siconos::algebra::SiconosVector6>& twist,
                 const Eigen::Ref<const siconos::algebra::SiconosMatrix66>& inertia,
                 Eigen::Ref<siconos::algebra::SiconosVector3> mGyr);
}  // namespace newton_euler

}  // namespace siconos::modeling
#endif  // NEWTONEULERNLDS_H
