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

#include <limits>

#include "DynamicalSystem.hpp"
#include "SecondOrderDS.hpp"

namespace siconos::modeling {

/** Pointer to function for plug-in. */
typedef void (*FInt_NE)(double t, double *q, double *v, double *f, unsigned int size_z,
                        double *z);

typedef void (*FExt_NE)(double t, double *f, unsigned int size_z, double *z);

void computeT(std::shared_ptr<siconos::algebra::SiconosVector> q,
              std::shared_ptr<SecondOrderDS::Matrix> T);

/** Compute the force and moment vectors applied to a body with state
 *  q from a force vector at a given position. */
void computeExtForceAtPos(std::shared_ptr<siconos::algebra::SiconosVector> q,
                          bool isMextExpressedInInertialFrame,
                          std::shared_ptr<siconos::algebra::SiconosVector> force,
                          bool forceAbsRef,
                          std::shared_ptr<siconos::algebra::SiconosVector> pos, bool posAbsRef,
                          std::shared_ptr<siconos::algebra::SiconosVector> fExt,
                          std::shared_ptr<siconos::algebra::SiconosVector> mExt,
                          bool accumulate);

/**
    NewtonEuler non linear dynamical systems

    The equations of motion in the Newton-Euler formalism can be stated as

    \f[
    \left\{\begin{array}{rcl}
    M \dot v +  F_{int}(q,v, \Omega, t)&=& F_{ext}(t), \\
    I \dot \Omega + \Omega \wedge I\Omega  + M_{int}(q,v, \Omega, t) &=&
    M_{ext}(t), \\
    \dot q &=& T(q) [ v, \Omega] \\
    \dot R &=& R \tilde \Omega,\quad R^{-1}=R^T,\quad  \det(R)=1 .
    \end{array}\right.
    \f]

    with

    - \f$ x_G,v_G \f$ position and velocity of the center of mass expressed in a
    inertial frame of reference (world frame)
    - \f$ \Omega \f$ angular velocity vector expressed in the body-fixed frame (frame attached
   to the object)
    - \f$ R \f$ rotation matrix form the inertial frame to the
    body-fixed frame \f$ R^{-1}=R^T, \det(R)=1 \f$, i.e \f$ R\in SO^+(3) \f$
    - \f$ M=m\,I_{3\times 3} \f$ diagonal mass matrix with  \f$ m \in \mathbb{R} \f$ the scalar
   mass
    - \f$ I \f$ constant inertia matrix
    - \f$ F_{ext} \f$ and \f$ M_{ext} \f$ are the external applied forces and moment

    In the current implementation, \f$ R \f$ is parametrized by a unit quaternion.

*/
class NewtonEulerDS : public SecondOrderDS {
 protected:
  ACCEPT_SERIALIZATION(NewtonEulerDS);

  /** _twist contains the twist of the Newton Euler dynamical system.
   *  _twist[0:2] : \f$ v_G \in \RR^3 \f$ velocity of the center of mass in
   *  the inertial frame of reference (world frame).
   *  _twist[3:5] : \f$ \Omega\in\RR^3 \f$ angular velocity expressed in the
   *  body-fixed frame
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _twist{nullptr};

  /** Initial twist */
  std::shared_ptr<siconos::algebra::MapVectorType> _twist0{nullptr};
  std::unique_ptr<std::vector<double>> twist0_internal_storage{nullptr};

  /** _q contains the representation of the system
   *   In the current implementation, we have
   *   _q[0:2] : the coordinates of the center of mass expressed
   *      in the inertial frame of reference (world frame)
   *   _q[3:6] : an unit quaternion representing the orientation of the solid.
   *      This unit quaternion encodes the rotation mapping from the inertial
   *   frame of reference to the body-fixed frame
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _q{nullptr};

  /** Dimension of _q, is not necessary equal to _ndof. In our case, _qDim = 7
   * and  _ndof =6*/
  unsigned int _qDim{7};

  /** The time derivative of  \f$ q \f$ ,  \f$ \dot q \f$ */
  std::shared_ptr<siconos::algebra::SiconosVector> _dotq{nullptr};

  /** _acceleration contains the time derivative of the twist
   */
  std::shared_ptr<siconos::algebra::SiconosVector> _acceleration{nullptr};

  /** Memory vectors that stores the values within the time--step */
  siconos::algebra::SiconosMemory _twistMemory;
  siconos::algebra::SiconosMemory _qMemory;
  siconos::algebra::SiconosMemory _forcesMemory;
  siconos::algebra::SiconosMemory _dotqMemory;

  
  /** Inertial matrix
   */
  std::shared_ptr<Matrix> _I{nullptr};

  /** Scalar mass of the system
   */
  double _scalarMass{1.};

  /** Matrix depending on the parametrization of the orientation
   * \f$ v = T(q) \dot q \f$
   */
  std::shared_ptr<Matrix> _T{nullptr};

  /** Time derivative of T.
   *
   * \f$ \dot v = \dot T(q) \dot q + T(q) \ddot q \f$
   */
  std::shared_ptr<Matrix> _Tdot{nullptr};

  /** external forces of the system */
  std::shared_ptr<siconos::algebra::SiconosVector> _fExt{nullptr};

  /** boolean if _fext is constant (set thanks to setFExtPtr for instance)
   * false by default */
  bool _hasConstantFExt{false};

  /** internal forces of the system */
  std::shared_ptr<siconos::algebra::SiconosVector> _fInt{nullptr};

  /** external moment expressed in the inertial frame */
  std::shared_ptr<siconos::algebra::SiconosVector> _mExt{nullptr};

  /** boolean if _mext is constant (set thanks to setMExtPtr for instance)
   * false by default */
  bool _hasConstantMExt{false};

  /** if true, we assume that mExt is given in inertial frame (default false) */
  bool _isMextExpressedInInertialFrame{false};

  ///** external moment expressed in the body-fixed frame  */
  // std::shared_ptr<siconos::algebra::SiconosVector> _mExtBodyFrame;

  /** internal moment of the forces */
  std::shared_ptr<siconos::algebra::SiconosVector> _mInt{nullptr};

  /** jacobian_q FInt  w.r.t q*/
  std::shared_ptr<Matrix> _jacobianFIntq{nullptr};

  /** jacobian_twist FInt  w.r.t the twist*/
  std::shared_ptr<Matrix> _jacobianFInttwist{nullptr};

  /** jacobian_q MInt w.r.t q */
  std::shared_ptr<Matrix> _jacobianMIntq{nullptr};

  /** jacobian_twist MInt  w.r.t the twist*/
  std::shared_ptr<Matrix> _jacobianMInttwist{nullptr};

  /** jacobian_q MExt w.r.t q*/
  std::shared_ptr<Matrix> _jacobianMExtq{nullptr};

  /** gyroscpical moment  */
  std::shared_ptr<siconos::algebra::SiconosVector> _mGyr{nullptr};

  /** jacobian_twist of mGyr w.r.t the twist*/
  std::shared_ptr<Matrix> _jacobianMGyrtwist{nullptr};

  /** wrench (q,twist,t)= [ fExt - fInt ; mExtBodyFrame - mGyr - mInt ]^T */
  std::shared_ptr<siconos::algebra::SiconosVector> _wrench{nullptr};

  /** jacobian_q forces*/
  std::shared_ptr<Matrix> _jacobianWrenchq{nullptr};

  /** jacobian_{twist} forces*/
  std::shared_ptr<Matrix> _jacobianWrenchTwist;

  /** if true, we set the gyroscopic forces equal to 0 (default false) **/
  bool _nullifyMGyr{false};

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianFIntqByFD{true};

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianFInttwistByFD{true};

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianMIntqByFD{true};

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianMInttwistByFD{true};

  /** value of the step in finite difference */
  double _epsilonFD{sqrt(std::numeric_limits<double>::epsilon())};

  /** Plugin to compute strength of external forces */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginFExt{nullptr};

  /** Plugin to compute the external moment expressed in the inertial frame  */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginMExt{nullptr};

  /** Plugin to compute strength of internal forces */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginFInt{nullptr};

  /** Plugin to compute moments of internal forces */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginMInt{nullptr};

  ///** The following code is commented because the jacobian of _mInt and _fInt
  // *  are not yet used by the numerical scheme.
  // *  Will be needed by a fully implicit scheme for instance.
  // */
  /* jacobian_q */
  //  std::shared_ptr<Matrix> _jacobianqmInt;
  /* jacobian_{qDot} */
  //  std::shared_ptr<Matrix> _jacobianqDotmInt;

  /** NewtonEulerDS plug-in to compute \f$ \nabla_qF_{Int}(\dot q, q, t) \f$, id =
   *  "jacobianFIntq"
   *  @param time : current time
   *  @param sizeOfq : size of vector q
   *  @param q : pointer to the first element of q
   *  @param twist : pointer to the first element of twist
   *  @param[in,out] jacob : pointer to the first element of the jacobian
   *  @param  size of vector z
   *  @param[in,out] z  : a vector of user-defined parameters
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJacqFInt{nullptr};

  /** NewtonEulerDS plug-in to compute \f$ \nabla_{\dot q}F_{Int}(\dot q, q,
   *  t) \f$, id = "jacobianFIntTwist"
   *  @param time : current time
   *  @param sizeOfq : size of vector q
   *  @param q : pointer to the first element of q
   *  @param twist : pointer to the first element of twist
   *  @param[in,out] jacob : pointer to the first element of the jacobian
   *  @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJactwistFInt{nullptr};

  /** NewtonEulerDS plug-in to compute \f$ \nabla_qM_{Int}(\dot q, q, t) \f$, id =
   *  "jacobianMInttwist"
   *  @param time : current time
   *  @param sizeOfq : size of vector q
   *  @param q : pointer to the first element of q
   *  @param twist : pointer to the first element of twist
   *  @param[in,out] jacob : pointer to the first element of the jacobian
   *  @param  size of vector z
   *  @param[in,out] z  : a vector of user-defined parameters
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJacqMInt{nullptr};

  /** NewtonEulerDS plug-in to compute \f$ \nabla_{\dot q}M_{Int}(\dot q, q,
   *  t) \f$, id = "jacobianMInttwist"
   *  @param time : current time
   *  @param sizeOfq : size of vector q
   *  @param q : pointer to the first element of q
   *  @param twist : pointer to the first element of twist
   *  @param[in,out] jacob : pointer to the first element of the jacobian
   *  @param  size of vector z
   *  @param[in,out] z  : a vector of user-defined parameters
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJactwistMInt{nullptr};

  // Internal constant, used to easily identify the different blocs in the rhs
  static constexpr auto jacobianXBloc00_ = 0;
  static constexpr auto jacobianXBloc01_ = 1;
  static constexpr auto jacobianXBloc10_ = 2;
  static constexpr auto jacobianXBloc11_ = 3;
  static constexpr auto zeroMatrix_ = 4;
  static constexpr auto zeroMatrixqDim_ = 5;
  static constexpr auto numberOfRhsMatrices_ = 6;

  /** A container of matrices to save matrices that are involed in first order
   *  from of NewtonEulerDS system values (jacobianXBloc10, jacobianXBloc11,
   *  zeroMatrix, idMatrix) No get-set functions at the time. Only used as a
   *  protected member.*/
  std::vector<std::shared_ptr<Matrix>> _rhsMatrices = {
      nullptr, nullptr, nullptr, nullptr};

  /** Default constructor
   */
  NewtonEulerDS();

  /** build all _plugin... PluggedObject */
  void _zeroPlugin() override;

 public:
  // === CONSTRUCTORS - DESTRUCTOR ===

  /** constructor from a minimum set of data
   *
   *  \param position initial coordinates of this DynamicalSystem
   *  \param twist initial twist of this DynamicalSystem
   *  \param mass the mass
   *  \param inertia the inertia matrix
   */
  NewtonEulerDS(std::shared_ptr<siconos::algebra::SiconosVector> position,
                std::shared_ptr<siconos::algebra::SiconosVector> twist, double mass,
                std::shared_ptr<Matrix> inertia);

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
  void computeJacobianRhsx(double time) override;

  /** reset non-smooth part of the rhs (i.e. p), for all 'levels' */
  void resetAllNonSmoothParts() override;

  /** set nonsmooth part of the rhs (i.e. p) to zero for a given level
   *
   *  \param level
   */
  void resetNonSmoothPart(unsigned int level) override;

  // -- forces --
  /** get forces
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> forces() const override {
    return _wrench;
  }

  // -- Jacobian Forces w.r.t q --

  /** \return the jacobian matrix of forces, with respect to q */
  inline std::shared_ptr<Matrix> jacobianqForces() const override {
    return _jacobianWrenchq;
  }

  /** \return the jacobian matrix  of forces with respect to velocity */
  inline std::shared_ptr<Matrix> jacobianvForces() const override {
    return _jacobianWrenchTwist;
  }

  /** Returns dimension of vector q */
  virtual inline unsigned int getqDim() const { return _qDim; }

  // -- q --

  /** get q (pointer link)
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> q() const override { return _q; }

  /** set value of generalized coordinates vector (copy)
   *
   *  \param newValue
   */
  void setQ(const siconos::algebra::SiconosVector &newValue) override;

  /** set value of generalized coordinates vector (pointer link)
   *
   *  \param newPtr
   */
  void setQPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) override;

  /** set initial state (copy)
   *
   *  \param newValue
   */
  void setQ0(const siconos::algebra::SiconosVector &newValue) override;

  /** set initial state (pointer link)
   *
   *  \param newPtr
   */
  void setQ0Ptr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) override;

  // -- twist --

  /** get twist
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> twist() const { return _twist; }

  /** get twist
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   *  this accessor is left to get a uniform access to velocity.
   *  This should be removed with MechanicalDS class
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> velocity() const override {
    return _twist;
  }

  inline std::shared_ptr<siconos::algebra::MapVectorType> twist0() const { return _twist0; }

  inline std::shared_ptr<siconos::algebra::MapVectorType> velocity0() const override {
    return _twist0;
  }
  /** set  velocity (copy)
   *
   *  \param newValue
   */
  void setVelocity(const siconos::algebra::SiconosVector &newValue) override;

  /** set  velocity (pointer link)
   *
   *  \param newPtr
   */
  void setVelocityPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) override;

  /** set initial velocity (copy)
   *
   *  \param newValue
   */
  void setVelocity0(const siconos::algebra::SiconosVector &newValue) override;

  /** set initial velocity (pointer link)
   *
   *  \param newPtr
   */
  void setVelocity0Ptr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) override;

  /** get acceleration (pointer link)
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  std::shared_ptr<siconos::algebra::SiconosVector> acceleration() const override {
    return _acceleration;
  };

  /** default function to compute the mass
   */
  void computeMass() override{};

  /** function to compute the mass
   *
   *  \param position value used to evaluate the mass matrix
   */
  void computeMass(std::shared_ptr<siconos::algebra::SiconosVector> position) override{};

  /** Get the linear velocity in the absolute (inertial) or relative
   *  (body) frame of reference.
   *
   *  \param absoluteRef If true, velocity is returned in the inertial
   *  frame, otherwise velocity is returned in the body frame.
   *  \return A siconos::algebra::SiconosVector of size 3 containing the linear velocity
   *  of this dynamical system.
   */
  std::shared_ptr<siconos::algebra::SiconosVector> linearVelocity(bool absoluteRef) const;

  /** Fill a siconos::algebra::SiconosVector with the linear velocity in the absolute
   *  (inertial) or relative (body) frame of reference.
   *
   *  \param absoluteRef If true, velocity is returned in the inertial
   *  frame, otherwise velocity is returned in the body frame.
   *  \param v A siconos::algebra::SiconosVector of size 3 to receive the linear velocity.
   */
  void linearVelocity(bool absoluteRef, siconos::algebra::SiconosVector &v) const;

  /** Get the angular velocity in the absolute (inertial) or relative
   *  (body) frame of reference.
   *
   *  \param absoluteRef If true, velocity is returned in the inertial
   *  frame, otherwise velocity is returned in the body frame.
   *  \return A siconos::algebra::SiconosVector of size 3 containing the angular velocity
   *  of this dynamical system.
   */
  std::shared_ptr<siconos::algebra::SiconosVector> angularVelocity(bool absoluteRef) const;

  /** Fill a siconos::algebra::SiconosVector with the angular velocity in the absolute
   *  (inertial) or relative (body) frame of reference.
   *
   *  \param absoluteRef If true, velocity is returned in the inertial
   *  frame, otherwise velocity is returned in the body frame.
   *  \param w A siconos::algebra::SiconosVector of size 3 to receive the angular velocity.
   */
  void angularVelocity(bool absoluteRef, siconos::algebra::SiconosVector &w) const;

  // -- p --

  // -- Mass --

  /** get mass value
   *
   *  \return a double
   */
  inline double scalarMass() const { return _scalarMass; };

  /** Modify the scalar mass */
  void setScalarMass(double mass) {
    _scalarMass = mass;
    updateMassMatrix();
  };

  /** \return the inertia matrix
   */
  std::shared_ptr<Matrix> inertia() { return _I; };

  /**
    Modify the inertia matrix (pointer link)

    \param newInertia the new inertia matrix
  */
  void setInertia(std::shared_ptr<Matrix> newInertia) {
    _I = newInertia;
    updateMassMatrix();
  }

  /**
    Modify the inertia matrix.

    \param ix x component
    \param iy y component
    \param iz z component
  */
  void setInertia(double ix, double iy, double iz);

  /** to be called after scalar mass or inertia matrix have changed */
  void updateMassMatrix();

  // -- Fext --
  /** get fExt
   *
   *  \return pointer on a plugged vector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> fExt() const { return _fExt; }

  /** set fExt to pointer newPtr
   *
   *  \param   newPtr a SP to a Simple vector
   */
  inline void setFExtPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
    _fExt = newPtr;
    _hasConstantFExt = true;
  }

  /** get mExt
   *
   *  \return pointer on a plugged vector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> mExt() const { return _mExt; }

  /** set mExt to pointer newPtr
   *
   *  \param newPtr a SP to a Simple vector
   */
  inline void setMExtPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
    _mExt = newPtr;
    _hasConstantMExt = true;
  }

  /** get mGyr
   *
   *  \return pointer on a plugged vector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> mGyr() const { return _mGyr; }

  inline std::shared_ptr<Matrix> T() { return _T; }
  inline std::shared_ptr<Matrix> Tdot() {
    assert(_Tdot);
    return _Tdot;
  }

  inline std::shared_ptr<siconos::algebra::SiconosVector> dotq() { return _dotq; }

  /** get all the values of the state vector q stored in memory
   *
   *  \return a memory
   */
  inline const siconos::algebra::SiconosMemory &qMemory() override { return _qMemory; }

  /** get all the values of the state vector twist stored in memory
   *
   *  \return a memory
   */
  inline const siconos::algebra::SiconosMemory &twistMemory() { return _twistMemory; }
  /** get all the values of the state vector twist stored in memory
   *
   *  \return a memory
   */
  inline const siconos::algebra::SiconosMemory &velocityMemory() override {
    return _twistMemory;
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

  inline const siconos::algebra::SiconosMemory &forcesMemory() override {
    return _forcesMemory;
  }

  inline const siconos::algebra::SiconosMemory &dotqMemory() { return _dotqMemory; }

  /** To compute the kinetic energy
   */
  double computeKineticEnergy();

  // --- miscellaneous ---

  /** print the data to the screen
   */
  void display(bool brief = true) const override;

  //  inline std::shared_ptr<Matrix> jacobianZFL() const { return
  //  jacobianZFL; }

  void setIsMextExpressedInInertialFrame(bool value);

  inline void setNullifyMGyr(bool value) { _nullifyMGyr = value; }

  virtual void normalizeq();

  /**
      Allocate memory for the lu factorization of the mass of the system.
      Useful for some integrators with system inversion involving the mass
  */
  void init_inverse_mass() override;

  /**
      Update the content of the lu factorization of the mass of the system,
      if required.
  */
  void update_inverse_mass() override;

  inline void setComputeJacobianFIntqByFD(bool value) { _computeJacobianFIntqByFD = value; }
  inline void setComputeJacobianFIntvByFD(bool value) {
    _computeJacobianFInttwistByFD = value;
  }
  inline void setComputeJacobianMIntqByFD(bool value) { _computeJacobianMIntqByFD = value; }
  inline void setComputeJacobianMIntvByFD(bool value) {
    _computeJacobianMInttwistByFD = value;
  }

  /** allow to set a specified function to compute _fExt
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeFExtFunction(const std::string &pluginPath, const std::string &functionName);

  /** allow to set a specified function to compute _mExt
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeMExtFunction(const std::string &pluginPath, const std::string &functionName);

  /** set a specified function to compute _fExt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeFExtFunction(FExt_NE fct);

  /** set a specified function to compute _mExt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeMExtFunction(FExt_NE fct);

  /** allow to set a specified function to compute _fInt
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeFIntFunction(const std::string &pluginPath, const std::string &functionName);

  /** allow to set a specified function to compute _mInt
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeMIntFunction(const std::string &pluginPath, const std::string &functionName);

  /** set a specified function to compute _fInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeFIntFunction(FInt_NE fct);

  /** set a specified function to compute _mInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeMIntFunction(FInt_NE fct);

  /** allow to set a specified function to compute the jacobian w.r.t q of the
   *  internal forces
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFIntqFunction(const std::string &pluginPath,
                                       const std::string &functionName);

  /** allow to set a specified function to compute the jacobian following v of
   *  the internal forces w.r.t.
   *
   *  \param pluginPath: the complete path to the plugin
   *  \param functionName: the name of the function to use in this plugin
   */
  void setComputeJacobianFIntvFunction(const std::string &pluginPath,
                                       const std::string &functionName);

  /** set a specified function to compute jacobian following q of the FInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFIntqFunction(FInt_NE fct);

  /** set a specified function to compute jacobian following v of the FInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFIntvFunction(FInt_NE fct);

  /** allow to set a specified function to compute the jacobian w.r.t q of the
   *  internal forces
   *
   *  \param pluginPath: the complete path to the plugin
   *  \param functionName: the name of the function to use in this plugin
   */
  void setComputeJacobianMIntqFunction(const std::string &pluginPath,
                                       const std::string &functionName);
  /** allow to set a specified function to compute the jacobian following v of
   *  the internal forces w.r.t.
   *
   *  \param pluginPath: the complete path to the plugin
   *  \param functionName: the name of the function to use in this plugin
   */
  void setComputeJacobianMIntvFunction(const std::string &pluginPath,
                                       const std::string &functionName);

  /** set a specified function to compute jacobian following q of the FInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianMIntqFunction(FInt_NE fct);

  /** set a specified function to compute jacobian following v of the FInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianMIntvFunction(FInt_NE fct);

  /** function to compute the external forces
   *
   *  \param time the current time
   */
  virtual void computeFExt(double time);

  /** default function to compute the external forces
   *
   *  \param time the current time
   *  \param fExt the computed external force (in-out param)
   */
  virtual void computeFExt(double time, std::shared_ptr<siconos::algebra::SiconosVector> fExt);

  /** function to compute the external moments
   *  The external moments are expressed by default in the body frame, since the
   *  Euler equation for Omega is written in the body--fixed frame. Nevertheless,
   *  if _isMextExpressedInInertialFrame) is set to true, we assume that the
   *  external moment is given in the inertial frame and we perform the rotation
   *  afterwards
   *
   *  \param time the current time
   *  \param mExt the computed external
   *  moment (in-out param)
   */
  virtual void computeMExt(double time, std::shared_ptr<siconos::algebra::SiconosVector> mExt);

  virtual void computeMExt(double time);

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
   *  \param pos A position at which force should be applied.  If nullptr,
   *  the center of mass is assumed.
   *  \param posAbsRef If true, pos is in inertial frame, otherwise it
   *  is in body frame.
   */
  void addExtForceAtPos(std::shared_ptr<siconos::algebra::SiconosVector> force,
                        bool forceAbsRef,
                        std::shared_ptr<siconos::algebra::SiconosVector> pos =
                            std::shared_ptr<siconos::algebra::SiconosVector>(),
                        bool posAbsRef = false);

  void computeJacobianMExtqExpressedInInertialFrameByFD(
      double time, std::shared_ptr<siconos::algebra::SiconosVector> q);
  void computeJacobianMExtqExpressedInInertialFrame(
      double time, std::shared_ptr<siconos::algebra::SiconosVector> q);

  /** default function to compute the internal forces
   *
   *  \param time the current time
   */
  // void computeFInt(double time);

  /** function to compute the internal forces
   *
   *  \param time the current time
   *  \param q
   *  \param v
   */
  void computeFInt(double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
                   std::shared_ptr<siconos::algebra::SiconosVector> v);

  /** default function to compute the internal forces
   *
   *  \param time the current time
   *  \param q
   *  \param v
   *  \param fInt the computed internal force vector
   */
  virtual void computeFInt(double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
                           std::shared_ptr<siconos::algebra::SiconosVector> v,
                           std::shared_ptr<siconos::algebra::SiconosVector> fInt);

  /** default function to compute the internal moments
   *
   *  \param time the current time
   *  \param q
   *  \param v
   */
  void computeMInt(double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
                   std::shared_ptr<siconos::algebra::SiconosVector> v);

  /** default function to compute the internal moments
   *
   *  \param time the current time
   *  \param q
   *  \param v
   *  \param mInt the computed internal moment vector
   */
  virtual void computeMInt(double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
                           std::shared_ptr<siconos::algebra::SiconosVector> v,
                           std::shared_ptr<siconos::algebra::SiconosVector> mInt);

  /** default function to update the plugins functions using a new time:
   *
   *  \param time  the current time
   */
  void updatePlugins(double time) override{};

  /** Allocate memory for forces and its jacobian.
   */
  void init_forces() override;

  /** Default function to compute forces
   *
   *  \param time double, the current time
   */
  virtual void computeForces(double time);

  /** function to compute forces with some specific values for q and twist (ie
   *  not those of the current state). \param time double : the current time
   *
   *  \param q std::shared_ptr<siconos::algebra::SiconosVector>: pointers on q
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>: pointers on twist
   */
  void computeForces(double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
                     std::shared_ptr<siconos::algebra::SiconosVector> twist) override;

  /** Default function to compute the jacobian w.r.t. q of forces
   *
   *  \param time double, the current time
   */
  void computeJacobianqForces(double time) override;

  /** Default function to compute the jacobian w.r.t. v of forces
   *
   *  \param time double, the current time
   */
  void computeJacobianvForces(double time) override;

  /** function to compute gyroscopic forces with some specific values for q and
   *  twist (ie not those of the current state).
   *
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>: pointers on  twist vector
   */
  virtual void computeMGyr(std::shared_ptr<siconos::algebra::SiconosVector> twist);

  /** function to compute gyroscopic forces with some specific values for q and
   *  twist (ie not those of the current state).
   *
   *  \param twist pointer to twist vector
   *  \param mGyr pointer to gyroscopic forces
   */
  virtual void computeMGyr(std::shared_ptr<siconos::algebra::SiconosVector> twist,
                           std::shared_ptr<siconos::algebra::SiconosVector> mGyr);

  /** Default function to compute the jacobian following q of mGyr
   *
   *  \param time the current time
   */
  virtual void computeJacobianMGyrtwist(double time);

  /** Default function to compute the jacobian following q of mGyr
   *  by forward finite difference
   *
   *  \param time the current time
   *  \param q current state
   *  \param twist pointer to twist vector
   */
  virtual void computeJacobianMGyrtwistByFD(
      double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
      std::shared_ptr<siconos::algebra::SiconosVector> twist);

  // /** Default function to compute the jacobian following v of mGyr
  //  *  \param time the current time
  //  */
  // virtual void computeJacobianvForces(double time);

  /** To compute the jacobian w.r.t q of the internal forces
   *
   *  \param time double : the current time
   */
  void computeJacobianFIntq(double time);

  /** To compute the jacobian w.r.t v of the internal forces
   *
   *  \param time double : the current time
   */
  void computeJacobianFIntv(double time);

  /** To compute the jacobian w.r.t q of the internal forces
   *
   *  \param time double
   *  \param position std::shared_ptr<siconos::algebra::SiconosVector>
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>
   */
  virtual void computeJacobianFIntq(double time,
                                    std::shared_ptr<siconos::algebra::SiconosVector> position,
                                    std::shared_ptr<siconos::algebra::SiconosVector> twist);
  /** To compute the jacobian w.r.t q of the internal forces
   *  by forward finite difference
   *
   *  \param time double
   *  \param position std::shared_ptr<siconos::algebra::SiconosVector>
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>
   */
  void computeJacobianFIntqByFD(double time,
                                std::shared_ptr<siconos::algebra::SiconosVector> position,
                                std::shared_ptr<siconos::algebra::SiconosVector> twist);

  /** To compute the jacobian w.r.t. v of the internal forces
   *
   *  \param time double: the current time
   *  \param position std::shared_ptr<siconos::algebra::SiconosVector>
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>
   */
  virtual void computeJacobianFIntv(double time,
                                    std::shared_ptr<siconos::algebra::SiconosVector> position,
                                    std::shared_ptr<siconos::algebra::SiconosVector> twist);

  /** To compute the jacobian w.r.t v of the internal forces
   *  by forward finite difference
   *
   *  \param time double
   *  \param position std::shared_ptr<siconos::algebra::SiconosVector>
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>
   */
  void computeJacobianFIntvByFD(double time,
                                std::shared_ptr<siconos::algebra::SiconosVector> position,
                                std::shared_ptr<siconos::algebra::SiconosVector> twist);

  /** To compute the jacobian w.r.t q of the internal forces
   *
   *  \param time double : the current time
   */
  virtual void computeJacobianMIntq(double time);

  /** To compute the jacobian w.r.t v of the internal forces
   *
   *  \param time double : the current time
   */
  virtual void computeJacobianMIntv(double time);

  /** To compute the jacobian w.r.t q of the internal forces
   *
   *  \param time double : the current time,
   *  \param position std::shared_ptr<siconos::algebra::SiconosVector>
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>
   */
  virtual void computeJacobianMIntq(double time,
                                    std::shared_ptr<siconos::algebra::SiconosVector> position,
                                    std::shared_ptr<siconos::algebra::SiconosVector> twist);

  /** To compute the jacobian w.r.t q of the internal moments
   *  by forward finite difference
   *
   *  \param time double
   *  \param position std::shared_ptr<siconos::algebra::SiconosVector>
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>
   */
  void computeJacobianMIntqByFD(double time,
                                std::shared_ptr<siconos::algebra::SiconosVector> position,
                                std::shared_ptr<siconos::algebra::SiconosVector> twist);

  /** To compute the jacobian w.r.t. v of the internal forces
   *
   *  \param time double: the current time
   *  \param position std::shared_ptr<siconos::algebra::SiconosVector>
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>
   */
  virtual void computeJacobianMIntv(double time,
                                    std::shared_ptr<siconos::algebra::SiconosVector> position,
                                    std::shared_ptr<siconos::algebra::SiconosVector> twist);

  /** To compute the jacobian w.r.t v of the internal moments
   *  by forward finite difference
   *
   *  \param time double
   *  \param position std::shared_ptr<siconos::algebra::SiconosVector>
   *  \param twist std::shared_ptr<siconos::algebra::SiconosVector>
   */
  void computeJacobianMIntvByFD(double time,
                                std::shared_ptr<siconos::algebra::SiconosVector> position,
                                std::shared_ptr<siconos::algebra::SiconosVector> twist);

  virtual void computeT();

  virtual void computeTdot();

  // visitors hook
  void acceptSP(std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const override;
  Type acceptType(types::FindType &ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // NEWTONEULERNLDS_H
