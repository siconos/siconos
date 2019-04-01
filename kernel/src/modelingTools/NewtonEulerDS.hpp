/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "DynamicalSystem.hpp"
#include "SecondOrderDS.hpp"
#include "BoundaryCondition.hpp"

//#define DEBUG_NOCOLOR
//#define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <debug.h>

/** Pointer to function for plug-in. */
typedef void (*FInt_NE)(double t, double* q, double* v, double *f, unsigned int size_z,  double* z);
typedef void (*FExt_NE)(double t, double* f, unsigned int size_z, double *z);

void computeT(SP::SiconosVector q, SP::SimpleMatrix T);

/** Compute the force and moment vectors applied to a body with state
 * q from a force vector at a given position. */
void computeExtForceAtPos(SP::SiconosVector q, bool isMextExpressedInInertialFrame,
                          SP::SiconosVector force, bool forceAbsRef,
                          SP::SiconosVector pos, bool posAbsRef,
                          SP::SiconosVector fExt, SP::SiconosVector mExt,
                          bool accumulate);

/** NewtonEuler non linear dynamical systems

  The equations of motion in the Newton-Euler formalism can be stated as

  \rst

  .. math::
     :nowrap:

      \left\{\begin{array}{rcl}
      M \\dot v +  F_{int}(q,v, \Omega, t)&=& F_{ext}(t), \             \
      I \dot \Omega + \Omega \wedge I\Omega  + M_{int}(q,v, \Omega, t) &=&  M_{ext}(t), \ \
      \dot q &=& T(q) [ v, \Omega] \                                    \
      \dot R &=& R \tilde \Omega,\quad R^{-1}=R^T,\quad  \det(R)=1 .
      \end{array}\right.
 \endrst

 with
  <ul>
  <li> \f$x_G,v_G\f$ position and velocity of the center of mass expressed in a inertial frame of
  reference (world frame) </li>
  <li> \f$\Omega\f$ angular velocity vector expressed in the body-fixed frame (frame attached to the object) </li>
  <li> \f$R\f$ rotation matrix form the inertial frame to the body-fixed frame \f$R^{-1}=R^T, \det(R)=1\f$, i.e \f$ R\in SO^+(3)\f$  </li>
  <li> \f$M=m\,I_{3\times 3}\f$ diagonal mass matrix with  \f$m \in \mathbb{R}\f$ the scalar mass  </li>
  <li> \f$I\f$ constant inertia matrix </li>
  <li> \f$F_{ext}\f$ and \f$ M_{ext}\f$ are the external applied forces and moment  </li>
  </ul>


  In the current implementation, \f$R\f$ is parametrized by a unit quaternion.

*/
class NewtonEulerDS : public SecondOrderDS
{
protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(NewtonEulerDS);

  /** Common code for constructors
   * should be replaced in C++11 by delegating constructors
   */
  void _init();

  // -- MEMBERS --

  /** _twist contains the twist of the Newton Euler dynamical system.
   *  _twist[0:2] : \f$v_G \in \RR^3 \f$ velocity of the center of mass in
   * the inertial frame of reference (world frame).
   *  _twist[3:5] : \f$\Omega\in\RR^3\f$ angular velocity expressed in the body-fixed frame
   */
  SP::SiconosVector _twist;

  /** Initial twist */
  SP::SiconosVector _twist0;

  /** _q contains the representation of the system
   * In the current implementation, we have
   *   _q[0:2] : the coordinates of the center of mass expressed
   *      in the inertial frame of reference (world frame)
   *   _q[3:6] : an unit quaternion representing the orientation of the solid.
   *      This unit quaternion encodes the rotation mapping from the inertial frame of reference
   *      to the body-fixed frame
   */
  SP::SiconosVector _q;

  /** Dimension of _q, is not necessary equal to _ndof. In our case, _qDim = 7 and  _ndof =6*/
  unsigned int _qDim;

  /** The time derivative of \f$q\f$, \f$\dot q\f$*/
  SP::SiconosVector _dotq;

  /** _acceleration contains the time derivative of the twist
   */
  SP::SiconosVector _acceleration;

  /** Memory vectors that stores the values within the time--step */
  SiconosMemory _twistMemory;
  SiconosMemory _qMemory;
  SiconosMemory _forcesMemory;
  SiconosMemory _dotqMemory;

  /** Inertial matrix
   */
  SP::SiconosMatrix _I;

  /** Scalar mass of the system
   */
  double _scalarMass;

  /** Matrix depending on the parametrization of the orientation
   * \f$v = T(q) \dot q\f$
   */
  SP::SimpleMatrix _T;

  /** Time derivative of T.
   *
   * \f$\dot v = \dot T(q) \dot q + T(q) \ddot q\f$
   */
  SP::SimpleMatrix _Tdot;

  /** external forces of the system */
  SP::SiconosVector _fExt;

  /** boolean if _fext is constant (set thanks to setFExtPtr for instance)
   * false by default */
  bool _hasConstantFExt;

  /** internal forces of the system */
  SP::SiconosVector _fInt;

  /** external moment expressed in the inertial frame */
  SP::SiconosVector _mExt;

  /** boolean if _mext is constant (set thanks to setMExtPtr for instance)
   * false by default */
  bool _hasConstantMExt;

  /** if true, we assume that mExt is given in inertial frame (default false)  */
  bool _isMextExpressedInInertialFrame;

  /** external moment expressed in the body-fixed frame  */
  // SP::SiconosVector _mExtBodyFrame;

  /** internal moment of the forces */
  SP::SiconosVector _mInt;

  /** jacobian_q FInt  w.r.t q*/
  SP::SimpleMatrix _jacobianFIntq;

  /** jacobian_twist FInt  w.r.t the twist*/
  SP::SimpleMatrix _jacobianFInttwist;

  /** jacobian_q MInt w.r.t q */
  SP::SimpleMatrix _jacobianMIntq;

  /** jacobian_twist MInt  w.r.t the twist*/
  SP::SimpleMatrix _jacobianMInttwist;

  /** jacobian_q MExt w.r.t q*/
  SP::SimpleMatrix _jacobianMExtq;

  /** gyroscpical moment  */
  SP::SiconosVector _mGyr;

  /** jacobian_twist of mGyr w.r.t the twist*/
  SP::SimpleMatrix _jacobianMGyrtwist;

  /** wrench (q,twist,t)= [ fExt - fInt ; mExtBodyFrame - mGyr - mInt ]^T */

  SP::SiconosVector _wrench;

  /** jacobian_q forces*/
  SP::SimpleMatrix _jacobianWrenchq;

  /** jacobian_{twist} forces*/
  SP::SimpleMatrix _jacobianWrenchTwist;


  /** if true, we set the gyroscopic forces equal to 0 (default false) **/
  bool _nullifyMGyr;

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianFIntqByFD;

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianFInttwistByFD;

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianMIntqByFD;

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianMInttwistByFD;

  /** value of the step in finite difference */
  double _epsilonFD;

  /** Plugin to compute strength of external forces */
  SP::PluggedObject _pluginFExt;

  /** Plugin to compute the external moment expressed in the inertial frame  */
  SP::PluggedObject _pluginMExt;

  /** Plugin to compute strength of internal forces */
  SP::PluggedObject _pluginFInt;

  /** Plugin to compute moments of internal forces */
  SP::PluggedObject _pluginMInt;

  /** The following code is commented because the jacobian of _mInt and _fInt
   *  are not yet used by the numerical scheme.
   *  Will be needed by a fully implicit scheme for instance.
   */
  /** jacobian_q */
  //  SP::SimpleMatrix _jacobianqmInt;
  /** jacobian_{qDot} */
  //  SP::SimpleMatrix _jacobianqDotmInt;

  /** NewtonEulerDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianFIntq"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param twist : pointer to the first element of twist
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacqFInt;

  /** NewtonEulerDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianFIntTwist"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param twist : pointer to the first element of twist
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJactwistFInt;

  /** NewtonEulerDS plug-in to compute \f$\nabla_qM_{Int}(\dot q, q, t)\f$, id = "jacobianMInttwist"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param twist : pointer to the first element of twist
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacqMInt;


  /** NewtonEulerDS plug-in to compute \f$\nabla_{\dot q}M_{Int}(\dot q, q, t)\f$, id = "jacobianMInttwist"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param twist : pointer to the first element of twist
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJactwistMInt;



  enum NewtonEulerDSRhsMatrices {jacobianXBloc00, jacobianXBloc01, jacobianXBloc10, jacobianXBloc11, zeroMatrix, zeroMatrixqDim, numberOfRhsMatrices};

  /** A container of matrices to save matrices that are involed in first order from of
   * NewtonEulerDS system values (jacobianXBloc10, jacobianXBloc11, zeroMatrix, idMatrix)
   * No get-set functions at the time. Only used as a protected member.*/
  VectorOfSMatrices _rhsMatrices;



  /** Default constructor
   */
  NewtonEulerDS();

  /** build all _plugin... PluggedObject */
  void _zeroPlugin();

public:

  // === CONSTRUCTORS - DESTRUCTOR ===

  /** constructor from a minimum set of data
   *  \param position initial coordinates of this DynamicalSystem
   *  \param twist initial twist of this DynamicalSystem
   *  \param mass the mass
   *  \param inertia the inertia matrix
   */
  NewtonEulerDS(SP::SiconosVector position,
                SP::SiconosVector twist,
                double mass,
                SP::SiconosMatrix inertia);

  /** destructor */
  virtual ~NewtonEulerDS();

  /*! @name Right-hand side computation */
  //@{

  /** reset the state to the initial state */
  void resetToInitialState();

  /** allocate (if needed)  and compute rhs and its jacobian.
   * \param time of initialization
   */
  void initRhs(double time) ;

  /** set nonsmooth input to zero
   *  \param level input-level to be initialized.
   */
  void initializeNonSmoothInput(unsigned int level) ;

  /** update right-hand side for the current state
   *  \param time of interest
   */
  virtual void computeRhs(double time);

  /** update \f$\nabla_x rhs\f$ for the current state
   *  \param time of interest
   */
  virtual void computeJacobianRhsx(double time);

  /** reset non-smooth part of the rhs (i.e. p), for all 'levels' */
  void resetAllNonSmoothParts();

  /** set nonsmooth part of the rhs (i.e. p) to zero for a given level
   * \param level
   */
  void resetNonSmoothPart(unsigned int level);

  // -- forces --
  /** get forces
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector forces() const
  {
    return _wrench;
  }

  // -- Jacobian Forces w.r.t q --


  /** get JacobianqForces
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianqForces() const
  {
    return _jacobianWrenchq;
  }

  /** get JacobianvForces
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianvForces() const
  {
    return _jacobianWrenchTwist;
  }

  ///@}

  /*! @name Attributes access
    @{ */


  /** Returns dimension of vector q */
  virtual inline unsigned int getqDim() const
  {
    return _qDim;
  }

  // -- q --

  /** get q (pointer link)
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q() const
  {
    return _q;
  }

  /** set value of generalized coordinates vector (copy)
   *  \param newValue
   */
  virtual void setQ(const SiconosVector& newValue);

  /** set value of generalized coordinates vector (pointer link)
   *  \param newPtr
   */
  virtual void setQPtr(SP::SiconosVector newPtr);

  /** set initial state (copy)
   *  \param newValue
   */
  void setQ0(const SiconosVector& newValue);

  /** set initial state (pointer link)
   *  \param newPtr
   */
  void setQ0Ptr(SP::SiconosVector newPtr);

   // -- twist --

  /** get twist
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector twist() const
  {
    return _twist;
  }

  /** get twist
   *  \return pointer on a SiconosVector
   * this accessor is left to get a uniform access to velocity.
   * This should be removed with MechanicalDS class
   */
  inline SP::SiconosVector velocity() const
  {
    return _twist;
  }

  inline SP::SiconosVector twist0() const
  {
    return _twist0;
  }

  inline SP::SiconosVector velocity0
  () const
  {
    return _twist0;
  }
  /** set  velocity (copy)
   *  \param newValue
   */
  void setVelocity(const SiconosVector& newValue);

  /** set  velocity (pointer link)
   *  \param newPtr
   */
  void setVelocityPtr(SP::SiconosVector newPtr) ;

  /** set initial velocity (copy)
   *  \param newValue
   */
  void setVelocity0(const SiconosVector& newValue);

  /** set initial velocity (pointer link)
   *  \param newPtr
   */
  void setVelocity0Ptr(SP::SiconosVector newPtr) ;

  /** get acceleration (pointer link)
   *  \return pointer on a SiconosVector
   */
  SP::SiconosVector acceleration() const
  {
    return _acceleration;
  };

  /** default function to compute the mass
   */
  virtual void computeMass() {};

  /** function to compute the mass
   *  \param position value used to evaluate the mass matrix
   */
  virtual void computeMass(SP::SiconosVector position) {};

  /** Get the linear velocity in the absolute (inertial) or relative
   * (body) frame of reference.
   * \param absoluteRef If true, velocity is returned in the inertial
   *                    frame, otherwise velocity is returned in the
   *                    body frame.
   * \return A SiconosVector of size 3 containing the linear velocity
   *         of this dynamical system.
   */
  SP::SiconosVector linearVelocity(bool absoluteRef) const;

  /** Fill a SiconosVector with the linear velocity in the absolute
   * (inertial) or relative (body) frame of reference.
   * \param absoluteRef If true, velocity is returned in the inertial
   *                    frame, otherwise velocity is returned in the
   *                    body frame.
   * \param v A SiconosVector of size 3 to receive the linear velocity.
   */
  void linearVelocity(bool absoluteRef, SiconosVector &v) const;

  /** Get the angular velocity in the absolute (inertial) or relative
   * (body) frame of reference.
   * \param absoluteRef If true, velocity is returned in the inertial
   *                    frame, otherwise velocity is returned in the
   *                    body frame.
   * \return A SiconosVector of size 3 containing the angular velocity
   *         of this dynamical system.
   */
  SP::SiconosVector angularVelocity(bool absoluteRef) const;

  /** Fill a SiconosVector with the angular velocity in the absolute
   * (inertial) or relative (body) frame of reference.
   * \param absoluteRef If true, velocity is returned in the inertial
   *                    frame, otherwise velocity is returned in the
   *                    body frame.
   * \param w A SiconosVector of size 3 to receive the angular velocity.
   */
  void angularVelocity(bool absoluteRef, SiconosVector &w) const;

    // -- p --




  // -- Mass --

  /** get mass value
   *  \return a double
   */
  inline double scalarMass() const
  {
    return _scalarMass;
  };

  /** Modify the scalar mass */
  void setScalarMass(double mass)
  {
    _scalarMass = mass;
    updateMassMatrix();
  };

  /* Get the inertia matrix
   * \return a SP::SimpleMatrix
   */
  SP::SiconosMatrix inertia()
  {
    return _I;
  };

  /* Modify the inertia matrix (pointer link)
     \param newInertia the new inertia matrix
  */
  void setInertia(SP::SiconosMatrix newInertia)
  {
    _I = newInertia;
    updateMassMatrix();
  }

  /* Modify the inertia matrix.
     \param ix x component
     \param iy y component
     \param iz z component
  */
  void setInertia(double ix, double iy, double iz);

  /** to be called after scalar mass or inertia matrix have changed */
  void updateMassMatrix();

  // -- Fext --
  /** get fExt
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector fExt() const
  {
    return _fExt;
  }

  /** set fExt to pointer newPtr
   *  \param   newPtr a SP to a Simple vector
   */
  inline void setFExtPtr(SP::SiconosVector newPtr)
  {
    _fExt = newPtr;
    _hasConstantFExt = true;
  }

  /** set mExt to pointer newPtr
    *  \param newPtr a SP to a Simple vector
    */
  inline void setMExtPtr(SP::SiconosVector newPtr)
  {
    _mExt = newPtr;
    _hasConstantMExt = true;
  }

  /** get mGyr
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector mGyr() const
  {
    return _mGyr;
  }

  inline SP::SimpleMatrix T()
  {
    return _T;
  }
  inline SP::SimpleMatrix Tdot()
  {
    assert(_Tdot);
    return _Tdot;
  }

  inline SP::SiconosVector dotq()
  {
    return _dotq;
  }

  /** @} end of members access group. */

  /*! @name Memory vectors management  */
  //@{

  /** get all the values of the state vector q stored in memory
   *  \return a memory
   */
  inline const SiconosMemory& qMemory()
  {
    return _qMemory;
  }


  /** get all the values of the state vector twist stored in memory
   *  \return a memory
   */
  inline const SiconosMemory& twistMemory()
  {
    return _twistMemory;
  }
  /** get all the values of the state vector twist stored in memory
   *  \return a memory
   */
  inline const SiconosMemory& velocityMemory()
  {
    return _twistMemory;
  }

    /** initialize the SiconosMemory objects with a positive size.
   * \param steps the size of the SiconosMemory (i)
   */
  void initMemory(unsigned int steps);

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory();

  inline const SiconosMemory& forcesMemory()
  {
    return _forcesMemory;
  }

  inline const SiconosMemory& dotqMemory()
  {
    return _dotqMemory;
  }


  /** @} end of memory group. */

  /*! @name Miscellaneous public methods */
  //@{

  /** To compute the kinetic energy
   */
  double computeKineticEnergy();

  // --- miscellaneous ---

  /** print the data to the screen
   */
  void display(bool brief = true) const;

  //  inline SP::SiconosMatrix jacobianZFL() const { return jacobianZFL; }

  inline void setIsMextExpressedInInertialFrame(bool value)
  {
    _isMextExpressedInInertialFrame= value;
    if(!_jacobianMExtq)
      _jacobianMExtq.reset(new SimpleMatrix(3, _qDim));
    if(!_jacobianWrenchq)
      _jacobianWrenchq.reset(new SimpleMatrix(_ndof, _qDim));
  }

  inline void setNullifyMGyr(bool value)
  {
    _nullifyMGyr = value;
  }

  virtual void normalizeq();

  /** Allocate memory for the lu factorization of the mass of the system.
      Useful for some integrators with system inversion involving the mass
  */
  void init_inverse_mass();

  /** Update the content of the lu factorization of the mass of the system,
      if required.
  */
  void update_inverse_mass();

  //@}


  /*! @name Plugins management  */

  //@{

  inline void setComputeJacobianFIntqByFD(bool value)
  {
    _computeJacobianFIntqByFD=value;
  }
  inline void setComputeJacobianFIntvByFD(bool value)
  {
    _computeJacobianFInttwistByFD=value;
  }
  inline void setComputeJacobianMIntqByFD(bool value)
  {
    _computeJacobianMIntqByFD=value;
  }
  inline void setComputeJacobianMIntvByFD(bool value)
  {
    _computeJacobianMInttwistByFD=value;
  }


  /** allow to set a specified function to compute _fExt
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeFExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginFExt->setComputeFunction(pluginPath, functionName);
    if(!_fExt)
      _fExt.reset(new SiconosVector(3, 0));
    _hasConstantFExt = false;
  }
  /** allow to set a specified function to compute _mExt
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeMExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginMExt->setComputeFunction(pluginPath, functionName);
    if(!_mExt)
      _mExt.reset(new SiconosVector(3, 0));
    _hasConstantMExt = false;
  }

  /** set a specified function to compute _fExt
   *  \param fct a pointer on the plugin function
   */
  void setComputeFExtFunction(FExt_NE fct)
  {
    _pluginFExt->setComputeFunction((void*)fct);
    if(!_fExt)
      _fExt.reset(new SiconosVector(3, 0));
    _hasConstantFExt = false;
  }

  /** set a specified function to compute _mExt
   *  \param fct a pointer on the plugin function
   */
  void setComputeMExtFunction(FExt_NE fct)
  {
    _pluginMExt->setComputeFunction((void*)fct);
    if(!_mExt)
      _mExt.reset(new SiconosVector(3, 0));
    _hasConstantMExt = false;
  }

  /** allow to set a specified function to compute _fInt
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeFIntFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginFInt->setComputeFunction(pluginPath, functionName);
    if(!_fInt)
      _fInt.reset(new SiconosVector(3, 0));
  }
  /** allow to set a specified function to compute _mInt
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeMIntFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginMInt->setComputeFunction(pluginPath, functionName);
    if(!_mInt)
      _mInt.reset(new SiconosVector(3, 0));
  }

  /** set a specified function to compute _fInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeFIntFunction(FInt_NE fct)
  {
    _pluginFInt->setComputeFunction((void*)fct);
    if(!_fInt)
      _fInt.reset(new SiconosVector(3, 0));
  }

  /** set a specified function to compute _mInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeMIntFunction(FInt_NE fct)
  {
    _pluginMInt->setComputeFunction((void*)fct);
    if(!_mInt)
      _mInt.reset(new SiconosVector(3, 0));
  }

  /** allow to set a specified function to compute the jacobian w.r.t q of the internal forces
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFIntqFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** allow to set a specified function to compute the jacobian following v of the internal forces w.r.t.
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFIntvFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute jacobian following q of the FInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFIntqFunction(FInt_NE fct);

  /** set a specified function to compute jacobian following v of the FInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFIntvFunction(FInt_NE fct);

  /** allow to set a specified function to compute the jacobian w.r.t q of the internal forces
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianMIntqFunction(const std::string&  pluginPath, const std::string&  functionName);
  /** allow to set a specified function to compute the jacobian following v of the internal forces w.r.t.
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianMIntvFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute jacobian following q of the FInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianMIntqFunction(FInt_NE fct);

  /** set a specified function to compute jacobian following v of the FInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianMIntvFunction(FInt_NE fct);

  /**  function to compute the external forces
   * \param time the current time
   */
  virtual void computeFExt(double time);

  /** default function to compute the external forces
    * \param time the current time
    * \param fExt the computed external force (in-out param)
    */
  virtual void computeFExt(double time, SP::SiconosVector fExt);

  /** function to compute the external moments
   * The external moments are expressed by default in the body frame, since the Euler equation for
   * Omega is written in the body--fixed frame.
   * Nevertheless, if _isMextExpressedInInertialFrame) is set to true, we assume that
   * the external moment is given in the inertial frame and we perform the rotation afterwards
   * \param time the current time
   * \param mExt the computed external moment (in-out param)
   */
  virtual void computeMExt(double time, SP::SiconosVector mExt);

  virtual void computeMExt(double time);

  /** Adds a force/torque impulse to a body's FExt and MExt vectors in
   * either absolute (inertial) or relative (body) frame.  Modifies
   * contents of _fExt and _mExt! Therefore these must have been set
   * as constant vectors using setFExtPtr and setMExtPtr prior to
   * calling this function.  Adjustments to _mExt will take into
   * account the value of _isMextExpressedInInertialFrame.
   * \param force A force vector to be added.
   * \param forceAbsRef If true, force is in inertial frame, otherwise
   *                    it is in body frame.
   * \param pos A position at which force should be applied.  If NULL,
   *            the center of mass is assumed.
   * \param posAbsRef If true, pos is in inertial frame, otherwise it
   *                  is in body frame.
   */
  void addExtForceAtPos(SP::SiconosVector force, bool forceAbsRef,
                        SP::SiconosVector pos = SP::SiconosVector(),
                        bool posAbsRef = false);

  void computeJacobianMExtqExpressedInInertialFrameByFD(double time, SP::SiconosVector q);
  void computeJacobianMExtqExpressedInInertialFrame(double time, SP::SiconosVector q);

  /** default function to compute the internal forces
   *  \param time the current time
   */
  // void computeFInt(double time);

  /** function to compute the internal forces
   * \param time the current time
   * \param q
   * \param v
   */
  void computeFInt(double time, SP::SiconosVector q, SP::SiconosVector v);

  /** default function to compute the internal forces
   * \param time the current time
   * \param q
   * \param v
   * \param fInt the computed internal force vector
   */
  virtual void computeFInt(double time, SP::SiconosVector q, SP::SiconosVector v, SP::SiconosVector fInt);

  /** default function to compute the internal moments
   * \param time the current time
   * \param q
   * \param v
   */
  void computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v);

  /** default function to compute the internal moments
   * \param time the current time
   * \param q
   * \param v
   * \param mInt the computed internal moment vector
   */
  virtual void computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v, SP::SiconosVector mInt);

  /**default function to update the plugins functions using a new time:
   * \param time  the current time
   */
  virtual void updatePlugins(double time) {};

  /** Allocate memory for forces and its jacobian.
  */
  void init_forces();

  /** Default function to compute forces
   *  \param time double, the current time
   */
  virtual void computeForces(double time);

  /** function to compute forces with some specific values for q and twist (ie not those of the current state).
   *  \param time double : the current time
   *  \param q SP::SiconosVector: pointers on q
   *  \param twist SP::SiconosVector: pointers on twist
   */
  virtual void computeForces(double time,
                             SP::SiconosVector q,
                             SP::SiconosVector twist);

  /** Default function to compute the jacobian w.r.t. q of forces
   *  \param time double, the current time
   */
  virtual void computeJacobianqForces(double time);

  /** Default function to compute the jacobian w.r.t. v of forces
   *  \param time double, the current time
   */
  virtual void computeJacobianvForces(double time);


  /** function to compute gyroscopic forces with some specific values for q and twist (ie not those of the current state).
   *  \param twist SP::SiconosVector: pointers on  twist vector
   */
  virtual void computeMGyr(SP::SiconosVector twist);

  /** function to compute gyroscopic forces with some specific values for q and twist (ie not those of the current state).
   *  \param twist pointer to twist vector
   *  \param mGyr pointer to gyroscopic forces
   */
  virtual void computeMGyr(SP::SiconosVector twist, SP::SiconosVector mGyr);


  /** Default function to compute the jacobian following q of mGyr
   *  \param time the current time
   */
  virtual void computeJacobianMGyrtwist(double time);

  /** Default function to compute the jacobian following q of mGyr
   * by forward finite difference
   * \param time the current time
   * \param q current state
   * \param twist pointer to twist vector
   */
  virtual void computeJacobianMGyrtwistByFD(double time, SP::SiconosVector q, SP::SiconosVector twist);

  // /** Default function to compute the jacobian following v of mGyr
  //  *  \param time the current time
  //  */
  // virtual void computeJacobianvForces(double time);

  /** To compute the jacobian w.r.t q of the internal forces
   *  \param time double : the current time
   */
  void computeJacobianFIntq(double time);

  /** To compute the jacobian w.r.t v of the internal forces
   *  \param time double : the current time
   */
  void computeJacobianFIntv(double time);

  /** To compute the jacobian w.r.t q of the internal forces
   * \param time double
   * \param position SP::SiconosVector
   * \param twist SP::SiconosVector
   */
  virtual void computeJacobianFIntq(double time,
                                    SP::SiconosVector position,
                                    SP::SiconosVector twist);
  /** To compute the jacobian w.r.t q of the internal forces
   * by forward finite difference
   * \param time double
   * \param position SP::SiconosVector
   * \param twist SP::SiconosVector
   */
  void computeJacobianFIntqByFD(double time,
                                SP::SiconosVector position,
                                SP::SiconosVector twist);


  /** To compute the jacobian w.r.t. v of the internal forces
   *  \param time double: the current time
   * \param position SP::SiconosVector
   * \param twist SP::SiconosVector
   */
  virtual void computeJacobianFIntv(double time,
                                    SP::SiconosVector position,
                                    SP::SiconosVector twist);

  /** To compute the jacobian w.r.t v of the internal forces
   * by forward finite difference
   * \param time double
   * \param position SP::SiconosVector
   * \param twist SP::SiconosVector
   */
  void computeJacobianFIntvByFD(double time,
                                SP::SiconosVector position,
                                SP::SiconosVector twist);


  /** To compute the jacobian w.r.t q of the internal forces
   *  \param time double : the current time
   */
  virtual void computeJacobianMIntq(double time);

  /** To compute the jacobian w.r.t v of the internal forces
   *  \param time double : the current time
   */
  virtual void computeJacobianMIntv(double time);

  /** To compute the jacobian w.r.t q of the internal forces
   *  \param time double : the current time,
   * \param position SP::SiconosVector
   * \param twist SP::SiconosVector
   */
  virtual void computeJacobianMIntq(double time,
                                    SP::SiconosVector position,
                                    SP::SiconosVector twist);

  /** To compute the jacobian w.r.t q of the internal moments
   * by forward finite difference
   * \param time double
   * \param position SP::SiconosVector
   * \param twist SP::SiconosVector
   */
  void computeJacobianMIntqByFD(double time,
                                SP::SiconosVector position,
                                SP::SiconosVector twist);




  /** To compute the jacobian w.r.t. v of the internal forces
   *  \param time double: the current time
   * \param position SP::SiconosVector
   * \param twist SP::SiconosVector
   */
  virtual void computeJacobianMIntv(double time,
                                    SP::SiconosVector position,
                                    SP::SiconosVector twist);

  /** To compute the jacobian w.r.t v of the internal moments
   * by forward finite difference
   * \param time double
   * \param position SP::SiconosVector
   * \param twist SP::SiconosVector
   */
  void computeJacobianMIntvByFD(double time,
                                SP::SiconosVector position,
                                SP::SiconosVector twist);

   virtual void computeT();

  virtual void computeTdot();

 //@}


  ACCEPT_STD_VISITORS();

};
#endif // NEWTONEULERNLDS_H
