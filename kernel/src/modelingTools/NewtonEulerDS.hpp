/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/** \file NewtonEulerDS.hpp
 */

#ifndef NEWTONEULERNLDS_H
#define NEWTONEULERNLDS_H

#include "DynamicalSystem.hpp"
#include "BoundaryCondition.hpp"

/** Pointer to function for plug-in. */
typedef void (*FInt_NE)(double t, double* q, double* v, double *f, unsigned int size_z,  double* z);
typedef void (*FExt_NE)(double t, double* f, unsigned int size_z, double *z);


void computeMObjToAbs(SP::SiconosVector q, SP::SimpleMatrix mObjToAbs);
void computeT(SP::SiconosVector q, SP::SimpleMatrix T);
/** \class NewtonEulerDS
 *  \brief NewtonEuler non linear dynamical systems - Second Order Non Linear Dynamical Systems.
 *   NewtonEuler non linear dynamical systems - Derived from DynamicalSystem -
 *
 * The equations of motion in the Newton-Euler formalism can be stated as
 * \f{equation}
 * \label{eq:NewtonEuler}
 * \left\{\begin{array}{rcl}
 *   M \dot v +  F_{int}(q,v, \Omega, t)&=& F_{ext}(t), \\
 *   I \dot \Omega + \Omega \wedge I\Omega  + M_{int}(q,v, \Omega, t) &=&  M_{ext}(t), \\
 *   \dot q &=& T(q) [ v, \Omega] \\
 *   \dot R &=& R \tilde \Omega,\quad R^{-1}=R^T,\quad  \det(R)=1 .
 * \end{array}\right.
 * \f}
 * with
 * <ul>
 * <li> \f$x_G,v_G\f$ position and velocity of the center of mass expressed in a inertial frame of
 * reference (world frame) </li>
 * <li> \f$\Omega\f$ angular velocity vector expressed in the body-fixed frame (frame attached to the object) </li>
 * <li> \f$R\f$ rotation matrix form the inertial frame to the bosy-fixed frame \f$R^{-1}=R^T, \det(R)=1\f$, i.e \f$ R\in SO^+(3)\f$  </li>
 * <li> \f$M=m\,I_{3\times 3}\f$ diagonal mass matrix with  \f$m \in \mathbb{R}\f$ the scalar mass  </li>
 * <li> \f$I\f$ constant inertia matrix </li>
 * <li> \f$F_{ext}\f$ and \f$ M_{ext}\f$ are the external applied forces and torques  </li>
 * </ul>
 *
 *
 * In the current implementation, \f$R\f$ is parametrized by a unit quaternion.
 *
 */
class NewtonEulerDS : public DynamicalSystem
{
public:
  enum WorkNames {xfree, sizeWorkVec};
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEulerDS);

  void internalInit(SP::SiconosVector Q0, SP::SiconosVector Velocity0, double mass , SP::SiconosMatrix inertialMatrix);

  // -- MEMBERS --

  /** _v contains the velocity of the Newton Euler dynamical system.
   *  _v[0:2] : \f$v_G \in \RR^3 \f$ velocity of the center of mass in
   * the inertial frame of reference (world frame).
   *  _v[3:5] : \f$\Omega\in\RR^3\f$ angular velocity expressed in the body-fixed frame
   */
  SP::SiconosVector _v;

  /** Initial velocity */
  SP::SiconosVector _v0;

  /** Memory vectors that stores the values within the time--step */
  SP::SiconosMemory _vMemory;
  SP::SiconosMemory _qMemory;
  SP::SiconosMemory _forcesMemory;
  SP::SiconosMemory _dotqMemory;

  /** _q dimension, is not necessary _n. In our case, _qDim = 7 and  _n =6*/
  unsigned int _qDim;

  /** _q contains the representation of the system
   * In the current implementation, we have
   *   _q[0:2] : the coordinates of the center of mass expressed
   *      in the inertial frame of reference (world frame)
   *   _q[3:6] : an unit quaternion representing the orientation of the solid.
   *      This unit quaternion encodes the rotation mapping from the inertial frame of reference
   *      to the body-fixed frame
   */
  SP::SiconosVector _q;

  //SP::SiconosVector _deltaq;

  /** Initial position */
  SP::SiconosVector _q0;

  /** The time derivative of \f$q\f$, \f$\dot q\f$*/
  SP::SiconosVector _dotq;

  /* the rotation matrix that converts a vector in body coordinates (in the body fixed frame)
   * in the absolute coordinates in the inertial frame of reference.
   */
  SP::SimpleMatrix _MObjToAbs;

  /** Inertial matrix
   */
  SP::SiconosMatrix _I;

  /** Scalar mass of the system
   */
  double _mass;

  /** used for concatenate _I and _mass.I_3 */
  SP::SimpleMatrix _massMatrix;

  /** Contains the LU factorization of the Mass (or the iteration matrix.).
   */
  SP::SimpleMatrix _luW;

  /** Matrix depending on the parametrization of the orientation
   * \f$v = T(q) \dot q\f$
   */
  SP::SimpleMatrix _T;

  /** Time derivative of T.
   *
   * \f$\dot v = \dot T(q) \dot q + T(q) \ddot q\f$
   */
  SP::SimpleMatrix _Tdot;

  /** "Reaction" due to the non smooth law - The index corresponds to the dynamic levels. */
  std::vector<SP::SiconosVector> _p;

  /** external forces of the system */
  SP::SiconosVector _fExt;

  /** internal forces of the system */
  SP::SiconosVector _fInt;

  /** external moment of the forces */
  SP::SiconosVector _mExt;

  /** internal moment of the forces */
  SP::SiconosVector _mInt;

  /** jacobian_q FInt*/
  SP::SimpleMatrix _jacobianFIntq;

  /** jacobian_{v} FInt*/
  SP::SimpleMatrix _jacobianFIntv;

  /** jacobian_q MInt*/
  SP::SimpleMatrix _jacobianMIntq;

  /** jacobian_{v} MInt*/
  SP::SimpleMatrix _jacobianMIntv;

  /** internal forces of the system */
  SP::SiconosVector _fGyr;

  /** jacobian_v FGyr*/
  SP::SimpleMatrix _jacobianFGyrv;

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianFIntqByFD;

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianFIntvByFD;

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianMIntqByFD;

  /** If true, we compute the missing Jacobian by forward finite difference */
  bool _computeJacobianMIntvByFD;

  /** value of the step in finite difference */
  double _epsilonFD;

  /** Plugin to compute strength of external forces */
  SP::PluggedObject _pluginFExt;

  /** Plugin to compute moments of external forces */
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
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacqFInt;

  /** NewtonEulerDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianFIntv"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacvFInt;

  /** NewtonEulerDS plug-in to compute \f$\nabla_qM_{Int}(\dot q, q, t)\f$, id = "jacobianMIntq"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacqMInt;


  /** NewtonEulerDS plug-in to compute \f$\nabla_{\dot q}M_{Int}(\dot q, q, t)\f$, id = "jacobianMIntv"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacvMInt;

  /** forces(q,v,t)= fExt - fInt - fGyr */
  SP::SiconosVector _forces;

  /** jacobian_q forces*/
  SP::SimpleMatrix _jacobianqForces;

  /** jacobian_{v} forces*/
  SP::SimpleMatrix _jacobianvForces;

  /** Boundary condition applied to a dynamical system*/
  SP::BoundaryCondition _boundaryConditions;

  /** Reaction to an applied  boundary condition */
  SP::SiconosVector _reactionToBoundaryConditions;



  /** set links with DS members
   */
  void connectToDS();

  /** Default constructor
   */
  NewtonEulerDS();

  void zeroPlugin();

public:

  // === CONSTRUCTORS - DESTRUCTOR ===


  /** constructor from a minimum set of data
   *  \param position initial coordinates of this DynamicalSystem
   *  \param velocity initial velocity of this DynamicalSystem
   *  \param mass the mass
   *  \param inertia the inertia matrix
   */
  NewtonEulerDS(SP::SiconosVector position,
                SP::SiconosVector velocity,
                double mass,
                SP::SiconosMatrix inertia);




  /** destructor */
  virtual ~NewtonEulerDS();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

  /** allocate memory for forces and its jacobians, if required.
   */
  void initForces();

  /** Initialization function for the rhs and its jacobian.
   *  \param time the time of initialization
   */
  void initRhs(double time) ;

  /** dynamical system initialization function except for _p:
   *  mainly set memory and compute plug-in for initial state values.
   *  \param time the time of initialization, default value = 0
   *  \param size the size of the memory, default size = 1.
   */
  void initialize(double time = 0, unsigned int size = 1) ;

  /** dynamical system initialization function for _p
   *  \param level for _p
   */
  void initializeNonSmoothInput(unsigned int level) ;

  // === GETTERS AND SETTERS ===

  /** return the dim. of the system (n for first order). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline unsigned int getDim() const
  {
    return _n;
  }
  virtual inline unsigned int getqDim() const
  {
    return _qDim;
  }

  // -- q --

  /** get q
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q() const
  {
    return _q;
  }
  // inline SP::SiconosVector deltaq() const
  // {
  //   return _deltaq;
  // }


  // -- q0 --

  /** get q0
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q0() const
  {
    return _q0;
  }
  inline SP::SiconosVector v0() const
  {
    return _v0;
  }

  // Q memory

  /** get all the values of the state vector q stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory qMemory() const
  {
    return _qMemory;
  }
  inline SP::SiconosMemory vMemory() const
  {
    return _vMemory;
  }

  // -- velocity --

  /** get velocity
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity() const
  {
    return _v;
  }

  // -- velocity0 --

  /** get velocity0
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity0() const
  {
    return _v0;
  }

  // Velocity memory

  /** get all the values of the state vector velocity stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory velocityMemory() const
  {
    return _vMemory;
  }

  // -- p --

  /** get p
   *  \param level unsigned int, required level for p, default = 2
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector p(unsigned int level = 2) const
  {
    return _p[level];
  }



  // -- Mass --

  /** get mass value
   *  \return a double
   */
  inline double massValue() const
  {
    return _mass;
  };


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
  }

  /** set mExt to pointer newPtr
    *  \param newPtr a SP to a Simple vector
    */
  inline void setMExtPtr(SP::SiconosVector newPtr)
  {
    _mExt = newPtr;
  }



  // -- forces --

  /** get the value of forces
   *  \return SiconosVector
   */
  inline const SiconosVector getForces() const
  {
    return *_forces;
  }

  /** get forces
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector forces() const
  {
    return _forces;
  }

  // -- Jacobian Forces w.r.t q --


  /** get JacobianqForces
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SimpleMatrix jacobianqForces() const
  {
    return _jacobianqForces;
  }


  /** get JacobianvForces
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SimpleMatrix jacobianvForces() const
  {
    return _jacobianvForces;
  }
  //  inline SP::SiconosMatrix jacobianZFL() const { return jacobianZFL; }


  inline void setComputeJacobianFIntqByFD(bool value)
  {
    _computeJacobianFIntqByFD=value;
  }
  inline void setComputeJacobianFIntvByFD(bool value)
  {
    _computeJacobianFIntvByFD=value;
  }
  inline void setComputeJacobianMIntqByFD(bool value)
  {
    _computeJacobianMIntqByFD=value;
  }
  inline void setComputeJacobianMIntvByFD(bool value)
  {
    _computeJacobianMIntvByFD=value;
  }




  // --- PLUGINS RELATED FUNCTIONS ---

  /** allow to set a specified function to compute _fExt
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeFExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginFExt->setComputeFunction(pluginPath, functionName);
  }
  /** allow to set a specified function to compute _mExt
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeMExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginMExt->setComputeFunction(pluginPath, functionName);
  }

  /** set a specified function to compute _fExt
   *  \param fct a pointer on the plugin function
   */
  void setComputeFExtFunction(FExt_NE fct)
  {
    _pluginFExt->setComputeFunction((void*)fct);
  }

  /** set a specified function to compute _mExt
   *  \param fct a pointer on the plugin function
   */
  void setComputeMExtFunction(FExt_NE fct)
  {
    _pluginMExt->setComputeFunction((void*)fct);
  }

  /** allow to set a specified function to compute _fInt
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeFIntFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginFInt->setComputeFunction(pluginPath, functionName);
  }
  /** allow to set a specified function to compute _mInt
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputeMIntFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginMInt->setComputeFunction(pluginPath, functionName);
  }

  /** set a specified function to compute _fInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeFIntFunction(FInt_NE fct)
  {
    _pluginFInt->setComputeFunction((void*)fct);
  }

  /** set a specified function to compute _mInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeMExtFunction(FInt_NE fct)
  {
    _pluginMInt->setComputeFunction((void*)fct);
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

  /** default function to compute the external forces
   *  \param time the current time
   */
  virtual void computeFExt(double time);

  /** default function to compute the external moments
   * \param time the current time
   */
  virtual void computeMExt(double time);

  /** default function to compute the internal forces
   *  \param time the current time
   */
  void computeFInt(double time);

  /** default function to compute the internal moments
   * \param time the current time
   */
  void computeMInt(double time);

  /** default function to compute the internal forces
   * \param time the current time
   * \param q
   * \param v
   */
  void computeFInt(double time, SP::SiconosVector q, SP::SiconosVector v);

  /** default function to compute the internal moments
   * \param time the current time
   * \param q
   * \param v
   */
  void computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v);

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
   * \param mInt the computed internal moment vector
   */
  virtual void computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v, SP::SiconosVector mInt);

  /** Default function to compute the right-hand side term
   *  \param time current time
   *  \param isDSup flag to avoid recomputation of operators
   */

  virtual void computeRhs(double time, bool isDSup = false);

  /** Default function to compute jacobian of the right-hand side term according to x
   *  \param time current time
   *  \param isDup flag to avoid recomputation of operators
   */
  virtual void computeJacobianRhsx(double time, bool isDup = false);

  /** Default function to compute forces
   *  \param time double, the current time
   */
  virtual void computeForces(double time);

  /** function to compute forces with some specific values for q and velocity (ie not those of the current state).
   *  \param time double : the current time
   *  \param q SP::SiconosVector: pointers on q
   *  \param velocity SP::SiconosVector: pointers on velocity
   */
  virtual void computeForces(double time,
                             SP::SiconosVector q,
                             SP::SiconosVector velocity);

  /** Default function to compute the jacobian w.r.t. q of forces
   *  \param time double, the current time
   */
  virtual void computeJacobianqForces(double time);

  /** Default function to compute the jacobian w.r.t. v of forces
   *  \param time double, the current time
   */
  virtual void computeJacobianvForces(double time);


  /** function to compute gyroscopic forces with some specific values for q and velocity (ie not those of the current state).
   *  \param velocity SP::SiconosVector: pointers on  velocity vector
   */
  virtual void computeFGyr(SP::SiconosVector velocity);

  /** function to compute gyroscopic forces with some specific values for q and velocity (ie not those of the current state).
   *  \param velocity SP::SiconosVector: pointers on  velocity vector
   *  \param SP::SiconosVector fGyr
   */
  virtual void computeFGyr(SP::SiconosVector velocity, SP::SiconosVector fGyr);


  /** Default function to compute the jacobian following q of fGyr
   *  \param time the current time
   */
  virtual void computeJacobianFGyrv(double time);

  /** Default function to compute the jacobian following q of fGyr
   * by forward finite difference
   *  \param time the current time
   */
  virtual void computeJacobianFGyrvByFD(double time, SP::SiconosVector q, SP::SiconosVector velocity);

  // /** Default function to compute the jacobian following v of fGyr
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
   * \param velocity SP::SiconosVector
   */
  virtual void computeJacobianFIntq(double time,
                                    SP::SiconosVector position,
                                    SP::SiconosVector velocity);
  /** To compute the jacobian w.r.t q of the internal forces
   * by forward finite difference
   * \param time double
   * \param position SP::SiconosVector
   * \param velocity SP::SiconosVector
   */
  void computeJacobianFIntqByFD(double time,
                                SP::SiconosVector position,
                                SP::SiconosVector velocity);


  /** To compute the jacobian w.r.t. v of the internal forces
   *  \param time double: the current time
   * \param position SP::SiconosVector
   * \param velocity SP::SiconosVector
   */
  virtual void computeJacobianFIntv(double time,
                                       SP::SiconosVector position,
                                       SP::SiconosVector velocity);

  /** To compute the jacobian w.r.t v of the internal forces
   * by forward finite difference
   * \param time double
   * \param position SP::SiconosVector
   * \param velocity SP::SiconosVector
   */
  void computeJacobianFIntvByFD(double time,
                                SP::SiconosVector position,
                                SP::SiconosVector velocity);


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
   * \param velocity SP::SiconosVector
   */
  virtual void computeJacobianMIntq(double time,
                                    SP::SiconosVector position,
                                    SP::SiconosVector velocity);

  /** To compute the jacobian w.r.t q of the internal moments
   * by forward finite difference
   * \param time double
   * \param position SP::SiconosVector
   * \param velocity SP::SiconosVector
   */
  void computeJacobianMIntqByFD(double time,
                                SP::SiconosVector position,
                                SP::SiconosVector velocity);




  /** To compute the jacobian w.r.t. v of the internal forces
   *  \param time double: the current time
   * \param position SP::SiconosVector
   * \param velocity SP::SiconosVector
   */
  virtual void computeJacobianMIntv(double time,
                                       SP::SiconosVector position,
                                       SP::SiconosVector velocity);

  /** To compute the jacobian w.r.t v of the internal moments
   * by forward finite difference
   * \param time double
   * \param position SP::SiconosVector
   * \param velocity SP::SiconosVector
   */
  void computeJacobianMIntvByFD(double time,
                                SP::SiconosVector position,
                                SP::SiconosVector velocity);



  // --- miscellaneous ---

  /** print the data to the screen
   */
  void display() const;

  /** initialize the SiconosMemory objects with a positive size.
   * \param steps the size of the SiconosMemory (i)
   */
  void initMemory(unsigned int steps);

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory();

  /** set p[...] to zero
   */
  void resetAllNonSmoothPart();

  /** set p[...] to zero for a given level
   * \param level
   */
  void resetNonSmoothPart(unsigned int level);


  virtual void computeT();

  virtual void computeTdot();



  virtual void normalizeq();

  inline SP::SimpleMatrix mass()
  {
    return _massMatrix;
  }
  inline SP::SimpleMatrix luW()
  {
    return _luW;
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
  inline SP::SiconosMemory forcesMemory()
  {
    return _forcesMemory;
  }
  inline SP::SiconosMemory dotqMemory()
  {
    return _dotqMemory;
  }
  inline SP::SiconosVector dotq()
  {
    return _dotq;
  }

  /** set Boundary Conditions
   *  \param newbd BoundaryConditions
   */
  inline void setBoundaryConditions(SP::BoundaryCondition newbd)
  {
    _boundaryConditions = newbd;
  };

  /** get Boundary Conditions
   *  \return SP::BoundaryCondition pointer on a BoundaryConditions
   */
  inline SP::BoundaryCondition boundaryConditions()
  {
    return _boundaryConditions;
  };

  /** set Reaction to Boundary Conditions
   *  \param newrbd BoundaryConditions pointer
   */
  inline void setReactionToBoundaryConditions(SP::SiconosVector newrbd)
  {
    _reactionToBoundaryConditions = newrbd;
  };

  /** get Reaction to  Boundary Conditions
   *  \return pointer on a BoundaryConditions
   */
  inline SP::SiconosVector reactionToBoundaryConditions()
  {
    return _reactionToBoundaryConditions;
  };


  /** get the matrix converting the object coordinates in the absolute coordinates.
      \return SP::SimpleMatrix
   */
  SP::SimpleMatrix MObjToAbs()
  {
    return _MObjToAbs;
  }
  /*update the _MObjToAbs from the current quaternion.*/
  void computeMObjToAbs();

  // /* update the _MObjToAbs from a given quaternion.
  //  * \param q
  //  */
  // void computeMObjToAbs(SP::SiconosVector q);

  ACCEPT_STD_VISITORS();

};
#endif // NEWTONEULERNLDS_H
