/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
  \brief LagrangianDS class - Second Order Non Linear Dynamical Systems.
*/

#ifndef LAGRANGIANDS_H
#define LAGRANGIANDS_H

#include "DynamicalSystem.hpp"
#include "BoundaryCondition.hpp"
#include "SiconosConst.hpp"

/** Lagrangian non linear dynamical systems - Derived from DynamicalSystem -
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 29, 2004
 *
 * \class LagrangianDS
 * The class LagrangianDS  defines  and computes a generic ndof-dimensional
 * Lagrangian Non Linear Dynamical System of the form :
 * \f[
 * \begin{cases}
 * M(q,z) \dot v + N(v, q, z) + F_{Int}(v , q , t, z) = F_{Ext}(t, z) + p \\
 * \dot q = v
 * \end{cases}
 * \f]
 * where
 * - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 * - \f$ \dot q =v \in R^{ndof} \f$ the velocity, i. e. the time
 *      derivative of the generalized coordinates (Lagrangian systems).
 * - \f$ \ddot q =\dot v \in R^{ndof} \f$ the acceleration, i. e. the second
 *       time derivative of the generalized coordinates.
 * - \f$ p \in R^{ndof} \f$ the reaction forces due to the Non Smooth
 *       Interaction.
 * - \f$ M(q) \in R^{ndof \times ndof} \f$ is the inertia term saved
 *       in the SiconosMatrix mass.
 * - \f$ N(\dot q, q) \in R^{ndof}\f$ is the non linear inertia term
 *       saved in the SiconosVector _fGyr.
 * - \f$ F_{Int}(\dot q , q , t) \in R^{ndof} \f$ are the internal
 *       forces saved in the SiconosVector fInt.
 * - \f$ F_{Ext}(t) \in R^{ndof} \f$ are the external forces saved in
 *       the SiconosVector fExt.
 * - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic
 *       variables, some sort of discrete state.
 *
 * The equation of motion is also shortly denoted as:
 * \f[
 * M(q,z) \dot v = F(v, q, t, z) + p
 * \f]
 *
 * where
 * - \f$F(v, q, t, z) \in R^{ndof} \f$ collects the total forces
 * acting on the system, that is
 * \f[ F(v, q, t, z) =  F_{Ext}(t, z) -  FGyr(v, q, z) + F_{Int}(v, q , t, z) \f]
 * This vector is stored in the  SiconosVector _Forces
 *
 * Links with first order DynamicalSystem top-class are:
 *
 * \f$ n= 2 ndof \f$
 * \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$
 *
 * The rhs is given by:
 * \f[
 * \dot x = \left[\begin{array}{c}
 *  \dot q  \                                           \
 * \ddot q = M^{-1}(q)\left[F(v, q , t, z) + p \right]\\
 * \end{array}\right]
 * \f]
 * Its jacobian is:
 * \f[
 * \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
 *  0  & I \\
 * \nabla_{q}(M^{-1}(q)F(v, q , t, z)) &  \nabla_{\dot q}(M^{-1}(q)F(v, q , t, z)) \\
 * \end{array}\right]
 * \f]
 *  The input due to the non smooth law is:
 * \f[
 * r = \left[\begin{array}{c}0 \\ p \end{array}\right]
 * \f]
 *
 *  Main functionalities to handle a LagrangianDS are:
 *
 * - Construction: the only required operator is M. All the operators
 *      can be set using the plug-in mechanism.
 * - Initialization: compute state members and operators for time=t0
 *     (usually done when calling simulation->initialize)
 * - Computation at time t, thanks to "compute" functions. Any call to
 *      one of the following functions requires that the plug-in
 *      has been set properly thanks to the corresponding setPluginFunction:
 *        => computeMass     (setComputeMassFunction)
 *        => computeFInt     (setComputeFIntFunction)
 *        => computeFExt     (setComputeFExtFunction)
 *        => computeFGyr      (setComputeFGyrFunction)
 *        => computeJacobianFIntq         (setComputeJacobianFIntqFunction)
 *        => computeJacobianFintVelocity  (setComputeJacobianFintVelocityFunction)
 *        => computeJacobianFGyrq          (setComputeJacobianFGyrqFunction)
 *        => computeJacobianFGyrVelocity   (setComputeJacobianFGyrVelocityFunction)
 *        => computeRhs            (no set function)
 *        => computeJacobianRhsx   (no set function)
 *
 * About notation:
 *    - q[i] is the derivative number i of q.
 * Thus: q[0]=\f$ q \f$, global coordinates, q[1]=\f$ \dot q\f$, velocity, q[2]=\f$ \ddot q \f$, acceleration.
 *
 *
 */
class LagrangianDS : public DynamicalSystem
{
public:

  /** List of indices used to save tmp work matrics (last one is the size of the present list) */
  enum WorkNames {xfree, sizeWorkVec};
  enum WorkMatrixNames {invMass, jacobianXBloc10, jacobianXBloc11, zeroMatrix, idMatrix, coeffs_denseoutput, sizeWorkMat};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianDS);
  // -- MEMBERS --

  /** number of degrees of freedom of the system */
  unsigned int _ndof;

  /** state of the system. See details on top of page. */
  VectorOfVectors _q;

  /** initial coordinates of the system */
  SP::SiconosVector _q0;

  /** initial velocity of the system */
  SP::SiconosVector _velocity0;

  /** memory of previous coordinates of the system */
  SP::SiconosMemory _qMemory;

  /** memory of previous velocities of the system */
  SP::SiconosMemory _velocityMemory;

  SP::BlockMatrix _jacxRhs;

  /** "Reaction", generalized forces or imuplses due to the non smooth law
   * The index corresponds to the kinematic
   * level of the corresponding constraints. It mainly depends on what the simulation
   * part want to store, but some rules have to be followed. For instance :
   *  - for the constraints at the acceleration level, _p[2] stores the reaction forces,
   *  - for the constraints at the veocity level,  _p[1] stores the (discrete) reaction impulse
   *  - for the constraints at the position level, _p[0] stores the multiplier for a constraint
   * in position
   */
  std::vector<SP::SiconosVector> _p;

  /** memory of previous generalized forces due to constraints */
  VectorOfMemories _pMemory;


  /** mass of the system */
  SP::SiconosMatrix _mass;

  /** internal forces applied to  the system */
  SP::SiconosVector _fInt;

  /** jacobian_q FInt*/
  SP::SiconosMatrix _jacobianFIntq;

  /** jacobian_{qDot} FInt*/
  SP::SiconosMatrix _jacobianFIntqDot;

  /** external forces applied to the system */
  SP::SiconosVector _fExt;

  /** non-linear inertia term of the system */
  SP::SiconosVector _fGyr;

  /** jacobian_q FGyrq*/
  SP::SiconosMatrix _jacobianFGyrq;
  /** jacobian_{qDot} FGyrq*/
  SP::SiconosMatrix _jacobianFGyrqDot;

  /** forces(q[0],q[1],t)= fExt - fInt -FGyr */
  SP::SiconosVector _forces;

  /** jacobian_q forces*/
  SP::SiconosMatrix _jacobianqForces;

  /** jacobian_{qDot} forces*/
  SP::SiconosMatrix _jacobianqDotForces;

  /** memory of previous forces of the system */
  SP::SiconosMemory _forcesMemory;

  /** Boundary condition applied to a dynamical system*/
  SP::BoundaryCondition _boundaryConditions;

  /** Reaction to an applied  boundary condition */
  SP::SiconosVector _reactionToBoundaryConditions;

  /** set links with DS members
   */
  void connectToDS();

  /** Default constructor
   */
  LagrangianDS();


  // pointers to functions member to compute plug-in functions

  /** LagrangianDS plug-in to compute mass(q,t) - id = "mass"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param[in,out] mass : pointer to the first element of mass
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  //  void (*computeMassPtr)(unsigned int, double*, double*, unsigned int, double*);
  SP::PluggedObject _pluginMass;


  /** LagrangianDS plug-in to compute internal forces \f$F_{int}(t,q,\dot q)\f$ - id = "fInt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] fInt : pointer to the first element of fInt
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginFInt;
  //  FPtr6 computeFIntPtr;

  /** LagrangianDS plug-in to compute external forces \f$F_{Ext}(t)\f$, id = "fExt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param[in,out] fExt : pointer to the first element of fExt
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  //  void (*computeFExtPtr)(double, unsigned int, double*, unsigned int, double* );
  SP::PluggedObject _pluginFExt;

  /** LagrangianDS plug-in to compute \f$FGyr(\dot q, q)\f$, id = "FGyr"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] FGyr : pointer to the first element of FGyr
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr5 computeFGyrPtr;
  SP::PluggedObject _pluginFGyr;

  /** LagrangianDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianFIntq"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr6 computeJacobianFIntqPtr;
  SP::PluggedObject _pluginJacqFInt;

  /** LagrangianDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianFIntqDot"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr6 computeJacobianFIntqDotPtr;
  SP::PluggedObject _pluginJacqDotFInt;

  /** LagrangianDS plug-in to compute \f$\nabla_qFGyr(\dot q, q)\f$, id = "jacobianFGyrq"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr5 computeJacobianFGyrqPtr;
  SP::PluggedObject _pluginJacqFGyr;
  /** LagrangianDS plug-in to compute \f$\nabla_{\dot q}FGyr(\dot q, q)\f$, id = "jacobianFGyrqDot"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr5 computeJacobianFGyrqDotPtr;
  SP::PluggedObject _pluginJacqDotFGyr;

  virtual void zeroPlugin();
public:

  // === CONSTRUCTORS - DESTRUCTOR ===

  /** constructor from a minimum set of data
   *  \param position SiconosVector : initial coordinates of this DynamicalSystem
   *  \param velocity SiconosVector : initial velocity of this DynamicalSystem
   */
  LagrangianDS(SP::SiconosVector position, SP::SiconosVector velocity);

  /** constructor from a minimum set of data
   *  \param position SiconosVector : initial coordinates of this DynamicalSystem
   *  \param velocity SiconosVector : initial velocity of this DynamicalSystem
   *  \param mass SiconosMatrix : mass matrix
   */
  LagrangianDS(SP::SiconosVector position,
               SP::SiconosVector velocity, SP::SiconosMatrix mass);

  /** constructor from a minimum set of data
   *  \param position SiconosVector : initial coordinates of this DynamicalSystem
   *  \param velocity SiconosVector : initial velocity of this DynamicalSystem
   *  \param plugin std::string: plugin path to compute mass matrix
   */
  LagrangianDS(SP::SiconosVector position, SP::SiconosVector velocity, const std::string& plugin);

  /** destructor */
  virtual ~LagrangianDS();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

  /** allocate memory for forces and its jacobians, if required.
   */
  void initForces();

  /** Initialization function for the rhs and its jacobian.
   *  \param time of initialization
   */
  void initRhs(double time) ;

  /** dynamical system initialization function except for _p:
   *  mainly set memory and compute plug-in for initial state values.
   *  \param time of initialisation, default value = 0
   *  \param size the size of the memory, default size = 1.
   */
  void initialize(double time = 0, unsigned int size = 1) ;

  /** dynamical system initialization function for _p
   *  \param level for _p
   */
  void initializeNonSmoothInput(unsigned int level) ;

  // === GETTERS AND SETTERS ===

  /** to get the value of ndof
   *  \return the value of ndof
   */
  inline unsigned int getNdof() const
  {
    return _ndof;
  };

  /** to set ndof
   *  \param newNdof unsigned int : the value to set ndof
   */
  inline void setNdof(unsigned int newNdof)
  {
    _ndof = newNdof;
  };

  /** return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline unsigned int dimension() const
  {
    return _ndof;
  }

  // -- q --

  /** get q
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q() const
  {
    return _q[0];
  }

  /** set the value of q to newValue
   *  \param newValue
   */
  void setQ(const SiconosVector& newValue);

  /** set Q to pointer newPtr
   *  \param newPtr
   */
  void setQPtr(SP::SiconosVector newPtr);

  // -- q0 --

  /** get q0
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q0() const
  {
    return _q0;
  }

  /** set the value of q0 to newValue
   *  \param newValue
   */
  void setQ0(const SiconosVector& newValue);

  /** set Q0 to pointer newPtr
   *  \param newPtr
   */
  void setQ0Ptr(SP::SiconosVector newPtr);

  // Q memory

  /** get all the values of the state vector q stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory qMemory() const
  {
    return _qMemory;
  }

  // -- velocity --

  /** get velocity
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity() const
  {
    return _q[1];
  }

  /** set the value of velocity to newValue
   *  \param newValue
   */
  void setVelocity(const SiconosVector& newValue);

  /** set Velocity to pointer newPtr
   *  \param newPtr
   */
  void setVelocityPtr(SP::SiconosVector newPtr);

  // -- velocity0 --

  /** get velocity0
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity0() const
  {
    return _velocity0;
  }

  /** set the value of velocity0 to newValue
   *  \param newValue
   */
  void setVelocity0(const SiconosVector& newValue);

  /** set Velocity0 to pointer newPtr
   *  \param newPtr
   */
  void setVelocity0Ptr(SP::SiconosVector newPtr) ;

  // -- acceleration --

  /** get acceleration
   *  \return pointer on a SiconosVector
   */
  SP::SiconosVector acceleration() const ;

  // Velocity memory

  /** get all the values of the state vector velocity stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory velocityMemory() const
  {
    return _velocityMemory;
  }

  // -- p --

  /** get p
   *  \param level unsigned int, required level for p
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector p(unsigned int level) const
  {
    return _p[level];
  }

  /** set the value of p to newValue
   *  \param newValue
   *  \param level required level for p
   */
  void setP(const SiconosVector& newValue, unsigned int level);

  /** set P to pointer newPtr
   *  \param newPtr
   *  \param level required level for p, default
   */
  void setPPtr(SP::SiconosVector newPtr, unsigned int level);

  /** get all the values of the state vector p stored in memory
   * \param level
   *  \return a memory
   */
  inline SP::SiconosMemory pMemory(unsigned int level) const
  {
    return _pMemory[level];
  }

  // -- Mass --

  /** get mass
   *  \return pointer on a plugged-matrix
   */
  inline SP::SiconosMatrix mass() const
  {
    return _mass;
  }

  /** set mass to pointer newPtr
   *  \param newPtr a plugged matrix SP
   */
  inline void setMassPtr(SP::SiconosMatrix newPtr)
  {
    _mass = newPtr;
  }

  /** get MassLU: a copy of the mass matrix which is LU-factorized. Temporary function?
   *  \return a pointer on a SiconosMatrix
   */
  inline SP::SimpleMatrix massLU() const
  {
    return (_workMatrix[invMass]);
  }

  // --- fInt ---

  /** get fInt
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector fInt() const
  {
    return _fInt;
  }


  /** set fInt to pointer newPtr
   *  \param newPtr a SP to plugged vector
   */
  inline void setFIntPtr(SP::SiconosVector newPtr)
  {
    _fInt = newPtr;
  }

  // -- Fext --

  /** get fExt
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector fExt() const
  {
    return _fExt;
  }


  /** set fExt to pointer newPtr
   *  \param newPtr a SP to a Simple vector
   */
  inline void setFExtPtr(SP::SiconosVector newPtr)
  {
    _fExt = newPtr;
  }

  // -- FGyr --

  /** get FGyr
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector FGyr() const
  {
    return _fGyr;
  }


  /** set FGyr to pointer newPtr
   *  \param newPtr a SP to plugged vector
   */
  inline void setFGyrPtr(SP::SiconosVector newPtr)
  {
    _fGyr = newPtr;
  }


  // -- Jacobian FInt --

  /** get jacobianFIntq
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianFIntq() const
  {
    return _jacobianFIntq;
  }
  /** get jacobianFIntqDot
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianFIntqDot() const
  {
    return _jacobianFIntqDot;
  }
  //  inline SP::SiconosMatrix jacobianZFInt() const { return jacobianZFInt; }


  /** set jacobianFIntq to pointer newPtr
   *  \param newPtr a SP SiconosMatrix
   */
  inline void setJacobianFIntqPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianFIntq = newPtr;
  }
  /** set jacobianFIntqDot to pointer newPtr
   *  \param newPtr a SP SiconosMatrix
   */
  inline void setJacobianFIntqDotPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianFIntqDot = newPtr;
  }

  // -- Jacobian FGyr --


  /** get jacobianFGyrq
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianFGyrq() const
  {
    return _jacobianFGyrq;
  }
  /** get jacobianFGyrqDot
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianFGyrqDot() const
  {
    return _jacobianFGyrqDot;
  }


  /** set jacobianFGyrq to pointer newPtr
   *  \param newPtr a SP SiconosMatrix
   */
  inline void setJacobianFGyrqPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianFGyrq = newPtr;
  }
  /** set jacobianFGyrqDot to pointer newPtr
   *  \param newPtr a SP SiconosMatrix
   */
  inline void setJacobianFGyrqDotPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianFGyrqDot = newPtr;
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
  /** get forces
   *  \return pointer on a SiconosMemory
   */
  inline SP::SiconosMemory forcesMemory()
  {
    return _forcesMemory;
  }
  // -- Jacobian forces --


  /** get JacobianForces
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianqForces() const
  {
    return _jacobianqForces;
  }

  /** get JacobianqDotForces
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianqDotForces() const
  {
    return _jacobianqDotForces;
  }
  //  inline SP::SiconosMatrix jacobianZFL() const { return jacobianZFL; }

  // --- PLUGINS RELATED FUNCTIONS ---

  /** allow to set a specified function to compute the mass
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeMassFunction(const std::string&  pluginPath, const std::string&  functionName)
  {
    _pluginMass->setComputeFunction(pluginPath, functionName);
    if (!_mass) _mass.reset(new SimpleMatrix(_ndof, _ndof));
  }

  /** set a specified function to compute Mass
   *  \param fct a pointer on the plugin function
   */
  void setComputeMassFunction(FPtr7 fct)
  {
    _pluginMass->setComputeFunction((void*)fct);
  }

  /** allow to set a specified function to compute FInt
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeFIntFunction(const std::string&  pluginPath, const std::string&  functionName)
  {
    _pluginFInt->setComputeFunction(pluginPath, functionName);
    if (!_fInt) _fInt.reset(new SiconosVector(_ndof));
    //    Plugin::setFunction(&computeFIntPtr, pluginPath,functionName);
  }

  /** set a specified function to compute fInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeFIntFunction(FPtr6 fct)
  {
    _pluginFInt->setComputeFunction((void*)fct);
    //    computeFIntPtr = fct;
  }

  /** allow to set a specified function to compute Fext
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeFExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginFExt->setComputeFunction(pluginPath, functionName);
    if (!_fExt) _fExt.reset(new SiconosVector(_ndof));
    //    Plugin::setFunction(&computeFExtPtr, pluginPath,functionName);
  }

  /** set a specified function to compute fExt
   *  \param fct a pointer on the plugin function
   */
  void setComputeFExtFunction(VectorFunctionOfTime fct)
  {
    _pluginFExt->setComputeFunction((void*)fct);
    //   computeFExtPtr = fct ;
  }

  /** allow to set a specified function to compute the inertia
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeFGyrFunction(const std::string& pluginPath, const std::string&  functionName);

  /** set a specified function to compute FGyr
   *  \param fct a pointer on the plugin function
   */
  void setComputeFGyrFunction(FPtr5 fct);

  /** allow to set a specified function to compute the jacobian w.r.t q of the internal forces
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFIntqFunction(const std::string&  pluginPath, const std::string&  functionName);
  /** allow to set a specified function to compute the jacobian following qDot of the internal forces w.r.t.
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFIntqDotFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute jacobian following q of the FInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFIntqFunction(FPtr6 fct);
  /** set a specified function to compute jacobian following qDot of the FInt
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFIntqDotFunction(FPtr6 fct);

  /** allow to set a specified function to compute the jacobian w.r.t q of the the external forces
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFGyrqFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** allow to set a specified function to compute the jacobian w.r.t qDot of the the external strength
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFGyrqDotFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute the jacobian following q of FGyr
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFGyrqFunction(FPtr5 fct);
  /** set a specified function to compute the jacobian following qDot of FGyr
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFGyrqDotFunction(FPtr5 fct);

  /** get the value of the gradient according to \f$ x \f$ of the right-hand side
   *  \return BlockMatrix&
   */
  inline const BlockMatrix& getJacobianRhsx() const
  {
    return *_jacxRhs;
  }

  /** get gradient according to \f$ x \f$ of the right-hand side (pointer)
   *  \return pointer on a SiconosMatrix
   */
  inline SP::BlockMatrix jacobianRhsx() const
  {
    return _jacxRhs;
  }

  /** default function to compute the mass
   */
  virtual void computeMass();

  /** function to compute the mass
   *  \param position value used to evaluate the mass matrix
   */
  virtual void computeMass(SP::SiconosVector position);

  /** default function to compute the internal strengths
   *  \param time the current time
   */
  virtual void computeFInt(double time);

  /** function to compute the internal strengths
   *  with some specific values for position and velocity (ie not those of the current state).
   *  \param time the current time,
   *  \param position value used to evaluate the internal forces
   *  \param velocity value used to evaluate the internal forces
   */
  virtual void computeFInt(double time,
                           SP::SiconosVector position,
                           SP::SiconosVector velocity);

  /** default function to compute the external strengths
   *  \param time the current time
   */
  virtual void computeFExt(double time);

  /** default function to compute the inertia
   */
  virtual void computeFGyr();

  /** function to compute the inertia
   *  with some specific values for q and velocity (ie not those of the current state).
   *  \param position value used to evaluate the inertia forces
   *  \param velocity value used to evaluate the inertia forces
   */
  virtual void computeFGyr(SP::SiconosVector position,
                          SP::SiconosVector velocity);

  /** To compute the jacobian w.r.t q of the internal forces
   *  \param time the current time
   */
  virtual void computeJacobianFIntq(double time);
  /** To compute the jacobian w.r.t qDot of the internal forces
   *  \param time the current time
   */
  virtual void computeJacobianFIntqDot(double time);

  /** To compute the jacobian w.r.t q of the internal forces
   *  \param time the current time
   *  \param position value used to evaluate the jacobian
   *  \param velocity value used to evaluate the jacobian
   */
  virtual void computeJacobianFIntq(double time,
                                    SP::SiconosVector position,
                                    SP::SiconosVector velocity);

  /** To compute the jacobian w.r.t. qDot of the internal forces
   *  \param time the current time
   *  \param position value used to evaluate the jacobian
   *  \param velocity value used to evaluate the jacobian
   */
  virtual void computeJacobianFIntqDot(double time,
                                       SP::SiconosVector position,
                                       SP::SiconosVector velocity);

  /** function to compute the jacobian w.r.t. q of the inertia forces
   */
  virtual void computeJacobianFGyrq();

  /** function to compute the jacobian w.r.t. qDot of the inertia forces
   */
  virtual void computeJacobianFGyrqDot();

  /** function to compute the jacobian w.r.t. q of the inertia forces
   *  \param position value used to evaluate the jacobian
   *  \param velocity value used to evaluate the jacobian
   */
  virtual void computeJacobianFGyrq(SP::SiconosVector position, SP::SiconosVector velocity);

  /** function to compute the jacobian w.r.t. qDot of the inertia forces
   *  \param position value used to evaluate the jacobian
   *  \param velocity value used to evaluate the jacobian
   */
  virtual void computeJacobianFGyrqDot(SP::SiconosVector position, SP::SiconosVector velocity);

  /** Default function to compute the right-hand side term
   *  \param time current time
   *  \param isDSup flag to avoid recomputation of operators
   */
  virtual void computeRhs(double time, bool isDSup = false);

  /** Default function to compute jacobian of the right-hand side term according to x
   *  \param time the current time
   *  \param isDSup flag to avoid recomputation of operators
   */
  virtual void computeJacobianRhsx(double time, bool isDSup = false);

  /** Default function to compute forces
   *  \param time the current time
   */
  virtual void computeForces(double time);

  /** function to compute forces with some specific values for q and velocity (ie not those of the current state).
   *  \param time the current time
   *  \param q SP::SiconosVector: pointers on q
   *  \param velocity SP::SiconosVector: pointers on velocity
   */
  virtual void computeForces(double time,
                             SP::SiconosVector q,
                             SP::SiconosVector velocity);

  /** Default function to compute the jacobian w.r.t. q of forces
   *  \param time the current time
   */
  virtual void computeJacobianqForces(double time);

  /** Default function to compute the jacobian w.r.t. qDot of forces
   *  \param time the current time
   */
  virtual void computeJacobianqDotForces(double time);


  /** To compute the kinetic energy
   */
  double computeKineticEnergy();

  // --- miscellaneous ---

  /** print the data to the screen
   */
  void display() const;

  /** initialize the SiconosMemory objects with a positive size.
   *  \param size the size of the SiconosMemory. must be >= 0
   */
  void initMemory(unsigned int size);

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory();

  /** To compute \f$\frac{|q_{i+1} - qi|}{|q_i|}\f$ where \f$ q_{i+1}\f$ represents the present state and \f$ q_i\f$ the previous one
   * \return a double
   */
  /*  double dsConvergenceIndicator(); */

  /** set p[...] to zero
   */
  void resetAllNonSmoothPart();

  /** set p[...] to zero for a given level
      \param level
   */
  void resetNonSmoothPart(unsigned int level);

  /** Computes post-impact velocity, using pre-impact velocity and impulse (p) value.
   * Used in EventDriven (LsodarOSI->updateState)
   */
  void computePostImpactVelocity();

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


  /** to allocate memory for a new tmp matrix
   *  \param id the id of the SimpleMatrix
   *  \param sizeOfRows an int to set the size of rows
   *  \param sizeOfCols an int to set the size of cols
   */
  inline void allocateWorkMatrix(const WorkMatrixNames & id, int sizeOfRows, int sizeOfCols)
  {
    _workMatrix[id].reset(new SimpleMatrix(sizeOfRows, sizeOfCols));
  }

  /** get a temporary saved matrix
   * \param id the id of the SimpleMatrix
   *  \return a SP::SimpleMatrix
   */
  inline SP::SimpleMatrix getWorkMatrix(const WorkMatrixNames & id) const
  {
    return  _workMatrix[id];
  }

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianDS)

#endif // LAGRANGIANDS_H
