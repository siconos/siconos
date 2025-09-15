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

/*! \file LagrangianDS.hpp
  LagrangianDS class - Second Order Non Linear Dynamical Systems.
*/

#ifndef LAGRANGIANDS_H
#define LAGRANGIANDS_H

#include "BoundaryCondition.hpp"
#include "DynamicalSystem.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosConst.hpp"

/**   
   Lagrangian non linear dynamical systems -  \f$ M(q,z) \dot v = F(v, q, t, z) +
   p \f$

   This class defines and computes a generic ndof-dimensional
   Lagrangian Non Linear Dynamical System of the form :
   
   \f[
     M(q,z) \dot v + F_{gyr}(v, q, z) + F_{int}(v , q , t, z) = F_{ext}(t, z) + p \\
     \dot q = v
   \f]

   where
   
   -  \f$ q \in R^{ndof} \f$ is the set of the generalized coordinates,
   - \f$ \dot q =v \in R^{ndof} \f$ the velocity, i. e. the time
   derivative of the generalized coordinates (Lagrangian systems).
   - \f$ \ddot q =\\dot v \in R^{ndof} \f$ the acceleration, i. e. the second
   time derivative of the generalized coordinates.
   - \f$ p \in R^{ndof} \f$ the reaction forces due to the Non Smooth
   Interaction.
   - \f$ M(q) \in R^{ndof \times ndof} \f$ is the inertia term (access : mass()
   method).
   - \f$ F_{gyr}(\dot q, q) \in R^{ndof} \f$ is the non linear inertia term
   (access fGyr() method).
   - \f$ F_{int}(\dot q , q , t) \in R^{ndof} \f$ are the internal
   forces (access fInt() method).
   - \f$ F_{ext}(t) \in R^{ndof} \f$ are the external forces (access fExt()
   method).
   - \f$ z \in R^{zSize} \f$  is a vector of arbitrary algebraic
   variables, some sort of discrete state.
   
   The equation of motion is also shortly denoted as  \f$ M(q,z) \dot v = F(v, q,
   t, z) + p \f$ 
   
   where  \f$ F(v, q, t, z) \in R^{ndof} \f$ collects the total forces acting on the
   system, that is \f$ F(v, q, t, z) =  F_{ext}(t, z) -  F_{gyr}(v, q, z) + F_{int}(v, q , t, z) \f$
   
   This vector is saved and may be accessed using forces() method.
   
   q[i] is the derivative number i of q.
   Thus: q[0]=\f$ q \f$, global coordinates, q[1]= \f$  \dot q \f$ , velocity,
   q[2]= \f$  \ddot q \f$, acceleration.
   
   The following operators (and their jacobians) can be plugged, in the usual way
   (see User Guide, 'User-defined plugins')
   
   -  \f$ M(q) \f$  (computeMass())
   -  \f$ F_{gyr}(v, q, z) \f$  (computeFGyr())
   -  \f$ F_{int}(v , q , t, z) \f$  (computeFInt())
   -  \f$ F_{ext}(t, z) \f$  (computeFExt())
   
   If required (e.g. for Event-Driven like simulation), formulation as a
   first-order system is also available, and writes:
   
   - \f$ n= 2 ndof \f$
   - \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right] \f$ 
   - rhs given by:
   
   \f[

      \dot x = \left[\begin{array}{c}
      \dot q\\
      \ddot q = M^{-1}(q)\left[F(v, q , t, z) + p \right]\\
      \end{array}\right]

   \f]   \endrst

   - jacobian of the rhs, with respect to x

   \f[
       \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
       0  & I \\
       \nabla_{q}(M^{-1}(q)F(v, q , t, z)) &  \nabla_{\dot q}(M^{-1}(q)F(v, q ,
       t, z)) \\ \end{array}\right]
   \f]
   
   with the input due to the non smooth law:

   \f[

      \left[\begin{array}{c}
      0 \\
      p \end{array}\right]
   \f]

   In that case, use the following methods:
   - initRhs() to allocate/initialize memory for these new operators,
   - rhs() to get the rhs vector
   - computeRhs(), computeJacobianRhsx() ..., to update the content of rhs, its
   jacobians ...
   
*/
class LagrangianDS : public SecondOrderDS {

protected:
  
  ACCEPT_SERIALIZATION(LagrangianDS);

  /**
     Common code for constructors
     should be replaced in C++11 by delegating constructors
     
     \param position vector of initial positions
     \param velocity vector of initial velocities
  */
  void _init(SP::SiconosVector position, SP::SiconosVector velocity);

  // -- MEMBERS --

  /** state of the system. See details on top of page. */
  VectorOfVectors _q;

  /** initial velocity of the system */
  SP::SiconosVector _velocity0;

  /** memory of previous coordinates of the system */
  SiconosMemory _qMemory;

  /** memory of previous velocities of the system */
  SiconosMemory _velocityMemory;

  /** internal forces applied to  the system */
  SP::SiconosVector _fInt;

  // Should we use this enum to clarify notations in LagrangianDS?
  // enum LagrangianDSJacobianId {Jacobian_FInt_wrt_q, Jacobian_FInt_wrt_qDot,
  //                              Jacobian_FGyr_wrt_q, Jacobian_FGyr_wrt_qdot,
  //                              Jacobian_Force_wrt_q,
  //                              Jacobian_Forces_wrt_qDot, numberOfJacobians};

  /** jacobian_q FInt*/
  SP::SiconosMatrix _jacobianFIntq;

  /** jacobian_{qDot} FInt*/
  SP::SiconosMatrix _jacobianFIntqDot;

  /** external forces applied to the system */
  SP::SiconosVector _fExt;

  /** boolean if _fext is constant (set thanks to setFExtPtr for instance)
   * false by default */
  bool _hasConstantFExt = false;

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
  SiconosMemory _forcesMemory;

  enum LagrangianDSRhsMatrices {
    jacobianXBloc10,
    jacobianXBloc11,
    zeroMatrix,
    idMatrix,
    numberOfRhsMatrices
  };
  /** A container of matrices to save matrices that are involed in first order
   * from of LagrangianDS system values (jacobianXBloc10, jacobianXBloc11,
   * zeroMatrix, idMatrix) No get-set functions at the time. Only used as a
   * protected member.*/
  VectorOfSMatrices _rhsMatrices;

  // pointers to functions member to compute plug-in functions

  /** LagrangianDS plug-in to compute mass(q,t) - id = "mass"
   *
   *  \param sizeOfq : size of vector q
   *  \param q : pointer to the first element of q
   *  \param[in,out] mass : pointer to the first element of mass
   *  \param  size of vector z
   *  \param[in,out] z : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginMass;

  /** LagrangianDS plug-in to compute internal forces  \f$ F_{int}(t,q,\dot q) \f$ 
   *  - id = "fInt"
   *
   *  \param time : current time
   *  \param sizeOfq : size of vector q
   *  \param q : pointer to the first element of q
   *  \param velocity : pointer to the first element of velocity
   *  \param[in,out] fInt : pointer to the first element of fInt
   *  \param  size of vector z
   *  \param[in,out] z : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginFInt;
  //  FPtr6 computeFIntPtr;

  /** LagrangianDS plug-in to compute external forces  \f$ F_{Ext}(t) \f$ , id =
   *  "fExt"
   *
   *  \param time : current time
   *  \param sizeOfq : size of vector q
   *  \param[in,out] fExt : pointer to the first element of fExt
   *  \param  size of vector z
   *  \param[in,out] z : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginFExt;

  /** LagrangianDS plug-in to compute  \f$ FGyr(\dot q, q) \f$ , id = "FGyr"
   *  \param sizeOfq : size of vector q
   *  \param q : pointer to the first element of q
   *  \param velocity : pointer to the first element of velocity
   *  \param[in,out] FGyr : pointer to the first element of FGyr
   *  \param  size of vector z
   *  \param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginFGyr;

  /** LagrangianDS plug-in to compute  \f$ \nabla_qF_{Int}(\dot q, q, t) \f$ , id =
   *  "jacobianFIntq"
   *
   *  \param time : current time
   *  \param sizeOfq : size of vector q
   *  \param q : pointer to the first element of q
   *  \param velocity : pointer to the first element of velocity
   *  \param[in,out] jacob : pointer to the first element of the jacobian
   *  \param  size of vector z
   *  \param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacqFInt;

  /** LagrangianDS plug-in to compute  \f$ \nabla_{\dot q}F_{Int}(\dot q, q,
   *  t) \f$ , id = "jacobianFIntqDot"
   *
   *  \param time : current time
   *  \param sizeOfq : size of vector q
   *  \param q : pointer to the first element of q
   *  \param velocity : pointer to the first element of velocity
   *  \param[in,out] jacob : pointer to the first element of the jacobian
   *  \param  size of vector z
   *  \param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacqDotFInt;

  /** LagrangianDS plug-in to compute  \f$ \nabla_qFGyr(\dot q, q) \f$ , id =
   *  "jacobianFGyrq"
   *
   *  \param sizeOfq : size of vector q
   *  \param q : pointer to the first element of q
   *  \param velocity : pointer to the first element of velocity
   *  \param[in,out] jacob : pointer to the first element of the jacobian
   *  \param  size of vector z
   *  \param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacqFGyr;

  /** LagrangianDS plug-in to compute  \f$ \nabla_{\dot q}FGyr(\dot q, q) \f$ , id =
   *  "jacobianFGyrqDot"
   *
   *  \param sizeOfq : size of vector q
   *  \param q : pointer to the first element of q
   *  \param velocity : pointer to the first element of velocity
   *  \param[in,out] jacob : pointer to the first element of the jacobian
   *  \param  size of vector z
   *  \param[in,out] z  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginJacqDotFGyr;

  /** build all _plugin... PluggedObject */
  void _zeroPlugin() override;

  /** Default constructor */
  LagrangianDS() : SecondOrderDS(){};

public:
  /** constructor from initial state only,  \f$ dv = p  \f$ 
   *
   *  \param position SiconosVector : initial coordinates of this
   *  DynamicalSystem \param velocity SiconosVector : initial velocity of this
   *  DynamicalSystem
   */
  LagrangianDS(SP::SiconosVector position, SP::SiconosVector velocity);

  /** constructor from initial state and mass,  \f$ Mdv = p \f$
   *
   *  \param position SiconosVector : initial coordinates of this
   *  DynamicalSystem \param velocity SiconosVector : initial velocity of this
   *  DynamicalSystem \param mass SiconosMatrix : mass matrix
   */
  LagrangianDS(SP::SiconosVector position, SP::SiconosVector velocity,
               SP::SiconosMatrix mass);

  /** constructor from initial state and mass (plugin)  \f$ Mdv = p \f$
   *
   *  \param position SiconosVector : initial coordinates of this
   *  DynamicalSystem \param velocity SiconosVector : initial velocity of this
   *  DynamicalSystem \param plugin std::string: plugin path to compute mass
   *  matrix
   */
  LagrangianDS(SP::SiconosVector position, SP::SiconosVector velocity,
               const std::string &plugin);

  /** destructor */
  virtual ~LagrangianDS(){};

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
  void computeJacobianRhsx(double time) override;

  /** reset non-smooth part of the rhs (i.e. p), for all 'levels' */
  void resetAllNonSmoothParts() override;

  /** set nonsmooth part of the rhs (i.e. p) to zero for a given level
   *
   *  \param level
   */
  void resetNonSmoothPart(unsigned int level) override;

  /** set the value of the right-hand side, \f$ \dot x \f$
   *
   *  \param newValue SiconosVector
   */
  void setRhs(const SiconosVector &newValue) override
  {
    THROW_EXCEPTION(
        "LagrangianDS - setRhs call is forbidden for 2nd order systems.");
  }

  /** set right-hand side, \f$ \dot x \f$ (pointer link)
   *
   *  \param newPtr SP::SiconosVector
   */
  void setRhsPtr(SP::SiconosVector newPtr) override
  {
    THROW_EXCEPTION(
        "LagrangianDS - setRhsPtr call is forbidden for 2nd order systems.");
  }

  /* function to compute  \f$ F(v,q,t,z) \f$  for the current state
   *
   *  \param time the current time
   */
  // virtual void computeForces(double time);

  /** Compute  \f$ F(v,q,t,z) \f$ 
   *
   *  \param time the current time
   *  \param q SP::SiconosVector: pointers on q
   *  \param velocity SP::SiconosVector: pointers on velocity
   */
  void computeForces(double time, SP::SiconosVector q,
                     SP::SiconosVector velocity) override;

  /** Compute  \f$ \nabla_qF(v,q,t,z) \f$  for current  \f$ q,v \f$ 
   *  Default function to compute forces
   *
   *  \param time the current time
   */
  void computeJacobianqForces(double time) override;

  /** Compute  \f$ \nabla_{\dot q}F(v,q,t,z) \f$  for current  \f$ q,v \f$ 
   *
   *  \param time the current time
   */
  virtual void computeJacobianqDotForces(double time)
  {
    computeJacobianvForces(time);
  };

  /** Compute  \f$ \nabla_{\dot q}F(v,q,t,z) \f$  for current  \f$ q,v \f$ 
   *
   *  \param time the current time
   */
  void computeJacobianvForces(double time) override;

  /** generalized coordinates of the system (vector of size dimension())
   *
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q() const override { return _q[0]; }

  /** set value of generalized coordinates vector (copy)
   *
   *  \param newValue
   */
  void setQ(const SiconosVector &newValue) override;

  /** set value of generalized coordinates vector (pointer link)
   *
   *  \param newPtr
   */
  void setQPtr(SP::SiconosVector newPtr) override;

  /** set initial state (copy)
   *
   *  \param newValue
   */
  void setQ0(const SiconosVector &newValue) override;

  /** set initial state (pointer link)
   *
   *  \param newPtr
   */
  void setQ0Ptr(SP::SiconosVector newPtr) override;

  /** get velocity vector (pointer link)
   *
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity() const override { return _q[1]; }

  /** set velocity vector (copy)
   *
   *  \param newValue
   */
  void setVelocity(const SiconosVector &newValue) override;

  /** set velocity vector (pointer link)
   *
   *  \param newPtr
   */
  void setVelocityPtr(SP::SiconosVector newPtr) override;

  /** get initial velocity (pointer)
   *
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity0() const override { return _velocity0; }

  /** set initial velocity (copy)
   *
   *  \param newValue
   */
  void setVelocity0(const SiconosVector &newValue) override;

  /** set initial velocity (pointer link)
   *
   *  \param newPtr
   */
  void setVelocity0Ptr(SP::SiconosVector newPtr) override;

  /** get acceleration (pointer link)
   *
   *  \return pointer on a SiconosVector
   */
  SP::SiconosVector acceleration() const override { return _q[2]; };

  /** get \$F_{int}\$ (pointer link)
   *
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector fInt() const { return _fInt; }

  /** set  \$F_{int}\$ (pointer link)
   *
   *  \param newPtr a SP to plugged vector
   */
  inline void setFIntPtr(SP::SiconosVector newPtr) { _fInt = newPtr; }

  /** get  \f$ F_{ext} \f$ , (pointer link)
   *
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector fExt() const { return _fExt; }

  /** set  \f$ F_{ext} \f$ , (pointer link)
   *
   *  \param newPtr a SP to a Simple vector
   */
  inline void setFExtPtr(SP::SiconosVector newPtr)
  {
    _fExt = newPtr;
    _hasConstantFExt = true;
  }

  /** get  \f$ F_{gyr} \f$ , (pointer link)
   *
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector fGyr() const { return _fGyr; }

  /** set  \f$ F_{gyr} \f$ , (pointer link)
   *
   *  \param newPtr a SP to plugged vector
   */
  inline void setFGyrPtr(SP::SiconosVector newPtr) { _fGyr = newPtr; }

  /** get  \f$ \nabla_qF_{int} \f$ , (pointer link)
   *
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianFIntq() const { return _jacobianFIntq; }

  /** get  \f$ \nabla_{\dot q}F_{int} \f$ , (pointer link)
   *
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianFIntqDot() const
  {
    return _jacobianFIntqDot;
  }
  void allocateMass();
  void allocateJacobianFIntq();
  void allocateJacobianFIntqDot();
  void allocateFInt();
  void allocateFExt();

  /** set  \f$ \nabla_{q}F_{int} \f$ , (pointer link)
   *
   *  \param newPtr a pointer to a SiconosMatrix
   */
  inline void setJacobianFIntqPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianFIntq = newPtr;
  }

  /** set  \f$ \nabla_{\dot q}F_{int} \f$ , (pointer link)
   *
   *  \param newPtr a pointer to a SiconosMatrix
   */
  inline void setJacobianFIntqDotPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianFIntqDot = newPtr;
  }

  /** get  \f$ \nabla_{q}F_{gyr} \f$ , (pointer link)
   *
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianFGyrq() const { return _jacobianFGyrq; }

  
  /** get  \f$ \nabla_{\dot q}F_{gyr} \f$ , (pointer link)
   *
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianFGyrqDot() const
  {
    return _jacobianFGyrqDot;
  }

  /** get  \f$ \nabla_{q}F_{gyr} \f$ , (pointer link)
   *
   *  \param newPtr a SP SiconosMatrix
   */
  inline void setJacobianFGyrqPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianFGyrq = newPtr;
  }

  /** get  \f$ \nabla_{\dot q}F_{gyr} \f$ , (pointer link)
   *
   *  \param newPtr a SP SiconosMatrix
   */
  inline void setJacobianFGyrqDotPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianFGyrqDot = newPtr;
  }

  /** get  \f$  F(v,q,t,z) \f$  (pointer  link)
   *
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector forces() const override { return _forces; }

  /** get  \f$  \nabla_qF(v,q,t,z) \f$  (pointer  link)
   *
   *  \return pointer on a SiconosMatrix
   */
  virtual inline SP::SiconosMatrix jacobianqForces() const override
  {
    return _jacobianqForces;
  }

  /** get \f$ \nabla_{\dot q}F(v,q,t,z) \f$  (pointer  link)
   *
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianvForces() const override
  {
    return _jacobianqDotForces;
  }

  /** get all the values of the state vector q stored in memory.
   *  note: not const due to SchatzmanPaoliOSI::initializeWorkVectorsForDS
   *
   *  \return a memory
   */
  inline const SiconosMemory &qMemory() override { return _qMemory; }

  /** get all the values of the state vector velocity stored in memory.
   *  note: not const due to SchatzmanPaoliOSI::initializeWorkVectorsForDS
   *
   *  \return a memory
   */
  inline const SiconosMemory &velocityMemory() override
  {
    return _velocityMemory;
  }

  /** get all the values of the state vector p stored in memory
   *
   *  \param level
   *  \return a memory
   */
  inline const SiconosMemory &pMemory(unsigned int level)
  {
    return _pMemory[level];
  }

  /** get forces in memory buff
   *
   *  \return pointer on a SiconosMemory
   */
  inline const SiconosMemory &forcesMemory() override { return _forcesMemory; }

  /** initialize the SiconosMemory objects with a positive size.
   *
   *  \param size the size of the SiconosMemory. must be >= 0
   */
  void initMemory(unsigned int size) override;

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   *  \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory() override;

  /** allow to set a specified function to compute the mass
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this
   *  plugin
   */
  void setComputeMassFunction(const std::string &pluginPath,
                              const std::string &functionName)
  {
    _pluginMass->setComputeFunction(pluginPath, functionName);
    if (!_mass)
      _mass.reset(new SimpleMatrix(_ndof, _ndof));
    _hasConstantMass = false;
  }

  /** set a specified function to compute Mass
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeMassFunction(FPtr7 fct)
  {
    _pluginMass->setComputeFunction((void *)fct);
    if (!_mass)
      _mass.reset(new SimpleMatrix(_ndof, _ndof));
    _hasConstantMass = false;
  }

  /** allow to set a specified function to compute FInt
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this
   *  plugin
   */
  void setComputeFIntFunction(const std::string &pluginPath,
                              const std::string &functionName);

  /** set a specified function to compute fInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeFIntFunction(FPtr6 fct);

  /** allow to set a specified function to compute Fext
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this
   * plugin
   */
  void setComputeFExtFunction(const std::string &pluginPath,
                              const std::string &functionName)
  {
    _pluginFExt->setComputeFunction(pluginPath, functionName);
    if (!_fExt)
      _fExt.reset(new SiconosVector(_ndof));
    _hasConstantFExt = false;
  }

  /** set a specified function to compute fExt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeFExtFunction(VectorFunctionOfTime fct)
  {
    _pluginFExt->setComputeFunction((void *)fct);
    if (!_fExt)
      _fExt.reset(new SiconosVector(_ndof));
    //   computeFExtPtr = fct ;
    _hasConstantFExt = false;
  }

  /** allow to set a specified function to compute the inertia
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this
   * plugin
   */
  void setComputeFGyrFunction(const std::string &pluginPath,
                              const std::string &functionName);

  /** set a specified function to compute FGyr
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeFGyrFunction(FPtr5 fct);

  /** allow to set a specified function to compute the jacobian w.r.t q of the
   *  internal forces 
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFIntqFunction(const std::string &pluginPath,
                                       const std::string &functionName);
  /** allow to set a specified function to compute the jacobian
   *  of the internal forces w.r.t. q
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFIntqDotFunction(const std::string &pluginPath,
                                          const std::string &functionName);

  /** set a specified function to compute jacobian following q of the FInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFIntqFunction(FPtr6 fct);
  /** set a specified function to compute jacobian following qDot of the FInt
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFIntqDotFunction(FPtr6 fct);

  /** allow to set a specified function to compute the jacobian w.r.t q of the
   *  the external forces
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFGyrqFunction(const std::string &pluginPath,
                                       const std::string &functionName);

  /** allow to set a specified function to compute the jacobian w.r.t qDot of
   *  the the external strength
   *
   *  \param pluginPath std::string : the complete path to the plugin
   *  \param functionName std::string : the name of the function to use in this plugin
   */
  void setComputeJacobianFGyrqDotFunction(const std::string &pluginPath,
                                          const std::string &functionName);

  /** set a specified function to compute the jacobian following q of FGyr
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFGyrqFunction(FPtr5 fct);
  /** set a specified function to compute the jacobian following qDot of FGyr
   *
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianFGyrqDotFunction(FPtr5 fct);

  /** default function to compute the mass
   */
  void computeMass() override;

  /** function to compute the mass
   *
   *  \param position value used to evaluate the mass matrix
   */
  void computeMass(SP::SiconosVector position) override;

  /** default function to compute the internal strengths
   *
   *  \param time the current time
   */
  virtual void computeFInt(double time);

  /** function to compute the internal strengths
   *  with some specific values for position and velocity (ie not those of the
   *  current state). 
   *
   *  \param time the current time,
   *  \param position value used to evaluate the internal forces
   *  \param velocity value used to evaluate the internal forces
   */
  virtual void computeFInt(double time, SP::SiconosVector position,
                           SP::SiconosVector velocity);

  /** default function to compute the external strengths
   *
   *  \param time the current time
   */
  virtual void computeFExt(double time);

  /** default function to compute the inertia
   */
  virtual void computeFGyr();

  /** function to compute the inertia
   *  with some specific values for q and velocity (ie not those of the current state). 
   *
   *  \param position value used to evaluate the inertia forces
   *  \param velocity value used to evaluate the inertia forces
   */
  virtual void computeFGyr(SP::SiconosVector position,
                           SP::SiconosVector velocity);

  /** To compute the jacobian w.r.t q of the internal forces
   *
   *  \param time the current time
   */
  virtual void computeJacobianFIntq(double time);
  /** To compute the jacobian w.r.t qDot of the internal forces
   *
   *  \param time the current time
   */
  virtual void computeJacobianFIntqDot(double time);

  /** To compute the jacobian w.r.t q of the internal forces
   *
   *  \param time the current time
   *  \param position value used to evaluate the jacobian
   *  \param velocity value used to evaluate the jacobian
   */
  virtual void computeJacobianFIntq(double time, SP::SiconosVector position,
                                    SP::SiconosVector velocity);

  /** To compute the jacobian w.r.t. qDot of the internal forces
   *
   *  \param time the current time
   *  \param position value used to evaluate the jacobian
   *  \param velocity value used to evaluate the jacobian
   */
  virtual void computeJacobianFIntqDot(double time, SP::SiconosVector position,
                                       SP::SiconosVector velocity);

  /** function to compute the jacobian w.r.t. q of the inertia forces
   */
  virtual void computeJacobianFGyrq();

  /** function to compute the jacobian w.r.t. qDot of the inertia forces
   */
  virtual void computeJacobianFGyrqDot();

  /** function to compute the jacobian w.r.t. q of the inertia forces
   *
   *  \param position value used to evaluate the jacobian
   *  \param velocity value used to evaluate the jacobian
   */
  virtual void computeJacobianFGyrq(SP::SiconosVector position,
                                    SP::SiconosVector velocity);

  /** function to compute the jacobian w.r.t. qDot of the inertia forces
   *
   *  \param position value used to evaluate the jacobian
   *  \param velocity value used to evaluate the jacobian
   */
  virtual void computeJacobianFGyrqDot(SP::SiconosVector position,
                                       SP::SiconosVector velocity);

  /** default function to update the plugins functions using a new time:
   *
   *  \param time  the current time
   */
  void updatePlugins(double time) override{};

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
  void init_generalized_coordinates(unsigned int level);

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

  /** Allocate memory for forces and its jacobian.
   */
  void init_forces() override;

  ACCEPT_STD_VISITORS();
};

TYPEDEF_SPTR(LagrangianDS)

#endif // LAGRANGIANDS_H
