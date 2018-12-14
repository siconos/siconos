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

/*! \file SecondOrder.hpp
  SecondOrderDS class - Second Order Non Linear Dynamical Systems.
*/

#ifndef SECONDORDERDS_H
#define SECONDORDERDS_H

#include "DynamicalSystem.hpp"
#include "BoundaryCondition.hpp"
#include "SiconosConst.hpp"

/** Second Order non linear dynamical systems - \f$M(q,z) \dot v = F(v, q, t, z) + p \f$

 This class defines and computes a generic ndof-dimensional
 second order Non Linear Dynamical System of the form :

 \rst
  .. math::

     M(q,z) \\dot v  = F(v, q, t,  z)  + p \\
     \\dot q = G(q,v)

 \endrst

 where

 - \f$q \in R^{ndof} \f$ is the set of the coordinates,
 - \f$ \\dot q =v \in R^{ndof} \f$ the velocity,
 - \f$ \ddot q =\\dot v \in R^{ndof} \f$ the acceleration, i. e. the second
 time derivative of the generalized coordinates.
 - \f$ p \in R^{ndof} \f$ the reaction forces due to the Non Smooth
 Interaction.
 - \f$ M(q) \in R^{ndof \times ndof} \f$ is the inertia term (access : mass() method).
 - \f$ F(\\dot q , q , t) \in R^{ndof} \f$ are the forces (access forces() method).
 - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discrete state.

 q[i] is the derivative number i of q.
 Thus: q[0]=\f$ q \f$, global coordinates, q[1]=\f$ \\dot q\f$, velocity, q[2]=\f$ \ddot q \f$, acceleration.

 The following operators (and their jacobians) can be plugged, in the usual way (see User Guide, 'User-defined plugins')

 - \f$M(q)\f$ (computeMass())
 - \f$F(v , q , t, z)\f$ (computeF())
3
 If required (e.g. for Event-Driven like simulation), formulation as a first-order system is also available, and writes:

 - \f$ n= 2 ndof \f$
 - \f$ x = \left[\begin{array}{c}q \\ \\dot q\end{array}\right]\f$
 - rhs given by:

   \rst

   .. math::

      \\dot x = \left[\begin{array}{c}
      \\dot q\\
      \ddot q = M^{-1}(q)\left[F(v, q , t, z) + p \right]\\
      \end{array}\right]

   \endrst

 - jacobian of the rhs, with respect to x

   \rst

    .. math::

       \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
       0  & I \\
       \nabla_{q}(M^{-1}(q)F(v, q , t, z)) &  \nabla_{\\dot q}(M^{-1}(q)F(v, q , t, z)) \\
       \end{array}\right]

     \endrst

     with the input due to the non smooth law:

     \rst

     .. math::

      \left[\begin{array}{c}
      0 \\
      p \end{array}\right]

      \endrst

  In that case, use the following methods:
    - initRhs() to allocate/initialize memory for these new operators,
    - rhs() to get the rhs vector
    - computeRhs(), computeJacobianRhsx() ..., to update the content of rhs, its jacobians ...

*/
class SecondOrderDS : public DynamicalSystem
{

protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(SecondOrderDS);

  // -- MEMBERS --

  /** number of degrees of freedom of the system */
  unsigned int _ndof;

  /** mass of the system */
  SP::SimpleMatrix _mass;

  /** true if the  mass matrix is constant */
  bool _hasConstantMass;

  /** inverse or factorization of the mass of the system */
  SP::SimpleMatrix _inverseMass;

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

  /** Boundary condition applied to a dynamical system*/
  SP::BoundaryCondition _boundaryConditions;

  /** Reaction to an applied  boundary condition */
  SP::SiconosVector _reactionToBoundaryConditions;

  // /** Default constructor */
  SecondOrderDS():DynamicalSystem(Type::SecondOrderDS) {};

  /** minimal constructor, from state dimension
      result in \f$ \dot x = r \f$
   *  \param dimension size of the system (n)
   */
  SecondOrderDS(unsigned int dimension, unsigned int ndof):DynamicalSystem(dimension),
                                                           _ndof(ndof),  _hasConstantMass(true)
  {};

public:

  /** destructor */
  virtual ~SecondOrderDS() {};

  /** get p
   *  \param level unsigned int, required level for p, default = 2
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector p(unsigned int level = 2) const
  {
    return _p[level];
  }

  /** get mass matrix (pointer link)
   *  \return SP::SiconosMatrix
   */
  inline SP::SimpleMatrix mass() const
  {
    return _mass;
  }

  /** get (pointer) inverse or LU-factorization of the mass,
   * used for LU-forward-backward computation
   *  \return pointer SP::SimpleMatrix
   */
  inline SP::SimpleMatrix inverseMass() const
  {
    return _inverseMass;
  }

  /** set mass to pointer newPtr
   *  \param newPtr a plugged matrix SP
   */
  inline void setMassPtr(SP::SimpleMatrix newPtr)
  {
    _mass = newPtr;
    _hasConstantMass = true;
  }


  /*! @name Right-hand side computation */
  //@{

  /** reset the state to the initial state */
  virtual void resetToInitialState()=0;

  /** allocate (if needed)  and compute rhs and its jacobian.
   * \param time of initialization
   */
  virtual void initRhs(double time)=0;

  /** set nonsmooth input to zero
   *  \param level input-level to be initialized.
   */
  virtual void initializeNonSmoothInput(unsigned int level) =0 ;

  /** update right-hand side for the current state
   *  \param time of interest
   */
  virtual void computeRhs(double time) =0;

  /** update \f$\nabla_x rhs\f$ for the current state
   *  \param time of interest
   */
  virtual void computeJacobianRhsx(double time) =0;

  /** reset non-smooth part of the rhs (i.e. p), for all 'levels' */
  void resetAllNonSmoothParts() =0;

  /** set nonsmooth part of the rhs (i.e. p) to zero for a given level
   * \param level
   */
  void resetNonSmoothPart(unsigned int level) =0;

  /** set the value of the right-hand side, \f$ \dot x \f$
   *  \param newValue SiconosVector
   */
  void setRhs(const SiconosVector& newValue)
  {
    RuntimeException::selfThrow("SecondOrderDS - setRhs call is forbidden for 2nd order systems.");
  }

  /** set right-hand side, \f$ \dot x \f$ (pointer link)
   *  \param newPtr SP::SiconosVector
   */
  void setRhsPtr(SP::SiconosVector newPtr)
  {
    RuntimeException::selfThrow("SecondOrderDS - setRhsPtr call is forbidden for 2nd order systems.");
  }

  /** function to compute \f$F(v,q,t,z)\f$ for the current state
   *  \param time the current time
   */
  //virtual void computeForces(double time);

  /** Compute \f$F(v,q,t,z)\f$
   *  \param time the current time
   *  \param q SP::SiconosVector: pointers on q
   *  \param velocity SP::SiconosVector: pointers on velocity
   */
  virtual void computeForces(double time,
                             SP::SiconosVector q,
                             SP::SiconosVector velocity) =0;

  /** Compute \f$\nabla_qF(v,q,t,z)\f$ for current \f$q,v\f$
      Default function to compute forces
   *  \param time the current time
   */
  virtual void computeJacobianqForces(double time) =0;

  /** Compute \f$\nabla_{\dot q}F(v,q,t,z)\f$ for current \f$q,v\f$
   *  \param time the current time
   */
  virtual void computeJacobianvForces(double time)=0;

  ///@}

  /*! @name Attributes access
    @{ */

  /** return the number of degrees of freedom of the system
   *  \return an unsigned int.
   */
  inline unsigned int dimension() const
  {
    return _ndof;
  }
  /** generalized coordinates of the system (vector of size dimension())
   *  \return pointer on a SiconosVector
   */
  virtual SP::SiconosVector q() const =0;

  /** set value of generalized coordinates vector (copy)
   *  \param newValue
   */
  virtual void setQ(const SiconosVector& newValue) =0;

  /** set value of generalized coordinates vector (pointer link)
   *  \param newPtr
   */
  virtual void setQPtr(SP::SiconosVector newPtr) =0;

  /** return initial state of the system
   *  \return pointer on a SiconosVector
   */
  virtual  SP::SiconosVector q0() const =0;

  /** set initial state (copy)
   *  \param newValue
   */
  virtual void setQ0(const SiconosVector& newValue) =0;

  /** set initial state (pointer link)
   *  \param newPtr
   */
  virtual void setQ0Ptr(SP::SiconosVector newPtr) =0;

  /** get velocity vector (pointer link)
   *  \return pointer on a SiconosVector
   */
  virtual inline SP::SiconosVector velocity() const =0;

  /** set velocity vector (copy)
   *  \param newValue
   */
  virtual void setVelocity(const SiconosVector& newValue) =0;

  /** set velocity vector (pointer link)
   *  \param newPtr
   */
  virtual void setVelocityPtr(SP::SiconosVector newPtr) =0;

  /** get initial velocity (pointer)
   *  \return pointer on a SiconosVector
   */
  virtual  SP::SiconosVector velocity0() const =0;
  /** set initial velocity (copy)
   *  \param newValue
   */
  virtual void setVelocity0(const SiconosVector& newValue) =0;

  /** set initial velocity (pointer link)
   *  \param newPtr
   */
  virtual void setVelocity0Ptr(SP::SiconosVector newPtr) =0;

  /** get acceleration (pointer link)
   *  \return pointer on a SiconosVector
   */
  virtual SP::SiconosVector acceleration() const =0;

  /** get \f$ F(v,q,t,z)\f$ (pointer  link)
   *  \return pointer on a SiconosVector
   */
  virtual  SP::SiconosVector forces() const =0;

  /** get \f$ \nabla_qF(v,q,t,z)\f$ (pointer  link)
   *  \return pointer on a SiconosMatrix
   */
  virtual  SP::SimpleMatrix jacobianqForces() const =0;

  /** get \f$ \nabla_{\dot q}F(v,q,t,z)\f$ (pointer  link)
   *  \return pointer on a SiconosMatrix
   */
  virtual  SP::SimpleMatrix jacobianvForces() const =0;

  ///@}

  /*! @name Memory vectors management
    @{ */

  /** get all the values of the state vector q stored in memory.
   * note: not const due to SchatzmanPaoliOSI::initializeWorkVectorsForDS
   *  \return a memory
   */
  virtual const SiconosMemory& qMemory() =0;

  /** get all the values of the state vector velocity stored in memory.
   * note: not const due to SchatzmanPaoliOSI::initializeWorkVectorsForDS
   *  \return a memory
   */
  virtual const SiconosMemory& velocityMemory() =0;


  /** get forces in memory buff
   *  \return pointer on a SiconosMemory
   */
  virtual  const SiconosMemory& forcesMemory() =0;

  /** initialize the SiconosMemory objects with a positive size.
   *  \param size the size of the SiconosMemory. must be >= 0
   */
  virtual void initMemory(unsigned int size) =0;

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  virtual void swapInMemory() =0;

  ///@}

  /** default function to compute the mass
   */
  virtual void computeMass() =0;

  /** function to compute the mass
   *  \param position value used to evaluate the mass matrix
   */
  virtual void computeMass(SP::SiconosVector position) =0;

  /**default function to update the plugins functions using a new time:
   * \param time  the current time
   */
  virtual void updatePlugins(double time) =0;

  ///@}

  /** print the data of the dynamical system on the standard output
   */
  virtual void display() const =0;

  /** set Boundary Conditions
   *  \param newbd BoundaryConditions
   */
  virtual void setBoundaryConditions(SP::BoundaryCondition newbd);

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


  /** Allocate memory for the lu factorization of the mass of the system.
      Useful for some integrators with system inversion involving the mass
  */
  virtual void init_inverse_mass() =0;

  /** Update the content of the lu factorization of the mass of the system,
      if required.
  */
  virtual void update_inverse_mass() =0;

  /** Allocate memory for forces and its jacobian.
  */
  virtual void init_forces()=0;

  ///@}

  ACCEPT_STD_VISITORS();

};

//TYPEDEF_SPTR(SecondOrderDS)

#endif // LAGRANGIANDS_H
