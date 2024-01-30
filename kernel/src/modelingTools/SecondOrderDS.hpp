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

/*! \file SecondOrder.hpp
  SecondOrderDS class - Second Order Non Linear Dynamical Systems.
*/

#ifndef SECONDORDERDS_H
#define SECONDORDERDS_H

#include <memory>

#include "DynamicalSystem.hpp"
#include "SiconosException.hpp"
#include "SimpleMatrix.hpp"

namespace siconos::algebra {

}

namespace siconos::modeling {

class BoundaryCondition;

/**

    Second Order non linear dynamical systems -  \f$ M(q,z) \dot v = F(v, q, t, z) + p \f$

    This class defines and computes a generic ndof-dimensional
    second order Non Linear Dynamical System of the form :

    \f[
     M(q,z) \dot v  = F(v, q, t,  z)  + p \\
     \dot q = G(q,v)
    \f]

    where

    - \f$ q \in R^{ndof} \f$ is the set of the coordinates,
    - \f$ \dot q =v \in R^{ndof} \f$ the velocity,
    - \f$ \ddot q = \dot v \in R^{ndof} \f$ the acceleration, i. e. the second
    time derivative of the generalized coordinates.
    - \f$ p \in R^{ndof} \f$ the reaction forces due to the Non Smooth Interaction.
    - \f$ M(q) \in R^{ndof \times ndof} \f$ is the inertia term (access : mass() method).
    - \f$ F( \dot q , q , t) \in R^{ndof} \f$ are the forces (access forces()
    method).
    - \f$ z \in R^{zSize} \f$ is a vector of arbitrary algebraic variables, some
    sort of discrete state.

    q[i] is the derivative number i of q.
    Thus: q[0]= \f$ q \f$, global coordinates, q[1]= \f$ \dot q \f$, velocity,
    q[2]= \f$ \ddot q \f$, acceleration.

    The following operators (and their jacobians) can be plugged, in the usual way
    (see User Guide, 'User-defined plugins')

    - \f$ M(q) \f$ (computeMass())
    - \f$ F(v , q , t, z) \f$ (computeF())

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

    \f]

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
class SecondOrderDS : public DynamicalSystem {
 protected:
  ACCEPT_SERIALIZATION(SecondOrderDS);

  // -- MEMBERS --

  /** number of degrees of freedom of the system */
  unsigned int _ndof{0};

  /** mass of the system */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _mass{nullptr};

  /** true if the  mass matrix is constant */
  bool _hasConstantMass = false;

  /** inverse or factorization of the mass of the system */
  std::shared_ptr<siconos::algebra::SimpleMatrix> _inverseMass{nullptr};

  /** "Reaction", generalized forces or imuplses due to the non smooth law
   * The index corresponds to the kinematic
   * level of the corresponding constraints. It mainly depends on what the
   * simulation part want to store, but some rules have to be followed. For
   * instance :
   *  - for the constraints at the acceleration level, _p[2] stores the reaction
   * forces,
   *  - for the constraints at the veocity level,  _p[1] stores the (discrete)
   * reaction impulse
   *  - for the constraints at the position level, _p[0] stores the multiplier
   * for a constraint in position
   */
  std::vector<std::shared_ptr<siconos::algebra::SiconosVector>> _p = {nullptr, nullptr,
                                                                      nullptr};

  /** Initial position */
  std::shared_ptr<siconos::algebra::SiconosVector> _q0{nullptr};

  /** Boundary condition applied to a dynamical system*/
  std::shared_ptr<siconos::modeling::BoundaryCondition> _boundaryConditions{nullptr};

  /** Reaction to an applied  boundary condition */
  std::shared_ptr<siconos::algebra::SiconosVector> _reactionToBoundaryConditions{nullptr};

  // Rule of five
  SecondOrderDS() = default;
  SecondOrderDS(const SecondOrderDS &) = delete;
  SecondOrderDS &operator=(const SecondOrderDS &) = delete;
  SecondOrderDS(SecondOrderDS &&) = delete;
  SecondOrderDS &operator=(SecondOrderDS &&) = delete;

  /** minimal constructor, from state dimension
   *  result in \f$ \dot x = r \f$
   *
   *  \param dimension size of the system (n)
   */
  SecondOrderDS(unsigned int dimension, unsigned int ndof)
      : DynamicalSystem(dimension), _ndof(ndof), _hasConstantMass(true){};

 public:
  /** destructor */
  virtual ~SecondOrderDS() noexcept = default;

  /** get p
   *
   *  \param level unsigned int, required level for p, default = 2
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> p(unsigned int level = 2) const
  {
    return _p[level];
  }

  /** get mass matrix (pointer link)
   *
   *  \return std::shared_ptr<siconos::algebra::SiconosMatrix>
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> mass() const { return _mass; }

  /** get (pointer) inverse or LU-factorization of the mass,
   *  used for LU-forward-backward computation
   *
   *  \return pointer std::shared_ptr<siconos::algebra::SimpleMatrix>
   */
  inline std::shared_ptr<siconos::algebra::SimpleMatrix> inverseMass() const
  {
    return _inverseMass;
  }

  /** set mass to pointer newPtr
   *
   *  \param newPtr a plugged matrix SP
   */
  void setMassPtr(std::shared_ptr<siconos::algebra::SimpleMatrix> newPtr);

  /** set the value of the right-hand side, \f$ \dot x \f$
   *
   *  \param newValue siconos::algebra::SiconosVector
   */
  void setRhs(const siconos::algebra::SiconosVector &newValue) override
  {
    THROW_EXCEPTION("SecondOrderDS - setRhs call is forbidden for 2nd order systems.");
  }

  /** set right-hand side, \f$ \dot x \f$ (pointer link)
   *
   *  \param newPtr std::shared_ptr<siconos::algebra::SiconosVector>
   */
  void setRhsPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) override
  {
    THROW_EXCEPTION("SecondOrderDS - setRhsPtr call is forbidden for 2nd order systems.");
  }

  /* function to compute \f$ F(v,q,t,z) \f$ for the current state
   *
   *  \param time the current time
   */
  // virtual void computeForces(double time);

  /** Compute \f$ F(v,q,t,z) \f$
   *
   *  \param time the current time
   *  \param q std::shared_ptr<siconos::algebra::SiconosVector>: pointers on q
   *  \param velocity std::shared_ptr<siconos::algebra::SiconosVector>: pointers on velocity
   */
  virtual void computeForces(double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
                             std::shared_ptr<siconos::algebra::SiconosVector> velocity) = 0;

  /** Compute \f$ \nabla_qF(v,q,t,z) \f$ for current \f$ q,v \f$
   *  Default function to compute forces
   *
   *  \param time the current time
   */
  virtual void computeJacobianqForces(double time) = 0;

  /** Compute \f$ \nabla_{\dot q}F(v,q,t,z) \f$ for current \f$ q,v \f$
   *
   *  \param time the current time
   */
  virtual void computeJacobianvForces(double time) = 0;

  /** return the number of degrees of freedom of the system
   *
   *  \return an unsigned int.
   */
  inline unsigned int dimension() const override { return _ndof; }

  /** generalized coordinates of the system (vector of size dimension())
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> q() const = 0;

  /** set value of generalized coordinates vector (copy)
   *
   *  \param newValue
   */
  virtual void setQ(const siconos::algebra::SiconosVector &newValue) = 0;

  /** set value of generalized coordinates vector (pointer link)
   *
   *  \param newPtr
   */
  virtual void setQPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) = 0;

  /** get initial state (pointer link)
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  std::shared_ptr<siconos::algebra::SiconosVector> q0() const { return _q0; }

  /** set initial state (copy)
   *
   *  \param newValue
   */
  virtual void setQ0(const siconos::algebra::SiconosVector &newValue) = 0;

  /** set initial state (pointer link)
   *
   *  \param newPtr
   */
  virtual void setQ0Ptr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) = 0;

  /** get velocity vector (pointer link)
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> velocity() const = 0;

  /** set velocity vector (copy)
   *
   *  \param newValue
   */
  virtual void setVelocity(const siconos::algebra::SiconosVector &newValue) = 0;

  /** set velocity vector (pointer link)
   *
   *  \param newPtr
   */
  virtual void setVelocityPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) = 0;

  /** get initial velocity (pointer)
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> velocity0() const = 0;

  /** set initial velocity (copy)
   *
   *  \param newValue
   */
  virtual void setVelocity0(const siconos::algebra::SiconosVector &newValue) = 0;

  /** set initial velocity (pointer link)
   *
   *  \param newPtr
   */
  virtual void setVelocity0Ptr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) = 0;

  /** get acceleration (pointer link)
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> acceleration() const = 0;

  /** get \f$ F(v,q,t,z) \f$ (pointer  link)
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> forces() const = 0;

  /** \return \f$ \nabla_qF(v,q,t,z) \f$ (pointer  link) */
  virtual std::shared_ptr<siconos::algebra::SiconosMatrix> jacobianqForces() const = 0;

  /** get \f$ \nabla_{\dot q}F(v,q,t,z) \f$ (pointer  link)
   *
   *  \return pointer on a SiconosMatrix
   */
  virtual std::shared_ptr<siconos::algebra::SiconosMatrix> jacobianvForces() const = 0;

  /** get all the values of the state vector q stored in memory.
   *  note: not const due to SchatzmanPaoliOSI::initializeWorkVectorsForDS
   *
   *  \return a memory
   */
  virtual const siconos::algebra::SiconosMemory &qMemory() = 0;

  /** get all the values of the state vector velocity stored in memory.
   *  note: not const due to SchatzmanPaoliOSI::initializeWorkVectorsForDS
   *
   *  \return a memory
   */
  virtual const siconos::algebra::SiconosMemory &velocityMemory() = 0;

  /** get forces in memory buff
   *
   *  \return pointer on a SiconosMemory
   */
  virtual const siconos::algebra::SiconosMemory &forcesMemory() = 0;

  /** initialize the SiconosMemory objects with a positive size.
   *
   *  \param size the size of the SiconosMemory. must be >= 0
   */
  void initMemory(unsigned int size) override = 0;

  /** default function to compute the mass
   */
  virtual void computeMass() = 0;

  /** function to compute the mass
   *
   *  \param position value used to evaluate the mass matrix
   */
  virtual void computeMass(std::shared_ptr<siconos::algebra::SiconosVector> position) = 0;

  /** set Boundary Conditions
   *
   *  \param newbd BoundaryConditions
   */
  virtual void setBoundaryConditions(
      std::shared_ptr<siconos::modeling::BoundaryCondition> newbd);

  /** get Boundary Conditions
   *
   *  \return std::shared_ptr<siconos::modeling::BoundaryCondition> pointer on a
   * BoundaryConditions
   */
  inline std::shared_ptr<siconos::modeling::BoundaryCondition> boundaryConditions()
  {
    return _boundaryConditions;
  };

  /** set Reaction to Boundary Conditions
   *
   *  \param newrbd BoundaryConditions pointer
   */
  inline void setReactionToBoundaryConditions(
      std::shared_ptr<siconos::algebra::SiconosVector> newrbd)
  {
    _reactionToBoundaryConditions = newrbd;
  };

  /** get Reaction to  Boundary Conditions
   *
   *  \return pointer on a BoundaryConditions
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> reactionToBoundaryConditions()
  {
    return _reactionToBoundaryConditions;
  };

  /**
      Allocate memory for the lu factorization of the mass of the system.
      Useful for some integrators with system inversion involving the mass
  */
  virtual void init_inverse_mass() = 0;

  /**
      Update the content of the lu factorization of the mass of the system,
      if required.
  */
  virtual void update_inverse_mass() = 0;

  /** Allocate memory for forces and its jacobian.
   */
  virtual void init_forces() = 0;

  //  ACCEPT_STD_VISITORS();
};
}  // namespace siconos::modeling

#endif  // LAGRANGIANDS_H
