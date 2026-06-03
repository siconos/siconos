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

namespace siconos::modeling {

class BoundaryCondition;

/**

    Second Order non linear dynamical systems -  \f$ M(q) \dot v = F_{total}(v, q, t) + p \f$

    This class defines and computes a generic ndof-dimensional
    second order Non Linear Dynamical System of the form :

    \f[
     M(q) \dot v  = F_{total}(v, q, t)  + p \\
     \dot q = G(q,v)
    \f]

    where

    - \f$ q \in R^{ndof} \f$ is the set of the coordinates,
    - \f$ \dot q =v \in R^{ndof} \f$ the velocity,
    - \f$ \ddot q = \dot v \in R^{ndof} \f$ the acceleration, i. e. the second
    time derivative of the generalized coordinates.
    - \f$ p \in R^{ndof} \f$ the reaction forces due to the Non Smooth Interaction.
    - \f$ M(q) \in R^{ndof \times ndof} \f$ is the inertia term (access : mass() method).
    - \f$ F_{total}( \dot q , q , t) \in R^{ndof} \f$ are the forces (access totalForces()
    method).
    - \f$ z \in R^{zSize} \f$ is a vector of arbitrary algebraic variables, some
    sort of discrete state.

    q[i] is the derivative number i of q.
    Thus: q[0]= \f$ q \f$, global coordinates, q[1]= \f$ \dot q \f$, velocity,
    q[2]= \f$ \ddot q \f$, acceleration.

    The following operators (and their jacobians) can be plugged, in the usual way
    (see User Guide, 'User-defined plugins')

    - \f$ M(q) \f$ (computeMass(...))
    - \f$ F_{total}(v , q , t) \f$ (computeTotalForces(...))

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

    \f]

    - jacobian of the rhs, with respect to x

    \f[
       \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
       0  & I \\
       \nabla_{q}(M^{-1}(q)F_{total}l}(v, q , t)) &  \nabla_{\dot q}(M^{-1}(q)F_{total}(v, q ,
       t)) \\ \end{array}\right]

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
    - computeRhs(), computeJacobianRhsOver_x() ..., to update the content of rhs, its
    jacobians ...

*/
class SecondOrderDS : public DynamicalSystem {
 protected:
  ACCEPT_SERIALIZATION(SecondOrderDS);

  // -- MEMBERS --

  /** number of degrees of freedom of the system */
  siconos::algebra::Index ndof_{0};

  /** "Reaction", generalized forces or impulses due to the non smooth law
   * The index corresponds to the kinematic
   * level of the corresponding constraints. It mainly depends on what the
   * simulation part want to store, but some rules have to be followed. For
   * instance :
   *  - for the constraints at the acceleration level, p_[2] stores the reaction
   * forces,
   *  - for the constraints at the veocity level,  p_[1] stores the (discrete)
   * reaction impulse
   *  - for the constraints at the position level, p_[0] stores the multiplier
   * for a constraint in position
   */
  siconos::algebra::blocks::SharedVector p_ = {nullptr, nullptr, nullptr};

  bool hasLUMass_{false};

  /** Boundary condition applied to a dynamical system*/
  std::shared_ptr<siconos::modeling::BoundaryCondition> boundaryConditions_{nullptr};

  /** Reaction to an applied  boundary condition */
  std::shared_ptr<siconos::algebra::SiconosVector> reactionToBoundaryConditions_{nullptr};

  // Rule of five
  SecondOrderDS() = default;
  SecondOrderDS(const SecondOrderDS&) = delete;
  SecondOrderDS& operator=(const SecondOrderDS&) = delete;
  SecondOrderDS(SecondOrderDS&&) = delete;
  SecondOrderDS& operator=(SecondOrderDS&&) = delete;

  /** minimal constructor, from state dimension
   *  result in \f$ \dot x = r \f$
   *
   *  \param dimension dimension of corresponding first order system
   *  \param ndof number of degrees of freedom
   */
  SecondOrderDS(siconos::algebra::Index dimension, siconos::algebra::Index ndof)
      : DynamicalSystem(dimension), ndof_(ndof){};

 public:
  /** destructor */
  virtual ~SecondOrderDS() noexcept = default;

  /** \return a read-only view onto the nonsmooth force or impulse
   *
   *  \param level required level for p, default = 2
   */
  inline auto p_read(unsigned int level = 2) const {
    return siconos::algebra::ConstMapVectorType(p_[level]->data(), p_[level]->size());
  }

  /** \return nonsmooth force or impulse (pointer link)
   *
   *  \param level unsigned int, required level for p, default = 2
   */
  std::shared_ptr<siconos::algebra::SiconosVector> p(unsigned int level = 2) const {
    return p_[level];
  }

  /** \return nonsmooth force or impulse (only for pybind11/python bindings)
   *
   *  \param level unsigned int, required level for p, default = 2
   */
  siconos::algebra::SiconosVector& p_python(unsigned int level = 2) const {
    return *(p_[level]);
  }

  /** \return the number of degrees of freedom of the system */
  inline siconos::algebra::Index dimension() const override { return ndof_; }

  /** Fake function to access the dimension used to allocate iteration matrix
   *  in the integrators
   * Usually equal to ndof but might be larger (see fem).
   * \return number of degrees of freedom of the system
   */
  inline virtual siconos::algebra::Index real_size() const override { return ndof_; }

  /** \return a read-only view on the generalized coordinates of the system */
  virtual const siconos::algebra::ConstMapVectorType q_read() const = 0;

  /** \return the generalized coordinates of the system (pointer link) */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> q() const = 0;

  // FP: override SecondOrderDS. Used only in visitors of MechanicsIO. To be reviewed ...
  virtual const siconos::algebra::ConstMapVectorType velocity_read() const = 0;

  /** \return a read-only view on acceleration vector */
  virtual const siconos::algebra::ConstMapVectorType acceleration_read() const = 0;

  /** \return the acceleration vector (pointer link) */
  virtual std::shared_ptr<siconos::algebra::SiconosVector> acceleration() const = 0;

  /** get all the values of the state vector q stored in memory.
   *  note: not const due to SchatzmanPaoliOSI::initializeWorkVectorsForDS
   *
   *  \return a memory
   */
  virtual const siconos::algebra::SiconosMemory& qMemory() = 0;

  /** get forces in memory buff
   *
   *  \return pointer on a SiconosMemory
   */
  virtual const siconos::algebra::SiconosMemory& forcesMemory() = 0;

  /** initialize the SiconosMemory objects with a positive size.
   *
   *  \param size the size of the SiconosMemory. must be >= 0
   */
  void initMemory(siconos::algebra::blocks::size_type size) override = 0;

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
  inline std::shared_ptr<siconos::modeling::BoundaryCondition> boundaryConditions() const {
    return boundaryConditions_;
  };

  /** set Reaction to Boundary Conditions
   *
   *  \param newrbd BoundaryConditions pointer
   */
  inline void setReactionToBoundaryConditions(
      std::shared_ptr<siconos::algebra::SiconosVector> newrbd) {
    reactionToBoundaryConditions_ = newrbd;
  };

  /** get Reaction to  Boundary Conditions
   *
   *  \return pointer on a BoundaryConditions
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> reactionToBoundaryConditions() {
    return reactionToBoundaryConditions_;
  };

  /** get Reaction to  Boundary Conditions
   *
   *  \return pointer on a BoundaryConditions
   */
  inline siconos::algebra::SiconosVector& reactionToBoundaryConditions_python() {
    return *reactionToBoundaryConditions_;
  };

  /**
      Allocate memory and lu-factorize the mass of the system.
      Useful for some integrators with system inversion involving the mass
  */
  virtual void init_lu_mass() = 0;

  bool hasLUMass() const { return hasLUMass_; }
};
}  // namespace siconos::modeling

#endif  // LAGRANGIANDS_H
