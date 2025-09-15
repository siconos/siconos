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

/*! \file RigidBody2dDS.hpp
  \brief Definition of an abstract 3D rigid body above NewtonEulerDS
*/


#ifndef RigidBody2dDS_h
#define RigidBody2dDS_h

#include <MechanicsFwd.hpp>
#include <LagrangianLinearTIDS.hpp>
#include <SiconosVisitor.hpp>
#include <SiconosContactor.hpp>

class RigidBody2dDS : public LagrangianLinearTIDS,
                      public std::enable_shared_from_this<RigidBody2dDS>
{
protected:
  ACCEPT_SERIALIZATION(RigidBody2dDS);

  RigidBody2dDS() : LagrangianLinearTIDS() {};

  /** a scalar mass in the case of RigidBody2dDS */
  double _scalarMass;

  SP::SiconosContactorSet _contactors;
  bool _useContactorInertia;

  /** If false, bodies connected to this body by a joint will not
   * collide. See also NewtonEulerJointR::_allowSelfCollide */
  bool _allowSelfCollide = true;

public:

  RigidBody2dDS(SP::SiconosVector position,
                SP::SiconosVector velocity,
                SP::SiconosMatrix mass = SP::SimpleMatrix());

  RigidBody2dDS(SP::SiconosVector position,
		SP::SiconosVector velocity,
		double mass,
		double inertia);

  virtual ~RigidBody2dDS();

  double scalarMass()
  {
    return _scalarMass;
  };


  void setUseContactorInertia(bool use) { _useContactorInertia = use; }

  bool useContactorInertia() { return _useContactorInertia; }

  /** Return the value of the _allowSelfCollide flag. */
  bool allowSelfCollide() { return _allowSelfCollide; }

  /** Set the value of the _allowSelfCollide flag. */
  void setAllowSelfCollide(bool x) { _allowSelfCollide = x; }

  /** \return the contactor set associated with this body */
  SP::SiconosContactorSet contactors() const { return _contactors; }

  /** Provide a set of contactors to the body.
   *
   *  \param c A SP::SiconosContactorSet */
  void setContactors(SP::SiconosContactorSet c) { _contactors = c; }

  /** Make the base position of the contactors equal to the DS q vector.
   *
   *  \return a SP::SiconosVector */
  virtual SP::SiconosVector base_position() { return q(); }

  ACCEPT_BASE_VISITORS(LagrangianLinearTIDS);
  
};

#endif /* RigidBody2dDS_h */
