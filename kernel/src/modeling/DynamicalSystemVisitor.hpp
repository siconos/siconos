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
 * Unless required by applicable law or agreed to in writing, softwa∏re
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/*! \file DynamicalSystemVisitor.hpp
  Abstract class - General interface for visitors of  Dynamical Systems.
*/

#ifndef DSVisitor_H
#define DSVisitor_H

#include <memory>

#include "SiconosException.hpp"
#include "SiconosVector.hpp"

// For any DS class that might be visited, visit functions must be declared here.
// Even those for DS in other components like mechanics.
//

#ifdef HAVE_SICONOS_MECHANICS
namespace siconos::collision {
class RigidBody2dDS;
class RigidBodyDS;

namespace native::bodies {
class Disk;
class ExternalBody;
class CircularDS;
class Circle;
class SphereLDS;
class SphereNEDS;
}  // namespace native::bodies

}  // namespace siconos::collision

namespace siconos::fem::cable {
class CableDS;
}

#ifdef SICONOS_HAS_OpenCASCADE
namespace siconos::mechanics::occ {
class OccBody;
}
#endif
#endif

namespace siconos::modeling {

class DynamicalSystem;
class LagrangianDS;
class NewtonEulerDS;
class LagrangianLinearTIDS;

namespace dynamical_systems {

struct Visitor {
  virtual ~Visitor() noexcept = default;

  // Define this function in your visitor if you need to handle
  // internal memory (view) without copy. See example in MechanicsIO.cpp
  virtual void setMap(Eigen::Ref<siconos::algebra::SiconosVector>){};

  virtual void visit(std::shared_ptr<siconos::modeling::DynamicalSystem>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "DynamicalSystem in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(siconos::modeling::DynamicalSystem &) {
    THROW_EXCEPTION(
        "you must define a visit function for DynamicalSystem in "
        "a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(const siconos::modeling::DynamicalSystem &) {
    THROW_EXCEPTION(
        "you must define a visit function for DynamicalSystem in "
        "a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::modeling::LagrangianDS>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: LagrangianDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(siconos::modeling::LagrangianDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for LagrangianDS in a "
        "derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::modeling::LagrangianDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for LagrangianDS in a "
        "derived class of dynamical_systems::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::modeling::LagrangianLinearTIDS>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: LagrangianLinearTIDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(siconos::modeling::LagrangianLinearTIDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for LagrangianLinearTIDS in a "
        "derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::modeling::LagrangianLinearTIDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for LagrangianLinearTIDS in a "
        "derived class of dynamical_systems::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::modeling::NewtonEulerDS>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: NewtonEulerDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(siconos::modeling::NewtonEulerDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for NewtonEulerDS in a "
        "derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::modeling::NewtonEulerDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for NewtonEulerDS in a "
        "derived class of dynamical_systems::Visitor");
  }

#ifdef HAVE_SICONOS_MECHANICS
  virtual void visit(std::shared_ptr<siconos::collision::RigidBody2dDS>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "RigidBody2dDS in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(siconos::collision::RigidBody2dDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for RigidBody2dDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::collision::RigidBody2dDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for RigidBody2dDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::RigidBodyDS>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "RigidBodyDS in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(siconos::collision::RigidBodyDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for RigidBodyDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::collision::RigidBodyDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for RigidBodyDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::SphereNEDS>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "SphereNEDS in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereNEDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereNEDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::ExternalBody>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: ExternalBody "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(siconos::collision::native::bodies::ExternalBody &) {
    THROW_EXCEPTION(
        "you must define a visit function for ExternalBody in a "
        "derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::ExternalBody &) {
    THROW_EXCEPTION(
        "you must define a visit function for ExternalBody in a "
        "derived class of dynamical_systems::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::Disk>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "Disk in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::Disk &) {
    THROW_EXCEPTION(
        "you must define a visit function for Disk "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::Disk &) {
    THROW_EXCEPTION(
        "you must define a visit function for Disk "
        "in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::Circle>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "Circle in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::Circle &) {
    THROW_EXCEPTION(
        "you must define a visit function for Circle "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::Circle &) {
    THROW_EXCEPTION(
        "you must define a visit function for Circle "
        "in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::SphereLDS>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "SphereLDS in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereLDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereLDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDS "
        "in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::fem::cable::CableDS>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "CableDS in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(siconos::fem::cable::CableDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for CableDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::fem::cable::CableDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for CableDS "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(siconos::collision::native::bodies::CircularDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for CircularDS in a "
        "derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::CircularDS &) {
    THROW_EXCEPTION(
        "you must define a visit function for CircularDS in a "
        "derived class of dynamical_systems::Visitor");
  }
#ifdef SICONOS_HAS_OpenCASCADE
  virtual void visit(std::shared_ptr<siconos::mechanics::occ::OccBody>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "OccBody in a derived class of dynamical_systems::Visitor");
  }

  virtual void visit(siconos::mechanics::occ::OccBody &) {
    THROW_EXCEPTION(
        "you must define a visit function for OccBody "
        "in a derived class of dynamical_systems::Visitor");
  }
  virtual void visit(const siconos::mechanics::occ::OccBody &) {
    THROW_EXCEPTION(
        "you must define a visit function for OccBody "
        "in a derived class of dynamical_systems::Visitor");
  }
#endif
#endif
};
}  // namespace dynamical_systems
}  // namespace siconos::modeling
#endif