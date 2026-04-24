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
/*! \file DynamicalSystem.hpp
  Abstract class - General interface for all Dynamical Systems.
*/

#ifndef RelVisitor_H
#define RelVisitor_H

#include <memory>

#include "SiconosException.hpp"
#include "SiconosVector.hpp"

// For any DS class that might be visited, visit functions must be declared here.
// Even those for DS in other components like mechanics.
//

#ifdef HAVE_SICONOS_MECHANICS
namespace siconos::collision {

namespace native::bodies {
class DiskDiskR;
class CircleCircleR;
class DiskMovingPlanR;
class DiskPlanR;
class SphereLDSSphereLDSR;
class SphereNEDSSphereNEDSR;
class SphereLDSPlanR;
class SphereNEDSPlanR;
}  // namespace native::bodies

class ContactR;
class Contact5DR;
class Contact2d3DR;
class Contact2dR;

#ifdef SICONOS_HAS_BULLET
namespace bullet {
class BulletR;
}
#endif

}  // namespace siconos::collision

namespace siconos::joints {
class CouplerJointR;
class CylindricalJointR;
class FixedJointR;
class JointFrictionR;
class JointStopR;
class KneeJointR;
class NewtonEulerJointR;
class PivotJointR;
class PrismaticJointR;
}  // namespace siconos::joints

#ifdef SICONOS_HAS_OpenCASCADE
namespace siconos::mechanics::occ {
class OccBody;
}
#endif
#endif

namespace siconos::modeling {
class Relation;
class LagrangianScleronomousR;

namespace relations {

struct Visitor {
  virtual ~Visitor() noexcept = default;

  // Define this function in your visitor if you need to handle
  // internal memory (view) without copy. See example in MechanicsIO.cpp
  virtual void setMap(Eigen::Ref<siconos::algebra::SiconosVector>){};

  virtual void visit(std::shared_ptr<siconos::modeling::Relation>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "Relation in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::modeling::Relation &) {
    THROW_EXCEPTION(
        "you must define a visit function for Relation in "
        "a derived class of relations::Visitor");
  }

  virtual void visit(const siconos::modeling::Relation &) {
    THROW_EXCEPTION(
        "you must define a visit function for Relation in "
        "a derived class of relations::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::modeling::LagrangianScleronomousR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: LagrangianScleronomousR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::modeling::LagrangianScleronomousR &) {
    THROW_EXCEPTION(
        "you must define a visit function for LagrangianScleronomousR in a "
        "derived class of relations::Visitor");
  }
  virtual void visit(const siconos::modeling::LagrangianScleronomousR &) {
    THROW_EXCEPTION(
        "you must define a visit function for LagrangianScleronomousR in a "
        "derived class of relations::Visitor");
  }

#ifdef HAVE_SICONOS_MECHANICS
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::DiskDiskR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "DiskDiskR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::DiskDiskR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskDiskR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::DiskDiskR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskDiskR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::CircleCircleR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "CircleCircleR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::CircleCircleR &) {
    THROW_EXCEPTION(
        "you must define a visit function for CircleCircleR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::CircleCircleR &) {
    THROW_EXCEPTION(
        "you must define a visit function for CircleCircleR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::DiskMovingPlanR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "DiskMovingPlanR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::DiskMovingPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskMovingPlanR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::DiskMovingPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskMovingPlanR "
        "in a derived class of relations::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::DiskPlanR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to DiskPlanR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::collision::native::bodies::DiskPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskPlanR in a "
        "derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::DiskPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskPlanR in a "
        "derived class of relations::Visitor");
  }

  virtual void visit(
      std::shared_ptr<siconos::collision::native::bodies::SphereLDSSphereLDSR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "SphereLDSSphereLDSR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereLDSSphereLDSR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDSSphereLDSR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereLDSSphereLDSR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDSSphereLDSR "
        "in a derived class of relations::Visitor");
  }

  virtual void visit(
      std::shared_ptr<siconos::collision::native::bodies::SphereNEDSSphereNEDSR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "SphereNEDSSphereNEDSR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereNEDSSphereNEDSR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDSSphereNEDSR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereNEDSSphereNEDSR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDSSphereNEDSR "
        "in a derived class of relations::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::SphereLDSPlanR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "SphereLDSPlanR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereLDSPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDSPlanR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereLDSPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDSPlanR "
        "in a derived class of relations::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::SphereNEDSPlanR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "SphereNEDSPlanR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereNEDSPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDSPlanR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereNEDSPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDSPlanR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::ContactR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "ContactR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::ContactR &) {
    THROW_EXCEPTION(
        "you must define a visit function for ContactR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::ContactR &) {
    THROW_EXCEPTION(
        "you must define a visit function for ContactR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::Contact5DR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "Contact5DR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::Contact5DR &) {
    THROW_EXCEPTION(
        "you must define a visit function for Contact5DR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::Contact5DR &) {
    THROW_EXCEPTION(
        "you must define a visit function for Contact5DR "
        "in a derived class of relations::Visitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::Contact2dR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "Contact2dR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::Contact2dR &) {
    THROW_EXCEPTION(
        "you must define a visit function for Contact2dR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::Contact2dR &) {
    THROW_EXCEPTION(
        "you must define a visit function for Contact2dR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::Contact2d3DR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "Contact2d3DR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::Contact2d3DR &) {
    THROW_EXCEPTION(
        "you must define a visit function for Contact2d3DR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::Contact2d3DR &) {
    THROW_EXCEPTION(
        "you must define a visit function for Contact2d3DR "
        "in a derived class of relations::Visitor");
  }
#ifdef SICONOS_HAS_BULLET

  virtual void visit(std::shared_ptr<siconos::collision::bullet::BulletR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "BulletR in a derived class of relations::Visitor");
  }

  virtual void visit(siconos::collision::bullet::BulletR &) {
    THROW_EXCEPTION(
        "you must define a visit function for BulletR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::collision::bullet::BulletR &) {
    THROW_EXCEPTION(
        "you must define a visit function for BulletR "
        "in a derived class of relations::Visitor");
  }
#endif
  virtual void visit(std::shared_ptr<siconos::joints::CouplerJointR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "CouplerJointR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::CouplerJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for CouplerJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::CouplerJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for CouplerJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::joints::CylindricalJointR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "CylindricalJointR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::CylindricalJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for CylindricalJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::CylindricalJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for CylindricalJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::joints::FixedJointR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "FixedJointR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::FixedJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for FixedJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::FixedJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for FixedJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::joints::JointFrictionR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "JointFrictionR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::JointFrictionR &) {
    THROW_EXCEPTION(
        "you must define a visit function for JointFrictionR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::JointFrictionR &) {
    THROW_EXCEPTION(
        "you must define a visit function for JointFrictionR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::joints::JointStopR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "JointStopR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::JointStopR &) {
    THROW_EXCEPTION(
        "you must define a visit function for JointStopR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::JointStopR &) {
    THROW_EXCEPTION(
        "you must define a visit function for JointStopR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::joints::KneeJointR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "KneeJointR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::KneeJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for KneeJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::KneeJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for KneeJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::joints::NewtonEulerJointR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "NewtonEulerJointR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::NewtonEulerJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for NewtonEulerJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::NewtonEulerJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for NewtonEulerJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::joints::PivotJointR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "PivotJointR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::PivotJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for PivotJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::PivotJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for PivotJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(std::shared_ptr<siconos::joints::PrismaticJointR>) {
    THROW_EXCEPTION(
        "you must define a visit function for shared ptr to "
        "PrismaticJointR in a derived class of relations::Visitor");
  }
  virtual void visit(siconos::joints::PrismaticJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for PrismaticJointR "
        "in a derived class of relations::Visitor");
  }
  virtual void visit(const siconos::joints::PrismaticJointR &) {
    THROW_EXCEPTION(
        "you must define a visit function for PrismaticJointR "
        "in a derived class of relations::Visitor");
  }

#endif
};
}  // namespace relations
}  // namespace siconos::modeling
#endif