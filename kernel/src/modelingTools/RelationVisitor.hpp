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

#include "SiconosConfig.h"  // For HAVE_SICONOS ...
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
}  // namespace siconos::collision

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
        "you must define a visit function for SP :: "
        "Relation in a derived class of SiconosVisitor");
  }
  virtual void visit(siconos::modeling::Relation &) {
    THROW_EXCEPTION(
        "you must define a visit function for Relation in "
        "a derived class of SiconosVisitor");
  }

  virtual void visit(const siconos::modeling::Relation &) {
    THROW_EXCEPTION(
        "you must define a visit function for Relation in "
        "a derived class of SiconosVisitor");
  }

  virtual void visit(std::shared_ptr<siconos::modeling::LagrangianScleronomousR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: LagrangianScleronomousR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(siconos::modeling::LagrangianScleronomousR &) {
    THROW_EXCEPTION(
        "you must define a visit function for LagrangianScleronomousR in a "
        "derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::modeling::LagrangianScleronomousR &) {
    THROW_EXCEPTION(
        "you must define a visit function for LagrangianScleronomousR in a "
        "derived class of SiconosVisitor");
  }

#ifdef HAVE_SICONOS_MECHANICS
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::DiskDiskR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: "
        "DiskDiskR in a derived class of SiconosVisitor");
  }

  virtual void visit(siconos::collision::native::bodies::DiskDiskR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskDiskR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::collision::native::bodies::DiskDiskR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskDiskR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::CircleCircleR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: "
        "CircleCircleR in a derived class of SiconosVisitor");
  }

  virtual void visit(siconos::collision::native::bodies::CircleCircleR &) {
    THROW_EXCEPTION(
        "you must define a visit function for CircleCircleR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::collision::native::bodies::CircleCircleR &) {
    THROW_EXCEPTION(
        "you must define a visit function for CircleCircleR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::DiskMovingPlanR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: "
        "DiskMovingPlanR in a derived class of SiconosVisitor");
  }

  virtual void visit(siconos::collision::native::bodies::DiskMovingPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskMovingPlanR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::collision::native::bodies::DiskMovingPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskMovingPlanR "
        "in a derived class of SiconosVisitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::DiskPlanR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: DiskPlanR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(siconos::collision::native::bodies::DiskPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskPlanR in a "
        "derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::collision::native::bodies::DiskPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for DiskPlanR in a "
        "derived class of SiconosVisitor");
  }

  virtual void visit(
      std::shared_ptr<siconos::collision::native::bodies::SphereLDSSphereLDSR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: "
        "SphereLDSSphereLDSR in a derived class of SiconosVisitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereLDSSphereLDSR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDSSphereLDSR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereLDSSphereLDSR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDSSphereLDSR "
        "in a derived class of SiconosVisitor");
  }

  virtual void visit(
      std::shared_ptr<siconos::collision::native::bodies::SphereNEDSSphereNEDSR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: "
        "SphereNEDSSphereNEDSR in a derived class of SiconosVisitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereNEDSSphereNEDSR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDSSphereNEDSR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereNEDSSphereNEDSR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDSSphereNEDSR "
        "in a derived class of SiconosVisitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::SphereLDSPlanR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: "
        "SphereLDSPlanR in a derived class of SiconosVisitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereLDSPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDSPlanR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereLDSPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereLDSPlanR "
        "in a derived class of SiconosVisitor");
  }

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::SphereNEDSPlanR>) {
    THROW_EXCEPTION(
        "you must define a visit function for SP :: "
        "SphereNEDSPlanR in a derived class of SiconosVisitor");
  }

  virtual void visit(siconos::collision::native::bodies::SphereNEDSPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDSPlanR "
        "in a derived class of SiconosVisitor");
  }
  virtual void visit(const siconos::collision::native::bodies::SphereNEDSPlanR &) {
    THROW_EXCEPTION(
        "you must define a visit function for SphereNEDSPlanR "
        "in a derived class of SiconosVisitor");
  }
#endif
};
}  // namespace relations
}  // namespace siconos::modeling
#endif