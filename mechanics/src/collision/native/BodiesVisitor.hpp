/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file BodiesVisitor.hpp
  \brief interface to visitors for siconos bodies
*/
#ifndef BodiesVisitor_hpp
#define BodiesVisitor_hpp

#include "BodiesAndRelationsDeclaration.hpp"
#include "SiconosVisitor.hpp"

namespace siconos::modeling {
class LagrangianScleronomousR;
}

namespace siconos::collision::native {
struct BodiesVisitor : public siconos::internal::SiconosVisitor {
  BodiesVisitor() = default;
  BodiesVisitor(const BodiesVisitor &) = delete;
  BodiesVisitor(BodiesVisitor &&) = delete;
  BodiesVisitor &operator=(BodiesVisitor &&) = delete;
  BodiesVisitor &operator=(const BodiesVisitor &) = delete;

  using siconos::internal::SiconosVisitor::visit;

  // Dynamical Systems
  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::Circle>)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::Disk>)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::SphereLDS>)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::SphereNEDS>)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(std::shared_ptr<siconos::collision::native::bodies::ExternalBody>)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::collision::native::bodies::DiskDiskR &)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::collision::native::bodies::CircleCircleR &)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::collision::native::bodies::DiskMovingPlanR &)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::modeling::LagrangianScleronomousR &)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::collision::native::bodies::DiskPlanR &rel)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::collision::native::bodies::SphereLDSSphereLDSR &)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::collision::native::bodies::SphereNEDSSphereNEDSR &)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::collision::native::bodies::SphereLDSPlanR &rel)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual void visit(const siconos::collision::native::bodies::SphereNEDSPlanR &rel)
  {
    THROW_EXCEPTION("you must define a visit function in a derived class of BodiesVisitor");
  };

  virtual ~BodiesVisitor() noexcept = default;
};
}  // namespace siconos::collision::native
#endif
