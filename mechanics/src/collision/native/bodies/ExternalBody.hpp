// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file ExternalBody.hpp
  \brief External visitable lagrangian class.
*/
#ifndef ExternalBody_hpp
#define ExternalBody_hpp

#include <LagrangianDS.hpp>

namespace siconos::collision::native {
class SpaceFilter;
}

namespace siconos::collision::native::bodies {

class ExternalBody : public siconos::modeling::LagrangianDS,
                     public std::enable_shared_from_this<ExternalBody> {
 public:
  virtual void selfHash(siconos::collision::native::SpaceFilter&) = 0;

  virtual void selfFindInteractions(
      std::shared_ptr<siconos::collision::native::SpaceFilter>) = 0;

  virtual void acceptSP(modeling::dynamical_systems::Visitor& tourist) override {
    tourist.visit(shared_from_this());
  }
  virtual void accept(modeling::dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }

 protected:
  ACCEPT_SERIALIZATION(ExternalBody);
};
}  // namespace siconos::collision::native::bodies
#endif
