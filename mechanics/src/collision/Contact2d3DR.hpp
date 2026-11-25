// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef Contact2d3DR_hpp
#define Contact2d3DR_hpp

#include "Lagrangian2dVelocityAngularVelocityR.hpp"

namespace siconos::collision {
class BodyShapeRecord;

class Contact2d3DR : public siconos::modeling::Lagrangian2dVelocityAngularVelocityR {
 private:
  ACCEPT_SERIALIZATION(Contact2d3DR);

 public:
  virtual ~Contact2d3DR() noexcept = default;

  /* For users that may require extra information about contacts. */
  std::shared_ptr<siconos::collision::BodyShapeRecord> bodyShapeRecordA{nullptr};
  std::shared_ptr<siconos::collision::BodyShapeRecord> bodyShapeRecordB{nullptr};

  /** Update this contact point information.
   *
   *  \param pos1 Position on ds1 in ds1 frame.
   *  \param pos2 Position on ds2 in ds2 frame (or world frame if ds2=null).
   *  \param normal Normal in ds2 frame (or world frame if ds2=null).
   */
  virtual void updateContactPointsInAbsoluteFrame(
      const siconos::algebra::SiconosVector2& pos1,
      const siconos::algebra::SiconosVector2& pos2,
      const siconos::algebra::SiconosVector2& normal);

  virtual void preDelete() {}
  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision
#endif
