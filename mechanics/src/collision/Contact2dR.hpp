// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef Contact2dR_hpp
#define Contact2dR_hpp

#include "Lagrangian2dVelocityR.hpp"

namespace siconos::collision {
class BodyShapeRecord;

class Contact2dR : public siconos::modeling::Lagrangian2dVelocityR {
 private:
  ACCEPT_SERIALIZATION(Contact2dR);

 public:
  /**/
  virtual ~Contact2dR() noexcept = default;

  /* For users that may require extra information about contacts. */
  std::shared_ptr<siconos::collision::BodyShapeRecord> bodyShapeRecordA{nullptr};
  std::shared_ptr<siconos::collision::BodyShapeRecord> bodyShapeRecordB{nullptr};

  /** Update this contact point information.
   *
   *  \param pos1 Position on ds1 in ds1 frame.
   *  \param pos2 Position on ds2 in ds2 frame (or world frame if ds2=null).
   *  \param normal Normal in ds2 frame (or world frame if ds2=null).
   */
  // virtual void updateContactPoints(const siconos::algebra::SiconosVector& pos1,
  //                                  const siconos::algebra::SiconosVector& pos2,
  //                                  const siconos::algebra::SiconosVector& normal);
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
