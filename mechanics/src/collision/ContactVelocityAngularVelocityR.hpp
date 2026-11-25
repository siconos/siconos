// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef ContactVelocityAngularVelocityR_hpp
#define ContactVelocityAngularVelocityR_hpp

#include "NewtonEulerVelocityAngularVelocityR.hpp"

namespace siconos::collision {
class BodyShapeRecord;

class ContactVelocityAngularVelocityR
    : public siconos::modeling::NewtonEulerVelocityAngularVelocityR {
 private:
  ACCEPT_SERIALIZATION(ContactVelocityAngularVelocityR);

 public:
  /* For users that may require extra information about contacts. */
  std::shared_ptr<siconos::collision::BodyShapeRecord> bodyShapeRecordA{nullptr};
  std::shared_ptr<siconos::collision::BodyShapeRecord> bodyShapeRecordB{nullptr};

  /**
     to compute the output y = h(q) of the Relation

      \param[in] q1 generalized coordinates vector of the fist dynamical system involved
      in the relation
      \param[in] q2 generalized coordinates vector of the second dynamical system
      involved in the relation
      \param[in,out] y the resulting vector
  */
  virtual void computeh(
      const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
      const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
      Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Update this contact point information.
   *
   *  \param pos1 Position on ds1 in ds1 frame.
   *  \param pos2 Position on ds2 in ds2 frame (or world frame if ds2=null).
   *  \param normal Normal in ds2 frame (or world frame if ds2=null).
   */
  virtual void updateContactPoints(const siconos::algebra::SiconosVector3& pos1,
                                   const siconos::algebra::SiconosVector3& pos2,
                                   const siconos::algebra::SiconosVector3& normal);

  virtual void preDelete() {}
  virtual void accept(modeling::relations::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision
#endif
