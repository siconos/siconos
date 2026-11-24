// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SphereNEDS.hpp

  \brief Definition of a 3D Sphere as a NewtonEulerDS (with
  quaternions).

*/
#ifndef SphereNEDS_h
#define SphereNEDS_h

#include "NewtonEulerDS.hpp"

namespace siconos::collision::native::bodies {
class SphereNEDS : public siconos::modeling::NewtonEulerDS,
                   public std::enable_shared_from_this<SphereNEDS> {
 protected:
  ACCEPT_SERIALIZATION(SphereNEDS);

  double radius{0.};

 public:
  SphereNEDS(double r, double m, Eigen::Ref<siconos::algebra::SiconosMatrix33> inertia,
             Eigen::Ref<siconos::algebra::SiconosVector7> position,
             Eigen::Ref<siconos::algebra::SiconosVector6> twist, siconos::algebra::AliasTag);

  SphereNEDS(double r, double m, const siconos::algebra::SiconosMatrix33& inertia,
             const siconos::algebra::SiconosVector7& position,
             const siconos::algebra::SiconosVector6& twist, siconos::algebra::CopyTag);

  ~SphereNEDS() noexcept = default;

  double getQ(siconos::algebra::Index pos);

  double getTwist(siconos::algebra::Index pos);

  inline double getRadius() const { return radius; };
  virtual void acceptSP(modeling::dynamical_systems::Visitor& tourist) override {
    tourist.visit(shared_from_this());
  }
  virtual void accept(modeling::dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};
}  // namespace siconos::collision::native::bodies
#endif /* SphereNEDS_h */
