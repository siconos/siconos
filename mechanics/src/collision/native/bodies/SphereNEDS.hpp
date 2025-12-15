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
