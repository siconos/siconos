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

/*! \file SphereLDS.hpp
  \brief Definition of a 3D sphere as a LagrangianDS (with Euler
         Angles)
*/

#ifndef SphereLDS_h
#define SphereLDS_h

#include "LagrangianDS.hpp"

namespace siconos::collision::native::bodies {

class SphereLDS : public siconos::modeling::LagrangianDS,
                  public std::enable_shared_from_this<SphereLDS> {
 protected:
  ACCEPT_SERIALIZATION(SphereLDS);

  double radius_{0.};
  double massValue_{0.};
  double inertia_{0.};

 public:
  /** constructor from initial state and velocity
   *
   *  \param R radius of the sphere
   *  \param m mass of the sphere
   *  \param position initial coordinates
   *  \param velocity initial velocity
   */

  SphereLDS(double, double, const siconos::algebra::SiconosVector6& position,
            const siconos::algebra::SiconosVector6& velocity);

  ~SphereLDS() noexcept = default;

  /** to compute the mass matrix operator \f$ M(q) \f$
   *
   *  \param position q vector
   */
  void computeMass(const Eigen::Ref<const siconos::algebra::SiconosVector>& position) override;

  double getQ(unsigned int pos);

  double getVelocity(unsigned int pos);

  inline double getMassValue() const { return massValue_; };

  inline double getRadius() const { return radius_; };
  virtual void acceptSP(modeling::dynamical_systems::Visitor& tourist) override {
    tourist.visit(shared_from_this());
  }
  virtual void accept(modeling::dynamical_systems::Visitor& tourist) const override {
    tourist.visit(*this);
  }
};

}  // namespace siconos::collision::native::bodies
#endif /* SphereLDS_h */
