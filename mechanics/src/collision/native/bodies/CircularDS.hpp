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

/** \file CircularDS.hpp
 */

/**
   Definition of a 2D circular shape - Inherits from LagrangianDS
*/

#ifndef CircularDS_h
#define CircularDS_h

#include "LagrangianDS.hpp"

namespace siconos::collision::native::bodies {
class CircularDS : public siconos::modeling::LagrangianDS {
 protected:
  ACCEPT_SERIALIZATION(CircularDS);

  double radius_{0.};
  double massValue_{0.};

 public:
  /** constructor from initial state and velocity
   *
   *  \param R radius
   *  \param m mass
   *  \param position initial coordinates
   *  \param velocity initial velocity
   */
  CircularDS(double, double, Eigen::Ref<siconos::algebra::SiconosVector> position,
             Eigen::Ref<siconos::algebra::SiconosVector> velocity);

  virtual ~CircularDS() noexcept = default;

  double getQ(siconos::algebra::Index pos);

  double getVelocity(siconos::algebra::Index pos);

  inline double getMassValue() const { return massValue_; };

  inline double getRadius() const { return radius_; };
};
}  // namespace siconos::collision::native::bodies
#endif /* CircularDS_h */
