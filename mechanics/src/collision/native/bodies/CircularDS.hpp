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

  double radius{0.};
  double massValue{0.};

 public:
  CircularDS(double, double, std::shared_ptr<siconos::algebra::SiconosVector>,
             std::shared_ptr<siconos::algebra::SiconosVector>);

  virtual ~CircularDS() noexcept = default;

  double getQ(unsigned int pos);

  double getVelocity(unsigned int pos);

  inline double getMassValue() const { return massValue; };

  inline double getRadius() const { return radius; };
};
}  // namespace siconos::collision::native::bodies
#endif /* CircularDS_h */
