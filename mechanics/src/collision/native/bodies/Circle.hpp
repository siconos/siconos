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

/*! \file Circle.hpp
 */
/** \class Circle
 * \brief Definition of a 2D Circle - Inherits from CircularDS
 */

#ifndef Circle_H
#define Circle_H

#include "CircularDS.hpp"

namespace siconos::collision::native::bodies {

class Circle : public CircularDS, public std::enable_shared_from_this<Circle> {
 private:
  ACCEPT_SERIALIZATION(Circle);

 public:
  /** constructor from initial state and velocity
   *
   *  \param R radius
   *  \param m mass
   *  \param position initial coordinates
   *  \param velocity initial velocity
   */
  Circle(double radius, double mass, Eigen::Ref<siconos::algebra::SiconosVector> position,
         Eigen::Ref<siconos::algebra::SiconosVector> velocity);

  /** destructor */
  virtual ~Circle() noexcept = default;
};
}  // namespace siconos::collision::native::bodies
#endif /* Circle_H */
