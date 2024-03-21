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

/*! \file SiconosCollisionQueryResult.hpp
\brief Holds one result of a query against the graph of body
contactors maintained by a SiconosCollisionManager.
*/

#ifndef SiconosCollisionQueryResult_h
#define SiconosCollisionQueryResult_h

#include "SiconosVector.hpp"
#include <memory>

#include "SiconosSerialization.hpp"


namespace siconos::modeling {
class SecondOrderDS;
}

namespace siconos::collision {

class SiconosShape;
class SiconosContactor;

/**
   Holds one result of a line segment intersection query
   against the graph of body contactors maintained by a
   SiconosCollisionManager
*/
class SiconosCollisionQueryResult {
 protected:
  ACCEPT_SERIALIZATION(SiconosCollisionQueryResult);

 public:
  /** Distance from reference point (start of line segment or query center) */
  double distance{0.};

  /** Body owning the contactor that was intersected, may be null for
   *  static contactors. */
  std::shared_ptr<siconos::modeling::SecondOrderDS> body{nullptr};

  /** The shape that was intersected. */
  std::shared_ptr<SiconosShape> shape{nullptr};

  /** The contactor that was intersected. */
  std::shared_ptr<SiconosContactor> contactor{nullptr};

  /** Closest point on contactor in world coordinates. */
  std::shared_ptr<siconos::algebra::SiconosVector> point{nullptr};
};

}  // namespace siconos::collision
#endif /* SiconosCollisionQueryResult.hpp */
